"""
3D real-time renderer using vispy.

Every pixel on screen is driven by per-particle state from the simulator —
nothing is painted on the background. The void is empty (as it is in real
space; the bright fixed stars you see at night are all *within* the volume
a cosmological sim of this scale would inhabit, so we don't paint extras
in for decoration).

Rendering pipeline
------------------
For each live particle we render:

1. **Core marker** at the particle's position.
   - In CPK mode (default): coloured by element identity (chemistry convention).
   - In astronomy mode (P key): coloured by the Planck blackbody RGB at the
     particle's effective temperature.
   - Size: covalent radius for atoms; cube-root-of-mass scaling for accreted
     bodies (real physical body radius).

2. **Inner halo** (additive blend) — represents the close-in radiative output.
   Radius is driven by Stefan-Boltzmann: apparent r ∝ √L ∝ T².

3. **Outer halo** (additive blend) — the dim radiative wings, 2.2× the inner
   halo. Where many halos overlap, the additive blending produces the
   volumetric brightening that real astronomy long-exposures show.

Halo *colour* is always the real Planck blackbody RGB at the particle's
temperature — radiation has no element identity. Halo *intensity* scales
with the same T² Stefan-Boltzmann factor.

Effective temperature in Kelvin is computed per particle as:
  - Free atoms / gas: kinetic-energy proxy (T ∝ KE), calibrated against
    a reference (hot_temperature_threshold ↔ reference_temperature_K).
  - Accreted bodies (planets, stars): stellar mass-luminosity relation,
    T ∝ √M (Eddington), calibrated so a body at the star_mass_threshold
    sits at the reference temperature.

Covalent bonds are drawn as line segments.
A minimal HUD in the window title shows live stats.

Controls (built into vispy TurntableCamera):
  Left-drag   → orbit
  Right-drag  → zoom
  Middle-drag → pan
  Scroll      → zoom
  Space       → pause / resume
  +/-         → speed up / slow down sim
  R           → reset camera (cancels cinematic mode)
  F           → cinematic camera (auto orbit + breathing elevation)
  P           → toggle astronomy palette (stellar colour-temperature)
  Q / Escape  → quit
  Left-click  → inspect particle (click again to deselect)
"""

from __future__ import annotations
import time
from collections import Counter
import numpy as np

from vispy import app, scene
from vispy.scene import visuals
from vispy.visuals.transforms import STTransform

from sim.diagnostics import kinetic_energy, gravitational_pe, linear_momentum
from sim.elements import ELEMENTS_LIST
from sim.extinction import compute_visibility, per_particle_optical_depth
from sim.illumination import compute_received_rgb_flux
from sim.molecules import identify_molecules, molecule_counts
from sim.polymers import polymer_stats, polymer_of
from sim.states import state_counts
from sim.world import World

_FPS_TARGET    = 60
_FPS_WINDOW    = 0.5     # smoothing window for FPS measurement (seconds)
_SPEED_UP      = 2.0     # speed multiplier per keypress
_SPEED_MAX     = 64.0
_SPEED_MIN     = 0.0625  # 1/16×
_CAM_ORBIT_STEP = 5.0    # degrees per arrow-key press
_CAM_ZOOM_STEP  = 1.15   # zoom factor per Page Up/Down press
_CAM_AZIMUTH_0  = 30.0   # default camera azimuth  (degrees)
_CAM_ELEVATION_0 = 20.0  # default camera elevation (degrees)
_PICK_PIXEL_RADIUS = 15.0  # click tolerance for particle picking
_ENERGY_REFRESH_FRAMES = 30   # how often to recompute O(N²) gravitational PE

# CPK colours for rendering classification
_COLOUR_PLANET = np.array([0.55, 0.50, 0.42], dtype=np.float32)   # rocky grey-brown
_COLOUR_STAR   = np.array([1.00, 0.92, 0.65], dtype=np.float32)   # warm yellow-white
_COLOUR_HOT    = np.array([1.00, 0.40, 0.10], dtype=np.float32)   # heat tint blended in
_COLOUR_PLASMA = np.array([0.70, 0.85, 1.00], dtype=np.float32)   # ionised: bluish-white

# ---------------------------------------------------------------------------
# Rendering pipeline — every visual element is driven by per-particle physical
# state. The halo around each particle represents its thermal radiation:
#
#   Luminosity ∝ T⁴       (Stefan-Boltzmann)
#   Apparent halo radius ∝ √Luminosity ∝ T²
#
# Where temperature T (Kelvin) is the kinetic temperature of the particle,
# computed from its kinetic energy and the eV→Kelvin calibration of the
# simulator's custom time units. The halo *colour* uses the real Planck-curve
# blackbody RGB at that temperature (Tanner Helland's polynomial fit to the
# CIE colour matching functions).
#
# Nothing is painted on the background. The void is empty.
# ---------------------------------------------------------------------------

# Halo size scaling — base halo is a small fixed fraction of the particle
# size; the bulk of the halo extent is driven by the T² Stefan-Boltzmann
# perceptual-radius scaling.
_HALO_BASE_SIZE_MUL  = 1.8     # cool particle halo as multiple of core size
_HALO_RADIATION_GAIN = 6.0     # extra halo radius per (T/T_ref)² unit
_HALO_INNER_ALPHA    = 0.45    # inner halo opacity (cool floor)
_HALO_OUTER_ALPHA    = 0.10    # outer halo opacity (cool floor)
_BODY_HALO_MUL       = 4.0     # accreted bodies get bigger halos (more emitting area)

# Cinematic camera
_CINEMATIC_ORBIT_DEG_PER_SEC = 6.0
_CINEMATIC_ELEV_AMPLITUDE    = 12.0
_CINEMATIC_ELEV_PERIOD_SEC   = 35.0


class Viewer:
    def __init__(self, world: World, cfg):
        self.world  = world
        self.cfg    = cfg
        self.paused = False
        self.speed  = 1.0   # time-scale multiplier

        # ---- canvas -------------------------------------------------------
        # Pure black background: a darker void makes additive glow + the
        # starfield pop with higher dynamic range than the original deep-blue.
        bg = cfg.renderer.background
        self.canvas = scene.SceneCanvas(
            title='Universe Simulator',
            size=(cfg.renderer.width, cfg.renderer.height),
            bgcolor=(0.005, 0.005, 0.01, 1.0),
            keys='interactive',
            show=True,
        )
        self.view = self.canvas.central_widget.add_view()
        self.view.camera           = 'turntable'
        self.view.camera.fov       = cfg.renderer.fov
        self.view.camera.distance  = cfg.renderer.camera_distance
        self.view.camera.azimuth   = _CAM_AZIMUTH_0
        self.view.camera.elevation = _CAM_ELEVATION_0

        # ---- radiative halos --------------------------------------------
        # Two superimposed Markers visuals render each particle's thermal
        # emission, sized by Stefan-Boltzmann (apparent radius ∝ T²) and
        # additively blended so overlapping halos brighten the scene the
        # way they would in a real long-exposure photo. No fake stars are
        # added — the void is empty.
        self.markers_nebula = visuals.Markers(parent=self.view.scene)
        self.markers_nebula.antialias = 1
        self.markers_nebula.set_gl_state(
            'translucent', blend=True,
            blend_func=('src_alpha', 'one'),
            depth_test=False,
        )
        self.markers_nebula.order = -2

        self.markers_inner_glow = visuals.Markers(parent=self.view.scene)
        self.markers_inner_glow.antialias = 1
        self.markers_inner_glow.set_gl_state(
            'translucent', blend=True,
            blend_func=('src_alpha', 'one'),
            depth_test=False,
        )
        self.markers_inner_glow.order = -1

        # ---- Phase-3 atmosphere halo ------------------------------------
        # For accreted bodies that have absorbed significant light-element
        # content (H, He), an extra halo is rendered between the body and
        # the radiative inner glow.  Colour is the mass-weighted CPK of the
        # body's light constituents — this is the gravitationally-bound gas
        # envelope, distinct from the body's thermal radiation.
        self.markers_atmosphere = visuals.Markers(parent=self.view.scene)
        self.markers_atmosphere.antialias = 1
        self.markers_atmosphere.set_gl_state(
            'translucent', blend=True,
            blend_func=('src_alpha', 'one'),
            depth_test=False,
        )
        self.markers_atmosphere.order = -3

        # ---- Phase-3 dense core (heavy-element fraction) ----------------
        # For bodies with substantial heavy-element content (Fe, Si, Ti, …),
        # a smaller, more saturated marker is rendered ON TOP of the body
        # to visualise the dense core — heavy elements naturally settle
        # inward in a differentiated body (real gravitational stratification).
        # The atmosphere halo + body + core layers together produce the
        # crust-mantle-core appearance of a real planet.
        self.markers_core = visuals.Markers(parent=self.view.scene)
        self.markers_core.antialias = 1
        self.markers_core.order = 1                # in front of the body marker

        # ---- primary markers (opaque cores) ------------------------------
        self.markers = visuals.Markers(parent=self.view.scene)
        self.markers.antialias = 1

        self.lines = visuals.Line(
            parent=self.view.scene,
            method='gl',
            connect='segments',
        )

        # ---- timer --------------------------------------------------------
        self._timer = app.Timer(
            interval=1.0 / _FPS_TARGET,
            connect=self._on_timer,
            start=True,
        )

        # ---- selection / inspect ------------------------------------------
        self._selected: int | None = None

        self._select_marker = visuals.Markers(parent=self.view.scene)
        self._select_marker.antialias = 0

        self._info_text = visuals.Text(
            '', color=(1.0, 1.0, 1.0, 0.9), font_size=9,
            parent=self.canvas.scene, anchor_x='left', anchor_y='top',
        )
        self._info_text.transform = STTransform(translate=(10, 60))

        # ---- key bindings -------------------------------------------------
        self.canvas.events.key_press.connect(self._on_key)
        self.canvas.events.mouse_press.connect(self._on_mouse_press)

        # ---- perf tracking ------------------------------------------------
        self._last_wall  = time.perf_counter()
        self._fps_acc    = 0.0
        self._fps_frames = 0
        self._fps        = 0.0

        # ---- conservation diagnostics -------------------------------------
        # Recompute the O(N²) PE every _ENERGY_REFRESH_FRAMES; cache between.
        self._ke           = 0.0
        self._pe           = 0.0
        self._momentum_mag = 0.0
        self._energy_ref   = None     # set once we have ≥2 particles
        self._energy_frame = 0

        # ---- cinematic mode (auto camera) ---------------------------------
        self._cinematic_mode: bool = False
        self._cinematic_t:    float = 0.0
        self._astronomy_palette: bool = False     # 'P' to toggle

    # ------------------------------------------------------------------
    # Timer callback — drives simulation + render each frame
    # ------------------------------------------------------------------

    def _on_timer(self, event) -> None:
        now  = time.perf_counter()
        wall = now - self._last_wall
        self._last_wall = now

        self._fps_acc    += wall
        self._fps_frames += 1
        if self._fps_acc >= _FPS_WINDOW:
            self._fps        = self._fps_frames / self._fps_acc
            self._fps_acc    = 0.0
            self._fps_frames = 0

        if not self.paused:
            dt      = self.cfg.simulation.time_step * self.speed
            n_steps = self.cfg.simulation.steps_per_frame
            for _ in range(n_steps):
                self.world.step(dt)

        # Cinematic camera — slow orbit + breathing elevation
        if self._cinematic_mode:
            self._cinematic_t += wall
            cam = self.view.camera
            cam.azimuth   = (_CAM_AZIMUTH_0 + _CINEMATIC_ORBIT_DEG_PER_SEC * self._cinematic_t) % 360.0
            cam.elevation = _CAM_ELEVATION_0 + _CINEMATIC_ELEV_AMPLITUDE * np.sin(
                2.0 * np.pi * self._cinematic_t / _CINEMATIC_ELEV_PERIOD_SEC
            )

        self._update_visuals()
        self._update_title()

    # ------------------------------------------------------------------
    # Update vispy visuals from world state
    # ------------------------------------------------------------------

    def _update_visuals(self) -> None:
        w   = self.world
        n   = w.n
        cfg = self.cfg.renderer

        if n == 0:
            self._selected = None
            zeros = np.zeros((1, 3), dtype=np.float32)
            transparent = (0.0, 0.0, 0.0, 0.0)
            self.markers.set_data(pos=zeros, face_color=transparent)
            self.markers_inner_glow.set_data(pos=zeros, face_color=transparent, size=1, edge_width=0)
            self.markers_nebula.set_data(pos=zeros, face_color=transparent, size=1, edge_width=0)
            self.markers_atmosphere.set_data(pos=zeros, face_color=transparent, size=1, edge_width=0)
            self.markers_core.set_data(pos=zeros, face_color=transparent, size=1, edge_width=0)
            self._select_marker.set_data(pos=zeros, face_color=transparent, size=1, edge_width=0)
            self._info_text.text = ''
            self.lines.set_data(pos=np.zeros((2, 3)))
            return

        if self._selected is not None and self._selected >= n:
            self._selected = None

        pos   = w.positions[:n].copy()
        e_ids = w.elem_ids[:n]

        # --- Base CPK colours and covalent-radius sizes --------------------
        colors = np.ones((n, 4), dtype=np.float32)
        sizes  = np.empty(n, dtype=np.float32)
        for idx in range(n):
            elem           = ELEMENTS_LIST[e_ids[idx]]
            colors[idx, :3] = elem.color
            sizes[idx]      = cfg.particle_size_base + cfg.particle_size_scale * elem.covalent_radius

        # --- Kinetic energy per particle (also used for halo intensity) ---
        ke    = 0.5 * w.masses[:n] * np.sum(w.velocities[:n] ** 2, axis=1)
        hot_T = float(getattr(cfg, 'hot_temperature_threshold', 0.0))

        # --- Temperature tint (cold → CPK, hot → orange-red) --------------
        if hot_T > 0.0:
            heat = np.clip(ke / hot_T, 0.0, 1.0).astype(np.float32)[:, np.newaxis]
            colors[:, :3] = colors[:, :3] * (1.0 - heat) + _COLOUR_HOT * heat

        # --- Plasma override (ionised particles) ---------------------------
        ion_mask = w.ionized[:n]
        if ion_mask.any():
            colors[ion_mask, :3] = _COLOUR_PLASMA

        # --- Override for accreted bodies (planets / stars) ----------------
        rend_cfg = cfg
        planet_threshold = float(getattr(rend_cfg, 'planet_mass_threshold', 1e9))
        star_threshold   = float(getattr(rend_cfg, 'star_mass_threshold',   1e9))
        body_size_scale  = float(getattr(rend_cfg, 'body_size_scale',       3.0))

        elem_masses   = np.array([ELEMENTS_LIST[e_ids[i]].mass for i in range(n)], dtype=np.float64)
        actual_masses = w.masses[:n]
        mass_ratios   = actual_masses / elem_masses

        planet_mask = (mass_ratios >= planet_threshold) & (mass_ratios < star_threshold)
        star_mask   = mass_ratios >= star_threshold

        # ---- Composition-derived body colouring (Phase 3) -----------------
        # Instead of a flat planet/star palette swatch, blend each body's
        # colour from the mass-weighted CPK colours of every element it has
        # accreted. A rocky world dominated by Fe + Si is iron-red-grey;
        # a primordial gas giant dominated by H + He is white-pale-cyan;
        # mixed-composition bodies take on intermediate tones. The data
        # comes from world.composition, which is summed through accretion.
        body_mask    = planet_mask | star_mask
        atmosphere_visible = np.zeros(n, dtype=bool)
        atmosphere_colors  = np.zeros((n, 3), dtype=np.float32)
        atmosphere_frac    = np.zeros(n, dtype=np.float32)
        core_visible       = np.zeros(n, dtype=bool)
        core_colors        = np.zeros((n, 3), dtype=np.float32)
        core_frac          = np.zeros(n, dtype=np.float32)
        if body_mask.any():
            elem_cpk = np.array(
                [e.color for e in ELEMENTS_LIST], dtype=np.float32,
            )                                              # (N_elem, 3)
            elem_masses_all = np.array(
                [e.mass for e in ELEMENTS_LIST], dtype=np.float32,
            )                                              # (N_elem,)

            comp     = w.composition[:n][body_mask]        # (B, N_elem)
            totals   = comp.sum(axis=1, keepdims=True)
            totals   = np.maximum(totals, 1e-12)
            fractions = comp / totals                      # (B, N_elem)
            body_colors = (fractions @ elem_cpk).astype(np.float32)  # (B, 3)

            body_sizes = (body_size_scale * np.cbrt(mass_ratios)).astype(np.float32)
            body_sizes[star_mask] *= 1.5

            colors[body_mask, :3] = body_colors
            sizes[body_mask]      = body_sizes[body_mask]

            # ---- Atmosphere extraction (light-element gas envelope) -----
            # Real planets and stars hold gravitationally-bound gas
            # envelopes of their lightest constituents (H, He). We project
            # the body's composition onto its light components and render
            # the result as a separate halo around the body.
            light_elem_mask = elem_masses_all < 5.0           # H, D, He, He3, Li, Be8
            light_comp      = comp[:, light_elem_mask]        # (B, N_light)
            light_total     = light_comp.sum(axis=1)          # (B,)
            light_frac_body = light_total / totals.flatten()  # (B,)
            has_atm         = light_frac_body > 0.05          # ≥5% light material

            if has_atm.any():
                light_cpk      = elem_cpk[light_elem_mask]    # (N_light, 3)
                light_totals_b = np.maximum(light_total[:, None], 1e-12)
                light_fractions = light_comp / light_totals_b
                atm_colors_body = (light_fractions @ light_cpk).astype(np.float32)

                # Scatter back to full-particle arrays
                body_indices = np.where(body_mask)[0]
                atm_indices  = body_indices[has_atm]
                atmosphere_visible[atm_indices] = True
                atmosphere_colors[atm_indices]  = atm_colors_body[has_atm]
                atmosphere_frac[atm_indices]    = light_frac_body[has_atm].astype(np.float32)

            # ---- Dense-core extraction (heavy-element rocky interior) ---
            # Identify "heavy" elements (mass > 20 amu): Fe, Si, Ti, Ca, etc.
            # If a body has ≥ 20% heavy content, render a smaller, denser
            # core marker on top of it — real gravitational differentiation
            # would settle these into the body's centre.
            heavy_elem_mask = elem_masses_all > 20.0
            heavy_comp      = comp[:, heavy_elem_mask]
            heavy_total     = heavy_comp.sum(axis=1)
            heavy_frac_body = heavy_total / totals.flatten()
            has_core        = heavy_frac_body > 0.20

            if has_core.any():
                heavy_cpk        = elem_cpk[heavy_elem_mask]
                heavy_totals_b   = np.maximum(heavy_total[:, None], 1e-12)
                heavy_fractions  = heavy_comp / heavy_totals_b
                core_colors_body = (heavy_fractions @ heavy_cpk).astype(np.float32)

                body_indices = np.where(body_mask)[0]
                core_indices = body_indices[has_core]
                core_visible[core_indices] = True
                core_colors[core_indices]  = core_colors_body[has_core]
                core_frac[core_indices]    = heavy_frac_body[has_core].astype(np.float32)

        self.markers.set_data(
            pos=pos,
            face_color=colors,
            size=sizes,
            edge_width=0,
        )

        # --- Physical temperature per particle (Kelvin) -------------------
        T_K = self._particle_temperatures_K(
            ke, mass_ratios, planet_mask, star_mask, cfg,
        )

        # --- Doppler + Hubble redshift on the rendered colour --------------
        # A source moving radially away from the observer is red-shifted; one
        # approaching is blue-shifted. Wien's displacement says λ·T = const,
        # so the apparent blackbody temperature is
        #     T_apparent = T_emitted / (1 + v_radial / c)
        # In our sim the particle's velocity already contains the Hubble
        # expansion term as well as any kinematic motion, so the same
        # formula handles cosmological redshift too. (At our regime the
        # GR distinction is sub-percent and absorbed into the calibrated
        # ``speed_of_light_sim`` constant.)
        cam_pos = self._camera_position_world()
        c_sim   = float(getattr(cfg, 'speed_of_light_sim', 0.0))
        if c_sim > 0.0 and cam_pos is not None:
            rel        = w.positions[:n] - cam_pos[None, :]
            dist       = np.linalg.norm(rel, axis=1)
            dist_safe  = np.maximum(dist, 1e-6)
            unit_out   = rel / dist_safe[:, None]
            v_radial   = np.einsum('ij,ij->i', w.velocities[:n], unit_out)
            shift_factor = 1.0 + v_radial / c_sim
            # Clip below tiny positive to avoid divide-by-zero / negative T
            T_K_seen   = T_K / np.clip(shift_factor, 1e-3, None)
        else:
            T_K_seen   = T_K

        # Blackbody RGB at the *apparent* (Doppler-shifted) temperature —
        # this is the real colour reaching the observer after the relativistic
        # frame shift.  Inter-particle illumination still uses the rest-frame
        # T_K below (the receiver doesn't sit in the camera's frame).
        bb_rgb = self._blackbody_rgb(T_K_seen)

        # --- Astronomy palette mode replaces the *core* colour with the
        # blackbody colour as well. The halo colour is always blackbody,
        # because radiation has no element identity.
        if self._astronomy_palette:
            astro = np.ones((n, 4), dtype=np.float32)
            astro[:, :3] = bb_rgb
            self.markers.set_data(pos=pos, face_color=astro, size=sizes, edge_width=0)

        # --- Self-emission RGB flux (real Stefan-Boltzmann) ---------------
        # L_self ∝ T⁴. Multiplied by the particle's own Planck blackbody RGB
        # gives an RGB-valued flux. This is what the particle radiates.
        ref_T_K  = float(getattr(cfg, 'reference_temperature_K', 20000.0))
        T_ratio  = (T_K / ref_T_K).astype(np.float32)
        L_self   = T_ratio ** 4                            # (N,) dimensionless luminosity
        rgb_self = bb_rgb * L_self[:, None]                # (N, 3) self-emission flux

        # --- Received illumination (inverse-square from all hot emitters) -
        # Real photon-transport-without-rays: each hot particle is treated
        # as a point source, flux Φ = L/(4πr²) is summed over every other
        # particle. Reflection-nebula colouring (gas takes on the star's
        # colour) emerges naturally because each source contributes its own
        # blackbody RGB to the receiver.
        albedo = float(getattr(cfg, 'illumination_albedo', 0.7))
        rgb_received = compute_received_rgb_flux(
            pos, T_K, self._blackbody_rgb, ref_T_K,
        )

        # --- Total flux = own emission + reflected illumination ----------
        # Bodies (planets / stars) get more emitting surface area → boosted self.
        body_mul = np.ones(n, dtype=np.float32)
        body_mask = planet_mask | star_mask
        if body_mask.any():
            body_mul[body_mask] = _BODY_HALO_MUL
        rgb_total = rgb_self * body_mul[:, None] + albedo * rgb_received

        # --- Beer-Lambert extinction along the line of sight to the camera --
        # Cold gas between an emitter and the observer absorbs / scatters
        # photons:  F_obs = F_emitted · exp(−τ).  Produces dark-nebula
        # silhouettes (Horsehead, Coalsack) when dense cold clouds sit in
        # front of bright sources.  Uses the cam_pos computed earlier for
        # the Doppler shift.
        #
        # Wavelength-dependent ("reddening"): blue is absorbed more than
        # red in real interstellar dust (Fitzpatrick 1999 extinction curve).
        # We apply per-RGB-channel optical depths via the extinction_rgb_weights
        # vector — same exp(−τ) but with different τ per colour.
        if cam_pos is not None:
            opacity_base = float(getattr(cfg, 'extinction_opacity', 0.05))
            kernel_sigma = float(getattr(cfg, 'extinction_kernel_sigma', 30.0))
            optical_d = per_particle_optical_depth(
                T_K, sizes,
                cold_threshold_K=3000.0,
                base_opacity=opacity_base,
            )
            visibility = compute_visibility(pos, optical_d, cam_pos, kernel_sigma)

            rgb_weights = getattr(cfg, 'extinction_rgb_weights', None)
            if rgb_weights is not None:
                # exp(−τ_channel · weight_channel) per RGB.  Equivalent to:
                #     visibility_channel = visibility ** weight_channel
                w_arr = np.array(rgb_weights, dtype=np.float32)
                if w_arr.shape == (3,):
                    vis_r = visibility ** w_arr[0]
                    vis_g = visibility ** w_arr[1]
                    vis_b = visibility ** w_arr[2]
                    visibility_rgb = np.stack([vis_r, vis_g, vis_b], axis=1)
                    rgb_total = rgb_total * visibility_rgb
                else:
                    rgb_total = rgb_total * visibility[:, None]
            else:
                rgb_total = rgb_total * visibility[:, None]

        # Raw scalar luminance (BT.601 perceptual weights) — drives halo size
        # before tone-mapping, because the *physical* radiated power sets the
        # apparent angular extent (Stefan-Boltzmann), not the camera response.
        luminance_raw = (0.299 * rgb_total[:, 0]
                         + 0.587 * rgb_total[:, 1]
                         + 0.114 * rgb_total[:, 2]).astype(np.float32)
        sqrt_lum_raw = np.sqrt(np.maximum(luminance_raw, 0.0))

        # --- HDR tone mapping (Reinhard) ---------------------------------
        # L_displayed = exposure · L_raw / (1 + exposure · L_raw)
        # Real-camera response: linear at low flux, asymptotic to 1 at high.
        # Bright fusing cores don't clip to white; dim gas isn't crushed black.
        # The single "exposure" knob is the same one a photographer dials.
        exposure = float(getattr(cfg, 'render_exposure', 5.0))
        e_rgb = exposure * rgb_total
        rgb_displayed = (e_rgb / (1.0 + e_rgb)).astype(np.float32)

        # Halo size: from raw (physical) luminance, not tone-mapped value
        inner_size = (sizes * (_HALO_BASE_SIZE_MUL
                               + sqrt_lum_raw * _HALO_RADIATION_GAIN)).astype(np.float32)
        outer_size = (inner_size * 2.2).astype(np.float32)

        # Halo colour: tone-mapped RGB (carries colour AND brightness)
        # Halo alpha: tone-mapped luminance with a soft floor
        disp_luminance = (0.299 * rgb_displayed[:, 0]
                          + 0.587 * rgb_displayed[:, 1]
                          + 0.114 * rgb_displayed[:, 2]).astype(np.float32)
        inner_alpha = np.clip(_HALO_INNER_ALPHA + 0.7 * disp_luminance, 0.0, 1.0)
        outer_alpha = np.clip(_HALO_OUTER_ALPHA + 0.3 * disp_luminance, 0.0, 0.7)

        inner_rgba = np.ones((n, 4), dtype=np.float32)
        inner_rgba[:, :3] = rgb_displayed
        inner_rgba[:,  3] = inner_alpha

        outer_rgba = np.ones((n, 4), dtype=np.float32)
        outer_rgba[:, :3] = rgb_displayed
        outer_rgba[:,  3] = outer_alpha

        self.markers_inner_glow.set_data(
            pos=pos, face_color=inner_rgba, size=inner_size, edge_width=0,
        )
        self.markers_nebula.set_data(
            pos=pos, face_color=outer_rgba, size=outer_size, edge_width=0,
        )

        # --- Dense core layer (Phase 3 differentiation) -------------------
        if core_visible.any():
            # Core is rendered at ~50% of body size, opaque, in the heavy
            # constituents' colour — visible "centre" of the body.
            core_size = (sizes * 0.55).astype(np.float32)
            core_rgba = np.zeros((n, 4), dtype=np.float32)
            core_rgba[:, :3] = core_colors
            core_rgba[:,  3] = np.where(core_visible, 1.0, 0.0).astype(np.float32)
            self.markers_core.set_data(
                pos=pos, face_color=core_rgba, size=core_size, edge_width=0,
            )
        else:
            self.markers_core.set_data(
                pos=np.zeros((1, 3), dtype=np.float32),
                face_color=(0.0, 0.0, 0.0, 0.0),
                size=1, edge_width=0,
            )

        # --- Atmosphere halo (Phase 3) -----------------------------------
        # Gas envelope around bodies that have absorbed substantial light
        # elements.  Rendered between the body and the thermal halo, in
        # the mass-weighted CPK colour of the body's light constituents.
        if atmosphere_visible.any():
            atm_alpha = np.clip(0.5 * atmosphere_frac, 0.0, 0.6).astype(np.float32)
            atm_alpha[~atmosphere_visible] = 0.0
            # Size: 1.5× body size + scale with light fraction
            atm_size  = (sizes * (1.5 + 1.8 * atmosphere_frac)).astype(np.float32)
            atm_rgba  = np.ones((n, 4), dtype=np.float32)
            atm_rgba[:, :3] = atmosphere_colors
            atm_rgba[:,  3] = atm_alpha
            self.markers_atmosphere.set_data(
                pos=pos, face_color=atm_rgba, size=atm_size, edge_width=0,
            )
        else:
            self.markers_atmosphere.set_data(
                pos=np.zeros((1, 3), dtype=np.float32),
                face_color=(0.0, 0.0, 0.0, 0.0),
                size=1, edge_width=0,
            )

        # --- Selection ring + info text ------------------------------------
        sel = self._selected
        if sel is not None:
            self._select_marker.set_data(
                pos=pos[sel:sel + 1],
                face_color=(0.0, 0.0, 0.0, 0.0),
                size=float(sizes[sel]) + 10.0,
                edge_width=2,
                edge_color=(1.0, 1.0, 1.0, 1.0),
            )
            self._info_text.text = self._particle_info(sel)
        else:
            self._select_marker.set_data(
                pos=np.zeros((1, 3), dtype=np.float32),
                face_color=(0.0, 0.0, 0.0, 0.0),
                size=1,
                edge_width=0,
            )
            self._info_text.text = ''

        # --- Bond lines ----------------------------------------------------
        if cfg.show_bonds and len(w.bonds) > 0:
            bond_pts = []
            for bond in w.bonds:
                bond_pts.append(w.positions[bond.i])
                bond_pts.append(w.positions[bond.j])
            self.lines.set_data(
                pos=np.array(bond_pts, dtype=np.float32),
                color=(0.9, 0.9, 0.9, cfg.bond_alpha),
                connect='segments',
            )
        else:
            self.lines.set_data(pos=np.zeros((2, 3)), color=(0, 0, 0, 0))

        self.canvas.update()

    def _update_title(self) -> None:
        w      = self.world
        counts = Counter(ELEMENTS_LIST[w.elem_ids[i]].symbol for i in range(w.n))
        top    = ' '.join(f'{s}:{c}' for s, c in counts.most_common(4))
        sc     = state_counts(w)
        status = 'PAUSED ' if self.paused else ''

        # Conservation diagnostics (throttled — gravitational PE is O(N²))
        self._energy_frame += 1
        if self._energy_frame % _ENERGY_REFRESH_FRAMES == 0 or self._energy_ref is None:
            self._ke = kinetic_energy(w)
            self._pe = gravitational_pe(w)
            self._momentum_mag = float(np.linalg.norm(linear_momentum(w)))
            if self._energy_ref is None and w.n >= 2:
                self._energy_ref = self._ke + self._pe   # set baseline once

        e_total = self._ke + self._pe
        drift   = ''
        if self._energy_ref is not None and self._energy_ref != 0:
            drift_pct = 100.0 * (e_total - self._energy_ref) / abs(self._energy_ref)
            drift     = f' drift={drift_pct:+.2f}%'

        # Molecule census — show top-N by count
        mol_str = ''
        poly_str = ''
        if w.bonds:
            mc = molecule_counts(w)
            if mc:
                top_mols = sorted(mc.items(), key=lambda kv: -kv[1])[:4]
                mol_str  = '  mols:[' + ' '.join(f'{f}:{c}' for f, c in top_mols) + ']'

            # Phase-5 polymer census: longest chain, mean chain, motif counts.
            # Only shown when at least one polymer chain reaches length ≥ 3,
            # so simple diatomics don't clutter the HUD.
            stats = polymer_stats(w)
            if stats['max_chain'] >= 3:
                poly_str = (
                    f'  poly:[max={stats["max_chain"]} '
                    f'mean={stats["mean_chain"]:.1f}]'
                )
                if stats['motif_counts']:
                    motif_bits = ' '.join(
                        f'{k}:{v}' for k, v in stats['motif_counts'].items()
                    )
                    poly_str += f' motifs:[{motif_bits}]'

        self.canvas.title = (
            f'Universe | {status}'
            f'n={w.n}  bonds={len(w.bonds)}  '
            f'fusions={w.total_fusions}  accreted={w.total_accretions}  '
            f'rxns={w.total_reactions}  '
            f'S:{sc["solid"]} L:{sc["liquid"]} G:{sc["gas"]} P:{sc["plasma"]}  '
            f'E={e_total:.2e}{drift}  clamps={w.total_velocity_clamps}  '
            f't={w.time:.1f}s  fps={self._fps:.0f}  '
            f'speed={self.speed:.1f}x  [{top}]'
            f'{mol_str}{poly_str}'
        )

    # ------------------------------------------------------------------
    # Keyboard
    # ------------------------------------------------------------------

    def _on_key(self, event) -> None:
        k   = event.key.name if event.key else ''
        cam = self.view.camera
        if k == 'Space':
            self.paused = not self.paused
        elif k in ('Equal', '+'):
            self.speed = min(self.speed * _SPEED_UP, _SPEED_MAX)
        elif k in ('Minus', '-'):
            self.speed = max(self.speed / _SPEED_UP, _SPEED_MIN)
        elif k == 'Left':
            cam.azimuth -= _CAM_ORBIT_STEP
        elif k == 'Right':
            cam.azimuth += _CAM_ORBIT_STEP
        elif k == 'Up':
            cam.elevation = min(cam.elevation + _CAM_ORBIT_STEP, 90.0)
        elif k == 'Down':
            cam.elevation = max(cam.elevation - _CAM_ORBIT_STEP, -90.0)
        elif k == 'PageUp':
            cam.distance = max(cam.distance / _CAM_ZOOM_STEP, 1.0)
        elif k == 'PageDown':
            cam.distance *= _CAM_ZOOM_STEP
        elif k == 'R':
            cam.azimuth   = _CAM_AZIMUTH_0
            cam.elevation = _CAM_ELEVATION_0
            cam.distance  = self.cfg.renderer.camera_distance
            cam.center    = (0.0, 0.0, 0.0)
            self._cinematic_mode = False
        elif k == 'F':
            self._cinematic_mode = not self._cinematic_mode
            self._cinematic_t = 0.0
        elif k == 'P':
            self._astronomy_palette = not self._astronomy_palette
        elif k in ('Q', 'Escape'):
            self.canvas.app.quit()

    # ------------------------------------------------------------------
    # Mouse — particle picking
    # ------------------------------------------------------------------

    def _on_mouse_press(self, event) -> None:
        if event.button != 1:
            return
        w = self.world
        if w.n == 0:
            return

        click_xy = np.array(event.pos[:2], dtype=np.float32)
        screen_xy, in_front = self._project_to_canvas(w.positions[:w.n])
        if not in_front.any():
            return

        dists = np.linalg.norm(screen_xy - click_xy, axis=1)
        dists[~in_front] = np.inf
        nearest = int(np.argmin(dists))

        if dists[nearest] < _PICK_PIXEL_RADIUS:
            self._selected = nearest if self._selected != nearest else None
        else:
            self._selected = None

    def _project_to_canvas(self, pos3d: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Project 3D world positions to 2D canvas pixel coords.

        Returns (screen_xy, in_front) where in_front masks particles
        not behind the camera. Handles the perspective divide that
        vispy's transform.map() does not do automatically.
        """
        pts = np.asarray(pos3d, dtype=np.float32)
        tr  = self.markers.transforms.get_transform('visual', 'canvas')
        out = tr.map(pts)

        if out.shape[1] == 4:
            w_h      = out[:, 3]
            in_front = w_h > 1e-6
            safe_w   = np.where(np.abs(w_h) > 1e-9, w_h, 1.0)
            screen_xy = out[:, :2] / safe_w[:, np.newaxis]
        else:
            screen_xy = out[:, :2]
            in_front  = np.ones(pts.shape[0], dtype=bool)

        return screen_xy, in_front

    # ------------------------------------------------------------------
    # Real-physics colour: Planck blackbody RGB
    # ------------------------------------------------------------------
    # Maps temperature in Kelvin to perceived RGB using Tanner Helland's
    # polynomial fit to the CIE colour-matching functions over the Planck
    # locus (the standard graphics conversion used in physically-based
    # rendering pipelines).
    #
    # Reference:
    #   Tanner Helland, "How to Convert Temperature (K) to RGB" (2012)
    #   The piecewise fit reproduces the CIE blackbody curve from ~1000 K
    #   (deep red) through ~5800 K (sunlight) to ~30000 K (hot O-class blue).

    def _camera_position_world(self) -> np.ndarray | None:
        """Camera position in world coordinates, derived from the
        TurntableCamera's azimuth / elevation / distance / center.

        vispy's TurntableCamera doesn't expose `position` directly; we
        reconstruct it from the spherical-coordinate parameters. Used by
        the Beer-Lambert extinction pass to set up line-of-sight rays.
        """
        cam = self.view.camera
        try:
            azimuth   = float(cam.azimuth)
            elevation = float(cam.elevation)
            distance  = float(cam.distance)
            center    = np.asarray(cam.center, dtype=np.float32)
        except (AttributeError, TypeError):
            return None
        # Spherical → Cartesian (vispy uses degrees + a specific convention)
        a = np.deg2rad(azimuth)
        e = np.deg2rad(elevation)
        offset = distance * np.array([
            np.cos(e) * np.cos(a),
            np.cos(e) * np.sin(a),
            np.sin(e),
        ], dtype=np.float32)
        return (center + offset).astype(np.float32)

    @staticmethod
    def _blackbody_rgb(T_kelvin: np.ndarray) -> np.ndarray:
        """Vectorised Tanner-Helland Planck-curve RGB. Output: (N, 3) in 0..1."""
        T = np.clip(T_kelvin, 1000.0, 40000.0) / 100.0

        # Red
        r = np.where(
            T <= 66.0,
            255.0,
            329.698727446 * np.power(np.maximum(T - 60.0, 1e-9), -0.1332047592),
        )
        # Green
        g_low  = 99.4708025861 * np.log(np.maximum(T, 1e-9)) - 161.1195681661
        g_high = 288.1221695283 * np.power(np.maximum(T - 60.0, 1e-9), -0.0755148492)
        g = np.where(T <= 66.0, g_low, g_high)
        # Blue
        b = np.where(
            T >= 66.0,
            255.0,
            np.where(
                T <= 19.0,
                0.0,
                138.5177312231 * np.log(np.maximum(T - 10.0, 1e-9)) - 305.0447927307,
            ),
        )

        rgb = np.clip(np.column_stack([r, g, b]) / 255.0, 0.0, 1.0)
        return rgb.astype(np.float32)

    def _particle_temperatures_K(
        self,
        ke: np.ndarray,
        mass_ratios: np.ndarray,
        planet_mask: np.ndarray,
        star_mask: np.ndarray,
        cfg,
    ) -> np.ndarray:
        """Effective Kelvin temperature per particle.

        Two regimes — each defensible against the standard textbook source:

        1. **Free atoms / gas**: kinetic temperature.
            T_K = (KE / hot_temperature_threshold) · reference_temperature_K
           This is a *linear* scaling because our sim's energy unit is
           itself a custom calibration; the reference point (KE = hot_T →
           T = ref_T_K) is the only calibrated constant.

        2. **Accreted bodies (planets / stars)**: stellar
           mass-luminosity relation. Main-sequence stars satisfy
            T ∝ M^(1/2)   (Eddington luminosity argument)
           giving small bodies (rocky planets) cool red colours and
           massive bodies (stellar masses) hot blue-white colours.
           Calibrated so a body at the star_mass_threshold sits at the
           reference temperature.

        For an accreted body we take the *maximum* of the two estimates
        so a fast-moving star is not mis-coloured by its kinetic
        proxy alone.
        """
        hot_T_sim = float(getattr(cfg, 'hot_temperature_threshold', 5.0e5))
        ref_T_K   = float(getattr(cfg, 'reference_temperature_K', 20000.0))

        # Gas / atomic regime
        T_K = (ke / hot_T_sim) * ref_T_K

        # Accreted bodies: stellar mass-luminosity relation T ∝ √M
        star_thresh = float(getattr(cfg, 'star_mass_threshold', 200.0))
        body_mask = planet_mask | star_mask
        if body_mask.any():
            T_K_body = ref_T_K * np.sqrt(np.maximum(
                mass_ratios[body_mask] / star_thresh, 0.0,
            ))
            T_K[body_mask] = np.maximum(T_K[body_mask], T_K_body)

        return T_K.astype(np.float32)

    def _particle_info(self, i: int) -> str:
        w    = self.world
        elem = ELEMENTS_LIST[w.elem_ids[i]]
        mass_ratio = w.masses[i] / elem.mass
        speed      = float(np.linalg.norm(w.velocities[i]))
        pos        = w.positions[i]

        bonded = [
            ELEMENTS_LIST[w.elem_ids[b.j if b.i == i else b.i]].symbol
            for b in w.bonds if b.i == i or b.j == i
        ]

        planet_threshold = float(getattr(self.cfg.renderer, 'planet_mass_threshold', 1e9))
        star_threshold   = float(getattr(self.cfg.renderer, 'star_mass_threshold',   1e9))
        if mass_ratio >= star_threshold:
            kind = 'STAR'
        elif mass_ratio >= planet_threshold:
            kind = 'PLANET'
        else:
            kind = 'atom'

        lines = [
            f'[{kind}] {elem.name} ({elem.symbol})  Z={elem.Z}',
            f'mass: {w.masses[i]:.3g}  ratio: {mass_ratio:.2f}x',
            f'pos:   ({pos[0]:.1f}, {pos[1]:.1f}, {pos[2]:.1f})',
            f'speed: {speed:.1f}',
            f'bonds: {w.bond_counts[i]}' + (f'  → {" ".join(bonded)}' if bonded else ''),
        ]

        # Composition breakdown for accreted bodies (Phase 3) — shows the
        # mass fractions of every absorbed element. A pristine atom has
        # 100% of its own element; an accreted body has a mix that drove
        # its rendered colour.
        comp = w.composition[i]
        total = comp.sum()
        if total > 0 and mass_ratio >= 1.5:                 # only worth showing for accreted bodies
            mass_fractions = comp / total
            top = sorted(
                ((ELEMENTS_LIST[k].symbol, float(mass_fractions[k]))
                 for k in range(len(ELEMENTS_LIST)) if mass_fractions[k] > 0.01),
                key=lambda kv: -kv[1],
            )[:5]
            if top:
                comp_str = ' '.join(f'{s}:{f*100:.0f}%' for s, f in top)
                lines.append(f'comp:  {comp_str}')

        # Molecule membership — if this atom is part of a bonded cluster,
        # show the formula and (where known) the common name.
        if w.bond_counts[i] > 0:
            from sim.molecules import molecule_of
            mol = molecule_of(w, i)
            if mol is not None:
                if mol.name:
                    lines.append(f'molecule: {mol.formula}  ({mol.name})')
                else:
                    lines.append(f'molecule: {mol.formula}')

            # Phase-5 polymer topology — chain length, branches, rings, motifs.
            # Only shown for components large enough to have non-trivial
            # structure (size ≥ 3); simple diatomics are already covered by
            # the molecule line above.
            poly = polymer_of(w, i)
            if poly is not None and poly.size >= 3:
                topo = f'chain={poly.chain_length}'
                if poly.branches:
                    topo += f' branches={poly.branches}'
                if poly.rings:
                    topo += f' rings={poly.rings}'
                lines.append(f'polymer: size={poly.size}  {topo}')
                if poly.motifs:
                    lines.append(f'motifs: {", ".join(poly.motifs)}')

        return '\n'.join(lines)

    # ------------------------------------------------------------------
    # Run
    # ------------------------------------------------------------------

    def run(self) -> None:
        app.run()
