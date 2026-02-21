"""
Distributions Module — Analiza rozkładów i topologii struktury
================================================================
Moduł obliczeniowy dla:
  - Rozkład grubości ścian (wall thickness distribution)
  - Rozkład szerokości kanałów (channel width distribution)
  - Connectivity / dead-end fraction / perkolacja
  - Accessible Surface Area (ASA, probe-based)
  - Throat / constriction analysis
  - Printability score

Metody bazują na voxelizacji + distance transform (EDT).

Refs:
  - Hildebrand, T. & Rüegsegger, P. (1997) JBMR 12(7):1167-1174
    (local thickness via distance transform)
  - Münch, B. & Holzer, L. (2008) J. Am. Ceram. Soc. 91(12):4059-4067
    (medial axis / throat detection)

Author: Claude (Anthropic)
Date: 2026-02-20
"""

import numpy as np
from scipy.ndimage import (
    distance_transform_edt, label, binary_erosion,
    generate_binary_structure, binary_dilation
)
from scipy.ndimage import maximum_filter
from typing import Dict, Optional, Tuple, List
import trimesh
import time
import logging

logger = logging.getLogger(__name__)


class VoxelAnalyzer:
    """
    Analiza struktury gyroidalnej na podstawie voxelizacji.

    Pipeline:
    1. Mesh → voxel grid (solid/void)
    2. Distance transform na solid i void
    3. Ekstrakcja metryk z distance transform
    """

    def __init__(self, mesh: trimesh.Trimesh,
                 container_geometry: str,
                 container_params: Dict,
                 unit_cell_mm: float,
                 printer_pixel_um: float = 22.0):
        """
        Args:
            mesh: Trimesh object
            container_geometry: 'cylinder', 'sphere', etc.
            container_params: {'diameter': float, 'height': float, ...}
            unit_cell_mm: Unit cell size [mm] — determines voxel resolution
            printer_pixel_um: Printer pixel size [µm] for printability checks
        """
        self.mesh = mesh
        self.container_geometry = container_geometry
        self.container_params = container_params
        self.unit_cell_mm = unit_cell_mm
        self.printer_pixel_um = printer_pixel_um

        # Cache
        self._voxel_grid = None
        self._pitch = None
        self._solid_edt = None  # Distance transform of solid phase
        self._void_edt = None   # Distance transform of void phase

    def _voxelize(self, voxels_per_unit_cell: int = 12):
        """
        Voxelize mesh with given resolution.

        Args:
            voxels_per_unit_cell: Voxels per unit cell (higher = better, slower)
        """
        if self._voxel_grid is not None:
            return

        pitch = self.unit_cell_mm / voxels_per_unit_cell
        logger.info(f"Voxelizing: pitch={pitch:.4f} mm, "
                    f"resolution={voxels_per_unit_cell} vox/unit_cell")

        t0 = time.time()
        voxel = self.mesh.voxelized(pitch=pitch)
        self._voxel_grid = voxel.matrix  # True = solid
        self._pitch = pitch

        logger.info(f"Voxel grid: {self._voxel_grid.shape}, "
                    f"time={time.time()-t0:.1f}s")

    def _ensure_edt(self, voxels_per_unit_cell: int = 12):
        """Compute distance transforms if not cached."""
        self._voxelize(voxels_per_unit_cell)

        if self._solid_edt is None:
            logger.info("Computing EDT for solid phase...")
            # EDT on solid: distance from each solid voxel to nearest void
            self._solid_edt = distance_transform_edt(
                self._voxel_grid, sampling=self._pitch
            )

        if self._void_edt is None:
            logger.info("Computing EDT for void phase...")
            # EDT on void: distance from each void voxel to nearest solid
            self._void_edt = distance_transform_edt(
                ~self._voxel_grid, sampling=self._pitch
            )

    # ==================== WALL THICKNESS DISTRIBUTION ====================

    def calc_wall_thickness_distribution(self, voxels_per_uc: int = 12) -> Dict:
        """
        Rozkład grubości ścian (local wall thickness).

        Metoda Hildebrand & Rüegsegger (1997):
        - EDT na fazie solid
        - Local maxima EDT = promień wpisanej sfery = połowa grubości lokalnej
        - Wall thickness = 2 × local_max(EDT_solid)

        Returns:
            Dict z histogramem, percentylami, statystykami
        """
        t0 = time.time()
        self._ensure_edt(voxels_per_uc)

        # Local maxima w solid EDT (szczyty = centra ścian)
        solid_mask = self._voxel_grid
        edt = self._solid_edt.copy()
        edt[~solid_mask] = 0

        # Filtr max 3x3x3 — local maxima to gdzie edt == local_max
        local_max = maximum_filter(edt, size=3)
        is_local_max = (edt == local_max) & solid_mask & (edt > 0)

        # Wall thickness = 2 × inscribed sphere radius
        wall_values = edt[is_local_max] * 2.0  # mm
        wall_values_um = wall_values * 1000  # µm

        if len(wall_values) == 0:
            return {'error': 'No wall thickness measurements found'}

        # Histogram
        n_bins = min(50, max(10, len(wall_values) // 20))
        hist_counts, hist_edges = np.histogram(wall_values_um, bins=n_bins)

        # Percentyle
        p10 = np.percentile(wall_values_um, 10)
        p50 = np.percentile(wall_values_um, 50)
        p90 = np.percentile(wall_values_um, 90)

        results = {
            'mean_um': float(np.mean(wall_values_um)),
            'std_um': float(np.std(wall_values_um)),
            'min_um': float(np.min(wall_values_um)),
            'max_um': float(np.max(wall_values_um)),
            'p10_um': float(p10),
            'p50_um': float(p50),
            'p90_um': float(p90),
            'n_measurements': len(wall_values),
            'hist_counts': hist_counts.tolist(),
            'hist_edges_um': hist_edges.tolist(),
            'uniformity': float(1.0 - np.std(wall_values_um) / np.mean(wall_values_um))
                          if np.mean(wall_values_um) > 0 else 0,
            'calculation_time': time.time() - t0,
        }

        return results

    # ==================== CHANNEL WIDTH DISTRIBUTION ====================

    def calc_channel_width_distribution(self, voxels_per_uc: int = 12) -> Dict:
        """
        Rozkład szerokości kanałów (local channel diameter).

        Analogiczna metoda do wall thickness, ale na void phase.
        EDT na void → local maxima → 2 × radius = channel width.

        Returns:
            Dict z histogramem, percentylami, statystykami
        """
        t0 = time.time()
        self._ensure_edt(voxels_per_uc)

        void_mask = ~self._voxel_grid
        edt = self._void_edt.copy()
        edt[~void_mask] = 0

        # Local maxima
        local_max = maximum_filter(edt, size=3)
        is_local_max = (edt == local_max) & void_mask & (edt > 0)

        channel_values = edt[is_local_max] * 2.0  # mm
        channel_values_um = channel_values * 1000  # µm

        if len(channel_values) == 0:
            return {'error': 'No channel width measurements found'}

        n_bins = min(50, max(10, len(channel_values) // 20))
        hist_counts, hist_edges = np.histogram(channel_values_um, bins=n_bins)

        p10 = np.percentile(channel_values_um, 10)
        p50 = np.percentile(channel_values_um, 50)
        p90 = np.percentile(channel_values_um, 90)

        results = {
            'mean_um': float(np.mean(channel_values_um)),
            'std_um': float(np.std(channel_values_um)),
            'min_um': float(np.min(channel_values_um)),
            'max_um': float(np.max(channel_values_um)),
            'p10_um': float(p10),
            'p50_um': float(p50),
            'p90_um': float(p90),
            'n_measurements': len(channel_values),
            'hist_counts': hist_counts.tolist(),
            'hist_edges_um': hist_edges.tolist(),
            'uniformity': float(1.0 - np.std(channel_values_um) / np.mean(channel_values_um))
                          if np.mean(channel_values_um) > 0 else 0,
            'calculation_time': time.time() - t0,
        }

        return results

    # ==================== CONNECTIVITY / DEAD-END ====================

    def calc_connectivity(self, voxels_per_uc: int = 12) -> Dict:
        """
        Analiza connectivity kanałów (void space).

        Metoda:
        1. Label connected components w void space
        2. Sprawdź które komponenty łączą inlet z outlet
           (Z_min → Z_max dla cylindra)
        3. Dead-end = void niepołączony z through-path

        Returns:
            Dict z connected porosity, dead-end fraction, etc.
        """
        t0 = time.time()
        self._voxelize(voxels_per_uc)

        void_mask = ~self._voxel_grid

        # Label connected components (26-connectivity for 3D)
        struct = generate_binary_structure(3, 3)  # 26-connectivity
        labeled, n_components = label(void_mask, structure=struct)

        total_void_voxels = np.sum(void_mask)

        if total_void_voxels == 0:
            return {'error': 'No void voxels found'}

        # Znajdź komponenty łączące inlet z outlet
        # Inlet = bottom slice (z=0), Outlet = top slice (z=max)
        nz = self._voxel_grid.shape[2]

        inlet_labels = set(np.unique(labeled[:, :, 0])) - {0}
        outlet_labels = set(np.unique(labeled[:, :, nz-1])) - {0}

        # Through-connected = present at both inlet AND outlet
        through_labels = inlet_labels & outlet_labels

        # Count voxels in each category
        through_voxels = 0
        inlet_only_voxels = 0
        outlet_only_voxels = 0
        isolated_voxels = 0

        for comp_label in range(1, n_components + 1):
            comp_size = np.sum(labeled == comp_label)
            if comp_label in through_labels:
                through_voxels += comp_size
            elif comp_label in inlet_labels:
                inlet_only_voxels += comp_size
            elif comp_label in outlet_labels:
                outlet_only_voxels += comp_size
            else:
                isolated_voxels += comp_size

        dead_end_voxels = inlet_only_voxels + outlet_only_voxels + isolated_voxels

        # Frakcje
        connected_fraction = through_voxels / total_void_voxels if total_void_voxels > 0 else 0
        dead_end_fraction = dead_end_voxels / total_void_voxels if total_void_voxels > 0 else 0
        isolated_fraction = isolated_voxels / total_void_voxels if total_void_voxels > 0 else 0

        # Connected porosity (only through-connected void)
        total_voxels = self._voxel_grid.size
        connected_porosity = through_voxels / total_voxels

        results = {
            'n_components': n_components,
            'n_through_connected': len(through_labels),
            'connected_fraction': float(connected_fraction),
            'connected_porosity': float(connected_porosity),
            'dead_end_fraction': float(dead_end_fraction),
            'isolated_fraction': float(isolated_fraction),
            'through_voxels': int(through_voxels),
            'dead_end_voxels': int(dead_end_voxels),
            'isolated_voxels': int(isolated_voxels),
            'total_void_voxels': int(total_void_voxels),
            'calculation_time': time.time() - t0,
        }

        # Warnings
        warnings = []
        if connected_fraction < 0.9:
            warnings.append(f'Only {connected_fraction:.0%} of void is through-connected!')
        if isolated_fraction > 0.05:
            warnings.append(f'{isolated_fraction:.0%} of void is completely isolated (trapped).')
        if dead_end_fraction > 0.2:
            warnings.append(f'{dead_end_fraction:.0%} of void is dead-end — may slow desorption.')
        results['warnings'] = warnings

        return results

    # ==================== ACCESSIBLE SURFACE AREA ====================

    def calc_accessible_surface_area(self, probe_radii_nm: List[float] = None,
                                     voxels_per_uc: int = 12) -> Dict:
        """
        Accessible Surface Area (ASA) — powierzchnia dostępna dla sondy.

        Metoda:
        1. Erozja void space o promień sondy (w voxelach)
        2. Powierzchnia po erozji = ASA
        3. ASA/total_SA = accessibility fraction

        Fizycznie: nie cała powierzchnia jest dostępna dla cząsteczki
        o skończonym rozmiarze. Wąskie kanały i szyjki są niedostępne
        dla dużych analitów.

        Args:
            probe_radii_nm: Lista promieni sondy [nm]. Default: [0, 1, 2, 5, 10]
            voxels_per_uc: Rozdzielczość voxelizacji

        Returns:
            Dict z ASA dla każdego promienia sondy
        """
        t0 = time.time()
        self._voxelize(voxels_per_uc)

        if probe_radii_nm is None:
            probe_radii_nm = [0, 1, 2, 5, 10]

        pitch_nm = self._pitch * 1e6  # mm → nm
        total_sa = self.mesh.area  # mm²

        profiles = []

        for probe_r_nm in probe_radii_nm:
            if probe_r_nm == 0:
                # No erosion — full surface
                profiles.append({
                    'probe_radius_nm': 0,
                    'accessible_SA_mm2': total_sa,
                    'fraction': 1.0,
                })
                continue

            # Convert probe radius to voxels
            probe_voxels = probe_r_nm / pitch_nm

            if probe_voxels < 0.5:
                # Probe smaller than voxel — no effect
                profiles.append({
                    'probe_radius_nm': probe_r_nm,
                    'accessible_SA_mm2': total_sa,
                    'fraction': 1.0,
                })
                continue

            # Erosion iterations (approximate)
            n_erosions = max(1, int(round(probe_voxels)))

            # Erode void space
            void_mask = ~self._voxel_grid
            eroded_void = void_mask.copy()
            for _ in range(n_erosions):
                eroded_void = binary_erosion(eroded_void)

            # Accessible surface = boundary between solid and eroded void
            # Count boundary voxels (solid that touch eroded void)
            dilated_void = binary_dilation(eroded_void)
            boundary = dilated_void & self._voxel_grid
            n_boundary = np.sum(boundary)

            # Original boundary (solid touching original void)
            dilated_orig = binary_dilation(void_mask)
            boundary_orig = dilated_orig & self._voxel_grid
            n_boundary_orig = max(1, np.sum(boundary_orig))

            # ASA fraction
            fraction = n_boundary / n_boundary_orig
            accessible_sa = total_sa * fraction

            profiles.append({
                'probe_radius_nm': probe_r_nm,
                'accessible_SA_mm2': float(accessible_sa),
                'fraction': float(fraction),
            })

        results = {
            'profiles': profiles,
            'total_SA_mm2': total_sa,
            'pitch_nm': pitch_nm,
            'calculation_time': time.time() - t0,
        }

        return results

    # ==================== THROAT ANALYSIS ====================

    def calc_throat_distribution(self, voxels_per_uc: int = 12) -> Dict:
        """
        Analiza przewężeń (throat / constriction analysis).

        Metoda:
        - EDT na void space
        - Throat = local minimum EDT wzdłuż medial axis
        - Szukamy saddle points (punkty siodłowe) w EDT

        Uproszczenie: szukamy local minima wzdłuż osi Z
        w obrębie medial axis kanałów.

        Returns:
            Dict z rozkładem throat sizes
        """
        t0 = time.time()
        self._ensure_edt(voxels_per_uc)

        void_mask = ~self._voxel_grid
        edt = self._void_edt.copy()
        edt[~void_mask] = 0

        # Medial axis (approximate): local maxima w EDT void
        local_max = maximum_filter(edt, size=3)
        medial_mask = (edt == local_max) & void_mask & (edt > 0)

        if np.sum(medial_mask) == 0:
            return {'error': 'No medial axis found'}

        # Na medial axis szukamy local minima (=throats)
        # Filtr minimum 3x3x3
        from scipy.ndimage import minimum_filter
        local_min = minimum_filter(edt, size=5)
        is_throat = (edt == local_min) & medial_mask

        throat_radii = edt[is_throat]  # mm (promienie)
        throat_diameters_um = throat_radii * 2.0 * 1000  # µm

        if len(throat_diameters_um) == 0:
            # Fallback: użyj P10 z channel distribution jako throat estimate
            return {
                'error': 'No clear throats detected (structure may be very uniform)',
                'note': 'For gyroid TPMS, throats are typically ~0.7× mean channel width'
            }

        n_bins = min(40, max(10, len(throat_diameters_um) // 10))
        hist_counts, hist_edges = np.histogram(throat_diameters_um, bins=n_bins)

        p10 = np.percentile(throat_diameters_um, 10)
        p50 = np.percentile(throat_diameters_um, 50)
        p90 = np.percentile(throat_diameters_um, 90)

        # Throat-to-pore ratio
        pore_values = edt[medial_mask] * 2.0 * 1000
        mean_pore = np.mean(pore_values) if len(pore_values) > 0 else 1
        throat_to_pore = float(np.mean(throat_diameters_um) / mean_pore) if mean_pore > 0 else 0

        results = {
            'mean_um': float(np.mean(throat_diameters_um)),
            'std_um': float(np.std(throat_diameters_um)),
            'min_um': float(np.min(throat_diameters_um)),
            'max_um': float(np.max(throat_diameters_um)),
            'p10_um': float(p10),
            'p50_um': float(p50),
            'p90_um': float(p90),
            'throat_to_pore_ratio': throat_to_pore,
            'n_throats': len(throat_diameters_um),
            'hist_counts': hist_counts.tolist(),
            'hist_edges_um': hist_edges.tolist(),
            'calculation_time': time.time() - t0,
        }

        # Warnings
        warnings = []
        min_throat = np.min(throat_diameters_um)
        if min_throat < 30:
            warnings.append(f'Throat min {min_throat:.0f} µm — risk of clogging!')
        if throat_to_pore < 0.3:
            warnings.append(f'Throat/pore ratio {throat_to_pore:.2f} — high flow resistance.')
        results['warnings'] = warnings

        return results

    # ==================== PRINTABILITY SCORE ====================

    def calc_printability(self, voxels_per_uc: int = 12) -> Dict:
        """
        Printability score — ocena realizowalności druku DLP.

        Sprawdza:
        1. Min wall thickness vs rozdzielczość drukarki
        2. Min channel width vs ryzyko bleed/zalania
        3. Trapped resin risk (zamknięte pory)
        4. Feature size statistics

        Args:
            voxels_per_uc: Rozdzielczość voxelizacji

        Returns:
            Dict z printability score i szczegółami
        """
        t0 = time.time()
        self._ensure_edt(voxels_per_uc)

        pixel_um = self.printer_pixel_um
        pixel_mm = pixel_um / 1000.0

        # ---- Wall thickness check ----
        solid_mask = self._voxel_grid
        edt_solid = self._solid_edt.copy()
        edt_solid[~solid_mask] = 0

        wall_values_um = edt_solid[solid_mask & (edt_solid > 0)] * 2.0 * 1000

        if len(wall_values_um) == 0:
            return {'error': 'No solid voxels found'}

        min_wall = float(np.min(wall_values_um))
        pct_thin_walls = float(np.mean(wall_values_um < pixel_um * 2) * 100)

        # ---- Channel width check ----
        void_mask = ~self._voxel_grid
        edt_void = self._void_edt.copy()
        edt_void[~void_mask] = 0

        channel_values_um = edt_void[void_mask & (edt_void > 0)] * 2.0 * 1000

        if len(channel_values_um) > 0:
            min_channel = float(np.min(channel_values_um))
            pct_narrow_channels = float(np.mean(channel_values_um < pixel_um * 3) * 100)
        else:
            min_channel = 0
            pct_narrow_channels = 100

        # ---- Trapped resin (isolated void) ----
        struct = generate_binary_structure(3, 3)
        labeled, n_components = label(void_mask, structure=struct)
        total_void = np.sum(void_mask)

        # Check which components touch the exterior
        # Exterior = any face of bounding box
        nz = self._voxel_grid.shape[2]
        ny = self._voxel_grid.shape[1]
        nx = self._voxel_grid.shape[0]

        exterior_labels = set()
        for face in [
            labeled[0, :, :], labeled[-1, :, :],
            labeled[:, 0, :], labeled[:, -1, :],
            labeled[:, :, 0], labeled[:, :, -1]
        ]:
            exterior_labels.update(set(np.unique(face)) - {0})

        trapped_voxels = 0
        for comp in range(1, n_components + 1):
            if comp not in exterior_labels:
                trapped_voxels += np.sum(labeled == comp)

        trapped_fraction = trapped_voxels / total_void if total_void > 0 else 0

        # ---- Score calculation ----
        issues = []
        score = 100  # Start at perfect

        # Wall thickness scoring
        if min_wall < pixel_um:
            score -= 40
            issues.append(f'FAIL: Min wall {min_wall:.0f} µm < pixel {pixel_um:.0f} µm')
        elif min_wall < pixel_um * 2:
            score -= 20
            issues.append(f'WARNING: Min wall {min_wall:.0f} µm < 2× pixel {pixel_um*2:.0f} µm')
        elif min_wall < pixel_um * 3:
            score -= 5
            issues.append(f'CAUTION: Min wall {min_wall:.0f} µm ≈ 3× pixel')

        if pct_thin_walls > 10:
            score -= 15
            issues.append(f'WARNING: {pct_thin_walls:.1f}% of walls thinner than 2 pixels')

        # Channel width scoring
        if min_channel < pixel_um * 2:
            score -= 30
            issues.append(f'FAIL: Min channel {min_channel:.0f} µm — will likely clog')
        elif min_channel < pixel_um * 3:
            score -= 15
            issues.append(f'WARNING: Min channel {min_channel:.0f} µm — bleed risk')

        if pct_narrow_channels > 20:
            score -= 10
            issues.append(f'WARNING: {pct_narrow_channels:.1f}% of channels < 3 pixels wide')

        # Trapped resin scoring
        if trapped_fraction > 0.1:
            score -= 20
            issues.append(f'FAIL: {trapped_fraction:.0%} trapped resin — cannot drain!')
        elif trapped_fraction > 0.02:
            score -= 10
            issues.append(f'WARNING: {trapped_fraction:.1%} trapped resin')

        score = max(0, score)

        # Rating
        if score >= 85:
            rating = 'EXCELLENT'
        elif score >= 70:
            rating = 'GOOD'
        elif score >= 50:
            rating = 'MARGINAL'
        elif score >= 30:
            rating = 'POOR'
        else:
            rating = 'FAIL'

        results = {
            'score': score,
            'rating': rating,
            'printer_pixel_um': pixel_um,

            'min_wall_um': min_wall,
            'pct_thin_walls': pct_thin_walls,
            'min_channel_um': min_channel,
            'pct_narrow_channels': pct_narrow_channels,
            'trapped_resin_fraction': float(trapped_fraction),
            'n_trapped_components': n_components - len(exterior_labels),

            'issues': issues,
            'calculation_time': time.time() - t0,
        }

        return results
