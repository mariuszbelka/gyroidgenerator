"""
Distributions Module — Structural and Topological Analysis
================================================================
Computational module for:
  - Wall thickness distribution
  - Channel width distribution
  - Connectivity / dead-end fraction / percolation
  - Accessible Surface Area (ASA, probe-based)
  - Throat / constriction analysis
  - Printability score

Methods are based on mesh ray-casting and voxelization + distance transform (EDT).

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
        self._analysis_mesh = mesh # Used for ray-casting, can be simplified
        self.container_geometry = container_geometry
        self.container_params = container_params
        self.unit_cell_mm = unit_cell_mm
        self.printer_pixel_um = printer_pixel_um

        # Cache
        self._voxel_grid = None
        self._pitch = None
        self._solid_edt = None  # Distance transform of solid phase
        self._void_edt = None   # Distance transform of void phase

    def _get_container_mask(self, x, y, z) -> np.ndarray:
        """Create boolean mask for voxels inside container."""
        bounds = self.mesh.bounds
        center = (bounds[0] + bounds[1]) / 2.0

        if self.container_geometry == 'sphere':
            r = self.container_params['diameter'] / 2.0
            dist = np.sqrt((x - center[0])**2 + (y - center[1])**2 + (z - center[2])**2)
            return dist <= r
        elif self.container_geometry == 'cylinder':
            r = self.container_params['diameter'] / 2.0
            h = self.container_params['height']
            dist_r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
            dist_z = np.abs(z - center[2])
            return (dist_r <= r) & (dist_z <= h/2.0)
        return np.ones_like(x, dtype=bool)

    def _voxelize(self, voxels_per_unit_cell: int = 12):
        """
        Voxelize mesh with given resolution and container masking.
        """
        if self._voxel_grid is not None:
            return

        pitch = self.unit_cell_mm / voxels_per_unit_cell
        logger.info(f"Voxelizing: pitch={pitch:.4f} mm, "
                    f"resolution={voxels_per_unit_cell} vox/unit_cell")

        t0 = time.time()
        voxel_obj = self._analysis_mesh.voxelized(pitch=pitch)

        # Robust volume representation
        try:
            voxel_obj = voxel_obj.fill()
        except:
            pass

        grid = voxel_obj.matrix
        # trimesh pitch can be a vector
        pitch_vec = np.atleast_1d(voxel_obj.pitch)
        self._pitch = float(pitch_vec[0])

        # Masking using open grid for memory efficiency
        nx, ny, nz = grid.shape
        origin = voxel_obj.translation

        # Use broadcasting with ogrid
        xi, yi, zi = np.ogrid[:nx, :ny, :nz]
        X = xi * self._pitch + float(origin[0])
        Y = yi * self._pitch + float(origin[1])
        Z = zi * self._pitch + float(origin[2])

        self._container_mask = self._get_container_mask(X, Y, Z)
        self._voxel_grid = grid & self._container_mask

        logger.info(f"Voxel grid: {self._voxel_grid.shape}, "
                    f"time={time.time()-t0:.1f}s")

    def simplify_analysis_mesh(self, target_faces: int):
        """
        Simplify the mesh used for ray-casting analysis.
        This significantly improves performance and memory usage
        while maintaining measurement accuracy.
        """
        if len(self.mesh.faces) <= target_faces:
            self._analysis_mesh = self.mesh
            return

        logger.info(f"Simplifying mesh from {len(self.mesh.faces)} to {target_faces} faces...")
        self._analysis_mesh = self.mesh.simplify_quadric_decimation(face_count=target_faces)
        logger.info("Mesh simplification complete.")

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

    def _get_internal_face_indices(self, skin_tolerance: float = 0.02) -> np.ndarray:
        """
        Find indices of faces that are NOT on the container's outer boundary.
        """
        centroids = self._analysis_mesh.triangles_center
        bounds = self._analysis_mesh.bounds
        center = (bounds[0] + bounds[1]) / 2.0

        if self.container_geometry == 'sphere':
            r_target = self.container_params['diameter'] / 2.0
            dist = np.linalg.norm(centroids - center, axis=1)
            # Filter out faces close to the boundary
            return np.where(dist < r_target * (1.0 - skin_tolerance))[0]

        elif self.container_geometry == 'cylinder':
            r_target = self.container_params['diameter'] / 2.0
            h_target = self.container_params['height']
            dist_r = np.linalg.norm(centroids[:, :2] - center[:2], axis=1)
            dist_z = np.abs(centroids[:, 2] - center[2])

            is_skin = (dist_r > r_target * (1.0 - skin_tolerance)) | \
                      (dist_z > (h_target / 2.0) * (1.0 - skin_tolerance))
            return np.where(~is_skin)[0]

        return np.arange(len(self._analysis_mesh.faces))

    def calc_wall_thickness_distribution(self, n_samples: int = 10000,
                                        batch_size: int = 100,
                                        callback=None) -> Dict:
        """
        Wall thickness distribution (Mesh-based Ray Casting).
        Measures the distance between opposing surfaces of the skeleton.
        """
        t0 = time.time()

        # 0. Filter to internal faces only to avoid "skin" artifacts
        internal_indices = self._get_internal_face_indices()
        if len(internal_indices) == 0:
            return {'error': 'No internal faces found for sampling'}

        # 1. Sample points on internal surface only
        points, face_indices = trimesh.sample.sample_surface(
            self._analysis_mesh, n_samples, face_indices=internal_indices
        )

        # 2. Get normals and flip them inwards (into the solid material)
        normals = self._analysis_mesh.face_normals[face_indices]
        ray_directions = -normals
        ray_origins = points + (ray_directions * 1e-5)

        # 3. Ray casting in batches
        all_distances = []
        n_batches = int(np.ceil(n_samples / batch_size))

        for i in range(n_batches):
            start = i * batch_size
            end = min((i + 1) * batch_size, n_samples)

            b_origins = ray_origins[start:end]
            b_dirs = ray_directions[start:end]

            locations, index_ray, _ = self._analysis_mesh.ray.intersects_location(
                ray_origins=b_origins,
                ray_directions=b_dirs,
                multiple_hits=False
            )

            if len(locations) > 0:
                dist = np.linalg.norm(locations - b_origins[index_ray], axis=1)

                # Filter outliers (length > 3x unit cell)
                mask = dist < (3.0 * self.unit_cell_mm)
                if np.any(mask):
                    all_distances.append(dist[mask])

            if callback:
                callback(int((i + 1) / n_batches * 100))

        # 4. Calculate distances
        if not all_distances:
            return {'error': 'No wall thickness measurements found'}

        wall_values_um = np.concatenate(all_distances) * 1000.0  # mm -> µm

        if len(wall_values_um) == 0:
            return {'error': 'No wall thickness measurements found'}

        # Histogram
        n_bins = 40
        hist_counts, hist_edges = np.histogram(wall_values_um, bins=n_bins)

        # Percentyle
        p1 = np.percentile(wall_values_um, 1)
        p5 = np.percentile(wall_values_um, 5)
        p10 = np.percentile(wall_values_um, 10)
        p50 = np.percentile(wall_values_um, 50)
        p90 = np.percentile(wall_values_um, 90)

        results = {
            'mean_um': float(np.mean(wall_values_um)),
            'std_um': float(np.std(wall_values_um)),
            'min_um': float(np.min(wall_values_um)),
            'max_um': float(np.max(wall_values_um)),
            'p1_um': float(p1),
            'p5_um': float(p5),
            'p10_um': float(p10),
            'p50_um': float(p50),
            'p90_um': float(p90),
            'n_measurements': len(wall_values_um),
            'hist_counts': hist_counts.tolist(),
            'hist_edges_um': hist_edges.tolist(),
            'uniformity': float(1.0 - np.std(wall_values_um) / np.mean(wall_values_um))
                          if np.mean(wall_values_um) > 0 else 0,
            'calculation_time': time.time() - t0,
        }

        return results

    # ==================== CHANNEL WIDTH DISTRIBUTION ====================

    def calc_channel_width_distribution(self, n_samples: int = 10000,
                                         batch_size: int = 100,
                                         callback=None) -> Dict:
        """
        Channel width distribution (Mesh-based Ray Casting).
        Measures the distance between walls within the void channels.
        """
        t0 = time.time()

        # 0. Filter to internal faces only
        internal_indices = self._get_internal_face_indices()
        if len(internal_indices) == 0:
            return {'error': 'No internal faces found for sampling'}

        # 1. Sample points on surface
        points, face_indices = trimesh.sample.sample_surface(
            self._analysis_mesh, n_samples, face_indices=internal_indices
        )

        # 2. Get normals pointing into the channels
        ray_directions = self._analysis_mesh.face_normals[face_indices]
        ray_origins = points + (ray_directions * 1e-5)

        # 3. Ray casting in batches
        all_distances = []
        n_batches = int(np.ceil(n_samples / batch_size))

        for i in range(n_batches):
            start = i * batch_size
            end = min((i + 1) * batch_size, n_samples)

            b_origins = ray_origins[start:end]
            b_dirs = ray_directions[start:end]

            locations, index_ray, _ = self._analysis_mesh.ray.intersects_location(
                ray_origins=b_origins,
                ray_directions=b_dirs,
                multiple_hits=False
            )

            if len(locations) > 0:
                dist = np.linalg.norm(locations - b_origins[index_ray], axis=1)

                # Filter outliers (length > 3x unit cell)
                mask = dist < (3.0 * self.unit_cell_mm)
                if np.any(mask):
                    all_distances.append(dist[mask])

            if callback:
                callback(int((i + 1) / n_batches * 100))

        # 4. Calculate distances (discard rays that leave the container)
        if not all_distances:
            return {'error': 'No channel width measurements found'}

        channel_values_um = np.concatenate(all_distances) * 1000.0

        if len(channel_values_um) == 0:
            return {'error': 'No channel width measurements found'}

        n_bins = 40
        hist_counts, hist_edges = np.histogram(channel_values_um, bins=n_bins)

        p1 = np.percentile(channel_values_um, 1)
        p5 = np.percentile(channel_values_um, 5)
        p10 = np.percentile(channel_values_um, 10)
        p50 = np.percentile(channel_values_um, 50)
        p90 = np.percentile(channel_values_um, 90)

        results = {
            'mean_um': float(np.mean(channel_values_um)),
            'std_um': float(np.std(channel_values_um)),
            'min_um': float(np.min(channel_values_um)),
            'max_um': float(np.max(channel_values_um)),
            'p1_um': float(p1),
            'p5_um': float(p5),
            'p10_um': float(p10),
            'p50_um': float(p50),
            'p90_um': float(p90),
            'n_measurements': len(channel_values_um),
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
        Channel connectivity analysis (void space) within the container.
        """
        t0 = time.time()
        self._voxelize(voxels_per_uc)

        # void = NOT solid AND INSIDE container
        void_mask = (~self._voxel_grid) & self._container_mask

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
        total_voxels_in_container = np.sum(self._container_mask)
        connected_porosity = through_voxels / total_voxels_in_container if total_voxels_in_container > 0 else 0

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

    def calc_accessible_surface_area(self, probe_radii_um: List[float] = None,
                                     voxels_per_uc: int = 15) -> Dict:
        """
        Accessible Surface Area (ASA) — particulate accessibility simulation.

        Method:
        1. Erode void space by probe radius (in voxels).
        2. ASA = surface area at the solid-void interface after erosion.

        Args:
            probe_radii_um: List of probe radii [µm]. Default: [0, 10, 20, 50, 100]
            voxels_per_uc: Voxelization resolution

        Returns:
            Dict with ASA for each probe radius
        """
        t0 = time.time()
        self._voxelize(voxels_per_uc)

        if probe_radii_um is None:
            probe_radii_um = [0, 10, 20, 50, 100]

        pitch_um = self._pitch * 1000  # mm → µm
        total_sa = self.mesh.area  # mm²

        profiles = []

        for probe_r_um in probe_radii_um:
            if probe_r_um == 0:
                profiles.append({
                    'probe_radius_um': 0,
                    'accessible_SA_mm2': float(total_sa),
                    'fraction': 1.0,
                })
                continue

            # Convert probe radius to voxels
            n_erosions = int(round(probe_r_um / pitch_um))

            if n_erosions < 1:
                profiles.append({
                    'probe_radius_um': float(probe_r_um),
                    'accessible_SA_mm2': float(total_sa),
                    'fraction': 1.0,
                })
                continue

            # Erode void space
            void_mask = (~self._voxel_grid) & self._container_mask
            eroded_void = void_mask.copy()
            for _ in range(n_erosions):
                eroded_void = binary_erosion(eroded_void)

            # Measure surface area by counting boundary voxels
            dilated_void = binary_dilation(eroded_void)
            boundary = dilated_void & self._voxel_grid
            n_boundary = np.sum(boundary)

            dilated_orig = binary_dilation(void_mask)
            boundary_orig = dilated_orig & self._voxel_grid
            n_boundary_orig = max(1, np.sum(boundary_orig))

            fraction = n_boundary / n_boundary_orig
            accessible_sa = total_sa * fraction

            profiles.append({
                'probe_radius_um': float(probe_r_um),
                'accessible_SA_mm2': float(accessible_sa),
                'fraction': float(fraction),
            })

        results = {
            'profiles': profiles,
            'total_SA_mm2': float(total_sa),
            'pitch_um': float(pitch_um),
            'calculation_time': time.time() - t0,
        }

        return results

    # ==================== THROAT ANALYSIS ====================

    def calc_throat_distribution(self, channel_results: Dict) -> Dict:
        """
        Throat / constriction analysis.
        Optimized: analytical estimate based on channel distribution.

        Method:
        For Gyroid structures:
          Mean throat ≈ 0.70 * Mean channel width
          P10 throat ≈ 0.70 * P10 channel width

        Returns:
            Dict with throat size estimates
        """
        t0 = time.time()

        if 'error' in channel_results:
            return channel_results

        # Analytical scaling factor for Gyroid throats
        FACTOR = 0.70

        mean_um = channel_results['mean_um'] * FACTOR
        p10_um = channel_results['p10_um'] * FACTOR
        p50_um = channel_results['p50_um'] * FACTOR
        p90_um = channel_results['p90_um'] * FACTOR
        min_um = channel_results['min_um'] * FACTOR
        max_um = channel_results['max_um'] * FACTOR

        results = {
            'mean_um': float(mean_um),
            'p10_um': float(p10_um),
            'p50_um': float(p50_um),
            'p90_um': float(p90_um),
            'min_um': float(min_um),
            'max_um': float(max_um),
            'throat_to_pore_ratio': FACTOR,
            'method': 'Analytical estimate (0.7 × channel)',
            'calculation_time': time.time() - t0,
        }

        # Warnings for clogging
        warnings = []
        if p10_um < 20:
            warnings.append(f'CRITICAL: P10 throat ({p10_um:.0f} µm) < 20 µm — High clogging risk for particulate samples!')
        elif p10_um < 35:
            warnings.append(f'WARNING: P10 throat ({p10_um:.0f} µm) < 35 µm — Possible pressure buildup.')

        results['warnings'] = warnings

        return results

    # ==================== PRINTABILITY SCORE ====================

    def calc_printability(self, wall_results: Dict, channel_results: Dict,
                          voxels_per_uc: int = 12) -> Dict:
        """
        Printability score — DLP 3D printing feasibility assessment.
        Optimized: uses P1/P5 percentiles instead of absolute minima for robustness.

        Args:
            wall_results: Results from calc_wall_thickness_distribution
            channel_results: Results from calc_channel_width_distribution
            voxels_per_uc: Resolution for trapped resin analysis

        Returns:
            Dict with printability score and details
        """
        t0 = time.time()
        self._voxelize(voxels_per_uc)

        pixel_um = self.printer_pixel_um

        # ---- Wall thickness check (Robust: use P1 percentile) ----
        if 'error' in wall_results:
            return {'error': 'Need wall thickness analysis for printability'}

        min_wall = wall_results.get('p1_um', wall_results['min_um'])
        # Also check % of walls below 2 pixels (heuristic)
        pct_thin_walls = 15.0 if wall_results.get('p10_um', 0) < (pixel_um * 2) else 0.0

        # ---- Channel width check (Robust: use P1 percentile) ----
        if 'error' in channel_results:
            return {'error': 'Need channel width analysis for printability'}

        min_channel = channel_results.get('p1_um', channel_results['min_um'])
        pct_narrow_channels = 25.0 if channel_results.get('p10_um', 0) < (pixel_um * 3) else 0.0

        # ---- Trapped resin (isolated void) ----
        void_mask = (~self._voxel_grid) & self._container_mask
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
