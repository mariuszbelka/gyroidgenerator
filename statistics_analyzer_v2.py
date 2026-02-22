"""
Gyroid Generator - Statistics Analyzer v2.0
Modular statistics system with user-selectable metrics

Author: Claude (Anthropic)
Date: 2026-02-19
"""

import numpy as np
import trimesh
from scipy.ndimage import distance_transform_edt
from typing import Dict, List, Tuple, Optional
import time


class StatisticsCalculator:
    """
    Modular statistics calculator for gyroid meshes.
    Each metric can be calculated independently.
    """

    def __init__(self, mesh: trimesh.Trimesh,
                 container_geometry: str,
                 container_params: Dict):
        """
        Initialize statistics calculator.

        Args:
            mesh: Trimesh object (gyroid structure)
            container_geometry: 'cylinder', 'sphere', etc.
            container_params: Dict with geometry parameters
                For cylinder: {'diameter': float, 'height': float}
                For sphere: {'diameter': float}
        """
        self.mesh = mesh
        self.container_geometry = container_geometry
        self.container_params = container_params

        # Cache for voxel grid (expensive to compute)
        self._voxel_grid = None
        self._voxel_resolution = None
        self._solid_voxels = None
        self._void_voxels = None


    # ==================== INSTANT METRICS ====================

    def calculate_basic_properties(self) -> Dict:
        """
        Calculate basic mesh properties (instant).

        Returns:
            Dict with: vertices, faces, watertight, volume, surface_area
        """
        start_time = time.time()

        results = {
            'vertices': len(self.mesh.vertices),
            'faces': len(self.mesh.faces),
            'watertight': self.mesh.is_watertight,
            'surface_area': self.mesh.area,  # mm²
        }

        # Calculate total container volume (always possible)
        total_volume = self._calculate_container_volume()
        results['total_volume'] = total_volume

        # Volume calculation (only if watertight)
        if self.mesh.is_watertight:
            results['gyroid_volume'] = self.mesh.volume  # mm³
            results['void_volume'] = total_volume - results['gyroid_volume'] if total_volume else None
            results['porosity'] = (results['void_volume'] / total_volume) * 100 if total_volume and results['void_volume'] is not None else None
        else:
            results['gyroid_volume'] = None
            results['void_volume'] = None
            results['porosity'] = None

        results['calculation_time'] = time.time() - start_time
        return results


    def calculate_specific_surface_area(self, basic_props: Dict,
                                       density: float = 1.05) -> Dict:
        """
        Calculate specific surface area (instant).

        Args:
            basic_props: Output from calculate_basic_properties()
            density: Material density in g/cm³

        Returns:
            Dict with SSA_volumetric (m²/mL), SSA_gravimetric (m²/g), mass (g)
        """
        start_time = time.time()

        surface_area_mm2 = basic_props['surface_area']
        surface_area_m2 = surface_area_mm2 / 1e6  # mm² → m²

        results = {}

        if basic_props['total_volume'] is not None:
            total_volume_mm3 = basic_props['total_volume']
            total_volume_mL = total_volume_mm3 / 1000  # mm³ → mL

            # Volumetric SSA
            results['SSA_volumetric'] = surface_area_m2 / total_volume_mL  # m²/mL

            # Mass calculation
            gyroid_volume_cm3 = basic_props['gyroid_volume'] / 1000  # mm³ → cm³
            mass_g = gyroid_volume_cm3 * density
            results['mass_g'] = mass_g

            # Gravimetric SSA
            results['SSA_gravimetric'] = surface_area_m2 / mass_g  # m²/g
        else:
            results['SSA_volumetric'] = None
            results['SSA_gravimetric'] = None
            results['mass_g'] = None

        results['density_g_cm3'] = density
        results['calculation_time'] = time.time() - start_time
        return results


    def calculate_external_surface_area(self, atol: float = 0.05) -> Dict:
        """
        Detect faces on the container boundary to calculate real external area.
        Takes into account "holes" in the gyroid at the surface.

        Args:
            atol: Absolute tolerance for boundary detection (mm)

        Returns:
            Dict with:
              - external_area_mm2: Actual area on boundary (with holes)
              - analytical_container_area_mm2: Full area of container without holes
              - internal_area_mm2: Total - actual external
              - external_fraction_percent
        """
        start_time = time.time()

        # 1. Analytical area
        analytical_area = 0
        if self.container_geometry == 'cylinder':
            r = self.container_params['diameter'] / 2
            h = self.container_params['height']
            analytical_area = 2 * np.pi * r * h + 2 * np.pi * r**2
        elif self.container_geometry == 'sphere':
            r = self.container_params['diameter'] / 2
            analytical_area = 4 * np.pi * r**2
        elif self.container_geometry == 'cube':
            s = self.container_params.get('size', 10.0)
            analytical_area = 6 * s**2

        # 2. Detect boundary faces
        face_centroids = self.mesh.triangles_center
        face_areas = self.mesh.area_faces

        is_external = np.zeros(len(face_centroids), dtype=bool)

        if self.container_geometry == 'cylinder':
            r_target = self.container_params['diameter'] / 2
            # Distance from Z axis (assuming centered at XY=0)
            r_actual = np.sqrt(face_centroids[:, 0]**2 + face_centroids[:, 1]**2)

            # Side wall
            on_side = np.abs(r_actual - r_target) < atol
            # Caps
            z_min, z_max = self.mesh.bounds[:, 2]
            on_caps = (np.abs(face_centroids[:, 2] - z_min) < atol) | \
                      (np.abs(face_centroids[:, 2] - z_max) < atol)

            is_external = on_side | on_caps

        elif self.container_geometry == 'sphere':
            r_target = self.container_params['diameter'] / 2
            center = self.mesh.centroid
            dist_from_center = np.linalg.norm(face_centroids - center, axis=1)
            is_external = np.abs(dist_from_center - r_target) < atol

        elif self.container_geometry == 'cube':
            bounds = self.mesh.bounds
            for i in range(3): # x, y, z
                on_bound = (np.abs(face_centroids[:, i] - bounds[0, i]) < atol) | \
                           (np.abs(face_centroids[:, i] - bounds[1, i]) < atol)
                is_external |= on_bound

        actual_external_area = np.sum(face_areas[is_external])
        total_area = self.mesh.area
        internal_area = total_area - actual_external_area

        results = {
            'external_area_mm2': float(actual_external_area),
            'analytical_container_area_mm2': float(analytical_area),
            'internal_area_mm2': float(internal_area),
            'total_area_mm2': float(total_area),
            'external_fraction_percent': float(actual_external_area / total_area * 100) if total_area > 0 else 0,
            'calculation_time': time.time() - start_time
        }
        return results


    # ==================== THEORETICAL ESTIMATES ====================

    def estimate_wall_thickness_theoretical(self, unit_cell: float,
                                           wall_design: float) -> Dict:
        """
        Theoretical wall thickness from design parameters (instant).

        Args:
            unit_cell: Unit cell size in mm
            wall_design: Design wall thickness in mm

        Returns:
            Dict with wall_thickness_mm, accuracy_note
        """
        start_time = time.time()

        results = {
            'wall_thickness_mm': wall_design,
            'wall_thickness_um': wall_design * 1000,
            'method': 'theoretical',
            'accuracy': 'Design parameter (actual may vary ±10%)',
            'calculation_time': time.time() - start_time
        }
        return results


    def estimate_channel_width_theoretical(self, unit_cell: float,
                                          wall_thickness: float) -> Dict:
        """
        Theoretical channel width (instant).

        Args:
            unit_cell: Unit cell size in mm
            wall_thickness: Wall thickness in mm

        Returns:
            Dict with channel_width_mm
        """
        start_time = time.time()

        # Simple approximation: channel = unit_cell - wall
        channel_width = unit_cell - wall_thickness

        results = {
            'channel_width_mm': channel_width,
            'channel_width_um': channel_width * 1000,
            'method': 'theoretical',
            'formula': 'unit_cell - wall_thickness',
            'calculation_time': time.time() - start_time
        }
        return results




    # ==================== VOXEL-BASED ANALYSIS ====================

    def _voxelize_mesh(self, voxels_per_unit_cell: int = 10):
        """
        Convert mesh to voxel grid (expensive - cached).

        Args:
            voxels_per_unit_cell: Resolution (higher = more accurate, slower)
        """
        if self._voxel_grid is not None and self._voxel_resolution == voxels_per_unit_cell:
            return  # Already computed

        print(f"INFO: Voxelizing mesh (resolution: {voxels_per_unit_cell} voxels/unit)...")

        # Determine voxel size based on container dimensions
        if self.container_geometry == 'cylinder':
            max_dim = max(self.container_params['diameter'],
                         self.container_params['height'])
        elif self.container_geometry == 'sphere':
            max_dim = self.container_params['diameter']
        else:
            max_dim = self.mesh.bounds[1].max() - self.mesh.bounds[0].min()

        # Estimate unit cell from dimensions (rough)
        estimated_unit_cells = max_dim / 0.2  # Assume ~0.2mm unit cell
        n_voxels = int(estimated_unit_cells * voxels_per_unit_cell)
        n_voxels = max(50, min(n_voxels, 500))  # Clamp to reasonable range

        # Voxelize using trimesh
        voxel_grid = self.mesh.voxelized(pitch=max_dim/n_voxels)

        # Convert to numpy boolean array (True = solid, False = void)
        self._voxel_grid = voxel_grid.matrix
        self._voxel_resolution = voxels_per_unit_cell
        self._voxel_pitch = max_dim / n_voxels

        # Separate phases
        self._solid_voxels = self._voxel_grid
        self._void_voxels = ~self._voxel_grid

        print(f"INFO: Voxel grid created: {self._voxel_grid.shape}, pitch={self._voxel_pitch:.4f}mm")


    def calculate_desorption_path_fast(self) -> Dict:
        """
        Fast desorption path: distance from each solid voxel to nearest void voxel.
        Reflects diffusion distance within the polymer wall.

        Returns:
            Dict with statistics and distribution
        """
        start_time = time.time()

        # Ensure voxel grid exists
        self._voxelize_mesh(voxels_per_unit_cell=12)

        print("INFO: Calculating desorption paths (solid-to-void EDT)...")

        # EDT on solid: foreground=solid (True), background=void (False)
        # scipy.ndimage.distance_transform_edt(input) calculates distance
        # of each non-zero (True) voxel to the nearest zero (False) voxel.
        edt = distance_transform_edt(self._solid_voxels, sampling=self._voxel_pitch)

        # Extract distances only for solid voxels
        distances = edt[self._solid_voxels]

        if len(distances) == 0:
            return {'error': 'No solid voxels found'}

        # Statistics
        mean_path = np.mean(distances)
        max_path = np.max(distances)
        min_path = np.min(distances)
        std_path = np.std(distances)
        median_path = np.median(distances)

        # Histogram
        hist, bins = np.histogram(distances, bins=50)

        results = {
            'mean_desorption_path_mm': float(mean_path),
            'median_desorption_path_mm': float(median_path),
            'max_desorption_path_mm': float(max_path),
            'min_desorption_path_mm': float(min_path),
            'std_desorption_path_mm': float(std_path),
            'n_samples': len(distances),
            'method': 'solid_to_void_edt',
            'histogram_counts': hist.tolist(),
            'histogram_bins': bins.tolist(),
            'calculation_time': time.time() - start_time
        }

        return results


    def calculate_desorption_path_sampling(self, n_samples: int = 1000) -> Dict:
        """
        Statistical sampling method for desorption path using mesh proximity.
        Directly measures distance to nearest phase boundary.

        Args:
            n_samples: Number of random points to sample from solid phase

        Returns:
            Dict with statistics and distribution
        """
        start_time = time.time()

        # 1. Sample points inside the solid phase
        # Use voxelization to find solid regions for efficient sampling
        self._voxelize_mesh(voxels_per_unit_cell=12)
        solid_coords = np.argwhere(self._solid_voxels)

        if len(solid_coords) == 0:
            return {'error': 'No solid voxels found'}

        if len(solid_coords) > n_samples:
            indices = np.random.choice(len(solid_coords), n_samples, replace=False)
            sampled_voxels = solid_coords[indices]
        else:
            sampled_voxels = solid_coords
            n_samples = len(solid_coords)

        # Convert to physical coordinates and add small random jitter
        # within voxel to avoid grid alignment artifacts
        jitter = (np.random.rand(n_samples, 3) - 0.5) * self._voxel_pitch
        sampled_points = (sampled_voxels * self._voxel_pitch) + jitter

        print(f"INFO: Querying mesh proximity for {n_samples} points...")
        from trimesh.proximity import closest_point

        _, distances, _ = closest_point(self.mesh, sampled_points)

        # Statistics
        mean_path = np.mean(distances)
        max_path = np.max(distances)
        min_path = np.min(distances)
        std_path = np.std(distances)
        median_path = np.median(distances)

        # Histogram
        hist, bins = np.histogram(distances, bins=50)

        results = {
            'mean_desorption_path_mm': float(mean_path),
            'median_desorption_path_mm': float(median_path),
            'max_desorption_path_mm': float(max_path),
            'min_desorption_path_mm': float(min_path),
            'std_desorption_path_mm': float(std_path),
            'n_samples': n_samples,
            'method': 'mesh_proximity_sampling',
            'histogram_counts': hist.tolist(),
            'histogram_bins': bins.tolist(),
            'calculation_time': time.time() - start_time
        }

        return results


    # ==================== HELPER METHODS ====================

    def _calculate_container_volume(self) -> float:
        """Calculate total volume of container."""
        if self.container_geometry == 'cylinder':
            r = self.container_params['diameter'] / 2
            h = self.container_params['height']
            return np.pi * r**2 * h

        elif self.container_geometry == 'sphere':
            r = self.container_params['diameter'] / 2
            return (4/3) * np.pi * r**3

        return None


# ==================== CONVENIENCE WRAPPER ====================

class StatisticsAnalyzer:
    """
    High-level interface for statistics calculation.
    Manages StatisticsCalculator and provides user-friendly API.
    """

    def __init__(self, mesh: trimesh.Trimesh,
                 container_geometry: str,
                 container_params: Dict,
                 unit_cell: float,
                 wall_thickness: float):
        """
        Initialize analyzer.

        Args:
            mesh: Trimesh object
            container_geometry: 'cylinder', 'sphere', etc.
            container_params: Container dimensions
            unit_cell: Unit cell size in mm (for theoretical calcs)
            wall_thickness: Wall thickness in mm (for theoretical calcs)
        """
        self.calculator = StatisticsCalculator(mesh, container_geometry, container_params)
        self.unit_cell = unit_cell
        self.wall_thickness = wall_thickness

        # Cache for results
        self.results = {}


    def calculate_selected_metrics(self, selections: Dict) -> Dict:
        """
        Calculate only selected metrics.

        Args:
            selections: Dict of boolean flags

        Returns:
            Dict with all calculated results
        """
        results = {}
        total_time_start = time.time()

        # Basic properties (usually needed for other metrics)
        if selections.get('basic_properties', True):
            print("Calculating: Basic Properties...")
            results['basic_properties'] = self.calculator.calculate_basic_properties()

        # SSA (depends on basic properties)
        if selections.get('specific_surface_area', True):
            if 'basic_properties' in results:
                print("Calculating: Specific Surface Area...")
                density = selections.get('density_g_cm3', 1.05)
                results['specific_surface_area'] = \
                    self.calculator.calculate_specific_surface_area(
                        results['basic_properties'], density
                    )

        # External surface
        if selections.get('external_surface', True):
            print("Calculating: External Surface Area...")
            results['external_surface'] = \
                self.calculator.calculate_external_surface_area()

        # Theoretical estimates
        if selections.get('wall_thickness_theoretical', True):
            print("Calculating: Wall Thickness (Theoretical)...")
            results['wall_thickness_theoretical'] = \
                self.calculator.estimate_wall_thickness_theoretical(
                    self.unit_cell, self.wall_thickness
                )

        if selections.get('channel_width_theoretical', True):
            print("Calculating: Channel Width (Theoretical)...")
            results['channel_width_theoretical'] = \
                self.calculator.estimate_channel_width_theoretical(
                    self.unit_cell, self.wall_thickness
                )

        # Voxel-based analysis (slower)
        if selections.get('desorption_path_fast', False):
            print("Calculating: Desorption Path (Fast)...")
            results['desorption_path_fast'] = \
                self.calculator.calculate_desorption_path_fast()

        if selections.get('desorption_path_sampling', False):
            print("Calculating: Desorption Path (Sampling)...")
            n_samples = selections.get('desorption_sampling_n', 1000)
            results['desorption_path_sampling'] = \
                self.calculator.calculate_desorption_path_sampling(n_samples)

        results['total_calculation_time'] = time.time() - total_time_start

        return results


    def export_to_csv(self, results: Dict, filename: str):
        """
        Export results to CSV file.

        Args:
            results: Output from calculate_selected_metrics()
            filename: Output file path
        """
        import csv

        with open(filename, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['Metric', 'Value', 'Unit', 'Method'])

            # Flatten results dict
            for category, data in results.items():
                if category == 'total_calculation_time':
                    writer.writerow([category, f"{data:.2f}", 's', ''])
                    continue

                if isinstance(data, dict):
                    for key, value in data.items():
                        if isinstance(value, (int, float)):
                            writer.writerow([f"{category}.{key}", f"{value:.4f}", '', ''])
                        elif isinstance(value, str):
                            writer.writerow([f"{category}.{key}", value, '', ''])

        print(f"Results exported to: {filename}")


# ============================================================
# Wrapper classes for gui_main.py / statistics_tab_v2.py API
# ============================================================

class GeometryParams:
    """Simple data holder for container geometry."""
    def __init__(self, shape: str, dimensions: Dict):
        self.shape = shape
        self.dimensions = dimensions


class GyroidParams:
    """Simple data holder for gyroid design parameters."""
    def __init__(self, unit_cell: float, wall_thickness: float, threshold: float = 0.0):
        self.unit_cell = unit_cell
        self.wall_thickness = wall_thickness
        self.threshold = threshold


class StatisticsV2:
    """
    Unified statistics interface for v2/v3 GUI.

    Wraps StatisticsCalculator and transforms return keys
    to match what statistics_tab_v2.py display_results() expects.
    """

    def __init__(self, mesh: trimesh.Trimesh,
                 geometry: GeometryParams,
                 gyroid_params: GyroidParams):
        self.mesh = mesh
        self.geometry = geometry
        self.gyroid_params = gyroid_params

        # Create underlying calculator
        self._calc = StatisticsCalculator(
            mesh=mesh,
            container_geometry=geometry.shape,
            container_params=geometry.dimensions
        )

        # Cache basic props
        self._basic_props = None

    def _ensure_basic_props(self):
        if self._basic_props is None:
            self._basic_props = self._calc.calculate_basic_properties()
        return self._basic_props

    def get_basic_stats(self) -> Dict:
        """Get basic mesh statistics (instant)."""
        bp = self._ensure_basic_props()
        return {
            'vertices': bp['vertices'],
            'faces': bp['faces'],
            'watertight': bp['watertight'],
            'surface_area': bp['surface_area'],
            'total_volume': bp.get('total_volume', 0),
            'gyroid_volume': bp.get('gyroid_volume', 0),
            'void_volume': bp.get('void_volume', 0),
            'porosity': bp.get('porosity'),
        }

    def calc_specific_surface_area(self, density: float = 1.05) -> Dict:
        """Calculate SSA. Transforms keys: SSA_volumetric → ssa_volumetric."""
        bp = self._ensure_basic_props()
        raw = self._calc.calculate_specific_surface_area(bp, density=density)
        return {
            'ssa_volumetric': raw.get('SSA_volumetric'),
            'ssa_gravimetric': raw.get('SSA_gravimetric'),
            'mass_g': raw.get('mass_g'),
            'calculation_time': raw.get('calculation_time', 0),
        }

    def calc_wall_thickness_theoretical(self) -> Dict:
        """Transforms: wall_thickness_um → mean_wall_thickness."""
        raw = self._calc.estimate_wall_thickness_theoretical(
            unit_cell=self.gyroid_params.unit_cell,
            wall_design=self.gyroid_params.wall_thickness
        )
        return {
            'mean_wall_thickness': raw.get('wall_thickness_um', 0),
            'method': raw.get('method', 'theoretical'),
            'accuracy': raw.get('accuracy', ''),
            'calculation_time': raw.get('calculation_time', 0),
        }

    def calc_channel_width_theoretical(self) -> Dict:
        """Transforms: channel_width_um → mean_channel_width."""
        raw = self._calc.estimate_channel_width_theoretical(
            unit_cell=self.gyroid_params.unit_cell,
            wall_thickness=self.gyroid_params.wall_thickness
        )
        return {
            'mean_channel_width': raw.get('channel_width_um', 0),
            'method': raw.get('method', 'theoretical'),
            'formula': raw.get('formula', ''),
            'calculation_time': raw.get('calculation_time', 0),
        }


    def calc_mass(self, density: float = 1.05) -> Dict:
        """Calculate mass from volume and density."""
        bp = self._ensure_basic_props()
        start_time = time.time()

        results = {'density_g_cm3': density}

        if bp.get('gyroid_volume') is not None:
            volume_cm3 = bp['gyroid_volume'] / 1000  # mm³ → cm³
            mass_g = volume_cm3 * density
            results['mass_g'] = mass_g
            results['mass_mg'] = mass_g * 1000
        else:
            results['mass_g'] = None
            results['mass_mg'] = None

        results['calculation_time'] = time.time() - start_time
        return results

    def calc_external_surface_area(self) -> Dict:
        """Adds internal/total/fraction keys for display."""
        raw = self._calc.calculate_external_surface_area()
        bp = self._ensure_basic_props()

        external = raw.get('external_area_mm2', 0) or 0
        analytical = raw.get('analytical_container_area_mm2', 0) or 0
        total_sa = bp.get('surface_area', 0)
        internal = max(0, total_sa - external)

        return {
            'external_area_mm2': external,
            'analytical_area_mm2': analytical,
            'internal_area_mm2': internal,
            'total_area_mm2': total_sa,
            'external_fraction_percent': (external / total_sa * 100) if total_sa > 0 else 0,
            'calculation_time': raw.get('calculation_time', 0),
        }

    def calc_desorption_path_fast(self) -> Dict:
        """Transforms desorption path keys: mean_desorption_path_mm → mean_path_mm."""
        raw = self._calc.calculate_desorption_path_fast()

        if 'error' in raw:
            return raw

        return {
            'mean_path_mm': raw.get('mean_desorption_path_mm', 0),
            'median_path_mm': raw.get('median_desorption_path_mm', 0),
            'std_path_mm': raw.get('std_desorption_path_mm', 0),
            'min_path_mm': max(0, raw.get('min_desorption_path_mm', 0)),
            'max_path_mm': raw.get('max_desorption_path_mm', 0),
            'samples': raw.get('n_samples', 0),
            'method': 'Voxel-based EDT (Solid-to-Void)',
            'accuracy': 'Full voxel coverage',
            'histogram_counts': raw.get('histogram_counts'),
            'histogram_bins': raw.get('histogram_bins'),
            'calculation_time': raw.get('calculation_time', 0),
        }

    def calc_desorption_path_sampling(self, n_samples: int = 1000) -> Dict:
        """Transforms desorption path keys for sampling method."""
        raw = self._calc.calculate_desorption_path_sampling(
            n_samples=n_samples
        )

        if 'error' in raw:
            return raw

        return {
            'mean_path_mm': raw.get('mean_desorption_path_mm', 0),
            'median_path_mm': raw.get('median_desorption_path_mm', 0),
            'std_path_mm': raw.get('std_desorption_path_mm', 0),
            'min_path_mm': max(0, raw.get('min_desorption_path_mm', 0)),
            'max_path_mm': raw.get('max_desorption_path_mm', 0),
            'samples': raw.get('n_samples', n_samples),
            'method': 'Mesh-based Proximity (Monte Carlo)',
            'accuracy': f'Direct mesh query ({n_samples} samples)',
            'histogram_counts': raw.get('histogram_counts'),
            'histogram_bins': raw.get('histogram_bins'),
            'calculation_time': raw.get('calculation_time', 0),
        }

    def export_summary_dict(self, results: Dict) -> Dict:
        """Flatten results into a simple key→value dict for CSV export."""
        summary = {}

        # Add design parameters
        summary['Shape'] = self.geometry.shape
        summary['Unit Cell [mm]'] = self.gyroid_params.unit_cell
        summary['Wall Thickness [mm]'] = self.gyroid_params.wall_thickness

        # Add basic stats
        bp = self._ensure_basic_props()
        summary['Vertices'] = bp['vertices']
        summary['Faces'] = bp['faces']
        summary['Watertight'] = bp['watertight']
        if bp.get('porosity') is not None:
            summary['Porosity [%]'] = f"{bp['porosity']:.1f}"
            summary['Total Volume [mm³]'] = f"{bp['total_volume']:.2f}"
            summary['Gyroid Volume [mm³]'] = f"{bp['gyroid_volume']:.2f}"
        summary['Surface Area [mm²]'] = f"{bp['surface_area']:.2f}"

        # Flatten calculation results
        for name, data in results.items():
            if isinstance(data, dict):
                for key, value in data.items():
                    if key in ('calculation_time', 'histogram_counts', 'histogram_bins'):
                        continue
                    if isinstance(value, float):
                        summary[f"{name}: {key}"] = f"{value:.4f}"
                    elif isinstance(value, (int, bool, str)):
                        summary[f"{name}: {key}"] = str(value)
            elif isinstance(data, (int, float, str)):
                summary[name] = str(data)

        return summary
