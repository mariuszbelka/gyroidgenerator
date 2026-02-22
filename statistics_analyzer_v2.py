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
        
        # Volume calculation (only if watertight)
        if self.mesh.is_watertight:
            results['gyroid_volume'] = self.mesh.volume  # mm³
            
            # Calculate total container volume
            total_volume = self._calculate_container_volume()
            results['total_volume'] = total_volume
            results['void_volume'] = total_volume - results['gyroid_volume']
            results['porosity'] = (results['void_volume'] / total_volume) * 100  # %
        else:
            results['gyroid_volume'] = None
            results['total_volume'] = None
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
    
    
    def calculate_external_surface_area(self) -> Dict:
        """
        Calculate external surface area of container (instant).
        
        Returns:
            Dict with external_area (mm²)
        """
        start_time = time.time()
        
        if self.container_geometry == 'cylinder':
            r = self.container_params['diameter'] / 2
            h = self.container_params['height']
            # Cylinder: 2πrh + 2πr²
            external_area = 2 * np.pi * r * h + 2 * np.pi * r**2
        
        elif self.container_geometry == 'sphere':
            r = self.container_params['diameter'] / 2
            # Sphere: 4πr²
            external_area = 4 * np.pi * r**2
        
        else:
            external_area = None
        
        results = {
            'external_area_mm2': external_area,
            'calculation_time': time.time() - start_time
        }
        return results
    
    
    def calculate_tortuosity(self) -> Dict:
        """
        Calculate theoretical tortuosity for gyroid (instant).
        
        Returns:
            Dict with tortuosity, effective_diffusivity_factor
        """
        start_time = time.time()
        
        # Gyroid theoretical tortuosity
        tau = 1.41
        
        # Effective diffusivity: D_eff = D_free / τ²
        D_eff_factor = 1.0 / (tau ** 2)
        
        results = {
            'tortuosity': tau,
            'effective_diffusivity_factor': D_eff_factor,
            'note': 'Theoretical value for gyroid TPMS',
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
    
    
    def estimate_desorption_path_theoretical(self) -> Dict:
        """
        Theoretical desorption path using geometric approximation (instant).
        
        Returns:
            Dict with mean_path_mm, max_path_mm
        """
        start_time = time.time()
        
        if self.container_geometry == 'cylinder':
            r = self.container_params['diameter'] / 2
            # Approximation: mean ≈ R/2, max ≈ R
            mean_path = r / 2
            max_path = r
        
        elif self.container_geometry == 'sphere':
            r = self.container_params['diameter'] / 2
            # For sphere: similar approximation
            mean_path = r / 2
            max_path = r
        
        else:
            mean_path = None
            max_path = None
        
        results = {
            'mean_desorption_path_mm': mean_path,
            'max_desorption_path_mm': max_path,
            'method': 'theoretical_geometric',
            'accuracy': 'Approximate (±20%)',
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
        self._voxel_origin = voxel_grid.translation
        self._voxel_pitch = voxel_grid.pitch
        self._voxel_resolution = voxels_per_unit_cell
        
        # Separate phases
        self._solid_voxels = self._voxel_grid
        self._void_voxels = ~self._voxel_grid
        
        print(f"INFO: Voxel grid created: {self._voxel_grid.shape}, pitch={self._voxel_pitch} mm")
    
    
    def calculate_desorption_path_fast(self, tau: float = 1.41) -> Dict:
        """
        Fast desorption path using Euclidean distance + tortuosity correction.
        
        Args:
            tau: Tortuosity factor (default 1.41 for gyroid)
        
        Returns:
            Dict with mean_path, max_path, std_dev, distribution
        """
        start_time = time.time()
        
        # Ensure voxel grid exists
        self._voxelize_mesh(voxels_per_unit_cell=10)
        
        print("INFO: Calculating desorption paths (fast method)...")
        
        # Get coordinates of solid voxels
        solid_coords = np.argwhere(self._solid_voxels)
        
        if len(solid_coords) == 0:
            return {'error': 'No solid voxels found'}
        
        # Convert voxel coords to physical coords (mm)
        # origin + indices * pitch
        solid_coords_mm = self._voxel_origin + solid_coords * self._voxel_pitch

        # Container parameters (using absolute coordinates from container_params)
        # Note: We assume the container is centered at the mesh center or as defined
        # In this platform, shapes are usually centered at (0,0,cz) or (0,0,0)
        # We should use the same logic as in mesh_generator.py

        # For sphere, we assume it's centered at (0,0,0) as per SphereGeometry default
        # For cylinder, centered at (0,0, height/2)
        
        if self.container_geometry == 'cylinder':
            center_x, center_y = 0.0, 0.0
            height = self.container_params['height']
            center_z = height / 2.0
            radius = self.container_params['diameter'] / 2
        elif self.container_geometry == 'sphere':
            center_x, center_y, center_z = 0.0, 0.0, 0.0
            radius = self.container_params['diameter'] / 2
        
        # Calculate Euclidean distance to edge for each solid voxel
        distances = []
        
        for coord in solid_coords_mm:
            if self.container_geometry == 'cylinder':
                # Distance to cylindrical wall
                dx = coord[0] - center_x
                dy = coord[1] - center_y
                radial_dist = np.sqrt(dx**2 + dy**2)
                dist_to_wall = radius - radial_dist
                
                # Distance to top/bottom caps
                dist_to_top = (self.container_params['height'] - coord[2])
                dist_to_bottom = coord[2]
                
                # Minimum distance to any boundary
                min_dist = min(dist_to_wall, dist_to_top, dist_to_bottom)
            
            elif self.container_geometry == 'sphere':
                # Distance to sphere surface
                dx = coord[0] - center_x
                dy = coord[1] - center_y
                dz = coord[2] - center_z
                dist_from_center = np.sqrt(dx**2 + dy**2 + dz**2)
                min_dist = radius - dist_from_center
            
            distances.append(min_dist)
        
        distances = np.array(distances)
        
        # Apply tortuosity correction
        distances_corrected = distances * tau
        
        # Statistics
        mean_path = np.mean(distances_corrected)
        max_path = np.max(distances_corrected)
        min_path = np.min(distances_corrected)
        std_path = np.std(distances_corrected)
        median_path = np.median(distances_corrected)
        
        # Histogram for distribution
        hist, bins = np.histogram(distances_corrected, bins=50)
        
        results = {
            'mean_desorption_path_mm': mean_path,
            'median_desorption_path_mm': median_path,
            'max_desorption_path_mm': max_path,
            'min_desorption_path_mm': min_path,
            'std_desorption_path_mm': std_path,
            'tortuosity_applied': tau,
            'n_samples': len(distances),
            'method': 'euclidean_corrected',
            'histogram_counts': hist.tolist(),
            'histogram_bins': bins.tolist(),
            'calculation_time': time.time() - start_time
        }
        
        return results
    
    
    def calculate_desorption_path_sampling(self, n_samples: int = 1000,
                                          tau: float = 1.41) -> Dict:
        """
        Statistical sampling method for desorption path (more accurate).
        
        Args:
            n_samples: Number of random points to sample from solid phase
            tau: Tortuosity factor
        
        Returns:
            Dict with statistics and distribution
        """
        start_time = time.time()
        
        # Ensure voxel grid exists
        self._voxelize_mesh(voxels_per_unit_cell=15)  # Higher res
        
        print(f"INFO: Calculating desorption paths (sampling method, n={n_samples})...")
        
        # Get coordinates of solid voxels
        solid_coords = np.argwhere(self._solid_voxels)
        
        if len(solid_coords) == 0:
            return {'error': 'No solid voxels found'}
        
        # Random sample
        if len(solid_coords) > n_samples:
            indices = np.random.choice(len(solid_coords), n_samples, replace=False)
            sampled_coords = solid_coords[indices]
        else:
            sampled_coords = solid_coords
            n_samples = len(solid_coords)
        
        # Convert to physical coords
        sampled_coords_mm = self._voxel_origin + sampled_coords * self._voxel_pitch
        
        # Calculate container parameters (centered at origin)
        if self.container_geometry == 'cylinder':
            center_x, center_y = 0.0, 0.0
            height = self.container_params['height']
            center_z = height / 2.0
            radius = self.container_params['diameter'] / 2
        elif self.container_geometry == 'sphere':
            center_x, center_y, center_z = 0.0, 0.0, 0.0
            radius = self.container_params['diameter'] / 2
        
        # For each sampled point, calculate path to edge
        distances = []
        
        for i, coord in enumerate(sampled_coords_mm):
            if i % 100 == 0:
                print(f"  Processing sample {i}/{n_samples}...")
            
            if self.container_geometry == 'cylinder':
                # Distance to cylindrical surface
                dx = coord[0] - center_x
                dy = coord[1] - center_y
                radial_dist = np.sqrt(dx**2 + dy**2)
                dist_to_wall = radius - radial_dist
                
                # Distance to caps
                dist_to_top = self.container_params['height'] - coord[2]
                dist_to_bottom = coord[2]
                
                min_dist = min(dist_to_wall, dist_to_top, dist_to_bottom)
            
            elif self.container_geometry == 'sphere':
                dx = coord[0] - center_x
                dy = coord[1] - center_y
                dz = coord[2] - center_z
                dist_from_center = np.sqrt(dx**2 + dy**2 + dz**2)
                min_dist = radius - dist_from_center
            
            distances.append(min_dist)
        
        distances = np.array(distances)
        
        # Apply tortuosity
        distances_corrected = distances * tau
        
        # Statistics
        mean_path = np.mean(distances_corrected)
        max_path = np.max(distances_corrected)
        min_path = np.min(distances_corrected)
        std_path = np.std(distances_corrected)
        median_path = np.median(distances_corrected)
        q25 = np.percentile(distances_corrected, 25)
        q75 = np.percentile(distances_corrected, 75)
        
        # Histogram
        hist, bins = np.histogram(distances_corrected, bins=50)
        
        results = {
            'mean_desorption_path_mm': mean_path,
            'median_desorption_path_mm': median_path,
            'max_desorption_path_mm': max_path,
            'min_desorption_path_mm': min_path,
            'std_desorption_path_mm': std_path,
            'q25_mm': q25,
            'q75_mm': q75,
            'tortuosity_applied': tau,
            'n_samples': n_samples,
            'method': 'statistical_sampling',
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
            selections: Dict of boolean flags:
                {
                    'basic_properties': True,
                    'specific_surface_area': True,
                    'wall_thickness_theoretical': True,
                    'channel_width_theoretical': True,
                    'desorption_path_theoretical': True,
                    'desorption_path_fast': False,
                    'desorption_path_sampling': False,
                    'desorption_sampling_n': 1000,
                    'external_surface': True,
                    'tortuosity': True,
                    'density_g_cm3': 1.05
                }
        
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
        
        # Tortuosity
        if selections.get('tortuosity', True):
            print("Calculating: Tortuosity...")
            results['tortuosity'] = self.calculator.calculate_tortuosity()
        
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
        
        if selections.get('desorption_path_theoretical', True):
            print("Calculating: Desorption Path (Theoretical)...")
            results['desorption_path_theoretical'] = \
                self.calculator.estimate_desorption_path_theoretical()
        
        # Voxel-based analysis (slower)
        if selections.get('desorption_path_fast', False):
            print("Calculating: Desorption Path (Fast)...")
            tau = results.get('tortuosity', {}).get('tortuosity', 1.41)
            results['desorption_path_fast'] = \
                self.calculator.calculate_desorption_path_fast(tau)
        
        if selections.get('desorption_path_sampling', False):
            print("Calculating: Desorption Path (Sampling)...")
            n_samples = selections.get('desorption_sampling_n', 1000)
            tau = results.get('tortuosity', {}).get('tortuosity', 1.41)
            results['desorption_path_sampling'] = \
                self.calculator.calculate_desorption_path_sampling(n_samples, tau)
        
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
    
    def calc_specific_surface_area(self) -> Dict:
        """Calculate SSA. Transforms keys: SSA_volumetric → ssa_volumetric."""
        bp = self._ensure_basic_props()
        raw = self._calc.calculate_specific_surface_area(bp)
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
    
    def calc_distance_to_edge_theoretical(self) -> Dict:
        """Transforms: mean_desorption_path_mm → mean_distance_mm."""
        raw = self._calc.estimate_desorption_path_theoretical()
        return {
            'mean_distance_mm': raw.get('mean_desorption_path_mm'),
            'max_distance_mm': raw.get('max_desorption_path_mm'),
            'method': raw.get('method', 'theoretical_geometric'),
            'accuracy': raw.get('accuracy', ''),
            'calculation_time': raw.get('calculation_time', 0),
        }
    
    def calc_tortuosity(self) -> Dict:
        """Adds 'method' key from 'note'."""
        raw = self._calc.calculate_tortuosity()
        return {
            'tortuosity': raw.get('tortuosity', 1.41),
            'effective_diffusivity_factor': raw.get('effective_diffusivity_factor'),
            'method': raw.get('note', 'Theoretical value for gyroid TPMS'),
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
        total_sa = bp.get('surface_area', 0)
        internal = max(0, total_sa - external)
        
        return {
            'external_area_mm2': external,
            'internal_area_mm2': internal,
            'total_area_mm2': total_sa,
            'external_fraction_percent': (external / total_sa * 100) if total_sa > 0 else 0,
            'calculation_time': raw.get('calculation_time', 0),
        }
    
    def calc_desorption_path_fast(self) -> Dict:
        """Transforms desorption path keys: mean_desorption_path_mm → mean_path_mm."""
        raw = self._calc.calculate_desorption_path_fast(tau=1.41)
        
        if 'error' in raw:
            return raw
        
        return {
            'mean_path_mm': raw.get('mean_desorption_path_mm', 0),
            'median_path_mm': raw.get('median_desorption_path_mm', 0),
            'std_path_mm': raw.get('std_desorption_path_mm', 0),
            'min_path_mm': max(0, raw.get('min_desorption_path_mm', 0)),
            'max_path_mm': raw.get('max_desorption_path_mm', 0),
            'samples': raw.get('n_samples', 0),
            'method': raw.get('method', 'fast_voxel'),
            'accuracy': 'Full voxel coverage',
            'histogram_counts': raw.get('histogram_counts'),
            'histogram_bins': raw.get('histogram_bins'),
            'calculation_time': raw.get('calculation_time', 0),
        }
    
    def calc_desorption_path_sampling(self, n_samples: int = 1000) -> Dict:
        """Transforms desorption path keys for sampling method."""
        raw = self._calc.calculate_desorption_path_sampling(
            n_samples=n_samples, tau=1.41
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
            'method': raw.get('method', 'monte_carlo_sampling'),
            'accuracy': f'Monte Carlo ({n_samples} samples)',
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
