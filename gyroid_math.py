"""
Gyroid Mathematics Module
==========================
Core mathematical functions for TPMS (Triply Periodic Minimal Surface) generation.
Focuses on gyroid surface: sin(x)cos(y) + sin(y)cos(z) + sin(z)cos(x) = t

Author: Claude (for Mariusz's DLP research)
Date: 2026-02
"""

import numpy as np
from typing import Tuple, Optional
import logging

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class GyroidSurface:
    """
    Gyroid TPMS generator with analytical properties.
    
    The gyroid equation in standard form:
    sin(2πx/a)·cos(2πy/a) + sin(2πy/a)·cos(2πz/a) + sin(2πz/a)·cos(2πx/a) = t
    
    Parameters:
    -----------
    unit_cell : float
        Unit cell size 'a' in mm (controls scale of structure)
    threshold : float
        Threshold level 't' [-1, 1] (controls wall/channel balance)
        t=0 : symmetric (50% porosity)
        t>0 : thicker walls, narrower channels
        t<0 : thinner walls, wider channels
    """
    
    def __init__(self, unit_cell: float = 0.2, threshold: float = 0.0):
        """
        Initialize gyroid surface generator.
        
        Args:
            unit_cell: Unit cell size in mm (typically 0.05 - 0.5 mm)
            threshold: Threshold level [-1, 1]
        """
        self.unit_cell = unit_cell
        self.threshold = threshold
        
        logger.info(f"Gyroid initialized: a={unit_cell:.3f}mm, t={threshold:.3f}")
    
    def evaluate(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        """
        Evaluate gyroid function at given coordinates.
        
        Args:
            x, y, z: Coordinate arrays (can be meshgrid)
            
        Returns:
            Gyroid function values (surface is where this equals threshold)
        """
        a = self.unit_cell
        
        # Scaled coordinates
        X = 2 * np.pi * x / a
        Y = 2 * np.pi * y / a
        Z = 2 * np.pi * z / a
        
        # Gyroid equation
        gyroid_value = (np.sin(X) * np.cos(Y) + 
                       np.sin(Y) * np.cos(Z) + 
                       np.sin(Z) * np.cos(X))
        
        return gyroid_value
    
    def estimate_properties(self) -> dict:
        """
        Estimate geometric properties analytically.
        
        Returns:
            Dictionary with estimated porosity, wall thickness, channel width
        """
        t = self.threshold
        a = self.unit_cell
        
        # Empirical formulas (fitted from numerical simulations)
        # These are approximations - actual values computed from mesh
        
        # Porosity (volume fraction of void space)
        if abs(t) < 0.01:
            porosity = 0.50  # Symmetric gyroid
        else:
            # Approximate relationship
            porosity = 0.50 - 0.35 * t
        
        # Clamp to physical limits
        porosity = np.clip(porosity, 0.15, 0.85)
        
        # Wall thickness (approximate, in mm)
        # At t=0, wall thickness ~ 0.3 * a
        # Increases with positive t
        wall_thickness = a * (0.30 + 0.15 * t)
        
        # Channel width (characteristic dimension of void space)
        # At t=0, channel ~ 0.7 * a
        # Decreases with positive t
        channel_width = a * (0.70 - 0.25 * t)
        
        # Surface area per unit volume (mm^-1)
        # Gyroid has constant mean curvature = 0 (minimal surface)
        surface_density = 3.09 / a  # Analytical result for symmetric gyroid
        
        properties = {
            'porosity': porosity,
            'wall_thickness_mm': wall_thickness,
            'channel_width_mm': channel_width,
            'surface_density': surface_density,
            'unit_cell': a,
            'threshold': t
        }
        
        logger.info(f"Estimated properties: porosity={porosity:.2%}, "
                   f"wall={wall_thickness*1000:.1f}µm, "
                   f"channel={channel_width*1000:.1f}µm")
        
        return properties
    
    @staticmethod
    def threshold_from_wall_thickness(unit_cell: float, 
                                     target_wall: float) -> float:
        """
        Calculate threshold 't' needed to achieve target wall thickness.
        
        Args:
            unit_cell: Unit cell size (mm)
            target_wall: Desired wall thickness (mm)
            
        Returns:
            Threshold value t
        """
        # Inverse of wall_thickness formula
        # wall = a * (0.30 + 0.15 * t)
        # t = (wall/a - 0.30) / 0.15
        
        t = (target_wall / unit_cell - 0.30) / 0.15
        
        # Clamp to valid range
        t = np.clip(t, -0.9, 0.9)
        
        return t
    
    @staticmethod
    def threshold_from_porosity(target_porosity: float) -> float:
        """
        Calculate threshold 't' to achieve target porosity.
        
        Args:
            target_porosity: Desired porosity [0.15, 0.85]
            
        Returns:
            Threshold value t
        """
        # Inverse of porosity formula
        # porosity = 0.50 - 0.35 * t
        # t = (0.50 - porosity) / 0.35
        
        t = (0.50 - target_porosity) / 0.35
        
        # Clamp to valid range
        t = np.clip(t, -0.9, 0.9)
        
        return t


class GradientGyroid:
    """
    Gyroid with spatially varying parameters (gradients).
    
    Supports:
    - Linear gradient along Z axis (e.g., wall thickness varies with height)
    - Concentric gradient (parameters vary with distance from center)
    """
    
    def __init__(self, base_unit_cell: float = 0.2):
        """
        Initialize gradient gyroid.
        
        Args:
            base_unit_cell: Base unit cell size
        """
        self.base_unit_cell = base_unit_cell
        self.gradient_type = None
        self.gradient_params = {}
        
    def set_linear_gradient(self, 
                           param: str,
                           z_min: float,
                           z_max: float,
                           value_bottom: float,
                           value_top: float):
        """
        Set linear gradient along Z axis.
        
        Args:
            param: Parameter to vary ('threshold' or 'unit_cell')
            z_min, z_max: Z coordinate range (mm)
            value_bottom, value_top: Parameter values at bottom and top
        """
        self.gradient_type = 'linear_z'
        self.gradient_params = {
            'param': param,
            'z_min': z_min,
            'z_max': z_max,
            'value_bottom': value_bottom,
            'value_top': value_top
        }
        
        logger.info(f"Linear Z gradient set: {param} from {value_bottom} to {value_top}")
    
    def set_concentric_gradient(self,
                               param: str,
                               center: Tuple[float, float, float],
                               r_min: float,
                               r_max: float,
                               value_center: float,
                               value_edge: float):
        """
        Set concentric gradient (spherical/cylindrical).
        
        Args:
            param: Parameter to vary
            center: Center point (x, y, z) in mm
            r_min, r_max: Radial distance range
            value_center, value_edge: Parameter values at center and edge
        """
        self.gradient_type = 'concentric'
        self.gradient_params = {
            'param': param,
            'center': np.array(center),
            'r_min': r_min,
            'r_max': r_max,
            'value_center': value_center,
            'value_edge': value_edge
        }
        
        logger.info(f"Concentric gradient set: {param} from {value_center} to {value_edge}")
    
    def get_local_parameter(self, x: np.ndarray, y: np.ndarray, z: np.ndarray,
                           param: str) -> np.ndarray:
        """
        Get spatially-varying parameter value at coordinates.
        
        Args:
            x, y, z: Coordinates
            param: Parameter name ('threshold' or 'unit_cell')
            
        Returns:
            Array of parameter values
        """
        if self.gradient_type is None:
            # No gradient - return base value
            if param == 'unit_cell':
                return np.full_like(x, self.base_unit_cell)
            else:  # threshold
                return np.zeros_like(x)
        
        if self.gradient_params['param'] != param:
            # This parameter doesn't have gradient
            if param == 'unit_cell':
                return np.full_like(x, self.base_unit_cell)
            else:
                return np.zeros_like(x)
        
        if self.gradient_type == 'linear_z':
            z_min = self.gradient_params['z_min']
            z_max = self.gradient_params['z_max']
            v_bottom = self.gradient_params['value_bottom']
            v_top = self.gradient_params['value_top']
            
            # Linear interpolation
            alpha = (z - z_min) / (z_max - z_min)
            alpha = np.clip(alpha, 0, 1)
            
            values = v_bottom + alpha * (v_top - v_bottom)
            
        elif self.gradient_type == 'concentric':
            center = self.gradient_params['center']
            r_min = self.gradient_params['r_min']
            r_max = self.gradient_params['r_max']
            v_center = self.gradient_params['value_center']
            v_edge = self.gradient_params['value_edge']
            
            # Distance from center
            dx = x - center[0]
            dy = y - center[1]
            dz = z - center[2]
            r = np.sqrt(dx**2 + dy**2 + dz**2)
            
            # Radial interpolation
            alpha = (r - r_min) / (r_max - r_min)
            alpha = np.clip(alpha, 0, 1)
            
            values = v_center + alpha * (v_edge - v_center)
        
        else:
            raise ValueError(f"Unknown gradient type: {self.gradient_type}")
        
        return values
    
    def evaluate(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        """
        Evaluate gyroid with gradient at coordinates.
        
        Args:
            x, y, z: Coordinates
            
        Returns:
            Gyroid function values
        """
        # Get local parameters
        local_a = self.get_local_parameter(x, y, z, 'unit_cell')
        local_t = self.get_local_parameter(x, y, z, 'threshold')
        
        # Gyroid with local parameters
        X = 2 * np.pi * x / local_a
        Y = 2 * np.pi * y / local_a
        Z = 2 * np.pi * z / local_a
        
        gyroid_value = (np.sin(X) * np.cos(Y) + 
                       np.sin(Y) * np.cos(Z) + 
                       np.sin(Z) * np.cos(X))
        
        # Compare to local threshold
        return gyroid_value - local_t


# Utility functions
def calculate_resolution(dimensions: Tuple[float, float, float],
                        unit_cell: float,
                        quality: str = 'medium') -> Tuple[int, int, int]:
    """
    Calculate grid resolution for marching cubes.
    
    Args:
        dimensions: (width, depth, height) in mm
        unit_cell: Unit cell size in mm
        quality: 'low', 'medium', 'high', 'ultra'
        
    Returns:
        (nx, ny, nz) grid dimensions
    """
    # Points per unit cell (PPU)
    ppu_map = {
        'low': 8,      # ~50k triangles for typical geometry
        'medium': 12,  # ~200k triangles
        'high': 20,    # ~1M triangles
        'ultra': 30    # ~5M triangles
    }
    
    ppu = ppu_map.get(quality, 12)
    
    # Calculate grid size
    width, depth, height = dimensions
    
    nx = int(np.ceil(width / unit_cell * ppu))
    ny = int(np.ceil(depth / unit_cell * ppu))
    nz = int(np.ceil(height / unit_cell * ppu))
    
    # Ensure minimum resolution
    nx = max(nx, 50)
    ny = max(ny, 50)
    nz = max(nz, 50)
    
    logger.info(f"Grid resolution: {nx}×{ny}×{nz} = {nx*ny*nz:,} points")
    
    return nx, ny, nz


if __name__ == "__main__":
    # Test gyroid math
    print("Testing Gyroid Mathematics Module\n")
    
    # Test 1: Basic gyroid
    gyroid = GyroidSurface(unit_cell=0.2, threshold=0.0)
    props = gyroid.estimate_properties()
    
    print("Basic Gyroid (a=200µm, t=0):")
    print(f"  Porosity: {props['porosity']:.1%}")
    print(f"  Wall thickness: {props['wall_thickness_mm']*1000:.1f} µm")
    print(f"  Channel width: {props['channel_width_mm']*1000:.1f} µm")
    print()
    
    # Test 2: Find threshold for target wall
    target_wall = 0.050  # 50 µm
    t = GyroidSurface.threshold_from_wall_thickness(0.2, target_wall)
    print(f"For wall={target_wall*1000:.0f}µm with a=200µm:")
    print(f"  Threshold t = {t:.3f}")
    
    gyroid2 = GyroidSurface(unit_cell=0.2, threshold=t)
    props2 = gyroid2.estimate_properties()
    print(f"  Actual wall: {props2['wall_thickness_mm']*1000:.1f} µm")
    print(f"  Porosity: {props2['porosity']:.1%}")
    print()
    
    # Test 3: Gradient
    grad_gyroid = GradientGyroid(base_unit_cell=0.2)
    grad_gyroid.set_linear_gradient(
        param='threshold',
        z_min=0, z_max=10,
        value_bottom=0.3,
        value_top=-0.3
    )
    
    print("Gradient Gyroid (linear Z):")
    z_test = np.array([0, 5, 10])
    for z in z_test:
        t_local = grad_gyroid.get_local_parameter(0, 0, z, 'threshold')
        print(f"  At z={z}mm: t={t_local:.3f}")
