"""
Mesh Generator Module
=====================
Converts implicit gyroid surface to triangle mesh using marching cubes.
Handles different base geometries (cylinder, sphere, cone, etc.)

Author: Claude
Date: 2026-02
"""

import numpy as np
from typing import Tuple, Optional, List
import logging
from skimage import measure
import trimesh

from gyroid_math import GyroidSurface, GradientGyroid, calculate_resolution

logger = logging.getLogger(__name__)


class BaseGeometry:
    """Base class for container geometries."""

    def contains(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        """
        Check if points are inside the geometry.

        Args:
            x, y, z: Coordinate arrays

        Returns:
            Boolean array: True where points are inside
        """
        raise NotImplementedError

    def get_bounds(self) -> Tuple[float, float, float]:
        """
        Get bounding box dimensions.

        Returns:
            (width, depth, height) in mm
        """
        raise NotImplementedError

    def get_center(self) -> Tuple[float, float, float]:
        """Get center point (x, y, z)."""
        raise NotImplementedError


class CylinderGeometry(BaseGeometry):
    """Cylindrical container."""

    def __init__(self, diameter: float, height: float,
                 center: Optional[Tuple[float, float, float]] = None):
        """
        Args:
            diameter: Cylinder diameter (mm)
            height: Cylinder height (mm)
            center: Center point (x, y, z). If None, centered at origin.
        """
        self.diameter = diameter
        self.radius = diameter / 2
        self.height = height

        if center is None:
            # Center at origin, extend symmetrically
            self.center = (0, 0, height / 2)
        else:
            self.center = center

        logger.info(f"Cylinder: D={diameter:.2f}mm, H={height:.2f}mm, "
                   f"center={self.center}")

    def contains(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        """Check if points inside cylinder."""
        cx, cy, cz = self.center

        # Radial distance from cylinder axis
        r = np.sqrt((x - cx)**2 + (y - cy)**2)

        # Height check
        z_min = cz - self.height / 2
        z_max = cz + self.height / 2

        inside = (r <= self.radius) & (z >= z_min) & (z <= z_max)

        return inside

    def get_bounds(self) -> Tuple[float, float, float]:
        """Return (width, depth, height)."""
        return (self.diameter, self.diameter, self.height)

    def get_center(self) -> Tuple[float, float, float]:
        return self.center


class SphereGeometry(BaseGeometry):
    """Spherical container."""

    def __init__(self, diameter: float,
                 center: Optional[Tuple[float, float, float]] = None):
        """
        Args:
            diameter: Sphere diameter (mm)
            center: Center point. If None, at origin.
        """
        self.diameter = diameter
        self.radius = diameter / 2
        self.center = center if center is not None else (0, 0, 0)

        logger.info(f"Sphere: D={diameter:.2f}mm, center={self.center}")

    def contains(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        cx, cy, cz = self.center
        r = np.sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2)
        return r <= self.radius

    def get_bounds(self) -> Tuple[float, float, float]:
        return (self.diameter, self.diameter, self.diameter)

    def get_center(self) -> Tuple[float, float, float]:
        return self.center


class SpindleGeometry(BaseGeometry):
    """
    Spindle (fusiform) shape - thicker in middle, tapered at ends.
    Like an ellipsoid of revolution.
    """

    def __init__(self, max_diameter: float, length: float,
                 center: Optional[Tuple[float, float, float]] = None):
        """
        Args:
            max_diameter: Maximum diameter at center (mm)
            length: Total length along Z axis (mm)
            center: Center point
        """
        self.max_diameter = max_diameter
        self.max_radius = max_diameter / 2
        self.length = length
        self.center = center if center is not None else (0, 0, length / 2)

        logger.info(f"Spindle: D_max={max_diameter:.2f}mm, L={length:.2f}mm")

    def contains(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        cx, cy, cz = self.center

        # Ellipsoid equation: (x/a)² + (y/a)² + (z/c)² <= 1
        a = self.max_radius  # Semi-axis in x, y
        c = self.length / 2  # Semi-axis in z

        term = ((x - cx) / a)**2 + ((y - cy) / a)**2 + ((z - cz) / c)**2

        return term <= 1.0

    def get_bounds(self) -> Tuple[float, float, float]:
        return (self.max_diameter, self.max_diameter, self.length)

    def get_center(self) -> Tuple[float, float, float]:
        return self.center


class ConeGeometry(BaseGeometry):
    """Conical container."""

    def __init__(self, base_diameter: float, height: float,
                 apex_at_top: bool = True,
                 center: Optional[Tuple[float, float, float]] = None):
        """
        Args:
            base_diameter: Diameter at base (mm)
            height: Cone height (mm)
            apex_at_top: If True, apex at top; if False, apex at bottom
            center: Center of base
        """
        self.base_diameter = base_diameter
        self.base_radius = base_diameter / 2
        self.height = height
        self.apex_at_top = apex_at_top

        if center is None:
            self.center = (0, 0, 0 if apex_at_top else height)
        else:
            self.center = center

        logger.info(f"Cone: D_base={base_diameter:.2f}mm, H={height:.2f}mm")

    def contains(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        cx, cy, cz = self.center

        if self.apex_at_top:
            # Base at z=cz, apex at z=cz+height
            z_base = cz
            z_apex = cz + self.height

            # Linear radius variation
            z_normalized = (z - z_base) / self.height
            z_normalized = np.clip(z_normalized, 0, 1)

            local_radius = self.base_radius * (1 - z_normalized)
        else:
            # Apex at z=cz-height, base at z=cz
            z_apex = cz - self.height
            z_base = cz

            z_normalized = (z - z_apex) / self.height
            z_normalized = np.clip(z_normalized, 0, 1)

            local_radius = self.base_radius * z_normalized

        # Radial check
        r = np.sqrt((x - cx)**2 + (y - cy)**2)

        # Height check
        if self.apex_at_top:
            in_height = (z >= z_base) & (z <= z_apex)
        else:
            in_height = (z >= z_apex) & (z <= z_base)

        inside = (r <= local_radius) & in_height

        return inside

    def get_bounds(self) -> Tuple[float, float, float]:
        return (self.base_diameter, self.base_diameter, self.height)

    def get_center(self) -> Tuple[float, float, float]:
        return self.center


class CubeGeometry(BaseGeometry):
    """Cubic container."""

    def __init__(self, size: float,
                 center: Optional[Tuple[float, float, float]] = None):
        """
        Args:
            size: Cube side length (mm)
            center: Center point
        """
        self.size = size
        self.center = center if center is not None else (0, 0, size / 2)

        logger.info(f"Cube: size={size:.2f}mm")

    def contains(self, x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        cx, cy, cz = self.center
        half = self.size / 2

        inside = ((np.abs(x - cx) <= half) &
                 (np.abs(y - cy) <= half) &
                 (np.abs(z - cz) <= half))

        return inside

    def get_bounds(self) -> Tuple[float, float, float]:
        return (self.size, self.size, self.size)

    def get_center(self) -> Tuple[float, float, float]:
        return self.center


class GyroidMeshGenerator:
    """
    Main mesh generator combining gyroid surface with container geometry.
    """

    def __init__(self, geometry: BaseGeometry,
                 gyroid: GyroidSurface,
                 quality: str = 'medium'):
        """
        Args:
            geometry: Container geometry (cylinder, sphere, etc.)
            gyroid: Gyroid surface generator
            quality: Mesh quality ('low', 'medium', 'high', 'ultra')
        """
        self.geometry = geometry
        self.gyroid = gyroid
        self.quality = quality

        self.mesh = None
        self.generation_time = None

    def generate(self) -> trimesh.Trimesh:
        """
        Generate gyroid mesh using marching cubes.

        Returns:
            Trimesh object
        """
        import time
        start_time = time.time()

        logger.info("Starting mesh generation...")

        # Get bounds
        width, depth, height = self.geometry.get_bounds()
        cx, cy, cz = self.geometry.get_center()

        # Calculate grid resolution
        nx, ny, nz = calculate_resolution(
            (width, depth, height),
            self.gyroid.unit_cell,
            self.quality
        )

        # Create coordinate grid
        x = np.linspace(cx - width/2, cx + width/2, nx)
        y = np.linspace(cy - depth/2, cy + depth/2, ny)
        z = np.linspace(cz - height/2, cz + height/2, nz)

        X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

        logger.info("Evaluating gyroid function...")

        # Evaluate gyroid
        gyroid_values = self.gyroid.evaluate(X, Y, Z)

        # Apply container geometry mask
        logger.info("Applying geometry mask...")
        inside_mask = self.geometry.contains(X, Y, Z)

        # Set outside values to large number (will be excluded by marching cubes)
        scalar_field = np.where(inside_mask, gyroid_values, 10.0)

        # Extract isosurface using marching cubes
        logger.info("Running marching cubes...")

        try:
            verts, faces, normals, values = measure.marching_cubes(
                scalar_field,
                level=self.gyroid.threshold,
                spacing=(
                    width / (nx - 1),
                    depth / (ny - 1),
                    height / (nz - 1)
                ),
                allow_degenerate=False
            )

            # Translate vertices to correct position
            verts[:, 0] += cx - width / 2
            verts[:, 1] += cy - depth / 2
            verts[:, 2] += cz - height / 2

            # Create trimesh object
            self.mesh = trimesh.Trimesh(vertices=verts, faces=faces,
                                       vertex_normals=normals)

            # Fix mesh issues



            # Try to make watertight (best effort)
            if not self.mesh.is_watertight:
                logger.warning("Mesh is not watertight - attempting to fix...")
                trimesh.repair.fill_holes(self.mesh)
                trimesh.repair.fix_normals(self.mesh)

            self.generation_time = time.time() - start_time

            logger.info(f"Mesh generated successfully in {self.generation_time:.1f}s")
            logger.info(f"  Vertices: {len(self.mesh.vertices):,}")
            logger.info(f"  Faces: {len(self.mesh.faces):,}")
            logger.info(f"  Watertight: {self.mesh.is_watertight}")

            return self.mesh

        except Exception as e:
            logger.error(f"Marching cubes failed: {str(e)}")
            raise

    def export_stl(self, filename: str):
        """
        Export mesh to STL file.

        Args:
            filename: Output STL path
        """
        if self.mesh is None:
            raise ValueError("No mesh generated yet. Call generate() first.")

        self.mesh.export(filename)

        # Get file size
        import os
        file_size = os.path.getsize(filename) / (1024 * 1024)  # MB

        logger.info(f"STL exported: {filename} ({file_size:.1f} MB)")

    def get_statistics(self) -> dict:
        """
        Calculate mesh statistics.

        Returns:
            Dictionary with mesh properties
        """
        if self.mesh is None:
            raise ValueError("No mesh generated yet.")

        stats = {
            'n_vertices': len(self.mesh.vertices),
            'n_faces': len(self.mesh.faces),
            'is_watertight': self.mesh.is_watertight,
            'volume_mm3': self.mesh.volume if self.mesh.is_watertight else None,
            'surface_area_mm2': self.mesh.area,
            'generation_time_s': self.generation_time
        }

        # Calculate container volume
        bounds = self.geometry.get_bounds()

        if isinstance(self.geometry, CylinderGeometry):
            container_volume = np.pi * (self.geometry.radius ** 2) * self.geometry.height
        elif isinstance(self.geometry, SphereGeometry):
            container_volume = (4/3) * np.pi * (self.geometry.radius ** 3)
        elif isinstance(self.geometry, CubeGeometry):
            container_volume = self.geometry.size ** 3
        elif isinstance(self.geometry, SpindleGeometry):
            # Ellipsoid volume
            a = self.geometry.max_radius
            c = self.geometry.length / 2
            container_volume = (4/3) * np.pi * a * a * c
        elif isinstance(self.geometry, ConeGeometry):
            container_volume = (1/3) * np.pi * (self.geometry.base_radius ** 2) * self.geometry.height
        else:
            container_volume = None

        stats['container_volume_mm3'] = container_volume

        if self.mesh.is_watertight and container_volume is not None:
            stats['porosity'] = 1 - (self.mesh.volume / container_volume)
            stats['gyroid_volume_mm3'] = self.mesh.volume
            stats['void_volume_mm3'] = container_volume - self.mesh.volume
        else:
            stats['porosity'] = None
            stats['gyroid_volume_mm3'] = None
            stats['void_volume_mm3'] = None

        return stats


# Test code
if __name__ == "__main__":
    print("Testing Mesh Generator\n")

    # Test: Generate simple cylinder gyroid
    geometry = CylinderGeometry(diameter=4.6, height=10.0)
    gyroid = GyroidSurface(unit_cell=0.2, threshold=0.0)

    generator = GyroidMeshGenerator(geometry, gyroid, quality='low')

    print("Generating mesh...")
    mesh = generator.generate()

    stats = generator.get_statistics()
    print("\nMesh Statistics:")
    for key, value in stats.items():
        if value is not None:
            if 'volume' in key or 'area' in key:
                print(f"  {key}: {value:.2f}")
            elif 'porosity' in key:
                print(f"  {key}: {value:.1%}")
            elif isinstance(value, (int, np.integer)):
                print(f"  {key}: {value:,}")
            else:
                print(f"  {key}: {value}")

    # Export
    # generator.export_stl("test_cylinder.stl")
    print("\nTest completed successfully!")
