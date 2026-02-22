"""
Cross-Section Module — Visualization and Export
==========================================================
Generates 2D cross-sections (XY, XZ, YZ) from the gyroid mesh.

Output:
  - Vector/Binary image (solid/void)
  - Local porosity from slice area
  - PNG/SVG export

Author: Claude (Anthropic)
Date: 2026-02-20
"""

import numpy as np
from typing import Dict, Optional, Tuple
import trimesh
import time
import logging

logger = logging.getLogger(__name__)


class CrossSectionAnalyzer:
    """Generate and analyze cross-sections of gyroid mesh."""

    def __init__(self, mesh: trimesh.Trimesh, unit_cell_mm: float):
        self.mesh = mesh
        self.unit_cell_mm = unit_cell_mm
        self._voxel_grid = None
        self._pitch = None
        self._origin = None

    def _voxelize(self, voxels_per_uc: int = 15):
        """Voxelize mesh for cross-section generation."""
        if self._voxel_grid is not None:
            return

        pitch = self.unit_cell_mm / voxels_per_uc
        voxel = self.mesh.voxelized(pitch=pitch)
        self._voxel_grid = voxel.matrix
        self._pitch = pitch
        self._origin = voxel.transform[:3, 3]

    def get_bounds(self) -> Dict:
        """Get mesh bounds for slider ranges."""
        bounds = self.mesh.bounds  # [[xmin,ymin,zmin], [xmax,ymax,zmax]]
        return {
            'x_min': float(bounds[0][0]),
            'x_max': float(bounds[1][0]),
            'y_min': float(bounds[0][1]),
            'y_max': float(bounds[1][1]),
            'z_min': float(bounds[0][2]),
            'z_max': float(bounds[1][2]),
        }

    def get_slice(self, plane: str = 'XY', position: float = None) -> Dict:
        """
        Get a 2D cross-section slice using mesh slicing (accurate).

        Args:
            plane: 'XY' (at given Z), 'XZ' (at given Y), 'YZ' (at given X)
            position: Position along the normal axis [mm]. None = center.

        Returns:
            Dict with:
              - 'path': trimesh.path.Path2D object
              - 'porosity': local porosity of this slice
              - 'extent': [xmin, xmax, ymin, ymax] for plotting
              - 'position': actual position used
        """
        t0 = time.time()
        bounds = self.get_bounds()

        if plane == 'XY':
            normal = [0, 0, 1]
            if position is None: position = (bounds['z_min'] + bounds['z_max']) / 2
            origin = [0, 0, position]
            extent = [bounds['x_min'], bounds['x_max'], bounds['y_min'], bounds['y_max']]
            axis_labels = ('X [mm]', 'Y [mm]')
        elif plane == 'XZ':
            normal = [0, 1, 0]
            if position is None: position = (bounds['y_min'] + bounds['y_max']) / 2
            origin = [0, position, 0]
            extent = [bounds['x_min'], bounds['x_max'], bounds['z_min'], bounds['z_max']]
            axis_labels = ('X [mm]', 'Z [mm]')
        elif plane == 'YZ':
            normal = [1, 0, 0]
            if position is None: position = (bounds['x_min'] + bounds['x_max']) / 2
            origin = [position, 0, 0]
            extent = [bounds['y_min'], bounds['y_max'], bounds['z_min'], bounds['z_max']]
            axis_labels = ('Y [mm]', 'Z [mm]')
        else:
            return {'error': f'Unknown plane: {plane}'}

        # Perform slicing
        section = self.mesh.section(plane_origin=origin, plane_normal=normal)

        if section is None:
            return {'error': 'No intersection found'}

        # Convert to 2D
        path_2d, _ = section.to_planar()

        # Calculate accurate area
        solid_area = path_2d.area

        # Determine container cross-section area
        container_area = 1.0
        if self.mesh.metadata.get('geometry_type') == 'cylinder' or 'diameter' in self.mesh.metadata:
            # Fallback to metadata or assume circular
            r = self.mesh.bounds[1, 0] # assuming centered
            container_area = np.pi * r**2
        else:
            # Rectangular bounding box of the slice
            container_area = (extent[1] - extent[0]) * (extent[3] - extent[2])

        local_porosity = 1.0 - (solid_area / container_area) if container_area > 0 else 0

        return {
            'path': path_2d,
            'porosity': float(local_porosity),
            'solid_area': float(solid_area),
            'extent': extent,
            'position': float(position),
            'plane': plane,
            'axis_labels': axis_labels,
            'title': f'{plane} plane at {position:.3f} mm',
            'calculation_time': time.time() - t0,
        }

    def export_slice_image(self, slice_data: Dict, filepath: str,
                           dpi: int = 300, show_grid: bool = False):
        """
        Export cross-section as PNG or SVG.

        Args:
            slice_data: Output from get_slice()
            filepath: Output path (.png or .svg)
            dpi: Resolution
            show_grid: Overlay unit cell grid
        """
        if filepath.lower().endswith('.svg'):
            # Direct vector export from trimesh path
            slice_data['path'].export(filepath)
            return filepath

        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(8, 8))

        # Plot vector path
        for entity in slice_data['path'].entities:
            disc = entity.discrete(slice_data['path'].vertices)
            ax.plot(disc[:, 0], disc[:, 1], 'k-', linewidth=1)
            # Fill solid area
            ax.fill(disc[:, 0], disc[:, 1], 'k', alpha=0.8)

        # Grid overlay
        if show_grid:
            bounds = slice_data['extent']
            uc = self.unit_cell_mm
            for x in np.arange(bounds[0], bounds[1], uc):
                ax.axvline(x, color='blue', alpha=0.2, linewidth=0.5)
            for y in np.arange(bounds[2], bounds[3], uc):
                ax.axhline(y, color='blue', alpha=0.2, linewidth=0.5)

        ax.set_xlabel(slice_data['axis_labels'][0])
        ax.set_ylabel(slice_data['axis_labels'][1])
        ax.set_title(f"{slice_data['title']}\n"
                     f"Local Porosity: {slice_data['porosity']:.1%}")
        ax.set_aspect('equal')
        ax.set_xlim(slice_data['extent'][0], slice_data['extent'][1])
        ax.set_ylim(slice_data['extent'][2], slice_data['extent'][3])

        fig.tight_layout()
        fig.savefig(filepath, dpi=dpi, bbox_inches='tight')
        plt.close(fig)

        return filepath
