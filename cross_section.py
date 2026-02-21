"""
Cross-Section Module — Wizualizacja i eksport przekrojów
==========================================================
Generuje przekroje 2D (XY, XZ, YZ) z mesha gyroidalnego.

Output:
  - Obraz binarny (solid/void)
  - Lokalna porowatość z przekroju
  - Eksport PNG/SVG

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

    def get_slice(self, plane: str = 'XY', position: float = None,
                  resolution: int = 15) -> Dict:
        """
        Get a 2D cross-section slice.

        Args:
            plane: 'XY' (at given Z), 'XZ' (at given Y), 'YZ' (at given X)
            position: Position along the normal axis [mm]. None = center.
            resolution: Voxels per unit cell

        Returns:
            Dict with:
              - 'image': 2D numpy array (True=solid, False=void)
              - 'porosity': local porosity of this slice
              - 'extent': [xmin, xmax, ymin, ymax] for plotting
              - 'position': actual position used
        """
        t0 = time.time()
        self._voxelize(resolution)

        grid = self._voxel_grid
        pitch = self._pitch
        origin = self._origin

        bounds = self.get_bounds()

        if plane == 'XY':
            axis = 2  # Z axis
            if position is None:
                position = (bounds['z_min'] + bounds['z_max']) / 2
            idx = int(round((position - origin[2]) / pitch))
            idx = np.clip(idx, 0, grid.shape[2] - 1)
            slice_2d = grid[:, :, idx].T  # Transpose for correct orientation
            extent = [bounds['x_min'], bounds['x_max'],
                      bounds['y_min'], bounds['y_max']]
            axis_labels = ('X [mm]', 'Y [mm]')
            title = f'XY plane at Z = {position:.3f} mm'

        elif plane == 'XZ':
            axis = 1  # Y axis
            if position is None:
                position = (bounds['y_min'] + bounds['y_max']) / 2
            idx = int(round((position - origin[1]) / pitch))
            idx = np.clip(idx, 0, grid.shape[1] - 1)
            slice_2d = grid[:, idx, :].T
            extent = [bounds['x_min'], bounds['x_max'],
                      bounds['z_min'], bounds['z_max']]
            axis_labels = ('X [mm]', 'Z [mm]')
            title = f'XZ plane at Y = {position:.3f} mm'

        elif plane == 'YZ':
            axis = 0  # X axis
            if position is None:
                position = (bounds['x_min'] + bounds['x_max']) / 2
            idx = int(round((position - origin[0]) / pitch))
            idx = np.clip(idx, 0, grid.shape[0] - 1)
            slice_2d = grid[idx, :, :].T
            extent = [bounds['y_min'], bounds['y_max'],
                      bounds['z_min'], bounds['z_max']]
            axis_labels = ('Y [mm]', 'Z [mm]')
            title = f'YZ plane at X = {position:.3f} mm'
        else:
            return {'error': f'Unknown plane: {plane}'}

        # Local porosity
        total_pixels = slice_2d.size
        solid_pixels = np.sum(slice_2d)
        local_porosity = 1.0 - (solid_pixels / total_pixels) if total_pixels > 0 else 0

        return {
            'image': slice_2d,
            'porosity': float(local_porosity),
            'solid_fraction': float(solid_pixels / total_pixels) if total_pixels > 0 else 0,
            'extent': extent,
            'position': float(position),
            'plane': plane,
            'axis_labels': axis_labels,
            'title': title,
            'resolution_px': slice_2d.shape,
            'pixel_size_um': pitch * 1000,
            'calculation_time': time.time() - t0,
        }

    def export_slice_image(self, slice_data: Dict, filepath: str,
                           dpi: int = 300, show_grid: bool = False):
        """
        Export cross-section as PNG image.

        Args:
            slice_data: Output from get_slice()
            filepath: Output path (.png or .svg)
            dpi: Resolution
            show_grid: Overlay unit cell grid
        """
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 1, figsize=(8, 8))

        # Plot binary image
        ax.imshow(
            slice_data['image'],
            extent=slice_data['extent'],
            cmap='gray_r',  # solid=black, void=white
            origin='lower',
            interpolation='nearest',
            aspect='equal'
        )

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
                     f"Porosity: {slice_data['porosity']:.1%}")

        fig.tight_layout()
        fig.savefig(filepath, dpi=dpi, bbox_inches='tight')
        plt.close(fig)

        return filepath
