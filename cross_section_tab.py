"""
Cross-Section Tab — Wizualizacja przekrojów 2D
=================================================

Author: Claude (Anthropic)
Date: 2026-02-20
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QPushButton, QComboBox, QSlider, QFileDialog, QCheckBox, QSpinBox
)
from PyQt6.QtCore import Qt
import numpy as np
import logging

logger = logging.getLogger(__name__)

import matplotlib
matplotlib.use('QtAgg')
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from cross_section import CrossSectionAnalyzer


class CrossSectionTab(QWidget):
    """Tab for cross-section visualization and export."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._analyzer = None
        self._current_slice = None
        self._init_ui()

    def _init_ui(self):
        layout = QHBoxLayout(self)

        # Left: Controls
        controls = QVBoxLayout()
        controls.setSpacing(10)

        # Plane selection
        plane_group = QGroupBox("Cross-Section Plane")
        pg_layout = QVBoxLayout(plane_group)
        self.plane_combo = QComboBox()
        self.plane_combo.addItems(['XY (at Z)', 'XZ (at Y)', 'YZ (at X)'])
        self.plane_combo.currentIndexChanged.connect(self._on_plane_changed)
        pg_layout.addWidget(self.plane_combo)
        controls.addWidget(plane_group)

        # Position slider
        pos_group = QGroupBox("Position")
        posl = QVBoxLayout(pos_group)

        self.pos_slider = QSlider(Qt.Orientation.Horizontal)
        self.pos_slider.setRange(0, 100)
        self.pos_slider.setValue(50)
        self.pos_slider.valueChanged.connect(self._on_slider_changed)
        posl.addWidget(self.pos_slider)

        self.pos_label = QLabel("Position: center")
        posl.addWidget(self.pos_label)

        controls.addWidget(pos_group)

        # Options
        opt_group = QGroupBox("Display Options")
        ol = QVBoxLayout(opt_group)

        self.cb_grid = QCheckBox("Show unit cell grid")
        ol.addWidget(self.cb_grid)

        h = QHBoxLayout()
        h.addWidget(QLabel("Resolution:"))
        self.res_spin = QSpinBox()
        self.res_spin.setRange(8, 30)
        self.res_spin.setValue(15)
        h.addWidget(self.res_spin)
        h.addWidget(QLabel("vox/UC"))
        ol.addLayout(h)

        controls.addWidget(opt_group)

        # Buttons
        self.btn_update = QPushButton("▶ Update View")
        self.btn_update.clicked.connect(self._update_view)
        controls.addWidget(self.btn_update)

        self.btn_export_png = QPushButton("Export PNG...")
        self.btn_export_png.clicked.connect(lambda: self._export('png'))
        controls.addWidget(self.btn_export_png)

        self.btn_export_svg = QPushButton("Export SVG...")
        self.btn_export_svg.clicked.connect(lambda: self._export('svg'))
        controls.addWidget(self.btn_export_svg)

        # Info
        self.info_label = QLabel("")
        self.info_label.setWordWrap(True)
        controls.addWidget(self.info_label)

        controls.addStretch()

        controls_widget = QWidget()
        controls_widget.setLayout(controls)
        controls_widget.setMaximumWidth(250)
        layout.addWidget(controls_widget)

        # Right: Figure
        self.figure = Figure(figsize=(8, 8))
        self.canvas = FigureCanvas(self.figure)
        layout.addWidget(self.canvas, stretch=1)

    def set_data(self, mesh, unit_cell_mm):
        """Set mesh for cross-section analysis."""
        self._analyzer = CrossSectionAnalyzer(mesh, unit_cell_mm)
        self._update_slider_range()
        self._update_view()

    def _get_plane(self) -> str:
        idx = self.plane_combo.currentIndex()
        return ['XY', 'XZ', 'YZ'][idx]

    def _update_slider_range(self):
        """Update slider range based on current plane."""
        if self._analyzer is None:
            return

        bounds = self._analyzer.get_bounds()
        plane = self._get_plane()

        if plane == 'XY':
            self._slider_min = bounds['z_min']
            self._slider_max = bounds['z_max']
        elif plane == 'XZ':
            self._slider_min = bounds['y_min']
            self._slider_max = bounds['y_max']
        else:
            self._slider_min = bounds['x_min']
            self._slider_max = bounds['x_max']

    def _slider_to_position(self) -> float:
        """Convert slider value (0-100) to physical position."""
        frac = self.pos_slider.value() / 100.0
        return self._slider_min + frac * (self._slider_max - self._slider_min)

    def _on_plane_changed(self):
        self._update_slider_range()
        self.pos_slider.setValue(50)

    def _on_slider_changed(self):
        pos = self._slider_to_position()
        self.pos_label.setText(f"Position: {pos:.3f} mm")

    def _update_view(self):
        """Generate and display cross-section."""
        if self._analyzer is None:
            self.info_label.setText("⚠ Generate a mesh first!")
            return

        plane = self._get_plane()
        position = self._slider_to_position()

        slice_data = self._analyzer.get_slice(
            plane=plane,
            position=position,
            resolution=self.res_spin.value()
        )

        if 'error' in slice_data:
            self.info_label.setText(f"⚠ {slice_data['error']}")
            return

        self._current_slice = slice_data

        # Display
        self.figure.clear()
        ax = self.figure.add_subplot(111)

        ax.imshow(
            slice_data['image'],
            extent=slice_data['extent'],
            cmap='gray_r',
            origin='lower',
            interpolation='nearest',
            aspect='equal'
        )

        # Grid overlay
        if self.cb_grid.isChecked() and self._analyzer is not None:
            uc = self._analyzer.unit_cell_mm
            bounds = slice_data['extent']
            for x in np.arange(bounds[0], bounds[1], uc):
                ax.axvline(x, color='cyan', alpha=0.3, linewidth=0.5)
            for y in np.arange(bounds[2], bounds[3], uc):
                ax.axhline(y, color='cyan', alpha=0.3, linewidth=0.5)

        ax.set_xlabel(slice_data['axis_labels'][0])
        ax.set_ylabel(slice_data['axis_labels'][1])
        ax.set_title(slice_data['title'])

        self.figure.tight_layout()
        self.canvas.draw()

        # Info
        self.info_label.setText(
            f"Porosity: {slice_data['porosity']:.1%}\n"
            f"Solid: {slice_data['solid_fraction']:.1%}\n"
            f"Resolution: {slice_data['resolution_px'][0]}×{slice_data['resolution_px'][1]} px\n"
            f"Pixel: {slice_data['pixel_size_um']:.1f} µm\n"
            f"Time: {slice_data['calculation_time']:.2f}s"
        )

    def _export(self, fmt):
        """Export current slice."""
        if self._current_slice is None or self._analyzer is None:
            return

        ext = 'PNG files (*.png)' if fmt == 'png' else 'SVG files (*.svg)'
        path, _ = QFileDialog.getSaveFileName(
            self, f"Export Cross-Section",
            f"cross_section.{fmt}", ext
        )
        if not path:
            return

        self._analyzer.export_slice_image(
            self._current_slice, path,
            dpi=300,
            show_grid=self.cb_grid.isChecked()
        )
        self.info_label.setText(f"Exported to {path}")
