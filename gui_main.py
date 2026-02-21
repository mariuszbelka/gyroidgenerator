"""
GUI Main Window — Gyroid Generator v3.0
=========================================
DLP Sorbent Design Platform

v3.0 adds:
  - Built-in 3D viewer (no pyglet dependency)
  - Predictions tab (ΔP, Van Deemter, desorption, capacity)
  - Distributions tab (wall/channel histograms, connectivity, ASA, throat, printability)
  - Cross-section tab (2D slices with export)
  - Compare tab (side-by-side design comparison)
  - "Save to Compare" workflow

Author: Claude (Anthropic)
Date: 2026-02-20
"""

import sys
import os
from pathlib import Path
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QLabel, QPushButton, QSlider, QComboBox, QGroupBox, QRadioButton,
    QButtonGroup, QSpinBox, QDoubleSpinBox, QTextEdit, QProgressBar,
    QFileDialog, QMessageBox, QTabWidget, QCheckBox, QScrollArea
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt6.QtGui import QFont, QIcon
import logging

# Import gyroid modules
from gyroid_math import GyroidSurface
from mesh_generator import (
    GyroidMeshGenerator, CylinderGeometry, SphereGeometry,
    SpindleGeometry, ConeGeometry, CubeGeometry
)
from statistics_analyzer_v2 import StatisticsV2, GeometryParams, GyroidParams
from statistics_tab_v2 import StatisticsTabV2

# v3 modules
from viewer_3d import MeshViewer3D
from predictions_tab import PredictionsTab
from distributions_tab import DistributionsTab
from cross_section_tab import CrossSectionTab
from compare_tab import CompareTab

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


class ParameterGroup(QGroupBox):
    """Group box for related parameters with sliders."""
    
    def __init__(self, title: str, parent=None):
        super().__init__(title, parent)
        self.layout = QVBoxLayout()
        self.setLayout(self.layout)
        self.parameters = {}
    
    def add_slider(self, name: str, label: str, 
                   min_val: float, max_val: float, default: float,
                   step: float = 0.01, suffix: str = '', decimals: int = 2):
        """Add a parameter with slider and spinbox."""
        container = QWidget()
        hlayout = QHBoxLayout()
        container.setLayout(hlayout)
        
        lbl = QLabel(label)
        lbl.setMinimumWidth(100)
        hlayout.addWidget(lbl)
        
        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setMinimum(0)
        slider.setMaximum(1000)
        default_pos = int((default - min_val) / (max_val - min_val) * 1000)
        slider.setValue(default_pos)
        hlayout.addWidget(slider, stretch=1)
        
        spinbox = QDoubleSpinBox()
        spinbox.setMinimum(min_val)
        spinbox.setMaximum(max_val)
        spinbox.setValue(default)
        spinbox.setSingleStep(step)
        spinbox.setDecimals(decimals)
        if suffix:
            spinbox.setSuffix(f' {suffix}')
        hlayout.addWidget(spinbox)
        
        # Connect slider <-> spinbox
        def slider_to_spin(val):
            value = min_val + (val / 1000.0) * (max_val - min_val)
            spinbox.setValue(value)
        
        def spin_to_slider(val):
            pos = int((val - min_val) / (max_val - min_val) * 1000)
            slider.blockSignals(True)
            slider.setValue(pos)
            slider.blockSignals(False)
        
        slider.valueChanged.connect(slider_to_spin)
        spinbox.valueChanged.connect(spin_to_slider)
        
        self.layout.addWidget(container)
        
        self.parameters[name] = {
            'slider': slider,
            'spinbox': spinbox,
            'min': min_val,
            'max': max_val,
        }
        
        return slider, spinbox
    
    def get_value(self, name: str) -> float:
        """Get parameter value."""
        return self.parameters[name]['spinbox'].value()
    
    def set_value(self, name: str, value: float):
        """Set parameter value."""
        self.parameters[name]['spinbox'].setValue(value)


class GeneratorThread(QThread):
    """Background thread for mesh generation."""
    
    progress = pyqtSignal(int, str)  # percentage, message
    finished = pyqtSignal(dict)      # result dict
    error = pyqtSignal(str)          # error message
    
    def __init__(self, parameters: dict):
        super().__init__()
        self.parameters = parameters
    
    def run(self):
        try:
            self.progress.emit(10, "Initializing gyroid math...")
            
            from mesh_generator import (
                GyroidMeshGenerator, CylinderGeometry, SphereGeometry,
                SpindleGeometry, ConeGeometry, CubeGeometry
            )
            from gyroid_math import GyroidSurface
            
            shape_type = self.parameters['shape_type']
            unit_cell = self.parameters['unit_cell']
            wall_thickness = self.parameters['wall_thickness']
            quality = self.parameters['quality']
            
            self.progress.emit(20, "Creating gyroid surface...")
            
            # Create gyroid surface
            threshold = GyroidSurface.threshold_from_wall_thickness(
                unit_cell, wall_thickness
            )
            gyroid = GyroidSurface(
                unit_cell=unit_cell,
                threshold=threshold
            )
            
            self.progress.emit(30, f"Creating {shape_type} geometry...")
            
            # Create geometry
            if shape_type == 'cylinder':
                geometry = CylinderGeometry(
                    diameter=self.parameters['diameter'],
                    height=self.parameters['height']
                )
            elif shape_type == 'sphere':
                geometry = SphereGeometry(
                    diameter=self.parameters['diameter']
                )
            elif shape_type == 'spindle':
                geometry = SpindleGeometry(
                    max_diameter=self.parameters['diameter'],
                    length=self.parameters['height']
                )
            elif shape_type == 'cone':
                geometry = ConeGeometry(
                    base_diameter=self.parameters['diameter'],
                    height=self.parameters['height']
                )
            elif shape_type == 'cube':
                geometry = CubeGeometry(
                    size=self.parameters['size']
                )
            
            self.progress.emit(40, "Generating mesh (marching cubes)...")
            
            # Generate mesh
            generator = GyroidMeshGenerator(geometry, gyroid, quality=quality)
            mesh = generator.generate()
            
            self.progress.emit(100, "Complete!")
            
            self.finished.emit({
                'mesh': mesh,
                'generator': generator,
                'unit_cell': unit_cell,
                'wall_thickness': wall_thickness,
                'threshold': threshold,
                'shape': shape_type,
                'quality': quality,
            })
            
        except Exception as e:
            logger.error(f"Generation error: {str(e)}", exc_info=True)
            self.error.emit(str(e))


class MainWindow(QMainWindow):
    """Main application window."""
    
    def __init__(self):
        super().__init__()
        self.current_mesh = None
        self.current_generator = None
        self.current_geometry_params = None
        self.current_gyroid_params = None
        
        self.init_ui()
        
        # Initial update
        self.update_estimated_properties()
    
    def init_ui(self):
        """Initialize the UI."""
        self.setWindowTitle("Gyroid Generator v3.0 — DLP Sorbent Design Platform")
        self.setGeometry(100, 100, 1500, 950)
        
        # Central widget
        central = QWidget()
        self.setCentralWidget(central)
        
        # Main horizontal layout
        main_layout = QHBoxLayout()
        central.setLayout(main_layout)
        
        left_panel = self.create_control_panel()
        main_layout.addWidget(left_panel, stretch=0)
        
        right_panel = self.create_preview_panel()
        main_layout.addWidget(right_panel, stretch=1)
        
        self.statusBar().showMessage("Ready")
    
    # ================================================================
    # CONTROL PANEL (left side)
    # ================================================================
    
    def create_control_panel(self) -> QWidget:
        """Create left control panel with scroll."""
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setMaximumWidth(380)
        
        panel = QWidget()
        layout = QVBoxLayout()
        panel.setLayout(layout)
        
        # === SHAPE SELECTION ===
        shape_group = QGroupBox("Shape Selection")
        shape_layout = QVBoxLayout()
        shape_group.setLayout(shape_layout)
        
        self.shape_selector = QComboBox()
        self.shape_selector.addItems([
            'Cylinder',
            'Sphere',
            'Spindle',
            'Cone',
            'Cube'
        ])
        self.shape_selector.currentTextChanged.connect(self.on_shape_changed)
        shape_layout.addWidget(self.shape_selector)
        
        layout.addWidget(shape_group)
        
        # === DIMENSIONS ===
        dim_group = QGroupBox("Dimensions")
        dim_layout = QVBoxLayout()
        dim_group.setLayout(dim_layout)
        
        self.diameter_spinbox = self.add_labeled_spinbox(
            dim_layout, "Diameter:", 1.0, 100.0, 4.6, "mm"
        )
        self.height_spinbox = self.add_labeled_spinbox(
            dim_layout, "Height/Length:", 1.0, 200.0, 50.0, "mm"
        )
        
        layout.addWidget(dim_group)
        
        # === GYROID PARAMETERS ===
        self.gyroid_group = ParameterGroup("Gyroid Parameters")
        
        self.gyroid_group.add_slider(
            'unit_cell', 
            'Unit Cell Size:',
            min_val=0.05, max_val=0.5, default=0.2, step=0.01,
            suffix='mm', decimals=3
        )
        
        wall_slider, wall_spinbox = self.gyroid_group.add_slider(
            'wall_thickness',
            'Wall Thickness:',
            min_val=0.02, max_val=0.15, default=0.06, step=0.001,
            suffix='mm', decimals=3
        )
        
        wall_spinbox.valueChanged.connect(self.update_estimated_properties)
        self.gyroid_group.parameters['unit_cell']['spinbox'].valueChanged.connect(
            self.update_estimated_properties
        )
        
        self.estimates_label = QLabel()
        self.estimates_label.setWordWrap(True)
        self.estimates_label.setStyleSheet(
            "background-color: #f0f0f0; padding: 10px; border-radius: 5px;"
        )
        self.gyroid_group.layout.addWidget(self.estimates_label)
        
        layout.addWidget(self.gyroid_group)
        
        # === MESH QUALITY ===
        quality_group = QGroupBox("Mesh Quality")
        quality_layout = QVBoxLayout()
        quality_group.setLayout(quality_layout)
        
        self.quality_selector = QComboBox()
        self.quality_selector.addItems(['Low (~50k tri)', 'Medium (~200k tri)', 
                                       'High (~1M tri)', 'Ultra (~5M tri)'])
        self.quality_selector.setCurrentIndex(1)
        quality_layout.addWidget(self.quality_selector)
        
        self.time_estimate_label = QLabel("⏱ Estimated time: ~45 seconds")
        self.time_estimate_label.setStyleSheet("color: #666;")
        quality_layout.addWidget(self.time_estimate_label)
        
        layout.addWidget(quality_group)
        
        # === GENERATE BUTTON ===
        self.generate_btn = QPushButton("Generate Gyroid")
        self.generate_btn.setStyleSheet("""
            QPushButton {
                background-color: #2563eb;
                color: white;
                padding: 12px;
                font-size: 14pt;
                font-weight: bold;
                border-radius: 8px;
            }
            QPushButton:hover { background-color: #1e40af; }
            QPushButton:pressed { background-color: #1e3a8a; }
        """)
        self.generate_btn.clicked.connect(self.on_generate_clicked)
        layout.addWidget(self.generate_btn)
        
        # Progress
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)
        
        self.progress_label = QLabel("")
        self.progress_label.setVisible(False)
        self.progress_label.setStyleSheet("color: #666; font-style: italic;")
        layout.addWidget(self.progress_label)
        
        # === ACTION BUTTONS ===
        self.export_btn = QPushButton("Export STL...")
        self.export_btn.setEnabled(False)
        self.export_btn.clicked.connect(self.on_export_clicked)
        layout.addWidget(self.export_btn)
        
        self.compare_btn = QPushButton("Save to Compare")
        self.compare_btn.setEnabled(False)
        self.compare_btn.clicked.connect(self.on_save_to_compare)
        self.compare_btn.setStyleSheet("""
            QPushButton { padding: 6px; }
            QPushButton:disabled { color: #aaa; }
        """)
        layout.addWidget(self.compare_btn)
        
        # Spacer
        layout.addStretch()
        
        scroll.setWidget(panel)
        return scroll
    
    # ================================================================
    # PREVIEW PANEL (right side — tabs)
    # ================================================================
    
    def create_preview_panel(self) -> QWidget:
        """Create right panel with all tabs."""
        panel = QWidget()
        layout = QVBoxLayout()
        panel.setLayout(layout)
        
        self.tabs = QTabWidget()
        
        # Tab 1: 3D Viewer (FIXED — built-in, no pyglet)
        self.viewer_3d = MeshViewer3D()
        self.tabs.addTab(self.viewer_3d, "3D View")
        
        # Tab 2: Statistics v2.0 (unchanged from v2)
        self.statistics_tab = StatisticsTabV2()
        self.tabs.addTab(self.statistics_tab, "Statistics")
        
        # Tab 3: Predictions (NEW — ΔP, Van Deemter, desorption, capacity)
        self.predictions_tab = PredictionsTab()
        self.tabs.addTab(self.predictions_tab, "Predictions")
        
        # Tab 4: Distributions (NEW — histograms, connectivity, ASA, etc.)
        self.distributions_tab = DistributionsTab()
        self.tabs.addTab(self.distributions_tab, "Distributions")
        
        # Tab 5: Cross-Section (NEW — 2D slices)
        self.cross_section_tab = CrossSectionTab()
        self.tabs.addTab(self.cross_section_tab, "Cross-Section")
        
        # Tab 6: Compare (NEW — side-by-side comparison)
        self.compare_tab = CompareTab()
        self.tabs.addTab(self.compare_tab, "Compare")
        
        layout.addWidget(self.tabs)
        
        return panel
    
    # ================================================================
    # HELPERS
    # ================================================================
    
    def add_labeled_spinbox(self, layout, label: str, 
                           min_val: float, max_val: float, 
                           default: float, suffix: str = "") -> QDoubleSpinBox:
        """Helper to add labeled spinbox."""
        container = QWidget()
        hlayout = QHBoxLayout()
        container.setLayout(hlayout)
        
        lbl = QLabel(label)
        lbl.setMinimumWidth(100)
        hlayout.addWidget(lbl)
        
        spinbox = QDoubleSpinBox()
        spinbox.setMinimum(min_val)
        spinbox.setMaximum(max_val)
        spinbox.setValue(default)
        spinbox.setSingleStep((max_val - min_val) / 100)
        spinbox.setDecimals(2)
        if suffix:
            spinbox.setSuffix(f" {suffix}")
        hlayout.addWidget(spinbox, stretch=1)
        
        layout.addWidget(container)
        
        return spinbox
    
    # ================================================================
    # EVENT HANDLERS
    # ================================================================
    
    def on_shape_changed(self, shape: str):
        """Handle shape selection change."""
        shape_lower = shape.lower()
        needs_height = shape_lower in ['cylinder', 'spindle', 'cone']
        self.height_spinbox.setEnabled(needs_height)
        logger.info(f"Shape changed to: {shape}")
    
    def update_estimated_properties(self):
        """Update estimated properties display."""
        from gyroid_math import GyroidSurface
        
        unit_cell = self.gyroid_group.get_value('unit_cell')
        wall_thickness = self.gyroid_group.get_value('wall_thickness')
        
        threshold = GyroidSurface.threshold_from_wall_thickness(
            unit_cell, wall_thickness
        )
        
        gyroid = GyroidSurface(unit_cell=unit_cell, threshold=threshold)
        props = gyroid.estimate_properties()
        
        text = (
            f"<b>Estimated Properties:</b><br>"
            f"Threshold (t): {threshold:.3f}<br>"
            f"Channel width: ~{props['channel_width_mm']*1000:.0f} µm<br>"
            f"Porosity: ~{props['porosity']:.1%}<br>"
            f"Surface density: {props['surface_density']:.2f} mm⁻¹"
        )
        
        self.estimates_label.setText(text)
    
    def on_generate_clicked(self):
        """Handle generate button click."""
        logger.info("Generate button clicked")
        
        params = {
            'shape_type': self.shape_selector.currentText().lower(),
            'diameter': self.diameter_spinbox.value(),
            'height': self.height_spinbox.value(),
            'size': max(self.diameter_spinbox.value(), self.height_spinbox.value()),
            'unit_cell': self.gyroid_group.get_value('unit_cell'),
            'wall_thickness': self.gyroid_group.get_value('wall_thickness'),
            'quality': ['low', 'medium', 'high', 'ultra'][
                self.quality_selector.currentIndex()
            ]
        }
        
        self.generate_btn.setEnabled(False)
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)
        self.progress_label.setVisible(True)
        self.progress_label.setText("Starting...")
        
        self.generator_thread = GeneratorThread(params)
        self.generator_thread.progress.connect(self.on_generation_progress)
        self.generator_thread.finished.connect(self.on_generation_finished)
        self.generator_thread.error.connect(self.on_generation_error)
        self.generator_thread.start()
    
    def on_generation_progress(self, percentage: int, message: str):
        """Handle generation progress updates."""
        self.progress_bar.setValue(percentage)
        self.progress_label.setText(message)
        self.statusBar().showMessage(f"{percentage}% - {message}")
    
    def on_generation_finished(self, result: dict):
        """Handle successful generation — distribute data to all tabs."""
        self.current_mesh = result['mesh']
        self.current_generator = result['generator']
        
        # Re-enable UI
        self.generate_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.progress_label.setVisible(False)
        self.export_btn.setEnabled(True)
        self.compare_btn.setEnabled(True)
        
        # Update 3D viewer
        try:
            self.viewer_3d.set_mesh(self.current_mesh)
        except Exception as e:
            logger.error(f"3D viewer error: {e}")
        
        # Get basic stats
        stats = result['generator'].get_statistics()
        
        # Parameters
        shape = result['shape']
        diameter = self.diameter_spinbox.value()
        height = self.height_spinbox.value()
        unit_cell = result['unit_cell']
        wall_thickness = result['wall_thickness']
        threshold = result['threshold']
        quality = result['quality']
        
        # Create geometry params
        if shape == 'cylinder':
            dimensions = {'diameter': diameter, 'height': height}
        elif shape == 'sphere':
            dimensions = {'diameter': diameter}
        elif shape == 'spindle':
            dimensions = {'diameter': diameter, 'height': height}
        elif shape == 'cone':
            dimensions = {'diameter': diameter, 'height': height}
        elif shape == 'cube':
            dimensions = {'diameter': diameter, 'height': diameter}
        else:
            dimensions = {'diameter': diameter, 'height': height}
        
        self.current_geometry_params = GeometryParams(shape=shape, dimensions=dimensions)
        self.current_gyroid_params = GyroidParams(
            unit_cell=unit_cell,
            wall_thickness=wall_thickness,
            threshold=threshold
        )
        
        # === Feed data to Statistics tab (v2 — unchanged) ===
        try:
            stats_analyzer = StatisticsV2(
                mesh=self.current_mesh,
                geometry=self.current_geometry_params,
                gyroid_params=self.current_gyroid_params
            )
            self.statistics_tab.set_statistics_analyzer(stats_analyzer)
            logger.info("Statistics v2.0 analyzer initialized")
        except Exception as e:
            logger.error(f"Statistics init error: {e}", exc_info=True)
        
        # === Feed data to Predictions tab (v3) ===
        try:
            basic_stats = {
                'porosity': stats.get('porosity', 0.5),
                'surface_area': stats.get('surface_area', self.current_mesh.area),
                'volume': stats.get('solid_volume', 
                          self.current_mesh.volume if self.current_mesh.is_watertight else 0),
                'vertices': stats.get('n_vertices', len(self.current_mesh.vertices)),
                'faces': stats.get('n_faces', len(self.current_mesh.faces)),
                'watertight': self.current_mesh.is_watertight,
            }
            
            # Try to get external SA from stats
            if 'external_surface_area' in stats:
                basic_stats['external_area_mm2'] = stats['external_surface_area']
            
            self.predictions_tab.set_data(
                mesh=self.current_mesh,
                geometry_params=self.current_geometry_params,
                gyroid_params=self.current_gyroid_params,
                basic_stats=basic_stats
            )
            logger.info("Predictions tab data set")
        except Exception as e:
            logger.error(f"Predictions init error: {e}", exc_info=True)
        
        # === Feed data to Distributions tab (v3) ===
        try:
            self.distributions_tab.set_data(
                mesh=self.current_mesh,
                container_geometry=shape,
                container_params=dimensions,
                unit_cell_mm=unit_cell,
            )
            logger.info("Distributions tab data set")
        except Exception as e:
            logger.error(f"Distributions init error: {e}", exc_info=True)
        
        # === Feed data to Cross-Section tab (v3) ===
        try:
            self.cross_section_tab.set_data(
                mesh=self.current_mesh,
                unit_cell_mm=unit_cell,
            )
            logger.info("Cross-section tab data set")
        except Exception as e:
            logger.error(f"Cross-section init error: {e}", exc_info=True)
        
        self.statusBar().showMessage("Generation complete!", 5000)
        logger.info(f"Mesh generated: {stats.get('n_vertices', '?'):,} vertices")
    
    def on_generation_error(self, error_msg: str):
        """Handle generation error."""
        self.generate_btn.setEnabled(True)
        self.progress_bar.setVisible(False)
        self.progress_label.setVisible(False)
        
        QMessageBox.critical(self, "Generation Error", 
                           f"Failed to generate mesh:\n\n{error_msg}")
        
        self.statusBar().showMessage("Generation failed", 5000)
        logger.error(f"Generation error: {error_msg}")
    
    def on_export_clicked(self):
        """Export mesh to STL file."""
        if self.current_mesh is None:
            return
        
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export STL", "gyroid.stl",
            "STL Files (*.stl);;All Files (*)"
        )
        
        if filename:
            try:
                self.current_generator.export_stl(filename)
                file_size = Path(filename).stat().st_size / (1024 * 1024)
                
                QMessageBox.information(
                    self, "Export Successful",
                    f"STL exported successfully!\n\n"
                    f"File: {filename}\n"
                    f"Size: {file_size:.1f} MB"
                )
                logger.info(f"STL exported: {filename}")
                
            except Exception as e:
                QMessageBox.critical(
                    self, "Export Error",
                    f"Failed to export STL:\n\n{str(e)}"
                )
                logger.error(f"Export error: {str(e)}")
    
    def on_save_to_compare(self):
        """Save current design to comparison table."""
        if self.current_mesh is None:
            return
        
        shape = self.shape_selector.currentText()
        uc = self.gyroid_group.get_value('unit_cell')
        wt = self.gyroid_group.get_value('wall_thickness')
        quality = ['Low', 'Medium', 'High', 'Ultra'][self.quality_selector.currentIndex()]
        
        name = f"{shape} UC={uc:.3f} WT={wt:.3f} ({quality})"
        
        params = {
            'shape': shape,
            'unit_cell': uc,
            'wall_thickness': wt,
            'quality': quality,
            'diameter': self.diameter_spinbox.value(),
            'height': self.height_spinbox.value(),
        }
        
        # Collect all available stats
        stats_dict = {}
        gen_stats = self.current_generator.get_statistics()
        if gen_stats:
            stats_dict['porosity'] = gen_stats.get('porosity')
            stats_dict['surface_area'] = gen_stats.get('surface_area', self.current_mesh.area)
            stats_dict['volume'] = gen_stats.get('solid_volume', 0)
            stats_dict['ssa_volumetric'] = gen_stats.get('ssa_volumetric')
        
        # Add distribution results if available
        dist_results = getattr(self.distributions_tab, '_results', {})
        if 'wall_thickness' in dist_results:
            r = dist_results['wall_thickness']
            stats_dict['wall_mean_um'] = r.get('mean_um')
            stats_dict['wall_p10_um'] = r.get('p10_um')
            stats_dict['wall_p90_um'] = r.get('p90_um')
            stats_dict['wall_uniformity'] = r.get('uniformity')
        
        if 'channel_width' in dist_results:
            r = dist_results['channel_width']
            stats_dict['channel_mean_um'] = r.get('mean_um')
            stats_dict['channel_p10_um'] = r.get('p10_um')
            stats_dict['channel_p90_um'] = r.get('p90_um')
        
        if 'connectivity' in dist_results:
            r = dist_results['connectivity']
            stats_dict['connected_fraction'] = r.get('connected_fraction')
            stats_dict['connected_porosity'] = r.get('connected_porosity')
            stats_dict['dead_end_fraction'] = r.get('dead_end_fraction')
        
        if 'throat' in dist_results:
            r = dist_results['throat']
            stats_dict['throat_mean_um'] = r.get('mean_um')
            stats_dict['throat_p10_um'] = r.get('p10_um')
            stats_dict['throat_to_pore_ratio'] = r.get('throat_to_pore_ratio')
        
        if 'printability' in dist_results:
            r = dist_results['printability']
            stats_dict['printability_score'] = r.get('score')
            stats_dict['printability_rating'] = r.get('rating')
        
        # Add prediction results
        pred_results = getattr(self.predictions_tab, '_results', {})
        if 'dp' in pred_results:
            stats_dict['delta_P_bar'] = pred_results['dp'].get('delta_P_bar')
        if 'vd' in pred_results:
            stats_dict['H_min_um'] = pred_results['vd'].get('H_min_um')
            stats_dict['N_per_m'] = pred_results['vd'].get('N_per_m')
            stats_dict['u_opt_mm_s'] = pred_results['vd'].get('u_opt_m_s', 0) * 1000
        if 'desorption' in pred_results:
            stats_dict['t_total_90_s'] = pred_results['desorption'].get('t_total_90_s')
        if 'capacity' in pred_results:
            stats_dict['Q_volumetric_mg_mL'] = pred_results['capacity'].get('Q_volumetric_mg_mL')
        
        self.compare_tab.add_design(name, params, stats_dict)
        self.tabs.setCurrentWidget(self.compare_tab)
        self.statusBar().showMessage(f"Saved '{name}' to comparison", 3000)


def main():
    """Main entry point."""
    app = QApplication(sys.argv)
    app.setStyle('Fusion')
    
    window = MainWindow()
    window.show()
    
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
