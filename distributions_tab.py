"""
Distributions Tab — GUI zakładka rozkładów i topologii
========================================================
Zawiera:
  - Histogram wall thickness
  - Histogram channel width
  - Connectivity / dead-end analysis
  - Accessible Surface Area (ASA)
  - Throat / constriction analysis
  - Printability score

Author: Claude (Anthropic)
Date: 2026-02-20
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QPushButton, QTextEdit, QSpinBox, QCheckBox, QProgressBar,
    QTabWidget, QFileDialog, QScrollArea, QDoubleSpinBox
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
import numpy as np
import time
import logging

logger = logging.getLogger(__name__)

import matplotlib
matplotlib.use('QtAgg')
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from distributions import VoxelAnalyzer


class AnalysisWorker(QThread):
    """Background thread for heavy voxel analysis."""
    finished = pyqtSignal(str, dict)  # (analysis_name, results)
    progress = pyqtSignal(str)  # status message
    sub_progress = pyqtSignal(str, int)  # (analysis_name, percentage)
    error = pyqtSignal(str, str)  # (analysis_name, error_message)

    def __init__(self, analyzer: VoxelAnalyzer, analyses: list, voxels_per_uc: int,
                 simplify: bool = False, target_faces: int = 200000):
        super().__init__()
        self.analyzer = analyzer
        self.analyses = analyses
        self.voxels_per_uc = voxels_per_uc
        self.simplify = simplify
        self.target_faces = target_faces

    def run(self):
        results_cache = {}

        if self.simplify:
            self.progress.emit("Simplifying mesh for analysis...")
            try:
                self.analyzer.simplify_analysis_mesh(self.target_faces)
            except Exception as e:
                logger.warning(f"Simplification failed, using original mesh: {e}")

        # Determine dynamic sample count based on mesh complexity
        # Target: ~10% of faces for sparse meshes, up to 50k for dense meshes
        n_faces = len(self.analyzer.mesh.faces)
        dynamic_n_samples = max(10000, min(50000, n_faces // 5))

        # Determine execution order (some depend on others)
        ordered_analyses = []
        for n in ['wall_thickness', 'channel_width', 'connectivity', 'asa', 'throat', 'printability']:
            if n in self.analyses:
                ordered_analyses.append(n)

        for name in ordered_analyses:
            try:
                self.progress.emit(f"Computing {name}...")

                def make_callback(n):
                    return lambda p: self.sub_progress.emit(n, p)

                if name == 'wall_thickness':
                    result = self.analyzer.calc_wall_thickness_distribution(
                        n_samples=dynamic_n_samples,
                        callback=make_callback(name)
                    )
                    results_cache['wall'] = result
                elif name == 'channel_width':
                    result = self.analyzer.calc_channel_width_distribution(
                        n_samples=dynamic_n_samples,
                        callback=make_callback(name)
                    )
                    results_cache['channel'] = result
                elif name == 'connectivity':
                    result = self.analyzer.calc_connectivity(self.voxels_per_uc)
                elif name == 'asa':
                    result = self.analyzer.calc_accessible_surface_area(
                        voxels_per_uc=self.voxels_per_uc
                    )
                elif name == 'throat':
                    # Depends on channel results
                    channel_res = results_cache.get('channel')
                    if not channel_res:
                        channel_res = self.analyzer.calc_channel_width_distribution(n_samples=10000)
                    result = self.analyzer.calc_throat_distribution(channel_res)
                elif name == 'printability':
                    # Depends on wall and channel
                    wall_res = results_cache.get('wall')
                    if not wall_res: wall_res = self.analyzer.calc_wall_thickness_distribution(n_samples=5000)
                    channel_res = results_cache.get('channel')
                    if not channel_res: channel_res = self.analyzer.calc_channel_width_distribution(n_samples=5000)

                    result = self.analyzer.calc_printability(
                        wall_results=wall_res,
                        channel_results=channel_res,
                        voxels_per_uc=self.voxels_per_uc
                    )
                else:
                    continue

                self.finished.emit(name, result)

            except Exception as e:
                logger.exception(f"Error in {name}")
                self.error.emit(name, str(e))


class DistributionsTab(QWidget):
    """Tab for distribution analysis and topology."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._mesh = None
        self._analyzer = None
        self._results = {}
        self._worker = None
        self._init_ui()

    def _init_ui(self):
        layout = QHBoxLayout(self)

        # Left: Controls
        controls_scroll = QScrollArea()
        controls_scroll.setWidgetResizable(True)
        controls_scroll.setMaximumWidth(300)
        controls_widget = QWidget()
        controls_layout = QVBoxLayout(controls_widget)

        # Resolution
        res_group = QGroupBox("Voxel Resolution")
        rg_layout = QVBoxLayout(res_group)
        h = QHBoxLayout()
        h.addWidget(QLabel("Voxels/unit cell:"))
        self.resolution_spin = QSpinBox()
        self.resolution_spin.setRange(6, 30)
        self.resolution_spin.setValue(12)
        h.addWidget(self.resolution_spin)
        rg_layout.addLayout(h)
        rg_layout.addWidget(QLabel("Higher = more accurate, slower"))
        controls_layout.addWidget(res_group)

        # Printer settings
        printer_group = QGroupBox("Printer Settings")
        pg_layout = QVBoxLayout(printer_group)
        h = QHBoxLayout()
        h.addWidget(QLabel("Pixel size [µm]:"))
        self.pixel_spin = QDoubleSpinBox()
        self.pixel_spin.setRange(1, 200)
        self.pixel_spin.setValue(22.0)
        self.pixel_spin.setDecimals(1)
        h.addWidget(self.pixel_spin)
        pg_layout.addLayout(h)
        controls_layout.addWidget(printer_group)

        # Mesh Optimization (New in v3.2)
        opt_group = QGroupBox("Mesh Optimization (Memory Safety)")
        og_layout = QVBoxLayout(opt_group)

        self.cb_simplify = QCheckBox("Simplify mesh for analysis")
        self.cb_simplify.setChecked(True)
        self.cb_simplify.setToolTip("Reduces triangle count before ray-casting. Recommended for meshes > 500k faces.")
        og_layout.addWidget(self.cb_simplify)

        h = QHBoxLayout()
        h.addWidget(QLabel("Target faces:"))
        self.target_faces_spin = QSpinBox()
        self.target_faces_spin.setRange(50000, 1000000)
        self.target_faces_spin.setSingleStep(50000)
        self.target_faces_spin.setValue(200000)
        h.addWidget(self.target_faces_spin)
        og_layout.addLayout(h)

        self.face_count_label = QLabel("Original: Unknown faces")
        og_layout.addWidget(self.face_count_label)

        controls_layout.addWidget(opt_group)

        # Analysis selection
        sel_group = QGroupBox("Select Analyses")
        sg_layout = QVBoxLayout(sel_group)

        self.cb_wall = QCheckBox("Wall Thickness (Ray Casting)")
        self.cb_wall.setChecked(True)
        self.cb_wall.setToolTip("Uses mesh-based ray casting to measure the distance from surface to surface through the solid phase.\n"
                                "High precision, independent of voxel resolution. Includes outlier filtering.")
        sg_layout.addWidget(self.cb_wall)

        self.cb_channel = QCheckBox("Channel Width (Ray Casting)")
        self.cb_channel.setChecked(True)
        self.cb_channel.setToolTip("Uses mesh-based ray casting to measure internal voids.\n"
                                   "Accurately captures channel dimensions without discretization artifacts.")
        sg_layout.addWidget(self.cb_channel)

        self.cb_connectivity = QCheckBox("Topological Connectivity")
        self.cb_connectivity.setChecked(True)
        self.cb_connectivity.setToolTip("Voxel-based analysis (scipy.ndimage) to detect connected void components.\n"
                                        "Restricted to the actual container volume. Checks for Z-axis through-paths.")
        sg_layout.addWidget(self.cb_connectivity)

        self.cb_asa = QCheckBox("Particulate Accessibility (ASA)")
        self.cb_asa.setChecked(True)
        self.cb_asa.setToolTip("Probe-based simulation of particle access using binary erosion of the masked void space.\n"
                               "Quantifies the fraction of surface area reachable by particles of specific radii (10-100 µm).")
        sg_layout.addWidget(self.cb_asa)

        self.cb_throat = QCheckBox("Throat / Clogging Analysis")
        self.cb_throat.setChecked(True)
        self.cb_throat.setToolTip("Analytical estimate of the narrowest constrictions (throats) in the channels.")
        sg_layout.addWidget(self.cb_throat)

        self.cb_printability = QCheckBox("Printability Score")
        self.cb_printability.setChecked(True)
        self.cb_printability.setToolTip("Evaluates if the design can be printed on the target DLP printer (robust percentile-based scoring).")
        sg_layout.addWidget(self.cb_printability)

        controls_layout.addWidget(sel_group)

        # Buttons
        self.btn_calculate = QPushButton("▶ Run Analysis")
        self.btn_calculate.setStyleSheet("QPushButton { font-weight: bold; padding: 8px; }")
        self.btn_calculate.clicked.connect(self._run_analysis)
        controls_layout.addWidget(self.btn_calculate)

        self.progress_label = QLabel("Ready")
        controls_layout.addWidget(self.progress_label)

        self.btn_export = QPushButton("Export Results CSV...")
        self.btn_export.clicked.connect(self._export_csv)
        controls_layout.addWidget(self.btn_export)

        controls_layout.addStretch()
        controls_scroll.setWidget(controls_widget)
        layout.addWidget(controls_scroll)

        # Right: Result tabs
        self.result_tabs = QTabWidget()

        # Wall thickness tab
        self.wall_widget = QWidget()
        wl = QVBoxLayout(self.wall_widget)
        self.wall_figure = Figure(figsize=(6, 4))
        self.wall_canvas = FigureCanvas(self.wall_figure)
        wl.addWidget(self.wall_canvas)
        self.wall_text = QTextEdit()
        self.wall_text.setReadOnly(True)
        self.wall_text.setMaximumHeight(120)
        wl.addWidget(self.wall_text)
        self.result_tabs.addTab(self.wall_widget, "Wall Thickness")

        # Channel width tab
        self.channel_widget = QWidget()
        cl = QVBoxLayout(self.channel_widget)
        self.channel_figure = Figure(figsize=(6, 4))
        self.channel_canvas = FigureCanvas(self.channel_figure)
        cl.addWidget(self.channel_canvas)
        self.channel_text = QTextEdit()
        self.channel_text.setReadOnly(True)
        self.channel_text.setMaximumHeight(120)
        cl.addWidget(self.channel_text)
        self.result_tabs.addTab(self.channel_widget, "Channel Width")

        # Connectivity tab
        self.conn_text = QTextEdit()
        self.conn_text.setReadOnly(True)
        self.result_tabs.addTab(self.conn_text, "Topology")

        # ASA tab
        self.asa_widget = QWidget()
        al = QVBoxLayout(self.asa_widget)
        self.asa_figure = Figure(figsize=(6, 4))
        self.asa_canvas = FigureCanvas(self.asa_figure)
        al.addWidget(self.asa_canvas)
        self.asa_text = QTextEdit()
        self.asa_text.setReadOnly(True)
        self.asa_text.setMaximumHeight(120)
        al.addWidget(self.asa_text)
        self.result_tabs.addTab(self.asa_widget, "Accessibility")

        # Throat tab
        self.throat_widget = QWidget()
        tl = QVBoxLayout(self.throat_widget)
        self.throat_figure = Figure(figsize=(6, 4))
        self.throat_canvas = FigureCanvas(self.throat_figure)
        tl.addWidget(self.throat_canvas)
        self.throat_text = QTextEdit()
        self.throat_text.setReadOnly(True)
        self.throat_text.setMaximumHeight(120)
        tl.addWidget(self.throat_text)
        self.result_tabs.addTab(self.throat_widget, "Throat")

        # Printability tab
        self.print_text = QTextEdit()
        self.print_text.setReadOnly(True)
        self.result_tabs.addTab(self.print_text, "Printability")

        layout.addWidget(self.result_tabs, stretch=1)

    def set_data(self, mesh, container_geometry, container_params,
                 unit_cell_mm, gyroid_params=None):
        """Set mesh data for analysis."""
        self._mesh = mesh
        n_faces = len(mesh.faces)
        self.face_count_label.setText(f"Original: {n_faces:,} faces")

        # Auto-suggest target faces (e.g., cap at 300k if mesh is huge)
        if n_faces > 500000:
            self.target_faces_spin.setValue(300000)
            self.cb_simplify.setChecked(True)
        else:
            self.cb_simplify.setChecked(False)

        self._analyzer = VoxelAnalyzer(
            mesh=mesh,
            container_geometry=container_geometry,
            container_params=container_params,
            unit_cell_mm=unit_cell_mm,
            printer_pixel_um=self.pixel_spin.value()
        )

    def _run_analysis(self):
        """Start analysis in background thread."""
        if self._analyzer is None:
            self.progress_label.setText("⚠ Generate a mesh first!")
            return

        # Update printer pixel
        self._analyzer.printer_pixel_um = self.pixel_spin.value()

        # Collect selected analyses
        analyses = []
        if self.cb_wall.isChecked():
            analyses.append('wall_thickness')
        if self.cb_channel.isChecked():
            analyses.append('channel_width')
        if self.cb_connectivity.isChecked():
            analyses.append('connectivity')
        if self.cb_asa.isChecked():
            analyses.append('asa')
        if self.cb_throat.isChecked():
            analyses.append('throat')
        if self.cb_printability.isChecked():
            analyses.append('printability')

        if not analyses:
            self.progress_label.setText("⚠ Select at least one analysis!")
            return

        self.btn_calculate.setEnabled(False)
        self.progress_label.setText("Running analysis...")

        self._worker = AnalysisWorker(
            self._analyzer, analyses, self.resolution_spin.value(),
            simplify=self.cb_simplify.isChecked(),
            target_faces=self.target_faces_spin.value()
        )
        self._worker.finished.connect(self._on_result)
        self._worker.progress.connect(self._on_progress)
        self._worker.sub_progress.connect(self._on_sub_progress)
        self._worker.error.connect(self._on_error)
        self._worker.finished.connect(self._check_done)
        self._pending = set(analyses)
        self._worker.start()

    def _on_progress(self, msg):
        self.progress_label.setText(msg)

    def _on_sub_progress(self, name, percentage):
        current_msg = self.progress_label.text().split(" (")[0]
        self.progress_label.setText(f"{current_msg} ({percentage}%)")

    def _on_error(self, name, error):
        self._pending.discard(name)
        logger.error(f"Analysis {name} failed: {error}")
        self.progress_label.setText(f"⚠ {name} failed: {error}")

    def _check_done(self, name, result):
        self._pending.discard(name)
        if not self._pending:
            self.btn_calculate.setEnabled(True)
            self.progress_label.setText("✅ All analyses complete!")

    def _on_result(self, name, result):
        """Handle analysis result."""
        self._results[name] = result

        if name == 'wall_thickness':
            self._display_wall_thickness(result)
        elif name == 'channel_width':
            self._display_channel_width(result)
        elif name == 'connectivity':
            self._display_connectivity(result)
        elif name == 'asa':
            self._display_asa(result)
        elif name == 'throat':
            self._display_throat(result)
        elif name == 'printability':
            self._display_printability(result)

    def _display_wall_thickness(self, r):
        if 'error' in r:
            self.wall_text.setText(f"⚠ {r['error']}")
            return

        # Plot histogram
        self.wall_figure.clear()
        ax = self.wall_figure.add_subplot(111)

        edges = r['hist_edges_um']
        centers = [(edges[i] + edges[i+1])/2 for i in range(len(edges)-1)]
        ax.bar(centers, r['hist_counts'], width=(edges[1]-edges[0])*0.9,
               color='#3498db', edgecolor='black', linewidth=0.5)

        # Mark percentiles
        ax.axvline(r['p5_um'], color='orange', linestyle=':', label=f"P5={r['p5_um']:.0f} µm")
        ax.axvline(r['p10_um'], color='red', linestyle='--', label=f"P10={r['p10_um']:.0f} µm")
        ax.axvline(r['p50_um'], color='green', linestyle='-', linewidth=2, label=f"P50={r['p50_um']:.0f} µm")
        ax.axvline(r['p90_um'], color='red', linestyle='--', label=f"P90={r['p90_um']:.0f} µm")

        ax.set_xlabel('Wall Thickness [µm]')
        ax.set_ylabel('Count')
        ax.set_title('Wall Thickness Distribution (Ray Casting)')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        self.wall_figure.tight_layout()
        self.wall_canvas.draw()

        # Text
        self.wall_text.setText(
            f"Mean: {r['mean_um']:.1f} µm | Std: {r['std_um']:.1f} µm\n"
            f"P1: {r['p1_um']:.1f} µm | P5: {r['p5_um']:.1f} µm | P10: {r['p10_um']:.1f} µm\n"
            f"P50: {r['p50_um']:.1f} µm | P90: {r['p90_um']:.1f} µm\n"
            f"Uniformity: {r['uniformity']:.2f} (1.0=perfect) | "
            f"Rays: {r['n_measurements']:,} | Time: {r['calculation_time']:.1f}s"
        )

    def _display_channel_width(self, r):
        if 'error' in r:
            self.channel_text.setText(f"⚠ {r['error']}")
            return

        self.channel_figure.clear()
        ax = self.channel_figure.add_subplot(111)

        edges = r['hist_edges_um']
        centers = [(edges[i] + edges[i+1])/2 for i in range(len(edges)-1)]
        ax.bar(centers, r['hist_counts'], width=(edges[1]-edges[0])*0.9,
               color='#e74c3c', edgecolor='black', linewidth=0.5)

        ax.axvline(r['p10_um'], color='blue', linestyle='--', label=f"P10={r['p10_um']:.0f} µm")
        ax.axvline(r['p50_um'], color='green', linestyle='-', linewidth=2, label=f"P50={r['p50_um']:.0f} µm")
        ax.axvline(r['p90_um'], color='blue', linestyle='--', label=f"P90={r['p90_um']:.0f} µm")

        ax.set_xlabel('Channel Width [µm]')
        ax.set_ylabel('Count')
        ax.set_title('Channel Width Distribution')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        self.channel_figure.tight_layout()
        self.channel_canvas.draw()

        self.channel_text.setText(
            f"Mean: {r['mean_um']:.1f} µm | Std: {r['std_um']:.1f} µm\n"
            f"Min: {r['min_um']:.1f} µm | Max: {r['max_um']:.1f} µm\n"
            f"P10: {r['p10_um']:.1f} µm | P50: {r['p50_um']:.1f} µm | P90: {r['p90_um']:.1f} µm\n"
            f"Uniformity: {r['uniformity']:.2f} | "
            f"Measurements: {r['n_measurements']:,} | Time: {r['calculation_time']:.1f}s"
        )

    def _display_connectivity(self, r):
        if 'error' in r:
            self.conn_text.setText(f"⚠ {r['error']}")
            return

        text = f"""Topological Connectivity (Voxel-based):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Note: This module checks for continuous topological paths (islands/bubbles).
  It does not calculate hydrodynamic dead volumes (stagnant flow zones).

  Connected components: {r['n_components']}
  Through-connected (inlet→outlet): {r['n_through_connected']}

  Through-connected void: {r['connected_fraction']:.1%}
  Dead-end void: {r['dead_end_fraction']:.1%}
  Isolated void (trapped): {r['isolated_fraction']:.1%}

  Connected porosity: {r['connected_porosity']:.3f}

  Voxel counts:
    Through-connected: {r['through_voxels']:,}
    Dead-end: {r['dead_end_voxels']:,}
    Isolated: {r['isolated_voxels']:,}
    Total void: {r['total_void_voxels']:,}

  Time: {r['calculation_time']:.1f}s
"""
        for w in r.get('warnings', []):
            text += f"\n  ⚠ {w}"

        self.conn_text.setText(text)

    def _display_asa(self, r):
        if 'error' in r:
            self.asa_text.setText(f"⚠ {r['error']}")
            return

        profiles = r['profiles']

        # Plot ASA vs probe radius
        self.asa_figure.clear()
        ax = self.asa_figure.add_subplot(111)

        radii = [p['probe_radius_um'] for p in profiles]
        fractions = [p['fraction'] * 100 for p in profiles]

        ax.plot(radii, fractions, 'bo-', linewidth=2, markersize=8)
        ax.fill_between(radii, fractions, alpha=0.2)

        ax.set_xlabel('Probe Radius [µm]')
        ax.set_ylabel('Accessible Surface Area [%]')
        ax.set_title('Particulate Accessibility vs Probe Size')
        ax.set_ylim(0, 105)
        ax.grid(True, alpha=0.3)
        self.asa_figure.tight_layout()
        self.asa_canvas.draw()

        # Text
        text = f"Total SA: {r['total_SA_mm2']:.2f} mm²\n\n"
        text += "Probe [µm]  │  ASA [mm²]     │  Fraction\n"
        text += "────────────┼────────────────┼──────────\n"
        for p in profiles:
            text += f"  {p['probe_radius_um']:6.1f}     │  {p['accessible_SA_mm2']:10.2f}   │  {p['fraction']:.1%}\n"
        text += f"\nTime: {r['calculation_time']:.1f}s"
        self.asa_text.setText(text)

    def _display_throat(self, r):
        if 'error' in r:
            self.throat_text.setText(f"⚠ {r.get('error', '')}  {r.get('note', '')}")
            return

        # No plot for analytical throat (optional)
        self.throat_figure.clear()
        ax = self.throat_figure.add_subplot(111)
        ax.text(0.5, 0.5, "Analytical Throat Estimation\n(0.7 × Channel Width)",
                ha='center', va='center', fontsize=12)
        self.throat_canvas.draw()

        text = (
            f"ESTIMATED THROAT SIZES:\n"
            f"━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
            f"  Mean Throat: {r['mean_um']:.1f} µm\n"
            f"  P10 Throat:  {r['p10_um']:.1f} µm (Clogging limit)\n"
            f"  P50 Throat:  {r['p50_um']:.1f} µm\n"
            f"  P90 Throat:  {r['p90_um']:.1f} µm\n\n"
            f"  Method: {r['method']}\n"
            f"  Time: {r['calculation_time']:.2f}s"
        )
        for w in r.get('warnings', []):
            text += f"\n  ⚠ {w}"
        self.throat_text.setText(text)

    def _display_printability(self, r):
        if 'error' in r:
            self.print_text.setText(f"⚠ {r['error']}")
            return

        # Color-coded score
        score = r['score']
        rating = r['rating']

        if score >= 85:
            color = "🟢"
        elif score >= 70:
            color = "🟡"
        elif score >= 50:
            color = "🟠"
        else:
            color = "🔴"

        text = f"""{color} PRINTABILITY SCORE: {score}/100 — {rating}
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Printer pixel: {r['printer_pixel_um']:.0f} µm

  Wall Thickness:
    Minimum: {r['min_wall_um']:.0f} µm
    Thin walls (<2 pixels): {r['pct_thin_walls']:.1f}%

  Channel Width:
    Minimum: {r['min_channel_um']:.0f} µm
    Narrow channels (<3 pixels): {r['pct_narrow_channels']:.1f}%

  Trapped Resin:
    Fraction: {r['trapped_resin_fraction']:.1%}
    Trapped components: {r['n_trapped_components']}

  Time: {r['calculation_time']:.1f}s
"""
        if r['issues']:
            text += "\nIssues:\n"
            for issue in r['issues']:
                text += f"  • {issue}\n"

        self.print_text.setText(text)

    def _export_csv(self):
        """Export all distribution results."""
        if not self._results:
            return

        path, _ = QFileDialog.getSaveFileName(
            self, "Export Distributions", "distributions.csv",
            "CSV files (*.csv)"
        )
        if not path:
            return

        try:
            with open(path, 'w', encoding='utf-8-sig') as f:
                f.write("Analysis,Metric,Value,Unit\n")

                for name, r in self._results.items():
                    if 'error' in r:
                        f.write(f"{name},error,{r['error']},\n")
                        continue

                    if name in ('wall_thickness', 'channel_width', 'throat'):
                        f.write(f"{name},mean,{r['mean_um']:.2f},um\n")
                        if 'std_um' in r:
                            f.write(f"{name},std,{r['std_um']:.2f},um\n")
                        f.write(f"{name},min,{r['min_um']:.2f},um\n")
                        f.write(f"{name},max,{r['max_um']:.2f},um\n")
                        f.write(f"{name},p10,{r['p10_um']:.2f},um\n")
                        f.write(f"{name},p50,{r['p50_um']:.2f},um\n")
                        f.write(f"{name},p90,{r['p90_um']:.2f},um\n")

                    elif name == 'connectivity':
                        f.write(f"connectivity,connected_fraction,{r['connected_fraction']:.4f},\n")
                        f.write(f"connectivity,dead_end_fraction,{r['dead_end_fraction']:.4f},\n")
                        f.write(f"connectivity,isolated_fraction,{r['isolated_fraction']:.4f},\n")
                        f.write(f"connectivity,connected_porosity,{r['connected_porosity']:.4f},\n")

                    elif name == 'printability':
                        f.write(f"printability,score,{r['score']},\n")
                        f.write(f"printability,rating,{r['rating']},\n")
                        f.write(f"printability,min_wall_um,{r['min_wall_um']:.1f},um\n")
                        f.write(f"printability,min_channel_um,{r['min_channel_um']:.1f},um\n")

            logger.info(f"Distributions exported to {path}")
        except Exception as e:
            logger.error(f"Export error: {e}")
