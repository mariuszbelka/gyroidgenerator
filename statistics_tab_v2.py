"""
Statistics Tab Widget v2.0
===========================
Modular statistics UI with checkbox selection.

User can select which statistics to calculate:
- INSTANT calculations (always fast)
- ACCURATE calculations (optional, slower)

Author: Claude
Date: 2026-02-19
Version: 2.0
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox,
    QCheckBox, QRadioButton, QButtonGroup, QLabel,
    QPushButton, QTextEdit, QComboBox, QSpinBox,
    QProgressBar, QScrollArea
)
from PyQt6.QtCore import Qt, pyqtSignal, QThread
from PyQt6.QtGui import QFont
import logging

logger = logging.getLogger(__name__)


class StatisticsCalculationThread(QThread):
    """Background thread for statistics calculations."""

    progress = pyqtSignal(str, int)  # (message, percentage)
    finished = pyqtSignal(dict)  # Results dictionary
    error = pyqtSignal(str)  # Error message

    def __init__(self, stats_analyzer, selected_calcs):
        super().__init__()
        self.stats_analyzer = stats_analyzer
        self.selected_calcs = selected_calcs  # List of (name, function, args)

    def run(self):
        """Execute selected calculations."""
        try:
            results = {}
            total = len(self.selected_calcs)

            for i, (name, func, args) in enumerate(self.selected_calcs):
                self.progress.emit(f"Calculating {name}...", int((i/total) * 100))

                # Call calculation function
                result = func(**args) if args else func()
                results[name] = result

                self.progress.emit(f"✓ {name} complete", int(((i+1)/total) * 100))

            self.finished.emit(results)

        except Exception as e:
            logger.error(f"Statistics calculation error: {e}", exc_info=True)
            self.error.emit(str(e))


class StatisticsTabV2(QWidget):
    """
    Modular statistics tab with checkbox selection.

    Groups:
    1. Basic Stats (always shown, instant)
    2. Instant Calculations (checkboxes)
    3. Accurate Analysis (checkboxes with options)
    """

    def __init__(self):
        super().__init__()

        self.stats_analyzer = None  # Will be set after mesh generation
        self.calculation_thread = None
        self.results = {}

        self._init_ui()

    def _init_ui(self):
        """Initialize UI layout."""
        layout = QVBoxLayout()

        # Scroll area for long content
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll_widget = QWidget()
        scroll_layout = QVBoxLayout()

        # ===== BASIC STATISTICS =====
        self.basic_group = self._create_basic_stats_group()
        scroll_layout.addWidget(self.basic_group)

        # ===== INSTANT CALCULATIONS =====
        self.instant_group = self._create_instant_calcs_group()
        scroll_layout.addWidget(self.instant_group)

        # ===== ACCURATE ANALYSIS =====
        self.accurate_group = self._create_accurate_calcs_group()
        scroll_layout.addWidget(self.accurate_group)

        # ===== CONTROL BUTTONS =====
        self.controls = self._create_controls()
        scroll_layout.addWidget(self.controls)

        # ===== RESULTS DISPLAY =====
        self.results_text = QTextEdit()
        self.results_text.setReadOnly(True)
        self.results_text.setFont(QFont("Courier", 9))
        self.results_text.setPlaceholderText("Results will appear here...")
        scroll_layout.addWidget(self.results_text)

        scroll_layout.addStretch()
        scroll_widget.setLayout(scroll_layout)
        scroll.setWidget(scroll_widget)

        layout.addWidget(scroll)
        self.setLayout(layout)

    def _create_basic_stats_group(self) -> QGroupBox:
        """Basic statistics (always shown, instant)."""
        group = QGroupBox("Basic Statistics (Instant)")
        layout = QVBoxLayout()

        self.basic_stats_label = QLabel("Generate mesh first...")
        self.basic_stats_label.setWordWrap(True)
        self.basic_stats_label.setFont(QFont("Courier", 9))

        layout.addWidget(self.basic_stats_label)
        group.setLayout(layout)
        return group

    def _create_instant_calcs_group(self) -> QGroupBox:
        """Instant calculations group with checkboxes."""
        group = QGroupBox("Instant Calculations (<1 second)")
        layout = QVBoxLayout()

        info = QLabel("✓ Select which instant calculations to perform:")
        layout.addWidget(info)

        # Material density input
        density_layout = QHBoxLayout()
        density_layout.addWidget(QLabel("Material Density [g/cm³]:"))
        self.density_spin = QDoubleSpinBox()
        self.density_spin.setRange(0.1, 5.0)
        self.density_spin.setValue(1.05)
        self.density_spin.setSingleStep(0.01)
        self.density_spin.setDecimals(2)
        density_layout.addWidget(self.density_spin)
        density_layout.addStretch()
        layout.addLayout(density_layout)

        # Checkboxes for each instant calculation
        self.cb_ssa = QCheckBox("Specific Surface Area (m²/mL, m²/g)")
        self.cb_ssa.setChecked(True)  # Default ON

        self.cb_wall_theoretical = QCheckBox("Wall Thickness (theoretical)")
        self.cb_wall_theoretical.setChecked(True)

        self.cb_channel_theoretical = QCheckBox("Channel Width (theoretical)")
        self.cb_channel_theoretical.setChecked(True)

        self.cb_mass = QCheckBox("Mass (g)")
        self.cb_mass.setChecked(True)

        self.cb_external_area = QCheckBox("External Surface Area (with holes detection)")
        self.cb_external_area.setChecked(True)

        layout.addWidget(self.cb_ssa)
        layout.addWidget(self.cb_wall_theoretical)
        layout.addWidget(self.cb_channel_theoretical)
        layout.addWidget(self.cb_mass)
        layout.addWidget(self.cb_external_area)

        # Select/deselect all
        btn_layout = QHBoxLayout()
        btn_select_all = QPushButton("Select All")
        btn_select_all.clicked.connect(lambda: self._set_all_instant(True))
        btn_deselect_all = QPushButton("Deselect All")
        btn_deselect_all.clicked.connect(lambda: self._set_all_instant(False))
        btn_layout.addWidget(btn_select_all)
        btn_layout.addWidget(btn_deselect_all)
        layout.addLayout(btn_layout)

        group.setLayout(layout)
        return group

    def _create_accurate_calcs_group(self) -> QGroupBox:
        """Accurate calculations group (slower, optional)."""
        group = QGroupBox("Accurate Analysis (Slower - Optional)")
        layout = QVBoxLayout()

        warning = QLabel("⚠ These calculations may take 30s-2min depending on mesh size")
        warning.setStyleSheet("color: #ff9800;")
        layout.addWidget(warning)

        # Desorption Path Analysis
        desorption_label = QLabel("<b>Desorption Path Analysis (Internal diffusion):</b>")
        layout.addWidget(desorption_label)

        self.cb_desorption_fast = QCheckBox("Voxel-based (EDT solid-to-void)")
        self.cb_desorption_fast.setChecked(False)  # Default OFF
        fast_info = QLabel("  → ~5-15 seconds, full coverage")
        fast_info.setStyleSheet("color: #666; margin-left: 20px;")

        self.cb_desorption_sampling = QCheckBox("Statistical Sampling")
        self.cb_desorption_sampling.setChecked(False)  # Default OFF

        # Sampling options
        sampling_layout = QHBoxLayout()
        sampling_layout.addSpacing(20)
        sampling_layout.addWidget(QLabel("Samples:"))
        self.desorption_samples = QComboBox()
        self.desorption_samples.addItems(["100", "500", "1000", "2000"])
        self.desorption_samples.setCurrentText("1000")
        sampling_layout.addWidget(self.desorption_samples)
        sampling_layout.addWidget(QLabel("(more = slower but accurate)"))
        sampling_layout.addStretch()

        sampling_info = QLabel("  → ~30s-2min, ±5-10% accuracy")
        sampling_info.setStyleSheet("color: #666; margin-left: 20px;")

        layout.addWidget(self.cb_desorption_fast)
        layout.addWidget(fast_info)
        layout.addWidget(self.cb_desorption_sampling)
        layout.addLayout(sampling_layout)
        layout.addWidget(sampling_info)

        group.setLayout(layout)
        return group

    def _create_controls(self) -> QWidget:
        """Control buttons."""
        widget = QWidget()
        layout = QHBoxLayout()

        self.btn_calculate = QPushButton("📊 Calculate Selected")
        self.btn_calculate.clicked.connect(self.calculate_selected)
        self.btn_calculate.setEnabled(False)

        self.btn_export_csv = QPushButton("💾 Export CSV...")
        self.btn_export_csv.clicked.connect(self.export_csv)
        self.btn_export_csv.setEnabled(False)

        self.btn_clear = QPushButton("🗑 Clear Results")
        self.btn_clear.clicked.connect(self.clear_results)

        layout.addWidget(self.btn_calculate)
        layout.addWidget(self.btn_export_csv)
        layout.addWidget(self.btn_clear)
        layout.addStretch()

        widget.setLayout(layout)
        return widget

    def _set_all_instant(self, checked: bool):
        """Select/deselect all instant calculations."""
        self.cb_ssa.setChecked(checked)
        self.cb_wall_theoretical.setChecked(checked)
        self.cb_channel_theoretical.setChecked(checked)
        self.cb_mass.setChecked(checked)
        self.cb_external_area.setChecked(checked)

    def set_statistics_analyzer(self, analyzer):
        """
        Set the statistics analyzer (after mesh generation).

        Args:
            analyzer: StatisticsV2 instance
        """
        self.stats_analyzer = analyzer
        self.btn_calculate.setEnabled(True)

        # Show basic stats immediately
        self.show_basic_stats()

    def show_basic_stats(self):
        """Display basic statistics (instant)."""
        if not self.stats_analyzer:
            return

        basic = self.stats_analyzer.get_basic_stats()

        text = "BASIC MESH STATISTICS\n"
        text += "=" * 50 + "\n"
        text += f"Vertices:    {basic['vertices']:,}\n"
        text += f"Faces:       {basic['faces']:,}\n"
        text += f"Watertight:  {'Yes' if basic['watertight'] else 'No'}\n"

        if basic['porosity'] is not None:
            text += f"\nVolume:      {basic['total_volume']:.2f} mm³\n"
            text += f"Gyroid:      {basic['gyroid_volume']:.2f} mm³\n"
            text += f"Void:        {basic['void_volume']:.2f} mm³\n"
            text += f"Porosity:    {basic['porosity']:.1f}%\n"
        else:
            text += f"\nVolume:      {basic['total_volume']:.2f} mm³\n"
            text += "Porosity:    N/A (not watertight)\n"

        text += f"\nSurface Area: {basic['surface_area']:.2f} mm²\n"

        self.basic_stats_label.setText(text)

    def calculate_selected(self):
        """Calculate selected statistics."""
        if not self.stats_analyzer:
            return

        density = self.density_spin.value()

        # Build list of calculations to perform
        calcs = []

        # Instant calculations
        if self.cb_ssa.isChecked():
            calcs.append(('Specific Surface Area',
                         self.stats_analyzer.calc_specific_surface_area,
                         {}))

        if self.cb_wall_theoretical.isChecked():
            calcs.append(('Wall Thickness (Theoretical)',
                         self.stats_analyzer.calc_wall_thickness_theoretical,
                         {}))

        if self.cb_channel_theoretical.isChecked():
            calcs.append(('Channel Width (Theoretical)',
                         self.stats_analyzer.calc_channel_width_theoretical,
                         {}))

        if self.cb_mass.isChecked():
            calcs.append(('Mass',
                         self.stats_analyzer.calc_mass,
                         {'density': density}))

        if self.cb_external_area.isChecked():
            calcs.append(('External Surface Area',
                         self.stats_analyzer.calc_external_surface_area,
                         {}))

        # Accurate calculations
        if self.cb_desorption_fast.isChecked():
            calcs.append(('Desorption Path (FAST)',
                         self.stats_analyzer.calc_desorption_path_fast,
                         {}))

        if self.cb_desorption_sampling.isChecked():
            n_samples = int(self.desorption_samples.currentText())
            calcs.append(('Desorption Path (SAMPLING)',
                         self.stats_analyzer.calc_desorption_path_sampling,
                         {'n_samples': n_samples}))

        if len(calcs) == 0:
            self.results_text.setText("⚠ No statistics selected!")
            return

        # Disable button during calculation
        self.btn_calculate.setEnabled(False)
        self.btn_export_csv.setEnabled(False)

        # Start calculation thread
        self.calculation_thread = StatisticsCalculationThread(
            self.stats_analyzer, calcs
        )
        self.calculation_thread.progress.connect(self.on_progress)
        self.calculation_thread.finished.connect(self.on_calculation_finished)
        self.calculation_thread.error.connect(self.on_calculation_error)
        self.calculation_thread.start()

    def on_progress(self, message: str, percentage: int):
        """Handle progress updates."""
        # Could add progress bar here
        logger.info(f"Statistics progress: {message} ({percentage}%)")

    def on_calculation_finished(self, results: dict):
        """Handle calculation completion."""
        self.results = results
        self.display_results()

        self.btn_calculate.setEnabled(True)
        self.btn_export_csv.setEnabled(True)

    def on_calculation_error(self, error_msg: str):
        """Handle calculation error."""
        self.results_text.setText(f"❌ ERROR:\n{error_msg}")
        self.btn_calculate.setEnabled(True)

    def display_results(self):
        """Display calculation results."""
        text = "\n" + "=" * 70 + "\n"
        text += "STATISTICS RESULTS\n"
        text += "=" * 70 + "\n\n"

        for name, result in self.results.items():
            text += f"{'─' * 70}\n"
            text += f"{name}:\n"
            text += f"{'─' * 70}\n"

            if 'error' in result:
                text += f"  ❌ Error: {result['error']}\n\n"
                continue

            # Format result based on type
            if 'ssa_volumetric' in result:
                text += f"  Volumetric:   {result['ssa_volumetric']:.2f} m²/mL\n"
                if result['ssa_gravimetric']:
                    text += f"  Gravimetric:  {result['ssa_gravimetric']:.2f} m²/g\n"
                text += f"  ⏱ Time: {result['calculation_time']:.3f}s\n\n"

            elif 'mean_wall_thickness' in result:
                text += f"  Mean: {result['mean_wall_thickness']:.1f} µm\n"
                text += f"  Method: {result['method']}\n"
                if 'accuracy' in result:
                    text += f"  Accuracy: {result['accuracy']}\n"
                text += f"  ⏱ Time: {result['calculation_time']:.3f}s\n\n"

            elif 'mean_channel_width' in result:
                text += f"  Mean: {result['mean_channel_width']:.1f} µm\n"
                text += f"  Method: {result['method']}\n"
                text += f"  ⏱ Time: {result['calculation_time']:.3f}s\n\n"

            elif 'mean_path_mm' in result:
                # Desorption path results
                text += f"  Mean path:   {result['mean_path_mm']:.2f} mm\n"
                text += f"  Median path: {result['median_path_mm']:.2f} mm\n"
                text += f"  Std dev:     {result['std_path_mm']:.2f} mm\n"
                text += f"  Min path:    {result['min_path_mm']:.2f} mm\n"
                text += f"  Max path:    {result['max_path_mm']:.2f} mm\n"
                text += f"  Samples:     {result['samples']}\n"
                text += f"  Method:      {result['method']}\n"
                text += f"  Accuracy:    {result['accuracy']}\n"
                text += f"  ⏱ Time: {result['calculation_time']:.1f}s\n\n"

            elif 'tortuosity' in result:
                text += f"  τ = {result['tortuosity']:.2f}\n"
                text += f"  Method: {result['method']}\n"
                text += f"  ⏱ Time: {result['calculation_time']:.3f}s\n\n"

            elif 'mass_g' in result:
                if result['mass_g']:
                    text += f"  Mass: {result['mass_g']:.4f} g\n"
                    text += f"  Density: {result['density_g_cm3']:.2f} g/cm³\n"
                else:
                    text += f"  Mass: N/A (mesh not watertight)\n"
                text += f"  ⏱ Time: {result['calculation_time']:.3f}s\n\n"

            elif 'external_area_mm2' in result:
                text += f"  External: {result['external_area_mm2']:.1f} mm²\n"
                text += f"  Internal: {result['internal_area_mm2']:.1f} mm²\n"
                text += f"  Total:    {result['total_area_mm2']:.1f} mm²\n"
                text += f"  External fraction: {result['external_fraction_percent']:.1f}%\n"
                text += f"  ⏱ Time: {result['calculation_time']:.3f}s\n\n"

        text += "=" * 70 + "\n"
        text += "✅ All calculations complete!\n"

        self.results_text.setText(text)

    def export_csv(self):
        """Export results to CSV."""
        if not self.results:
            return

        from PyQt6.QtWidgets import QFileDialog

        filename, _ = QFileDialog.getSaveFileName(
            self, "Export Statistics", "", "CSV Files (*.csv)"
        )

        if filename:
            try:
                summary = self.stats_analyzer.export_summary_dict(self.results)

                import csv
                with open(filename, 'w', newline='', encoding='utf-8-sig') as f:
                    writer = csv.writer(f)
                    writer.writerow(['Parameter', 'Value'])
                    for key, value in summary.items():
                        writer.writerow([key, value])

                logger.info(f"Statistics exported to {filename}")
                self.results_text.append(f"\n✅ Exported to: {filename}")

            except Exception as e:
                logger.error(f"Export failed: {e}")
                self.results_text.append(f"\n❌ Export failed: {e}")

    def clear_results(self):
        """Clear results display."""
        self.results_text.clear()
        self.results = {}
