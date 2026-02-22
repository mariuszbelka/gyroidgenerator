"""
Predictions Tab — GUI zakładka predykcji fizykochemicznych
============================================================
Zawiera:
  - Ciśnienie wsteczne ΔP (wykres + wartość)
  - Van Deemter / HETP (wykres + optimum)
  - Czas równowagi desorpcji (bar chart)
  - Pojemność wiązania

Author: Claude (Anthropic)
Date: 2026-02-20
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QPushButton, QComboBox, QDoubleSpinBox, QTextEdit,
    QSplitter, QFileDialog, QScrollArea, QTabWidget
)
from PyQt6.QtCore import Qt, QThread, pyqtSignal
import numpy as np
import logging

logger = logging.getLogger(__name__)

# Import matplotlib with Qt backend
import matplotlib
matplotlib.use('QtAgg')
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from predictions import (
    SOLVENT_PRESETS, ANALYTE_PRESETS, POLYMER_PRESETS,
    calc_pressure_drop, calc_pressure_drop_curve,
    calc_van_deemter, calc_desorption_time,
    calc_binding_capacity, calc_hydraulic_diameter
)


class PredictionsTab(QWidget):
    """Tab for physicochemical predictions."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._mesh = None
        self._geometry_params = None
        self._gyroid_params = None
        self._basic_stats = None
        self._init_ui()

    def _init_ui(self):
        layout = QHBoxLayout(self)

        # Left: Controls
        controls_scroll = QScrollArea()
        controls_scroll.setWidgetResizable(True)
        controls_scroll.setMaximumWidth(320)
        controls_widget = QWidget()
        controls_layout = QVBoxLayout(controls_widget)

        # === Solvent ===
        solvent_group = QGroupBox("Solvent (Mobile Phase)")
        sg_layout = QVBoxLayout(solvent_group)

        self.solvent_combo = QComboBox()
        for name in SOLVENT_PRESETS:
            self.solvent_combo.addItem(name)
        sg_layout.addWidget(self.solvent_combo)

        self.viscosity_label = QLabel("η = 0.890 mPa·s")
        sg_layout.addWidget(self.viscosity_label)
        self.solvent_combo.currentTextChanged.connect(self._update_viscosity_label)

        controls_layout.addWidget(solvent_group)

        # === Analyte ===
        analyte_group = QGroupBox("Analyte")
        ag_layout = QVBoxLayout(analyte_group)

        self.analyte_combo = QComboBox()
        for name in ANALYTE_PRESETS:
            self.analyte_combo.addItem(name)
        self.analyte_combo.addItem("Manual")
        ag_layout.addWidget(self.analyte_combo)

        h = QHBoxLayout()
        h.addWidget(QLabel("D_mobile [m²/s]:"))
        self.d_mobile_spin = QDoubleSpinBox()
        self.d_mobile_spin.setDecimals(2)
        self.d_mobile_spin.setRange(0.01, 100.0)
        self.d_mobile_spin.setValue(1.0)
        self.d_mobile_spin.setSuffix(" ×10⁻⁹")
        h.addWidget(self.d_mobile_spin)
        ag_layout.addLayout(h)

        self.analyte_combo.currentTextChanged.connect(self._update_d_mobile)
        self._update_d_mobile(self.analyte_combo.currentText())

        controls_layout.addWidget(analyte_group)

        # === Material ===
        material_group = QGroupBox("Polymer / Material")
        mg_layout = QVBoxLayout(material_group)

        self.polymer_combo = QComboBox()
        for name in POLYMER_PRESETS:
            self.polymer_combo.addItem(name)
        self.polymer_combo.addItem("Custom")
        mg_layout.addWidget(self.polymer_combo)

        h = QHBoxLayout()
        h.addWidget(QLabel("D_polymer [m²/s]:"))
        self.d_polymer_spin = QDoubleSpinBox()
        self.d_polymer_spin.setDecimals(2)
        self.d_polymer_spin.setRange(0.000001, 1000.0)
        self.d_polymer_spin.setValue(1.0)
        self.d_polymer_spin.setSuffix(" ×10⁻¹³")
        h.addWidget(self.d_polymer_spin)
        mg_layout.addLayout(h)

        self.polymer_combo.currentTextChanged.connect(self._update_d_polymer)
        self._update_d_polymer(self.polymer_combo.currentText())

        h = QHBoxLayout()
        h.addWidget(QLabel("Density [g/cm³]:"))
        self.density_spin = QDoubleSpinBox()
        self.density_spin.setDecimals(3)
        self.density_spin.setRange(0.5, 3.0)
        self.density_spin.setValue(1.05)
        h.addWidget(self.density_spin)
        mg_layout.addLayout(h)

        controls_layout.addWidget(material_group)

        # === Flow ===
        flow_group = QGroupBox("Flow Conditions")
        fg_layout = QVBoxLayout(flow_group)

        h = QHBoxLayout()
        h.addWidget(QLabel("Flow rate [µL/min]:"))
        self.flow_spin = QDoubleSpinBox()
        self.flow_spin.setDecimals(1)
        self.flow_spin.setRange(0.1, 10000.0)
        self.flow_spin.setValue(100.0)
        h.addWidget(self.flow_spin)
        fg_layout.addLayout(h)

        h = QHBoxLayout()
        h.addWidget(QLabel("K_geom:"))
        self.k_geom_spin = QDoubleSpinBox()
        self.k_geom_spin.setRange(1.0, 1000.0)
        self.k_geom_spin.setValue(32.0)
        self.k_geom_spin.setToolTip("Geometry Resistance Factor. 32 = straight capillaries (Hagen-Poiseuille).\nTypical for gyroids: 30-50.")
        h.addWidget(self.k_geom_spin)
        fg_layout.addLayout(h)

        controls_layout.addWidget(flow_group)

        # === Binding ===
        binding_group = QGroupBox("Binding Capacity")
        bg_layout = QVBoxLayout(binding_group)

        h = QHBoxLayout()
        h.addWidget(QLabel("q_surface [mg/m²]:"))
        self.q_surface_spin = QDoubleSpinBox()
        self.q_surface_spin.setDecimals(2)
        self.q_surface_spin.setRange(0.01, 100.0)
        self.q_surface_spin.setValue(2.0)
        self.q_surface_spin.setToolTip("Surface binding density. Typical monolayer capacity for small molecules is ~0.5 - 2.0 mg/m².")
        h.addWidget(self.q_surface_spin)
        bg_layout.addLayout(h)

        h = QHBoxLayout()
        h.addWidget(QLabel("Accessibility [0-1]:"))
        self.accessibility_spin = QDoubleSpinBox()
        self.accessibility_spin.setDecimals(2)
        self.accessibility_spin.setRange(0.0, 1.0)
        self.accessibility_spin.setValue(0.8)
        self.accessibility_spin.setSingleStep(0.05)
        h.addWidget(self.accessibility_spin)
        bg_layout.addLayout(h)

        controls_layout.addWidget(binding_group)

        # === Buttons ===
        self.btn_calculate = QPushButton("▶ Calculate All Predictions")
        self.btn_calculate.setStyleSheet("QPushButton { font-weight: bold; padding: 8px; }")
        self.btn_calculate.clicked.connect(self._calculate_all)
        controls_layout.addWidget(self.btn_calculate)

        self.btn_export = QPushButton("Export Results CSV...")
        self.btn_export.clicked.connect(self._export_csv)
        controls_layout.addWidget(self.btn_export)

        controls_layout.addStretch()
        controls_scroll.setWidget(controls_widget)
        layout.addWidget(controls_scroll)

        # Right: Results (tabs with plots)
        self.result_tabs = QTabWidget()

        # ΔP tab
        self.dp_widget = QWidget()
        dp_layout = QVBoxLayout(self.dp_widget)
        self.dp_figure = Figure(figsize=(6, 4))
        self.dp_canvas = FigureCanvas(self.dp_figure)
        dp_layout.addWidget(self.dp_canvas)
        self.dp_text = QTextEdit()
        self.dp_text.setReadOnly(True)
        self.dp_text.setMaximumHeight(150)
        dp_layout.addWidget(self.dp_text)
        self.result_tabs.addTab(self.dp_widget, "ΔP")

        # Van Deemter tab
        self.vd_widget = QWidget()
        vd_layout = QVBoxLayout(self.vd_widget)
        self.vd_figure = Figure(figsize=(6, 4))
        self.vd_canvas = FigureCanvas(self.vd_figure)
        vd_layout.addWidget(self.vd_canvas)
        self.vd_text = QTextEdit()
        self.vd_text.setReadOnly(True)
        self.vd_text.setMaximumHeight(180)
        vd_layout.addWidget(self.vd_text)
        self.result_tabs.addTab(self.vd_widget, "Van Deemter")

        # Desorption tab
        self.des_widget = QWidget()
        des_layout = QVBoxLayout(self.des_widget)
        self.des_figure = Figure(figsize=(6, 4))
        self.des_canvas = FigureCanvas(self.des_figure)
        des_layout.addWidget(self.des_canvas)
        self.des_text = QTextEdit()
        self.des_text.setReadOnly(True)
        self.des_text.setMaximumHeight(150)
        des_layout.addWidget(self.des_text)
        self.result_tabs.addTab(self.des_widget, "Desorption")

        # Capacity tab
        self.cap_text = QTextEdit()
        self.cap_text.setReadOnly(True)
        self.result_tabs.addTab(self.cap_text, "Capacity")

        layout.addWidget(self.result_tabs, stretch=1)

        # Store results for export
        self._results = {}

    def set_data(self, mesh, geometry_params, gyroid_params, basic_stats):
        """Set mesh and parameters for calculations."""
        self._mesh = mesh
        self._geometry_params = geometry_params
        self._gyroid_params = gyroid_params
        self._basic_stats = basic_stats

    def _get_viscosity(self) -> float:
        name = self.solvent_combo.currentText()
        return SOLVENT_PRESETS[name]['viscosity_Pa_s']

    def _get_d_mobile(self) -> float:
        return self.d_mobile_spin.value() * 1e-9

    def _get_d_polymer(self) -> float:
        return self.d_polymer_spin.value() * 1e-13

    def _update_viscosity_label(self, name):
        if name in SOLVENT_PRESETS:
            eta = SOLVENT_PRESETS[name]['viscosity_Pa_s']
            self.viscosity_label.setText(f"η = {eta*1000:.3f} mPa·s")

    def _update_d_mobile(self, name):
        if name in ANALYTE_PRESETS:
            d = ANALYTE_PRESETS[name]['D_mobile']
            self.d_mobile_spin.setValue(d / 1e-9)
            self.d_mobile_spin.setEnabled(False)
        else:
            self.d_mobile_spin.setEnabled(True)

    def _update_d_polymer(self, name):
        if name in POLYMER_PRESETS:
            d = POLYMER_PRESETS[name]['D_polymer']
            self.d_polymer_spin.setValue(d / 1e-13)
            self.d_polymer_spin.setEnabled(False)
        else:
            self.d_polymer_spin.setEnabled(True)

    def _get_column_dims(self):
        """Extract column dimensions from geometry params."""
        if self._geometry_params is None:
            return None, None

        dims = self._geometry_params.dimensions
        shape = self._geometry_params.shape

        if shape == 'cylinder':
            return dims.get('diameter', 4.6) / 1000.0, dims.get('height', 50) / 1000.0
        elif shape == 'sphere':
            d = dims.get('diameter', 2.0) / 1000.0
            return d, d
        else:
            # Approximate
            d = dims.get('diameter', dims.get('width', 4.6)) / 1000.0
            h = dims.get('height', 50) / 1000.0
            return d, h

    def _calculate_all(self):
        """Run all predictions."""
        if self._mesh is None or self._basic_stats is None:
            self.dp_text.setText("⚠ Generate a mesh first!")
            return

        self.btn_calculate.setEnabled(False)
        self.btn_calculate.setText("Calculating...")

        try:
            self._calc_pressure_drop()
            self._calc_van_deemter()
            self._calc_desorption()
            self._calc_capacity()
        except Exception as e:
            logger.exception("Prediction error")
            self.dp_text.setText(f"⚠ Error: {e}")
        finally:
            self.btn_calculate.setEnabled(True)
            self.btn_calculate.setText("▶ Calculate All Predictions")

    def _calc_pressure_drop(self):
        """Calculate and plot ΔP."""
        stats = self._basic_stats
        col_d, col_L = self._get_column_dims()

        if col_d is None:
            self.dp_text.setText("⚠ Cannot determine column dimensions")
            return

        porosity = stats.get('porosity')
        if porosity is not None:
            porosity /= 100.0 # Convert % to fraction
        else:
            porosity = 0.5

        total_sa = stats.get('surface_area', 1.0)
        external_sa = stats.get('external_area_mm2', 0)
        if external_sa == 0:
            # Estimate from geometry
            if self._geometry_params.shape == 'cylinder':
                d = self._geometry_params.dimensions.get('diameter', 4.6)
                h = self._geometry_params.dimensions.get('height', 50)
                external_sa = np.pi * d * h + 2 * np.pi * (d/2)**2
            elif self._geometry_params.shape == 'sphere':
                d = self._geometry_params.dimensions.get('diameter', 2.0)
                external_sa = np.pi * d**2

        internal_sa = total_sa - external_sa
        if internal_sa <= 0:
            internal_sa = total_sa * 0.85  # fallback

        # Container volume
        if self._geometry_params.shape == 'cylinder':
            d_mm = self._geometry_params.dimensions.get('diameter', 4.6)
            h_mm = self._geometry_params.dimensions.get('height', 50)
            total_vol = np.pi * (d_mm/2)**2 * h_mm
        elif self._geometry_params.shape == 'sphere':
            d_mm = self._geometry_params.dimensions.get('diameter', 2.0)
            total_vol = (4/3) * np.pi * (d_mm/2)**3
        else:
            solid_vol = stats.get('volume')
            if solid_vol is not None and porosity < 1:
                total_vol = solid_vol / (1 - porosity)
            else:
                total_vol = 1.0

        d_h_mm = calc_hydraulic_diameter(porosity, total_vol, internal_sa)
        d_h_m = d_h_mm / 1000.0

        viscosity = self._get_viscosity()
        flow = self.flow_spin.value()
        k_geom = self.k_geom_spin.value()

        # Single point
        result = calc_pressure_drop(
            porosity, d_h_m, col_L, viscosity, flow, col_d, k_geom=k_geom
        )
        self._results['dp'] = result

        # Curve
        curve = calc_pressure_drop_curve(
            porosity, d_h_m, col_L, viscosity, col_d, k_geom=k_geom,
            flow_range_uL_min=(1.0, max(flow * 5, 1000.0))
        )

        # Plot
        self.dp_figure.clear()
        ax = self.dp_figure.add_subplot(111)
        ax.plot(curve['flow_rates_uL_min'], curve['delta_P_bar'], 'b-', linewidth=2)
        ax.axvline(flow, color='r', linestyle='--', alpha=0.7, label=f'Current: {flow} µL/min')

        # Only show limits if within range
        max_p = max(curve['delta_P_bar'])
        if max_p > 400:
            ax.axhline(400, color='orange', linestyle=':', alpha=0.5, label='HPLC limit (400 bar)')
        if max_p > 1000:
            ax.axhline(1000, color='red', linestyle=':', alpha=0.5, label='UHPLC limit (1000 bar)')

        ax.set_xlabel('Flow rate [µL/min]')
        ax.set_ylabel('ΔP [bar]')
        ax.set_title('Pressure Drop vs Flow Rate (Modified Darcy)')
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
        self.dp_figure.tight_layout()
        self.dp_canvas.draw()

        # Text
        if 'error' in result:
            self.dp_text.setText(f"⚠ {result['error']}")
            return

        p_display = f"{result['delta_P_bar']:.4f} bar"
        if result['delta_P_bar'] < 0.1:
            p_display = f"{result['delta_P_mbar']:.2f} mbar"

        text = f"""Pressure Drop Results:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  ΔP = {p_display}  ({result['delta_P_psi']:.2f} PSI)
  Flow: {flow:.1f} µL/min
  u_superficial: {result['u_superficial_m_s']*1000:.3f} mm/s
  u_interstitial: {result['u_interstitial_m_s']*1000:.3f} mm/s
  d_h: {d_h_mm*1000:.1f} µm
  Permeability: {result['permeability_m2']:.2e} m²
  Porosity: {porosity:.3f} | K_geom: {k_geom}
"""
        for w in result.get('warnings', []):
            text += f"\n  ⚠ {w}"

        self.dp_text.setText(text)

    def _calc_van_deemter(self):
        """Calculate and plot Van Deemter curve."""
        gp = self._gyroid_params
        if gp is None:
            self.vd_text.setText("⚠ No gyroid parameters")
            return

        # Use hydraulic diameter for Van Deemter as well for consistency
        stats = self._basic_stats
        porosity = stats.get('porosity')
        if porosity is not None:
            porosity /= 100.0
        else:
            porosity = 0.5
        total_sa = stats.get('surface_area', 1.0)
        # External area estimation
        if self._geometry_params.shape == 'cylinder':
            d = self._geometry_params.dimensions.get('diameter', 4.6)
            h = self._geometry_params.dimensions.get('height', 50)
            external_sa = np.pi * d * h + 2 * np.pi * (d/2)**2
        elif self._geometry_params.shape == 'sphere':
            d = self._geometry_params.dimensions.get('diameter', 2.0)
            external_sa = np.pi * d**2
        else:
            external_sa = total_sa * 0.15

        internal_sa = total_sa - external_sa
        if internal_sa <= 0: internal_sa = total_sa * 0.85

        if self._geometry_params.shape == 'cylinder':
            d_mm = self._geometry_params.dimensions.get('diameter', 4.6)
            h_mm = self._geometry_params.dimensions.get('height', 50)
            total_vol = np.pi * (d_mm/2)**2 * h_mm
        elif self._geometry_params.shape == 'sphere':
            d_mm = self._geometry_params.dimensions.get('diameter', 2.0)
            total_vol = (4/3) * np.pi * (d_mm/2)**3
        else:
            total_vol = stats.get('volume', 1.0) / (1 - porosity) if porosity < 1 else 1.0

        d_h_mm = calc_hydraulic_diameter(porosity, total_vol, internal_sa)
        d_h_m = d_h_mm / 1000.0

        unit_cell_m = gp.unit_cell / 1000.0
        wall_m = gp.wall_thickness / 1000.0

        result = calc_van_deemter(
            unit_cell_m=unit_cell_m,
            wall_thickness_m=wall_m,
            hydraulic_diameter_m=d_h_m,
            porosity=porosity,
            D_mobile=self._get_d_mobile(),
            D_polymer=self._get_d_polymer(),
        )
        self._results['vd'] = result

        # Plot
        self.vd_figure.clear()
        ax = self.vd_figure.add_subplot(111)

        u = np.array(result['u_values_m_s']) * 1000  # → mm/s
        H = np.array(result['H_total_m']) * 1e6  # → µm

        ax.plot(u, H, 'b-', linewidth=2, label='H total')
        ax.plot(u, np.array(result['H_A_m']) * 1e6, 'g--', alpha=0.6, label='A (eddy)')
        ax.plot(u, np.array(result['H_B_m']) * 1e6, 'r--', alpha=0.6, label='B (axial)')
        ax.plot(u, np.array(result['H_Cs_m']) * 1e6, 'm--', alpha=0.6, label='Cs (stat. phase)')
        ax.plot(u, np.array(result['H_Cm_m']) * 1e6, 'c--', alpha=0.6, label='Cm (mob. phase)')

        # Mark optimum
        u_opt = result['u_opt_m_s'] * 1000
        H_min = result['H_min_um']
        ax.plot(u_opt, H_min, 'ro', markersize=10, zorder=5,
                label=f'Optimum: H={H_min:.1f} µm @ u={u_opt:.2f} mm/s')

        ax.set_xlabel('Linear velocity u [mm/s]')
        ax.set_ylabel('HETP H [µm]')
        ax.set_title('Van Deemter Curve')
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(bottom=0)
        self.vd_figure.tight_layout()
        self.vd_canvas.draw()

        # Text
        col_d, col_L = self._get_column_dims()
        L_mm = col_L * 1000 if col_L else 50

        text = f"""Van Deemter Results:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  H_min = {result['H_min_um']:.1f} µm
  u_opt = {result['u_opt_m_s']*1000:.3f} mm/s
  N/m = {result['N_per_m']:.0f} plates/m (Note: For extraction, low N is common)

  For column L = {L_mm:.1f} mm:
    N = {result['N_per_m'] * col_L:.0f} theoretical plates

  Terms:
    A (eddy diffusion) = {result['A_m']*1e6:.1f} µm
    B (axial diffusion) = {result['B_m2_s']:.2e} m²/s
    C_s (stationary) = {result['C_s_s']*1000:.4f} ms
    C_m (mobile) = {result['C_m_s']*1000:.4f} ms
    Dominant: {result['dominant_term']}
    C_s fraction: {result['C_s_fraction']:.0%}

  D_mobile: {result['D_mobile']:.2e} m²/s
  D_polymer: {result['D_polymer']:.2e} m²/s
  k (retention): {result['k_retention']}
"""
        self.vd_text.setText(text)

    def _calc_desorption(self):
        """Calculate and plot desorption time."""
        gp = self._gyroid_params
        if gp is None:
            self.des_text.setText("⚠ No gyroid parameters")
            return

        wall_m = gp.wall_thickness / 1000.0
        stats = self._basic_stats

        # Void volume in mm3 (uL)
        porosity = stats.get('porosity')
        if porosity is not None:
            porosity /= 100.0
        else:
            porosity = 0.5
        if self._geometry_params.shape == 'cylinder':
            d_mm = self._geometry_params.dimensions.get('diameter', 4.6)
            h_mm = self._geometry_params.dimensions.get('height', 50)
            total_vol = np.pi * (d_mm/2)**2 * h_mm
        elif self._geometry_params.shape == 'sphere':
            d_mm = self._geometry_params.dimensions.get('diameter', 2.0)
            total_vol = (4/3) * np.pi * (d_mm/2)**3
        else:
            total_vol = stats.get('volume', 1.0) / (1 - porosity) if porosity < 1 else 1.0

        void_vol = total_vol * porosity
        flow = self.flow_spin.value()

        result = calc_desorption_time(
            wall_thickness_m=wall_m,
            D_polymer=self._get_d_polymer(),
            void_volume_mm3=void_vol,
            flow_rate_uL_min=flow
        )
        self._results['desorption'] = result

        # Plot: bar chart comparing wall vs channel time
        self.des_figure.clear()
        ax = self.des_figure.add_subplot(111)

        labels = ['Wall\n(90%)', 'Convective\nSweep', 'Total\n(90%)', 'Total\n(99%)']
        times = [
            result['t_wall_90_s'],
            result['t_sweep_s'],
            result['t_total_90_s'],
            result['t_total_99_s'],
        ]
        colors = ['#e74c3c', '#3498db', '#2ecc71', '#f39c12']

        bars = ax.bar(labels, times, color=colors, edgecolor='black', linewidth=0.5)

        # Add value labels on bars
        for bar, t in zip(bars, times):
            label = _format_time_short(t)
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height(),
                    label, ha='center', va='bottom', fontsize=9, fontweight='bold')

        ax.set_ylabel('Time')
        ax.set_title(f'Desorption Time Breakdown\n{result["rate_limiting"]}')

        # Smart Y-axis label
        max_t = max(times)
        if max_t < 1:
            ax.set_ylabel('Time [ms]')
            ax.set_yticklabels([f'{t*1000:.0f}' for t in ax.get_yticks()])
        elif max_t < 60:
            ax.set_ylabel('Time [s]')
        elif max_t < 3600:
            ax.set_ylabel('Time [min]')
            ax.set_yticklabels([f'{t/60:.1f}' for t in ax.get_yticks()])

        ax.grid(True, alpha=0.3, axis='y')
        self.des_figure.tight_layout()
        self.des_canvas.draw()

        # Text
        text = f"""Desorption Time Results:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Wall diffusion (Rate-limiting):
    t_50%  = {_format_time_short(result['t_wall_50_s'])}
    t_90%  = {_format_time_short(result['t_wall_90_s'])}
    t_95%  = {_format_time_short(result['t_wall_95_s'])}
    t_99%  = {_format_time_short(result['t_wall_99_s'])}

  Convective sweep (Channel):
    t_sweep = {_format_time_short(result['t_sweep_s'])}

  Total (wall + sweep):
    t_90%  = {result['t_total_90_formatted']}
    t_99%  = {result['t_total_99_formatted']}

  Rate-limiting step: {result['rate_limiting']}
  Wall/Sweep ratio: {result['t_wall_to_sweep_ratio']:.1f}×
"""
        self.des_text.setText(text)

    def _calc_capacity(self):
        """Calculate binding capacity."""
        stats = self._basic_stats
        if stats is None:
            self.cap_text.setText("⚠ No mesh data")
            return

        total_sa = stats.get('surface_area', 0)
        external_sa = stats.get('external_area_mm2', 0)
        if external_sa == 0:
            external_sa = total_sa * 0.15
        internal_sa = total_sa - external_sa

        porosity = stats.get('porosity')
        if porosity is not None:
            porosity /= 100.0
        else:
            porosity = 0.5

        solid_vol = stats.get('volume') or 0

        # Container volume
        if self._geometry_params.shape == 'cylinder':
            d_mm = self._geometry_params.dimensions.get('diameter', 4.6)
            h_mm = self._geometry_params.dimensions.get('height', 50)
            total_vol = np.pi * (d_mm/2)**2 * h_mm
        elif self._geometry_params.shape == 'sphere':
            d_mm = self._geometry_params.dimensions.get('diameter', 2.0)
            total_vol = (4/3) * np.pi * (d_mm/2)**3
        else:
            total_vol = solid_vol / (1 - porosity) if porosity < 1 else solid_vol * 2

        result = calc_binding_capacity(
            internal_surface_area_mm2=internal_sa,
            total_volume_mm3=total_vol,
            solid_volume_mm3=solid_vol,
            polymer_density_g_cm3=self.density_spin.value(),
            q_surface_mg_m2=self.q_surface_spin.value(),
            accessibility=self.accessibility_spin.value(),
        )
        self._results['capacity'] = result

        # Auto-scaling units for capacity
        cap_mg = result['total_capacity_mg']
        if cap_mg < 0.001:
            cap_display = f"{cap_mg*1e6:.2f} ng ({cap_mg*1000:.4f} µg)"
        elif cap_mg < 1.0:
            cap_display = f"{cap_mg*1000:.2f} µg ({cap_mg:.4f} mg)"
        else:
            cap_display = f"{cap_mg:.4f} mg"

        text = f"""Binding Capacity Results:
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Total capacity: {cap_display}

  Volumetric:   {result['Q_volumetric_mg_mL']:.4f} mg/mL
  Gravimetric:  {result['Q_gravimetric_mg_g']:.4f} mg/g

  Verdict: {result['verdict']}

  Parameters:
    Effective SA: {result['effective_SA_cm2']:.2f} cm²
    q_surface: {result['q_surface_mg_m2']:.2f} mg/m²
    Accessibility: {result['accessibility']:.0%}
    Mass: {result['mass_g']*1000:.2f} mg

  Reference values (commercial sorbents):
    C18 silica SPE:  5-20 mg/mL
    HLB polymer SPE: 10-40 mg/mL
    MIP:             1-10 mg/mL
"""
        self.cap_text.setText(text)

    def _export_csv(self):
        """Export all results to CSV."""
        if not self._results:
            return

        path, _ = QFileDialog.getSaveFileName(
            self, "Export Predictions", "predictions.csv",
            "CSV files (*.csv)"
        )
        if not path:
            return

        try:
            with open(path, 'w', encoding='utf-8-sig') as f:
                f.write("Metric,Value,Unit\n")

                if 'dp' in self._results:
                    r = self._results['dp']
                    if 'error' not in r:
                        f.write(f"delta_P_bar,{r['delta_P_bar']:.4f},bar\n")
                        f.write(f"delta_P_psi,{r['delta_P_psi']:.1f},PSI\n")
                        f.write(f"permeability,{r['permeability_m2']:.4e},m2\n")

                if 'vd' in self._results:
                    r = self._results['vd']
                    f.write(f"H_min,{r['H_min_um']:.2f},um\n")
                    f.write(f"u_opt,{r['u_opt_m_s']*1000:.4f},mm/s\n")
                    f.write(f"N_per_m,{r['N_per_m']:.0f},plates/m\n")
                    f.write(f"A_term,{r['A_m']*1e6:.2f},um\n")
                    f.write(f"C_s,{r['C_s_s']:.6e},s\n")
                    f.write(f"C_m,{r['C_m_s']:.6e},s\n")

                if 'desorption' in self._results:
                    r = self._results['desorption']
                    f.write(f"t_wall_90,{r['t_wall_90_s']:.4f},s\n")
                    f.write(f"t_channel,{r['t_channel_s']:.4f},s\n")
                    f.write(f"t_total_90,{r['t_total_90_s']:.4f},s\n")
                    f.write(f"t_total_99,{r['t_total_99_s']:.4f},s\n")

                if 'capacity' in self._results:
                    r = self._results['capacity']
                    f.write(f"Q_volumetric,{r['Q_volumetric_mg_mL']:.4f},mg/mL\n")
                    f.write(f"Q_gravimetric,{r['Q_gravimetric_mg_g']:.4f},mg/g\n")
                    f.write(f"total_capacity,{r['total_capacity_mg']:.6f},mg\n")

            logger.info(f"Predictions exported to {path}")
        except Exception as e:
            logger.error(f"Export error: {e}")


def _format_time_short(seconds: float) -> str:
    """Format time concisely."""
    if seconds < 0.001:
        return f"{seconds*1e6:.0f} µs"
    elif seconds < 1:
        return f"{seconds*1000:.1f} ms"
    elif seconds < 60:
        return f"{seconds:.1f} s"
    elif seconds < 3600:
        return f"{seconds/60:.1f} min"
    else:
        return f"{seconds/3600:.1f} h"
