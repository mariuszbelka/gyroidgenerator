"""
Compare Tab — Porównywarka projektów
=======================================
Pozwala zapisywać i porównywać do 4 geometrii obok siebie.

Author: Claude (Anthropic)
Date: 2026-02-20
"""

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QPushButton, QTextEdit, QFileDialog, QTableWidget,
    QTableWidgetItem, QHeaderView, QMessageBox
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor
import logging

logger = logging.getLogger(__name__)


class CompareTab(QWidget):
    """Tab for comparing multiple gyroid designs."""

    MAX_DESIGNS = 6

    def __init__(self, parent=None):
        super().__init__(parent)
        self._designs = []  # List of {'name': str, 'params': dict, 'stats': dict}
        self._init_ui()

    def _init_ui(self):
        layout = QVBoxLayout(self)

        # Top bar
        top = QHBoxLayout()
        top.addWidget(QLabel(f"Compare up to {self.MAX_DESIGNS} designs:"))

        self.btn_clear = QPushButton("Clear All")
        self.btn_clear.clicked.connect(self._clear_all)
        top.addWidget(self.btn_clear)

        self.btn_export = QPushButton("Export CSV...")
        self.btn_export.clicked.connect(self._export_csv)
        top.addWidget(self.btn_export)

        top.addStretch()
        self.count_label = QLabel("0 designs saved")
        top.addWidget(self.count_label)

        layout.addLayout(top)

        # Table
        self.table = QTableWidget()
        self.table.setAlternatingRowColors(True)
        layout.addWidget(self.table)

    def add_design(self, name: str, params: dict, stats: dict):
        """
        Add a design to comparison.

        Args:
            name: Display name (e.g. "Sphere 2mm UC=0.2 WT=0.06")
            params: {'shape': str, 'unit_cell': float, 'wall_thickness': float, ...}
            stats: Dict with all available metrics
        """
        if len(self._designs) >= self.MAX_DESIGNS:
            QMessageBox.warning(
                self, "Compare", f"Maximum {self.MAX_DESIGNS} designs. Clear some first."
            )
            return

        self._designs.append({
            'name': name,
            'params': params,
            'stats': stats,
        })

        self.count_label.setText(f"{len(self._designs)} designs saved")
        self._update_table()

    def _clear_all(self):
        self._designs.clear()
        self.count_label.setText("0 designs saved")
        self.table.clear()
        self.table.setRowCount(0)
        self.table.setColumnCount(0)

    def _update_table(self):
        """Rebuild comparison table."""
        if not self._designs:
            return

        # Collect all metric keys across all designs
        all_keys = []
        key_labels = {}

        # Parameter rows
        param_keys = [
            ('shape', 'Shape', ''),
            ('unit_cell', 'Unit Cell', 'mm'),
            ('wall_thickness', 'Wall Thickness', 'mm'),
            ('quality', 'Quality', ''),
        ]

        # Stat rows — ordered by importance
        stat_keys = [
            ('porosity', 'Porosity', ''),
            ('connected_porosity', 'Connected Porosity', ''),
            ('surface_area', 'Surface Area', 'mm²'),
            ('ssa_volumetric', 'SSA Volumetric', 'm²/mL'),
            ('volume', 'Solid Volume', 'mm³'),
            ('mass_g', 'Mass', 'g'),

            # Distributions
            ('wall_mean_um', 'Wall Thickness Mean', 'µm'),
            ('wall_p10_um', 'Wall Thickness P10', 'µm'),
            ('wall_p90_um', 'Wall Thickness P90', 'µm'),
            ('wall_uniformity', 'Wall Uniformity', ''),
            ('channel_mean_um', 'Channel Width Mean', 'µm'),
            ('channel_p10_um', 'Channel Width P10', 'µm'),
            ('channel_p90_um', 'Channel Width P90', 'µm'),

            # Connectivity
            ('connected_fraction', 'Through-Connected Void', ''),
            ('dead_end_fraction', 'Dead-End Void', ''),

            # Throat
            ('throat_mean_um', 'Throat Mean', 'µm'),
            ('throat_p10_um', 'Throat P10 (critical)', 'µm'),
            ('throat_to_pore_ratio', 'Throat/Pore Ratio', ''),

            # Printability
            ('printability_score', 'Printability Score', '/100'),
            ('printability_rating', 'Printability Rating', ''),

            # Predictions
            ('delta_P_bar', 'ΔP', 'bar'),
            ('H_min_um', 'HETP (H_min)', 'µm'),
            ('N_per_m', 'Plates/m (N/m)', ''),
            ('u_opt_mm_s', 'u_opt', 'mm/s'),
            ('t_total_90_s', 'Desorption t90', 's'),
            ('Q_volumetric_mg_mL', 'Capacity (vol)', 'mg/mL'),
        ]

        # Setup table
        n_designs = len(self._designs)
        n_rows = len(param_keys) + 1 + len(stat_keys)  # +1 for separator

        self.table.setRowCount(n_rows)
        self.table.setColumnCount(n_designs + 2)  # Metric, Unit, designs...

        headers = ['Metric', 'Unit'] + [d['name'] for d in self._designs]
        self.table.setHorizontalHeaderLabels(headers)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)

        row = 0

        # Parameters section
        for key, label, unit in param_keys:
            self.table.setItem(row, 0, QTableWidgetItem(label))
            self.table.setItem(row, 1, QTableWidgetItem(unit))
            for col, design in enumerate(self._designs):
                val = design['params'].get(key, '')
                self.table.setItem(row, col + 2, QTableWidgetItem(str(val)))
            row += 1

        # Separator
        for col in range(n_designs + 2):
            item = QTableWidgetItem("─" * 20)
            item.setBackground(QColor(230, 230, 230))
            self.table.setItem(row, col, item)
        row += 1

        # Stats section
        for key, label, unit in stat_keys:
            self.table.setItem(row, 0, QTableWidgetItem(label))
            self.table.setItem(row, 1, QTableWidgetItem(unit))

            values = []
            for design in self._designs:
                val = design['stats'].get(key, None)
                values.append(val)

            # Format and highlight best/worst
            for col, val in enumerate(values):
                if val is None:
                    self.table.setItem(row, col + 2, QTableWidgetItem("—"))
                elif isinstance(val, float):
                    item = QTableWidgetItem(f"{val:.4g}")
                    self.table.setItem(row, col + 2, item)
                else:
                    self.table.setItem(row, col + 2, QTableWidgetItem(str(val)))

            # Color-code: green for best, red for worst
            numeric_vals = [(i, v) for i, v in enumerate(values)
                           if v is not None and isinstance(v, (int, float))]
            if len(numeric_vals) >= 2:
                # Determine if higher or lower is better
                higher_better = key in (
                    'porosity', 'connected_porosity', 'surface_area',
                    'ssa_volumetric', 'connected_fraction',
                    'wall_uniformity', 'throat_to_pore_ratio',
                    'printability_score', 'N_per_m', 'Q_volumetric_mg_mL'
                )
                lower_better = key in (
                    'dead_end_fraction', 'delta_P_bar', 'H_min_um',
                    't_total_90_s'
                )

                sorted_vals = sorted(numeric_vals, key=lambda x: x[1])
                if higher_better:
                    best_idx = sorted_vals[-1][0]
                    worst_idx = sorted_vals[0][0]
                elif lower_better:
                    best_idx = sorted_vals[0][0]
                    worst_idx = sorted_vals[-1][0]
                else:
                    best_idx = worst_idx = None

                if best_idx is not None:
                    self.table.item(row, best_idx + 2).setBackground(
                        QColor(200, 255, 200))  # Light green
                if worst_idx is not None and worst_idx != best_idx:
                    self.table.item(row, worst_idx + 2).setBackground(
                        QColor(255, 200, 200))  # Light red

            row += 1

    def _export_csv(self):
        """Export comparison table as CSV."""
        if not self._designs:
            return

        path, _ = QFileDialog.getSaveFileName(
            self, "Export Comparison", "comparison.csv", "CSV files (*.csv)"
        )
        if not path:
            return

        try:
            with open(path, 'w', encoding='utf-8-sig') as f:
                # Header
                names = [d['name'] for d in self._designs]
                f.write("Metric,Unit," + ",".join(names) + "\n")

                # Rows from table
                for row in range(self.table.rowCount()):
                    cells = []
                    for col in range(self.table.columnCount()):
                        item = self.table.item(row, col)
                        cells.append(item.text() if item else "")
                    f.write(",".join(cells) + "\n")

            logger.info(f"Comparison exported to {path}")
        except Exception as e:
            logger.error(f"Export error: {e}")
