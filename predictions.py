"""
Predictions Module — Fizykochemiczne predykcje transportowe
=============================================================
Moduł obliczeniowy dla:
  - Ciśnienie wsteczne ΔP (Kozeny-Carman)
  - Van Deemter / HETP
  - Czas równowagi desorpcji (t90, t99)
  - Pojemność wiązania

Wzory i stałe z:
  - Guiochon, G. "Fundamentals of Preparative and Nonlinear Chromatography" (2006)
  - Knox, J.H. J. Chromatogr. A (1999) 831, 3-15
  - Schure, M.R. et al. J. Chromatogr. A (2004) 1031, 79-86
  - Crank, J. "Mathematics of Diffusion" (1975)

Author: Claude (Anthropic) for Mariusz's DLP research
Date: 2026-02-20
"""

import numpy as np
from typing import Dict, Optional, List, Tuple
import logging

logger = logging.getLogger(__name__)


# ==================== SOLVENT PRESETS ====================

SOLVENT_PRESETS = {
    'Water (25°C)': {
        'viscosity_Pa_s': 0.890e-3,
        'description': 'Pure water at 25°C'
    },
    'ACN/H₂O 20:80': {
        'viscosity_Pa_s': 0.807e-3,
        'description': 'Acetonitrile/water 20:80 v/v at 25°C'
    },
    'ACN/H₂O 50:50': {
        'viscosity_Pa_s': 0.573e-3,
        'description': 'Acetonitrile/water 50:50 v/v at 25°C'
    },
    'ACN/H₂O 80:20': {
        'viscosity_Pa_s': 0.398e-3,
        'description': 'Acetonitrile/water 80:20 v/v at 25°C'
    },
    'ACN 100%': {
        'viscosity_Pa_s': 0.343e-3,
        'description': 'Pure acetonitrile at 25°C'
    },
    'MeOH/H₂O 50:50': {
        'viscosity_Pa_s': 1.090e-3,
        'description': 'Methanol/water 50:50 v/v at 25°C'
    },
    'MeOH 100%': {
        'viscosity_Pa_s': 0.544e-3,
        'description': 'Pure methanol at 25°C'
    },
}

# ==================== ANALYTE PRESETS ====================

ANALYTE_PRESETS = {
    'Small molecule (~300 Da)': {
        'D_mobile': 1.0e-9,
        'description': 'e.g. drugs, pesticides, MW ~300'
    },
    'Small molecule (~500 Da)': {
        'D_mobile': 7.0e-10,
        'description': 'e.g. steroids, MW ~500'
    },
    'Peptide (~2 kDa)': {
        'D_mobile': 3.0e-10,
        'description': 'e.g. angiotensin, MW ~2000'
    },
    'Protein (~15 kDa)': {
        'D_mobile': 1.2e-10,
        'description': 'e.g. lysozyme, MW ~14300'
    },
    'Protein (~66 kDa)': {
        'D_mobile': 6.0e-11,
        'description': 'e.g. BSA, MW ~66000'
    },
}

# ==================== MATERIAL PRESETS ====================

POLYMER_PRESETS = {
    'Porous / Hydrogel (fast)': {
        'D_polymer': 1.0e-11,
        'description': 'Soft hydrogels or high nanoporosity materials'
    },
    'Swollen Resin (medium)': {
        'D_polymer': 1.0e-13,
        'description': 'Resin swollen by solvent (e.g. in MeOH)'
    },
    'Solid Crosslinked (slow)': {
        'D_polymer': 1.0e-15,
        'description': 'Standard solid DLP print, highly crosslinked'
    },
}


# ==================== PRESSURE DROP ====================

def calc_pressure_drop(
    porosity: float,
    hydraulic_diameter_m: float,
    column_length_m: float,
    viscosity_Pa_s: float,
    flow_rate_uL_min: float,
    column_diameter_m: float,
    k_geom: float = 32.0
) -> Dict:
    """
    Ciśnienie wsteczne z uogólnionego równania dla monolitów (Darcy-Weisbach / Hagen-Poiseuille).

    ΔP = (K_geom × η × u_superficial × L) / (d_h² × ε)

    Gdzie:
      K_geom: Geometryczna stała oporu (domyślnie 32 dla kapilar,
              dla gyroidów zazwyczaj 30-50).
      d_h: Średnica hydrauliczna = 4 × ε / S_volumetric

    Args:
        porosity: Porowatość [-] (0-1)
        hydraulic_diameter_m: Średnica hydrauliczna [m]
        column_length_m: Długość kolumny [m]
        viscosity_Pa_s: Lepkość fazy ruchomej [Pa·s]
        flow_rate_uL_min: Przepływ [µL/min]
        column_diameter_m: Średnica kolumny [m]
        k_geom: Stała geometryczna (Geometry Resistance Factor)

    Returns:
        Dict z ΔP w różnych jednostkach + parametry pośrednie
    """
    if porosity <= 0 or porosity >= 1:
        return {'error': 'Porosity must be between 0 and 1'}
    if hydraulic_diameter_m <= 0:
        return {'error': 'Hydraulic diameter must be positive'}

    # Konwersja przepływu: µL/min → m³/s
    flow_rate_m3_s = flow_rate_uL_min * 1e-9 / 60.0

    # Przekrój kolumny [m²]
    column_area_m2 = np.pi * (column_diameter_m / 2) ** 2

    # Prędkość powierzchniowa (superficial velocity) [m/s]
    u_superficial = flow_rate_m3_s / column_area_m2

    # Prędkość liniowa w porach (interstitial) [m/s]
    u_interstitial = u_superficial / porosity

    # Nowy wzór: ΔP = (K_geom × η × u_superficial × L) / (d_h² × ε)
    delta_P_Pa = (k_geom * viscosity_Pa_s * u_superficial * column_length_m) / \
                 (hydraulic_diameter_m ** 2 * porosity)

    # Konwersje jednostek
    delta_P_bar = delta_P_Pa / 1e5
    delta_P_mbar = delta_P_Pa / 100.0
    delta_P_MPa = delta_P_Pa / 1e6
    delta_P_psi = delta_P_Pa / 6894.76

    # Permeability (Darcy) [m²]
    # K = η × u × L / ΔP
    if delta_P_Pa > 0:
        permeability_m2 = viscosity_Pa_s * u_superficial * column_length_m / delta_P_Pa
    else:
        permeability_m2 = None

    results = {
        'delta_P_Pa': delta_P_Pa,
        'delta_P_bar': delta_P_bar,
        'delta_P_mbar': delta_P_mbar,
        'delta_P_MPa': delta_P_MPa,
        'delta_P_psi': delta_P_psi,
        'u_superficial_m_s': u_superficial,
        'u_interstitial_m_s': u_interstitial,
        'permeability_m2': permeability_m2,
        'k_geom': k_geom,
        'flow_rate_uL_min': flow_rate_uL_min,
        'viscosity_Pa_s': viscosity_Pa_s,
    }

    # Ostrzeżenia (tylko dla wysokich ciśnień)
    warnings = []
    if delta_P_bar > 400:
        warnings.append('ΔP > 400 bar — exceeds typical HPLC limit!')
    if delta_P_bar > 1000:
        warnings.append('ΔP > 1000 bar — exceeds UHPLC limit!')
    results['warnings'] = warnings

    return results


def calc_pressure_drop_curve(
    porosity: float,
    hydraulic_diameter_m: float,
    column_length_m: float,
    viscosity_Pa_s: float,
    column_diameter_m: float,
    k_geom: float = 32.0,
    flow_range_uL_min: Tuple[float, float] = (1.0, 5000.0),
    n_points: int = 50
) -> Dict:
    """
    Krzywa ΔP vs przepływ.

    Returns:
        Dict z tablicami flow_rates i delta_P_bar
    """
    flows = np.linspace(flow_range_uL_min[0], flow_range_uL_min[1], n_points)
    pressures = []

    for f in flows:
        result = calc_pressure_drop(
            porosity, hydraulic_diameter_m, column_length_m,
            viscosity_Pa_s, f, column_diameter_m
        )
        if 'error' in result:
            pressures.append(0)
        else:
            pressures.append(result['delta_P_bar'])

    return {
        'flow_rates_uL_min': flows.tolist(),
        'delta_P_bar': pressures,
    }


# ==================== VAN DEEMTER ====================

def calc_van_deemter(
    unit_cell_m: float,
    wall_thickness_m: float,
    hydraulic_diameter_m: float,
    porosity: float,
    D_mobile: float,
    D_polymer: float,
    tortuosity: float = 1.41,
    u_range: Optional[Tuple[float, float]] = None,
    n_points: int = 100
) -> Dict:
    """
    Krzywa Van Deemter: H = A + B/u + (C_s + C_m) × u

    Model dla struktur TPMS:
      A = 2·λ·d_p     (eddy diffusion)
      B = 2·γ·D_m/τ²   (longitudinal diffusion)
      C_s = (2/3)·k/(1+k)² · d_f²/D_s  (stationary phase mass transfer)
      C_m = ω · d_h²/D_m  (mobile phase mass transfer)

    Parametry:
      - λ ≈ 0.1 (regularna struktura TPMS, bardzo niskie rozmycie wirowe)
      - γ ≈ 0.6 (obstruction factor)
      - ω ≈ 0.1 (współczynnik strukturalny dla monolitów)

    Args:
        unit_cell_m: Rozmiar komórki [m]
        wall_thickness_m: Grubość ścianki [m]
        hydraulic_diameter_m: Średnica hydrauliczna kanałów [m]
        porosity: Porowatość [-]
        D_mobile: Współczynnik dyfuzji w fazie ruchomej [m²/s]
        D_polymer: Współczynnik dyfuzji w polimerze [m²/s]
        tortuosity: Tortuosity [-]
        u_range: Zakres prędkości liniowej (m/s). None = auto.
        n_points: Liczba punktów na krzywej

    Returns:
        Dict z krzywą Van Deemter i parametrami
    """
    # Parametry modelu
    lambda_eddy = 0.1    # Eddy diffusion parameter (regularna struktura TPMS)
    gamma_obstr = 0.6    # Obstruction factor for longitudinal diffusion
    omega_mobile = 0.1   # Structural factor for monoliths (mobile phase)

    # Efektywny rozmiar „ziarna" dla gyroidu ≈ unit cell
    d_p = unit_cell_m

    # Grubość filmu fazy stacjonarnej ≈ połowa grubości ścianki
    d_f = wall_thickness_m / 2.0

    # Średnica hydrauliczna kanałów
    d_h = hydraulic_diameter_m

    # Współczynnik retencji (retention factor)
    k_retention = 2.0

    # ---- Termy Van Deemter ----

    # A: Eddy diffusion
    A = 2 * lambda_eddy * d_p

    # B: Longitudinal (axial) diffusion
    D_eff_longitudinal = (gamma_obstr / tortuosity**2) * D_mobile
    B = 2 * D_eff_longitudinal

    # C_s: Mass transfer resistance in stationary phase
    C_s = (2.0/3.0) * (k_retention / (1 + k_retention)**2) * d_f**2 / D_polymer

    # C_m: Mass transfer resistance in mobile phase
    # Uproszczony wzór dla monolitów oparty na d_h
    C_m = omega_mobile * (d_h**2 / D_mobile)

    # C total
    C_total = C_s + C_m

    # ---- Automatyczny zakres prędkości ----
    if u_range is None:
        # u_opt ~ sqrt(B/C)
        if C_total > 0:
            u_opt_est = np.sqrt(B / C_total)
            u_min = u_opt_est * 0.1
            u_max = u_opt_est * 10.0
        else:
            u_min = 1e-5
            u_max = 1e-2
    else:
        u_min, u_max = u_range

    # Zapewnij fizycznie sensowny zakres
    u_min = max(u_min, 1e-7)
    u_max = min(u_max, 1.0)

    # ---- Oblicz krzywą ----
    u_values = np.linspace(u_min, u_max, n_points)

    H_A = np.full_like(u_values, A)
    H_B = B / u_values
    H_Cs = C_s * u_values
    H_Cm = C_m * u_values
    H_total = H_A + H_B + H_Cs + H_Cm

    # ---- Optimum ----
    if C_total > 0:
        u_opt = np.sqrt(B / C_total)
        H_min = A + 2 * np.sqrt(B * C_total)
    else:
        u_opt = u_values[np.argmin(H_total)]
        H_min = np.min(H_total)

    # Liczba półek teoretycznych na metr
    if H_min > 0:
        N_per_m = 1.0 / H_min
    else:
        N_per_m = None

    results = {
        # Parametry modelu
        'A_m': A,
        'B_m2_s': B,
        'C_s_s': C_s,
        'C_m_s': C_m,
        'C_total_s': C_total,

        # Optimum
        'u_opt_m_s': u_opt,
        'H_min_m': H_min,
        'H_min_um': H_min * 1e6,
        'N_per_m': N_per_m,

        # Krzywa
        'u_values_m_s': u_values.tolist(),
        'H_total_m': H_total.tolist(),
        'H_A_m': H_A.tolist(),
        'H_B_m': H_B.tolist(),
        'H_Cs_m': H_Cs.tolist(),
        'H_Cm_m': H_Cm.tolist(),

        # Parametry wejściowe
        'D_mobile': D_mobile,
        'D_polymer': D_polymer,
        'k_retention': k_retention,
        'lambda_eddy': lambda_eddy,
        'gamma_obstruction': gamma_obstr,
        'tortuosity': tortuosity,

        # Dominujący term
        'dominant_term': 'C_s (stationary)' if C_s > C_m else 'C_m (mobile)',
        'C_s_fraction': C_s / C_total if C_total > 0 else 0,
    }

    return results


# ==================== DESORPTION TIME ====================

def calc_desorption_time(
    wall_thickness_m: float,
    D_polymer: float,
    void_volume_mm3: float,
    flow_rate_uL_min: float
) -> Dict:
    """
    Czas równowagi desorpcji (model SPE).

    Model dwuetapowy:
    1) Dyfuzja przez ściankę polimeru (rate-limiting):
       t_x = (wall/2)² × α_x / D_polymer
       Gdzie α_x to Fourier number (Fo) dla desorpcji z połowy grubości ścianki.

    2) Wymywanie konwekcyjne kanałów (Convective Sweep):
       t_sweep = Void_Volume / Flow_Rate

    Args:
        wall_thickness_m: Grubość ścianki [m]
        D_polymer: Współczynnik dyfuzji w polimerze [m²/s]
        void_volume_mm3: Objętość kanałów [mm³]
        flow_rate_uL_min: Przepływ [µL/min]

    Returns:
        Dict z czasami desorpcji
    """
    # 1. Dyfuzja przez ściankę (Crank, half-thickness model)
    L_half = wall_thickness_m / 2.0

    # Zaktualizowane współczynniki alfa (Fo) dla połowy grubości (l)
    alpha = {
        '50%': 0.196,
        '90%': 0.848,
        '95%': 1.129,
        '99%': 1.781,
        '99.9%': 2.715,
    }

    t_wall = {}
    for pct, a in alpha.items():
        t_wall[pct] = L_half**2 * a / D_polymer

    # 2. Wymywanie kanałów (t_sweep)
    # 1 mm³ = 1 µL
    t_sweep_min = void_volume_mm3 / flow_rate_uL_min if flow_rate_uL_min > 0 else 0
    t_sweep_s = t_sweep_min * 60.0

    # Total
    t_total = {}
    for pct in alpha:
        t_total[pct] = t_wall[pct] + t_sweep_s

    # Rate-limiting step
    ratio = t_wall['90%'] / t_sweep_s if t_sweep_s > 0 else float('inf')

    if ratio > 10:
        rate_limiting = 'Wall diffusion (strongly dominant)'
    elif ratio > 2:
        rate_limiting = 'Wall diffusion (dominant)'
    elif ratio > 0.5:
        rate_limiting = 'Mixed (wall ≈ sweep)'
    else:
        rate_limiting = 'Convective sweep (dominant)'

    results = {
        't_wall_50_s': t_wall['50%'],
        't_wall_90_s': t_wall['90%'],
        't_wall_95_s': t_wall['95%'],
        't_wall_99_s': t_wall['99%'],
        't_sweep_s': t_sweep_s,
        't_total_90_s': t_total['90%'],
        't_total_95_s': t_total['95%'],
        't_total_99_s': t_total['99%'],
        'rate_limiting': rate_limiting,
        't_wall_to_sweep_ratio': ratio,

        # Formatted for display
        't_wall_90_formatted': _format_time(t_wall['90%']),
        't_total_90_formatted': _format_time(t_total['90%']),
        't_total_99_formatted': _format_time(t_total['99%']),
    }

    return results


def _format_time(seconds: float) -> str:
    """Format time in human-readable units."""
    if seconds < 1:
        return f"{seconds*1000:.1f} ms"
    elif seconds < 60:
        return f"{seconds:.1f} s"
    elif seconds < 3600:
        return f"{seconds/60:.1f} min"
    else:
        return f"{seconds/3600:.1f} h"


# ==================== BINDING CAPACITY ====================

def calc_binding_capacity(
    internal_surface_area_mm2: float,
    total_volume_mm3: float,
    solid_volume_mm3: float,
    polymer_density_g_cm3: float,
    q_surface_mg_m2: float,
    accessibility: float = 0.8,
    probe_accessible_fraction: float = 1.0
) -> Dict:
    """
    Pojemność wiązania analitu.

    Model:
      Q_vol = ASA × q_surface × accessibility / V_total
      Q_mass = Q_vol × V_total / mass

    Gdzie:
      ASA = internal_SA × probe_accessible_fraction
      q_surface = powierzchniowa gęstość wiązania [mg/m²]
        Typowe wartości:
          MIP (imprinted polymer): 1-5 mg/m²
          C18 silica: 2-4 mg/m²
          Ion exchange: 0.5-2 mg/m²
      accessibility = frakcja powierzchni realnie dostępna (0.5-1.0)
      probe_accessible_fraction = z analizy ASA (probe-based)

    Args:
        internal_surface_area_mm2: Wewnętrzna powierzchnia [mm²]
        total_volume_mm3: Objętość kontenera [mm³]
        solid_volume_mm3: Objętość fazy stałej [mm³]
        polymer_density_g_cm3: Gęstość polimeru [g/cm³]
        q_surface_mg_m2: Powierzchniowa pojemność [mg/m²]
        accessibility: Współczynnik dostępności [0-1]
        probe_accessible_fraction: Z ASA analysis [0-1]

    Returns:
        Dict z pojemnościami w różnych jednostkach
    """
    # Konwersje
    internal_SA_m2 = internal_surface_area_mm2 / 1e6
    total_volume_mL = total_volume_mm3 / 1000.0
    solid_volume_cm3 = solid_volume_mm3 / 1000.0
    mass_g = solid_volume_cm3 * polymer_density_g_cm3

    # Dostępna powierzchnia
    effective_SA_m2 = internal_SA_m2 * accessibility * probe_accessible_fraction

    # Pojemności
    total_capacity_mg = effective_SA_m2 * q_surface_mg_m2

    Q_volumetric = total_capacity_mg / total_volume_mL if total_volume_mL > 0 else 0  # mg/mL
    Q_gravimetric = total_capacity_mg / mass_g if mass_g > 0 else 0  # mg/g

    # Verdict
    if total_capacity_mg < 0.001:  # < 1 µg
        verdict = "Ultra-trace analysis / Micro-extraction only. Capacity is extremely low due to non-porous solid walls."
    elif total_capacity_mg < 0.1:   # 1 µg - 100 µg
        verdict = "Analytical scale. Suitable for highly sensitive detectors (MS/MS, fluorescence)."
    else:
        verdict = "Preparative / Standard SPE scale."

    results = {
        'total_capacity_mg': total_capacity_mg,
        'Q_volumetric_mg_mL': Q_volumetric,
        'Q_gravimetric_mg_g': Q_gravimetric,
        'effective_SA_m2': effective_SA_m2,
        'effective_SA_cm2': effective_SA_m2 * 1e4,
        'mass_g': mass_g,
        'q_surface_mg_m2': q_surface_mg_m2,
        'accessibility': accessibility,
        'probe_accessible_fraction': probe_accessible_fraction,
        'verdict': verdict
    }

    return results


# ==================== HYDRAULIC DIAMETER ====================

def calc_hydraulic_diameter(
    porosity: float,
    total_volume_mm3: float,
    internal_surface_area_mm2: float
) -> float:
    """
    Średnica hydrauliczna kanałów.

    d_h = 4 × ε × V_total / S_internal

    Ref: Standard hydraulic diameter definition for porous media.

    Args:
        porosity: Porowatość (0-1)
        total_volume_mm3: Objętość kontenera [mm³]
        internal_surface_area_mm2: Wewnętrzna powierzchnia [mm²]

    Returns:
        d_h w mm
    """
    if internal_surface_area_mm2 <= 0:
        return 0

    void_volume = porosity * total_volume_mm3
    d_h = 4.0 * void_volume / internal_surface_area_mm2

    return d_h
