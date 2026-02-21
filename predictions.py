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


# ==================== PRESSURE DROP ====================

def calc_pressure_drop(
    porosity: float,
    hydraulic_diameter_m: float,
    column_length_m: float,
    viscosity_Pa_s: float,
    flow_rate_uL_min: float,
    column_diameter_m: float
) -> Dict:
    """
    Ciśnienie wsteczne z równania Kozeny-Carman.

    ΔP = (150 × η × u × L) / (d_h² × ε³)

    Dla struktur gyroidalnych stosujemy zmodyfikowaną wersję
    z hydraulicznym średnikiem kanału d_h obliczanym z geometrii:
        d_h = 4 × ε × V_total / S_internal

    Refs:
      - Kozeny (1927), Carman (1937)
      - Guiochon, "Fundamentals of Preparative and Nonlinear Chromatography"

    Args:
        porosity: Porowatość [-] (0-1)
        hydraulic_diameter_m: Średnica hydrauliczna [m]
        column_length_m: Długość kolumny [m]
        viscosity_Pa_s: Lepkość fazy ruchomej [Pa·s]
        flow_rate_uL_min: Przepływ [µL/min]
        column_diameter_m: Średnica kolumny [m]

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

    # Kozeny-Carman: ΔP = K_KC × η × u × L / (d_h² × ε³)
    # K_KC = 150 dla packed beds (Blake-Kozeny)
    # Dla monolitycznych struktur TPMS: K_KC ≈ 120-180
    # Używamy 150 jako standardową wartość
    K_KC = 150.0

    delta_P_Pa = (K_KC * viscosity_Pa_s * u_superficial * column_length_m) / \
                 (hydraulic_diameter_m ** 2 * porosity ** 3)

    # Konwersje jednostek
    delta_P_bar = delta_P_Pa / 1e5
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
        'delta_P_MPa': delta_P_MPa,
        'delta_P_psi': delta_P_psi,
        'u_superficial_m_s': u_superficial,
        'u_interstitial_m_s': u_interstitial,
        'permeability_m2': permeability_m2,
        'kozeny_carman_const': K_KC,
        'flow_rate_uL_min': flow_rate_uL_min,
        'viscosity_Pa_s': viscosity_Pa_s,
    }

    # Ostrzeżenia
    warnings = []
    if delta_P_bar > 400:
        warnings.append('ΔP > 400 bar — exceeds typical HPLC limit!')
    if delta_P_bar > 1000:
        warnings.append('ΔP > 1000 bar — exceeds UHPLC limit!')
    if delta_P_bar > 1500:
        warnings.append('ΔP > 1500 bar — EXTREME! Not feasible.')
    results['warnings'] = warnings

    return results


def calc_pressure_drop_curve(
    porosity: float,
    hydraulic_diameter_m: float,
    column_length_m: float,
    viscosity_Pa_s: float,
    column_diameter_m: float,
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
    channel_width_m: float,
    porosity: float,
    D_mobile: float,
    D_polymer: float,
    tortuosity: float = 1.41,
    u_range: Optional[Tuple[float, float]] = None,
    n_points: int = 100
) -> Dict:
    """
    Krzywa Van Deemter: H = A + B/u + (C_s + C_m) × u

    Model:
      A = 2·λ·d_p     (eddy diffusion / multiple paths)
      B = 2·γ·D_m      (longitudinal/axial diffusion)
      C_s = f_s × d_f² / D_s   (mass transfer in stationary phase)
      C_m = f_m × d_c² / D_m   (mass transfer in mobile phase)

    Parametry dla gyroidu TPMS:
      - λ ≈ 0.5 (regularna struktura, niski eddy diffusion)
        Ref: Knox (1999), dla monolitów regularnych
      - γ ≈ 0.6 (obstructive factor)
        Ref: Torquato (2002), dla porowatych mediów ε~0.5
      - f_s = 1/12 (geometria: dyfuzja z płaskiej płytki)
        Ref: Crank (1975), rozwiązanie dla płytki grubości d_f
      - f_m = 1/96 (dyfuzja w kanale cylindrycznym)
        Ref: Golay (1958) / Aris (1956)

    Args:
        unit_cell_m: Rozmiar komórki [m]
        wall_thickness_m: Grubość ścianki [m]
        channel_width_m: Szerokość kanału [m]
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
    lambda_eddy = 0.5    # Eddy diffusion parameter (regularna struktura)
    gamma_obstr = 0.6    # Obstruction factor for longitudinal diffusion
    f_stationary = 1/12  # Film mass transfer factor (flat plate geometry)
    f_mobile = 1/96      # Mobile phase mass transfer factor (cylindrical channel)

    # Efektywny rozmiar „ziarna" dla gyroidu ≈ unit cell
    d_p = unit_cell_m

    # Grubość filmu fazy stacjonarnej ≈ połowa grubości ścianki
    d_f = wall_thickness_m / 2.0

    # Charakterystyczny wymiar kanału
    d_c = channel_width_m

    # Współczynnik retencji (retention factor)
    # k = (1-ε)/ε × K_distribution
    # Dla ogólnej analizy zakładamy k = 2 (typowa wartość)
    k_retention = 2.0

    # ---- Termy Van Deemter ----

    # A: Eddy diffusion (multiple paths)
    # Dla regularnej struktury gyroidu λ jest niskie
    A = 2 * lambda_eddy * d_p

    # B: Longitudinal (axial) diffusion
    # D_eff = (γ/τ²) × D_m (w porach) + obstructed solid contribution
    D_eff_longitudinal = (gamma_obstr / tortuosity**2) * D_mobile
    B = 2 * D_eff_longitudinal

    # C_s: Mass transfer resistance in stationary phase
    # C_s = (2/3) × k/(1+k)² × d_f²/(D_s)
    # Ref: van Deemter (1956), model płaskiej płytki
    C_s = (2.0/3.0) * (k_retention / (1 + k_retention)**2) * d_f**2 / D_polymer

    # C_m: Mass transfer resistance in mobile phase
    # C_m = f(k) × d_c² / D_m
    # f(k) = (1 + 6k + 11k²) / (96 × (1+k)²) — Golay equation
    f_k = (1 + 6*k_retention + 11*k_retention**2) / (96 * (1 + k_retention)**2)
    C_m = f_k * d_c**2 / D_mobile

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
    D_mobile: float,
    mean_desorption_path_m: float,
    tortuosity: float = 1.41
) -> Dict:
    """
    Czas równowagi desorpcji.

    Model dwuetapowy:
    1) Dyfuzja przez ściankę polimeru (rate-limiting):
       Rozwiązanie dla płaskiej płytki grubości L = wall_thickness
       Ref: Crank (1975), "Mathematics of Diffusion", Chapter 4

       Czas do x% desorpcji z płytki:
         t_x = (L/2)² × α_x / D_polymer

       Gdzie α_x:
         α_50  = 0.049 (50% desorbed)
         α_90  = 0.66  (90% desorbed)
         α_95  = 0.90  (95% desorbed)
         α_99  = 1.27  (99% desorbed)
         α_99.9 = 1.88 (99.9% desorbed)

       Wyprowadzenie: z serii Fouriera Mt/M∞ = 1 - Σ 8/((2n+1)²π²) × exp(-(2n+1)²π²Dt/L²)

    2) Transport kanałami do krawędzi kontenera:
         t_channel = path² / (2 × D_eff)
         D_eff = D_mobile / τ²

    Args:
        wall_thickness_m: Grubość ścianki [m]
        D_polymer: Współczynnik dyfuzji w polimerze [m²/s]
        D_mobile: Współczynnik dyfuzji w fazie ruchomej [m²/s]
        mean_desorption_path_m: Średnia ścieżka desorpcji [m]
        tortuosity: Krętość [-]

    Returns:
        Dict z czasami desorpcji
    """
    # Half-thickness (dyfuzja z obu stron ścianki)
    L_half = wall_thickness_m / 2.0

    # Współczynniki Cranka (z serii Fouriera, plane sheet)
    alpha = {
        '50%': 0.049,
        '90%': 0.660,
        '95%': 0.897,
        '99%': 1.270,
        '99.9%': 1.878,
    }

    # Czasy dyfuzji przez ściankę
    t_wall = {}
    for pct, a in alpha.items():
        t_wall[pct] = L_half**2 * a / D_polymer

    # Transport kanałami
    D_eff_channel = D_mobile / tortuosity**2
    t_channel = mean_desorption_path_m**2 / (2 * D_eff_channel)

    # Total (suma — bo procesy są sekwencyjne)
    t_total = {}
    for pct in alpha:
        t_total[pct] = t_wall[pct] + t_channel

    # Rate-limiting step
    t_wall_90 = t_wall['90%']
    ratio = t_wall_90 / t_channel if t_channel > 0 else float('inf')

    if ratio > 10:
        rate_limiting = 'Wall diffusion (strongly dominant)'
    elif ratio > 2:
        rate_limiting = 'Wall diffusion (dominant)'
    elif ratio > 0.5:
        rate_limiting = 'Mixed (wall ≈ channel)'
    elif ratio > 0.1:
        rate_limiting = 'Channel transport (dominant)'
    else:
        rate_limiting = 'Channel transport (strongly dominant)'

    results = {
        't_wall_50_s': t_wall['50%'],
        't_wall_90_s': t_wall['90%'],
        't_wall_95_s': t_wall['95%'],
        't_wall_99_s': t_wall['99%'],
        't_channel_s': t_channel,
        't_total_90_s': t_total['90%'],
        't_total_95_s': t_total['95%'],
        't_total_99_s': t_total['99%'],
        'rate_limiting': rate_limiting,
        't_wall_to_channel_ratio': ratio,
        'D_eff_channel': D_eff_channel,

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
