"""
Predictions Module — Physicochemical Transport Predictions
=============================================================
Computational module for:
  - Backpressure ΔP (Darcy-based)
  - Van Deemter / HETP (band broadening)
  - Desorption equilibrium time (Crank diffusion + Convective Sweep)
  - Binding capacity

Formulas and constants from:
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
    Calculate backpressure using a generalized monolith equation (Darcy-based).

    ΔP = (K_geom * viscosity * u_superficial * L) / (d_h² * epsilon)

    Where:
      K_geom: Geometry resistance factor (default 32 for capillaries,
              typically 30-50 for TPMS/Gyroids).
      d_h: Hydraulic diameter = 4 * epsilon / S_volumetric

    Args:
        porosity: Porosity (fraction 0-1)
        hydraulic_diameter_m: Hydraulic diameter [m]
        column_length_m: Column/Bed length [m]
        viscosity_Pa_s: Mobile phase viscosity [Pa·s]
        flow_rate_uL_min: Flow rate [µL/min]
        column_diameter_m: Column diameter [m]
        k_geom: Geometry Resistance Factor

    Returns:
        Dict with pressure drop in various units and intermediate velocities
    """
    if porosity <= 0 or porosity >= 1:
        return {'error': 'Porosity must be between 0 and 1'}
    if hydraulic_diameter_m <= 0:
        return {'error': 'Hydraulic diameter must be positive'}

    # Convert flow rate: µL/min -> m³/s
    flow_rate_m3_s = flow_rate_uL_min * 1e-9 / 60.0

    # Column cross-sectional area [m²]
    column_area_m2 = np.pi * (column_diameter_m / 2) ** 2

    # Superficial velocity (u) [m/s]
    superficial_velocity_m_s = flow_rate_m3_s / column_area_m2

    # Interstitial velocity (u_e) [m/s]
    interstitial_velocity_m_s = superficial_velocity_m_s / porosity

    # Equation: ΔP = (K_geom * eta * u * L) / (d_h² * porosity)
    pressure_drop_pascal = (k_geom * viscosity_Pa_s * superficial_velocity_m_s * column_length_m) / \
                           (hydraulic_diameter_m ** 2 * porosity)

    # Unit conversions
    pressure_drop_bar = pressure_drop_pascal / 1e5
    pressure_drop_mbar = pressure_drop_pascal / 100.0
    pressure_drop_megapascal = pressure_drop_pascal / 1e6
    pressure_drop_psi = pressure_drop_pascal / 6894.76

    # Permeability (Darcy) [m²]
    if pressure_drop_pascal > 0:
        permeability_m2 = viscosity_Pa_s * superficial_velocity_m_s * column_length_m / pressure_drop_pascal
    else:
        permeability_m2 = None

    results = {
        'pressure_drop_pascal': pressure_drop_pascal,
        'pressure_drop_bar': pressure_drop_bar,
        'pressure_drop_mbar': pressure_drop_mbar,
        'pressure_drop_megapascal': pressure_drop_megapascal,
        'pressure_drop_psi': pressure_drop_psi,
        'superficial_velocity_m_s': superficial_velocity_m_s,
        'interstitial_velocity_m_s': interstitial_velocity_m_s,
        'permeability_m2': permeability_m2,
        'k_geom': k_geom,
        'flow_rate_uL_min': flow_rate_uL_min,
        'viscosity_Pa_s': viscosity_Pa_s,
    }

    # Warnings for high pressure
    warnings = []
    if pressure_drop_bar > 400:
        warnings.append('ΔP > 400 bar — exceeds typical HPLC limit!')
    if pressure_drop_bar > 1000:
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
    Generate a pressure drop vs flow rate curve.

    Returns:
        Dict with flow_rates and pressure_drop_bar arrays
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
            pressures.append(result['pressure_drop_bar'])

    return {
        'flow_rates_uL_min': flows.tolist(),
        'pressure_drop_bar': pressures,
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
    Calculate Van Deemter curve: H = A + B/u + (C_s + C_m) * u

    Model for TPMS structures:
      A = 2 * lambda * d_p       (eddy diffusion)
      B = 2 * gamma * D_m / tau² (longitudinal diffusion)
      C_s = (2/3) * k/(1+k)² * d_f² / D_s (stationary phase mass transfer)
      C_m = omega * d_h² / D_m   (mobile phase mass transfer)

    Parameters for ordered TPMS:
      - lambda ≈ 0.1 (ordered structure, extremely low eddy dispersion)
      - gamma ≈ 0.6 (obstruction factor)
      - omega ≈ 0.1 (structural factor for monoliths)

    Args:
        unit_cell_m: Unit cell size [m]
        wall_thickness_m: Wall thickness [m]
        hydraulic_diameter_m: Hydraulic diameter of channels [m]
        porosity: Porosity [fraction]
        D_mobile: Diffusion coefficient in mobile phase [m²/s]
        D_polymer: Diffusion coefficient in polymer [m²/s]
        tortuosity: Tortuosity [-]
        u_range: Linear velocity range [m/s]. None = auto.
        n_points: Number of points on curve

    Returns:
        Dict with Van Deemter curve and optimal parameters
    """
    # Model parameters
    lambda_eddy = 0.1    # Eddy diffusion parameter (regular TPMS structure)
    gamma_obstr = 0.6    # Obstruction factor for longitudinal diffusion
    omega_mobile = 0.1   # Structural factor for monoliths (mobile phase)

    # Effective "particle" size for gyroid ≈ unit cell
    d_p = unit_cell_m

    # Stationary phase film thickness ≈ half of wall thickness
    d_f = wall_thickness_m / 2.0

    # Hydraulic diameter of channels
    d_h = hydraulic_diameter_m

    # Retention factor (k)
    k_retention = 2.0

    # ---- Van Deemter Terms ----

    # A: Eddy diffusion
    A = 2 * lambda_eddy * d_p

    # B: Longitudinal (axial) diffusion
    axial_diffusion_coeff = (gamma_obstr / tortuosity**2) * D_mobile
    B = 2 * axial_diffusion_coeff

    # C_s: Mass transfer resistance in stationary phase
    C_s = (2.0/3.0) * (k_retention / (1 + k_retention)**2) * d_f**2 / D_polymer

    # C_m: Mass transfer resistance in mobile phase
    C_m = omega_mobile * (d_h**2 / D_mobile)

    # C total
    C_total = C_s + C_m

    # ---- Automatic Velocity Range ----
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

    # Ensure physically sensible range
    u_min = max(u_min, 1e-7)
    u_max = min(u_max, 1.0)

    # ---- Calculate Curve ----
    u_values = np.linspace(u_min, u_max, n_points)

    H_A = np.full_like(u_values, A)
    H_B = B / u_values
    H_Cs = C_s * u_values
    H_Cm = C_m * u_values
    H_total = H_A + H_B + H_Cs + H_Cm

    # ---- Optimum ----
    if C_total > 0:
        optimal_velocity_m_s = np.sqrt(B / C_total)
        minimal_hetp_m = A + 2 * np.sqrt(B * C_total)
    else:
        optimal_velocity_m_s = u_values[np.argmin(H_total)]
        minimal_hetp_m = np.min(H_total)

    # Number of theoretical plates per meter
    if minimal_hetp_m > 0:
        theoretical_plates_per_meter = 1.0 / minimal_hetp_m
    else:
        theoretical_plates_per_meter = None

    results = {
        # Term parameters
        'eddy_diffusion_term_m': A,
        'longitudinal_diffusion_term_m2_s': B,
        'stationary_phase_mass_transfer_s': C_s,
        'mobile_phase_mass_transfer_s': C_m,
        'total_mass_transfer_s': C_total,

        # Optimum
        'optimal_velocity_m_s': optimal_velocity_m_s,
        'minimal_hetp_m': minimal_hetp_m,
        'minimal_hetp_um': minimal_hetp_m * 1e6,
        'theoretical_plates_per_meter': theoretical_plates_per_meter,

        # Curve data
        'u_values_m_s': u_values.tolist(),
        'H_total_m': H_total.tolist(),
        'H_A_m': H_A.tolist(),
        'H_B_m': H_B.tolist(),
        'H_Cs_m': H_Cs.tolist(),
        'H_Cm_m': H_Cm.tolist(),

        # Input parameters
        'D_mobile': D_mobile,
        'D_polymer': D_polymer,
        'k_retention': k_retention,
        'lambda_eddy': lambda_eddy,
        'gamma_obstruction': gamma_obstr,
        'tortuosity': tortuosity,

        # Analysis
        'dominant_term': 'C_s (stationary)' if C_s > C_m else 'C_m (mobile)',
        'stationary_phase_fraction': C_s / C_total if C_total > 0 else 0,
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
    Calculate desorption equilibrium time (SPE model).

    Two-stage model:
    1) Diffusion through the polymer wall (often rate-limiting):
       t_x = (wall/2)² * alpha_x / D_polymer
       Where alpha_x is the Fourier number (Fo) for desorption from a slab of half-thickness.

    2) Convective sweep of the channels:
       t_sweep = Void_Volume / Flow_Rate

    Args:
        wall_thickness_m: Wall thickness [m]
        D_polymer: Diffusion coefficient in polymer [m²/s]
        void_volume_mm3: Void (channel) volume [mm³]
        flow_rate_uL_min: Flow rate [µL/min]

    Returns:
        Dict with desorption times and rate-limiting analysis
    """
    # 1. Wall diffusion (Crank, half-thickness slab model)
    L_half = wall_thickness_m / 2.0

    # Updated alpha coefficients (Fo) for half-thickness (l)
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
    Calculate analyte binding capacity.

    Model:
      Q_vol = ASA * q_surface * accessibility / V_total
      Q_mass = Q_vol * V_total / mass

    Where:
      ASA = internal_SA * probe_accessible_fraction
      q_surface = surface binding density [mg/m²]
        Typical values:
          MIP (imprinted polymer): 1-5 mg/m²
          C18 silica: 2-4 mg/m²
          Ion exchange: 0.5-2 mg/m²
      accessibility = fraction of surface physically accessible (0.5-1.0)
      probe_accessible_fraction = from ASA analysis (probe-based)

    Args:
        internal_surface_area_mm2: Internal surface area [mm²]
        total_volume_mm3: Container total volume [mm³]
        solid_volume_mm3: Solid phase volume [mm³]
        polymer_density_g_cm3: Material density [g/cm³]
        q_surface_mg_m2: Surface capacity [mg/m²]
        accessibility: Global accessibility factor [0-1]
        probe_accessible_fraction: Size-exclusion fraction from ASA [0-1]

    Returns:
        Dict with capacities in various units
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
    Calculate hydraulic diameter of the channels.

    d_h = 4 * epsilon * V_total / S_internal

    Ref: Standard hydraulic diameter definition for porous media.

    Args:
        porosity: Porosity [0-1]
        total_volume_mm3: Container volume [mm³]
        internal_surface_area_mm2: Internal surface area [mm²]

    Returns:
        d_h in mm
    """
    if internal_surface_area_mm2 <= 0:
        return 0

    void_volume = porosity * total_volume_mm3
    d_h = 4.0 * void_volume / internal_surface_area_mm2

    return d_h
