# Gyroid Generator — FULL DEVELOPMENT ROADMAP

**Date:** 2026-02-20
**Target:** Platform for design, prediction, and validation of 3D-printed sorbents/LC phases.

---

## DESIGN PHILOSOPHY

### Three Model Levels (Architectural Decision)

| Level | Representation | Data Source |
|--------|----------------|---------------|
| **Nominal (CAD/STL)** | What you design | Gyroid Generator v3.3 |
| **As-printed (predicted)** | What will likely print | Correction Models (v5.0) |
| **As-measured** | Reality | SEM, micro-CT, Mass, ΔP |

Metrics are calculated separately for each level. This decouples design issues from printing issues and chemical surface issues.

### Two Application Profiles

| Profile | Priority KPIs |
|--------|-----------------|
| **SPE / Batch Extraction** | Desorption t90, Capacity, Dead-end fraction, Accessibility |
| **LC Column** | ΔP, HETP (Van Deemter), Dispersion, Channel uniformity |

---

## CURRENT STATUS — v3.3: "Precision & Physics Intelligence"
**Status: COMPLETED**
**Key Achievements:**
- **Mesh-based Ray Casting**: High-precision wall/channel thickness distributions.
- **Scientific Fidelity**: Darcy, Van Deemter, and Crank kinetics models updated and verified.
- **Container Masking**: Physically consistent porosity and connectivity metrics (restricted to object volume).
- **Skin Filtering**: Removal of container boundary artifacts from distributions.
- **STL Import**: Support for analyzing external files.
- **Memory Safety**: Optimized masking and resolution capping for high-density meshes.

---

## PHASE 2 — v4.0: "Transport Intelligence"
**Target: AFTER v3 VALIDATION**

### 2.1 True Tortuosity (Beyond constant 1.41)
- Geodesic tortuosity on voxel grid (pathfinding A*).
- Hydraulic tortuosity proxy.
- Path length distribution (anisotropy X, Y, Z).

### 2.2 Structural Anisotropy
- Anisotropy tensor from geometry.
- Metrics for flow-axis vs. perpendicular axes.

### 2.3 RTD Proxy (Residence Time Distribution)
- Simplified transit time distribution from channel graph.
- Fast/slow path contribution analysis.

### 2.4 Axial Dispersion Proxy
- Geometric Péclet number estimate.
- Dispersion index based on tortuosity and constrictions.

### 2.5 Gradient UI Integration
- Direct GUI control for Linear Z and Concentric gradients (already in math engine).

---

## PHASE 3 — v5.0: "As-Printed Model"
**Target: AFTER EXPERIMENTAL CALIBRATION**

### 3.1 Post-Print Geometry Prediction
- Inputs: Exposure, Dp, Ec, Bleed proxy.
- Model: Overcure/Bleed compensation, narrowing of channels, shrinkage.

### 3.2 Post-Print Metric Recalculation
- Nominal vs. As-Printed comparison for all v3.3 metrics.

---

## PHASE 4 — v6.0: "Validation + Calibration"
**Target: WITH EXPERIMENTAL DATASETS**

### 4.1 Experimental Data Import
- BET Porosity, Mass, measured ΔP, SEM dimensions.

### 4.2 Dashboard: Predicted vs. Measured
- Residual plots, R², RMSE, Bias correction.

### 4.3 Bias Correction (Calibration)
- Linear correction fits for all physical models.

---

## PHASE 5 — v7.0+: "High Fidelity"
**Target: FUTURE / RESEARCH**

### 5.1 Flow Solver (Stokes / LBM)
- Direct CFD simulation on the voxel grid.
- Shear stress and true ΔP calculation.

### 5.2 Pore Network Model (PoreSpy Integration)
- Extraction of pore networks for advanced drainage analysis.

---

## TOP 10 METRICS (Implementation Priority)

| # | Metric | Status | Implementation Level |
|---|---------|------|-----------|
| 1 | Porosity (Connected Only) | ✅ | High (Masked Voxel) |
| 2 | Accessible Surface Area (ASA) | ✅ | High (Probe-based) |
| 3 | Throat Size Distribution | ✅ | High (Analytical-Scale) |
| 4 | Dead-end Void Fraction | ✅ | High (Masked Voxel) |
| 5 | Backpressure (ΔP) | ✅ | High (Darcy-based) |
| 6 | Van Deemter (HETP) | ✅ | High (TPMS-optimized) |
| 7 | Desorption Kinetics (t90) | ✅ | High (Crank Model) |
| 8 | Printability Score | ✅ | Medium (Percentile-based) |
| 9 | True Tortuosity | ❌ | Planned (v4.0) |
| 10 | Flow Anisotropy | ❌ | Planned (v4.0) |

---

## SCIENTIFIC IMPACT OF v3.3

The current version provides a complete **Methodological Publication Framework**:
1. **Pre-print Screening**: Eliminating resin waste via printability scores.
2. **LC Candidate Filtering**: Early ΔP and HETP selection.
3. **Connectivity Verification**: Detecting "False SSA" in dead-end zones.
4. **Realistic Capacity**: Size-exclusion (ASA) based load estimation.

This framework is sufficient for a high-impact paper: *"Computational design and performance prediction of 3D-printed gyroid sorbents."*
