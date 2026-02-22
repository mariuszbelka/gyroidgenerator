# Gyroid Generator v3.3 — DLP Sorbent Design Platform

A specialized design and analysis tool for Triply Periodic Minimal Surface (TPMS) structures, specifically optimized for DLP 3D-printed sorbents in Liquid Chromatography (LC) and Solid Phase Extraction (SPE).

## Quick Start

```cmd
# Run the application
python gui_main.py
```

Requirements:
```cmd
pip install -r requirements.txt
```

---

## Key Features (v3.3)

### 🔬 High-Precision Geometric Analysis
Unlike standard voxel-based tools, Gyroid Generator v3.3 uses **Mesh Ray-Casting** and **Masked Voxelization** for superior accuracy:
- **Ray-Casting Distributions**: Wall and channel thicknesses are measured by casting thousands of rays against the exact mesh triangulation, avoiding the pixelation errors of Euclidean Distance Transforms (EDT).
- **Skin Filtering**: Automatically removes container boundary artifacts to ensure distributions reflect the internal gyroid structure.
- **Container Masking**: All topological metrics (connectivity, ASA) are restricted to the actual design volume (sphere, cylinder, etc.), ensuring physically consistent porosity values.

### 🧪 Physicochemical Predictions
The platform translates geometric data into chromatographic performance indicators using established physical models:

#### 1. Darcy-Based Pressure Drop (ΔP)
Uses a generalized monolith equation to predict backpressure:
`ΔP = (K_geom * η * u * L) / (d_h² * ε)`
- **K_geom**: Geometry resistance factor (default 32 for capillaries, ~35-50 for gyroids).
- **η**: Viscosity (presets for common LC solvents).
- **d_h**: Hydraulic diameter derived from internal wetted area.

#### 2. Van Deemter / HETP Analysis
Predicts band broadening (Efficiency) for TPMS structures:
`H = A + B/u + (C_s + C_m) * u`
- **A (Eddy Diffusion)**: λ=0.1 (extremely low due to ordered TPMS geometry).
- **B (Longitudinal)**: Includes obstruction factors for porous media.
- **C (Mass Transfer)**: Separates stationary phase diffusion (C_s) and mobile phase resistance (C_m).

#### 3. Desorption Equilibrium Time
A two-stage kinetics model for SPE optimization:
- **Wall Diffusion**: Based on Crank's slab model using half-thickness.
- **Convective Sweep**: Calculation of the time required to flush the void volume at a given flow rate.

#### 4. Binding Capacity
Estimates total analyte load based on **Accessible Surface Area (ASA)**:
`Q = ASA * Surface_Density * Accessibility_Factor`

---

## GUI Tabs

- **3D View**: Interactive OpenGL viewer (Top/Front/Side views).
- **Statistics**: Instant mesh properties, external/internal surface area split, and wall diffusion path distribution.
- **Predictions**: Performance charts (ΔP vs Flow, Van Deemter curve, Desorption breakdown).
- **Distributions**: Wall/Channel thickness histograms and printability scoring (DLP-specific).
- **Cross-Section**: 2D slice visualization and SVG/PNG export.
- **Compare**: Side-by-side comparison of up to 6 designs with best/worst metric highlighting.

---

## File Architecture

```
gui_main.py              ← Main entry point and tab integration.
gyroid_math.py           ← Core TPMS equations and resolution safety logic.
mesh_generator.py        ← Isosurface extraction (Marching Cubes) and container logic.
statistics_analyzer_v2.py← Voxel/Mesh-based geometric metrics engine.
predictions.py           ← Physics engine (Darcy, Van Deemter, Crank models).
distributions.py         ← Topological engine (Ray-casting, Connectivity, ASA).
viewer_3d.py             ← High-performance OpenGL mesh viewer.
```

---

## Scientific Transparency

Gyroid Generator is an **open tool**. All calculations are based on published chromatographic theory. Users are encouraged to adjust constants (K_geom, λ, ω) in the GUI to align the models with their experimental findings.

**Version:** 3.3 (Reviewer Edition)
**Author:** Jules (Anthropic) for DLP Research
**Date:** 2026-02-20
