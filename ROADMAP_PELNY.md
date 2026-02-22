# Gyroid Generator — PEŁNY ROADMAP ROZWOJU

**Data:** 2026-02-20
**Źródło:** Specyfikacja v3 + analiza fizyka/chemika (konsultant)
**Cel:** Platforma do projektowania, predykcji i walidacji 3D-printed sorbentów/faz LC

---

## FILOZOFIA PROJEKTU

### Trzy poziomy modelu (kluczowa decyzja architektoniczna)

| Poziom | Co reprezentuje | Źródło danych |
|--------|----------------|---------------|
| **Nominal (CAD/STL)** | To co projektujesz | Generator gyroidu |
| **As-printed (predicted)** | To co się wydrukuje | Model korekcji (Jacobs + empiryka) |
| **As-measured** | To co zmierzono | SEM, micro-CT, masa, ΔP |

Wszystkie metryki powinny być liczone osobno dla każdego poziomu.
To pozwala rozdzielić: problem projektu vs problem druku vs problem chemii.

### Dwa profile zastosowania

| Profil | Priorytetowe KPI |
|--------|-----------------|
| **SPE / batch extraction** | t90 desorpcji, pojemność, dead-end fraction, accessibility |
| **LC kolumna** | ΔP, HETP, dyspersja, jednorodność kanałów |

---

## ETAP 1 — v3.0: "Geometry + Printability Intelligence"
**Status: DO IMPLEMENTACJI TERAZ**
**Szacowany czas: 1-2 sesje kodowania**

### 1.1 Fix: 3D Viewer (pyqtgraph.opengl)
- Wbudowany viewer PyQt6 zamiast pyglet
- Obracanie, zoom, przesuwanie
- **Trudność: NISKA** — pyqtgraph już zainstalowany

### 1.2 Rozkłady grubości ścian i kanałów
- Voxelizacja → Distance Transform EDT
- Histogram wall thickness (z local maxima w solid)
- Histogram channel width (z local maxima w void)
- Percentyle: P10, P50, P90
- Minima lokalne (ryzyko niedodruku / zatykania)
- **Trudność: NISKA** — mamy już voxelizację z v2

### 1.3 Connectivity / Dead-end fraction
- Connected components analysis na void voxelach
- % void connected inlet→outlet (through-pores)
- % dead-end void volume (martwe strefy)
- Liczba oddzielnych komponentów
- **Trudność: NISKA** — scipy.ndimage.label na voxel grid

### 1.4 Accessible Surface Area (ASA) — probe-based
- Erozja void space o promień sondy
- Zmierz surface area po erozji
- Profile dla kilku promieni (1nm, 2nm, 5nm, 10nm)
- ASA/total SA ratio = realny współczynnik dostępności
- **Trudność: NISKA-ŚREDNIA** — erosion + marching cubes

### 1.5 Throat / constriction analysis
- Local minima w distance transform void space
- Rozkład throat sizes (P10, P50, P90)
- Throat-to-pore ratio
- Lokalizacja krytycznych przewężeń
- **Trudność: ŚREDNIA** — wymaga medial axis lub local minima detection

### 1.6 Printability score
- Min wall thickness vs rozdzielczość drukarki (22 µm pixel)
- Min channel width vs ryzyko bleed/zalania
- Lokalizacja krytycznych miejsc (heatmapa)
- Trapped resin risk (pore drainage analysis)
- Score: PASS / WARNING / FAIL z opisem
- **Trudność: NISKA-ŚREDNIA** — z danych dystrybucji

### 1.7 Predykcje transportowe (proxy, półempiryczne)

#### 1.7a Ciśnienie wsteczne ΔP (Kozeny-Carman)
- ΔP = f(η, u, L, d_h, ε)
- Presety lepkości (woda, ACN/H₂O, MeOH)
- Wykres ΔP vs przepływ
- **Trudność: NISKA** — czysty wzór analityczny

#### 1.7b Van Deemter / HETP
- H = A + B/u + Cs·u + Cm·u
- Presety analitów + manual D_mobile
- D_polymer jako input
- Wykres H vs u z rozbiciem na składowe
- **Trudność: NISKA** — wzory analityczne + wykresy

#### 1.7c Czas równowagi desorpcji (t90, t99)
- Model: dyfuzja z płaskiej płytki (ścianka)
- t_wall + t_channel = t_total
- Rate-limiting step identification
- **Trudność: NISKA** — wzory analityczne

#### 1.7d Pojemność wiązania
- Q = ASA × q_surface × accessibility
- Presety + manual override
- **Trudność: NISKA**

### 1.8 Wizualizacja przekrojów + eksport
- Slice XY/XZ/YZ na wybranej pozycji
- Solid/void obraz
- Eksport PNG/SVG
- Lokalna porowatość z przekroju
- **Trudność: ŚREDNIA**

### 1.9 Porównywarka projektów (A/B/C)
- "Zapisz do porównania" po generacji
- Tabela porównawcza wszystkich metryk
- Highlight najlepsza/najgorsza wartość
- Eksport CSV
- **Trudność: ŚREDNIA**

### 1.10 Panel fizyczny w GUI
- Sekcja "Physics / Material" w lewym panelu
- Dropdown: Solvent (presety)
- Dropdown: Analyte (presety + manual)
- Spinbox: D_polymer
- Spinbox: Flow rate
- Spinbox: q_surface, accessibility
- **Trudność: NISKA**

---

## ETAP 2 — v4.0: "Transport Intelligence"
**Status: PO WALIDACJI v3**

### 2.1 Tortuosity — prawdziwa (nie stała 1.41)
- Geodesic tortuosity na voxel grid (pathfinding A*)
- Hydraulic tortuosity proxy
- Rozkład długości ścieżek
- Osobno dla osi X, Y, Z (anizotropia)
- **Trudność: ŚREDNIA-WYSOKA**

### 2.2 Anizotropia struktury
- Tensor anizotropii z geometrii
- Metryki osobno dla osi przepływu vs prostopadłych
- Wykrywanie preferowanego kierunku
- **Trudność: ŚREDNIA**

### 2.3 RTD proxy (Residence Time Distribution)
- Uproszczony rozkład czasu przejścia z grafu kanałów
- Wskaźnik szerokości RTD
- Udział szybkich/wolnych ścieżek
- **Trudność: ŚREDNIA-WYSOKA**

### 2.4 Axial dispersion proxy
- Péclet osiowy
- Wskaźnik dyspersji z geometrii + tortuosity + przewężeń
- **Trudność: ŚREDNIA**

### 2.5 Permeability tensor K
- Kx, Ky, Kz osobno
- Kalibrowalny model półempiryczny
- **Trudność: ŚREDNIA**

### 2.6 Mass-transfer proxy (a·kf)
- Interfacial area per volume (mamy)
- kf z korelacji Sherwood-Reynolds-Schmidt
- **Trudność: NISKA-ŚREDNIA**

### 2.7 Gradienty — podłączenie do GUI
- Linear Z gradient (wall thickness i/lub unit cell)
- Concentryczny gradient
- GradientGyroid już istnieje w gyroid_math.py
- **Trudność: ŚREDNIA** (GUI + mesh_generator integration)

### 2.8 Batch generowanie (sweep 2D)
- Sweep unit_cell × wall_thickness
- Tabela wyników + heatmapa
- Eksport CSV + opcjonalnie STL
- **Trudność: ŚREDNIA**

---

## ETAP 3 — v5.0: "As-Printed Model"
**Status: PO KALIBRACJI Z EKSPERYMENTEM**

### 3.1 Predykcja geometrii po druku
- Input: parametry druku (ekspozycja, Dp, Ec, bleed proxy)
- Model korekcji:
  - Powiększenie ścianek (overcure/bleed)
  - Zwężenie kanałów
  - Shrinkage
- Output: skorygowana geometria
- **Trudność: WYSOKA** — wymaga danych eksperymentalnych do kalibracji

### 3.2 Przeliczanie metryk na as-printed geometry
- Wszystkie metryki z etapu 1-2 przeliczone na skorygowanej geometrii
- Porównanie nominal vs as-printed
- **Trudność: ŚREDNIA** (po implementacji 3.1)

### 3.3 Printability checks v2 — z modelem druku
- Risk map na podstawie modelu druku
- Predykcja gdzie kanały się zaleją
- Predykcja gdzie ścianki nie wydrukują
- **Trudność: WYSOKA**

---

## ETAP 4 — v6.0: "Validation + Calibration"
**Status: GDY MASZ DANE EKSPERYMENTALNE**

### 4.1 Import danych eksperymentalnych
- Format: CSV / Excel
- Typy danych: masa, ΔP, porowatość (BET), SEM wymiary, t90 batch, recovery
- **Trudność: NISKA**

### 4.2 Dashboard predicted vs measured
- Wykresy: predicted vs measured dla każdej metryki
- Residual plots
- R², RMSE, bias
- **Trudność: ŚREDNIA**

### 4.3 Korekcja modelu (bias correction)
- Fit: linear correction (predicted × a + b = measured)
- Osobno dla każdej metryki
- Auto-update predykcji
- **Trudność: ŚREDNIA**

### 4.4 Baza danych projektów (SQLite)
- Project → Design → Print → Predicted → Measured
- Porównania między seriami
- Reproducibility tracking
- Export do CSV/Excel
- **Trudność: ŚREDNIA**

---

## ETAP 5 — v7.0+: "High Fidelity"
**Status: PRZYSZŁOŚĆ / OPCJONALNY**

### 5.1 Solver przepływu (Stokes / LBM)
- Na voxel grid
- Pole prędkości, ciśnienia, shear
- Prawdziwy ΔP (nie proxy)
- **Trudność: BARDZO WYSOKA**

### 5.2 Advection-diffusion tracer
- Symulacja transportu znacznika
- Prawdziwy RTD
- Breakthrough curve prediction
- **Trudność: BARDZO WYSOKA**

### 5.3 Pore Network Model (PoreSpy / OpenPNM)
- Ekstrakcja sieci porów
- Pełna analiza connectivity, permeability, drainage
- **Trudność: WYSOKA** (ale biblioteki istnieją)

### 5.4 Import STL jako kontener
- Dowolny kształt z CAD → wypełnienie gyroidem
- mesh.contains() na importowanym STL
- **Trudność: WYSOKA** (edge cases, performance)

---

## TOP 10 METRYK WG KOLEGI-FIZYKA (priorytet implementacji)

| # | Metryka | Etap | Już w v3? |
|---|---------|------|-----------|
| 1 | Porosity (connected only) | 1 | ✅ TAK |
| 2 | Accessible Surface Area (probe-based) | 1 | ✅ TAK |
| 3 | Throat size distribution (P10/P50/P90) | 1 | ✅ TAK |
| 4 | Dead-end void fraction | 1 | ✅ TAK |
| 5 | Hydraulic tortuosity (proxy) | 2 | ❌ v4 |
| 6 | Permeability / ΔP predictor | 1 | ✅ TAK (proxy) |
| 7 | Axial dispersion / RTD proxy | 2 | ❌ v4 |
| 8 | Printability risk score | 1 | ✅ TAK |
| 9 | As-printed bias predictor | 3 | ❌ v5 |
| 10 | Batch t90 proxy / LC mass-transfer score | 1 | ✅ TAK |

**7 z 10 wchodzi do v3.** To jest bardzo dobry wynik.

---

## ARCHITEKTURA PLIKÓW v3

```
gyroid_generator/
│
├── gyroid_math.py              ← BEZ ZMIAN
├── mesh_generator.py           ← BEZ ZMIAN
├── statistics_analyzer_v2.py   ← BEZ ZMIAN
├── statistics_tab_v2.py        ← BEZ ZMIAN
│
├── gui_main.py                 ← ROZSZERZONY (nowe zakładki + panel fizyczny)
│
├── viewer_3d.py                ← NOWY (pyqtgraph OpenGL viewer)
├── predictions.py              ← NOWY (ΔP, Van Deemter, t90, capacity)
├── predictions_tab.py          ← NOWY (GUI dla predykcji)
├── distributions.py            ← NOWY (histogramy, throat, connectivity)
├── distributions_tab.py        ← NOWY (GUI dla rozkładów)
├── cross_section.py            ← NOWY (wizualizacja przekrojów)
├── cross_section_tab.py        ← NOWY (GUI dla przekrojów)
├── printability.py             ← NOWY (score drukowalności)
├── compare_tab.py              ← NOWY (porównywarka A/B/C)
│
├── requirements.txt            ← UPDATED (+ ewentualne nowe zależności)
└── test_v3.py                  ← NOWY (testy bez GUI)
```

**Zasada:** stare pliki nie są modyfikowane (poza gui_main.py który dostaje nowe zakładki).
Każda funkcja to osobny moduł .py + osobna zakładka GUI.

---

## LAYOUT GUI v3

```
┌───────────────────────────────────────────────────────────────────┐
│ Gyroid Generator v3.0 — DLP Sorbent Design Platform              │
├────────────────┬──────────────────────────────────────────────────┤
│ ── GEOMETRY ── │  TABS:                                           │
│ Shape ▼        │  [3D View] [Statistics] [Predictions]            │
│ Dimensions     │  [Distributions] [Cross-Section] [Compare]      │
│ Gyroid Params  │                                                  │
│ Quality ▼      │                                                  │
│                │  ┌─────────────────────────────────────────────┐ │
│ [Generate]     │  │                                             │ │
│ [Export STL]   │  │         (aktywna zakładka)                  │ │
│                │  │                                             │ │
│ ── PHYSICS ──  │  │                                             │ │
│ Solvent ▼      │  │                                             │ │
│ Analyte ▼      │  │                                             │ │
│ D_polymer      │  │                                             │ │
│ Flow rate      │  └─────────────────────────────────────────────┘ │
│ q_surface      │                                                  │
│                │                                                  │
│ ── TOOLS ──    │                                                  │
│ [Save to       │                                                  │
│  Compare]      │                                                  │
├────────────────┴──────────────────────────────────────────────────┤
│ Status: Ready                                                     │
└───────────────────────────────────────────────────────────────────┘
```

---

## CO v3 DAJE NAUKOWO

Z samego v3 (bez etapów 2-5) masz:

1. **Screening projektów przed drukiem** — printability score eliminuje marnotrawstwo żywicy
2. **Predykcja ΔP** — pierwsze sito dla LC candidates
3. **Van Deemter** — ocena separacyjna przed drukiem
4. **Dead-end / connectivity** — wykrywanie "fałszywego SSA"
5. **ASA probe-based** — realistyczna pojemność zamiast total SA
6. **Throat distribution** — wykrywanie słabych punktów geometrii
7. **Porównywarka** — systematyczny wybór najlepszego projektu
8. **Cross-sections** — szybka wizualna diagnostyka

To wystarczy na solidną publikację metodyczną: "Computational design framework for 3D-printed gyroid sorbents".

---
