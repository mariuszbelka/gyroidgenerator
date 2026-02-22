# Gyroid Generator v3.0 — DLP Sorbent Design Platform

## Uruchomienie

```cmd
python gui_main.py
```

Wymagania (już zainstalowane z v2):
```cmd
pip install -r requirements.txt
```

---

## Nowe funkcje w v3 (w porównaniu z v2)

### Zakładka: 3D View
**Wbudowany viewer 3D** — zastępuje trimesh.show() (który wymagał pyglet).

- Obracanie: lewy przycisk myszy
- Zoom: scroll
- Przesuwanie: prawy przycisk myszy
- Przyciski widoku: Top / Front / Side / Reset
- Automatyczny decimation dla dużych meshów (>200k ścian)

**Technologia:** pyqtgraph OpenGL (już zainstalowany).

---

### Zakładka: Statistics (bez zmian z v2)
Wszystkie instant statistics i desorption path z v2 działają bez zmian.

---

### Zakładka: Predictions
Predykcje fizykochemiczne bazujące na wzorach analitycznych.

#### Panel fizyczny (wewnątrz zakładki Predictions)
Parametry materiałowe ustawia się bezpośrednio w zakładce Predictions (obok predykcji które ich używają):
- **Solvent** — dropdown z presetami lepkości:
  - Water 25°C: η = 0.890 mPa·s
  - ACN/H₂O 20:80 do 80:20
  - MeOH/H₂O 50:50, MeOH 100%
- **Analyte** — dropdown z presetami D_mobile + manual override:
  - Mała cząsteczka 300 Da: D = 1.0 × 10⁻⁹ m²/s
  - Mała cząsteczka 500 Da: D = 7.0 × 10⁻¹⁰ m²/s
  - Peptyd 2 kDa: D = 3.0 × 10⁻¹⁰ m²/s
  - Białko 15 kDa: D = 1.2 × 10⁻¹⁰ m²/s
  - Białko 66 kDa: D = 6.0 × 10⁻¹¹ m²/s
- **D_polymer** — współczynnik dyfuzji w polimerze [m²/s], domyślnie 10⁻¹¹
- **Flow rate** — przepływ [µL/min]

#### 1. Ciśnienie wsteczne ΔP
**Wzór:** Kozeny-Carman (1937)
```
ΔP = (150 × η × u × L) / (d_h² × ε³)
```
- η = lepkość [Pa·s] z presetu
- u = prędkość powierzchniowa [m/s] = F / A_column
- L = długość kolumny [m]
- d_h = średnica hydrauliczna = 4εV/S_internal [m]
- ε = porowatość z mesha

**Co pokazuje:**
- ΔP w barach i PSI
- Prędkość powierzchniowa i interstycjalna
- Permeability (Darcy) [m²]
- Wykres ΔP vs przepływ
- Ostrzeżenia: >400 bar (HPLC limit), >1000 bar (UHPLC limit)

**Stała K_KC = 150** — standard Blake-Kozeny dla packed beds i monolitów.

#### 2. Van Deemter / HETP
**Wzór:**
```
H = A + B/u + C_s·u + C_m·u
```

Termy:
- **A = 2λd_p** — dyspersja wirowa (multiple paths)
  - λ = 0.5 (regularna struktura TPMS, niższe niż packed bed λ~1.5)
  - d_p = unit_cell (efektywny rozmiar „ziarna")
- **B = 2γD_m/τ²** — dyfuzja podłużna (longitudinal diffusion)
  - γ = 0.6 (obstruction factor)
  - τ = tortuosity
- **C_s = (2/3) × k/(1+k)² × (d_f²/D_s)** — transfer masy w fazie stacjonarnej
  - d_f = wall_thickness/2 (połowa grubości ścianki)
  - D_s = D_polymer
  - k = 2.0 (domyślny retention factor)
- **C_m = f(k) × d_c²/D_m** — transfer masy w fazie ruchomej
  - f(k) = (1+6k+11k²)/(96(1+k)²) — równanie Golaya (1958)
  - d_c = channel_width

**Co pokazuje:**
- H_min [µm] i u_opt [mm/s]
- N/m (liczba półek na metr)
- Wykres H vs u z rozbiciem na A, B, C_s, C_m
- Który term dominuje (C_s czy C_m)
- Frakcja C_s w C_total

**Refs:** van Deemter (1956), Knox (1999), Golay (1958), Aris (1956)

#### 3. Czas równowagi desorpcji
**Model:** Dyfuzja z płaskiej płytki (Crank, 1975)

```
t_wall = (wall/2)² × α / D_polymer
t_channel = path² / (2 × D_mobile/τ²)
t_total = t_wall + t_channel
```

Współczynniki Cranka (z serii Fouriera):
- α_50% = 0.049
- α_90% = 0.660
- α_95% = 0.897
- α_99% = 1.270

**Co pokazuje:**
- t_wall (czas dyfuzji przez ściankę) dla 50%, 90%, 95%, 99%
- t_channel (czas transportu kanałami)
- t_total (suma)
- Rate-limiting step (co dominuje: ścianka czy kanał?)
- Wykres słupkowy: t_wall vs t_channel

#### 4. Pojemność wiązania
**Wzór:**
```
Q = ASA × q_surface × accessibility / V_total
```

- ASA = dostępna powierzchnia wewnętrzna [m²]
- q_surface = gęstość wiązania [mg/m²] (podajesz w GUI)
- accessibility = współczynnik dostępności (0-1)

**Co pokazuje:**
- Q_vol [mg/mL] — pojemność objętościowa
- Q_mass [mg/g] — pojemność masowa
- Total capacity [mg] — dla danego kontenera

---

### Zakładka: Distributions
Analiza rozkładów i topologii struktury. **Bazuje na voxelizacji + Distance Transform (EDT).**

Wszystkie obliczenia uruchamiane jednym kliknięciem — zaznacz checkboxy i kliknij „Calculate".
Voxelizacja jest cachowana (liczona raz, reużywana).

#### 1. Wall Thickness Distribution
**Metoda:** Hildebrand & Rüegsegger (1997)
- EDT na fazie solid → local maxima = promienie wpisanych sfer
- Wall thickness = 2 × radius w każdym punkcie
- Histogram + percentyle P10, P50, P90

**Co pokazuje:**
- Histogram grubości ścian [µm]
- Mean, Std, Min, Max, P10, P50, P90
- Uniformity score (1.0 = idealnie jednorodne)
- **Min wall** — krytyczna dla drukowalności

#### 2. Channel Width Distribution
**Metoda:** Analogiczna do wall thickness, ale na fazie void.
- EDT na void → local maxima → channel width = 2 × radius

**Co pokazuje:**
- Histogram szerokości kanałów [µm]
- Mean, Std, Min, Max, P10, P50, P90
- **Min channel** — ryzyko bleed/zalania żywicą

#### 3. Connectivity / Dead-End Analysis
**Metoda:**
- Connected components (scipy.ndimage.label, 26-connectivity)
- Sprawdzenie które void-komponenty łączą Z_min z Z_max (through-pores)
- Dead-end = void niepołączony z obu stronami

**Co pokazuje:**
- Liczba komponentów void
- Through-connected fraction (% void łączącego inlet z outlet)
- Dead-end fraction (% martwych stref)
- Isolated fraction (% zamkniętego void)
- Connected porosity (porowatość tylko z through-pores)
- **Ostrzeżenia** jeśli dead-end > 20% lub isolated > 5%

**Dlaczego to ważne:** Wysoki SSA w dead-end zones to „fałszywa" korzyść.
W LC: dead-end = tailing i peak broadening.
W batch SPE: dead-end = wolna desorpcja i niski recovery.

#### 4. Accessible Surface Area (ASA) — probe-based
**Metoda:**
- Erozja void space o promień sondy (probe radius)
- Powierzchnia po erozji = ASA (dostępna dla cząsteczki o danym rozmiarze)
- Profile dla sond (probe radii): 0, 1, 2, 5, 10 nm

**Co pokazuje:**
- Wykres ASA vs promień sondy
- ASA/total_SA fraction dla każdego promienia
- Ile powierzchni jest realnie dostępne dla Twojego analitu

**Dlaczego to ważne:** Nie cała powierzchnia jest dostępna dla dużych cząsteczek.
Wąskie kanały i szyjki blokują dostęp → realna pojemność jest niższa.

#### 5. Throat / Constriction Analysis
**Metoda:**
- Local minima na medial axis kanałów (saddle points w EDT)
- Throat diameter = 2 × minimum radius w przewężeniu

**Co pokazuje:**
- Histogram throat sizes [µm]
- P10, P50, P90
- Throat-to-pore ratio (niższe = wyższy opór przepływu)
- **Ostrzeżenia** jeśli min throat < 30 µm lub ratio < 0.3

**Dlaczego to ważne:** Przewężenia dominują opór przepływu i ΔP.
Jeden wąski throat może zdominować cały element.

#### 6. Printability Score
**Metoda:** Sprawdzenie realizowalności druku DLP:
1. Min wall thickness vs pixel drukarki (22 µm)
2. Min channel width vs ryzyko bleed
3. Trapped resin fraction (zamknięty void)

**Scoring:**
- 85-100: EXCELLENT
- 70-84: GOOD
- 50-69: MARGINAL
- 30-49: POOR
- 0-29: FAIL

**Co sprawdza:**
- Wall < 1 pixel (22 µm) → FAIL (-40 pkt)
- Wall < 2 pixels (44 µm) → WARNING (-20 pkt)
- Channel < 2 pixels → FAIL (-30 pkt)
- Channel < 3 pixels → WARNING (-15 pkt)
- >10% trapped resin → FAIL (-20 pkt)
- >2% trapped resin → WARNING (-10 pkt)
- >10% thin walls → WARNING (-15 pkt)

---

### Zakładka: Cross-Section
Wizualizacja przekrojów 2D mesha.

**Jak używać:**
1. Wybierz płaszczyznę: XY, XZ, lub YZ
2. Slider: pozycja przekroju
3. Podgląd: solid=czarny, void=biały
4. Lokalna porowatość z przekroju

**Eksport:**
- PNG (obraz rastrowy)
- SVG (obraz wektorowy)

---

### Zakładka: Compare
Porównywanie do 6 projektów obok siebie.

**Workflow:**
1. Wygeneruj geometrię A → kliknij „Save to Compare"
2. Zmień parametry, wygeneruj B → „Save to Compare"
3. Przejdź do zakładki Compare → tabela porównawcza

**Co porównuje:**
- Wszystkie parametry (unit cell, wall, shape, quality)
- Statystyki podstawowe (vertices, faces, porosity, SSA)
- Eksport do CSV

---

## Architektura plików

```
gui_main.py              ← Główne okno, integracja zakładek (ROZSZERZONY z v2)
gyroid_math.py           ← Matematyka gyroidu (BEZ ZMIAN z v2)
mesh_generator.py        ← Generowanie meshów (BEZ ZMIAN z v2)
statistics_analyzer_v2.py← Statystyki instant + desorption (BEZ ZMIAN z v2)
statistics_tab_v2.py     ← GUI statystyk v2 (BEZ ZMIAN z v2)

viewer_3d.py             ← NOWY: wbudowany viewer 3D (pyqtgraph OpenGL)
predictions.py           ← NOWY: obliczenia ΔP, Van Deemter, t_desorption, capacity
predictions_tab.py       ← NOWY: GUI zakładki Predictions (wykresy matplotlib)
distributions.py         ← NOWY: voxel analysis (wall/channel/connectivity/ASA/throat/printability)
distributions_tab.py     ← NOWY: GUI zakładki Distributions (histogramy matplotlib)
cross_section.py         ← NOWY: generowanie przekrojów 2D
cross_section_tab.py     ← NOWY: GUI zakładki Cross-Section
compare_tab.py           ← NOWY: GUI porównywarka projektów

requirements.txt         ← Zaktualizowany
ROADMAP_PELNY.md         ← Plan rozwoju v4-v7
README.md                ← Ten plik
```

---

## Co zostaje bez zmian z v2

- ✅ Generowanie STL (wszystkie kształty)
- ✅ Export STL
- ✅ Instant statistics
- ✅ Desorption path (fast + sampling)
- ✅ gyroid_math.py, mesh_generator.py, statistics_analyzer_v2.py, statistics_tab_v2.py

---

## Znane ograniczenia v3

1. **Tortuosity = stała 1.41** — wartość literaturowa, nie liczona z mesha.
   Planowane: prawdziwa tortuosity geodezyjna w v4.
2. **Van Deemter: k = 2.0** — retention factor domyślny.
   Planowane: input w GUI w v4.
3. **Printability: pixel = 22 µm** — hardcoded dla drukarki Mariusza.
   Planowane: input w GUI.
4. **Connectivity: analiza Z-axis** — inlet=Z_min, outlet=Z_max.
   Dla sfer analiza we wszystkich kierunkach w v4.

---

## Kolejne wersje (ROADMAP)

Patrz: `ROADMAP_PELNY.md`

- **v4: Transport Intelligence** — prawdziwa tortuosity, RTD proxy, anizotropia, gradienty, batch sweep
- **v5: As-Printed Model** — korekcja geometrii po druku (bleed, shrinkage)
- **v6: Validation** — import danych eksperymentalnych, predicted vs measured
- **v7: High Fidelity** — solver przepływu LBM/Stokes, pore network model
