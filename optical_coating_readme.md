# Optical Coating Technologies Lab Analysis

## Overview

This project contains Python scripts used for analyzing experimental data from the "Optical Coating Technologies" practice lab. The objective is to gain insight into optical coating technologies by studying and applying fundamental characterization methods.

The scripts focus on determining the optical properties, specifically the **refractive index (n)** and **extinction coefficient (k)**, of common substrates (BK7, fused silica) and a single-layer Nb₂O₅ film.

The analysis is performed using data from two primary non-destructive methods:

1. **Spectrophotometry** – measuring Reflectance (R) and Transmittance (T)
2. **Ellipsometry** – measuring changes in light polarization (Ψ and Δ)

The results are then compared with established literature values.

### Course Information

- **Course:** Optical Coating Technologies
- **Institution:** ELI ALPS
- **Program:** LAScALA Erasmus Mundus Master
- **Semester:** 2025 Autumn Semester
- **Supervisor:** Veronika Hanyecz PhD

---

## 📁 Project Contents

- `Q3_coating.py` – Spectrophotometry analysis (Substrate n/k calculation and Swanepoel method for layer analysis)
- `Q4_coating.py` – Ellipsometry analysis (Bulk material n/k calculation)
- `syllabus_Opt_Coat_Tech_2025 (2).pdf` – Original lab manual and syllabus
- `*.csv` – Experimental data files (listed below)

---

## 🛠️ Dependencies

The scripts require the following Python libraries:

```bash
pip install pandas numpy matplotlib scipy
```

**Library purposes:**
- **pandas** – Loading and managing CSV data
- **numpy** – Numerical calculations and array manipulation
- **matplotlib** – Plotting results
- **scipy** – Peak finding functionality (`find_peaks`)

---

## 📂 Required Data Files

Ensure the following CSV files are in the same directory as the Python scripts:

**Substrate Data:**
- `1mm_T.Sample.Raw.csv` – Transmission data for 1mm Fused Silica
- `1mm_thick_R.Sample.Raw.csv` – Reflectance data for 1mm Fused Silica
- `6mm_T.Sample.Raw.csv` – Transmission data for 6mm BK7
- `6mm_thick_R.Sample.Raw.csv` – Reflectance data for 6mm BK7

**Layer Data:**
- `Nb2O4_T_TL.csv` – Transmission (T) and Lower Envelope (TL) data for Nb₂O₅ layer
- `Nb2O5_std.csv` – Literature (standard) refractive index values for Nb₂O₅

**Ellipsometry Data:**
- `50deg.csv` – Ellipsometry (Ψ, Δ) data at 50° incidence
- `55deg.csv` – Ellipsometry (Ψ, Δ) data at 55° incidence
- `60deg.csv` – Ellipsometry (Ψ, Δ) data at 60° incidence

---

## 🚀 Usage

Execute the scripts from your terminal. The scripts will automatically load the required data, perform calculations, and display plots.

**Run the Spectrophotometry (Swanepoel) analysis:**

```bash
python Q3_coating.py
```

**Run the Ellipsometry analysis:**

```bash
python Q4_coating.py
```

---

## 📊 Script Descriptions

### `Q3_coating.py` (Spectrophotometry Analysis)

This script performs analysis for Tasks 3, 4, 5, and 6 of the syllabus.

#### Part 1: Substrate Characterization (BK7 & Fused Silica)

- Loads transmission (T) and one-side reflectance (R) data for 1mm (Fused Silica) and 6mm (BK7) samples
- Converts T and R from percentage (0-100) to fractional (0-1.0) values
- Calculates extinction coefficient (k) using:
  - T = (1-R)² · e^(-αl) where α = 4πk/λ
- Calculates refractive index (n) by solving the quadratic equation derived from:
  - R = [(n-1)² + k²] / [(n+1)² + k²]
- Generates plots for T/R vs. wavelength and calculated n and k vs. wavelength for both substrates

#### Part 2: Single Layer Characterization (Swanepoel Method)

- Characterizes the Nb₂O₅ layer using the **Swanepoel method**
- Loads layer transmission (`T_swan`) and lower envelope (`TL_swan`) from CSV
- Uses substrate refractive index (`n_sub`) from Part 1 (6mm BK7 sample)
- Implements calculation functions for the transparent region based on syllabus formulas
- Identifies transmission peaks (maxima) using `scipy.signal.find_peaks`
- Calculates layer thickness (d) using wavelengths (λ₁, λ₂) and refractive indices (n₁, n₂) at two adjacent maxima
- Plots calculated layer refractive index and compares to literature data

### `Q4_coating.py` (Ellipsometry Analysis)

This script performs analysis for Tasks 8 and 9 of the syllabus.

- Loads ellipsometric angles Ψ (psi) and Δ (delta) for BK7 sample at 50°, 55°, and 60° incidence angles
- Converts wavelength data from eV to nm
- Calculates complex reflection ratio: ρ = tan(Ψ)e^(iΔ)
- Determines optical constants (n and k) by:
  - Calculating complex relative dielectric constant ⟨ε̃⟩
  - Solving relations: ε₁ = n² - k² and ε₂ = 2nk
- Generates plots for:
  - Ψ vs. wavelength
  - Δ vs. wavelength
  - Calculated n and k vs. wavelength for all three incidence angles

---

## 📈 Output

Both scripts generate matplotlib figures showing:
- Raw experimental data (T, R, Ψ, Δ)
- Calculated optical constants (n, k) vs. wavelength
- Comparison with literature values (where applicable)

---

## 📝 Notes

- Ensure all data files are in the correct format (CSV with appropriate headers)
- Wavelength units should be consistent across all input files
- The scripts assume standard CSV formatting with comma separators

---

## 👨‍💻 Author

**Code written by:** Isuru Withanawasam