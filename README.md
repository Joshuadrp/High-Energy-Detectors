# 🔬 Gamma-Ray Detector Lab

> **A Python project for characterizing gamma-ray detectors in high-energy astrophysics**

This project analyzes the performance of three detector types used in space based gamma-ray astronomy missions:
- **NAITI** (Thallium-doped Sodium Iodide) - Scintillator detector
- **BGO** (Bismuth Germanate) - Scintillator detector  
- **CdTe** (Cadmium Telluride) - Solid-state detector

Built for the **PHYC40870 Space Detector Laboratory** course, it provides a complete analysis pipeline for energy calibration, resolution characterization, efficiency measurements, and off angle testing.

---

## Features

#### **Data Loading & Preprocessing**
- Load spectral data from `.Spe` (NAITI/BGO) and `.mca` (CDTE) files
- Extract metadata (live time, real time, detector info)
- Background subtraction

#### **Peak Fitting & Identification**
- Gaussian peak fitting with background subtraction
- Full uncertainty propagation

#### **Energy Calibration**
- Calibration (channel → energy conversion)
- Linear calibration: E = c1 × channel + c0
- Validation and sorting of calibration points
- YAML-based configuration for detector specific calibrations

#### **Energy Resolution**
- Re-fit peaks in calibrated energy space
- Energy resolution calculation: R = FWHM / E × 100%
- Resolution curve fitting: R² = aE⁻² + bE⁻¹ + c
- Both linear and log visualization

#### **Source Activity & Half-Life**
- Calculate current source activity from calibration date
- Account for radioactive decay over elapsed time
- Uncertainty propagation for activity measurements
- Support for Cs-137, Co-60, Am-241, Ba-133

#### **Efficiency Analysis**
- **Intrinsic efficiency**: Fraction of incident photons detected
- **Absolute efficiency**: Fraction of emitted photons detected
- Logarithmic polynomial fitting: ln ε = a + b ln E + c(ln E)²
- Full error propagation through all calculations

#### **Angular Response**
- Off-axis response characterization (peak amplitude vs angle)
- FWHM variation with detector angle
- Normalization to on-axis measurements

*Note: Co-60 excluded from CdTe analysis due to low detection efficiency at high energies*

### 📈 Outputs & Visualizations (in order)

1. **Gaussian Fit (Channel Space)** - Raw spectrum with fitted peaks
2. **Energy Calibration Curve** - Channel vs Energy
3. **Gaussian Fit (Energy)** - Calibrated spectrum with fitted peaks
4. **Energy Resolution vs Energy** - Linear scale
5. **Energy Resolution vs Energy** - Log scale with fitted curve
6. **Intrinsic Peak Efficiency vs Energy** - Log with polynomial fit
7. **Angular Response (Peak Amplitude)** - Count rate vs angle
8. **Angular Response (FWHM)** - Peak width vs angle

---