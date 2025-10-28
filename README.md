## Overview

This laboratory involves calibrating and characterizing space radiation detectors to analyze and compare their performance for different missions. You will develop scripts, analyze spectral data, and report your findings. The project is highly focused on collaborative data analysis pipelines and reproducible scientific reporting.

***

## Primary Goals

- **Write a min 5-page brief** comparing your results to others.
- **Discuss** which detectors are suitable for each type of mission.

***

## Experimental Tasks

### 1. Energy Calibration

- Write code to calibrate each detector, relating the digital channel number to absolute energy:  
  $$
  E = c_1 \times \text{channel} + c_0
  $$
- Test the linearity of this calibration curve and, if needed, fit higher-order terms.

### 2. Source Activity

- Determine the activities of laboratory sources based on their measured activity at a fixed date. Use half-lives for decay corrections.
- Note which sources you are using and their corresponding measurement dates.

### 3. Energy Resolution

- Calculate the energy resolution $$ R $$ at each energy:
  $$
  R = \frac{\delta E}{E}
  $$
  where $$ \delta E $$ is the FWHM of the peak.
- Fit resolution as a function of energy using:
  $$
  R^2 = \left(\frac{\delta E}{E}\right)^2 = aE^{-2} + bE^{-1} + c
  $$
  or express using logarithmically scaled axes for clarity.

### 4. Detector Efficiency

- Compute **absolute efficiency**:
  $$
  \epsilon = \frac{\text{count rate}}{\text{source activity}}
  $$
- Compute **intrinsic efficiency**:
  $$
  \epsilon = \frac{\text{count rate}}{\text{rate that photon passes through the detector}}
  $$
- Fit intrinsic efficiency (log-log axes recommended):
  $$
  \ln \epsilon_p = a + b \ln E + c (\ln E)^2
  $$

### 5. Off-Axis Response

- Measure the detector response as a function of the source angle.
- Collect and analyze spectra at multiple angles, plot (e.g., FWHM vs. angle).

***

## Data Analysis Pipeline

### General Guidelines

- Build modular Python scripts for each analysis step (e.g., calibration, resolution, efficiency).
- Use separate files/functions, e.g. `calibration.py`, `resolution.py`.
- Prefer argument parsing (e.g., Python’s `argparse`) for workflow flexibility:
  ```
  python pipeline.py --detector BGO --calibrate
  ```
- Avoid absolute paths in your code; use relative paths with Python modules like `os` or `pathlib`.

### Required Codes (For Each Detector)

- Calibration curve (channel ↔ energy/wavelength)
- Energy resolution vs. energy (plot FWHM of photo-peaks)
- Absolute and intrinsic efficiencies vs. energy (with plots)
- Off-axis response: change in FWHM/peak as a function of angle

Ensure all scripts work as a single, integrated pipeline.

***

## Data Collection

- Collect data from specified detectors in groups.
- Maintain clear records, consistent naming conventions, and document all procedures.
- Save spectra as `.mca` or `.Spe` files (plain text).

***

## Data Extraction and Processing

- Use Python scripts to import spectral data, fit peaks, perform background subtraction, and extract relevant parameters (energy, FWHM, count rate).
- Fit Gaussian or compound models to spectral peaks:
  $$
  f(x; \mu, \sigma^2, A) = \frac{A}{\sigma \sqrt{2\pi}} \, \exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)
  $$
- Report uncertainties and propagate errors in all results.

***

## Comparing Results

- After analysis, compare your group’s results to those of other groups.
- Produce summary plots (e.g., energy resolution vs. energy) combining all data.
- Comment on similarities, differences, strengths, weaknesses, and mission suitability for each detector.

***

## Reporting & Submission

- **Final report**: Concise, contains discussion, code, and comparison.
- Include code and analysis pipeline (ideally via GitHub).
- Discuss which detectors are best suited for specific mission types and why.

***

## Collaboration and Best Practices

- Use version control and shared repositories.
- Agreement on directory structure and file naming is essential.
- Communicate regularly to avoid missing or duplicated data.

***

