from high_detectors_functions import *
# ============================================================================
# MAIN
# ============================================================================
base_path = "Foreigners/NaITI"

# Load data
cs137 = load_spe(f"{base_path}/CS137_aligned.Spe")
co60 = load_spe(f"{base_path}/CO60_aligned.Spe")
am241 = load_spe(f"{base_path}/AM_aligned.Spe")
ba133 = load_spe(f"{base_path}/BA_aligned.Spe")
background = load_spe(f"{base_path}/BACKGROUND.Spe")


# Subtract background
print("Subtracting background...")
cs137_no_bg = subtract_background(cs137, background)
co60_no_bg = subtract_background(co60, background)
am241_no_bg = subtract_background(am241, background)
ba133_no_bg = subtract_background(ba133, background)

# Fit and plot peaks
cs137_fit = plot_spectrum_with_fit(cs137_no_bg, "CS-137", peak_channel=300)
co60_fit = plot_spectrum_with_fit(co60_no_bg, "Co-60", peak_channel=490)
am241_fit = plot_spectrum_with_fit(am241_no_bg, "Am-241", peak_channel=50)
ba133_fit = plot_spectrum_with_fit(ba133_no_bg, "Ba-133", peak_channel=150)

# UNCERTAINTIES USING CHANNELS ----------------------------------------
pcov_cs137 = cs137_fit['pcov']
pcov_co60 = co60_fit['pcov']
pcov_am241 = am241_fit['pcov']
pcov_ba133 = ba133_fit['pcov']
errors = []

pcovs_isotopes = [pcov_co60, pcov_am241, pcov_ba133,pcov_cs137]
for i in pcovs_isotopes:
    errors.append(np.sqrt(np.diag(i)))

# Print uncertainties (CHANNEL UNCERTAINTIES)
print(f"CO60:\nMean (mu): {co60_fit['mu']:.2f} ± {co60_fit['mu_err']:.2f}")
print(f"AM247:\nMean (mu): {am241_fit['mu']:.2f} ± {am241_fit['mu_err']:.2f}")
print(f"BA133:\nMean (mu): {ba133_fit['mu']:.2f} ± {ba133_fit['mu_err']:.2f}")
print(f"CS137:\nMean (mu): {cs137_fit['mu']:.2f} ± {cs137_fit['mu_err']:.2f}")

"""
END OF CHANNEL UNCERTAINTIES
"""

# CALIBRATION EQUATION RESULTS FROM CALIBRATION FUNCTION FILE
naiti_file = "Foreigners/NaITI/NaITI.yaml"
x,y,coeffs, m,b, m_uncert,b_uncert = calibrate(naiti_file)
# print(m_uncert,b_uncert)
plot_calibration(x, y, coeffs)

cs137_energy = energy_calibration_equation(m,cs137_no_bg['channels'],b)
am241_energy = energy_calibration_equation(m,am241_no_bg['channels'],b)
ba133_energy = energy_calibration_equation(m,ba133_no_bg['channels'],b)
co60_energy = energy_calibration_equation(m,co60_no_bg['channels'],b)

cs137_fit_energy = plot_spectrum_with_fit_energy(energy=cs137_energy, y=cs137_no_bg['counts'], title="CS-137", peak_channel=700)
cs137_prop_error=propagate_energy_uncertainty(cs137_fit['mu'], errors[3][1], m, b, m_uncert, b_uncert)

co60_fit_energy = plot_spectrum_with_fit_energy(energy=co60_energy,y=co60_no_bg['counts'], title="Co-60", peak_channel=1332 )
co60_prop_error = propagate_energy_uncertainty(co60_fit['mu'], errors[0][1], m, b, m_uncert, b_uncert)

am241_fit_energy = plot_spectrum_with_fit_energy(energy=am241_energy,y=am241_no_bg['counts'], title="AM-241", peak_channel=59 )
am241_prop_error = propagate_energy_uncertainty(am241_fit['mu'], errors[1][1], m, b, m_uncert, b_uncert)

ba133_fit_energy = plot_spectrum_with_fit_energy(energy=ba133_energy,y=ba133_no_bg['counts'], title="BA-133", peak_channel=356 )
ba133_prop_error=propagate_energy_uncertainty(ba133_fit['mu'], errors[2][1], m, b, m_uncert, b_uncert)

# Print propagated uncertainties
print(f"\nCO60:\nMean (mu): {co60_fit_energy['mu']:.2f} ± {co60_prop_error:.2f}")
print(f"CS137:\nMean (mu): {cs137_fit_energy['mu']:.2f} ± {cs137_prop_error:.2f}")
print(f"AM241:\nMean (mu): {am241_fit_energy['mu']:.2f} ± {am241_prop_error:.2f}")
print(f"BA133:\nMean (mu): {ba133_fit_energy['mu']:.2f} ± {ba133_prop_error:.2f}")






