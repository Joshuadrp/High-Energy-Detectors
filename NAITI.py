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
print(cs137['channels'])

# Subtract background
print("Subtracting background...")
cs137 = subtract_background(cs137, background)
co60 = subtract_background(co60, background)
am241 = subtract_background(am241, background)
ba133 = subtract_background(ba133, background)

# Fit and plot peaks
cs137_fit = plot_spectrum_with_fit(cs137, title="CS-137", peak_channel=300)
co60_fit = plot_spectrum_with_fit(co60, "Co-60", peak_channel=490)
am241_fit = plot_spectrum_with_fit(am241, "Am-241", peak_channel=50)
ba133_fit = plot_spectrum_with_fit(ba133, "Ba-133", peak_channel=150)

# UNCERTAINTIES USING CHANNLES ----------------------------------------
pcov_cs137 = cs137_fit['pcov']
pcov_co60 = co60_fit['pcov']
pcov_am241 = am241_fit['pcov']
pcov_ba133 = ba133_fit['pcov']
errors = []

pcovs_isotopes = [pcov_co60, pcov_am241, pcov_ba133,pcov_cs137]
for i in pcovs_isotopes:
    errors.append(np.sqrt(np.diag(i)))

# Print uncertainties (CHANNEL UNCERTAINTIES)
print(f"\nCO60:\nMean (mu): {co60_fit['mu']:.2f} ± {errors[0][1]:.2f}")
print(f"\nAM247:\nMean (mu): {am241_fit['mu']:.2f} ± {errors[1][1]:.2f}")
print(f"\nBA133:\nMean (mu): {ba133_fit['mu']:.2f} ± {errors[2][1]:.2f}")
print(f"\nCS137:\nMean (mu): {cs137_fit['mu']:.2f} ± {errors[3][1]:.2f}")

"""
END OF CHANNEL UNCERTAINTIES
"""

# CALIBRATION EQUATION RESULTS FROM CALIBRATION FUNCTION FILE
naiti_file = "Foreigners/NaITI/NaITI.yaml"
x,y,coeffs, m,b, m_uncert,b_uncert = calibrate(naiti_file)
print(m_uncert,b_uncert)
plot_calibration(x, y, coeffs)

cs137_energy = energy_calibration_equation(m,cs137['channels'],b)
am241_energy = energy_calibration_equation(m,am241['channels'],b)
ba133_energy = energy_calibration_equation(m,ba133['channels'],b)
co60_energy = energy_calibration_equation(m,co60['channels'],b)

cs137_fit_energy = plot_spectrum_with_fit_energy(energy=cs137_energy, y=cs137['counts'], title="CS-137", peak_channel=700)
cs137_prop_error=propagate_energy_uncertainty(cs137_fit['mu'], errors[3][1], m, b, m_uncert, b_uncert)
print(cs137_prop_error)

co60_fit_energy = plot_spectrum_with_fit_energy(energy=co60_energy,y=co60['counts'], title="Co-60", peak_channel=1332 )
co60_prop_error = propagate_energy_uncertainty(co60_fit['mu'], errors[0][1], m, b, m_uncert, b_uncert)
print(co60_prop_error)

am241_fit_energy = plot_spectrum_with_fit_energy(energy=am241_energy,y=am241['counts'], title="AM-241", peak_channel=50 )
am241_prop_error = propagate_energy_uncertainty(am241_fit['mu'], errors[1][1], m, b, m_uncert, b_uncert)
print(am241_prop_error)

ba133_fit_energy = plot_spectrum_with_fit_energy(energy=ba133_energy,y=ba133['counts'], title="BA-133", peak_channel=356 )
ba133_prop_error=propagate_energy_uncertainty(ba133_fit['mu'], errors[2][1], m, b, m_uncert, b_uncert)
print(ba133_prop_error)

"""
Efficiency functions
"""
photon_amt_Cs, cesium = calc_half_life('137-Cs')
photon_amt_Am, amer = calc_half_life('241-Am')
photon_amt_Ba, bar = calc_half_life('133-Ba')
photon_amt_Co, cob = calc_half_life('60-Co')

cs_intrinsic = intrinsic()
am_instrinsic = intrinsic()
ba_intrinsic = intrinsic()
co_intrinsic = intrinsic()

Cs_eff = efficiency(cesium, cs137['counts'],photon_amt_Cs,cs137_energy,288.44,cs_intrinsic)
Am_eff = efficiency(amer, am241['counts'],photon_amt_Am,am241_energy,28.84,am_instrinsic)
Ba_eff = efficiency(bar, ba133['counts'],photon_amt_Ba,ba133_energy,157.25,ba_intrinsic)
Co_eff = efficiency(cob, co60['counts'],photon_amt_Co,co60_energy,493,co_intrinsic)

