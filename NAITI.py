from high_detectors_functions import *

base_path = "Foreigners/NaITI"
isotope_names = ['CS137', 'CO60', 'AM', 'BA']
peak_channels = [300, 490, 50, 150]
known_energies = [662, 1332, 59.5, 356]  #KeV

# Load data
data = [load_spe(f"{base_path}/{name}_aligned.Spe") for name in isotope_names]
background = load_spe(f"{base_path}/BACKGROUND.Spe")

# Subtract background
data_no_bg = [subtract_background(d, background) for d in data]

# Fit peaks
fits = []
for i, name in enumerate(isotope_names):
    fit = plot_spectrum_with_fit(data_no_bg[i], name, peak_channels[i])
    fits.append(fit)
    print(f"{name}: {fit['mu']:.2f} ± {fit['mu_err']:.2f}")

# Calibration
x, y, coeffs, m, b, m_uncert, b_uncert = calibrate("Foreigners/NaITI/NaITI.yaml")
plot_calibration(x, y, coeffs)

# Convert to energy and plot
for i, name in enumerate(isotope_names):
    energy = energy_calibration_equation(m, data_no_bg[i]['channels'], b)
    energy_fit = plot_spectrum_with_fit_energy(energy, data_no_bg[i]['counts'],
                                                name, known_energies[i])
    error = propagate_energy_uncertainty(fits[i]['mu'], fits[i]['mu_err'],
                                         m, b, m_uncert, b_uncert)
    print(f"{name}: {energy_fit['mu']:.2f} ± {error:.2f} keV")

# """
# Efficiency functions
# """
photon_amt_Cs, cesium, cs_act_err = calc_half_life('137-Cs')
# photon_amt_Am, amer, am_act_err = calc_half_life('241-Am')
# photon_amt_Ba, bar, bar_act_err = calc_half_life('133-Ba')
# photon_amt_Co, cob, cob_act_err = calc_half_life('60-Co')
#
cs_intrinsic, cs_intrinsic_err = intrinsic(photon_amt_Cs,cs_act_err,5.08,16)
print(cs_intrinsic)
#am_instrinsic = intrinsic()
#ba_intrinsic = intrinsic()
#co_intrinsic = intrinsic()

#Cs_eff = efficiency_uncertainty(cesium, cs137['counts'],photon_amt_Cs,cs137_energy,288.44,cs_intrinsic,cs137['real_time'])
#Am_eff = efficiency(amer, am241['counts'],photon_amt_Am,am241_energy,28.84,am_instrinsic)
#Ba_eff = efficiency(bar, ba133['counts'],photon_amt_Ba,ba133_energy,157.25,ba_intrinsic)
#Co_eff = efficiency(cob, co60['counts'],photon_amt_Co,co60_energy,493,co_intrinsic)

abs_eff, int_eff, abs_eff_err, int_eff_err = efficiency_uncertainty(
    nuclide='137-Cs',
    counts=CS137['counts'],
    energy=cs137_energy,
    peak_energy=288.44,                    # e.g., 15000 counts
    peak_counts_err=errors[3][0],            # e.g., 150 counts
    time=cs137['real_time'],                    # e.g., 3600 seconds
    emitted_counts=photon_amt_Cs,               # e.g., 50000 photons/s (total emission)
    emitted_counts_err=cs_act_err,       # e.g., 500 photons/s
    incident_counts=cs_intrinsic,               # e.g., 800 photons/s (hitting detector)
    incident_counts_err=cs_intrinsic_err        # e.g., 8 photons/s
)
