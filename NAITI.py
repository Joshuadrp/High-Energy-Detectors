from high_detectors_functions import *

base_path = "Foreigners/NaITI"
isotope_names = ['CS137', 'CO60', 'AM', 'BA']
peak_channels = [300, 490, 50, 150]
known_energies = [662, 1173, 59.5, 356]  #KeV

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
    print(f"{name}: {fit['A']}")
# Calibration
x, y, coeffs, m, b, m_uncert, b_uncert = calibrate("Foreigners/NaITI/NaITI.yaml")
plot_calibration(x, y, coeffs)

# Convert to energy and plot
total_energy = []
energy_fits = []
for i, name in enumerate(isotope_names):
    energy = energy_calibration_equation(m, data_no_bg[i]['channels'], b)
    total_energy.append(energy)
    
    peak_energy_estimate = m * fits[i]['mu'] + b
    
    energy_fit = plot_spectrum_with_fit_energy(energy, data_no_bg[i]['counts'],
                                                name, peak_energy_estimate)
    energy_fits.append(energy_fit)
    error = propagate_energy_uncertainty(fits[i]['mu'], fits[i]['mu_err'],
                                         m, b, m_uncert, b_uncert)
    print(f"{name}: {energy_fit['mu']:.2f} ± {error:.2f} keV")

# """
# Efficiency functions
# """
photon_amt_Cs, cesium, cs_act_err = calc_half_life('137-Cs')
photon_amt_Am, amer, am_act_err = calc_half_life('241-Am')
photon_amt_Ba, bar, bar_act_err = calc_half_life('133-Ba')
photon_amt_Co, cob, cob_act_err = calc_half_life('60-Co')
#intrinsic efficiency rates and errors
NAITI_detector_parameters = (5.08, 16)
cs_intrinsic, cs_intrinsic_err = intrinsic(photon_amt_Cs,cs_act_err,*NAITI_detector_parameters)
am_instrinsic, am_intrinsic_err = intrinsic(photon_amt_Am,am_act_err,*NAITI_detector_parameters)
ba_intrinsic, ba_intrinsic_err= intrinsic(photon_amt_Ba,bar_act_err,*NAITI_detector_parameters)
co_intrinsic, co_intrinsic_err = intrinsic(photon_amt_Co,cob_act_err,*NAITI_detector_parameters)

#absolute and intrinisc efficiencies and errors
cs_abs_eff, cs_int_eff, cs_abs_eff_err, cs_int_eff_err = efficiency_uncertainty(
    nuclide='137-Cs',
    peak_counts = fits[0]['A'],                    # e.g., 15000 counts
    peak_counts_err=fits[0]['amp_err'],            # e.g., 150 counts
    time=data_no_bg[0]['real_time'],                    # e.g., 3600 seconds
    emitted_counts=photon_amt_Cs,               # e.g., 50000 photons/s (total emission)
    emitted_counts_err=cs_act_err,       # e.g., 500 photons/s
    incident_counts=cs_intrinsic,               # e.g., 800 photons/s (hitting detector)
    incident_counts_err=cs_intrinsic_err        # e.g., 8 photons/s
)

co_abs_eff, co_int_eff, co_abs_eff_err, co_int_eff_err = efficiency_uncertainty(
    nuclide='60-Co',
    peak_counts = fits[1]['A'],                    # e.g., 15000 counts
    peak_counts_err=fits[1]['amp_err'],            # e.g., 150 counts
    time=data_no_bg[1]['real_time'],                    # e.g., 3600 seconds
    emitted_counts=photon_amt_Co,               # e.g., 50000 photons/s (total emission)
    emitted_counts_err=cob_act_err,       # e.g., 500 photons/s
    incident_counts=co_intrinsic,               # e.g., 800 photons/s (hitting detector)
    incident_counts_err=co_intrinsic_err        # e.g., 8 photons/s
)

am_abs_eff, am_int_eff, am_abs_eff_err, am_int_eff_err = efficiency_uncertainty(
    nuclide='241-Am',
    peak_counts = fits[2]['A'],                    # e.g., 15000 counts
    peak_counts_err=fits[2]['amp_err'],            # e.g., 150 counts
    time=data_no_bg[2]['real_time'],                    # e.g., 3600 seconds
    emitted_counts=photon_amt_Am,               # e.g., 50000 photons/s (total emission)
    emitted_counts_err=am_act_err,       # e.g., 500 photons/s
    incident_counts=am_instrinsic,               # e.g., 800 photons/s (hitting detector)
    incident_counts_err=am_intrinsic_err 
)

ba_abs_eff, ba_int_eff, ba_abs_eff_err, ba_int_eff_err = efficiency_uncertainty(
    nuclide='133-Ba',
    peak_counts = fits[3]['A'],                    # e.g., 15000 counts
    peak_counts_err=fits[3]['amp_err'],            # e.g., 150 counts
    time=data_no_bg[3]['real_time'],                    # e.g., 3600 seconds
    emitted_counts=photon_amt_Ba,               # e.g., 50000 photons/s (total emission)
    emitted_counts_err=bar_act_err,       # e.g., 500 photons/s
    incident_counts=ba_intrinsic,               # e.g., 800 photons/s (hitting detector)
    incident_counts_err=ba_intrinsic_err        # e.g., 8 photons/s
)

#off axis response
cs_angular = fit_off_axis_response("Foreigners/NaITI/CS_offaxis", "CS", background, 300,on_axis_file='Foreigners/NaITI/CS137_aligned.Spe')
plot_off_axis_response((cs_angular), 'cs137',661.57)

Am_angular = fit_off_axis_response("Foreigners/NaITI/AM_offaxis", "AM", background, 50,on_axis_file='Foreigners/NaITI/AM_aligned.Spe')
plot_off_axis_response((Am_angular), 'Am241',59.54)

cs_angular_FWHM = fit_off_axis_response_FWHM("Foreigners/NaITI/CS_offaxis", "CS", background, 300,on_axis_file='Foreigners/NaITI/CS137_aligned.Spe')
plot_off_axis_response_FWHM((cs_angular_FWHM), 'cs137',661.57)

AM_angular_FWHM = fit_off_axis_response_FWHM("Foreigners/NaITI/AM_offaxis", "AM", background, 50,on_axis_file='Foreigners/NaITI/AM_aligned.Spe')
plot_off_axis_response_FWHM((AM_angular_FWHM), 'AM241',59.54)