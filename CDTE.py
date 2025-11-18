from high_detectors_functions import *

base_path = "Foreigners/CDTE"
isotope_names = ['CS137', 'BA', 'AM']
peak_channels = [639, 213, 420]
known_energies = [662, 276, 59.54]

# Load data
data = [load_spe(f"{base_path}/{name}_aligned.mca") for name in isotope_names]
background = load_spe(f"{base_path}/BACKGROUND.mca")

# Subtract background
data_no_bg = [subtract_background(d, background) for d in data]

# Fit peaks
fits = []
print("Fitting...")
for i, name in enumerate(isotope_names):
    fit = plot_spectrum_with_fit(data_no_bg[i], name, peak_channels[i])
    fits.append(fit)
    print(f"{name}: {fit['mu']:.2f} ± {fit['mu_err']:.2f}")

# Calibration
x, y, coeffs, m, b, m_uncert, b_uncert = calibrate("Foreigners/CDTE/CDTE.yaml")
plot_calibration(x, y, coeffs)

print("\nEnergy fit...")
energy_fits = []
total_energy = []
for i, name in enumerate(isotope_names):
    energy = energy_calibration_equation(m, data_no_bg[i]['channels'], b)
    total_energy.append(energy)

    peak_energy_estimate = m * fits[i]['mu'] + b

    energy_fit = plot_spectrum_with_fit_energy(energy, data_no_bg[i]['counts'],
                                               name, peak_energy_estimate,
                                               window=100)
    energy_fits.append(energy_fit)
    error = propagate_energy_uncertainty(fits[i]['mu'], fits[i]['mu_err'],
                                         m, b, m_uncert, b_uncert)

    print(f"{name}: {energy_fit['mu']:.2f} ± {error:.2f} keV")

# ENERGY RESOLUTION
print("\nENERGY RESOLUTION")
resolutions = calculate_resolution(
    isotope_names,
    fits,
    energy_fits,
    m, b,
    known_energies
)
# Plot Resolution vs Energy
plot_resolution_vs_energy(resolutions, detector_name='CDTE')
a, b_coeff, c, fit_errors = fit_resolution_curve(resolutions)

# """
# Efficiency functions
# """
photon_amt_Cs, cesium, cs_act_err = calc_half_life('137-Cs')
photon_amt_Am, amer, am_act_err = calc_half_life('241-Am')
photon_amt_Ba, bar, bar_act_err = calc_half_life('133-Ba')
# photon_amt_Co, cob, cob_act_err = calc_half_life('60-Co')

#intrinsic efficiency rates and errors
CdTE_detector_parameters = (5.08, 16)
cs_intrinsic, cs_intrinsic_err = intrinsic(photon_amt_Cs,cs_act_err,*CdTE_detector_parameters)
am_instrinsic, am_intrinsic_err = intrinsic(photon_amt_Am,am_act_err,*CdTE_detector_parameters)
ba_intrinsic, ba_intrinsic_err= intrinsic(photon_amt_Ba,bar_act_err,*CdTE_detector_parameters)
# co_intrinsic, co_intrinsic_err = intrinsic(photon_amt_Co,cob_act_err,*CdTE_detector_parameters)

#absolute and intrinisc efficiencies and errors
cs_abs_eff, cs_int_eff, cs_abs_eff_err, cs_int_eff_err = efficiency_uncertainty(
    nuclide='137-Cs',
    peak_counts = fits[0]['A'],                   
    peak_counts_err=fits[0]['amp_err'],           
    time=data_no_bg[0]['real_time'],                   
    emitted_counts=photon_amt_Cs,               
    emitted_counts_err=cs_act_err,       
    incident_counts=cs_intrinsic,               
    incident_counts_err=cs_intrinsic_err        
)

# co_abs_eff, co_int_eff, co_abs_eff_err, co_int_eff_err = efficiency_uncertainty(
#     nuclide='60-Co',
#     peak_counts = fits[1]['A'],                    
#     peak_counts_err=fits[1]['amp_err'],            
#     time=data_no_bg[1]['real_time'],                   
#     emitted_counts=photon_amt_Co,               
#     emitted_counts_err=cob_act_err,      
#     incident_counts=co_intrinsic,               
#     incident_counts_err=co_intrinsic_err        
# )

am_abs_eff, am_int_eff, am_abs_eff_err, am_int_eff_err = efficiency_uncertainty(
    nuclide='241-Am',
    peak_counts = fits[2]['A'],                   
    peak_counts_err=fits[2]['amp_err'],          
    time=data_no_bg[2]['real_time'],                  
    emitted_counts=photon_amt_Am,               
    emitted_counts_err=am_act_err,      
    incident_counts=am_instrinsic,               
    incident_counts_err=am_intrinsic_err 
)

ba_abs_eff, ba_int_eff, ba_abs_eff_err, ba_int_eff_err = efficiency_uncertainty(
    nuclide='133-Ba',
    peak_counts = fits[1]['A'],                    
    peak_counts_err=fits[1]['amp_err'],            
    time=data_no_bg[1]['real_time'],                    
    emitted_counts=photon_amt_Ba,               
    emitted_counts_err=bar_act_err,       
    incident_counts=ba_intrinsic,               
    incident_counts_err=ba_intrinsic_err        
)
energies = [356, 662, 59.54]
intrinsic_efficiencies = [ba_int_eff, cs_int_eff, am_int_eff]
intrinsic_eff_errors = [ba_int_eff_err, cs_int_eff_err,am_int_eff_err]

# Fit the logarithmic polynomial (Equation 20)
intrinsic_result = fit_intrinsic_efficiency(energies, intrinsic_efficiencies, degree=1)

# Plot the data with the fit on log-log axes
plot_intrinsic_efficiency(energies, intrinsic_efficiencies, intrinsic_eff_errors, intrinsic_result, title='CDTE Intrinsic Peak Efficiency vs Energy')
#off axis response
ba_angular = fit_off_axis_response("Foreigners/CDTE/Ba_offaxis", "Ba", background, 213,on_axis_file='Foreigners/CDTE/BA_aligned.mca')
plot_off_axis_response((ba_angular), 'Ba133',276)

ba_angular_FWHM = fit_off_axis_response_FWHM("Foreigners/CDTE/Ba_offaxis", "Ba", background, 213,on_axis_file='Foreigners/CDTE/BA_aligned.mca')
plot_off_axis_response_FWHM((ba_angular_FWHM), 'Ba133',276)
