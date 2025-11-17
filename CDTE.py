from high_detectors_functions import *

base_path = "Foreigners/CDTE"
isotope_names = ['CS137', 'AM', 'BA'] #CdTe can't reliably detect Co-60
peak_channels = [222, 32, 83]
known_energies = [662,59.5,81]

# Load data
data = [load_spe(f"{base_path}/{name}_aligned.mca") for name in isotope_names]
# background = load_spe(f"{base_path}/BACKGROUND.mca")

# Subtract background
# data_no_bg = [subtract_background(d, background) for d in data]

# Fit peaks
fits = []
print("Fitting...")
for i, name in enumerate(isotope_names):
    fit = plot_spectrum_with_fit(data[i], name, peak_channels[i])
    fits.append(fit)
    print(f"{name}: {fit['mu']:.2f} ± {fit['mu_err']:.2f}")

# Calibration
x, y, coeffs, m, b, m_uncert, b_uncert = calibrate("Foreigners/CDTE/CDTE.yaml")
plot_calibration(x, y, coeffs)

# print("\nEnergy fit...")
# energy_fits = []
# total_energy = []
# for i, name in enumerate(isotope_names):
#     energy = energy_calibration_equation(m, data[i]['channels'], b)
#     total_energy.append(energy)
#
#     peak_energy_estimate = m * fits[i]['mu'] + b
#
#     energy_fit = plot_spectrum_with_fit_energy(energy, data[i]['counts'],
#                                                name, peak_energy_estimate,
#                                                window=50)
#     energy_fits.append(energy_fit)
#     error = propagate_energy_uncertainty(fits[i]['mu'], fits[i]['mu_err'],
#                                          m, b, m_uncert, b_uncert)
#
#     print(f"{name}: {energy_fit['mu']:.2f} ± {error:.2f} keV")
#
# # ENERGY RESOLUTION
# print("\nENERGY RESOLUTION")
# resolutions = calculate_resolution(
#     isotope_names,
#     fits,
#     energy_fits,
#     m, b,
#     known_energies
# )
# # Plot Resolution vs Energy
# plot_resolution_vs_energy(resolutions, detector_name='CDTE')
#
# # Fit the resolution curve: R² = aE⁻² + bE⁻¹ + c
# a, b_coeff, c, fit_errors = fit_resolution_curve(resolutions)



