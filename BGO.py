from high_detectors_functions import *

base_path = "Foreigners/BGO"
isotope_names = ['CS137', 'CO60', 'AM', 'BA']
peak_channels = [290, 520, 50, 180]
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