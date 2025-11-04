from high_detectors_functions import *
from calibration_func import *

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

# CALIBRATION EQUATION RESULTS FROM CALIBRATION FUNCTION FILE
naiti_file = "Foreigners/NaITI/NaITI.yaml"
x,y,coeffs, m,b = calibrate(naiti_file)
plot_calibration(x, y, coeffs)

test_channels = cs137['channels']
energy = energy_calibration_equation(m,cs137['channels'],b)
print(energy)

# Fit and plot peaks
print("\nCS-137")
cs137_fit = plot_spectrum_with_fit(cs137, energy=energy, title="CS-137", peak_channel=300)
print("\nCo-60")
co60_fit = plot_spectrum_with_fit(co60, "Co-60", peak_channel=490)
print("\nAm-241")
am241_fit = plot_spectrum_with_fit(am241, "Am-241", peak_channel=50)
print("\nBa-133")
ba133_fit = plot_spectrum_with_fit(ba133, "Ba-133", peak_channel=150)

# UNCERTAINTIES
pcov_cs137 = cs137_fit['pcov']
pcov_co60 = co60_fit['pcov']
pcov_am241 = am241_fit['pcov']
pcov_ba133 = ba133_fit['pcov']
errors = []

pcovs_isotopes = [pcov_co60, pcov_am241, pcov_ba133,pcov_cs137]
for i in pcovs_isotopes:
    errors.append(np.sqrt(np.diag(i)))

# Print uncertainties
print(f"\nCS137:\nMean (mu): {cs137_fit['mu']:.2f} ± {errors[3][1]:.2f}")
print(f"\nCO60:\nMean (mu): {co60_fit['mu']:.2f} ± {errors[0][1]:.2f}")
print(f"\nAM247:\nMean (mu): {am241_fit['mu']:.2f} ± {errors[1][1]:.2f}")
print(f"\nBA133:\nMean (mu): {ba133_fit['mu']:.2f} ± {errors[2][1]:.2f}\n")


