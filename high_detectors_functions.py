import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import yaml

"""
DATA LOADING FUNCTION
"""
def load_spe(filepath):
    """Load .spe (NaI/BGO) or .mca (CdTe) file and return spectrum with metadata"""

    # Check file extension
    if filepath.endswith('.mca'):
        # Binary format - CdTe detector
        try:
            counts = np.fromfile(filepath, dtype=np.uint32)
            counts = counts[128:]  # Skip header
        except:
            counts = np.fromfile(filepath, dtype=np.uint16)
            counts = counts[256:]

        channels = np.arange(len(counts))

        return {
            'channels': channels,
            'counts': counts,
            'filename': filepath.split('/')[-1],
            'detector': 'CdTe',
            'live_time': None,
            'real_time': None,
        }

    # ASCII format - .spe files (NaI/BGO)
    with open(filepath, 'r') as f:
        lines = f.readlines()

    metadata = {
        'filename': filepath.split('/')[-1],
        'detector': 'NaI',
        'live_time': None,
        'real_time': None,
    }

    for i, line in enumerate(lines):
        if '$MEAS_TIM:' in line and i + 1 < len(lines):
            times = lines[i + 1].strip().split()
            if len(times) >= 2:
                metadata['live_time'] = float(times[0])
                metadata['real_time'] = float(times[1])

    data_start = False
    skip_next = False
    counts = []

    for line in lines:
        line = line.strip()
        if '$DATA:' in line:
            data_start = True
            skip_next = True
            continue
        if skip_next:
            skip_next = False
            continue
        if '$' in line and data_start:
            break
        if data_start and line and line[0].isdigit():
            counts.append(float(line))

    channels = np.arange(len(counts))
    counts = np.array(counts)

    return {
        'channels': np.array(channels),
        'counts': counts,
        'filename': metadata['filename'],
        'detector': metadata['detector'],
        'live_time': metadata['live_time'],
        'real_time': metadata['real_time'],
    }

"""
PEAK FITTING
"""


def subtract_background(signal_data, background_data):
    """Subtract background spectrum from signal spectrum"""
    signal_counts = signal_data['counts'].copy()
    background_counts = background_data['counts']

    # Handle different lengths - min length
    min_length = min(len(signal_counts), len(background_counts))

    signal_counts = signal_counts[:min_length]
    background_counts = background_counts[:min_length]

    # Subtract and ensure no negative counts
    corrected_counts = signal_counts - background_counts
    corrected_counts = np.maximum(corrected_counts, 0)

    return {
        'channels': signal_data['channels'][:min_length],
        'counts': corrected_counts,
        'filename': signal_data['filename'],
        'detector': signal_data['detector'],
        'live_time': signal_data['live_time'],
        'real_time': signal_data['real_time'],
    }


def gaussian(x, A, mu, sigma):
    """Gaussian function"""
    return A * np.exp(-(x - mu) ** 2 / (2 * sigma ** 2))

def gaussian_with_background(x, A, mu, sigma, bg):
    """Gaussian with background"""
    return gaussian(x, A, mu, sigma) + bg

def fit_peak(channels, counts, peak_channel, window=50):
    """Fit Gaussian to a peak"""
    center_idx = np.argmin(np.abs(channels - peak_channel))
    start = max(0, center_idx - window)
    end = min(len(channels), center_idx + window)

    region_ch = channels[start:end]
    region_counts = counts[start:end]

    A_guess = np.sum(region_counts)
    mu_guess = region_ch[np.argmax(region_counts)]

    # ADAPTIVE SIGMA GUESS - scales with peak position
    # For channels: ~5, for energy (keV): ~peak_energy/100
    sigma_guess = max(5.0, mu_guess / 100)  # ← FIX HERE

    bg_guess = np.min(region_counts)

    popt, pcov = curve_fit(gaussian_with_background, region_ch, region_counts,
                           p0=[A_guess, mu_guess, sigma_guess, bg_guess],
                           maxfev=10000)

    A, mu, sigma, bg = popt
    errors = np.sqrt(np.diag(pcov))
    FWHM = 2.355 * sigma

    return {
        'mu': mu,
        'sigma': sigma,
        'FWHM': FWHM,
        'mu_err': errors[1],
        'sigma_err': errors[2],
        'amp_err': errors[0],
        'A': A,
        'bg': bg,
        'pcov': pcov,
        'region_ch': region_ch,
        'region_counts': region_counts,
        'fitted_curve': gaussian_with_background(region_ch, A, mu, sigma, bg)
    }


"""
UNCERTAINTIES
"""

def propagate_energy_uncertainty(channel, channel_err, m, b, m_err, b_err):
    """Propagate uncertainty from channel to energy"""
    energy_err = np.sqrt(
        (m * channel_err) ** 2 +
        (channel * m_err) ** 2 +
        b_err ** 2
    )
    return energy_err

"""
CALIBRATION
"""


def calibrate(yaml_file):
    with open(yaml_file) as f:
        content = f.read()
    peak_dict = yaml.safe_load(content)

    x = np.array(list(peak_dict['Peaks'].values())).flatten()
    y = np.array(list(peak_dict['Energies'].values())).flatten()

    coeffs, pcov = np.polyfit(x, y, 1, cov=True)
    m, b = coeffs[0], coeffs[1]
    m_uncert, b_uncert = np.sqrt(np.diag(pcov))

    return x, y, coeffs, m, b, m_uncert, b_uncert


def energy_calibration_equation(c1, channel, c0):
    return c1 * channel + c0

"""
PLOT FUNCTIONS
"""

def plot_spectrum_with_fit(data, title, peak_channel, window=50):
    """Plot spectrum with Gaussian fit"""
    fit_result = fit_peak(data['channels'], data['counts'], peak_channel, window)
    plt.figure(figsize=(12, 7))

    plt.plot(data['channels'], data['counts'], 'b-', linewidth=1, label='Data', alpha=0.7)
    plt.plot(fit_result['region_ch'], fit_result['fitted_curve'], 'r-',
             linewidth=2.5, label='Gaussian Fit')
    plt.axvspan(fit_result['region_ch'][0], fit_result['region_ch'][-1],
                alpha=0.1, color='yellow', label='Fit Region')

    plt.xlabel('Channel', fontsize=12)
    plt.ylabel('Counts', fontsize=12)
    plt.title(title, fontsize=14, fontweight='bold')
    plt.legend(loc='upper right', fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    return fit_result


def plot_spectrum_with_fit_energy(energy, y, title, peak_channel, window=50):
    """Plot spectrum with Gaussian fit"""
    fit_result = fit_peak(energy, y, peak_channel, window)
    plt.figure(figsize=(12, 7))

    plt.plot(energy, y, 'b-', linewidth=1, label='Data', alpha=0.7)
    plt.plot(fit_result['region_ch'], fit_result['fitted_curve'], 'r-',
             linewidth=2.5, label='Gaussian Fit')
    plt.axvspan(fit_result['region_ch'][0], fit_result['region_ch'][-1],
                alpha=0.1, color='yellow', label='Fit Region')

    plt.xlabel('Energy(Kev)', fontsize=12)
    plt.ylabel('Counts', fontsize=12)
    plt.title(title, fontsize=14, fontweight='bold')
    plt.legend(loc='upper right', fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

    return fit_result

def plot_calibration(x, y, coeffs):
    x_fit = np.linspace(x.min(), x.max(), 100)
    y_fit = np.polyval(coeffs, x_fit)
    plt.scatter(x, y, label='Data')
    plt.plot(x_fit, y_fit, 'r--', label='Fit')
    plt.xlabel('Channel')
    plt.ylabel('Energy (keV)')
    plt.title('Energy vs Channel')
    plt.legend()
    plt.show()

"""
ENERGY RESOLUTION
"""

def calculate_resolution(isotope_names, fits, energy_fits, m, b, known_energies, use_channel_conversion=None):
    """
    Calculate energy resolution for multiple isotopes

    Args:
        isotope_names: list of isotope names
        fits: list of channel-space fit results
        energy_fits: list of energy-space fit results
        m, b: calibration parameters
        known_energies: list of known energies for each isotope
        use_channel_conversion: list of isotope names to use channel conversion (e.g., ['CO60'])

    Returns:
        list of dicts with FWHM_energy, resolution, energy
    """
    if use_channel_conversion is None:
        use_channel_conversion = []

    resolutions = []

    for i, name in enumerate(isotope_names):
        if name in use_channel_conversion:
            # Use channel-space fit (more reliable for noisy peaks)
            FWHM_energy = m * fits[i]['FWHM']
            fitted_energy = m * fits[i]['mu'] + b
            print(f"  {name}: Using channel-space conversion")
        else:
            # Use energy-space fit
            FWHM_energy = energy_fits[i]['FWHM']
            fitted_energy = energy_fits[i]['mu']

        resolution = (FWHM_energy / fitted_energy) * 100

        resolutions.append({
            'FWHM_energy': FWHM_energy,
            'resolution': resolution,
            'energy': known_energies[i]
        })

        print(f"{name}: FWHM = {FWHM_energy:.2f} keV, Resolution = {resolution:.2f}%")

    return resolutions


def plot_resolution_vs_energy(resolutions, detector_name='NaI'):
    """
    Plot energy resolution as a function of energy

    Args:
        resolutions: list of dicts with 'energy', 'resolution', and optionally 'resolution_err'
        detector_name: name of detector
    """
    energies = [r['energy'] for r in resolutions]
    resolution_vals = [r['resolution'] for r in resolutions]

    # Check if resolution_err exists, if not, don't plot error bars
    has_errors = 'resolution_err' in resolutions[0]

    plt.figure(figsize=(10, 6))

    if has_errors:
        resolution_errs = [r['resolution_err'] for r in resolutions]
        plt.errorbar(energies, resolution_vals, yerr=resolution_errs,
                     fmt='o', markersize=8, capsize=5, label='Data')
    else:
        plt.plot(energies, resolution_vals, 'o', markersize=8, label='Data')

    plt.xlabel('Energy (keV)', fontsize=12)
    plt.ylabel('Energy Resolution (%)', fontsize=12)
    plt.title(f'{detector_name} Energy Resolution vs Energy', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()


def fit_resolution_curve(resolutions):
    """
    Fit the resolution curve: R² = aE⁻² + bE⁻¹ + c
    Plots on both linear and logarithmic scales

    Args:
        resolutions: list of dicts with 'energy' and 'resolution'

    Returns:
        fitted parameters a, b, c, errors
    """
    energies = np.array([r['energy'] for r in resolutions])
    resolution_vals = np.array([r['resolution'] for r in resolutions])

    # R² = a/E² + b/E + c
    def res_squared_model(E, a, b, c):
        return a / E ** 2 + b / E + c

    # Fit R² vs E
    R_squared = (resolution_vals / 100) ** 2  # Convert % to fraction

    popt, pcov = curve_fit(res_squared_model, energies, R_squared)
    a, b, c = popt
    errors = np.sqrt(np.diag(pcov))

    # print("\n=== Resolution Curve Fit ===")
    # print(f"R² = {a:.6f}/E² + {b:.6f}/E + {c:.6f}")
    # print(f"Uncertainties: a={errors[0]:.6f}, b={errors[1]:.6f}, c={errors[2]:.6f}")

    # Generate fit curve (use logspace for smooth log plot)
    E_fit = np.logspace(np.log10(min(energies)), np.log10(max(energies)), 100)
    R_squared_fit = res_squared_model(E_fit, a, b, c)
    R_fit = np.sqrt(np.abs(R_squared_fit)) * 100

    # Check if resolution_err exists
    has_errors = 'resolution_err' in resolutions[0]

    #Log-Log Scale (as recommended)

    plt.figure(figsize=(10, 6))

    if has_errors:
        resolution_errs = [r['resolution_err'] for r in resolutions]
        plt.errorbar(energies, resolution_vals, yerr=resolution_errs,
                     fmt='o', markersize=10, capsize=5, label='Data')
    else:
        plt.loglog(energies, resolution_vals, 'o', markersize=10, label='Data')

    plt.loglog(E_fit, R_fit, 'r-', linewidth=2, label='Fitted Curve')

    plt.xlabel('Energy (keV)', fontsize=12)
    plt.ylabel('Energy Resolution (%)', fontsize=12)
    plt.title('Energy Resolution vs Energy (Log-Log Scale)', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3, which='both')
    plt.legend()
    plt.tight_layout()
    plt.show()

    return a, b, c, errors

def intrinsic(activity,activity_err,diameter ,distance):

    intrinsic_rate = activity * (np.pi*(diameter/2)**2)/(4*np.pi*distance)
    intrinsic_rate_err = activity_err * (np.pi*(diameter/2)**2)/(4*np.pi*distance**2)
    return intrinsic_rate, intrinsic_rate_err



def efficiency_uncertainty(nuclide,counts, energy,peak_energy,peak_counts_err, time,
                          emitted_counts, emitted_counts_err,
                          incident_counts, incident_counts_err):

    idx = np.argmin(np.abs(energy - peak_energy))

# Get the counts at that peak
    peak_counts = counts[idx]
    # Count rate and its uncertainty
    count_rate = peak_counts / time
    count_rate_err = peak_counts_err / time

    # Absolute efficiency error: σ(ε_abs) = ε_abs * sqrt((σ_rate/rate)^2 + (σ_emitted/emitted)^2)
    abs_eff = count_rate / emitted_counts
    abs_eff_err = abs_eff * np.sqrt(
        (count_rate_err / count_rate)**2 +
        (emitted_counts_err / emitted_counts)**2
    )
    # Intrinsic efficiency error: σ(ε_int) = ε_int * sqrt((σ_rate/rate)^2 + (σ_incident/incident)^2)
    int_eff = count_rate / incident_counts
    int_eff_err = int_eff * np.sqrt(
        (count_rate_err / count_rate)**2 +
        (incident_counts_err / incident_counts)**2
    )
    print(f'The Absolute Efficiency of {nuclide} is {100*abs_eff:.4f}±{abs_eff_err*100:.4f}%')
    print(f'The Intrinsic Efficiency of {nuclide} is {100*int_eff:.4f}±{int_eff_err*100:.4f}%')
    return abs_eff, int_eff, abs_eff_err, int_eff_err