import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import yaml
import pandas as pd
"""
DATA LOADING FUNCTION
"""


def load_spe(filepath):
    """Load .spe file and return spectrum with metadata"""
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
    corrected_counts = signal_counts - background_counts
    corrected_counts = np.maximum(corrected_counts, 0)

    return {
        'channels': signal_data['channels'],
        'counts': corrected_counts,
        'filename': signal_data['filename'],
        'detector': signal_data['detector'],
        'live_time': signal_data['live_time'],
        'real_time': signal_data['real_time'],
    }


def gaussian(x, A, mu, sigma):
    """Gaussian function"""
    return (A / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))


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
    sigma_guess = 5.0
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
        'A': A,
        'bg': bg,
        'mu_err': errors[1],
        'sigma_err': errors[2],
        'A_err': errors[0],
        'bg_err': errors[3],
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
Efficiency Functions
"""

def calc_half_life(nuclide, elap_time=45.88, 
                   isotope_file='efficiency_files/isotope_data(1).yaml',
                   activity_file='efficiency_files/set1184_downstairs(1).dat'):
    
    with open(isotope_file, 'r') as f:
        content = f.read()
    isotopes = [yaml.safe_load('Isotope:' + section) 
                for section in content.split('Isotope:')[1:]]
    
    for isotope in isotopes:
        half_life = isotope['Half-life']
        if isinstance(half_life, str):
            half_life = half_life.replace('years', '').replace('year', '').strip()
            isotope['Half-life'] = float(half_life)
    
    isotope_dict = {isotope['Isotope']: isotope for isotope in isotopes}
    half_life_dict = {iso['Isotope']: iso['Half-life'] for iso in isotopes}
    
    # Also create reversed format dict (60-Co -> Co-60)
    reversed_half_life = {}
    for key in half_life_dict:
        parts = key.split('-')
        reversed_key = f"{parts[1]}-{parts[0]}"
        reversed_half_life[reversed_key] = half_life_dict[key]
    
    df = pd.read_csv(activity_file)
    nuclide_activity = {}
    for index, row in df.iterrows():
        key = str(row.iloc[0]).strip()
        nuclide_activity[key] = row.iloc[3]
    
    # Use reversed format for half-life lookup
    curie_amt = nuclide_activity[nuclide]*(1/2)**(elap_time/reversed_half_life[nuclide])
    photon_amt = curie_amt*37000
    
    print(f'--- {nuclide} ACTIVITY AND HALFLIFE ---')
    print(f'{curie_amt:.2f} is the current activity in uCi')
    print(f'{photon_amt:.2f} is the photons per second given the activity')
    print(f'Half life of {nuclide} is {reversed_half_life[nuclide]} years')

    return photon_amt, nuclide

def intrinsic(activity,diameter,distance):

    intrinsic_rate = activity * (np.pi*(diameter/2)**2)/(4*np.pi*distance)
    
    return intrinsic_rate

def efficiency(nuclide,detected_counts, emiited_counts,peak_energy,isotope_energy, incident_counts):
    peak_energy = peak_energy
    idx = np.argmin(np.abs(isotope_energy - peak_energy))

# Get the counts at that index
    peak_counts = detected_counts[idx]
    
    abs_eff = peak_counts/emiited_counts
    int_eff = peak_counts/incident_counts
    print(f'The Absolute Efficiency of {nuclide} is {100*abs_eff}%')
    print(f'The Intrinsic Efficiency of {nuclide} is {100*int_eff}%')
    return abs_eff, int_eff

