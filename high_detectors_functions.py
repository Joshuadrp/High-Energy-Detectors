import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

"""
DATA LOADING FUNCTION
"""

def load_spe(filepath):
    """Load .spe file and return spectrum with metadata"""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # Initialize metadata
    metadata = {
        'filename': filepath.split('/')[-1],
        'detector': 'NaI',
        'live_time': None,
        'real_time': None,
    }

    # Parse metadata
    for i, line in enumerate(lines):
        if '$MEAS_TIM:' in line and i + 1 < len(lines):
            times = lines[i + 1].strip().split()
            if len(times) >= 2:
                metadata['live_time'] = float(times[0])
                metadata['real_time'] = float(times[1])

    # Parse spectrum data
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

    # Subtract and ensure no negative counts
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
    # Get region around peak
    center_idx = np.argmin(np.abs(channels - peak_channel))
    start = max(0, center_idx - window)
    end = min(len(channels), center_idx + window)

    region_ch = channels[start:end]
    region_counts = counts[start:end]

    # Initial guesses
    A_guess = np.sum(region_counts)
    mu_guess = region_ch[np.argmax(region_counts)]
    sigma_guess = 5.0
    bg_guess = np.min(region_counts)

    # Fit
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
        'mu_error': errors[1],
        'pcov': pcov,  # Add covariance matrix
        'region_ch': region_ch,
        'region_counts': region_counts,
        'fitted_curve': gaussian_with_background(region_ch, A, mu, sigma, bg)
    }

"""
PLOT FIT SPECTRUM
"""


def plot_spectrum_with_fit(x,y, title, peak, window=50):
    """Plot spectrum with improved Gaussian fit overlaid on peak region"""
    fit_result = fit_peak(x, y, peak, window)
    plt.figure(figsize=(12, 7))
    # Plot data
    plt.plot(x,y, 'b-', linewidth=1, label='Data', alpha=0.7)
    # Plot Gaussian fit
    plt.plot(fit_result['region_ch'], fit_result['fitted_curve'], 'r-',
             linewidth=2.5, label='Gaussian Fit')
    # Highlight the fit region
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
