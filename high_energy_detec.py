import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


# ============================================================================
# DATA LOADING
# ============================================================================

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
        'channels': channels,
        'counts': counts,
        'filename': metadata['filename'],
        'detector': metadata['detector'],
        'live_time': metadata['live_time'],
        'real_time': metadata['real_time'],
    }


# ============================================================================
# PEAK FITTING
# ============================================================================

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
        'mu_err': errors[1],
        'region_ch': region_ch,
        'region_counts': region_counts,
        'fitted_curve': gaussian_with_background(region_ch, A, mu, sigma, bg)
    }


def plot_spectrum_with_fit(data, title, peak_channel, window=50):
    """Plot spectrum with Gaussian fit overlaid only on peak region"""
    fit_result = fit_peak(data['channels'], data['counts'], peak_channel, window)

    plt.figure(figsize=(10, 6))
    plt.plot(data['channels'], data['counts'], 'b-', linewidth=1, label='Data')
    plt.plot(fit_result['region_ch'], fit_result['fitted_curve'], 'r-', linewidth=2, label='Gaussian Fit')

    text = f"Peak: {fit_result['mu']:.2f}\nFWHM: {fit_result['FWHM']:.2f}\nσ: {fit_result['sigma']:.2f}"
    plt.text(0.7, 0.95, text, transform=plt.gca().transAxes,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.xlabel('Channel')
    plt.ylabel('Counts')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.show()

    return fit_result


# ============================================================================
# MAIN
# ============================================================================

base_path = "Foreigners/NaITI"

# Load data
cs137 = load_spe(f"{base_path}/CS137_aligned.Spe")
co60 = load_spe(f"{base_path}/CO60_align.Spe")
am241 = load_spe(f"{base_path}/AM_ALIGN.Spe")
ba133 = load_spe(f"{base_path}/BA_ALIGNED.Spe")
background = load_spe(f"{base_path}/BACKGROUND.Spe")

# Subtract background
print("Subtracting background...")
cs137 = subtract_background(cs137, background)
co60 = subtract_background(co60, background)
am241 = subtract_background(am241, background)
ba133 = subtract_background(ba133, background)

# Fit and plot peaks
print("\nCS-137")
cs137_fit = plot_spectrum_with_fit(cs137, title="CS-137", peak_channel=300)

print("\nCo-60")
co60_fit = plot_spectrum_with_fit(co60, "Co-60", peak_channel=100)

print("\nAm-241")
am241_fit = plot_spectrum_with_fit(am241, "Am-241", peak_channel=50)

print("\nBa-133")
ba133_fit = plot_spectrum_with_fit(ba133, "Ba-133", peak_channel=150)