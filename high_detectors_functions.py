import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import yaml
import pandas as pd
import glob
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
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def gaussian_with_background(x, A, mu, sigma,bg):
    """Gaussian with background"""
    return gaussian(x, A, mu, sigma) + bg

def fit_peak(channels, counts, peak_channel, window=50):
    """Fit Gaussian to a peak"""
    center_idx = np.argmin(np.abs(channels - peak_channel))
    start = max(0, center_idx - window)
    end = min(len(channels), center_idx + window)

    region_ch = channels[start:end]
    region_counts = counts[start:end]

    A_guess = A_guess = np.max(region_counts) - np.min(region_counts)
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
        'FWHM_err' : errors[2]*2.355,
        'mu_err': errors[1],
        'A': A,
        'amp_err': errors[0],
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
    nuclide_activity_uncert = {}
    for index, row in df.iterrows():
        key = str(row.iloc[0]).strip()
        nuclide_activity_uncert[key] = row.iloc[4]
    
    curie_amt = nuclide_activity[nuclide]*(1/2)**(elap_time/reversed_half_life[nuclide])
    photon_amt = curie_amt*37000
    
    curie_amt_err = nuclide_activity_uncert[nuclide] *(1/2)**(elap_time/reversed_half_life[nuclide])
    photon_amt_err = curie_amt_err * 37000
    
    print(f'--- {nuclide} ACTIVITY AND HALFLIFE ---')
    print(f'{curie_amt:.2f} is the current activity in uCi')
    print(f'{photon_amt:.2f} is the photons per second given the activity')
    print(f'Photon emission uncertainty: ±{photon_amt_err:.2f} photons/s')
    print(f'Half life of {nuclide} is {reversed_half_life[nuclide]} years')

    return photon_amt, nuclide, photon_amt_err

def intrinsic(activity,activity_err,diameter,distance):

    intrinsic_rate = activity * (np.pi*(diameter/2)**2)/(4*np.pi*distance**2)
    intrinsic_rate_err = activity_err * (np.pi*(diameter/2)**2)/(4*np.pi*distance**2)
    
    return intrinsic_rate, intrinsic_rate_err



def efficiency_uncertainty(nuclide, peak_counts,peak_counts_err, time,
                          emitted_counts, emitted_counts_err,
                          incident_counts, incident_counts_err):

    #idx = np.argmin(np.abs(energy - peak_energy))

# Get the counts at that peak
    #peak_counts = counts[idx]
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


def fit_angular_response(base_path, isotope_prefix, background_data, peak_channel,
                         window=50, on_axis_file=None):

    # Get files
    files = glob.glob(f"{base_path}/{isotope_prefix}*.Spe")

    # Parse angles
    angle_list = []
    for f in files:
        filename = f.split('/')[-1]  # Get just the filename
        # Remove prefix and .Spe to get angle string
        angle_str = filename.replace(isotope_prefix, '').replace('.Spe', '')
        angle = int(angle_str)  # Negative sign is already in the string
        angle_list.append((angle, f))

    angle_list.sort()

    # Fit peaks
    angles = []
    amplitudes = []
    errors = []

    # Add on-axis measurement if provided
    if on_axis_file:
        data = load_spe(on_axis_file)
        data_no_bg = subtract_background(data, background_data)
        try:
            fit_result = fit_peak(data_no_bg['channels'], data_no_bg['counts'],
                                 peak_channel, window)
            time = data_no_bg['live_time']
            count_rate = fit_result['A'] / time
            count_rate_err = fit_result['amp_err'] / time
            angles.append(0)
            amplitudes.append(count_rate)
            errors.append(count_rate_err)
        except Exception as e:
            print(f"Angle    0°: Fit failed - {e}")

    # Fit off-axis measurements
    for angle, filepath in angle_list:
        data = load_spe(filepath)
        data_no_bg = subtract_background(data, background_data)

        try:
            fit_result = fit_peak(data_no_bg['channels'], data_no_bg['counts'],
                                 peak_channel, window)

            # Normalize to counts per second using live_time
            time = data_no_bg['live_time']
            count_rate = fit_result['A'] / time
            count_rate_err = fit_result['amp_err'] / time

            angles.append(angle)
            amplitudes.append(count_rate)
            errors.append(count_rate_err)
        except Exception as e:
            print(f"Angle {angle:+4d}°: Fit failed - {e}")

    # Sort by angle so the plot line connects correctly
    sorted_data = sorted(zip(angles, amplitudes, errors))
    angles = [x[0] for x in sorted_data]
    amplitudes = [x[1] for x in sorted_data]
    errors = [x[2] for x in sorted_data]

    return {'angles': angles, 'amplitudes': amplitudes, 'errors': errors}


def plot_angular_response(angular_data, isotope_name, energy_kev):
    """
    Plot peak amplitude vs angle

    Parameters:
    -----------
    angular_data : dict
        Output from fit_angular_response() with 'angles', 'amplitudes', 'errors'
    isotope_name : str
        Name for plot label (e.g., 'Cs-137')
    energy_kev : float
        Energy in keV for plot label
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.errorbar(list(map(float, angular_data['angles'])), angular_data['amplitudes'],
                yerr=angular_data['errors'],
                fmt='o-', capsize=5, markersize=8, linewidth=2)

    ax.set_xlabel('Angle (degrees)', fontsize=13)
    ax.set_ylabel('Peak Count Rate (counts/s)', fontsize=13)
    ax.set_title(f'{isotope_name} ({energy_kev} keV) Angular Response',
                 fontsize=15, fontweight='bold')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

def fit_angular_response_FWHM(base_path, isotope_prefix, background_data, peak_channel,
                         window=50, on_axis_file=None):

    # Get files
    files = glob.glob(f"{base_path}/{isotope_prefix}*.Spe")

    # Parse angles
    angle_list = []
    for f in files:
        filename = f.split('/')[-1]  # Get just the filename
        # Remove prefix and .Spe to get angle string
        angle_str = filename.replace(isotope_prefix, '').replace('.Spe', '')
        angle = int(angle_str)  # Negative sign is already in the string
        angle_list.append((angle, f))

    angle_list.sort()

    # Fit peaks
    angles = []
    FWHMs = []
    errors = []

    # Add on-axis measurement if provided
    if on_axis_file:
        data = load_spe(on_axis_file)
        data_no_bg = subtract_background(data, background_data)
        try:
            fit_result = fit_peak(data_no_bg['channels'], data_no_bg['counts'],
                                 peak_channel, window)
            time = data_no_bg['live_time']
            FWHM = fit_result['FWHM']
            FWHM_err = fit_result['FWHM_err'] / time
            angles.append(0)
            FWHMs.append(FWHM)
            errors.append(FWHM_err)
        except Exception as e:
            print(f"Angle    0°: Fit failed - {e}")

    # Fit off-axis measurements
    for angle, filepath in angle_list:
        data = load_spe(filepath)
        data_no_bg = subtract_background(data, background_data)

        try:
            fit_result = fit_peak(data_no_bg['channels'], data_no_bg['counts'],
                                 peak_channel, window)

            # Normalize to counts per second using live_time
            time = data_no_bg['live_time']
            FWHM = fit_result['FWHM']
            FWHM_err = fit_result['FWHM_err']

            angles.append(angle)
            FWHMs.append(FWHM)
            errors.append(FWHM_err)
        except Exception as e:
            print(f"Angle {angle:+4d}°: Fit failed - {e}")

    # Sort by angle so the plot line connects correctly
    sorted_data = sorted(zip(angles, FWHMs, errors))
    angles = [x[0] for x in sorted_data]
    FWHMs = [x[1] for x in sorted_data]
    errors = [x[2] for x in sorted_data]

    return {'angles': angles, 'FWHMs': FWHMs, 'errors': errors}


def plot_angular_response_FWHM(angular_data, isotope_name, energy_kev):
    """
    Plot peak amplitude vs angle

    Parameters:
    -----------
    angular_data : dict
        Output from fit_angular_response() with 'angles', 'amplitudes', 'errors'
    isotope_name : str
        Name for plot label (e.g., 'Cs-137')
    energy_kev : float
        Energy in keV for plot label
    """
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.errorbar(list(map(float, angular_data['angles'])), angular_data['FWHMs'],
                yerr=angular_data['errors'],
                fmt='o-', capsize=5, markersize=8, linewidth=2)

    ax.set_xlabel('Angle (degrees)', fontsize=13)
    ax.set_ylabel('FWHM', fontsize=13)
    ax.set_title(f'{isotope_name} ({energy_kev} keV) Angular Response of FWHM',
                 fontsize=15, fontweight='bold')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()