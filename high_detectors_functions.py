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
        # PMCA text format - CdTe detector
        metadata = {
            'filename': filepath.split('/')[-1],
            'detector': 'CdTe',
            'live_time': None,
            'real_time': None,
        }

        counts = []
        IS_DATA = False

        # Open with latin-1 encoding (or just 'r' like file_parser)
        with open(filepath, "r", encoding='latin-1') as file:
            for line in file:
                line = line.strip()

                # Data section
                if line.startswith('<<DATA>>'):
                    IS_DATA = True
                    continue
                elif line.startswith('REAL_TIME'):
                    # Use the file_parser approach
                    metadata['real_time'] = float(line[12:])  # Skip "REAL_TIME - "
                    IS_DATA = False
                elif line.startswith('LIVE_TIME'):
                    # file_parser doesn't have this, but keep it
                    metadata['live_time'] = float(line[12:])  # Skip "LIVE_TIME - "
                    IS_DATA = False
                elif line.startswith('START_TIME'):
                    # Extract date if needed (file_parser uses this for DATE_MEAS)
                    # metadata['date'] = line[13:]  # Uncomment if you want date
                    IS_DATA = False
                elif line.startswith("<<END>>"):
                    IS_DATA = False
                    break

                # Collect data (file_parser uses int, not float)
                if IS_DATA:
                    try:
                        counts.append(int(line))  # Changed from float to int
                    except Exception:
                        continue

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

    # ASCII format - .spe files (NaI/BGO) - KEEP YOUR EXISTING CODE
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
        'FWHM_err' : errors[2]*2.355,
        'mu_err': errors[1],
        'sigma_err': errors[2],
        'amp_err': errors[0],
        'A': A,
        'amp_err': errors[0],
        'bg': bg,
        'pcov': pcov,
        'region_ch': region_ch,
        'region_counts': region_counts,
        'fitted_curve': gaussian_with_background(region_ch, A, mu, sigma, bg)
    }
"""
CALIBRATION
"""
def calibrate(yaml_file):
    with open(yaml_file) as f:
        content = f.read()
    peak_dict = yaml.safe_load(content)

    x = np.array(list(peak_dict['Peaks'].values())).flatten()
    y = np.array(list(peak_dict['Energies'].values())).flatten()

    # Sort to ensure proper order
    sort_idx = np.argsort(x)
    x = x[sort_idx]
    y = y[sort_idx]

    # Handle 2-point fit specially
    if len(x) == 2:
        # For exactly 2 points, calculate line manually, for cdte
        m = (y[1] - y[0]) / (x[1] - x[0])
        b = y[0] - m * x[0]
        coeffs = np.array([m, b])
        # Rough uncertainty estimates
        m_uncert = 0.01
        b_uncert = 1.0
    else:
        # 3+ points: use polyfit with covariance
        coeffs, pcov = np.polyfit(x, y, 1, cov=True)
        m, b = coeffs[0], coeffs[1]
        m_uncert, b_uncert = np.sqrt(np.diag(pcov))

    return x, y, coeffs, m, b, m_uncert, b_uncert


def energy_calibration_equation(c1, channel, c0):
    return c1 * channel + c0

def propagate_energy_uncertainty(channel, channel_err, m, b, m_err, b_err):
    """Propagate uncertainty from channel to energy"""
    energy_err = np.sqrt(
        (m * channel_err) ** 2 +
        (channel * m_err) ** 2 +
        b_err ** 2
    )
    return energy_err
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
SOURCE ACTIVITY
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


def plot_resolution_vs_energy(resolutions, detector_name):
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

def fit_intrinsic_efficiency(energies, intrinsic_efficiencies, degree=1):
    # Convert to numpy arrays
    energies = np.array(energies)
    intrinsic_efficiencies = np.array(intrinsic_efficiencies)

    # Take natural logarithms
    log_E = np.log(energies)
    log_eff = np.log(intrinsic_efficiencies)

    # Calculate weights if errors are provided
    
    coeffs, pcov = np.polyfit(log_E, log_eff, degree, cov=True)

    # Create fit function
    def fit_func(E):
        """Evaluate fitted efficiency at energy E (in keV)"""
        ln_E = np.log(E)
        ln_eff = np.polyval(coeffs, ln_E)
        return np.exp(ln_eff)

    # Print fit parameters
    if degree == 2:
        c, b, a = coeffs
       
    elif degree == 1:
        b, a = coeffs
        

    return {
        'coeffs': coeffs,
        'pcov': pcov,
        'log_E': log_E,
        'log_eff': log_eff,
        'fit_func': fit_func
    }

def plot_intrinsic_efficiency(energies, intrinsic_efficiencies, intrinsic_eff_errors=None,
                               fit_result=None, title='Intrinsic Peak Efficiency vs Energy'):
    
    fig, ax = plt.subplots(figsize=(10, 7))

    # Plot data points
    if intrinsic_eff_errors is not None:
        ax.errorbar(energies, intrinsic_efficiencies, yerr=intrinsic_eff_errors,
                   fmt='o', markersize=8, capsize=5, label='Data', linewidth=2)
    else:
        ax.plot(energies, intrinsic_efficiencies, 'o', markersize=8, label='Data')

    # Plot fit if provided
    if fit_result is not None:
        E_fit = np.logspace(np.log10(min(energies)*0.8), np.log10(max(energies)*1.2), 200)
        eff_fit = fit_result['fit_func'](E_fit)
        ax.plot(E_fit, eff_fit, 'r-', linewidth=2.5, label='Logarithmic Fit')

    # Set log scale on both axes
    ax.set_xscale('log')
    ax.set_yscale('log')

    # Labels and formatting
    ax.set_xlabel('Energy (keV)', fontsize=13)
    ax.set_ylabel('Intrinsic Peak Efficiency', fontsize=13)
    ax.set_title(title, fontsize=15, fontweight='bold')
    ax.legend(loc='best', fontsize=11)
    ax.grid(True, which='both', alpha=0.3)

    plt.tight_layout()
    plt.show()

"""
Off-Axis Response
"""
def fit_off_axis_response(base_path, isotope_prefix, background_data, peak_channel,
                         window=50, on_axis_file=None):

    # Get files - check for both .Spe and .mca extensions
    files = glob.glob(f"{base_path}/{isotope_prefix}*.Spe")
    files += glob.glob(f"{base_path}/{isotope_prefix}*.mca")

    # Parse angles
    angle_list = []
    for f in files:
        filename = f.split('/')[-1]  # Get just the filename
        # Remove prefix and extension (.Spe or .mca) to get angle string
        angle_str = filename.replace(isotope_prefix, '').replace('.Spe', '').replace('.mca', '')
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
            time = data_no_bg['real_time']
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
            time = data_no_bg['real_time']
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


def plot_off_axis_response(angular_data, isotope_name, energy_kev):
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

def fit_off_axis_response_FWHM(base_path, isotope_prefix, background_data, peak_channel,
                         window=50, on_axis_file=None):

    # Get files - check for both .Spe and .mca extensions
    files = glob.glob(f"{base_path}/{isotope_prefix}*.Spe")
    files += glob.glob(f"{base_path}/{isotope_prefix}*.mca")

    # Parse angles
    angle_list = []
    for f in files:
        filename = f.split('/')[-1]  # Get just the filename
        # Remove prefix and extension (.Spe or .mca) to get angle string
        angle_str = filename.replace(isotope_prefix, '').replace('.Spe', '').replace('.mca', '')
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
            time = data_no_bg['real_time']
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
            time = data_no_bg['real_time']
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


def plot_off_axis_response_FWHM(angular_data, isotope_name, energy_kev):
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