import numpy as np
import matplotlib.pyplot as plt

# ============================================================================
# STEP 1: DATA LOADING
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
# STEP 2: PEAK FITTING
# ============================================================================

def gaussian(x, A, mu, sigma):
    """Gaussian function"""
    return (A / (sigma * np.sqrt(2 * np.pi))) * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))


def gaussian_with_background(x, A, mu, sigma, bg):
    """Gaussian function with constant background"""
    return gaussian(x, A, mu, sigma) + bg


def plot_spectrum(data, title):
    """Plot a spectrum"""
    plt.figure(figsize=(10, 6))
    plt.plot(data['channels'], data['counts'])
    plt.xlabel('Channel')
    plt.ylabel('Counts')
    plt.title(title)
    plt.grid(True)
    plt.show()

# ============================================================================
# MAIN: FIT ALL PEAKS
# ============================================================================

# Load your data
base_path = "Foreigners/NaITI"

cs137 = load_spe(f"{base_path}/CS137_aligned.Spe")
co60 = load_spe(f"{base_path}/CO60_align.Spe")
am241 = load_spe(f"{base_path}/AM_ALIGN.Spe")
ba133 = load_spe(f"{base_path}/BA_ALIGNED.Spe")

# Plot CS-137 to find peak
plot_spectrum(cs137, f"CS137 aligned")
# Plot Co-60 to find peaks
plot_spectrum(co60, f"C060 aligned")
# Plot Am-241
plot_spectrum(am241, f"AM241 aligned")
# Plot Ba-133
plot_spectrum(ba133, f"BA133 aligned")


