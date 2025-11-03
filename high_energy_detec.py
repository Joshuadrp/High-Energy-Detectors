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
plt.figure(figsize=(10, 6))
plt.plot(cs137['channels'], cs137['counts'])
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.title('CS-137 Spectrum')
plt.grid(True)
plt.show()


# Plot Co-60 to find peaks
plt.figure(figsize=(10, 6))
plt.plot(co60['channels'], co60['counts'])
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.title('Co-60 Spectrum')
plt.grid(True)
plt.show()

# Plot Am-241
plt.figure(figsize=(10, 6))
plt.plot(am241['channels'], am241['counts'])
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.title('Am-241 Spectrum')
plt.grid(True)
plt.show()

# Plot Ba-133
plt.figure(figsize=(10, 6))
plt.plot(ba133['channels'], ba133['counts'])
plt.xlabel('Channel')
plt.ylabel('Counts')
plt.title('Ba-133 Spectrum')
plt.grid(True)
plt.show()

