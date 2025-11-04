import yaml
import numpy as np
import matplotlib.pyplot as plt

def calibrate(yaml_file):
    with open(yaml_file) as f:
        content = f.read()
    peak_dict = yaml.safe_load(content)

    x = np.array(list(peak_dict['Peaks'].values())).flatten()
    y = np.array(list(peak_dict['Energies'].values())).flatten()
    coeffs = np.polyfit(x, y, 1)
    m, b = coeffs[0], coeffs[1]

    print(f"Energy = {m:.2f} × Channel + {b:.2f}")
    print(f"Conversion factor: {m:.4f} keV/channel")

    return x, y, coeffs, m, b

def energy_calibration_equation(c1,channel, c0):
    return c1 * channel + c0

# def plot_calibration(x, y, coeffs):
#     x_fit = np.linspace(x.min(), x.max(), 100)
#     y_fit = np.polyval(coeffs, x_fit)
#     plt.scatter(x, y, label='Data')
#     plt.plot(x_fit, y_fit, 'r--', label='Fit')
#     plt.xlabel('Channel')
#     plt.ylabel('Energy (keV)')
#     plt.title('Energy vs Channel')
#     plt.legend()
#     plt.show()