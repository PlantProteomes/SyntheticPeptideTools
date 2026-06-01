from models import Scan
from models import Peptide, Modification
import numpy as np
from constants import ppm
from matplotlib import pyplot as plt

def initialize_peaks(spectrum: Scan, sequence: Peptide, modification: Modification, tolerance):
    sorted_mz_array = np.sort(spectrum.mz_array)
    sorted_intensity_array = spectrum.intensity_array[np.argsort(spectrum.mz_array)]
    testing_sequence = Peptide(sequence.raw_sequence, [Modification(m.position, m.delta, m.name, m.is_labile) for m in sequence.modifications])
    testing_modification = Modification(modification.position, modification.delta, modification.name, modification.is_labile)
    testing_sequence.modifications.append(testing_modification)
    fragments = testing_sequence.fragments()
    all_ions = [ion for ion_list in fragments.values() for ion in ion_list]
    found_peaks = []
    relevant_ions = []
    intensities = []

    for ion in all_ions:
        left = np.searchsorted(sorted_mz_array, ion.mz - ppm(ion.mz, tolerance), side="left")
        right = np.searchsorted(sorted_mz_array, ion.mz + ppm(ion.mz, tolerance), side="right")
        mz_slice = sorted_mz_array[left:right]
        intensity_slice = sorted_intensity_array[left:right]
        if len(mz_slice) == 0:
            continue
        best_index = np.argmax(intensity_slice) # registers tallest peak in window
        found_peaks.append(mz_slice[best_index])
        relevant_ions.append(ion)
        intensities.append(intensity_slice[best_index])

    found_peaks = np.array(found_peaks)
    intensities = np.array(intensities)
    weights = intensities / intensities.sum()
    return testing_sequence, testing_modification, found_peaks, relevant_ions, weights

def calculate_stdev(spectrum: Scan, sequence: Peptide, modification: Modification, tolerance):
    testing_sequence, testing_modification, found_peaks, relevant_ions, weights = initialize_peaks(spectrum, sequence, modification, tolerance)
    if len(found_peaks) == 0:
        raise ValueError("no peaks found")

    base_delta = testing_modification.delta
    results = []

    for i in range(-2000, 2000):
        testing_modification.delta = base_delta + i / 1000000
        expected_ion_masses = np.array([ion.mz for ion in relevant_ions])
        residuals = found_peaks - expected_ion_masses
        weighted_mean = np.dot(weights, residuals)
        weighted_variance = np.dot(weights, np.square(residuals - weighted_mean))
        stdev = np.sqrt(weighted_variance)
        results.append((base_delta + i / 1000000, stdev))

    return results

def plot_stdev(output, mass_delta, results, scan_number):
    plt.xlabel("Mass delta difference (Da)")
    plt.ylabel("Standard deviation of residuals")
    best_delta, best_stdev = min(results, key=lambda x: x[1])
    xlist = []
    ylist = []
    bounds = [(None, None), (None, None)]
    prev_y = 0
    prev_x = 0
    for xval, yval in results:
        xlist.append(xval - mass_delta)
        ylist.append(yval)
        if yval < 1.1*best_stdev <= prev_y:
            bounds[0] = (prev_x, prev_y)
        if yval >= 1.1*best_stdev > prev_y:
            bounds[1] = (xval - mass_delta, yval)
        prev_y = yval
        prev_x = xval - mass_delta
    plt.plot(xlist, ylist)
    plt.plot(best_delta - mass_delta, best_stdev, marker='*', color='green')
    plt.annotate(text=f"{best_delta:.6f}", xy=(best_delta - mass_delta, best_stdev), xytext=(0,-3), textcoords='offset points', ha='center', va='top')
    print(bounds)
    plt.annotate(text=f"{bounds[0][0]:.4f}", xy=(bounds[0]), xytext=(-3,-3), textcoords='offset points', ha='right', va='top')
    plt.annotate(text=f"{bounds[1][0]:.4f}", xy=(bounds[1]), xytext=(3,-3), textcoords='offset points', ha='left', va='top')
    plt.axhline(1.10 * best_stdev, linestyle="--")
    print("best stdev is " + str(best_stdev))
    plt.axvline(0, linestyle="--", lw = 0.5)
    plt.xlim(-0.002, 0.002)
    plt.title(f"Scan {scan_number}")
    plt.tight_layout()
    plt.savefig(f"C:\\Users\\carol\\OneDrive - Bellevue School District\\2025-2026\\6-7 Internship\\python_testing\\ncORF-89\\temp_directory\\new_stdev_{output}1")
