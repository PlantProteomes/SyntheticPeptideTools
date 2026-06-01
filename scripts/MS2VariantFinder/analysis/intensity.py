from models import Scan, MSRun
import numpy as np
from constants import ppm

def calculate_precursor_intensity(mz, spectrum: Scan, run: MSRun, tolerance):
    max_precursor_intensity = 0.0
    max_precursor_mz = 0.0
    total_precursor_intensity = 0.0
    intensity_row = []
    mz_row = []

    ms1_idx = run.ms1_spectra.index(spectrum)
    start = max(0, ms1_idx - 10)
    end = min(len(run.ms1_spectra), ms1_idx + 11)

    for i in range(start, end):
        curr_scan = run.ms1_spectra[i]
        sorted_mz_array = np.sort(curr_scan.mz_array)
        sorted_intensity_array = curr_scan.intensity_array[np.argsort(curr_scan.mz_array)]
        left = np.searchsorted(sorted_mz_array, mz - ppm(mz, tolerance), side="left")
        right = np.searchsorted(sorted_mz_array, mz + ppm(mz, tolerance), side="right")
        mz_slice = sorted_mz_array[left:right]
        intensity_slice = sorted_intensity_array[left:right]

        if len(intensity_slice) == 0:
            mz_row.append(None)
            intensity_row.append(None)
            continue

        bp_intensity = np.max(intensity_slice)
        bp_mz = mz_slice[np.argmax(intensity_slice)]
        if bp_intensity > max_precursor_intensity:
            max_precursor_intensity = bp_intensity
            max_precursor_mz = bp_mz
        total_precursor_intensity += bp_intensity
        intensity_row.append(bp_intensity)
        mz_row.append(bp_mz)

    return max_precursor_intensity, total_precursor_intensity, max_precursor_mz, mz_row, intensity_row

def find_max_ms1(run: MSRun, mz, tolerance):
    all_ms1_intensities = []
    all_ms1_mzs = []
    mz_row = []
    intensity_row = []
    best_sn_ratio = 0.0

    for ms1_scan in run.ms1_spectra:
        sorted_mz_array = np.sort(ms1_scan.mz_array)
        sorted_intensity_array = ms1_scan.intensity_array[np.argsort(ms1_scan.mz_array)]
        left = np.searchsorted(sorted_mz_array, mz - ppm(mz, tolerance), side="left")
        right = np.searchsorted(sorted_mz_array, mz + ppm(mz, tolerance), side="right")
        mz_slice = sorted_mz_array[left:right]
        intensity_slice = sorted_intensity_array[left:right]

        signal = np.sum(sorted_intensity_array[-3:-1]) / 2
        noise = np.sum(sorted_intensity_array[0:2]) / 2
        if signal / noise > best_sn_ratio:
            best_sn_ratio = signal / noise

        if len(intensity_slice) == 0: # just puts 0 if nothing found.
            all_ms1_intensities.append(0.0)
            all_ms1_mzs.append(0.0)
            continue

        bp_intensity = np.max(intensity_slice)
        bp_mz = mz_slice[np.argmax(intensity_slice)]
        all_ms1_intensities.append(bp_intensity)
        all_ms1_mzs.append(bp_mz)

    all_ms1_intensities = np.array(all_ms1_intensities)
    all_ms1_mzs = np.array(all_ms1_mzs)
    max_ms1_intensity = np.max(all_ms1_intensities)
    max_ms1_idx = np.argmax(all_ms1_intensities)
    max_ms1_mz = all_ms1_mzs[max_ms1_idx]
    best_ms1_spectrum = run.ms1_spectra[max_ms1_idx]

    start = max(0, int(max_ms1_idx) - 10)
    end = min(len(run.ms1_spectra), int(max_ms1_idx) + 11)
    for i in range(start, end):
        mz_row.append(all_ms1_mzs[i])
        intensity_row.append(all_ms1_intensities[i])

    return best_ms1_spectrum, max_ms1_intensity, max_ms1_mz, mz_row, intensity_row, best_sn_ratio