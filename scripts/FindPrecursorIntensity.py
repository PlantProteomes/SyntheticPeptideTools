# py find_precursor_intensity.py --mzml_file --window_size 10 --output 
# py FindPrecursorIntensity.py --mzml_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\251103_mEclipse_ncORF89-S1.mzML" --window_size 10 --output "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\peptide_089\089_s1\estimated_precursors_089_s1_TEST.csv"

import os
import argparse
import pandas as pd
from pyteomics import mzml
import csv
import numpy as np

TOLERANCE = 0.003

def find_peak_in_scan(scan, guess_mz, tolerance=TOLERANCE):
    mz_array = scan['m/z array']
    intensity_array = scan['intensity array']
    closest_mz = None
    closest_intensity = 0.0

    for mz, intensity in zip(mz_array, intensity_array):
        if abs(mz - guess_mz) <= tolerance:
            if closest_mz is None or abs(mz - guess_mz) < abs(closest_mz - guess_mz):
                closest_mz = mz
                closest_intensity = intensity
    return closest_mz, closest_intensity

def get_precursor_intensity(ms1_scans,ms2_scan_time,precursor_mz,window_size=10):
    ms2_inx = next((i for i, s in enumerate(ms1_scans) if s['ScanTime'] >= ms2_scan_time), None)
    if ms2_inx is None:
        return None, None, None, None
    
    start = max(ms2_inx - window_size, 0)
    end = min(ms2_inx + window_size + 1, len(ms1_scans))
    ms1_window = ms1_scans[start:end]

    mz_values = []
    intensity_values = []
    summed_intensity = 0.0

    for scan in ms1_window:
        mz, intensity = find_peak_in_scan(scan, precursor_mz)
        mz_values.append(mz if mz is not None else "NA")
        intensity_values.append(intensity)
        summed_intensity += intensity

    while len(mz_values) < 20:
        mz_values.append("NA")
        intensity_values.append(0.0)
    
    mz_values = mz_values[:20]
    intensity_values = intensity_values[:20]

    max_intensity = max(intensity_values)

    return {
        'mz_values': mz_values,
        'intensity_values': intensity_values,
        'summed_intensity': summed_intensity,
        'max_intensity': max_intensity
    }

def process_ms1_window(ms1_scans, center_idx, guess_mz, window_size):
    start = max(center_idx - window_size, 0)
    end = min(center_idx + window_size + 1, len(ms1_scans))
    ms1_window = ms1_scans[start:end]

    mz_values = []
    intensity_values = []
    summed_intensity = 0.0

    for scan in ms1_window:
        mz, intensity = find_peak_in_scan(scan, guess_mz)
        mz_values.append(mz if mz is not None else "NA")
        intensity_values.append(intensity)
        summed_intensity += intensity

    # Ensure exactly 20 points
    while len(mz_values) < 20:
        mz_values.append("NA")
        intensity_values.append(0.0)

    mz_values = mz_values[:20]
    intensity_values = intensity_values[:20]

    max_intensity = max(intensity_values) if intensity_values else 0.0

    return mz_values, intensity_values, summed_intensity, max_intensity

def main():
    parser = argparse.ArgumentParser(description='Estimate precursor m/z from MS1 scans around MS2 scans')
    parser.add_argument('--mzml_file', required=True, help='Input mzML file')
    parser.add_argument('--window_size', type=int, default=10, help='Number of MS1 scans before/after MS2')
    parser.add_argument('--output', default='estimated_precursor_values.csv', help='Output CSV file')
    args = parser.parse_args()

    if not os.path.isfile(args.mzml_file):
        print(f"ERROR: mzML file {args.mzml_file} not found")
        return

    # Step 1: Read all spectra and store MS1 and MS2 scans
    ms1_scans = []
    ms2_scans = []

    with mzml.read(args.mzml_file) as reader:
        for spectrum in reader:
            if spectrum['ms level'] == 1:
                ms1_scans.append({
                    'ScanNumber': int(spectrum['id'].split('=')[-1]) if 'id' in spectrum else None,
                    'ScanTime': spectrum['scanList']['scan'][0]['scan start time'] if 'scanList' in spectrum else None,
                    'm/z array': spectrum['m/z array'],
                    'intensity array': spectrum['intensity array']
                })
            elif spectrum['ms level'] == 2:
                # Extract precursor m/z from mzML
                precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                scan_number = int(spectrum['id'].split('=')[-1]) if 'id' in spectrum else None
                ms2_scans.append({
                    'ScanNumber': scan_number,
                    'ScanTime': spectrum['scanList']['scan'][0]['scan start time'] if 'scanList' in spectrum else None,
                    'Precursor_m/z': precursor_mz,
                    'Spectrum': spectrum
                })

    ms1_count = len(ms1_scans)
    results = []

    # Process MS1 scans to find apex
    apex_scan = max(ms1_scans, key=lambda s: sum(s['intensity array']))
    apex_idx = ms1_scans.index(apex_scan)
    apex_mz_array = apex_scan['m/z array']
    apex_int_array = apex_scan['intensity array']
    guess_mz_apex = apex_mz_array[np.argmax(apex_int_array)]

    mz_values, intensity_values, summed_intensity, max_intensity = process_ms1_window(
        ms1_scans, apex_idx, guess_mz_apex, args.window_size
    )

    apex_row = {
        'MS2Scan': apex_scan['ScanNumber'],
        'SummedIntensity': summed_intensity,
        'Maximum Precursor Intensity': max_intensity,
        'MS2 TIC': 0
    }
    for i in range(20):
        apex_row[f'mz_{i+1}'] = mz_values[i]
        apex_row[f'int_{i+1}'] = intensity_values[i]

    results.append(apex_row)

    # Process each MS2 scan
    for ms2 in ms2_scans:
        ms2_scan = ms2['ScanNumber']
        spectrum = ms2['Spectrum']
        guess_mz = ms2['Precursor_m/z']
        total_ion_current = np.sum(spectrum['intensity array'])

        features = get_precursor_intensity(
            ms1_scans,
            ms2['ScanTime'],
            guess_mz,
            args.window_size
        )

        if features is None:
            continue

        row = {
            'MS2Scan': ms2_scan,    
            'SummedIntensity': features['summed_intensity'],
            'Maximum Precursor Intensity': features['max_intensity'],
            'MS2 TIC': total_ion_current
        }

        for i in range(20):
            row[f'mz_{i+1}'] = features['mz_values'][i]
            row[f'int_{i+1}'] = features['intensity_values'][i]

        results.append(row)

    # Step 3: Write CSV
    apex = results[0]
    real_ms2_rows = results[1:]
    real_ms2_rows.sort(key=lambda x: x['MS2 TIC'], reverse=True)
    results = [apex] + real_ms2_rows

    fieldnames = ['MS2Scan', 'SummedIntensity', 'MS2 TIC', 'Maximum Precursor Intensity',''] + \
                 [f'mz_{i+1}' for i in range(20)] + [f'int_{i+1}' for i in range(20)]
    
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"INFO: Wrote estimated precursor m/z and intensities for {len(results)} MS2 scans to {args.output}")

if __name__ == "__main__":
    main()