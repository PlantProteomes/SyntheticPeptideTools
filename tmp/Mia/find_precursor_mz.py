# example usage: py find_precursor_mz.py --mzml_file input.mzML --ms2_scans 2873,2875,2801 --precursor_mz 676.781037,669.2960205,685.7739788 --tolerance_mz 0.01 --window_size 10 --output summed_precursor.csv
# example usage: py find_precursor_mz.py --mzml_file C:\\Users\\miawc\\OneDrive\\Documents\\ISB_INTERNSHIP\\mia_data\\NEW_250402_mEclipse_QC_ncORF-089.mzML --ms2_scans 2873,2875,2801 --precursor_mz 676.781037,669.2960205,685.7739788 --tolerance_mz 0.01 --window_size 10 --output C:\\Users\\miawc\\OneDrive\\Documents\\ISB_INTERNSHIP\\mia_data\\summed_precursor.csv


import os
import argparse
from pyteomics import mzml
import csv
import numpy as np

def find_estimated_precursor(ms1_window, guess_mz, tolerance):
    # Create a dict to sum intensities per m/z
    summed_intensity_dict = {}
    
    for scan in ms1_window:
        mz_array = scan['m/z array']
        intensity_array = scan['intensity array']
        for mz, intensity in zip(mz_array, intensity_array):
            if abs(mz - guess_mz) <= tolerance:
                mz_key = round(mz, 4)  # round to 4 decimals to bin similar m/z
                if mz_key in summed_intensity_dict:
                    summed_intensity_dict[mz_key] += intensity
                else:
                    summed_intensity_dict[mz_key] = intensity
    
    if not summed_intensity_dict:
        return None, 0.0
    
    # Find m/z with the highest summed intensity
    estimated_mz = max(summed_intensity_dict, key=summed_intensity_dict.get)
    total_intensity = summed_intensity_dict[estimated_mz]
    
    return estimated_mz, total_intensity

def main():
    parser = argparse.ArgumentParser(description='Estimate precursor m/z from MS1 scans around MS2 scans')
    parser.add_argument('--mzml_file', required=True, help='Input mzML file')
    parser.add_argument('--ms2_scans', type=str, required=True, help='Comma-separated list of MS2 scan numbers')
    parser.add_argument('--precursor_mz', type=str, required=True, help='Comma-separated list of guessed precursor m/z values')
    parser.add_argument('--tolerance_mz', type=float, default=0.01, help='Tolerance around guessed m/z')
    parser.add_argument('--window_size', type=int, default=10, help='Number of MS1 scans before/after MS2')
    parser.add_argument('--output', default='estimated_precursor.csv', help='Output CSV file')
    args = parser.parse_args()

    # Convert inputs to lists
    ms2_list = [int(x) for x in args.ms2_scans.split(',')]
    precursor_guess_list = [float(x) for x in args.precursor_mz.split(',')]

    if len(ms2_list) != len(precursor_guess_list):
        raise ValueError("Number of MS2 scans must match number of guessed precursor m/z values")

    if not os.path.isfile(args.mzml_file):
        print(f"ERROR: File {args.mzml_file} not found")
        return

    # Step 1: Read all spectra and store MS1 scans
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
                scan_number = int(spectrum['id'].split('=')[-1]) if 'id' in spectrum else None
                ms2_scans.append({
                    'ScanNumber': scan_number,
                    'ScanTime': spectrum['scanList']['scan'][0]['scan start time'] if 'scanList' in spectrum else None
                })

    ms1_count = len(ms1_scans)
    results = []

    # Step 2: Process each MS2 scan and corresponding guessed m/z
    for ms2_scan, guess_mz in zip(ms2_list, precursor_guess_list):
        ms2 = next((s for s in ms2_scans if s['ScanNumber'] == ms2_scan), None)
        if ms2 is None:
            print(f"WARNING: MS2 scan {ms2_scan} not found, skipping")
            continue

        ms2_idx = next((i for i, s in enumerate(ms1_scans) if s['ScanTime'] >= ms2['ScanTime']), None)
        if ms2_idx is None:
            print(f"WARNING: Could not find MS1 scans around MS2 scan {ms2_scan}, skipping")
            continue

        start = max(ms2_idx - args.window_size, 0)
        end = min(ms2_idx + args.window_size + 1, ms1_count)
        ms1_window = ms1_scans[start:end]

        estimated_mz, summed_intensity = find_estimated_precursor(ms1_window, guess_mz, args.tolerance_mz)

        results.append({
            'MS2Scan': ms2_scan,
            'Instrument_m/z': guess_mz,
            'Calculated_m/z': estimated_mz if estimated_mz else "NA",
            'SummedIntensity': summed_intensity
        })

    # Step 3: Write CSV
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['MS2Scan', 'Instrument_m/z', 'Calculated_m/z', 'SummedIntensity'])
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"INFO: Wrote estimated precursor m/z for {len(results)} MS2 scans to {args.output}")

if __name__ == "__main__":
    main()
