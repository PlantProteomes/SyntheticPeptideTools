# py find_precursor_intensity.py --mzml_file --excel_file --window_size 10 --output 
# py find_precursor_intensity.py --mzml_file C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\NEW_250402_mEclipse_QC_ncORF-089.mzML --excel_file C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\ms2_scans_089.xlsx --window_size 10 --output C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\scan_2771.csv
# py find_precursor_intensity.py --mzml_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\250402_mEclipse_QC_ncORF-055.mzML" --excel_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\new_peptide\scan_precursor_055.xlsx" --window_size 10 --output C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\new_peptide\estimated_precursors_055(NEW).csv
# excel file must have headers "scannumber" and "precursor"


import os
import argparse
import pandas as pd
from pyteomics import mzml
import csv

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

def main():
    parser = argparse.ArgumentParser(description='Estimate precursor m/z from MS1 scans around MS2 scans')
    parser.add_argument('--mzml_file', required=True, help='Input mzML file')
    parser.add_argument('--excel_file', required=True, help='Excel file containing columns "scannumber" and "precursor"')
    parser.add_argument('--window_size', type=int, default=10, help='Number of MS1 scans before/after MS2')
    parser.add_argument('--output', default='estimated_precursor_values.csv', help='Output CSV file')
    args = parser.parse_args()

    # Load Excel file
    if not os.path.isfile(args.excel_file):
        print(f"ERROR: Excel file {args.excel_file} not found")
        return
    df = pd.read_excel(args.excel_file)

    if 'scannumber' not in df.columns or 'precursor' not in df.columns:
        raise ValueError("Excel file must contain columns named 'scannumber' and 'precursor'")

    ms2_list = df['scannumber'].dropna().astype(int).tolist()
    precursor_guess_list = df['precursor'].dropna().astype(float).tolist()

    if len(ms2_list) != len(precursor_guess_list):
        raise ValueError("Number of scannumber and precursor entries must match")

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
                scan_number = int(spectrum['id'].split('=')[-1]) if 'id' in spectrum else None
                ms2_scans.append({
                    'ScanNumber': scan_number,
                    'ScanTime': spectrum['scanList']['scan'][0]['scan start time'] if 'scanList' in spectrum else None
                })

    ms1_count = len(ms1_scans)
    results = []

    # Step 2: Process each MS2 scan
    for ms2_scan, guess_mz in zip(ms2_list, precursor_guess_list):
        ms2 = next((s for s in ms2_scans if s['ScanNumber'] == ms2_scan), None)

        if ms2 is not None:
            ms2_idx = next((i for i, s in enumerate(ms1_scans) if s['ScanTime'] >= ms2['ScanTime']), None)
            if ms2_idx is None:
                print(f"WARNING: Could not find MS1 scans around MS2 scan {ms2_scan}, skipping")
                continue
        else:
            # NEW: If not in MS2, treat as MS1 apex scan
            print(f"INFO: Scan {ms2_scan} not found in MS2, treating as apex MS1 scan")
            ms2_idx = next((i for i, s in enumerate(ms1_scans) if s['ScanNumber'] == ms2_scan), None)
            if ms2_idx is None:
                print(f"WARNING: Scan {ms2_scan} not found in MS1, skipping")
                continue
        start = max(ms2_idx - args.window_size, 0)
        end = min(ms2_idx + args.window_size + 1, ms1_count)
        ms1_window = ms1_scans[start:end]

        mz_values = []
        intensity_values = []
        summed_intensity = 0.0

        for scan in ms1_window:
            mz, intensity = find_peak_in_scan(scan, guess_mz)
            mz_values.append(mz if mz is not None else "NA")
            intensity_values.append(intensity)
            summed_intensity += intensity

        while len(mz_values) < 20:
            mz_values.append("NA")
            intensity_values.append(0.0)
        mz_values = mz_values[:20]
        intensity_values = intensity_values[:20]

        row = {
            'MS2Scan': ms2_scan,
            'Instrument_m/z': guess_mz,
            'SummedIntensity': summed_intensity
        }
        for i in range(20):
            row[f'mz_{i+1}'] = mz_values[i]
            row[f'int_{i+1}'] = intensity_values[i]
        results.append(row)

    # Step 3: Write CSV
    fieldnames = ['MS2Scan', 'Instrument_m/z', 'SummedIntensity'] + \
                 [f'mz_{i+1}' for i in range(20)] + [f'int_{i+1}' for i in range(20)]
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)

    print(f"INFO: Wrote estimated precursor m/z and intensities for {len(results)} MS2 scans to {args.output}")

if __name__ == "__main__":
    main()