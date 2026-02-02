# py find_precursor_intensity.py --mzml_file --scan_number --precursor_mz --window_size 10 --output 
# py find_ms1_intensity.py --mzml_file C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\NEW_250402_mEclipse_QC_ncORF-089.mzML --scan_number 2771 --precursor_mz 657.31400280985 --window_size 10 --output C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\estimate_peptide_intensity

import argparse
from pyteomics import mzml
import csv

TOLERANCE = 0.003 

def find_peak_in_scan(scan, target_mz, tolerance=TOLERANCE):
    """Find closest m/z to target_mz within tolerance and return intensity"""
    mz_array = scan['m/z array']
    intensity_array = scan['intensity array']
    closest_intensity = 0.0
    closest_mz = None
    for mz, intensity in zip(mz_array, intensity_array):
        if abs(mz - target_mz) <= tolerance:
            if closest_mz is None or abs(mz - target_mz) < abs(closest_mz - target_mz):
                closest_mz = mz
                closest_intensity = intensity
    return closest_mz, closest_intensity

def main():
    parser = argparse.ArgumentParser(description='Estimate precursor intensity at apex MS1 scan')
    parser.add_argument('--mzml_file', required=True, help='Path to mzML file')
    parser.add_argument('--scan_number', type=int, required=True, help='Apex MS1 scan number')
    parser.add_argument('--precursor_mz', type=float, required=True, help='Precursor m/z of peptide')
    parser.add_argument('--window_size', type=int, default=10, help='Number of MS1 scans before/after apex')
    parser.add_argument('--output', required=True, help='Output CSV file')
    args = parser.parse_args()

    # Read all MS1 scans
    ms1_scans = []
    with mzml.read(args.mzml_file) as reader:
        for spectrum in reader:
            if spectrum['ms level'] == 1:
                scan_num = int(spectrum['id'].split('=')[-1]) if 'id' in spectrum else None
                ms1_scans.append({
                    'ScanNumber': scan_num,
                    'm/z array': spectrum['m/z array'],
                    'intensity array': spectrum['intensity array']
                })

    # Find index of apex scan
    apex_idx = next((i for i, s in enumerate(ms1_scans) if s['ScanNumber'] == args.scan_number), None)
    if apex_idx is None:
        print(f"ERROR: Scan {args.scan_number} not found in MS1 scans")
        return

    # Take ±window_size scans around apex
    start = max(apex_idx - args.window_size, 0)
    end = min(apex_idx + args.window_size + 1, len(ms1_scans))
    ms1_window = ms1_scans[start:end]

    mz_values = []
    intensity_values = []
    summed_intensity = 0.0

    for scan in ms1_window:
        mz, intensity = find_peak_in_scan(scan, args.precursor_mz)
        mz_values.append(mz if mz is not None else "NA")
        intensity_values.append(intensity)
        summed_intensity += intensity

    # Fill to 20 values for backward compatibility
    while len(mz_values) < 20:
        mz_values.append("NA")
        intensity_values.append(0.0)
    mz_values = mz_values[:20]
    intensity_values = intensity_values[:20]

    # Prepare row
    row = {
        'MS2Scan': args.scan_number,
        'Instrument_m/z': args.precursor_mz,
        'SummedIntensity': summed_intensity
    }
    for i in range(20):
        row[f'mz_{i+1}'] = mz_values[i]
        row[f'int_{i+1}'] = intensity_values[i]

    # Write CSV
    fieldnames = ['MS2Scan', 'Instrument_m/z', 'SummedIntensity'] + \
                 [f'mz_{i+1}' for i in range(20)] + [f'int_{i+1}' for i in range(20)]
    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow(row)

if __name__ == "__main__":
    main()