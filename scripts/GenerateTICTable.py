import os
import argparse
import timeit
import csv
from pyteomics import mzml

def main():
    # sets up a module with a description
    parser = argparse.ArgumentParser(description='Generate a TIC table from an mzML file')
    parser.add_argument('--mzml_file', required=True, help='Name of the mzML file to read')
    parser.add_argument('--output', default=None, help='Output CSV/TSV/XLSX file (default: same root as mzML)')
    args = parser.parse_args()

    # Check file exists
    if not os.path.isfile(args.mzml_file):
        print(f"ERROR: File '{args.mzml_file}' not found")
        return

    # Determine file root and output file
    file_root = os.path.splitext(os.path.basename(args.mzml_file))[0]
    output_file = args.output if args.output else f"{file_root}_TIC.csv"

    # Read MS1 spectra and collect TIC info
    t0 = timeit.default_timer()
    tic_data = []

    with mzml.read(args.mzml_file) as reader:
        for spectrum in reader:
            if spectrum['ms level'] == 1:
                scan_number = int(spectrum['id'].split('=')[-1]) if 'id' in spectrum else None
                scan_time = spectrum['scanList']['scan'][0]['scan start time'] if 'scanList' in spectrum else None
                tic = sum(spectrum['intensity array']) if 'intensity array' in spectrum else 0
                tic_data.append({
                    'FileRoot': file_root,
                    'ScanNumber': scan_number,
                    'ScanTime': scan_time,
                    'TIC': tic
                })

    # Write out CSV
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['FileRoot', 'ScanNumber', 'ScanTime', 'TIC'])
        writer.writeheader()
        for row in tic_data:
            writer.writerow(row)

    t1 = timeit.default_timer()
    print(f"INFO: Wrote {len(tic_data)} MS1 spectra to {output_file}")
    print(f"INFO: Elapsed time: {t1 - t0:.2f} seconds")

if __name__ == "__main__":
    main()
