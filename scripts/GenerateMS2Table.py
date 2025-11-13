# Program heavily drawn from read_mzML1.py and find_spectra.py made by Nathan Zhang and Eric Deutsch.
# ChatGPT also gave substantial coding assistance.
# py GenerateMS2Table.py --mzml_file --previous_list --precursor_mz
# py "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\repository\SyntheticPeptideTools\scripts\GenerateMS2Table.py" --mzml_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\250402_mEclipse_QC_ncORF-097.mzML" --precursor_mz 708.3294
# py "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\repository\SyntheticPeptideTools\scripts\GenerateMS2Table.py" --mzml_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\251103_mEclipse_ncORF89-S3.mzML" --precursor_mz 657.31400280985

import os
import argparse
import os.path
import timeit

import csv
import numpy as np

from pyteomics import mzml, auxiliary

class GenerateMS2Table:
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog="GenerateMS2Table",
            description="Writes out CSV table with MS2 spectra with columns including file root, scan no., scan time, and total ion current (TIC)")
        parser.add_argument("--mzml_file", help="Name of mzML file")
        parser.add_argument("--previous_list", help="Previous annotated list (accepts CSV, TSV files)")
        parser.add_argument("--precursor_mz", help="Precursor mz of peptide")

        args = parser.parse_args()

        if args.mzml_file is None or args.mzml_file == "":
            print('ERROR: Parameter --mzml_file must be provided. See --help for more information')
            return

        # tries to open the params file specified
        if not os.path.isfile(args.mzml_file):
            print(f"ERROR: File '{args.mzml_file}' not found or not a file")
            return

        self.mzml_file = args.mzml_file
        self.previous_list = args.previous_list
        self.precursor_mz = float(args.precursor_mz)
        self.spectra = []
        self.start = timeit.default_timer()
        self.stats = { 'counter': 0, 'ms1spectra': 0, 'ms2spectra': 0 }

    def read_mzml(self):
        with mzml.read(self.mzml_file) as reader:
            # iterates through mzml file and collects precursor mz, charge, scan no., time, precursor mass delta, tic,
            for spectrum in reader:
                self.stats['counter'] += 1

                if spectrum['ms level'] == 1:
                    self.stats['ms1spectra'] += 1
                elif spectrum['ms level'] == 2:
                    self.stats['ms2spectra'] += 1

                    precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                    charge = int(spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])
                    scan_number = 1 + spectrum['index']
                    scan_time = spectrum['scanList']['scan'][0]['scan start time']
                    mass_delta = precursor_mz * charge - self.precursor_mz * 2 - 1.00727 * (charge - 2)
                    total_ion_current = np.sum(spectrum['intensity array'])
                    injection_time = spectrum.get('scanList', {}).get('scan', [{}])[0].get('ion injection time', None)

                    spectrum_data = {'file root' : self.mzml_file,
                                     'scan number' : scan_number,
                                     'scan time' : scan_time,
                                     'total ion current' : total_ion_current,
                                     'precursor m/z' : precursor_mz,
                                     'precursor charge' : charge,
                                     'precursor mass delta' : mass_delta,
                                     'injection time' : injection_time
                                     }
                    self.spectra.append(spectrum_data)

    def read_annotation(self):
        if not self.previous_list:
            print("INFO: No previous annotation list provided; skipping annotation merge.")
            for row in self.spectra:
                row.update({"confidence": "", "modification": "", "usi": "", "comments": ""})
            return

        ext = os.path.splitext(self.previous_list)[1].lower()
        d = ',' if ext == '.csv' else '\t'

        annotations = {}
        with open(self.previous_list, newline='') as file:
            reader = csv.DictReader(file, delimiter=d)
            for row in reader:
                if not any(row.values()):
                    continue
                try:
                    scan_number = int(row.get('scan number', '').strip())
                except ValueError:
                    continue
                annotations[scan_number] = {
                    "confidence": row.get("confidence", ""),
                    "modification": row.get("modification", ""),
                    "usi": row.get("usi", ""),
                    "comments": row.get("comments", "")
                }

        blank_annotation_data = {"confidence": "", "modification": "", "usi": "", "comments": ""}
        for row in self.spectra:
            scan_number = row['scan number']
            row.update(annotations.get(scan_number, blank_annotation_data))

    # creates new merged csv file
    def write_csv(self):
        fieldnames = ["file root", "scan number", "injection time", "scan time", "total ion current", "precursor m/z", "precursor charge", "precursor mass delta", "confidence", "modification", "usi", "comments"]
        with open("ms2_table.csv", "w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            for row in self.spectra:
                writer.writerow(row)

def main():

    generate_table = GenerateMS2Table()
    generate_table.read_mzml()
    print(f"INFO: Read {generate_table.stats['counter']} spectra from {generate_table.mzml_file}")
    print(f"The number of ms1spectra is {generate_table.stats['ms1spectra']}")
    print(f"The number of ms2spectra is {generate_table.stats['ms2spectra']}")
    end = timeit.default_timer()
    print(f"Elapsed time to process spectra: {end - generate_table.start}")
    print(f"INFO: Processed {generate_table.stats['counter'] / (end - generate_table.start)} spectra per second")

    generate_table.read_annotation()
    print(f"INFO: Merged {generate_table.previous_list} and {generate_table.mzml_file} data.")
    generate_table.write_csv()
    print(f"INFO: Generated MS2 table. Length of file: {len(generate_table.spectra)}. Output file: 'ms2_table.csv'")

    final_end = timeit.default_timer()
    print(f"Total elapsed time: {final_end - generate_table.start}")

if __name__ == '__main__':
    main()
