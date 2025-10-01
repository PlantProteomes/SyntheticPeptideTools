# Program heavily drawn from read_mzML1.py and find_spectra.py made by Nathan Zhang and Eric Deutsch.
# ChatGPT also gave substantial coding assistance.

import os
import argparse
import os.path
import timeit
import matplotlib.pyplot as plt
import spectrum_utils.spectrum as sus
import spectrum_utils.plot as sup
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
                    scan_time = spectrum['scanList']['scan'][0]['scan start time'] # does scan time refer to ion injection time or scan start time?
                    mass_delta = precursor_mz * charge - 657.3140 * 2
                    total_ion_current = np.sum(spectrum['intensity array'])

                    spectrum_data = {'file root' : self.mzml_file,
                                     'scan number' : scan_number,
                                     'scan time' : scan_time,
                                     'total ion current' : total_ion_current,
                                     'precursor m/z' : precursor_mz,
                                     'precursor charge' : charge,
                                     'precursor mass delta' : mass_delta}
                    self.spectra.append(spectrum_data)

    def read_annotation(self):
        # defines delimiter
        ext = os.path.splitext(self.previous_list)[1].lower()
        d = ""

        # checks if extension is tsv or csv
        if ext == '.csv':
            d = ','
        elif ext == '.tsv':
            d = '\t'
        else:
            print(f"ERROR: Unsupported file extension '{ext}'")

        # collects confidence, annotations, modifications, and usi from original file, and keys them to scan no.
        annotations = {}
        with open(self.previous_list, newline = '') as file:
            reader = csv.DictReader(file, delimiter = d)
            for row in reader:
                if not any(row.values()):
                    continue

                scan_number = int(row['scan number'])
                annotation_data = {"confidence" : row["confidence"],
                                   "modification" : row["modification"],
                                   "usi" : row["usi"],
                                   "comments" : row["annotation"]}
                annotations[scan_number] = annotation_data

        print("Length of annotated file:", len(annotations))
        print("Length of spectra before merging:", len(self.spectra))

        # merges annotations to final spectra list
        for row in self.spectra:
            scan_number = row['scan number']
            if scan_number in annotations:
                row.update(annotations[scan_number])

    def write_csv(self):
        fieldnames = ["file root", "scan number", "scan time", "total ion current", "precursor m/z", "precursor charge", "precursor mass delta", "confidence", "modification", "usi", "comments"]
        with open("ms2_table.csv", "w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            for row in self.spectra:
                writer.writerow(row)

    # def write_tsv(self):
    # def write_xlsx(self):

def main():
    print("Generating MS2 table. Output file: 'ms2_table.csv'")

    generate_table = GenerateMS2Table()
    generate_table.read_mzml()
    generate_table.read_annotation()

    print("Writing csv. Length of file:", len(generate_table.spectra))
    generate_table.write_csv()
    end = timeit.default_timer()
    print("Elapsed time:", end - generate_table.start)
    print(f"INFO: Read {generate_table.stats['counter']} spectra from {generate_table.mzml_file}")
    print(f"The number of ms1spectra is {generate_table.stats['ms1spectra']}")
    print(f"The number of ms2spectra is {generate_table.stats['ms2spectra']}")
    print(f"INFO: Processed {generate_table.stats['counter'] / (end - generate_table.start)} spectra per second")

if __name__ == '__main__':
    main()