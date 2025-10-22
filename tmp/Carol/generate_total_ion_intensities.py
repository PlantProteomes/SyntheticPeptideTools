# Program references find_spectra.py made by Nathan Zhang and Eric Deutsch.
# Some reference taken from find_precursor_mz.py made by Mia Chen. Some help from ChatGPT was referenced.
# Sample command: --mzml_file "C:\Users\carol\Downloads\250402_mEclipse_QC_ncORF-089.mzML\250402_mEclipse_QC_ncORF-089.mzML" --precursor_mz 676.781037,1313.622732 --scan_numbers 2873,2775 --tolerance 0.01 --window 10

import os
import argparse
import os.path
import re
import timeit
import csv
import numpy as np
from pyteomics import mzml, auxiliary

class TotalIonIntensities:
    """
    Outputs precursor ion intensities from ms1 spectra given a list of target scan numbers and target precursor m/z.
    """
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog="GenerateTotalIonIntensities",
            description="Iterates through MS1 spectra to find scan intensities with specified precursor mz and sums them.")
        parser.add_argument("--mzml_file", help="Name of mzML file")
        parser.add_argument("--precursor_mz", help="Comma-separated list of precursor m/z to select, corresponding to the scan number list")
        parser.add_argument("--tolerance", help="Tolerance of m/z (default 0.01)", default=0.01, type=float)
        parser.add_argument("--scan_numbers", help="Comma-separated list of scan numbers of MS2 scan, corresponding to the precursor m/z list")
        parser.add_argument("--window", help="Number of ms1 scans to iterate through (before and after) (default 10)", default=10, type=int)
        args = parser.parse_args()

        if args.mzml_file is None or args.mzml_file == "":
            print('ERROR: Parameter --mzml_file must be provided. See --help for more information')
            return
        if not os.path.isfile(args.mzml_file):
            print(f"ERROR: File '{args.mzml_file}' not found or not a file")
            return

        if args.precursor_mz is None or args.precursor_mz == "":
            print('ERROR: Parameter --precursor_mz must be provided. See --help for more information')
            return
        if args.scan_numbers is None or args.scan_numbers == "":
            print('ERROR: Parameter --scan_numbers must be provided. See --help for more information')
            return

        self.mzml_file = args.mzml_file
        self.precursor_mz = list(map(float, args.precursor_mz.split(",")))
        self.scan_numbers = list(map(int, args.scan_numbers.split(",")))
        self.tolerance = args.tolerance
        self.start = timeit.default_timer()
        self.window = args.window

    def read_mzml(self, target_scan_number, target_precursor_mz): # writes singular row
        ms1_data = []
        summed_intensity = 0
        index = 0
        mz_dict = {}
        intensity_dict = {}

        with mzml.read(self.mzml_file) as reader:
            for spectrum in reader:
                if spectrum['ms level'] == 1:
                    scan_number = 1 + spectrum['index']
                    ms1_mz_array = spectrum['m/z array']
                    ms1_intensity_array = spectrum['intensity array']
                    ms1_data.append({'scan number' : scan_number, 'm/z array' : ms1_mz_array, 'intensity array' : ms1_intensity_array})

                elif spectrum['ms level'] == 2:
                    if 1 + spectrum['index'] == target_scan_number:
                        precursor_scan_number = int(re.search(r"scan=(\d+)", spectrum['precursorList']['precursor'][0]['spectrumRef']).group(1))

        for spectrum in ms1_data:
            if spectrum['scan number'] == precursor_scan_number:
                index = ms1_data.index(spectrum)
                break

        mz_current_array = ms1_data[index]['m/z array']
        peak_index = np.argmin(abs(np.abs(mz_current_array - target_precursor_mz)))
        summed_intensity += ms1_data[index]['intensity array'][peak_index]
        intensity_dict[0] = float(ms1_data[index]['intensity array'][peak_index])
        mz_dict["mz0"] = ms1_data[index]['m/z array'][peak_index]

        for i in range(self.window):
            mz_next_array = ms1_data[index + i + 1]['m/z array']
            mz_prev_array = ms1_data[index - i - 1]['m/z array']

            if np.min(np.abs(mz_next_array - target_precursor_mz)) <= self.tolerance:
                peak_index = np.argmin(np.abs(mz_next_array - target_precursor_mz))
                summed_intensity += ms1_data[index + i + 1]['intensity array'][peak_index]
                intensity_dict[i + 1] = float(ms1_data[index + i + 1]['intensity array'][peak_index])
                mz_dict["mz" + str(i + 1)] = ms1_data[index + i + 1]['m/z array'][peak_index]
            else:
                intensity_dict[i + 1] = 0
                mz_dict["mz" + str(i + 1)] = "Not Found"

            if np.min(np.abs(mz_prev_array - target_precursor_mz)) <= self.tolerance:
                peak_index = np.argmin(np.abs(mz_prev_array - target_precursor_mz))
                summed_intensity += ms1_data[index - i - 1]['intensity array'][peak_index]
                intensity_dict[-i - 1] = float(ms1_data[index - i - 1]['intensity array'][peak_index])
                mz_dict["mz" + str(-i - 1)] = ms1_data[index - i - 1]['m/z array'][peak_index]
            else:
                intensity_dict[-i - 1] = 0
                mz_dict["mz" + str(-i - 1)] = "Not Found"

        csv_data_sorted = {"scan number" : target_scan_number,
                           "precursor m/z" : target_precursor_mz,
                           "tolerance" : self.tolerance,
                           "total ion intensity" : summed_intensity,
                           **dict(sorted(mz_dict.items())),
                           **dict(sorted(intensity_dict.items()))}
        csv_data_sorted = {str(key) : val for key, val in csv_data_sorted.items()}
        return csv_data_sorted

    def write_file(self):
        with open('ion_intensities.csv', "w", newline="") as file:
            fieldnames = ["scan number", "precursor m/z", "tolerance", "total ion intensity",
                          "mz-10", "mz-9", "mz-8", "mz-7", "mz-6", "mz-5", "mz-4", "mz-3", "mz-2", "mz-1",
                          "mz0", "mz1", "mz2", "mz3", "mz4", "mz5", "mz6", "mz7", "mz8", "mz9", "mz10",
                          "-10", "-9", "-8", "-7", "-6", "-5", "-4", "-3", "-2", "-1",
                          "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10"]
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            for i in range(len(self.scan_numbers)):
                writer.writerow(self.read_mzml(self.scan_numbers[i], self.precursor_mz[i]))

def main():
    intensity_calculator = TotalIonIntensities()
    intensity_calculator.write_file()
    end = timeit.default_timer()
    print(f"INFO: Read spectra from {intensity_calculator.mzml_file}. Output file ion_intensities.csv. Elapsed time {end - intensity_calculator.start}")

if __name__ == '__main__':
    main()