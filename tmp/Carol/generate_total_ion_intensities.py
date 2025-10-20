# Program references find_spectra.py made by Nathan Zhang and Eric Deutsch.
# Some reference taken from find_precursor_mz.py made by Mia Chen. Some help from ChatGPT was referenced.



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
    TODO: Input comma-separated list of scan no. and precursor mz.
    """
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog="GenerateTotalIonIntensities",
            description="Iterates through MS1 spectra to find scan intensities with specified precursor mz and sums them.")
        parser.add_argument("--mzml_file", help="Name of mzML file")
        parser.add_argument("--precursor_mz", help="Precursor m/z to select")
        parser.add_argument("--tolerance", help="Tolerance of m/z")
        parser.add_argument("--scan_number", help="Scan number of MS2 scan")

        args = parser.parse_args()

        if args.mzml_file is None or args.mzml_file == "":
            print('ERROR: Parameter --mzml_file must be provided. See --help for more information')
            return
        if not os.path.isfile(args.mzml_file):
            print(f"ERROR: File '{args.mzml_file}' not found or not a file")
            return

        if args.tolerance is None:
            self.tolerance = 0.01
        else:
            self.tolerance = float(args.tolerance)

        self.mzml_file = args.mzml_file
        self.precursor_mz = float(args.precursor_mz)
        self.scan_number = int(args.scan_number)
        self.lower_bound = self.precursor_mz - self.tolerance
        self.upper_bound = self.precursor_mz + self.tolerance
        self.stats = {'counter': 0, 'ms1spectra': 0, 'ms2spectra': 0}

        # self.start = timeit.default_timer()

    def read_mzml(self):
        ms1_data = []
        summed_intensity = 0
        csv_data = [{}]
        index = 0

        with mzml.read(self.mzml_file) as reader:
            for spectrum in reader:
                self.stats['counter'] += 1

                if spectrum['ms level'] == 1:
                    self.stats['ms1spectra'] += 1
                    scan_number = 1 + spectrum['index']
                    ms1_mz_array = spectrum['m/z array']
                    ms1_intensity_array = spectrum['intensity array']
                    ms1_data.append({'scan number' : scan_number, 'm/z array' : ms1_mz_array, 'intensity array' : ms1_intensity_array})

                elif spectrum['ms level'] == 2:
                    self.stats['ms2spectra'] += 1
                    if 1 + spectrum['index'] == self.scan_number:
                        precursor_scan_number = int(re.search(r"scan=(\d+)", spectrum['precursorList']['precursor'][0]['spectrumRef']).group(1))

        for spectrum in ms1_data:
            if spectrum['scan number'] == precursor_scan_number:
                index = ms1_data.index(spectrum)
                print("index is", index)
                break

        mz_current_array = ms1_data[index]['m/z array']
        peak_index = np.argmin(abs(np.abs(mz_current_array - self.precursor_mz)))

        summed_intensity += ms1_data[index]['intensity array'][peak_index]
        csv_data[0][0] = float(ms1_data[index]['intensity array'][peak_index])

        for i in range(10):
            mz_next_array = ms1_data[index + i + 1]['m/z array']
            mz_prev_array = ms1_data[index - i - 1]['m/z array']

            # print("Closest value in next array no", index + i + 1, np.min(abs(np.abs(mz_next_array - self.precursor_mz))))
            # print("Closest value in previous no", index - i - 1, np.min(abs(np.abs(mz_prev_array - self.precursor_mz))))

            if np.min(np.abs(mz_next_array - self.precursor_mz)) <= self.tolerance:
                peak_index = np.argmin(np.abs(mz_next_array - self.precursor_mz))
                # print(peak_index, ms1_data[index + i + 1]['m/z array'][peak_index])
                summed_intensity += ms1_data[index + i + 1]['intensity array'][peak_index]
                # print("Intensity +=", ms1_data[index + i + 1]['intensity array'][peak_index])
                csv_data[0][i + 1] = float(ms1_data[index + i + 1]['intensity array'][peak_index])
            else:
                csv_data[0][i + 1] = "Not Found"

            if np.min(np.abs(mz_prev_array - self.precursor_mz)) <= self.tolerance:
                peak_index = np.argmin(np.abs(mz_prev_array - self.precursor_mz))
                # print(peak_index, ms1_data[index - i - 1]['m/z array'][peak_index])
                summed_intensity += ms1_data[index - i - 1]['intensity array'][peak_index]
                # print("Intensity +=", ms1_data[index - i - 1]['intensity array'][peak_index])
                csv_data[0][-i - 1] = float(ms1_data[index - i - 1]['intensity array'][peak_index])
            else:
                csv_data[0][-i - 1] = "Not Found"


        csv_data_sorted = {"scan number" : self.scan_number, "precursor m/z" : self.precursor_mz, "tolerance" : self.tolerance, "total ion intensity" : summed_intensity, **dict(sorted(csv_data[0].items()))}

        csv_data_sorted = {str(key) : val for key, val in csv_data_sorted.items()}
        print(csv_data_sorted)

        with open('ion_intensities.csv', "w", newline="") as file:
            fieldnames = csv_data_sorted.keys()
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(csv_data_sorted)
        print(summed_intensity)

    def write_file(self):
        return


def main():
    intensity_calculator = TotalIonIntensities()
    intensity_calculator.read_mzml()

if __name__ == '__main__':
    main()