# Content largely adapted from plot_ms2.py and generate_total_ion_intensities.py. ChatGPT was used for brainstorming program structure.

import os
import argparse
import os.path
from pyteomics import mzml, auxiliary
import pandas as pd
import matplotlib.pyplot as plt
import re

class MS1Plot:
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog="PlotMS1",
            description="Plots full and zoomed MS1 spectra of precursors given MS2 scan number: intensity vs. m/z.")
        parser.add_argument("--mzml_file", help="Name of mzML file", required=True)
        parser.add_argument("--xmin", help="Minimum x-axis (scan number) value for zoomed graph. Exponential notation accepted (ex. 1e8)", required=False)
        parser.add_argument("--xmax", help="Maximum x-axis (scan number) value for zoomed graph. Exponential notation accepted (ex. 1e8)", required=False)
        parser.add_argument("--ymin", help="Minimum y-axis (intensity) value. Exponential notation accepted (ex. 1e8)", required=False)
        parser.add_argument("--ymax", help="Maximum y-axis (intensity) value. Exponential notation accepted (ex. 1e8)", required=False)
        parser.add_argument("--output", help="Name of output files. Do not add extension", required=True)
        parser.add_argument("--ms2_scan_number", help="Scan number of MS2 precursor ion if MS1 scan is unknown", type=int)
        parser.add_argument("--ms1_scan_number", help="Scan number of desired MS1 scan if known (no need to input MS2)", type=int)
        args = parser.parse_args()

        if not os.path.isfile(args.mzml_file):
            raise FileNotFoundError

        self.mzml_file = args.mzml_file
        self.xmin = float(args.xmin) if args.xmin else None
        self.xmax = float(args.xmax) if args.xmax else None
        self.ymin = float(args.ymin) if args.ymin else None
        self.ymax = float(args.ymax) if args.ymax else None
        self.output = args.output
        self.ms2_scan_number = args.ms2_scan_number if args.ms2_scan_number else None
        self.ms1_scan_number = args.ms1_scan_number if args.ms1_scan_number else None

    def read_file(self):
        ms1_data = []

        with mzml.read(self.mzml_file) as reader:
            for spectrum in reader:
                precursor_scan_number = 0

                if spectrum['ms level'] == 1:
                    scan_number = 1 + spectrum['index']
                    ms1_mz_array = spectrum['m/z array']
                    ms1_intensity_array = spectrum['intensity array']
                    ms1_data.append(
                        {'scan number': scan_number, 'm/z array': ms1_mz_array, 'intensity array': ms1_intensity_array})

                if spectrum['ms level'] == 2:
                    if 1 + spectrum['index'] == self.ms2_scan_number:
                        print(f"scan number {self.ms2_scan_number} found")
                        precursor_scan_number = int(
                            re.search(r"scan=(\d+)", spectrum['precursorList']['precursor'][0]['spectrumRef']).group(1))
                        print(precursor_scan_number)
                        break

            for ms1_spectrum in ms1_data:
                if ms1_spectrum['scan number'] == precursor_scan_number:
                    print("spectrum found")
                    print(ms1_spectrum['m/z array'])
                    ms1_mz_array = ms1_spectrum['m/z array']
                    ms1_intensity_array = ms1_spectrum['intensity array']

        return ms1_mz_array, ms1_intensity_array, precursor_scan_number

    def get_ms1(self):
        with mzml.read(self.mzml_file) as reader:
            for spectrum in reader:
                if spectrum['ms level'] == 1 and 1 + spectrum['index'] == self.ms1_scan_number:
                    scan_number = 1 + spectrum['index']
                    ms1_mz_array = spectrum['m/z array']
                    ms1_intensity_array = spectrum['intensity array']
        return ms1_mz_array, ms1_intensity_array, scan_number

    def plot_ms1(self):
        ms1_mz_array, ms1_intensity_array, precursor_scan_number = self.read_file() if self.ms2_scan_number else self.get_ms1()
        df = pd.DataFrame({"mz" : ms1_mz_array, "intensity" : ms1_intensity_array})
        markers, stems, base = plt.stem(df["mz"], df["intensity"], markerfmt=" ")

        plt.xlabel("m/z")
        plt.ylabel("Intensity")
        plt.suptitle("Intensity vs. m/z of MS1 scan number " + str(precursor_scan_number))
        plt.savefig(self.output + '_plot.pdf')

        ax = plt.gca()
        ax.set_xlim(self.xmin, self.xmax)
        ax.set_ylim(self.ymin, self.ymax)
        plt.suptitle("Intensity vs. m/z of MS1 scan number " + str(precursor_scan_number) + " (Zoomed)")
        plt.savefig(self.output + '_plot_zoomed.pdf')

def main():
    ms1_plot = MS1Plot()
    print(f"INFO: Plotting to {ms1_plot.output}_plot.pdf and {ms1_plot.output}_plot_zoomed.pdf")
    ms1_plot.plot_ms1()

if __name__ == "__main__":
    main()