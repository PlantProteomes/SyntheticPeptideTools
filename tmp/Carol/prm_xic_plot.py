import re
import matplotlib.pyplot as plt
from pyteomics import mzml
import os
import argparse
import os.path
from matplotlib.backends.backend_pdf import PdfPages

class PRMPlotXIC:
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog="PlotXIC",
            description="Plots XIC graph with specified ppm tolerance.")
        parser.add_argument("--mzml_file", required=True)
        parser.add_argument("--ppm_tolerance", type=int, required=True)
        parser.add_argument("--precursor_mz", type=float, required=True)
        parser.add_argument("--output", help="Name of output file")
        parser.add_argument("--xmin", type=float)
        parser.add_argument("--xmax", type=float)
        args = parser.parse_args()

        if not os.path.isfile(args.mzml_file):
            raise FileNotFoundError

        self.mzml_file = args.mzml_file
        self.ppm_tolerance = args.ppm_tolerance
        self.precursor_mz = args.precursor_mz
        self.xmin = args.xmin
        self.xmax = args.xmax
        if args.output.endswith(".pdf"):
            self.output = args.output
        else:
            self.output = args.output + ".pdf"

    def plot_prm_xic(self):
        with mzml.read(self.mzml_file) as reader:
            retention_times = []
            intensities = []
            injection_times = []

            for spectrum in reader:
                filter_string = spectrum['scanList']['scan'][0]['filter string']
                match = re.search(r'NSI (\S+) (\S+).* \[([\d\.]+)\-([\d\.]+)\]', filter_string)
                if match:
                    scan_type = match.group(1)
                    ms_level = match.group(2)
                    start = float(match.group(3))
                    end = float(match.group(4))

                    if scan_type == "SIM" and ms_level == "ms" and start <= self.precursor_mz <= end:
                        ms1_mz_array = spectrum['m/z array']
                        ms1_intensity_array = spectrum['intensity array']
                        intensity_sum = 0.0
                        for i in range(len(ms1_mz_array)):
                            mz = ms1_mz_array[i]
                            if abs(mz - self.precursor_mz)/self.precursor_mz * 1000000 <= self.ppm_tolerance:
                                intensity_sum += ms1_intensity_array[i]
                        injection_time = spectrum['scanList']['scan'][0]['ion injection time']
                        injection_times.append(injection_time)
                        # normalized_intensity = intensity_sum
                        normalized_intensity = intensity_sum / injection_time * 100
                        intensities.append(normalized_intensity)
                        retention_times.append(spectrum['scanList']['scan'][0]['scan start time'])
                else:
                    print(f"ERROR: unable to interpret {filter_string}")
                    break

        with PdfPages(self.output) as pdf:
            plt.plot(retention_times, intensities)
            plt.xlabel("Retention time (s)")
            plt.ylabel("Intensity")
            plt.title("Extracted Ion Chromatogram")
            pdf.savefig()

            plt.xlim(self.xmin, self.xmax)
            pdf.savefig()
            plt.close()

def main():
    prm_plot = PRMPlotXIC()
    prm_plot.plot_prm_xic()

if __name__ == "__main__":
    main()