# python PlotIonInjectionTime.py --mzml_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\prm\260113_mEclipse_PRM_ncORF89-AlCl3.mzML" --output_file C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\peptide_089\ion_injection_plots\AlCl_injection_plot.pdf --modifications "Aluminum:669.2930" --scan_range 14000,18000

import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pyteomics import mzml
import numpy as np
import re


def ppm_to_da(mz, ppm):
    return mz * ppm / 1e6


def read_injection_time_sim(mzml_file, target_mz, ppm=4.0, scan_range=None):
    scan_numbers = []
    injection_times_ms = []

    tol_da = ppm_to_da(target_mz, ppm)

    with mzml.read(mzml_file) as reader:
        for spectrum in reader:
            if spectrum.get('ms level', 1) != 1:
                continue

            scan_number = int(spectrum.get('id', '0').split('=')[-1])
            if scan_range and not (scan_range[0] <= scan_number <= scan_range[1]):
                continue

            try:
                scan = spectrum['scanList']['scan'][0]
                filterstring = scan.get('filter string', '')
            except (KeyError, IndexError):
                continue

            match = re.search(r'SIM.*\[(\d+\.?\d*)-(\d+\.?\d*)\]', filterstring)
            if not match:
                continue

            start = float(match.group(1))
            end   = float(match.group(2))

            if not (start <= target_mz <= end):
                continue

            mz_array = spectrum['m/z array']
            intensity_array = spectrum['intensity array']

            mask = (mz_array >= target_mz - tol_da) & (mz_array <= target_mz + tol_da)
            xic_intensity = intensity_array[mask].sum()

            if xic_intensity == 0:
                continue

            injection_time = spectrum.get('ion injection time')
            if injection_time is None:
                continue

            injection_time_ms = float(injection_time)

            scan_numbers.append(scan_number)
            injection_times_ms.append(injection_time_ms)

    if not scan_numbers:
        print(f"No SIM scans found for m/z {target_mz}")

    return scan_numbers, injection_times_ms


def generate_injection_time_pdf(mzml_file, mod_dict, output_file,
                                scan_range=None, ppm=4.0):
    with PdfPages(output_file) as pdf:
        for mod_name, target_mz in mod_dict.items():
            print(f"Processing {mod_name} (m/z={target_mz})...")

            scans, inj_times = read_injection_time_sim(
                mzml_file, target_mz, ppm, scan_range
            )

            if not scans:
                continue

            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(scans, inj_times, label=f"{mod_name} Injection Time")

            if scan_range:
                ax.set_xlim(scan_range)

            ax.set_xlabel("Scan Number")
            ax.set_ylabel("Ion Injection Time (ms)")
            ax.set_title(f"Ion Injection Time for {mod_name} SIM scans")
            ax.grid(True)
            ax.legend()

            pdf.savefig()
            plt.close()

    print(f"PDF saved to: {output_file}")


def parse_range(range_str):
    return tuple(map(int, range_str.split(',')))


def main():
    parser = argparse.ArgumentParser(
        description="Plot ion injection time for Aluminum-containing SIM scans"
    )
    parser.add_argument("--mzml_file", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--modifications", required=True, help="Comma-separated list Name:m/z (e.g. Aluminum:669.2930)")
    parser.add_argument("--scan_range", default=None)
    parser.add_argument("--ppm", type=float, default=4.0)

    args = parser.parse_args()

    mod_dict = {}
    for mod in args.modifications.split(','):
        name, mz = mod.split(':')
        mod_dict[name.strip()] = float(mz.strip())

    scan_range = parse_range(args.scan_range) if args.scan_range else None

    generate_injection_time_pdf(
        args.mzml_file,
        mod_dict,
        args.output_file,
        scan_range,
        args.ppm
    )


if __name__ == "__main__":
    main()
