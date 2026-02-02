# example python PRM_XICPlots.py --mzml_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\prm\260113_mEclipse_PRM_ncORF89-AlCl3.mzML" --output_file C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\peptide_089\sim_plots\AlCl_SIM_Plots.pdf --modifications "Aluminum:669.2930" --scan_range 14000,18000
# 

import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pyteomics import mzml
import numpy as np
import re 

def ppm_to_da(mz, ppm):
    return mz * ppm / 1e6

def read_xic_sim(mzml_file, target_mz, ppm=4.0, scan_range=None, normalized=True):
    scan_numbers, norm_intensities = [], []
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
                scan_numbers.append(scan_number)
                norm_intensities.append(0)
                continue

            injection_time = scan.get('ion injection time', 1)
            if injection_time == 0:
                injection_time = 1

            if normalized:
                injection_time = scan.get('ion injection time', 1)
                if injection_time == 0:
                    injection_time = 1
                xic_intensity = (xic_intensity / injection_time) * 100

            scan_numbers.append(scan_number)
            norm_intensities.append(xic_intensity)

    if not scan_numbers:
        print(f"No SIM scans found for target m/z {target_mz} in the specified scan range.")
    return scan_numbers, norm_intensities

def generate_xic_pdf(mzml_file, mod_dict, pdf, scan_range=None, ppm=4.0, max_y=None, normalized=True):
    for mod_name, target_mz in mod_dict.items():
        print(f"Processing {mod_name} (m/z={target_mz}) for scan range {scan_range}...")
        xic_scans, xic_intensity = read_xic_sim(mzml_file, target_mz, ppm, scan_range, normalized)

        if not xic_scans:
            continue

        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(xic_scans, xic_intensity, color='purple', label=f'XIC {mod_name}')
                # , marker='*')

        if scan_range:
            ax.set_xlim(scan_range)
        if max_y is not None:
            ax.set_ylim(0, max_y)

        ax.set_xlabel('Scan Number')
        ax.set_ylabel('Normalized Intensity (XIC / ms * 100)' if normalized else 'Intensity (XIC)')
        ax.set_title(f'SIM XIC for {mod_name}: {target_mz} ± {ppm} ppm')
        ax.grid(True)
        ax.legend()
        pdf.savefig()
        plt.close()

    print(f"saved")

def parse_ranges(range_str):
    ranges = []
    for r in range_str.split(';'):
        start, end = map(int, r.split(','))
        ranges.append((start, end))
    return ranges


def main():
    parser = argparse.ArgumentParser(description="Generate XIC plots from mzML using ppm tolerance")
    parser.add_argument("--mzml_file", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--modifications", required=True, help="Comma-separated list of modifications Name:m/z")
    parser.add_argument("--scan_range", default=None, help="Optional scan range, e.g. 3000,6000")
    parser.add_argument("--ppm", type=float, default=4.0, help="Mass tolerance in ppm")
    parser.add_argument("--max_y", type=float, default=None, help="Optional maximum Y-axis value (intensity)")
    parser.add_argument("--no_normalize", action='store_true', help="Disable normalization of XIC intensities")

    args = parser.parse_args()

    mod_dict = {}
    for mod in args.modifications.split(','):
        name, mz = mod.split(':')
        mod_dict[name.strip()] = float(mz.strip())

    if args.scan_range:
        scan_ranges = parse_ranges(args.scan_range)
    else:
        scan_ranges = [None]

    normalize = not args.no_normalize

    with PdfPages(args.output_file) as pdf:
        for scan_range in scan_ranges:
            print(f"\nProcessing scan range: {scan_range}")
            generate_xic_pdf(args.mzml_file, mod_dict, pdf, scan_range, args.ppm, args.max_y, normalize)

if __name__ == "__main__":
    main()