# example:
# py PRM_SIM_MassDeltaPlots.py --mzml_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\prm\260113_mEclipse_PRM_ncORF89-AlCl3.mzML" --output_file C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\peptide_089\sim_plots\AlCldeltas.pdf --modifications "Aluminum:669.2930" --scan_range 15000,16000 --ppm 4

import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pyteomics import mzml
import numpy as np
import re


def ppm_to_da(mz, ppm):
    return mz * ppm / 1e6


def compute_mass_deltas(mzml_file, target_mz, ppm=4.0, scan_range=None):
    scan_numbers = []
    delta_ppms = []

    tol_da = ppm_to_da(target_mz, ppm)

    with mzml.read(mzml_file) as reader:
        for spectrum in reader:

            # MS1 only (SIM)
            if spectrum.get("ms level", 1) != 1:
                continue

            scan_number = int(spectrum.get("id", "0").split("=")[-1])

            if scan_range and not (scan_range[0] <= scan_number <= scan_range[1]):
                continue

            try:
                scan = spectrum["scanList"]["scan"][0]
                filterstring = scan.get("filter string", "")
            except (KeyError, IndexError):
                continue

            # Require SIM window
            match = re.search(r"\[(\d+\.?\d*)-(\d+\.?\d*)\]", filterstring)
            if not match:
                continue

            start = float(match.group(1))
            end   = float(match.group(2))

            # Only SIM windows containing target
            if not (start <= target_mz <= end):
                continue

            mz_array = spectrum["m/z array"]
            intensity_array = spectrum["intensity array"]

            mask = (mz_array >= target_mz - tol_da) & (mz_array <= target_mz + tol_da)

            if not np.any(mask):
                continue

            mz_vals = mz_array[mask]
            int_vals = intensity_array[mask]

            if int_vals.sum() == 0:
                continue

            # Intensity-weighted centroid m/z
            observed_mz = np.sum(mz_vals * int_vals) / np.sum(int_vals)

            delta_ppm = (observed_mz - target_mz) / target_mz * 1e6

            scan_numbers.append(scan_number)
            delta_ppms.append(delta_ppm)

    return scan_numbers, delta_ppms


def parse_ranges(range_str):
    ranges = []
    for r in range_str.split(";"):
        start, end = map(int, r.split(","))
        ranges.append((start, end))
    return ranges


def generate_delta_pdf(mzml_file, mod_dict, pdf, scan_range=None, ppm=4.0):
    for mod_name, target_mz in mod_dict.items():
        print(f"Processing {mod_name} (m/z={target_mz})...")

        scans, deltas = compute_mass_deltas(
            mzml_file,
            target_mz,
            ppm=ppm,
            scan_range=scan_range
        )

        if not scans:
            print(f"No data for {mod_name}")
            continue

        fig, ax = plt.subplots(figsize=(10, 6))

        ax.scatter(scans, deltas, s=18)
        ax.axhline(0, linestyle="--", linewidth=1)

        if scan_range:
            ax.set_xlim(scan_range)

        ax.set_xlabel("Scan Number")
        ax.set_ylabel("Mass Delta (ppm)")
        ax.set_title(f"Mass Error for {mod_name} ({target_mz:.4f}) ± {ppm} ppm")
        ax.grid(True)

        pdf.savefig(fig)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot mass deltas (ppm) vs scan number for SIM scans")
    parser.add_argument("--mzml_file", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--modifications", required=True, help="Name:m/z,Name:m/z")
    parser.add_argument("--scan_range", default=None)
    parser.add_argument("--ppm", type=float, default=4.0)

    args = parser.parse_args()

    mod_dict = {}
    for mod in args.modifications.split(","):
        name, mz = mod.split(":")
        mod_dict[name.strip()] = float(mz.strip())

    if args.scan_range:
        scan_ranges = parse_ranges(args.scan_range)
    else:
        scan_ranges = [None]

    with PdfPages(args.output_file) as pdf:
        for scan_range in scan_ranges:
            print(f"\nScan range: {scan_range}")
            generate_delta_pdf(
                args.mzml_file,
                mod_dict,
                pdf,
                scan_range=scan_range,
                ppm=args.ppm
            )

    print(f"Saved mass-delta plots to {args.output_file}")


if __name__ == "__main__":
    main()
