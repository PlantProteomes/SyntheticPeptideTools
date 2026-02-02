# py Overlay_PRM_XICGraphs.py --mzml_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\prm\260113_mEclipse_PRM_ncORF89.mzML" --output_file C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\peptide_089\SIM_plots\no_aluminum_SIMoverlay.pdf --modifications Oxidation:673.305995,TargetPeptide:657.3140 --ppm 10 --no_normalize --scan_range 15000,16000


import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pyteomics import mzml
import re


def ppm_to_da(mz, ppm):
    return mz * ppm / 1e6


def read_xic_sim(mzml_file, target_mz, ppm=4.0, scan_range=None, normalized=True):
    scan_numbers = []
    intensities = []
    tol_da = ppm_to_da(target_mz, ppm)

    with mzml.read(mzml_file) as reader:
        for spectrum in reader:

            if spectrum.get("ms level", 1) != 1:
                continue

            scan_number = int(spectrum["id"].split("=")[-1])
            if scan_range and not (scan_range[0] <= scan_number <= scan_range[1]):
                continue

            try:
                scan = spectrum["scanList"]["scan"][0]
                filterstring = scan.get("filter string", "")
            except (KeyError, IndexError):
                continue

            if "SIM" not in filterstring:
                continue

            match = re.search(r"\[(\d+\.?\d*)-(\d+\.?\d*)\]", filterstring)
            if not match:
                continue

            start = float(match.group(1))
            end = float(match.group(2))

            if not (start <= target_mz <= end):
                continue

            mz_array = spectrum["m/z array"]
            intensity_array = spectrum["intensity array"]

            mask = (mz_array >= target_mz - tol_da) & (mz_array <= target_mz + tol_da)
            xic_intensity = intensity_array[mask].sum()

            if xic_intensity == 0:
                scan_numbers.append(scan_number)
                intensities.append(0)
                continue

            if normalized:
                inj_time = scan.get("ion injection time", 1) or 1
                xic_intensity = (xic_intensity / inj_time) * 100

            scan_numbers.append(scan_number)
            intensities.append(xic_intensity)

    return scan_numbers, intensities


def parse_ranges(range_str):
    ranges = []
    for r in range_str.split(";"):
        start, end = map(int, r.split(","))
        ranges.append((start, end))
    return ranges


def main():
    parser = argparse.ArgumentParser(description="Overlay SIM XICs for multiple targets")
    parser.add_argument("--mzml_file", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument(
        "--modifications",
        required=True,
        help="Name:m/z,Name:m/z,Name:m/z..."
    )
    parser.add_argument("--scan_range", default=None)
    parser.add_argument("--ppm", type=float, default=4.0)
    parser.add_argument("--max_y", type=float, default=None)
    parser.add_argument("--no_normalize", action="store_true")

    args = parser.parse_args()

    mods = []
    for mod in args.modifications.split(","):
        name, mz = mod.split(":")
        mods.append((name.strip(), float(mz.strip())))

    scan_ranges = parse_ranges(args.scan_range) if args.scan_range else [None]
    normalize = not args.no_normalize

    cmap = plt.get_cmap("tab10")

    with PdfPages(args.output_file) as pdf:
        for scan_range in scan_ranges:
            fig, ax = plt.subplots(figsize=(10, 6))

            for i, (name, mz) in enumerate(mods):
                color = cmap(i % cmap.N)

                scans, intensities = read_xic_sim(
                    args.mzml_file,
                    mz,
                    ppm=args.ppm,
                    scan_range=scan_range,
                    normalized=normalize,
                )

                if name == "TargetPeptide":
                    intensities = [x / 500 for x in intensities]

                if not scans:
                    continue

                label = f"{name} ({mz:.4f})"
                if name == "TargetPeptide":
                    label += " /500"

                ax.plot(
                    scans,
                    intensities,
                    label=label,
                    color=color,
                    linewidth=2,
                )

                # Mark max
                max_idx = intensities.index(max(intensities))
                max_scan = scans[max_idx]
                max_intensity = intensities[max_idx]

                ax.scatter(
                    max_scan,
                    max_intensity,
                    marker="*",
                    s=200,
                    color=color,
                    zorder=6,
                )

                ax.annotate(
                    f"Scan {max_scan} Intensity {max_intensity:.2f}",
                    xy=(max_scan, max_intensity),
                    xytext=(5, 5),
                    textcoords="offset points",
                    fontsize=9,
                    color=color,
                    weight="bold",
                )

            if scan_range:
                ax.set_xlim(scan_range)
            if args.max_y is not None:
                ax.set_ylim(0, args.max_y)

            ax.set_xlabel("Scan Number")
            ax.set_ylabel(
                "Normalized Intensity (XIC / ms × 100)"
                if normalize
                else "Intensity (XIC)"
            )
            ax.set_title(f"SIM XIC Overlay (±{args.ppm} ppm)")
            ax.legend()
            ax.grid(True)

            pdf.savefig(fig)
            plt.close(fig)

    print(f"Saved overlay SIM XICs to {args.output_file}")


if __name__ == "__main__":
    main()