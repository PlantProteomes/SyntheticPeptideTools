# example usage: py MS1XICExtractor.py --mzml_file [file path] --output_file output.pdf --modifications "TargetPeptide:657.314, Aluminum:669.293" --scan_range "1420,1520" --xic_ppm 4 --delta_ppm 15

import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pyteomics import mzml
import numpy as np
import re
from matplotlib.patches import Rectangle

def ppm_to_da(mz, ppm):
    return mz * ppm / 1e6


def parse_ranges(range_str):
    ranges = []

    for r in range_str.split(";"):
        start, end = map(int, r.split(","))
        ranges.append((start, end))

    return ranges

# detect whether the file is DDA or PRM based on presence of SIM filter strings in MS1 spectra
def detect_acquisition_mode(mzml_file):

    with mzml.read(mzml_file) as reader:

        for spectrum in reader:

            if spectrum.get("ms level", 1) != 1:
                continue

            try:
                scan = spectrum["scanList"]["scan"][0]
                filterstring = scan.get("filter string", "")
            except Exception:
                filterstring = ""

            if "SIM" in filterstring.upper():
                return "PRM"

    return "DDA"

# check if target m/z falls within SIM window defined in filter string
def target_in_sim_window(spectrum, target_mz):

    try:
        scan = spectrum["scanList"]["scan"][0]
        filterstring = scan.get("filter string", "")
    except Exception:
        return False

    if "SIM" not in filterstring.upper():
        return False

    match = re.search(r"\[(\d+\.?\d*)-(\d+\.?\d*)\]", filterstring)

    if not match:
        return False

    start = float(match.group(1))
    end = float(match.group(2))

    return start <= target_mz <= end


def read_xic(
    mzml_file,
    target_mz,
    mode,
    ppm=4.0,
    scan_range=None,
    normalized=True
):

    scan_numbers = []
    intensities = []

    tol_da = ppm_to_da(target_mz, ppm)

    with mzml.read(mzml_file) as reader:

        for spectrum in reader:

            if spectrum.get("ms level", 1) != 1:
                continue

            scan_number = int(spectrum["id"].split("=")[-1])

            if scan_range and not (
                scan_range[0] <= scan_number <= scan_range[1]
            ):
                continue

            if mode == "PRM":

                if not target_in_sim_window(spectrum, target_mz):
                    continue

            mz_array = spectrum["m/z array"]
            intensity_array = spectrum["intensity array"]

            mask = (
                (mz_array >= target_mz - tol_da)
                & (mz_array <= target_mz + tol_da)
            )

            xic_intensity = intensity_array[mask].sum()

            if normalized:

                try:
                    scan = spectrum["scanList"]["scan"][0]
                    inj_time = scan.get("ion injection time", 1) or 1
                except Exception:
                    inj_time = 1

                xic_intensity /= inj_time

            scan_numbers.append(scan_number)
            intensities.append(xic_intensity)

    return scan_numbers, np.array(intensities)


def compute_mass_deltas(
    mzml_file,
    target_mz,
    mode,
    ppm=4.0,
    scan_range=None
):

    scan_numbers = []
    delta_ppms = []
    signals = []

    tol_da = ppm_to_da(target_mz, ppm)

    with mzml.read(mzml_file) as reader:

        for spectrum in reader:

            if spectrum.get("ms level", 1) != 1:
                continue

            scan_number = int(spectrum["id"].split("=")[-1])

            if scan_range and not (
                scan_range[0] <= scan_number <= scan_range[1]
            ):
                continue

            if mode == "PRM":

                if not target_in_sim_window(spectrum, target_mz):
                    continue

            mz_array = spectrum["m/z array"]
            intensity_array = spectrum["intensity array"]

            mask = (
                (mz_array >= target_mz - tol_da)
                & (mz_array <= target_mz + tol_da)
            )

            if not np.any(mask):
                continue

            mz_vals = mz_array[mask]
            int_vals = intensity_array[mask]

            total_intensity = int_vals.sum()

            if total_intensity == 0:
                continue

            observed_mz = np.sum(mz_vals * int_vals) / total_intensity

            delta_ppm = (
                (observed_mz - target_mz)
                / target_mz
                * 1e6
            )

            scan_numbers.append(scan_number)
            delta_ppms.append(delta_ppm)
            signals.append(total_intensity)

    return scan_numbers, delta_ppms, signals


def average_delta_window(scans, deltas, start_scan):
    for i, scan in enumerate(scans):
        if scan >= start_scan:
            window_scans = scans[i:i+5]
            window_deltas = deltas[i:i+5]

            if len(window_deltas) < 5:
                print(f"Warning: only found {len(window_deltas)} points after scan {start_scan}")

            avg = np.mean(window_deltas)

            return window_scans, window_deltas, avg

    return None, None, None


def read_tic(
    mzml_file,
    scan_range=None,
    normalized=True
):

    scan_numbers = []
    intensities = []

    with mzml.read(mzml_file) as reader:

        for spectrum in reader:

            if spectrum.get("ms level", 1) != 1:
                continue

            scan_number = int(spectrum["id"].split("=")[-1])

            if scan_range and not (
                scan_range[0] <= scan_number <= scan_range[1]
            ):
                continue

            intensity_array = spectrum["intensity array"]

            tic_intensity = intensity_array.sum()

            if normalized:

                try:
                    scan = spectrum["scanList"]["scan"][0]
                    inj_time = scan.get("ion injection time", 1) or 1
                except Exception:
                    inj_time = 1

                tic_intensity /= inj_time

            scan_numbers.append(scan_number)
            intensities.append(tic_intensity)

    return scan_numbers, np.array(intensities)


def read_injection_time(
    mzml_file,
    target_mz,
    mode,
    ppm=4.0,
    scan_range=None
):

    scan_numbers = []
    injection_times = []

    tol_da = ppm_to_da(target_mz, ppm)

    with mzml.read(mzml_file) as reader:

        for spectrum in reader:

            if spectrum.get("ms level", 1) != 1:
                continue

            scan_number = int(spectrum["id"].split("=")[-1])

            if scan_range and not (
                scan_range[0] <= scan_number <= scan_range[1]
            ):
                continue

            if mode == "PRM":

                if not target_in_sim_window(spectrum, target_mz):
                    continue

            mz_array = spectrum["m/z array"]
            intensity_array = spectrum["intensity array"]

            mask = (
                (mz_array >= target_mz - tol_da)
                & (mz_array <= target_mz + tol_da)
            )

            if intensity_array[mask].sum() == 0:
                continue

            try:
                scan = spectrum["scanList"]["scan"][0]
                inj_time = scan.get("ion injection time")
            except Exception:
                continue

            if inj_time is None:
                continue

            scan_numbers.append(scan_number)
            injection_times.append(float(inj_time))

    return scan_numbers, injection_times


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument("--mzml_file", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--modifications", required=True)

    parser.add_argument("--scan_range", default=None)

    parser.add_argument("--xic_ppm", type=float, default=4.0)
    parser.add_argument("--delta_ppm", type=float, default=4.0)

    parser.add_argument("--max_y", type=float, default=None)

    parser.add_argument("--no_normalize", action="store_true")

    args = parser.parse_args()

    mods = []

    for mod in args.modifications.split(","):

        name, mz = mod.split(":")
        mods.append((name.strip(), float(mz.strip())))

    scan_ranges = (
        parse_ranges(args.scan_range)
        if args.scan_range
        else [None]
    )

    normalize = not args.no_normalize

    cmap = plt.get_cmap("tab10")

    mode = detect_acquisition_mode(args.mzml_file)

    print(f"\nDetected acquisition mode: {mode}\n")

    with PdfPages(args.output_file) as pdf:

        for scan_range in scan_ranges:

            fig, ax = plt.subplots(figsize=(10, 6))

            xic_data = {}

            for name, mz in mods:

                scans, intensities = read_xic(
                    args.mzml_file,
                    mz,
                    mode,
                    args.xic_ppm,
                    scan_range,
                    normalize
                )

                if len(scans):
                    xic_data[name] = (scans, intensities, mz)

            fig, ax = plt.subplots(figsize=(10, 6))

            xbar_values = []

            for i, (name, mz) in enumerate(mods):

                color = cmap(i % cmap.N)

                scans, deltas, signals = compute_mass_deltas(
                    args.mzml_file,
                    mz,
                    mode,
                    args.delta_ppm,
                    scan_range
                )

                if not scans:
                    continue

                signals = np.array(signals)

                sizes = 18 + (signals / signals.max()) * 100

                ax.scatter(scans, deltas, s=sizes, color=color,
                           label=f"{name} ({mz:.4f})")

                if len(signals) >= 5:

                    best_idx = np.argmax(
                        [np.sum(signals[j:j+5]) for j in range(len(signals)-4)]
                    )

                    w_scans = scans[best_idx:best_idx+5]
                    w_deltas = deltas[best_idx:best_idx+5]

                    xbar = np.mean(w_deltas)
                    xbar_values.append(xbar)

                    # draw rectangle around the 5 points and annotate with xbar
                    rect = Rectangle(
                        (min(w_scans), min(w_deltas)),
                        max(w_scans) - min(w_scans),
                        max(w_deltas) - min(w_deltas),
                        fill=False,
                        edgecolor=color,
                        linewidth=1,
                    )
                    ax.add_patch(rect)

                    ax.annotate(
                        f"x\u0304{i+1} = {xbar:.2f}",
                        xy=(np.mean(w_scans), np.mean(w_deltas)),
                        xytext=(0, 18),
                        textcoords="offset points",
                        ha="center",
                        va="bottom",
                        fontsize=11,
                        color=color,
                        bbox=dict(
                            facecolor="white",
                            edgecolor="none",
                            alpha=0.75
                        ),
                        zorder=10
                    )

            ax.axhline(0, linestyle="--", linewidth=1)
            ax.set_ylim(-6, 7)

            if scan_range:
                ax.set_xlim(scan_range)

            # bottom-right difference (clean + readable)
            if len(xbar_values) >= 2:
                diff = xbar_values[0] - xbar_values[1]

                ax.annotate(
                    f"x\u03041 - x\u03042 = {diff:.2f}",
                    xy=(1, 0),
                    xycoords="axes fraction",
                    xytext=(-10, 10),
                    textcoords="offset points",
                    ha="right",
                    va="bottom",
                    bbox=dict(
                        facecolor="white",
                        edgecolor="none",
                        alpha=0.75
                    ),
                    fontsize=12,
                    zorder=10
                )

            ax.set_xlabel("Scan Number")
            ax.set_ylabel("Mass Delta (ppm)")
            ax.set_title(f"{mode} Mass Delta vs Scan Number")

            ax.legend(
                loc="upper left",
                bbox_to_anchor=(0.01, 0.99),
                frameon=True
            )

            ax.grid(True)

            pdf.savefig(fig)
            plt.close(fig)

    print(f"\nSaved PDF to {args.output_file}\n")


if __name__ == "__main__":
    main()