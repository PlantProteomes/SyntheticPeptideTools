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

# read XIC for given target m/z and mode, applying ppm tolerance and optional scan range filter
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

            # PRM/SIM restriction
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

# compute mass deltas in ppm for observed m/z vs target m/z across scans, with point size scaled by signal intensity
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

            # PRM restriction
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

# compute average delta ppm in a window of 5 scans starting from the first scan >= specified start scan, to reveal any trends in mass accuracy that may correlate with presence of target peptide or other modifications
# take the five most intense points in the window to compute the average, to focus on scans where target peptide is likely present and reduce noise from low-intensity scans where mass accuracy may be less reliable
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

# read TIC across scans, applying optional scan range filter and normalization by injection time
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

# read injection time for MS1 spectra where target m/z is detected, applying ppm tolerance and optional scan range filter
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

# main function to parse arguments, detect acquisition mode, and generate PDF with XIC, mass delta, TIC, and injection time plots for specified modifications and scan ranges

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

            # page 1: XIC overlay with new scaling logic to make TargetPeptide more visible when present, while preserving relative intensities of other modifications

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

                    xic_data[name] = (
                        scans,
                        intensities,
                        mz
                    )
            
            scale_factor = 1.0

            if "TargetPeptide" in xic_data:

                target_max = xic_data["TargetPeptide"][1].max()

                other_maxes = []

                for key in xic_data:

                    if key == "TargetPeptide":
                        continue

                    other_maxes.append(
                        xic_data[key][1].max()
                    )

                if len(other_maxes):

                    second_tallest = max(other_maxes)

                    desired_target_height = (
                        second_tallest * 1.5
                    )

                    if desired_target_height > 0:

                        scale_factor = (
                            target_max
                            / desired_target_height
                        )

            for i, (name, _) in enumerate(mods):

                if name not in xic_data:
                    continue

                scans, intensities, mz = xic_data[name]

                color = cmap(i % cmap.N)

                if (
                    name == "TargetPeptide"
                    and scale_factor != 1
                ):

                    intensities = intensities / 1e5

                    label = (
                        f"{name} "
                        f"(/1e5)"
                        f"(m/z {mz:.4f})"
                    )

                else:

                    label = f"{name} (m/z {mz:.4f})"

                ax.plot(
                    scans,
                    intensities,
                    linewidth=2,
                    color=color,
                    label=label
                )

            if scan_range:
                ax.set_xlim(scan_range)

            if args.max_y:
                ax.set_ylim(0, args.max_y)

            ax.set_xlabel("Scan Number")

            ax.set_ylabel(
                "Normalized XIC (Intensity / ms)"
                if normalize
                else "XIC Intensity"
            )

            ax.set_title(
                f"{mode} XIC Overlay "
                f"(± {args.xic_ppm} ppm)"
            )

            ax.legend(
                loc="upper left",
                bbox_to_anchor=(0.01, 0.99),
                frameon=True
            )

            ax.grid(True)

            pdf.savefig(fig)
            plt.close(fig)

            # page 2: mass delta vs scan number, with point size scaled by signal intensity, and horizontal line at 0 ppm to indicate perfect match

            fig, ax = plt.subplots(figsize=(10, 6))

            xbar_values = []

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

                        # rectangle
                        rect = Rectangle(
                            (min(w_scans), min(w_deltas)),
                            max(w_scans) - min(w_scans),
                            max(w_deltas) - min(w_deltas),
                            fill=False,
                            edgecolor=color,
                            linewidth=1
                        )
                        ax.add_patch(rect)

                        ax.annotate(
                            f"x\u0304{i+1} = {xbar:.3f}",
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

                if len(xbar_values) >= 2:
                    diff = xbar_values[0] - xbar_values[1]

                    ax.annotate(
                        f"x\u03041 - x\u03042 = {diff:.3f}",
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


            # page 3: TIC vs scan number, with annotation of scan with highest TIC, and optional normalization by injection time to reveal trends that may be obscured by varying injection times

            fig, ax = plt.subplots(figsize=(10, 6))

            scans, tic_intensity = read_tic(
                args.mzml_file,
                scan_range,
                normalize
            )

            if len(scans):

                ax.plot(
                    scans,
                    tic_intensity,
                    color = cmap(i % cmap.N),
                    linewidth=2,
                    label="Total Ion Chromatogram"
                )

                peak_idx = np.argmax(tic_intensity)

                peak_scan = scans[peak_idx]

                peak_intensity = tic_intensity[peak_idx]

                ax.annotate(
                    f"Scan {peak_scan}",
                    xy=(peak_scan, peak_intensity),
                    xytext=(10, 10),
                    textcoords="offset points",
                    arrowprops=dict(
                        arrowstyle="->"
                    ),
                    fontsize=10
                )

            if scan_range:
                ax.set_xlim(scan_range)

            if args.max_y:
                ax.set_ylim(0, args.max_y)

            ax.set_xlabel("Scan Number")

            ax.set_ylabel(
                "Normalized TIC (Intensity / ms)"
                if normalize
                else "TIC Intensity"
            )

            ax.set_title(
                f"{mode} Total Ion Chromatogram"
            )

            ax.legend(
                loc="upper left",
                bbox_to_anchor=(0.01, 0.99),
                frameon=True
            )

            ax.grid(True)

            pdf.savefig(fig)
            plt.close(fig)

            # page 4: injection time vs scan number for MS1 spectra where target m/z is detected, to reveal any trends in injection time that may correlate with presence of target peptide or other modifications

            fig, ax = plt.subplots(figsize=(10, 6))

            for i, (name, mz) in enumerate(mods):

                color = cmap(i % cmap.N)

                scans, inj = read_injection_time(
                    args.mzml_file,
                    mz,
                    mode,
                    args.xic_ppm,
                    scan_range
                )

                if not scans:
                    continue

                ax.plot(
                    scans,
                    inj,
                    marker="*",
                    linewidth=2,
                    color=color,
                    label=f"{name} ({mz:.4f})"
                )

            if scan_range:
                ax.set_xlim(scan_range)

            ax.set_xlabel("Scan Number")

            ax.set_ylabel(
                "Ion Injection Time (ms)"
            )

            ax.set_title(
                f"{mode} Injection Time"
            )

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
