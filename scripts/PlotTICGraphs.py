# py generate_tic_plots.py --mzml_file --annotation_file --output_file --ms1_zoom --overlay_zoom --ms2_ranges --precursor_mz --tolerance
# py TICplots.py --mzml_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\250402_mEclipse_QC_ncORF-055.mzML" --annotation_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\new_peptide\scan_identifications_055.xlsx" --tolerance 0.003 --precursor_mz 562.3260 --output_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\new_peptide\TIC_plots_055.pdf" --ms1_zoom 2000,2900 --overlay_zoom 2000,2900 --ms2_ranges 2000,2200;2200,2250;2250,2300 



import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from pyteomics import mzml
import numpy as np


def read_ms1_tic(mzml_file):
    ms1_scan_numbers, ms1_tic_values = [], []
    with mzml.read(mzml_file) as reader:
        for spectrum in reader:
            if spectrum['ms level'] == 1:
                scan_number = int(spectrum['id'].split('=')[-1])
                tic = sum(spectrum['intensity array'])
                ms1_scan_numbers.append(scan_number)
                ms1_tic_values.append(tic)
    return ms1_scan_numbers, ms1_tic_values


def read_ms2_tic(mzml_file):
    ms2_scan_numbers, ms2_tic_values = [], []
    with mzml.read(mzml_file) as reader:
        for spectrum in reader:
            if spectrum['ms level'] == 2:
                scan_number = int(spectrum['id'].split('=')[-1])
                tic = sum(spectrum['intensity array'])
                ms2_scan_numbers.append(scan_number)
                ms2_tic_values.append(tic)
    return ms2_scan_numbers, ms2_tic_values


def read_ms2_precursors(mzml_file, target_mz=None, tol=0.003):
    """Return a list of MS1 scan numbers to mark based on MS2 precursor m/z matching target_mz"""
    ms1_mark_scans = []
    with mzml.read(mzml_file) as reader:
        prev_ms1_scan = None
        for spectrum in reader:
            if spectrum['ms level'] == 1:
                prev_ms1_scan = int(spectrum['id'].split('=')[-1])
            elif spectrum['ms level'] == 2 and 'precursorList' in spectrum:
                precursor_mz = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z']
                if target_mz is not None and abs(precursor_mz - target_mz) <= tol:
                    if prev_ms1_scan is not None:
                        ms1_mark_scans.append(prev_ms1_scan)
    return ms1_mark_scans


def read_annotations(annotation_file):
    df = pd.read_excel(annotation_file)
    if not {'ScanNumber', 'Identification'}.issubset(df.columns):
        raise ValueError("Annotation file must contain columns 'ScanNumber' and 'Identification'")
    return dict(zip(df['ScanNumber'], df['Identification']))


def parse_range(range_str, multiple=False):
    if multiple:
        return [tuple(map(int, r.split(','))) for r in range_str.split(';')]
    else:
        return tuple(map(int, range_str.split(',')))


def generate_plots(ms1_scans, ms1_tics, ms2_scans, ms2_tics, annotations, output_file,
                   ms1_zoom=(2000, 4000), overlay_zoom=(2000, 4000), ms2_ranges=None, scans_to_mark=None):
    if scans_to_mark is None:
        scans_to_mark = list(annotations.keys())

    scaled_ms2_tics = [v * 80 for v in ms2_tics]

    def get_tic(scan_list, tic_list, scan):
        idx = next((i for i, s in enumerate(scan_list) if s > scan), len(tic_list) - 1)
        return tic_list[idx]

    y_list = [get_tic(ms1_scans, ms1_tics, s) for s in scans_to_mark]

    if ms2_ranges is None:
        ms2_ranges = [(2350, 2750), (2750, 2920), (3200, 3730)]

    with PdfPages(output_file) as pdf:
        # --- GRAPH 1: Full MS1 TIC ---
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(ms1_scans, ms1_tics, color='#ff00c1', label='MS1 TIC')
        ax.scatter(scans_to_mark, [y + 2e7 for y in y_list], color='green', marker='*', s=50)
        ax.set_xlabel('Scan Number')
        ax.set_ylabel('Total Ion Current (TIC)')
        ax.set_title('MS1 TIC vs Scan Number')
        ax.grid(True)
        pdf.savefig()
        plt.close()

        # --- GRAPH 2: Zoomed MS1 TIC ---
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(ms1_scans, ms1_tics, color='#ff00c1')
        ax.scatter(scans_to_mark, [y + 2e7 for y in y_list], color='green', marker='*', s=50)
        ax.set_xlim(ms1_zoom)
        ax.set_ylim(0, max(ms1_tics) * 1.1)
        ax.set_xlabel('Scan Number')
        ax.set_ylabel('Total Ion Current (TIC)')
        ax.set_title(f"MS1 TIC (Zoomed {ms1_zoom[0]}–{ms1_zoom[1]})")
        ax.grid(True)
        pdf.savefig()
        plt.close()

        # --- GRAPH 3: MS1 + MS2 Overlay ---
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(ms1_scans, ms1_tics, color='#ff00c1', label='MS1')
        ax.plot(ms2_scans, scaled_ms2_tics, color='blue', label='MS2 (x80)')
        ax.scatter(scans_to_mark, [y + 2e7 for y in y_list], color='green', marker='*', s=50)
        ax.set_xlim(overlay_zoom)
        ax.set_ylim(0, max(ms1_tics) * 1.1)
        ax.set_xlabel('Scan Number')
        ax.set_ylabel('Total Ion Current (TIC)')
        ax.set_title(f"MS1 and MS2 Overlay ({overlay_zoom[0]}–{overlay_zoom[1]})")
        ax.legend()
        ax.grid(True)
        pdf.savefig()
        plt.close()

        # --- GRAPH 4–n: Zoomed MS2 with annotations ---
        for x_min, x_max in ms2_ranges:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(ms2_scans, ms2_tics, color='blue')
            used_y = []

            for scan, label in annotations.items():
                if x_min <= scan <= x_max and scan in ms2_scans:
                    idx = ms2_scans.index(scan)
                    tic_val = ms2_tics[idx]
                    label_y = tic_val * 1.05

                    for prev_y in used_y:
                        if abs(prev_y - label_y) < max(ms2_tics) * 0.03:
                            label_y += max(ms2_tics) * 0.03
                    used_y.append(label_y)

                    ax.plot([scan, scan], [tic_val, label_y], color='black', linestyle='--', linewidth=0.7)
                    ax.text(scan, label_y, f"{scan}", color='black', fontsize=7, rotation=20, ha='center', va='bottom')

            legend_text = [f"{scan}: {label}" for scan, label in sorted(annotations.items()) if x_min <= scan <= x_max]
            legend_labels = "\n".join(legend_text)
            plt.gcf().text(1.02, 0.5, legend_labels, va='center', fontsize=8, transform=plt.gca().transAxes)
            plt.subplots_adjust(right=0.75)

            ax.set_xlabel('Scan Number')
            ax.set_ylabel('Total Ion Current (TIC)')
            ax.set_title(f'MS2 TIC {x_min}-{x_max}')
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(0, max(used_y) * 1.1 if used_y else max(ms2_tics))
            ax.grid(True)
            pdf.savefig()
            plt.close()


def main():
    parser = argparse.ArgumentParser(description="Generate MS1/MS2 TIC-over-ScanNumber plots with annotations.")
    parser.add_argument("--mzml_file", required=True, help="Path to mzML file")
    parser.add_argument("--annotation_file", required=True, help="Excel file with columns ScanNumber and Identification")
    parser.add_argument("--output_file", required=True, help="Output PDF file path")
    parser.add_argument("--ms1_zoom", default="2000,4000", help="Zoom range for MS1 TIC plot (e.g. 2000,4000)")
    parser.add_argument("--overlay_zoom", default="2000,4000", help="Zoom range for overlay plot (e.g. 2000,4000)")
    parser.add_argument("--ms2_ranges", default="2350,2750;2750,2920;3200,3730",
                        help="Semicolon-separated zoom ranges for MS2 (e.g. 2350,2750;2750,2920;3200,3730)")
    parser.add_argument("--precursor_mz", type=float, help="Peptide precursor m/z to mark on MS1 TIC")
    parser.add_argument("--tolerance", type=float, default=0.003, help="Tolerance for precursor m/z matching")
    
    args = parser.parse_args()

    ms1_zoom = parse_range(args.ms1_zoom)
    overlay_zoom = parse_range(args.overlay_zoom)
    ms2_ranges = parse_range(args.ms2_ranges, multiple=True)

    print("Reading MS1 data...")
    ms1_scans, ms1_tics = read_ms1_tic(args.mzml_file)
    print(f"  → {len(ms1_scans)} MS1 scans read")

    print("Reading MS2 data...")
    ms2_scans, ms2_tics = read_ms2_tic(args.mzml_file)
    print(f"  → {len(ms2_scans)} MS2 scans read")

    print("Reading annotations...")
    annotations = read_annotations(args.annotation_file)
    print(f"  → {len(annotations)} annotations loaded")

    # determine which MS1 scans to mark with stars
    if args.precursor_mz is not None:
        scans_to_mark = read_ms2_precursors(args.mzml_file, args.precursor_mz, args.tolerance)
    else:
        scans_to_mark = list(annotations.keys())

    print("Generating plots...")
    generate_plots(ms1_scans, ms1_tics, ms2_scans, ms2_tics, annotations,
                   args.output_file, ms1_zoom, overlay_zoom, ms2_ranges, scans_to_mark=scans_to_mark)
    print(f"PDF successfully written to: {args.output_file}")


if __name__ == "__main__":
    main()
