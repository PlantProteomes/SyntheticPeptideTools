#example py GenerateXICGraphs.py --mzml_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\251103_mEclipse_ncORF89-S3.mzML" --output_file C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\peptide_089\xic_plots\251103_mEclipse_ncORF89-S3_XICPlots_calcium.pdf --modifications "Calcium:676.2875" --scan_range 4000,6000 --tol 0.002
# py GenerateXICGraphs.py --mzml_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\251203_mEclipse_ncORF89-AlK(S04)2.mzML" --output_file C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\peptide_089\xic_plots\251203_mEclipse_ncORF89-AlK(S04)2_XICPlots.pdf --modifications "TargetPeptide:657.3140028,Deamidation:657.8060,Oxidation:665.3115,Calcium+Deamidated:676.7795,Aluminum:669.2930,Nickel+Deamidation:685.7659,Sodium:668.3050,Nickel:685.2738,Iron:683.7697,Zinc:688.2707,Propionamide:692.8326,Diethyl:685.3453,Methyl:664.3218,Dioxidation:673.3089,Formyl:671.3115,Dehydration:657.3140,Phospho:697.2972,Aluminum+Deamidation:669.7850,Calcium:676.2875" --scan_range 4000,6000 --tol 0.002


import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pyteomics import mzml
import numpy as np

def read_xic(mzml_file, target_mz, tol=0.002, scan_range=None):
    """Extract MS1 intensities for a given precursor m/z ± tol"""
    scan_numbers, intensities = [], []
    with mzml.read(mzml_file) as reader:
        for spectrum in reader:
            if spectrum['ms level'] == 1:
                scan_number = int(spectrum['id'].split('=')[-1])
                if scan_range and not (scan_range[0] <= scan_number <= scan_range[1]):
                    continue
                mz_array = spectrum['m/z array']
                intensity_array = spectrum['intensity array']
                mask = (mz_array >= target_mz - tol) & (mz_array <= target_mz + tol)
                xic_intensity = intensity_array[mask].sum()
                scan_numbers.append(scan_number)
                intensities.append(xic_intensity)
    return scan_numbers, intensities

def read_ms2(mzml_file):
    marks = []
    with mzml.read(mzml_file) as reader:
        prev_ms1_scan = None
        for spectrum in reader:
            if spectrum['ms level'] == 1:
                prev_ms1_scan = int(spectrum['id'].split('=')[-1])

            elif spectrum['ms level'] == 2 and 'precursorList' in spectrum:
                precursor_mz = spectrum['precursorList']['precursor'][0]\
                                            ['selectedIonList']['selectedIon'][0]\
                                            ['selected ion m/z']
                marks.append((precursor_mz, prev_ms1_scan))
    return marks

def generate_xic_pdf(mzml_file, mod_dict, output_file, scan_range=None, tol=0.002):
    all_ms2_marks = read_ms2(mzml_file)

    with PdfPages(output_file) as pdf:
        for mod_name, target_mz in mod_dict.items():
            print(f"Processing {mod_name} (m/z={target_mz})...")

            xic_scans, xic_intensity = read_xic(mzml_file, target_mz, tol, scan_range)
            scan_to_intensity = dict(zip(xic_scans, xic_intensity))
            marks = [scan for mz, scan in all_ms2_marks if abs(mz - target_mz) <= tol]
            offset = 1.05
            mark_intensities = [scan_to_intensity.get(s, 0) * offset for s in marks]

            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(xic_scans, xic_intensity, color='purple', label=f'XIC {mod_name}')

            if scan_range:
                ax.set_xlim(scan_range)

            ax.scatter(marks, mark_intensities, color='green', marker='*', s=50)

            ax.set_xlabel('Scan Number')
            ax.set_ylabel('Intensity')
            ax.set_title(f'XIC for {mod_name}. {target_mz} ± {tol} m/z')
            ax.grid(True)

            pdf.savefig()
            plt.close()

    print(f"PDF saved to: {output_file}")

def parse_range(range_str):
    return tuple(map(int, range_str.split(',')))

def main():
    parser = argparse.ArgumentParser(description="Generate XIC plots from mzML with modification names")
    parser.add_argument("--mzml_file", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--modifications", required=True,
                        help="Comma-separated list of modifications and m/z in the form Name:m/z")
    parser.add_argument("--scan_range", default=None,
                        help="Optional scan range to plot, e.g., 3000,6000")
    parser.add_argument("--tol", type=float, default=0.002,
                        help="Tolerance for m/z matching")
    args = parser.parse_args()

    mod_dict = {}
    for mod in args.modifications.split(','):
        name, mz = mod.split(':')
        mod_dict[name.strip()] = float(mz.strip())

    scan_range = parse_range(args.scan_range) if args.scan_range else None

    generate_xic_pdf(args.mzml_file, mod_dict, args.output_file, scan_range, args.tol)

if __name__ == "__main__":
    main()
