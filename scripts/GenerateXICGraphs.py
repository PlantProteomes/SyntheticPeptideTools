# py GenerateXICGraphs.py --mzml_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\251112_mEclipse_ncORF89-S4_MC.mzML" --output_file C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\peptide_089\xic_plots_4ppm\251103_mEclipse_ncORF89-MC_XICPlots.pdf --modifications "TargetPeptide:657.3140028,Aluminum:669.2930,Sodium:668.3050,Iron:683.7697,Aluminum+Deamidation:669.7850,Calcium:676.2875,sodium+deamidation:668.7970" --scan_range 4000,6000 
# py GenerateXICGraphs.py --mzml_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\NEW_250402_mEclipse_QC_ncORF-089.mzML" --output_file C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\peptide_089\xic_plots\251203_mEclipse_ncORF89-AlK(S04)2_XICPlots.pdf --modifications "TargetPeptide:657.3140028,Deamidation:657.8060,Oxidation:665.3115,Calcium+Deamidated:676.7795,Aluminum:669.2930,Nickel+Deamidation:685.7659,Sodium:668.3050,Nickel:685.2738,Iron:683.7697,Zinc:688.2707,Propionamide:692.8326,Diethyl:685.3453,Methyl:664.3218,Dioxidation:673.3089,Formyl:671.3115,Dehydration:657.3140,Phospho:697.2972,Aluminum+Deamidation:669.7850,Calcium:676.2875" --scan_range 4000,6000
# py GenerateXICGraphs.py --mzml_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\NEW_250402_mEclipse_QC_ncORF-089.mzML" --output_file C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\peptide_089\xic_plots_4ppm\250402_mEclipse_QC_ncORF-089_XICPlots.pdf --modifications "TargetPeptide:657.3140028,Deamidation:657.8060,Oxidation:665.3115,Calcium+Deamidated:676.7795,Aluminum:669.2930,Nickel+Deamidation:685.7659,Sodium:668.3050,Nickel:685.2738,Iron:683.7697,Zinc:688.2707,Propionamide:692.8326,Diethyl:685.3453,Methyl:664.3218,Dioxidation:673.3089,Formyl:671.3115,Dehydration:657.3140,Phospho:697.2972,Aluminum+Deamidation:669.7850,Calcium:676.2875,sodium+deamidation:668.7970" --scan_range 2000,4000 


import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pyteomics import mzml
import numpy as np

def ppm_to_da(mz, ppm):
    """Convert ppm tolerance to Daltons"""
    return mz * ppm / 1e6

def read_xic(mzml_file, target_mz, ppm=4.0, scan_range=None):
    """Extract MS1 intensities for a given precursor m/z ± ppm"""
    scan_numbers, intensities = [], []
    tol_da = ppm_to_da(target_mz, ppm)

    with mzml.read(mzml_file) as reader:
        for spectrum in reader:
            if spectrum['ms level'] == 1:
                scan_number = int(spectrum['id'].split('=')[-1])
                if scan_range and not (scan_range[0] <= scan_number <= scan_range[1]):
                    continue

                mz_array = spectrum['m/z array']
                intensity_array = spectrum['intensity array']

                mask = (mz_array >= target_mz - tol_da) & (mz_array <= target_mz + tol_da)
                xic_intensity = intensity_array[mask].sum()

                scan_numbers.append(scan_number)
                intensities.append(xic_intensity)

    return scan_numbers, intensities

def read_ms2(mzml_file):
    """
    Extract MS2 precursor m/z and the actual MS2 scan number.
    Returns a list of tuples: (precursor_mz, ms1_scan, ms2_scan)
    """
    ms2_marks = []
    with mzml.read(mzml_file) as reader:
        prev_ms1_scan = None
        for spectrum in reader:
            scan_number = int(spectrum['id'].split('=')[-1])
            if spectrum['ms level'] == 1:
                prev_ms1_scan = scan_number
            elif spectrum['ms level'] == 2 and 'precursorList' in spectrum:
                precursor_mz = spectrum['precursorList']['precursor'][0] \
                                            ['selectedIonList']['selectedIon'][0] \
                                            ['selected ion m/z']
                ms2_marks.append((precursor_mz, prev_ms1_scan, scan_number))
    return ms2_marks

def generate_xic_pdf(mzml_file, mod_dict, output_file, scan_range=None, ppm=4.0):
    all_ms2_marks = read_ms2(mzml_file)

    with PdfPages(output_file) as pdf:
        for mod_name, target_mz in mod_dict.items():
            print(f"Processing {mod_name} (m/z={target_mz})...")

            tol_da = ppm_to_da(target_mz, ppm)

            xic_scans, xic_intensity = read_xic(
                mzml_file, target_mz, ppm, scan_range
            )
            scan_to_intensity = dict(zip(xic_scans, xic_intensity))

            # Select MS2 events matching this target m/z
            ms2_events = [
                (ms2_scan, prev_ms1_scan)
                for mz, prev_ms1_scan, ms2_scan in all_ms2_marks
                if abs(mz - target_mz) <= tol_da
            ]

            # X and Y for stars
            offset = 1.05
            stars_x = [ms2_scan for ms2_scan, _ in ms2_events]
            stars_y = [scan_to_intensity.get(prev_ms1_scan, 0) * offset for _, prev_ms1_scan in ms2_events]

            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(xic_scans, xic_intensity, color='purple', label=f'XIC {mod_name}')

            if scan_range:
                ax.set_xlim(scan_range)

            # Plot MS2 stars
            ax.scatter(stars_x, stars_y, color='green', marker='*', s=50)

            # Annotate stars with MS2 scan numbers
            for x, y in zip(stars_x, stars_y):
                ax.annotate(
                    str(x),
                    xy=(x, y),
                    xytext=(0, 5),
                    textcoords='offset points',
                    ha='center',
                    va='bottom',
                    fontsize=8,
                    color='green',
                    rotation=90
                )

            ax.set_xlabel('Scan Number')
            ax.set_ylabel('Intensity')
            ax.set_title(f'XIC for {mod_name}: {target_mz} ± {ppm} ppm')
            ax.grid(True)

            pdf.savefig()
            plt.close()

    print(f"PDF saved to: {output_file}")

def parse_range(range_str):
    return tuple(map(int, range_str.split(',')))

def main():
    parser = argparse.ArgumentParser(
        description="Generate XIC plots from mzML using ppm tolerance"
    )
    parser.add_argument("--mzml_file", required=True)
    parser.add_argument("--output_file", required=True)
    parser.add_argument("--modifications", required=True,
                        help="Comma-separated list of modifications Name:m/z")
    parser.add_argument("--scan_range", default=None,
                        help="Optional scan range, e.g. 3000,6000")
    parser.add_argument("--ppm", type=float, default=4.0,
                        help="Mass tolerance in ppm")

    args = parser.parse_args()

    mod_dict = {}
    for mod in args.modifications.split(','):
        name, mz = mod.split(':')
        mod_dict[name.strip()] = float(mz.strip())

    scan_range = parse_range(args.scan_range) if args.scan_range else None

    generate_xic_pdf(
        args.mzml_file,
        mod_dict,
        args.output_file,
        scan_range,
        args.ppm
    )

if __name__ == "__main__":
    main()
