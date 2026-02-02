# example python CalculateCycleTime.py --mzml_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\mzml_files\prm\260113_mEclipse_PRM_ncORF89-AlK(S04)2.mzML" --output_file "C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\mia_data\peptide_089\cycle_times\AlK_cycle_time.xlsx
# 

import argparse
from pyteomics import mzml
import pandas as pd
import numpy as np

def compute_full_ms1_cycles(mzml_file):
    ms1_rts = []
    ms1_scan_numbers = []

    with mzml.read(mzml_file) as reader:
        for spectrum in reader:
            ms_level = spectrum.get('ms level', 1)
            if ms_level != 1:
                continue

            try:
                scan = spectrum['scanList']['scan'][0]
                filterstring = scan.get('filter string', '')
            except (KeyError, IndexError):
                continue

            if "SIM" in filterstring:
                continue

            rt = scan.get('scan start time', None)
            if rt is None:
                continue
            rt_sec = rt * 60 if rt < 100 else rt 

            ms1_rts.append(rt_sec)
            scan_number = int(spectrum.get('id', '0').split('=')[-1])
            ms1_scan_numbers.append(scan_number)

    ms1_rts = np.array(ms1_rts)
    cycle_times = np.diff(ms1_rts)  # delta RT between consecutive full MS1 scans

    return ms1_scan_numbers, ms1_rts, cycle_times

def save_to_excel(scan_numbers, ms1_rts, cycle_times, output_file):
    data = {
        "MS1 Scan #": scan_numbers[:-1],  # last scan has no next MS1
        "MS1 Retention Time (s)": ms1_rts[:-1],
        "Full Cycle Time (s)": cycle_times
    }
    df = pd.DataFrame(data)
    df.to_excel(output_file, index=False)
    print(f"Excel file saved to {output_file}")
    print(f"Average full cycle time: {np.mean(cycle_times):.2f} s")

def main():
    parser = argparse.ArgumentParser(description="Compute true full MS1 cycle times and export to Excel")
    parser.add_argument("--mzml_file", required=True, help="Input mzML file")
    parser.add_argument("--output_file", required=True, help="Output Excel file (.xlsx)")
    args = parser.parse_args()

    scan_numbers, ms1_rts, cycle_times = compute_full_ms1_cycles(args.mzml_file)
    save_to_excel(scan_numbers, ms1_rts, cycle_times, args.output_file)

if __name__ == "__main__":
    main()
