import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

csv_file = r"C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\data\NEW_250402_mEclipse_QC_ncORF-089_TIC.csv"
ms2_csv_file = r"C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\data\ms2_table.csv"
annotations_file = r"C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\repository\SyntheticPeptideTools\data\scan_annotations.xlsx"

scan_numbers, tic_values = [], []
with open(csv_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        scan_numbers.append(int(row['ScanNumber']))
        tic_values.append(float(row['TIC']))

ms2_scan_numbers, ms2_tic_values = [], []
with open(ms2_csv_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        ms2_scan_numbers.append(int(row['scan number']))
        ms2_tic_values.append(float(row['total ion current']))

scan_labels_df = pd.read_excel(annotations_file)
scan_labels = dict(zip(scan_labels_df['ScanNumber'], scan_labels_df['Identification']))

scans_to_mark = [3750, 4637, 5715, 5883, 6098, 3528, 5554, 2181, 5252, 5099, 5401, 
                 4310, 4135, 4477, 4945, 3954, 2649, 2859, 4793, 3077, 3294]
x_list, y_list = [], []
for scan_to_mark in scans_to_mark:
    idx = next((i for i, scan in enumerate(scan_numbers) if scan > scan_to_mark), len(tic_values)-1)
    x_list.append(scan_to_mark)
    y_list.append(tic_values[idx])

scaled_ms2_tic_values = [v * 10 for v in ms2_tic_values]

def annotate_scans(ax, scan_labels, scan_numbers, tic_values, offset=2e7, fontsize=6):
    for scan, label in scan_labels.items():
        idx = next((i for i, s in enumerate(scan_numbers) if s > scan), len(tic_values)-1)
        ax.text(scan, tic_values[idx]+offset, f"{label}", color='black', fontsize=fontsize, rotation=20, ha='center')

with PdfPages(r"C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\repository\SyntheticPeptideTools\data\TIC_over_ScanNumber_Graphs.pdf") as pdf:

    # GRAPH 1: Full MS1 TIC with stars
    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(scan_numbers, tic_values, color='#ff00c1')
    ax.scatter(x_list, [val+2e7 for val in y_list], color='green', marker='*', s=50)
    ax.set_xlabel('Scan Number')
    ax.set_ylabel('Total Ion Current (TIC)')
    ax.set_title('TIC vs Scan Number')
    ax.grid(True)
    pdf.savefig()
    plt.close()

    # GRAPH 2: Zoomed MS1 TIC
    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(scan_numbers, tic_values, color='#ff00c1')
    ax.scatter(x_list, [val+2e7 for val in y_list], color='green', marker='*', s=50)
    ax.set_xlabel('Scan Number')
    ax.set_ylabel('Total Ion Current (TIC)')
    ax.set_title('Zoomed In')
    ax.set_xlim(2000,4000)
    ax.set_ylim(0,2e9)
    ax.grid(True)
    pdf.savefig()
    plt.close()

    # GRAPH 3: MS1 + MS2 overlay
    fig, ax = plt.subplots(figsize=(10,6))
    ax.plot(scan_numbers, tic_values, color='#ff00c1', label='MS1')
    ax.plot(ms2_scan_numbers, scaled_ms2_tic_values, color='blue', label='MS2 x10')
    ax.scatter(x_list, [val+2e7 for val in y_list], color='green', marker='*', s=50)
    ax.set_xlabel('Scan Number')
    ax.set_ylabel('Total Ion Current (TIC)')
    ax.set_title('MS1 and MS2 Overlay')
    ax.set_xlim(2000,4000)
    ax.set_ylim(0,2e9)
    ax.grid(True)
    pdf.savefig()
    plt.close()

    # GRAPH 4-6: MS2 Zoomed with scan number annotations and side legend
    ms2_ranges = [(2350,2750),(2750,2920),(3200,3730)]
    for x_min, x_max in ms2_ranges:
        fig, ax = plt.subplots(figsize=(10,6))
        ax.plot(ms2_scan_numbers, ms2_tic_values, color='blue')
        
        used_y = []

        for scan, identification in scan_labels.items():
            if x_min <= scan <= x_max and scan in ms2_scan_numbers:
                idx = ms2_scan_numbers.index(scan)
                tic_val = ms2_tic_values[idx]
                
                label_y = tic_val * 1.05
                for prev_y in used_y:
                    if abs(prev_y - label_y) < max(ms2_tic_values)*0.03:
                        label_y += max(ms2_tic_values)*0.03
                used_y.append(label_y)
                
                ax.plot([scan, scan], [tic_val, label_y], color='black', linestyle='--', linewidth=0.7)
                ax.text(scan, label_y, f"{scan}", color='black', fontsize=7, rotation=20, ha='center', va='bottom')
        
        legend_text = [f"{scan}: {label}" for scan, label in sorted(scan_labels.items()) if x_min <= scan <= x_max]
        legend_labels = "\n".join(legend_text)
        plt.gcf().text(1.02, 0.5, legend_labels, va='center', fontsize=8, transform=plt.gca().transAxes)
        plt.subplots_adjust(right=0.75)  # Make space for the side legend
        
        ax.set_xlabel('Scan Number')
        ax.set_ylabel('Total Ion Current (TIC)')
        ax.set_title(f'MS2 TIC {x_min}-{x_max}')
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(0, max(used_y)*1.1 if used_y else max(ms2_tic_values))
        ax.grid(True)
        pdf.savefig()
        plt.close()



