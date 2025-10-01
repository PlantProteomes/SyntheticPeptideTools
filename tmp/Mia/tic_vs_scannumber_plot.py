import csv
from pyteomics import mzml
import matplotlib.pyplot as plt

csv_file = (r"C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\data\NEW_250402_mEclipse_QC_ncORF-089_TIC.csv")

scan_numbers = []
tic_values = []

with open(csv_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        scan_numbers.append(int(row['ScanNumber']))
        tic_values.append(float(row['TIC']))

plt.figure(figsize=(10, 6))
plt.plot(scan_numbers, tic_values, color='#ff00c1')

scans_to_annotate = [2880, 2859]
for scan_to_annotate in scans_to_annotate:
    counter = 0
    for scan in scan_numbers:
        if scan > scan_to_annotate:
            break
        counter += 1
    plt.text(scan_to_annotate, tic_values[counter], f'Scan {scan}', color='black')

plt.xlabel('Scan Number')
plt.ylabel('Total Ion Current (TIC)')
plt.title('TIC vs Scan Number')
plt.grid(True)
plt.savefig(r"C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\data\TIC_over_scannumber_plot(2).pdf")
