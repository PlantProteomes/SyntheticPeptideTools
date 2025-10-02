import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

csv_file = (r"C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\data\NEW_250402_mEclipse_QC_ncORF-089_TIC.csv")

scan_numbers = []
tic_values = []

with open(csv_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        scan_numbers.append(int(row['ScanNumber']))
        tic_values.append(float(row['TIC']))

"""
scans_to_annotate = [2880]

for scan_to_annotate in scans_to_annotate:
    counter = 0
    for scan in scan_numbers:
        if scan > scan_to_annotate:
            break
        counter += 1
    plt.text(scan_to_annotate, tic_values[counter], f'Scan {scan}', color='black')
"""

scans_to_mark = [3750, 4637, 5715, 5883, 6098, 3528, 5554, 2181, 5252, 5099, 5401, 4310, 4135, 4477, 4945, 3954, 2649, 2859, 4793, 3077, 3294]
x_list = [] 
y_list = []

for scan_to_mark in scans_to_mark:  
    counter = 0
    for scan in scan_numbers:
        if scan > scan_to_mark:
            break
        counter += 1
    x_list.append(scan_to_mark)
    y_list.append(tic_values[counter])

x = x_list
y = y_list


with PdfPages(r"C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\data\NEW_tic_over_scan_graph(2).pdf") as pdf:
    plt.figure(figsize=(10, 6))
    plt.plot(scan_numbers, tic_values, color='#ff00c1')
    plt.scatter(x, [val + 2e7 for val in y], color='green', marker='*', s=50)
    plt.xlabel('Scan Number')
    plt.ylabel('Total Ion Current (TIC)')
    plt.title('TIC vs Scan Number')
    plt.grid(True)
    pdf.savefig()
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.plot(scan_numbers, tic_values, color='#ff00c1')
    plt.scatter(x, [val + 2e7 for val in y], color='green', marker='*', s=50)
    plt.xlabel('Scan Number')
    plt.ylabel('Total Ion Current (TIC)')
    plt.title('Zoomed In')
    plt.grid(True)
    plt.xlim(2000, 4000)
    plt.ylim(0,2e9) 
    pdf.savefig()  
    plt.close()

    # plt.savefig(r"C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\data\TIC_over_scannumber_plot(2).pdf")
