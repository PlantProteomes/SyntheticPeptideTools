import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

csv_file = (r"C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\data\NEW_250402_mEclipse_QC_ncORF-089_TIC.csv")
ms2_csv_file = (r"C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\data\ms2_table.csv")

scan_numbers = []
tic_values = []

with open(csv_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        scan_numbers.append(int(row['ScanNumber']))
        tic_values.append(float(row['TIC']))

ms2_scan_numbers = []
ms2_tic_values = []

with open(ms2_csv_file, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        ms2_scan_numbers.append(int(row['scan number']))
        ms2_tic_values.append(float(row['total ion current']))

scan_labels = {
    2880: "deamidation",
    2748: "missing Q",
    2770: "dipeptide",
    2766: "dipeptide",
    2767: "loss of S",
    2703: "Gln->Gly",
    2472: "missing EE",
    2762: "dipeptide",
    2525: "missing AQD + Propionamide",
    2812: "missing AQ",
    2810: "Propionamide",
    3616: "missing AQD + biotin",
    2701: "missing AQ + SMA",
    3690: "Butene",
    3335: "missing R",
    2709: "Glu->Ala",
    2366: "missing AQDS",
    3704: "oxidation and carboxylation",
    2804: "Propionamide",
    2684: "missing AQ & Ethylphosphate",
    2873: "deamidation + Cation:K", 
    2463: "missing EE",
    2799: "acetald+26",
    2704: "missing AQD + nipcam", 
    2895: "extra E",
    2441: "missing AQD", 
    2820: "methyl", 
    3434: "missing AQD + diethyl",
    2616: "missing AE",
    2801: "deamidation + replacement of two protons with nickel",
    2828: "missing AQD + Carboxymethyl", 
    3294: "deamidation + amidation",
    2695: "missing AQD; addition of phosphate + hexose",
    2875: "Al cation", 

}

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

scaled_ms2_tic_values = [v * 10 for v in ms2_tic_values]

with PdfPages(r"C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\data\NEW_tic_over_scan_graph(5).pdf") as pdf:
# GRAPH 1: Full mS1 TIC vs Scan Number with Stars
    plt.figure(figsize=(10, 6))
    plt.plot(scan_numbers, tic_values, color='#ff00c1')
    plt.scatter(x, [val + 2e7 for val in y], color='green', marker='*', s=50)
    plt.xlabel('Scan Number')
    plt.ylabel('Total Ion Current (TIC)')
    plt.title('TIC vs Scan Number')
    for scan_to_annotate, label in scan_labels.items():
        counter = 0
        for scan in scan_numbers:
            if scan > scan_to_annotate:
                break
            counter += 1
        if counter >= len(tic_values):  
            counter = len(tic_values) - 1
    plt.grid(True)
    pdf.savefig()
    plt.close()

# GRAPH 2: Zoomed in MS1 TIC vs Scan Number with Stars
    plt.figure(figsize=(10, 6))
    plt.plot(scan_numbers, tic_values, color='#ff00c1')
    plt.scatter(x, [val + 2e7 for val in y], color='green', marker='*', s=50)
    plt.xlabel('Scan Number')
    plt.ylabel('Total Ion Current (TIC)')
    plt.title('Zoomed In')
    plt.grid(True)
    for scan_to_annotate, label in scan_labels.items():
        counter = 0
        for scan in scan_numbers:
            if scan > scan_to_annotate:
                break
            counter += 1
        if counter >= len(tic_values):  
            counter = len(tic_values) - 1
        plt.text(scan_to_annotate, tic_values[counter] + 2e7, f"{scan_to_annotate}", color='black', fontsize=6, rotation = 20)
    plt.xlim(2000, 4000)
    plt.ylim(0,2e9) 
    legend_text = [f"{scan}: {label}" for scan, label in sorted(scan_labels.items())]
    legend_labels = "\n".join(legend_text)
    plt.gcf().text(1.02, 0.5, legend_labels, va='center', fontsize=5, transform=plt.gca().transAxes)
    plt.subplots_adjust(right=0.75)
    pdf.savefig()  
    plt.close()

# GRAPH 3: MS1 and MS2 graph overlayed with Annotations 
    plt.figure(figsize=(10, 6))
    plt.plot(scan_numbers, tic_values, color='#ff00c1')
    plt.plot(ms2_scan_numbers, scaled_ms2_tic_values, color='blue')
    plt.scatter(x, [val + 2e7 for val in y], color='green', marker='*', s=50)
    plt.xlabel('Scan Number')
    plt.ylabel('Total Ion Current (TIC)')
    plt.title('TIC vs Scan Number')
    plt.grid(True)
    for scan_to_annotate, label in scan_labels.items():
        counter = 0
        for scan in scan_numbers:
            if scan > scan_to_annotate:
                break
            counter += 1
        if counter >= len(tic_values):  
            counter = len(tic_values) - 1
        plt.text(scan_to_annotate, tic_values[counter] + 2e7, f"{scan_to_annotate}", color='black', fontsize=6, rotation = 20)
    plt.xlim(2000, 4000)
    plt.ylim(0,2e9) 
    legend_text = [f"{scan}: {label}" for scan, label in sorted(scan_labels.items())]
    legend_labels = "\n".join(legend_text)
    plt.gcf().text(1.02, 0.5, legend_labels, va='center', fontsize=5, transform=plt.gca().transAxes)
    plt.subplots_adjust(right=0.75)
    pdf.savefig()
    plt.close()

# GRAPH 4: MS2 TIC vs Scan Number Graph 2350 to 2700
    plt.figure(figsize=(10, 6))
    plt.plot(ms2_scan_numbers, ms2_tic_values, color='blue')
    used_label_positions = []
    x_min, x_max = 2350, 2750
    plt.xlim(x_min, x_max)
    for scan_to_annotate, label in scan_labels.items():
        if x_min <= scan_to_annotate <= x_max:
            index = ms2_scan_numbers.index(scan_to_annotate)
            tic_val = ms2_tic_values[index]
            label_y = tic_val * 1.05
            for prev_y in used_label_positions:
                if abs(prev_y - label_y) < max(ms2_tic_values) * 0.03: 
                    label_y += max(ms2_tic_values) * 0.03
            used_label_positions.append(label_y)
            plt.plot([scan_to_annotate, scan_to_annotate], [tic_val, label_y],color='black', linestyle='--', linewidth=0.7)
            plt.text(scan_to_annotate, label_y, f"{scan_to_annotate}", color='black', fontsize=7, rotation=20, ha='center', va='bottom')

    plt.xlabel('Scan Number')
    plt.ylabel('Total Ion Current (TIC)')
    plt.title('MS2 TIC')
    plt.grid(True)
    plt.xlim(x_min, x_max)
    plt.ylim(0, max(used_label_positions) * 1.1) 

    legend_text = [f"{scan}: {label}" for scan, label in sorted(scan_labels.items())]
    legend_labels = "\n".join(legend_text)
    plt.gcf().text(1.02, 0.5, legend_labels, va='center', fontsize=5, transform=plt.gca().transAxes)
    plt.subplots_adjust(right=0.75)
    pdf.savefig()
    plt.close()

# GRAPH 5: MS2 TIC vs Scan Number Graph 2750 to 2920
    plt.figure(figsize=(10, 6))
    plt.plot(ms2_scan_numbers, ms2_tic_values, color='blue')
    used_label_positions = []
    x_min, x_max = 2750, 2920
    plt.xlim(x_min, x_max)
    for scan_to_annotate, label in scan_labels.items():
        if x_min <= scan_to_annotate <= x_max:
            index = ms2_scan_numbers.index(scan_to_annotate)
            tic_val = ms2_tic_values[index]
            label_y = tic_val * 1.05
            for prev_y in used_label_positions:
                if abs(prev_y - label_y) < max(ms2_tic_values) * 0.03: 
                    label_y += max(ms2_tic_values) * 0.03
            used_label_positions.append(label_y)
            plt.plot([scan_to_annotate, scan_to_annotate], [tic_val, label_y],color='black', linestyle='--', linewidth=0.7)
            plt.text(scan_to_annotate, label_y, f"{scan_to_annotate}", color='black', fontsize=7, rotation=20, ha='center', va='bottom')
    plt.xlabel('Scan Number')
    plt.ylabel('Total Ion Current (TIC)')
    plt.title('MS2 TIC')
    plt.grid(True)
    plt.xlim(2700, 2900)
    plt.ylim(0, max(used_label_positions) * 1.1) 

    legend_text = [f"{scan}: {label}" for scan, label in sorted(scan_labels.items())]
    legend_labels = "\n".join(legend_text)
    plt.gcf().text(1.02, 0.5, legend_labels, va='center', fontsize=5, transform=plt.gca().transAxes)
    plt.subplots_adjust(right=0.75)
    pdf.savefig()
    plt.close()

# GRAPH 6: MS2 TIC vs Scan Number Graph 3200 to 3730
    plt.figure(figsize=(10, 6))
    plt.plot(ms2_scan_numbers, ms2_tic_values, color='blue')
    used_label_positions = []
    x_min, x_max = 3200, 3730
    plt.xlim(x_min, x_max)
    for scan_to_annotate, label in scan_labels.items():
        if x_min <= scan_to_annotate <= x_max:
            index = ms2_scan_numbers.index(scan_to_annotate)
            tic_val = ms2_tic_values[index]
            label_y = tic_val * 1.05
            for prev_y in used_label_positions:
                if abs(prev_y - label_y) < max(ms2_tic_values) * 0.03: 
                    label_y += max(ms2_tic_values) * 0.03
            used_label_positions.append(label_y)
            plt.plot([scan_to_annotate, scan_to_annotate], [tic_val, label_y],color='black', linestyle='--', linewidth=0.7)
            plt.text(scan_to_annotate, label_y, f"{scan_to_annotate}", color='black', fontsize=7, rotation=20, ha='center', va='bottom')
    plt.xlabel('Scan Number')
    plt.ylabel('Total Ion Current (TIC)')
    plt.title('MS2 TIC')
    plt.grid(True)
    plt.xlim(3200, 3710)
    plt.ylim(0, max(used_label_positions) * 1.1) 

    legend_text = [f"{scan}: {label}" for scan, label in sorted(scan_labels.items())]
    legend_labels = "\n".join(legend_text)
    plt.gcf().text(1.02, 0.5, legend_labels, va='center', fontsize=5, transform=plt.gca().transAxes)
    plt.subplots_adjust(right=0.75)
    pdf.savefig()
    plt.close()
