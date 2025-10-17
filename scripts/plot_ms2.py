# Content largely adapted from generate_ms2_table.py.

import os
import argparse
import os.path
import timeit

import pandas as pd
import csv
import numpy as np
import matplotlib.pyplot as plt

class MS2Plot:
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog="PlotMS2",
            description="Plots original and zoomed MS2 spectra from annotated CSV file: total ion current (TIC) vs. scan number.")
        parser.add_argument("--file", help="Name of CSV file")
        parser.add_argument("--xmin", help="Minimum x-axis (scan number) value for zoomed graph. Exponential notation accepted (ex. 1e8)")
        parser.add_argument("--xmax", help="Maximum x-axis (scan number) value for zoomed graph. Exponential notation accepted (ex. 1e8)")
        parser.add_argument("--ymin", help="Minimum y-axis (TIC) value. Exponential notation accepted (ex. 1e8)")
        parser.add_argument("--ymax", help="Maximum y-axis (TIC) value. Exponential notation accepted (ex. 1e8)")

        args = parser.parse_args()

        if args.file is None or args.file == "":
            print('ERROR: Parameter --file must be provided. See --help for more information')
            return

        # tries to open the params file specified
        if not os.path.isfile(args.file):
            print(f"ERROR: File '{args.file}' not found or not a file")
            return

        if args.xmin is None or args.xmin == "" or args.ymin is None or args.ymin == "" or args.xmax is None or args.xmax == "" or args.ymax is None or args.ymax == "":
            print('ERROR: Parameters --xmin, --xmax, --ymin, --ymax must be provided. See --help for more information')
            return

        self.file = args.file
        self.xmin = float(args.xmin)
        self.xmax = float(args.xmax)
        self.ymin = float(args.ymin)
        self.ymax = float(args.ymax)
        self.spectra = []

    def read_file(self):
        ext = os.path.splitext(self.file)[1].lower()
        d = ","

        if ext == '.csv':
            d = ','
        else:
            print(f"ERROR: Unsupported file extension '{ext}'")

        with open(self.file, newline = '') as file:
            reader = csv.DictReader(file, delimiter = d)
            for row in reader:
                if not any(row.values()):
                    continue

                annotation_data = {'scan number' : int(row["scan number"]),
                                   'scan time' : float(row["scan time"]),
                                   'total ion current' : float(row["total ion current"]),
                                   'precursor m/z' : float(row["precursor m/z"]),
                                   "confidence" : row["confidence"],
                                   "modification" : row["modification"],
                                   "usi" : row["usi"],
                                   "comments" : row["comments"]}
                self.spectra.append(annotation_data)

    def plot_tic(self):
        df = pd.DataFrame(self.spectra)
        colors = ["blue"] * len(df)
        markers, stems, base = plt.stem(df["scan number"], df["total ion current"], markerfmt=" ")
        print(df)
        # highlights 20 tallest peaks
        tallest = df.nlargest(20, columns="total ion current")
        for i in range(len(df)):
            if i in tallest.index:
                colors[i] = "red"
        stems.set_colors(colors)

        plt.figure(figsize = (12,9))
        plt.xlabel("Scan number")
        plt.ylabel("Total ion current")
        plt.suptitle("Total Ion Current vs. Scan Number")
        plt.savefig('plot.pdf')

        ax = plt.gca()
        ax.set_xlim(self.xmin, self.xmax)
        ax.set_ylim(self.ymin, self.ymax)
        labels = []
        for i in tallest.index:
            x = tallest['scan number'][i]
            if self.spectra[i]["modification"] != "":
                label = str(x) + "\n" + self.spectra[i]["modification"]
                labels.append(str(x) + ": " + self.spectra[i]["modification"])
            else:
                label = str(x)
            if tallest['total ion current'][i] >= self.ymax:
                if self.spectra[i]["modification"] != "":
                    y = 0.9*self.ymax
                else:
                    y = self.ymax
            else:
                y = tallest['total ion current'][i]
            ax.annotate(label, xy=(x, y), xycoords="data", textcoords="data", size=7, rotation=45, rotation_mode='anchor')
            # rotation = 45, rotation_mode = 'anchor')
        plt.suptitle("Total Ion Current vs. Scan Number (Top 20 Labelled and Annotated)")
        plt.savefig('plot_zoomed_annotated.pdf')

        for i in tallest.index:
            x = tallest['scan number'][i]
            label = str(x)
            if tallest['total ion current'][i] >= self.ymax:
                if self.spectra[i]["modification"] != "":
                    y = 0.9*self.ymax
                else:
                    y = self.ymax
            else:
                y = tallest['total ion current'][i]
            ax.annotate(label, xy=(x, y), xycoords="data", textcoords="data", size=7, rotation=45, rotation_mode='anchor')

        text = "\n".join(labels)
        plt.gcf().text(0.98,0.5, text, fontsize=7)

        plt.savefig('plot_zoomed_legend.pdf')


def main():
    ms2_plot = MS2Plot()
    print(f"INFO: Reading CSV File {ms2_plot.file}")
    ms2_plot.read_file()
    print(f"INFO: Plotting to plot.png and plot_zoomed.png")
    ms2_plot.plot_tic()

if __name__ == "__main__":
    main()