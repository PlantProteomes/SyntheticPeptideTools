# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import pandas as pd
from pandas import DataFrame
import numpy as np

acids = {
	"A" : 71.0371,
	"R" : 156.1011,
	"N" : 114.0429,
	"D" : 115.0269,
	"C" : 103.0477,
	"Q" : 128.0586,
	"E" : 129.0426,
	"G" : 57.0215,
	"H" : 137.0589,
	"I" : 113.0841,
	"L" : 113.0841,
	"K" : 128.0950,
	"M" : 131.0405,
	"F" : 147.0684,
	"P" : 97.0528,
	"S" : 87.0320,
	"T" : 101.0477,
	"W" : 186.0793,
	"Y" : 163.0633,
	"V" : 99.0684
}

modifications = {
    "267" : 10.008269
}

def mass_to_charge_ratio(sequence, charge):
    mass = 18.006			# extra H - OH
    for letter in sequence:
        if letter == "R":
            mass += acids[letter] + 10.008269
        else:
            mass += acids[letter]
    mass += charge * 1.007		# extra 2 protons
    mz = mass / charge
    return f"{mz:0,.4f}"


def b_y_ions(sequence, charge):
    b_ions = []
    y_ions = []
    n = 0
    ion_masses = {"b-ions" : 1.007, "y-ions" : 1.007 + 18.006}
    for i in range(len(sequence) - 1):
        for ion_type in ["b", "y"]:
            if ion_type == "b":
                offset = i
            else:
                offset = len(sequence) - 1 - i
            amino_acid = sequence[offset]
            ion_type_name = f"{ion_type}-ions"
            mass = acids[amino_acid]
            if amino_acid == "R":
                mass += modifications["267"]

            ion_masses[ion_type_name] += mass

        b_ions.append(f"{ion_masses["b-ions"]:0,.4f}")
        y_ions.append(f"{ion_masses["y-ions"]:0,.4f}")

    by_dict = {"b" : b_ions, "y" : y_ions}
    by_df = pd.DataFrame(by_dict)
    return by_df

print(mass_to_charge_ratio("LTLPAKWER", 2))
print(b_y_ions("LTLPAKWER", 2))