import generate_amino_acid_masses
import os
import argparse
import os.path
import csv
import unimod_reader
import usi_checker
import sequence_parser


class MS2TableReader:
    """
    TODO: missing and extra should also check USIs
    """
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog="MS2TableReader",
            description="Reads and automatically filled MS2 CSV file.")
        parser.add_argument("--ms_run", help="Name of MS run. Example: 250402_mEclipse_QC_ncORF-055")
        parser.add_argument("--sequence", help="Sequence of unmodified peptide (include any tags enclosed in brackets). Example: LTLPAKWER[Label:13C(6)15N(4)]")
        parser.add_argument("--modifications_file", help="Name of modifications CSV file")
        parser.add_argument("--ms2_file", help="Name of CSV file")
        parser.add_argument("--tolerance", help="Tolerance of m/z delta (default 0.01)", default=0.01, type=float)
        args = parser.parse_args()

        self.ms_run = args.ms_run
        self.modifications_file = args.modifications_file
        self.sequence = args.sequence
        self.ms2_file = args.ms2_file
        self.tolerance = args.tolerance
        self.spectra = []

    def predict_modification(self, mass_delta):
        ext = os.path.splitext(self.modifications_file)[1].lower()
        d = ""
        if ext == '.csv':
            d = ','
        else:
            print(f"ERROR: Unsupported file extension '{ext}'")

        lowest_delta = self.tolerance
        most_likely_combination = None
        most_likely_mass_delta = 0
        for i in range(1,6):
            combinations = generate_amino_acid_masses.generate_masses(self.sequence, i)
            for combination in combinations:
                if abs(-mass_delta - combinations[combination]) <= lowest_delta:
                    lowest_delta = abs(-mass_delta - combinations[combination])
                    most_likely_mass_delta = -combinations[combination]
                    most_likely_combination = f"missing {combination}"
                if abs(mass_delta - combinations[combination]) <= lowest_delta:
                    lowest_delta = abs(mass_delta - combinations[combination])
                    most_likely_mass_delta = combinations[combination]
                    most_likely_combination = f"extra {combination}"

        # reads modifications file to check if there is anything
        modifications = {}
        with open(self.modifications_file, newline='') as file:
            reader = csv.DictReader(file, delimiter=',')
            for row in reader:
                if not any(row.values()):
                    continue
                modifications[row["modification"]] = float(row["mass delta"])
        closest_mod, closest_mod_val = min(modifications.items(), key=lambda x: abs(mass_delta - x[1]))

        if abs(mass_delta - closest_mod_val) <= self.tolerance and most_likely_combination is not None:
            return closest_mod if min(abs(mass_delta - closest_mod_val), abs(mass_delta - most_likely_mass_delta)) == closest_mod_val else most_likely_combination
        if abs(mass_delta - closest_mod_val) <= self.tolerance:
            return closest_mod
        if abs(mass_delta - most_likely_mass_delta) <= self.tolerance:
            return most_likely_combination
        return unimod_reader.search_unimod_by_mass(mass_delta, self.tolerance)

    def write_usi(self, scan_number, modification, charge):
        if modification == "":
            return ""
        ms_run = self.ms_run
        new_sequence_list = sequence_parser.parse_sequence(self.sequence)
        if "missing" in modification: # removes first seen instance
            loss = sequence_parser.parse_sequence(modification.split(" ")[-1])
            for letter in loss:
                new_sequence_list.remove(letter)
            new_sequence = "".join(new_sequence_list)
        elif "extra" in modification: #adds to n-term
            gain = sequence_parser.parse_sequence(modification.split(" ")[-1])
            for letter in gain:
                new_sequence_list.insert(0,letter)
            new_sequence = "".join(new_sequence_list)
        elif modification == "No mod":
            return f"mzspec:PXD999005:{ms_run}:scan:{scan_number}:LTLPAKWER[Label:13C(6)15N(4)]/{charge}"
        else:
            new_sequence = usi_checker.check_location(ms_run, scan_number, self.sequence, modification, charge, 0.002)
        if modification != "":
            usi = f"mzspec:PXD999005:{ms_run}:scan:{scan_number}:{new_sequence}/{charge}"
            return usi
        return ""

    def read_file(self):
        ext = os.path.splitext(self.ms2_file)[1].lower()
        d = ","
        if ext != ".csv":
            print(f"ERROR: Unsupported file extension '{ext}'")
        with open(self.ms2_file, newline='') as file:
            reader = csv.DictReader(file, delimiter=d)
            for row in reader:
                if not any(row.values()):
                    continue
                scan_number = int(row['scan number'])
                spectrum_data = {"file root" : row["file root"],
                    "scan number" : scan_number,
                                 "scan time" : row["scan time"],
                                 "total ion current" : row["total ion current"],
                                 "precursor m/z" : row["precursor m/z"],
                                 "precursor charge" : row["precursor charge"],
                                   "precursor mass delta" : float(row["precursor mass delta"]),
                                   "confidence": row["confidence"],
                                   "modification": row["modification"],
                                   "usi": row["usi"],
                                   "comments": row["comments"]}
                self.spectra.append(spectrum_data)

    def edit_file(self):
        for spectrum in self.spectra:
            if spectrum["modification"] == "":
                print(f"Reading {spectrum["scan number"]}")
                spectrum["modification"] = self.predict_modification(float(spectrum["precursor mass delta"]))
                spectrum["usi"] = self.write_usi(int(spectrum["scan number"]), spectrum["modification"], int(spectrum["precursor charge"]))
        return

    def write_file(self):
        fieldnames = ["file root", "scan number", "scan time", "total ion current", "precursor m/z", "precursor charge",
                      "precursor mass delta", "confidence", "modification", "usi", "comments"]
        with open("autofilled_ms2_table.csv", "w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            for row in self.spectra:
                writer.writerow(row)

def main():
    ms2_reader = MS2TableReader()
    ms2_reader.read_file()
    ms2_reader.edit_file()
    ms2_reader.write_file()

if __name__ == "__main__":
    main()