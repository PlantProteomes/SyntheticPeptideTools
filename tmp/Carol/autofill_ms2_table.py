# Some assistance was given by Eric Deutsch and ChatGPT in the making of this program.

import generate_amino_acid_masses
import os
import argparse
import os.path
import csv
import unimod_reader
import usi_checker
import sequence_parser
import timeit
from pyteomics import mzml
from collections import defaultdict

class MS2TableReader:
    """
    TODO: fix find_extra, make find_missing more efficient, instead of finding usi online search in the actual file
    """
    def __init__(self):
        parser = argparse.ArgumentParser(
            prog="MS2TableReader",
            description="Reads and automatically filled MS2 CSV file.")
        parser.add_argument("--ms_run", help="Name of MS run. Example: 250402_mEclipse_QC_ncORF-055")
        parser.add_argument("--sequence", help="Sequence of unmodified peptide (include any tags enclosed in brackets). Example: LTLPAKWER[Label:13C(6)15N(4)]")
        parser.add_argument("--modifications_file", help="Name of modifications CSV file")
        parser.add_argument("--ms2_file", help="Name of CSV file")
        parser.add_argument("--tolerance", help="Tolerance of m/z delta (default 0.002)", default=0.002, type=float)
        parser.add_argument("--pxd", help="Number coming after 'PXD:' in the scan", type=int)
        parser.add_argument("--mzml_file", help="Name of MZML file")
        parser.add_argument("--test_scan_number", help="Scan number to activate verbose mode", type=int)
        parser.add_argument("--test_modification", help="Modification to activate verbose mode")
        args = parser.parse_args()

        if args.mzml_file is None or args.mzml_file == "":
            print('ERROR: Parameter --mzml_file must be provided. See --help for more information')
            return
        if not os.path.isfile(args.mzml_file):
            print(f"ERROR: File '{args.mzml_file}' not found or not a file")
            return

        self.ms_run = args.ms_run
        self.modifications_file = args.modifications_file
        self.sequence = args.sequence
        self.ms2_file = args.ms2_file
        self.tolerance = args.tolerance
        self.pxd = args.pxd
        self.mzml_file = args.mzml_file
        self.test_scan_number = args.test_scan_number
        self.test_modification = args.test_modification
        self.raw_spectra = {}
        self.spectra = []
        self.scores = ""
        self.mass_delta_difference = 0.0
        self.start = timeit.default_timer()
        self.predict_time = 0.0
        self.localization_time = 0.0

    def read_mzml(self):
        with mzml.read(self.mzml_file) as reader:
            for spectrum in reader:
                if spectrum['ms level'] == 2:
                    scan_number = 1 + spectrum['index']
                    mz_array = spectrum['m/z array']
                    mz_dict = defaultdict(list)
                    intensity_array = spectrum['intensity array']
                    intensity_dict = defaultdict(list)
                    for mz in mz_array:
                        key = int(float(mz))
                        mz_dict[key].append(float(mz))
                    for intensity in intensity_array:
                        key = int(float(intensity))
                        intensity_dict[key].append(float(intensity))
                    mz_int_dict = defaultdict(dict)
                    for i in range(mz_array.size):
                        key = int(float(mz_array[i]))
                        mz_int_dict[key][float(mz_array[i])] = float(intensity_array[i])
                    self.raw_spectra[scan_number] = {"scan number" : scan_number,
                                                     "indexed m/z dictionary" : mz_dict,
                                                     "indexed intensity dictionary" : intensity_dict,
                                                     "indexed combined dictionary" : mz_int_dict}

    def predict_modification(self, mass_delta):
        ext = os.path.splitext(self.modifications_file)[1].lower()
        if ext != '.csv':
            print(f"ERROR: Unsupported file extension '{ext}'")

        dipeptide_mass = mass_delta * 2 - 18.0106
        if abs(dipeptide_mass - mass_delta) <= self.tolerance:
            self.mass_delta_difference = dipeptide_mass - mass_delta
            return "dipeptide"

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
        self.mass_delta_difference = closest_mod_val - mass_delta

        if abs(mass_delta - closest_mod_val) <= self.tolerance and most_likely_combination is not None:
            self.mass_delta_difference = closest_mod_val - mass_delta if min(abs(mass_delta - closest_mod_val), abs(mass_delta - most_likely_mass_delta)) == closest_mod_val else most_likely_mass_delta - mass_delta
            return closest_mod if min(abs(mass_delta - closest_mod_val), abs(mass_delta - most_likely_mass_delta)) == closest_mod_val else most_likely_combination
        if abs(mass_delta - closest_mod_val) <= self.tolerance:
            self.mass_delta_difference = closest_mod_val - mass_delta
            return closest_mod
        if abs(mass_delta - most_likely_mass_delta) <= self.tolerance:
            self.mass_delta_difference = most_likely_mass_delta - mass_delta
            return most_likely_combination
        unimod_mod = unimod_reader.search_unimod_by_mass(mass_delta, self.tolerance)
        if unimod_mod != "":
            unimod_val = unimod_reader.search_unimod_by_name(unimod_mod) if unimod_mod != "" else None
            self.mass_delta_difference = unimod_val - mass_delta
        return unimod_mod

    def write_usi(self, scan_number, modification, charge):
        if modification == "":
            return ""
        ms_run = self.ms_run
        new_sequence_list = sequence_parser.parse_sequence(self.sequence)
        if "missing" in modification:
            new_sequence, scores = usi_checker.check_missing(self.sequence, modification, 0.002, self.raw_spectra[scan_number])
            self.scores = scores
            return f"mzspec:PXD{self.pxd}:{ms_run}:scan:{scan_number}:{new_sequence}/{charge}"
        elif "extra" in modification: #adds to n-term
            gain = sequence_parser.parse_sequence(modification.split(" ")[-1])
            for letter in gain:
                new_sequence_list.insert(0,letter)
            new_sequence = "".join(new_sequence_list)
            return f"mzspec:PXD{self.pxd}:{ms_run}:scan:{scan_number}:{new_sequence}/{charge}"
        elif modification == "No mod":
            return f"mzspec:PXD{self.pxd}:{ms_run}:scan:{scan_number}:{self.sequence}/{charge}"
        else:
            new_sequence, scores = usi_checker.check_location(self.sequence, modification, 0.002, self.raw_spectra[scan_number], False)
        if modification != "":
            usi = f"mzspec:PXD{self.pxd}:{ms_run}:scan:{scan_number}:{new_sequence}/{charge}"
            self.scores = scores
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
                                 "ion injection time" : row["injection time"],
                                 "total ion current" : row["total ion current"],
                                 "maximum precursor intensity" : row["maximum precursor intensity"],
                                 "relative intensity" : row["relative intensity"],
                                 "precursor m/z" : row["precursor m/z"],
                                 "precursor charge" : row["precursor charge"],
                                 "precursor mass delta" : float(row["precursor mass delta"]),
                                 "mass delta difference" : "" if "mass delta difference" not in row.keys() else row["mass delta difference"],
                                 "localization scores" : "" if "localization scores" not in row.keys() else row["localization scores"],
                                 "confidence": row["confidence"],
                                 "modification": row["modification"],
                                 "type" : row["type"],
                                 "usi": row["usi"],
                                 "comments": row["comments"]}
                self.spectra.append(spectrum_data)

    def edit_file(self):
        counter = 0
        types = {"missing" : "synthesis error", "extra" : "synthesis error", "Cation:" : "cation", "dipeptide" : "synthesis error"}
        for spectrum in self.spectra:
            counter += 1
            if counter % 1000 == 0:
                print(counter)
            # print("Reading", spectrum["scan number"])
            if spectrum["confidence"] and spectrum["confidence"] != "predicted":
                continue
            t1_start = timeit.default_timer()
            spectrum["modification"] = self.predict_modification(float(spectrum["precursor mass delta"]))
            t1_end = timeit.default_timer()
            self.predict_time += t1_end - t1_start
            if not spectrum["modification"]:
                continue
            t2_start = timeit.default_timer()
            spectrum["usi"] = self.write_usi(int(spectrum["scan number"]), spectrum["modification"], int(spectrum["precursor charge"]))
            t2_end = timeit.default_timer()
            spectrum["mass delta difference"] = f"{self.mass_delta_difference:.4f}" if spectrum["modification"] else ""
            self.mass_delta_difference = 0.0
            spectrum["localization scores"] = self.scores
            self.scores = ""
            self.localization_time += t2_end - t2_start
            spectrum["confidence"] = "predicted"
            for key, value in types.items():
                if key in spectrum["modification"]:
                    spectrum["type"] = value
                    break
        return

    def write_file(self):
        fieldnames = ["file root", "scan number", "scan time", "ion injection time", "total ion current",
                      "maximum precursor intensity", "relative intensity", "precursor m/z", "precursor charge",
                      "precursor mass delta", "mass delta difference", "localization scores", "confidence", "type",
                      "modification", "usi", "comments"]
        with open("autofilled_ms2_table.csv", "w", newline="") as file:
            writer = csv.DictWriter(file, fieldnames=fieldnames)
            writer.writeheader()
            for row in self.spectra:
                writer.writerow(row)

    def test_scan(self, scan_number, modification):
        new_sequence, scores = usi_checker.check_location(self.sequence, modification, 0.002, self.raw_spectra[scan_number],
                                                      True)
        print(new_sequence)
        print(scores)

def main():
    ms2_reader = MS2TableReader()
    print(f"INFO: Reading MZML file {ms2_reader.mzml_file}")
    ms2_reader.read_mzml()
    print(f"INFO: Reading MS2 table {ms2_reader.ms2_file}")
    ms2_reader.read_file()
    print(f"INFO: Predicting modifications and USIs")
    ms2_reader.edit_file()
    print(f"INFO: Writing new MS2 table autofilled_ms2_table.csv")
    ms2_reader.write_file()
    end = timeit.default_timer()
    print(f"Elapsed time: {end - ms2_reader.start}. Time spent finding modifications: {ms2_reader.predict_time}. Time spent localizing: {ms2_reader.localization_time}")
    print(ms2_reader.test_scan(ms2_reader.test_scan_number, ms2_reader.test_modification)) if ms2_reader.test_scan_number and ms2_reader.test_modification else None

if __name__ == "__main__":
    main()