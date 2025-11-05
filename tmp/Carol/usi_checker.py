import requests
import json
import pandas as pd
import unimod_reader
import sequence_parser
import math

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
    "No modification" : 0,
    "R tag" : 10.008269,
    "Formyl" : 27.994915,
    "Deamidation" : 0.984016,
    "Acetald +26" : 26.1565,
    "Acetald +28" : 28.0313,
    "Amidation" : -0.984016,
    "Methyl" : 14.01565,
    "Oxidation" : 15.994915,
    "Propionamide" : 71.037114,
    "Butene" : 56.06253364,
    "Carboxy" : 43.989829,
    "Dehydration" : -18.010565
}

def calculate_b_y_ions(modified_sequence):
    """
    TODO: Account for nested brackets.
    :param modified_sequence:
    :return b_y_ions:
    """
    sequence_list = sequence_parser.parse_sequence(modified_sequence)

    b_ions = []
    y_ions = []
    ion_masses = {"b-ions": 1.0073, "y-ions": 1.0073 + 18.0105}
    all_mods = sequence_parser.get_mods(modified_sequence)
    mod_dict = {}
    for mod in all_mods:
        mod_dict[mod] = float(unimod_reader.search_unimod_by_name(mod)) if unimod_reader.search_unimod_by_name(mod) != "" else 0.0
    for i in range(len(sequence_list) - 1):
        for ion_type in ["b", "y"]:
            if ion_type == "b":
                offset = i
            else:
                offset = len(sequence_list) - 1 - i
            amino_acid = sequence_list[offset]
            ion_type_name = f"{ion_type}-ions"
            acid_mods = sequence_parser.get_mods(amino_acid) # extracts anything bracketed
            actual_acid = sequence_parser.get_acids(amino_acid) # extracts anything not bracketed
            mass = acids[actual_acid]
            for mod in acid_mods:
                if mod in mod_dict.keys():
                    mass += mod_dict[mod]
            ion_masses[ion_type_name] += mass

        b_ions.append(float(f"{ion_masses["b-ions"]:.4f}"))
        y_ions.append(float(f"{ion_masses["y-ions"]:.4f}"))

    by_dict = {"b": b_ions, "y": y_ions}
    by_df = pd.DataFrame(by_dict)
    return by_df

def check_location(ms_run, scan_number, sequence, modification, charge, tolerance):
    candidates = {}
    modified_sequence_list = sequence_parser.parse_sequence(sequence)
    modified_sequence_list.insert(0, "[" + modification + "]")
    modified_sequence = "".join(modified_sequence_list)
    url = f"https://proteomecentral.proteomexchange.org/api/proxi/v0.1/spectra?resultType=full&usi=mzspec:PXD999005:{ms_run}:scan:{scan_number}:{modified_sequence}/{charge}"
    response = requests.get(url)
    analysis = json.loads(response.text)
    intensity_array = analysis[0]["intensities"]
    intensity_dict = {}
    for intensity in intensity_array:
        int_intensity = math.floor(float(intensity))
        if int_intensity not in intensity_dict.keys():
            intensity_dict[int_intensity] = [float(intensity)]
        else:
            intensity_dict[int_intensity].append(float(intensity))
    mz_array = analysis[0]["mzs"]
    mz_dict = {}
    for mz in mz_array:
        int_mz = math.floor(float(mz))
        if int_mz not in mz_dict.keys():
            mz_dict[int_mz] = [float(mz)]
        else:
            mz_dict[int_mz].append(float(mz))

    length = len(sequence_parser.parse_sequence(sequence))

    for i in range(length + 1):
        modified_sequence_list = sequence_parser.parse_sequence(sequence)
        modified_sequence_list.insert(i, "[" + modification + "]")
        modified_sequence = "".join(modified_sequence_list)
        by_df = calculate_b_y_ions(modified_sequence)
        match_count = 0
        for b in by_df["b"]:
            b_int = math.floor(b)
            if b_int not in mz_dict.keys():
                continue
            for match in mz_dict[b_int]:
                if abs(match - b) <= tolerance:
                    match_count += 1
        for y in by_df["y"]:
            y_int = math.floor(y)
            if y_int not in mz_dict.keys():
                continue
            for match in mz_dict[y_int]:
                if abs(match - y) <= tolerance:
                    match_count += 1

        candidates[modified_sequence] = match_count
    best_sequence = max(candidates, key=candidates.get)
    return f"{best_sequence}"
