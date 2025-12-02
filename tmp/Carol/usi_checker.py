import requests
import json
import pandas as pd
import unimod_reader
import sequence_parser
import itertools
import timeit
from collections import defaultdict

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
    Returns DataFrame of b- and y-ions given a sequence of amino acids including modifications.
    :param str modified_sequence:
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

def check_location(sequence, modification, tolerance, spectrum, verbosity):
    # TODO: Weight by intensity
    unmodded_by_df = calculate_b_y_ions(sequence)
    candidates = {}
    approval_list = []
    modified_sequence_list = sequence_parser.parse_sequence(sequence)
    possible_locations = unimod_reader.localize(modification)
    mz_dict = spectrum['indexed m/z dictionary']
    intensity_dict = spectrum['indexed intensity dictionary']
    max_intensity = max(value for list in intensity_dict.values() for value in list)
    combined_dict = spectrum['indexed combined dictionary']
    code_string = ""

    length = len(sequence_parser.parse_sequence(sequence))
    mod_mass = float(unimod_reader.search_unimod_by_name(modification)) if unimod_reader.search_unimod_by_name(modification) != "" else 0.0
    is_approved = False

    for i in range(length + 1):
        intensity_sum = 0.0
        is_approved = False
        current_sequence = modified_sequence_list.copy()
        current_sequence.insert(i, "[" + modification + "]")
        current_by_df = unmodded_by_df.copy()
        match_score = 0.0
        match_count = 0
        print("current sequence", current_sequence) if verbosity else None

        if possible_locations is not None:
            if i == 0:
                if "N-term" in possible_locations:
                    is_approved = True
                    code_string += "NH2-"
                    print(f"modification {modification} on N-term approved by unimod") if verbosity else None
                else:
                    code_string += "nh2-"
            elif i == length + 1:
                if "C-term" in possible_locations:
                    is_approved = True
                    code_string += "COOH-"
                    print(f"modification {modification} on C-term approved by unimod") if verbosity else None
                else:
                    code_string += "cooh-"
            else:
                if modified_sequence_list[i-1] in possible_locations:
                    is_approved = True
                    code_string += f"{modified_sequence_list[i - 1]}-"
                    print(f"modification {modification} on {modified_sequence_list[i-1]} approved by unimod") if verbosity else None
                else:
                    code_string += f"{modified_sequence_list[i - 1].lower()}-"

        for row in current_by_df.itertuples():
            index = row.Index
            b = float(row.b)
            y = float(row.y)
            if index + 1 >= i:
                b += mod_mass
            b_int = int(b)
            if b_int in combined_dict:
                for match in combined_dict[b_int].keys():
                    if abs(match - b) <= tolerance:
                        print(f"found matching b-ion with mass {match} and intensity {combined_dict[b_int][match]}") if verbosity else None
                        intensity_sum += combined_dict[b_int][match]
                        match_count += 1
            if length - index <= i:
                y += mod_mass
            y_int = int(y)
            if y_int not in combined_dict:
                continue
            for match in combined_dict[y_int].keys():
                if abs(match - y) <= tolerance:
                    print(f"found matching y-ion with mass {match} and intensity {combined_dict[y_int][match]}") if verbosity else None
                    intensity_sum += combined_dict[y_int][match]
                    match_count += 1
        match_score = 5 * (intensity_sum / (match_count * max_intensity)) + 5 * (match_count / (2 * length - 2)) if match_count * max_intensity != 0 else 5 * (match_count / (2 * length - 2))
        code_string += f"{match_score:.2f}," if i < length else f"{match_score:.2f}"
        candidates["".join(current_sequence)] = match_score
        approval_list.append(is_approved)

    best_sequence = max(candidates, key=candidates.get)
    return f"{best_sequence}", code_string

def check_missing(sequence, modification, tolerance, spectrum):
    loss = sequence_parser.parse_sequence(modification.split(" ")[-1])
    sequence_list = sequence_parser.parse_sequence(sequence)
    mz_dict = spectrum['indexed m/z dictionary']
    intensity_dict = spectrum['indexed intensity dictionary']
    max_intensity = max(value for list in intensity_dict.values() for value in list)
    combined_dict = spectrum['indexed combined dictionary']
    code_string = ""

    letter_indices = {}
    for index, chara in enumerate(sequence_list):
        letter_indices.setdefault(chara, []).append(index)

    removal_indices = []
    for chara in loss:
        removal_indices.append(letter_indices[chara])

    candidates = {}
    known_products = set()
    for product in itertools.product(*removal_indices): # (2, 3, 1)
        # This can be optimized better for speed because itertools.product generates a lot of unnecessary stuff
        if len(set(product)) != len(product):
            continue
        product_set = frozenset(product)
        if product_set in known_products:
            continue
        known_products.add(product_set)
        modified_sequence_list = sequence_parser.parse_sequence(sequence)
        length = len(modified_sequence_list)
        sorted_product = sorted(product, reverse=True)
        modified_sequence_display = modified_sequence_list.copy()
        for i in range(len(sorted_product)):
            modified_sequence_display[sorted_product[i]] = "_"
            modified_sequence_list.pop(sorted_product[i])
        modified_sequence = "".join(modified_sequence_list)
        code_string += f"{"".join(modified_sequence_display)}-"
        by_df = calculate_b_y_ions(modified_sequence)
        match_score = 0.0
        match_count = 0
        intensity_sum = 0.0

        for row in by_df.itertuples():
            b = float(row.b)
            y = float(row.y)
            b_int = int(b)
            y_int = int(y)
            if b_int in combined_dict:
                for match in combined_dict[b_int].keys():
                    if abs(match - b) <= tolerance:
                        intensity_sum += combined_dict[b_int][match]
                        match_count += 1
            if y_int not in combined_dict:
                continue
            for match in combined_dict[y_int].keys():
                if abs(match - y) <= tolerance:
                    intensity_sum += combined_dict[y_int][match]
                    match_count += 1

        match_score = 5 * (intensity_sum / (match_count * max_intensity)) + 5 * (match_count / (2 * length - 2)) if match_count * max_intensity != 0 else 5 * (match_count / (2 * length - 2))
        code_string += f"{match_score:.2f},"
        candidates[modified_sequence] = match_score
    code_string = code_string[:-1]
    best_sequence = max(candidates, key=candidates.get)
    return f"{best_sequence}", code_string

def check_extra(ms_run, scan_number, sequence, modification, pxd, tolerance):
    gain = sequence_parser.parse_sequence(modification.split(" ")[-1])
    sequence_list = sequence_parser.parse_sequence(sequence)

    url = f"https://proteomecentral.proteomexchange.org/api/proxi/v0.1/spectra?resultType=full&usi=mzspec:PXD{pxd}:{ms_run}:scan:{scan_number}"
    response = requests.get(url)
    analysis = json.loads(response.text)
    intensity_array = analysis[0]["intensities"]
    intensity_dict = {}
    for intensity in intensity_array:
        int_intensity = int(float(intensity))
        if int_intensity not in intensity_dict.keys():
            intensity_dict[int_intensity] = [float(intensity)]
        else:
            intensity_dict[int_intensity].append(float(intensity))
    mz_array = analysis[0]["mzs"]
    mz_dict = {}
    for mz in mz_array:
        int_mz = int(float(mz))
        if int_mz not in mz_dict.keys():
            mz_dict[int_mz] = [float(mz)]
        else:
            mz_dict[int_mz].append(float(mz))

    removal_indices = {}
    for i in range(len(sequence_list)):
        if sequence_list[i] in gain:
            if sequence_list[i] not in removal_indices.keys():
                removal_indices[sequence_list[i]] = [i]
            else:
                removal_indices[sequence_list[i]].append(i)

    candidates = {}
    for product in itertools.product(*removal_indices.values()):  # (2, 3, 1)
        modified_sequence_list = sequence_parser.parse_sequence(sequence)
        sorted_product = sorted(product, reverse=True)
        for i in range(len(sorted_product)):
            modified_sequence_list.pop(sorted_product[i])
        modified_sequence = "".join(modified_sequence_list)
        by_df = calculate_b_y_ions(modified_sequence)
        match_count = 0
        for b in by_df["b"]:
            b_int = int(b)
            if b_int not in mz_dict.keys():
                continue
            for match in mz_dict[b_int]:
                if abs(match - b) <= tolerance:
                    print(f"Found match at {match}")
                    match_count += 1
        for y in by_df["y"]:
            y_int = int(y)
            if y_int not in mz_dict.keys():
                continue
            for match in mz_dict[y_int]:
                if abs(match - y) <= tolerance:
                    match_count += 1
        candidates[modified_sequence] = match_count

    best_sequence = max(candidates, key=candidates.get)
    return f"{best_sequence}"

# snippet = '''
# url = f"https://proteomecentral.proteomexchange.org/api/proxi/v0.1/spectra?resultType=full&usi=mzspec:PXD999007:251103_mEclipse_ncORF89-S1:scan:5001"
# response = requests.get(url)
# analysis = json.loads(response.text)
# '''


# if __name__ == "__main__":
#     print(timeit.timeit(stmt="check_missing('251103_mEclipse_ncORF89-S1', 4343, 'AQDSQVLEEER[Label:13C(6)15N(4)]', 'missing QV', 999007, 0.002)",setup="from __main__ import check_missing", number=10))
#     print(timeit.timeit(stmt="check_location('251103_mEclipse_ncORF89-S1', 5001, 'AQDSQVLEEER[Label:13C(6)15N(4)]', 'Asp->Gly', 999007, 0.002)",setup="from __main__ import check_location", number=10))
#     print(timeit.timeit(stmt=snippet, globals=globals(), number=10))

# print(check_location('AQDSQVLEEER[Label:13C(6)15N(4)]', 'Oxidation', 999007, 0.002))