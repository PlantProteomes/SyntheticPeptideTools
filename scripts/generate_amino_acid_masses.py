import itertools
import sequence_parser
import unimod_reader

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

def generate_masses(peptide, r):
    sequence = sequence_parser.parse_sequence(peptide)
    all_mods = sequence_parser.get_mods(peptide)
    mod_dict = {}
    for mod in all_mods:
        mod_dict[mod] = float(unimod_reader.search_unimod_by_name(mod)) if unimod_reader.search_unimod_by_name(mod) else 0.0
    combinations_list = sorted(list(set(itertools.combinations(sorted(sequence), r))))
    output = {}

    for i in range(len(combinations_list)):
        mass = 0
        string = ""
        for j in range(len(combinations_list[i])):
            mods = sequence_parser.get_mods(combinations_list[i][j])
            for mod in mods:
                mass += mod_dict[mod] if mod else 0.0
            mass += acids[sequence_parser.get_acids(combinations_list[i][j])]
            string += combinations_list[i][j]
        output[string] = float(f"{mass:0,.4f}")
    return output

