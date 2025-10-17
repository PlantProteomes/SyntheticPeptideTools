import itertools

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

def generate_masses(peptide, r, modification):
    print("Peptide sequence:", peptide)
    print("Length of subsequence:", r)
    print("Modification:", modification)

    sequence = []
    for letter in peptide:
        sequence.append(letter)

    combinations = set(itertools.combinations(sorted(sequence), r))
    combinations_list = sorted(list(combinations))

    for i in range(len(combinations_list)):
        mass = modifications[modification]
        for j in range(len(combinations_list[i])):
            if combinations_list[i][j] == "R":
                mass += modifications["R tag"]
            mass += acids[combinations_list[i][j]]
            print(combinations_list[i][j], end = "")
        print(f": {mass:0,.4f}")

generate_masses("AQDSQVLEEER", 2, "Acetald +26")