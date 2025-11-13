aa_masses = {
    "A": 71.0371,
    "R": 156.1011,
    "N": 114.0429,
    "D": 115.0269,
    "C": 103.0477,
    "Q": 128.0586,
    "E": 129.0426,
    "G": 57.0215,
    "H": 137.0589,
    "I": 113.0841,
    "L": 113.0841,
    "K": 128.0950,
    "M": 131.0405,
    "F": 147.0684,
    "P": 97.0528,
    "S": 87.0320,
    "T": 101.0477,
    "W": 186.0793,
    "Y": 163.0633,
    "V": 99.0684
}

ptm_mass_change = {
    'R': 10.008269,
    'F': 10.027228
}

peptide = "HEEHAHNVNTAF"
water_mass = (2*1.007825035 + 15.99491463)
proton_mass = 1.00727646688
charge = 2

def calc_mass(seq, include_water=False):
    total_mass = 0
    for aa in seq:
        base_mass = aa_masses[aa]
        ptm_mass = ptm_mass_change.get(aa, 0)
        total_mass += base_mass + ptm_mass
    if include_water:
        total_mass += water_mass
    return total_mass + proton_mass

# Calculate full peptide mass
full_mass = calc_mass(peptide, include_water=True)
print(f"Full peptide mass: {full_mass:.4f}")

# Calculate b-ions
for i in range(1, len(peptide)):
    b_fragment = peptide[:i]
    b_mass = calc_mass(b_fragment, include_water=False)
    print(f"b{i:<6} {b_fragment:<12} {b_mass:>14.4f}")

# Calculate y-ions
for i in range(1, len(peptide)):
    y_fragment = peptide[-i:]
    y_mass = calc_mass(y_fragment, include_water=True)
    print(f"y{i:<6} {y_fragment:<12} {y_mass:>14.4f}")
