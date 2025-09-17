aa_masses = {
    'A': 71.037113805,
    'Q': 128.058577540,
    'D': 115.026943065,
    'S': 87.032028435,
    'V': 99.068413945,
    'L': 113.084064015,
    'E': 129.042593135,
    'R': 156.101111050
}

ptm_mass_change = {
    'R': 10.008269
    }

peptide = "AQDSQVLEEER"
water_mass = (2*1.007825035 + 15.99491463)
proton_mass = 1.00727646688
charge = 2

total_mass = 0

def calc_mass(seq, include_water=False):
    total_mass = 0
    for aa in seq:
        base_mass = aa_masses[aa]
        ptm_mass = ptm_mass_change.get(aa, 0)
        total_mass += base_mass + ptm_mass
    if include_water:
        total_mass += water_mass
    return total_mass + proton_mass

for i in range(1, len(peptide)):
    b_fragment = peptide[:i]
    b_mass = calc_mass(b_fragment, include_water=False)
    print(f"b{i:<6} {b_fragment:<12} {b_mass:>14.4f}")

for i in range(1, len(peptide)):
    y_fragment = peptide[-i:]
    y_mass = calc_mass(y_fragment, include_water=True)
    print(f"y{i:<6} {y_fragment:<12} {y_mass:>14.4f}")

