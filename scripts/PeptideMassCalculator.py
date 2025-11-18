import pandas as pd

acids = {
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

ISOTOPE_SHIFT = {
    "13C": 1.003355,
    "15N": 0.997035,
    "18O": 2.004245,
    "2H": 1.006277
}

def parse_sequence_with_labels(sequence: str):
    import re

    pattern = re.compile(r"([A-Z])(?:\[Label:([^\]]+)\])?")
    parsed = []
    for match in pattern.finditer(sequence):
        aa = match.group(1)
        label = match.group(2)
        shift = 0.0
        if label:
            for iso, count in re.findall(r"(13C|15N|18O|2H)\((\d+)\)", label):
                shift += int(count) * ISOTOPE_SHIFT[iso]
        parsed.append((aa, shift))
    return parsed


def calculate_peptide_mass(sequence: str, charge: int = 1):
    parsed = parse_sequence_with_labels(sequence)
    mass = 18.01056  # Add H2O mass for termini
    for aa, shift in parsed:
        mass += acids[aa] + shift
    mz = (mass + charge * 1.007276) / charge
    return mass, mz


def b_y_ions(sequence: str):
    parsed = parse_sequence_with_labels(sequence)
    n = len(parsed)
    b_ions, y_ions = [], []

    b_mass = 1.007276  # proton
    y_mass = 19.01784  # proton + H2O

    for i in range(n - 1):
        aa, shift = parsed[i]
        b_mass += acids[aa] + shift
        b_ions.append(b_mass)

        aa_y, shift_y = parsed[-(i + 1)]
        y_mass += acids[aa_y] + shift_y
        y_ions.append(y_mass)

    y_ions = y_ions[::-1]  # reverse for matching b-ion order
    return pd.DataFrame({"b-ions": b_ions, "y-ions": y_ions})


# === Example Usage ===
if __name__ == "__main__":
    peptide = "HEEHAHNVNTAF[Label:13C(9)15N(1)]"

    mass, mz = calculate_peptide_mass(peptide, charge=2)
    print(f"Peptide: {peptide}")
    print(f"Monoisotopic mass: {mass:.4f} Da")
    print(f"m/z (z=2): {mz:.4f}")

    ions = b_y_ions(peptide)
    print("\nFragment ions (Da):")
    print(ions)