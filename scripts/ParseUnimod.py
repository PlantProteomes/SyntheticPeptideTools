# mod_dict = parse_unimod("path/to/unimod.obo")

import re

def parse_unimod(obo_file):

    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
    termini = ["N-term", "C-term", "Protein N-term", "Protein C-term"]
    positions = amino_acids + termini

    mod_dict = {}

    with open(obo_file, "r", encoding="utf-8") as f:
        lines = f.readlines()

    lines_iter = iter(lines)
    current = None

    for line in lines_iter:
        line = line.strip()

        if line == "[Term]":
            if current and current["Modification"] and current["Monoisotopic Mass"] is not None:
                mass=current["Monoisotopic Mass"]
                if mass not in mod_dict:
                    mod_dict[mass] = []
                mod_dict[mass].append(current)
            current = {"Modification": None, "Monoisotopic Mass": None}
            for pos in positions:
                current[pos] = ""  # initialize all positions as empty

        elif current is not None:
            # Get modification name
            if line.startswith("name: "):
                current["Modification"] = line[6:]

            # Get monoisotopic mass
            elif line.startswith("xref: delta_mono_mass "):
                m = re.match(r'xref: delta_mono_mass "(.+)"', line)
                if m:
                    try:
                        current["Monoisotopic Mass"] = float(m.group(1))
                    except:
                        pass
                    
            # Get sites
            elif line.startswith("xref: spec_") and "_site" in line:
                m = re.match(r'xref: spec_(\d+)_site "(.+)"', line)
                if m:
                    spec_num = m.group(1)
                    site = m.group(2)

                    # Look ahead for the corresponding position
                    position = None
                    prefix = f"xref: spec_{spec_num}_position "
                    for peek_line in lines_iter:
                        peek_line = peek_line.strip()
                        if peek_line.startswith(prefix):
                            pm = re.match(rf'{prefix}"(.+)"', peek_line)
                            if pm:
                                position = pm.group(1)
                            break

                    # Map the site + position to dictionary columns
                    if position:
                        if site in amino_acids and position == "Anywhere":
                            current[site] = "yes"
                        elif site == "N-term":
                            if position == "Any N-term":
                                current["N-term"] = "yes"
                            elif position == "Protein N-term":
                                current["Protein N-term"] = "yes"
                        elif site == "C-term":
                            if position == "Any C-term":
                                current["C-term"] = "yes"
                            elif position == "Protein C-term":
                                current["Protein C-term"] = "yes"

    # Save the last modification
    if current and current["Modification"] and current["Monoisotopic Mass"] is not None:
        mass = current["Monoisotopic Mass"]
        if mass not in mod_dict:
            mod_dict[mass] = []
        mod_dict[mass].append(current)
    return mod_dict



if __name__ == "__main__":
    obo_path = r"C:\Users\miawc\OneDrive\Documents\ISB_INTERNSHIP\projects\unimod.obo.txt" # replace with path to your unimod.obo file
    mod_dict = parse_unimod(obo_path)

    print("Done parsing")

# test
mass_to_check = 42.010565
if mass_to_check in mod_dict:
    print(mod_dict[mass_to_check])