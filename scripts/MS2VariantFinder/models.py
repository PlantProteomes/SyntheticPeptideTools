import itertools
from dataclasses import dataclass
import unimod
from constants import *
import numpy as np
from collections import defaultdict

class Scan:
    def __init__(
            self,
            scan_number: int,
            scan_type: str,
            ms_level: int,
            isolation_window: tuple,
            mz_array: np.ndarray,
            intensity_array: np.ndarray,
            rt: float,
            iit: float,
            tic: float,
            precursor_mz: float | None = None,
            precursor_charge: int | None = None,
            last_ms1_scan: int | None = None
    ):
        self.scan_number = scan_number
        self.scan_type = scan_type
        self.ms_level = ms_level
        self.isolation_window = isolation_window
        self.mz_array = mz_array
        self.intensity_array = intensity_array
        self.rt = rt
        self.iit = iit
        self.tic = tic

        if self.ms_level == 2:
            self.precursor_mz = precursor_mz
            self.precursor_charge = precursor_charge
            self.last_ms1_scan = last_ms1_scan

class MSRun:
    def __init__(self, scans: list[Scan], run_type):
        self.scans = scans
        self.run_type = run_type
        self.indexed_scans = {scan.scan_number: scan for scan in scans}
        self.ms1_spectra = [scan for scan in scans if scan.ms_level == 1]
        self.ms2_spectra = [scan for scan in scans if scan.ms_level == 2]

    def get_scan(self, scan_number):
        if scan_number in self.indexed_scans:
            return self.indexed_scans[scan_number]
        raise KeyError(f"Scan {scan_number} not found")

    def get_sim_scans(self):
        ms1_sim_spectra = defaultdict(list)
        for scan in self.ms1_spectra:
            if scan.scan_type == "SIM":
                ms1_sim_spectra[scan.isolation_window].append(scan)
        return ms1_sim_spectra

    def get_precursor(self, scan):
        if scan.ms_level == 1:
            raise AttributeError("MS1 scan has no precursor.")

        if self.run_type == "PRM":
            sim = self.get_sim_scans()
            for key in sim.keys():
                if key[0] < scan.precursor_mz < key[1]:
                    precursor = None
                    for sim_scan in sim[key]:
                        if sim_scan.scan_number >= scan.scan_number:
                            break
                        precursor = sim_scan
                    return precursor

        if self.run_type == "DDA":
            return self.indexed_scans[scan.last_ms1_scan]

        return None

@dataclass
class Modification:
    position: int
    delta: float
    name: str
    is_labile: bool = False

    @classmethod
    def from_string(cls, pos, mod_name, is_labile=False):
        return cls(
            position=pos,
            delta=float(unimod.get_mod(mod_name)["delta_mono_mass"]),
            name=mod_name,
            is_labile=is_labile
        )

    @property
    def is_n_term(self):
        return self.position == -1

class Peptide:
    def __init__(self, raw_sequence: str, modifications: list[Modification]):
        self.raw_sequence = raw_sequence
        self.modifications = modifications

    @classmethod
    def from_string(cls, peptide):
        bracket_depth = 0
        modification = ""
        raw_sequence = ""
        previous_chara = ""
        modifications = []
        counter = -1

        for chara in peptide:
            if chara == "[":
                if bracket_depth == 0:
                    previous_chara = raw_sequence[counter] if raw_sequence != "" else None
                    bracket_depth += 1
                    continue
                bracket_depth += 1
            if chara == "]":
                bracket_depth -= 1
            if chara == "{":
                if bracket_depth == 0:
                    previous_chara = "Labile"
                    bracket_depth += 1
                    continue
                bracket_depth += 1
            if chara == "}":
                bracket_depth -= 1
            if bracket_depth != 0:
                modification += chara
            if bracket_depth == 0 and modification != "":
                if previous_chara == "Labile":
                    modifications.append(Modification.from_string(counter, modification, is_labile=True))
                elif counter == -1:
                    modifications.append(Modification.from_string(counter, modification))
                else:
                    modifications.append(Modification.from_string(counter, modification))
                modification = ""
            if chara.isupper() and modification == "":
                counter += 1
                raw_sequence += chara
        return cls(raw_sequence, modifications)

    def remove_residues(self, indices):
        indices_set = set(indices)
        modded_sequence = "".join([char for idx, char in enumerate(self.raw_sequence) if idx not in indices_set])
        new_modifications = []
        for mod in self.modifications:
            if mod.position in indices_set:
                continue
            if mod.position == -1 and 0 in indices_set:
                continue
            offset = sum(1 for i in indices if i < mod.position)
            new_modifications.append(Modification(
                position=mod.position - offset,
                delta=mod.delta,
                name=mod.name,
                is_labile=mod.is_labile
            ))
        return Peptide(modded_sequence, new_modifications)

    def insert_residues(self, residue_index, index):
        raw = list(self.raw_sequence)
        raw.insert(index, self.get_residue(residue_index))
        modded_sequence = "".join(raw)
        new_modifications = []
        for mod in self.modifications:
            new_pos = mod.position + 1 if mod.position >= index else mod.position
            new_modifications.append(Modification(
                position=new_pos,
                delta=mod.delta,
                name=mod.name,
                is_labile=mod.is_labile
            ))
            if mod.position == residue_index:
                new_modifications.append(Modification(
                    position=index,
                    delta=mod.delta,
                    name=mod.name,
                    is_labile=mod.is_labile
                ))
        return Peptide(modded_sequence, new_modifications)

    def mass(self, charge):
        mass = WATER_MASS
        mass += sum(AA_MASSES[acid] for acid in self.raw_sequence)
        mass += sum(mod.delta for mod in self.modifications)
        mass += PROTON_MASS * charge
        return mass

    def mz(self, charge):
        return self.mass(charge) / charge

    # returns fragment ion masses
    def fragments(self):
        fragments = {"a": self.a_ions(), "b": self.b_ions(1), "b2": self.b_ions(2), "y": self.y_ions(1), "y2": self.y_ions(2)}
        return fragments

    # returns a-ion masses
    def a_ions(self):
        a_ions = []
        for i in range(1, len(self)):
            a_ions.append(Fragment(self, "a", i, 1))
        return a_ions

    # returns b- and b2-ion masses
    def b_ions(self, charge):
        b_ions = []
        for i in range(1, len(self)):
            b_ions.append(Fragment(self, "b", i, charge))
        return b_ions

    # returns y- and y2-ion masses
    def y_ions(self, charge):
        y_ions = []
        for i in range(1, len(self)):
            y_ions.append(Fragment(self, "y", i, charge))
        return y_ions

    def get_residue(self, position):
        if position == -1:
            return "N-term"
        return self.raw_sequence[position]

    def generate_error_masses(self, r, deduplicate=False):
        combinations_list = list(itertools.combinations(range(len(self)), r))
        output = {}
        for combination in combinations_list:
            mass = 0
            indices = []
            parts = []
            for j in combination:
                part = self.get_residue(j)
                indices.append(j)
                mass += AA_MASSES[self.get_residue(j)]
                for mod in self.modifications:
                    if mod.position == j:
                        mass += mod.delta
                        part += f"[{mod.name}]"
                parts.append(part)
            key = tuple(frozenset(sorted(indices))) if deduplicate else tuple(sorted(indices))
            output[key] = float(mass)
        return output

    def __str__(self):
        output = ""
        for i in range(-1, len(self.raw_sequence)):
            output += self.raw_sequence[i] if i != -1 else ""
            for modification in self.modifications:
                if modification.position == i:
                    if modification.is_labile:
                        output += f"{{{modification.name}}}"
                    else:
                        output += f"[{modification.name}]"
        return output

    def __len__(self):
        return len(self.raw_sequence)


class Fragment:
    def __init__(self, full_peptide, ion_type, position, charge):
        self.full_sequence = full_peptide.raw_sequence
        self.modifications = full_peptide.modifications
        self.ion_type = ion_type
        self.position = position
        if self.position < 0 or self.position > len(self.full_sequence):
            raise ValueError("Invalid fragment ion position.")
        self.charge = charge
        if ion_type == "b" or ion_type == "a":
            self.valid_mods = [mod for mod in self.modifications if (mod.position < position and not mod.is_labile)]
        elif ion_type == "y":
            self.valid_mods = [mod for mod in self.modifications
                               if (mod.position >= len(self.full_sequence) - position and not mod.is_labile)]
        else:
            raise ValueError("Invalid ion type. Accept only a, b, y ions.")

    @property
    def raw_fragment_sequence(self):
        if self.ion_type == "b" or self.ion_type == "a":
            return self.full_sequence[:self.position]
        elif self.ion_type == "y":
            return self.full_sequence[-self.position:]
        return None

    @property
    def mz(self):
        mass = 0
        for aa in self.raw_fragment_sequence:
            mass += AA_MASSES[aa]
        if self.ion_type == "b":
            mass += 0
        elif self.ion_type == "a":
            mass -= (CARBON_MASS + OXYGEN_MASS)
        elif self.ion_type == "y":
            mass += WATER_MASS
        mass += sum(mod.delta for mod in self.valid_mods)
        mass += PROTON_MASS * self.charge
        mz = mass / self.charge
        return mz

    def __str__(self):
        self.valid_mods.sort(key=lambda x: x.position)
        if self.ion_type == "b" or self.ion_type == "a":
            raw = list(self.full_sequence[:self.position])
            for mod in reversed(self.valid_mods):
                raw.insert(mod.position + 1, f"[{mod.name}]")
            return "".join(raw)
        elif self.ion_type == "y":
            raw = list(self.full_sequence[-self.position:])
            for mod in reversed(self.valid_mods):
                raw.insert(mod.position - (len(self.full_sequence) - self.position) + 1, f"[{mod.name}]")
            return "".join(raw)
        return None

    __repr__ = __str__
