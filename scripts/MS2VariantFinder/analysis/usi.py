from models import Peptide, Modification, Scan
import unimod
from constants import ppm
import numpy as np

def score_ions(ions, sorted_mz_array, sorted_intensity_array, total_intensity, tolerance):
    intensity_sum = 0
    match_count = 0
    for ion in ions:
        left = np.searchsorted(sorted_mz_array, ion.mz - ppm(ion.mz, tolerance), side="left")
        right = np.searchsorted(sorted_mz_array, ion.mz + ppm(ion.mz, tolerance), side="right")
        intensity_slice = sorted_intensity_array[left:right]
        matched_peak_intensity = np.sum(intensity_slice)
        match_count += 1 if matched_peak_intensity > 0 else 0
        intensity_sum += matched_peak_intensity

    score = 5 * intensity_sum / total_intensity + 5 * match_count / len(ions) # prioritizes intensity and number of matches
    return score

def localize(sequence: Peptide, mod_name: str, tolerance, spectrum: Scan):
    # Returns likely location of specific modification.
    sorted_mz_array = np.sort(spectrum.mz_array)
    sorted_intensity_array = spectrum.intensity_array[np.argsort(spectrum.mz_array)]
    total_intensity = np.sum(sorted_intensity_array)
    code_string_list = []
    best_score = 0.0
    best_sequence = ""
    testing_sequence = Peptide(sequence.raw_sequence, [Modification(m.position, m.delta, m.name, m.is_labile) for m in sequence.modifications])

    modification = Modification.from_string(-1, mod_name, False)
    testing_sequence.modifications.append(modification)

    for i in range(-2, len(testing_sequence)):
        code_string = ""
        if i == -2:
            modification.is_labile = True
        else:
            modification.is_labile = False
        modification.position = i if i > -2 else -1
        fragments = testing_sequence.fragments()
        all_ions = [ion for ion_list in fragments.values() for ion in ion_list]
        if modification.is_labile:
            code_string += "Labile-"
        else:
            is_approved = unimod.is_approved(mod_name, testing_sequence.get_residue(modification.position))
            code_string += (testing_sequence.get_residue(modification.position).upper() if is_approved else testing_sequence.get_residue(modification.position).lower()) + "-"

        score = score_ions(all_ions, sorted_mz_array, sorted_intensity_array, total_intensity, tolerance)
        code_string += f"{score:.2f}"
        code_string_list.append(code_string)

        if score > best_score:
            best_score = score
            best_sequence = str(testing_sequence)

    return best_sequence, best_score, ", ".join(code_string_list)

def synthesis_error(sequence: Peptide, mass_delta, tolerance):
    # Returns likely candidates for synthesis errors along with the window of masses that fall within desired tolerance.
    if mass_delta < 0:
        is_negative = True
        mass_delta *= -1
    else:
        is_negative = False
    for r in range(1, 5):
        if r > 1 and not is_negative:
            return None
        combinations = sequence.generate_error_masses(r, deduplicate=not is_negative)
        items = sorted(combinations.items(), key=lambda x: x[1])
        combination_list = [key for key, value in items]
        mass_list = np.array([value for key, value in items])
        left = np.searchsorted(mass_list, mass_delta - ppm(mass_delta, tolerance), side="left")
        right = np.searchsorted(mass_list, mass_delta + ppm(mass_delta, tolerance), side="right")
        mass_slice = mass_list[left:right]
        if len(mass_slice) != 0:
            errors = {}
            for i in range(len(mass_slice)):
                errors[combination_list[left + i]] = mass_slice[i]
            return errors
    return None

def localize_synthesis_error(sequence: Peptide, errors, mass_delta, tolerance, spectrum: Scan):
    sorted_mz_array = np.sort(spectrum.mz_array)
    sorted_intensity_array = spectrum.intensity_array[np.argsort(spectrum.mz_array)]
    total_intensity = np.sum(sorted_intensity_array)
    code_string_list = []
    best_score = 0.0
    best_sequence = ""
    best_error_str = ""
    if errors is None or len(errors) == 0:
        return None
    if mass_delta < 0:
        for error in errors:
            sequence_with_blanks = list(sequence.raw_sequence)
            for i in error:
                sequence_with_blanks[i] = "_"
            code_string = "".join(sequence_with_blanks) + "-"
            error_str = "".join(sequence.raw_sequence[i] for i in error)
            test_peptide = sequence.remove_residues(error)
            fragments = test_peptide.fragments()
            all_ions = [ion for ion_list in fragments.values() for ion in ion_list]

            score = score_ions(all_ions, sorted_mz_array, sorted_intensity_array, total_intensity, tolerance)
            code_string += f"{score:.2f}"
            code_string_list.append(code_string)

            if score > best_score:
                best_score = score
                best_sequence = str(test_peptide)
                best_error_str = "missing " + error_str
    elif mass_delta > 0:
        seen = set()
        for i in range(len(sequence)+1):
            for error in errors:
                residue_char = sequence.get_residue(error[0])
                mods_at_pos = frozenset((m.name, m.is_labile) for m in sequence.modifications if m.position == error[0])
                signature = (i, residue_char, mods_at_pos)
                if signature in seen:
                    continue
                seen.add(signature)
                test_peptide = sequence.insert_residues(error[0], i)
                fragments = test_peptide.fragments()
                all_ions = [ion for ion_list in fragments.values() for ion in ion_list]
                code_string = str(test_peptide) + "-"

                score = score_ions(all_ions, sorted_mz_array, sorted_intensity_array, total_intensity, tolerance)
                code_string += f"{score:.2f}"
                code_string_list.append(code_string)

                if score > best_score:
                    best_score = score
                    best_sequence = str(test_peptide)
                    best_error_str = "extra " + residue_char
    return best_sequence, best_score, ", ".join(code_string_list), best_error_str

def generate_usi(spectrum: Scan, sequence: Peptide, mass_delta, tolerance):
    final_candidates = []

    if abs(mass_delta) <= ppm(sequence.mz, tolerance):
        fragments = sequence.fragments()
        sorted_mz_array = np.sort(spectrum.mz_array)
        sorted_intensity_array = spectrum.intensity_array[np.argsort(spectrum.mz_array)]
        total_intensity = np.sum(sorted_intensity_array)
        all_ions = [ion for ion_list in fragments.values() for ion in ion_list]
        no_mod_score = score_ions(all_ions, sorted_mz_array, sorted_intensity_array, total_intensity, tolerance)
        final_candidates.append((no_mod_score, str(sequence), None, "No mod", 0.0, ""))

    candidate_mods = unimod.get_candidate_mods(mass_delta, tolerance, spectrum.precursor_mz)
    if candidate_mods is not None:
        # tiebreaking for mods with identical mass (score will not show any difference)
        closest_mod_mass_diff = min(abs(x["delta_mono_mass"] - mass_delta) for x in candidate_mods)
        tied_mods = [mod for mod in candidate_mods if
                     abs(mod["delta_mono_mass"] - spectrum.precursor_mz) == closest_mod_mass_diff]
        approved_mods = []
        for mod in tied_mods:
            for locale in sequence.raw_sequence:
                if unimod.is_approved(mod["name"], locale):
                    approved_mods.append(mod)  # check if unimod is ok
                    break
            tied_mods = [m for m in tied_mods if "->" not in m["name"]]
        if len(approved_mods) > 0:  # if at least some of them are approved, return the first of these; if none are, keep all of them and return the first.
            best_mod = approved_mods[0]
        else:
            best_mod = tied_mods[0]
        best_mod_sequence, best_mod_score, mod_code_string = localize(sequence, best_mod["name"], tolerance, spectrum)
        mod_type = "cation" if "Cation" in best_mod["name"] else ""
        final_candidates.append((best_mod_score, str(best_mod_sequence), mod_code_string, best_mod, best_mod["delta_mono_mass"], mod_type))

    # tiebreaking for synthesis errors with identical mass (score will show difference)
    candidate_synthesis_errors = synthesis_error(sequence, mass_delta, tolerance)
    if candidate_synthesis_errors is not None:
        closest_error_mass_diff = min(candidate_synthesis_errors.values(), key=lambda x: abs(x - mass_delta))
        tied_errors = [error for error in candidate_synthesis_errors.items() if error[1] == closest_error_mass_diff]
        best_error_sequence, best_error_score, error_code_string, best_error_str = localize_synthesis_error(sequence, tied_errors, mass_delta, tolerance, spectrum)
        final_candidates.append((best_error_score, str(best_error_sequence), error_code_string, best_error_str, closest_error_mass_diff, "synthesis error"))

    # comparing mods and synthesis errors if all plausible.
    if len(final_candidates) == 0:
        return None, None, None, None, None
    return max(final_candidates, key=lambda x: x[0])[1:]

