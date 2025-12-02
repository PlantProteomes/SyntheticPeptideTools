# Parses sequence (ex. AQDSQVLEEER[Label:13C(6)15N(4)]) as string. Can return list of amino acids, list of modifications (unbracketed), and list of unmodified amino acids
# Some assistance given by ChatGPT
import timeit

def parse_sequence(modified_sequence): # converts everything in the sequence to a list of modified amino acids
    bracket_depth = 0
    bracketed_sequence = ""
    all_units = []
    first_mod = ""
    for chara in modified_sequence:
        if chara.isupper() and bracket_depth == 0:
            if first_mod != "":
                all_units.append(first_mod + chara)
                first_mod = ""
            else:
                all_units.append(chara)
        if chara == "[":
            bracket_depth += 1
        if bracket_depth != 0:
            bracketed_sequence += chara
        if chara == "]":
            bracket_depth -= 1
        if bracket_depth == 0 and bracketed_sequence != "":
            if not all_units:
                first_mod = bracketed_sequence
            else:
                all_units[-1] += bracketed_sequence
            bracketed_sequence = ""

    return all_units

def get_mods(modified_sequence): # extracts only mods from sequence
    bracket_depth = 0
    bracketed_sequence = ""
    all_mods = []
    for chara in modified_sequence:
        if chara == "[":
            bracket_depth += 1
        if bracket_depth != 0:
            bracketed_sequence += chara
        if chara == "]":
            bracket_depth -= 1
        if bracket_depth == 0 and bracketed_sequence != "":
            all_mods.append(bracketed_sequence[1:-1])
            bracketed_sequence = ""

    return all_mods

def get_acids(modified_sequence): # extracts only unmodified amino acids
    bracket_depth = 0
    bracketed_sequence = ""
    unmodified_sequence = ""
    for chara in modified_sequence:
        if chara == "[":
            bracket_depth += 1
        if bracket_depth != 0:
            bracketed_sequence += chara
        if chara == "]":
            bracket_depth -= 1
        if bracketed_sequence == "":
            unmodified_sequence += chara
        if bracket_depth == 0 and bracketed_sequence != "":
            bracketed_sequence = ""

    return unmodified_sequence


# if __name__ == "__main__":
# #     print(timeit.timeit(stmt="check_missing('251103_mEclipse_ncORF89-S1', 4343, 'AQDSQVLEEER[Label:13C(6)15N(4)', 'missing QV', 999007, 0.002)",setup="from __main__ import check_missing", number=10))
#     print(timeit.timeit(stmt="parse_sequence('AQDSQVLEEER[Label:13C(6)15N(4)]')",setup="from __main__ import parse_sequence", number=10))
#     # print(timeit.timeit(stmt=snippet, globals=globals(), number=10))
