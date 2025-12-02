import fastobo
from urllib.request import urlopen
import math
import timeit

url = "https://proteomecentral.proteomexchange.org/extern/CVs/unimod.obo"
unimod = fastobo.load(urlopen(url))

# TODO: Add locations of modifications

unimod_dict = []

for i in range(1, len(unimod)):
    temp_dict = {}
    for clause in unimod[i]:
        if isinstance(clause, fastobo.term.NameClause):
            term_name = clause.name
            temp_dict["name"] = term_name
        elif isinstance(clause, fastobo.term.XrefClause):
            xref = clause.xref
            temp_dict[str(xref.id)] = xref.desc
            # if str(xref.id) == "delta_mono_mass":
            #     delta_mono_mass = float(xref.desc)
            #     temp_dict["delta mono mass"] = delta_mono_mass
    unimod_dict.append(temp_dict)


number_indexed_unimod_dict = {}
for entry in unimod_dict:
    indexed_mass = math.floor(float(entry["delta_mono_mass"]))
    if indexed_mass not in number_indexed_unimod_dict.keys():
        number_indexed_unimod_dict[indexed_mass] = [entry]
    else:
        number_indexed_unimod_dict[indexed_mass].append(entry)

def search_unimod_by_mass(mass_delta, tolerance):
    if math.floor(mass_delta) not in number_indexed_unimod_dict:
        return ""
    for i in number_indexed_unimod_dict[math.floor(mass_delta)]:
        if abs(float(i["delta_mono_mass"])-mass_delta) < tolerance:
            return i["name"]

    # for i in range(1, len(unimod)):
    #     term_name = ""
    #     for clause in unimod[i]:
    #         if isinstance(clause, fastobo.term.NameClause):
    #             term_name = clause.name
    #         elif isinstance(clause, fastobo.term.XrefClause):
    #             xref = clause.xref
    #             if str(xref.id) == "delta_mono_mass":
    #                 if abs(float(xref.desc) - mass_delta) <= tolerance:
    #                     return term_name
    return ""

name_indexed_unimod_dict = {}
for entry in unimod_dict:
    first_chara = entry["name"][0]
    if first_chara not in name_indexed_unimod_dict.keys():
        name_indexed_unimod_dict[first_chara] = [entry]
    else:
        name_indexed_unimod_dict[first_chara].append(entry)

def search_unimod_by_name(name):
    if name[0] not in name_indexed_unimod_dict:
        return ""
    for i in name_indexed_unimod_dict[name[0]]:
        if i["name"] == name:
            return float(i["delta_mono_mass"])

    # for i in range(1, len(unimod)):
    #     term_name = ""
    #     for clause in unimod[i]:
    #         if isinstance(clause, fastobo.term.NameClause):
    #             if str(clause.name) == name:
    #                 term_name = name
    #         elif isinstance(clause, fastobo.term.XrefClause):
    #             xref = clause.xref
    #             if str(xref.id) == "delta_mono_mass":
    #                 if term_name == name:
    #                     return xref.desc
    return ""

# if __name__ == "__main__":
# #     print(timeit.timeit(stmt="check_missing('251103_mEclipse_ncORF89-S1', 4343, 'AQDSQVLEEER[Label:13C(6)15N(4)', 'missing QV', 999007, 0.002)",setup="from __main__ import check_missing", number=10))
#     print(timeit.timeit(stmt="search_unimod_by_name('Deamidation')",setup="from __main__ import search_unimod_by_name", number=10))

def localize(mod_name):
    if mod_name[0] not in name_indexed_unimod_dict:
        return
    for i in name_indexed_unimod_dict[mod_name[0]]:
        if i["name"] == mod_name:
            locations = []
            for entry in i.keys():
                locations.append(i[entry]) if "site" in entry else None
            return locations
    return
