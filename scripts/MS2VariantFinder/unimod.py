import fastobo
from urllib.request import urlopen
from constants import ppm

url = "https://proteomecentral.proteomexchange.org/extern/CVs/unimod.obo"
_unimod_list = None
_number_index = None
_name_index = None

def _load():
    global _unimod_list, _number_index, _name_index
    if _unimod_list is not None:
        return

    _unimod = fastobo.load(urlopen(url))
    _unimod_list = []
    _number_index = {}
    _name_index = {}

    for i in range(1, len(_unimod)):
        temp_dict = {"locales": []}
        for clause in _unimod[i]:
            if isinstance(clause, fastobo.term.NameClause):
                term_name = clause.name
                temp_dict["name"] = term_name
            elif isinstance(clause, fastobo.term.XrefClause):
                xref = clause.xref
                if "site" in str(xref.id):
                    temp_dict["locales"].append(xref.desc)
                temp_dict[str(xref.id)] = xref.desc
        _unimod_list.append(temp_dict)

    for entry in _unimod_list:
        indexed_mass = int(float(entry["delta_mono_mass"]))
        _number_index.setdefault(indexed_mass, []).append(entry)
        first_chara = entry["name"][0]
        _name_index.setdefault(first_chara, []).append(entry)

def get_mod(name):
    _load()
    if name[0] not in _name_index:
        raise KeyError(f"Name {name} not found in Unimod.")
    for i in _name_index[name[0]]:
        if i["name"] == name:
            return i
    return None

# list of all mods that have masses within a specific mass delta tolerance
def get_candidate_mods(mass_delta, tolerance, precursor_mz):
    _load()
    candidates = []
    if int(mass_delta) not in _number_index:
        return None
    for i in _number_index[int(mass_delta)]:
        if abs(float(i["delta_mono_mass"]) - mass_delta) < ppm((precursor_mz + mass_delta), tolerance):
            candidates.append(i)
    if len(candidates) > 0:
        return candidates
    else:
        return None

def is_approved(name, locale):
    _load()
    if name[0] not in _name_index:
        raise KeyError(f"Name {name} not found in Unimod.")
    for i in _name_index[name[0]]:
        if i["name"] == name:
            if locale in i["locales"]:
                return True
    return False
