import fastobo
from urllib.request import urlopen

url = "https://proteomecentral.proteomexchange.org/extern/CVs/unimod.obo"
unimod = fastobo.load(urlopen(url))

# TODO: Implement dict-based search based on Mia's code

def search_unimod_by_mass(mass_delta, tolerance):
    for i in range(1, len(unimod)):
        term_name = ""
        for clause in unimod[i]:
            if isinstance(clause, fastobo.term.NameClause):
                term_name = clause.name
            elif isinstance(clause, fastobo.term.XrefClause):
                xref = clause.xref
                if str(xref.id) == "delta_mono_mass":
                    if abs(float(xref.desc) - mass_delta) <= tolerance:
                        return term_name
    return ""

def search_unimod_by_name(name):
    for i in range(1, len(unimod)):
        term_name = ""
        for clause in unimod[i]:
            if isinstance(clause, fastobo.term.NameClause):
                if str(clause.name) == name:
                    term_name = name
            elif isinstance(clause, fastobo.term.XrefClause):
                xref = clause.xref
                if str(xref.id) == "delta_mono_mass":
                    if term_name == name:
                        return xref.desc
    return ""
