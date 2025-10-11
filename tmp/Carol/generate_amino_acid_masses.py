import pandas as pd

acids = {
	"A" : 71.0371,
	"R" : 156.1011,
	"N" : 114.0429,
	"D" : 115.0269,
	"C" : 103.0477,
	"Q" : 128.0586,
	"E" : 129.0426,
	"G" : 57.0215,
	"H" : 137.0589,
	"I" : 113.0841,
	"L" : 113.0841,
	"K" : 128.0950,
	"M" : 131.0405,
	"F" : 147.0684,
	"P" : 97.0528,
	"S" : 87.0320,
	"T" : 101.0477,
	"W" : 186.0793,
	"Y" : 163.0633,
	"V" : 99.0684
}

modifications = {
    "267" : 10.008269
}

def print_pairs(sequence):
    sequence_list = []
    for letter in sequence:
        sequence_list.append(letter)
    print(f"Pairs in sequence {sequence}:")
    pairs = set()

    for i in range(len(sequence_list)):
        for j in range(i+1, len(sequence_list)):
            pair = tuple(sorted([sequence_list[i], sequence_list[j]]))
            pairs.add(pair)
    pairs_list = sorted(list(pairs))

    for pair in pairs_list:
        mass = acids[pair[0]] + acids[pair[1]]
        if pair[0] == "R":
            mass += modifications["267"]
        if pair[1] == "R":
            mass += modifications["267"]
        print(pair, f"{mass:0,.4f}")

def print_triples(sequence):
    sequence_list = []
    for letter in sequence:
        sequence_list.append(letter)
    print(f"Triples in sequence {sequence}:")
    triples = set()

    for i in range(len(sequence_list)):
        for j in range(i + 1, len(sequence_list)):
            for k in range (j + 1, len(sequence_list)):
                triple = tuple(sorted([sequence_list[i], sequence_list[j], sequence_list[k]]))
                triples.add(triple)
    triples_list = sorted(list(triples))

    for triple in triples_list:
        mass = acids[triple[0]] + acids[triple[1]] + acids[triple[2]]
        if triple[0] == "R":
            mass += modifications["267"]
        if triple[1] == "R":
            mass += modifications["267"]
        if triple[2] == "R":
            mass += modifications["267"]
        print(triple, f"{mass:0,.4f}")

print_pairs("AQDSQVLEEER")
print_triples("AQDSQVLEEER")