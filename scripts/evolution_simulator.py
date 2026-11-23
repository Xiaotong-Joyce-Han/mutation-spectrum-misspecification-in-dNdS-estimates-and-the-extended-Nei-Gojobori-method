# This script uses python to conduct evolution simulation.
import sys
import numpy as np
import pandas as pd
from collections import Counter
import os
import random

DIVERGENCE=0.01 # variable input
REPLICATE=100
# Folder structure:
# Root:
# ├── sequence/
# │   └── <index>.txt
# ├── codon_table/
# │   └── codon_table_1.csv
# |   └── codon_table_4.csv
# ├── Results/ # to store simulated sequence pairs
# │   └── <index>/
# │       └── <sele>/
# └── loading_data.csv

# each index corresponds to one mutation spectrum.
# they are indexed by csv row number in data/summary_of_mutation_spectrum_and_omega_ratio.csv
# for variable names, I usually use "species" for "index",
# although one species can have mutiple mutation spectra.

# loading_data.csv restores mutation spectrum parameters for each index.

# This script simulates 100 replicates of sequence evolution for each mutation spectrum.

ROOT_DIR="some path"
RESULTS_DIR="some path"

def count_base(sequence):
    counts = Counter(sequence)
    return [counts['A'], counts['T'], counts['C'], counts['G']]


def contain_unexpected_base(string):
    arr = np.array(list(string))
    valid_bases = np.array(['A', 'T', 'C', 'G'])
    return not np.all(np.isin(arr, valid_bases))


def precompute_base_positions(sequence, rank):
    positions = {'A': [], 'T': [], 'C': [], 'G': []}
    counts = {'A':0, 'T':0, 'C':0, 'G':0} 
    for i, base in enumerate(sequence):
        if base in positions:
            positions[base].append(i)
            rank.append(counts[base])
            counts[base] += 1
    return positions


def update_base_positions(positions, rank, base_type, new_base_type, index):
    if positions[base_type][-1] == index:
        positions[base_type].pop()
    else:
        positions[base_type][rank[index]] = positions[base_type].pop()
        rank[positions[base_type][rank[index]]] = rank[index]
    positions[new_base_type].append(index)
    rank[index] = len(positions[new_base_type])-1


def simulate_sequence(species, omega):
    
    root = ROOT_DIR
    sele = "nega" if omega=="0.3" else ("neu" if omega=="1" else "posi" )
    folder = RESULTS_DIR + species + "/" + sele + "/"

    sequence_path = root+"sequence/"+species+".txt"

    species = int(species)
    omega = float(omega)
    seed = int(species*1000 + omega*100)
    np.random.seed(seed)
    random.seed(seed)
    
    data_file_path = ROOT_DIR + "loading_data.csv"
    data = pd.read_csv(data_file_path)
    data = data[data['Row_index'] == species] 

    code = data.iloc[0, 8]
    r = data.iloc[0, 7]
    bias1, bias2, bias3, bias4, bias5, bias6 = data.iloc[0, 1:7]

    codon_table = pd.read_csv(ROOT_DIR + "codon_table/codon_table_" + str(code) + ".csv", header = None)

    propensity = [1, 1, r, r]
    bias = [[0, bias1, bias2, bias3],
            [bias1, 0, bias3, bias2],
            [bias5, bias6, 0, bias4],
            [bias6, bias5, bias4, 0]]

    codon_table_dict = {codon_table.at[i,0]: codon_table.at[i,1] for i in range(64)}

    def is_stop_codon(codon):
        return codon_table_dict[codon] == "stop"

    def compare_aa(codon, new_codon):
        return codon_table_dict[codon] == codon_table_dict[new_codon]

    genes = []
    current_gene = ""
    pseudo_indicator = False
    codons = []
    with open(sequence_path, 'r') as file:
        for l in file:
            if l.startswith(">"):
                if current_gene and (not pseudo_indicator) and len(current_gene)%3==0 and (not contain_unexpected_base(current_gene)):
                    genes.append(current_gene[:-3])
                current_gene = ""
                pseudo_indicator = "pseudo=true" in l
            else:
                current_gene += l.strip()
        if current_gene and (not pseudo_indicator) and len(current_gene)%3==0 and (not contain_unexpected_base(current_gene)):
            genes.append(current_gene[:-3])

    for gene in genes:
        for i in range(0, len(gene), 3):
            codon = gene[i: i+3]
            if not is_stop_codon(codon):
                codons.append(codon)

    original_genome = "".join(codons)
    seq_length = len(original_genome)
    prop_substitutions = DIVERGENCE
    num_of_substitutions = round(seq_length*prop_substitutions)

    for i in range(REPLICATE):

        sequence = list(original_genome)

        rank = []
        base_positions = precompute_base_positions(sequence, rank)
        base_counts = count_base(sequence)

        proposed_syns = 0
        proposed_nons = 0
        fixed_syns = 0
        fixed_nons = 0

        if omega<=1:
            fix_prob_syn = 1
            fix_prob_non = fix_prob_syn*omega
        else:
            fix_prob_syn = 0.33
            fix_prob_non = 0.99

        seq_index = str(i+1)

        while fixed_syns+fixed_nons < num_of_substitutions:
            probabilities = np.array(base_counts)*np.array(propensity)
            probabilities /= probabilities.sum()
            base_type_idx = np.random.choice(4, p=probabilities)
            base_type = np.array(['A', 'T', 'C', 'G'])[base_type_idx]

            if not base_positions[base_type]:
                continue

            index = random.choice(base_positions[base_type])
            new_base_type_idx = np.random.choice(4, p=bias[base_type_idx])
            new_base_type = np.array(['A', 'T', 'C', 'G'])[new_base_type_idx]

            codon = sequence[index-index%3:index-index%3+3]
            new_codon = codon[:index%3] + list(new_base_type) + codon[index%3+1:]
            if is_stop_codon(''.join(codon)) or is_stop_codon(''.join(new_codon)):
                continue

            if compare_aa(''.join(codon), ''.join(new_codon)):
                proposed_syns += 1
                if random.random() < fix_prob_syn:
                    sequence[index] = new_base_type
                    base_counts[base_type_idx] -= 1
                    base_counts[new_base_type_idx] += 1
                    update_base_positions(base_positions, rank, base_type, new_base_type, index)
                    fixed_syns += 1

            else:
                proposed_nons += 1
                if random.random() < fix_prob_non:
                    sequence[index] = new_base_type
                    base_counts[base_type_idx] -= 1
                    base_counts[new_base_type_idx] += 1
                    update_base_positions(base_positions, rank, base_type, new_base_type, index)
                    fixed_nons += 1                   

        sequence = ''.join(sequence)

        mainfile = folder + seq_index + ".txt"
        nsfile = folder + seq_index + "ns.txt"
        NSfile = folder + seq_index + "NS.txt"

        with open(mainfile, 'w') as f:
            f.write(f"2 {seq_length}\n")
            f.write("Initial\n")
            f.write(original_genome+"\n")
            f.write("Evolved\n")
            f.write(sequence)

        with open(nsfile, "w") as f:
            f.write(f"{fixed_nons}\t{fixed_syns}\t{fixed_nons/fixed_syns}")
        
        with open(NSfile, "w") as f:
            f.write(f"{proposed_nons}\t{proposed_syns}\t{proposed_nons/proposed_syns}")

if __name__ == "__main__":

    species = sys.argv[1]
    omega = sys.argv[2]

    simulate_sequence(species, omega)
