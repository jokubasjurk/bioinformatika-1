import numpy as np
from Bio import SeqIO
from Bio.Data.CodonTable import unambiguous_dna_by_name

files = ["bacterial1.fasta", "bacterial2.fasta", "bacterial3.fasta", "bacterial4.fasta", "mamalian1.fasta", "mamalian2.fasta", "mamalian3.fasta", "mamalian4.fasta"]
start_codon = 'ATG'
stop_codons = ['TAA', 'TAG', 'TGA']
alphabet_list = []

def find_orf(seq_record):
    i = 0
    codon_list = []
    while i < len(seq_record):
        if seq_record[i] == start_codon:
            start_pos = i
            j = i
            while j < len(seq_record):
                if seq_record[j] in stop_codons:
                    end_pos = j
                    codon_list.append(''.join(str(e) for e in seq_record[start_pos:end_pos + 1]))
                    i = j
                    break
                j += 1
        i += 1
    return codon_list

def find_orfs(triplets):
    orfs = []
    for triplet in triplets:
        orfs.append(find_orf(triplet))
    return orfs

def find_farthest(seq_record):
    i = 0
    codon_list = []
    while i < len(seq_record):
        if seq_record[i] in stop_codons:
            start_pos = i
            j = i + 1
            end_pos = -1
            while j < len(seq_record):
                if seq_record[j] == start_codon:
                    end_pos = j
                elif seq_record[j] in stop_codons:
                    if end_pos != -1:
                        codon_list.append(''.join(str(e) for e in seq_record[start_pos:end_pos + 1]))
                    i = j
                    break
                j += 1
        i += 1
    return codon_list


def find_farthest_starts(triplets):
    starts = []
    for triplet in triplets:
        starts.append(find_farthest(triplet))
    return starts


def find_longest_fragments(orfs):
    longest_orfs = []
    for orf in orfs:
        if len(orf) > 100:
            longest_orfs.append(orf)
    return longest_orfs


def find_codon_frequency(list):
    standard_table = unambiguous_dna_by_name["Bacterial"]
    letters = standard_table.nucleotide_alphabet
    alphabet_list = []
    for l1 in letters:
        for l2 in letters:
            for l3 in letters:
                codon = l1 + l2 + l3
                alphabet_list.append(codon)
    codon_alphabet = dict.fromkeys(alphabet_list, 0.0)

    total = 0
    for orf in list:
        codons_in_orf = [orf[i:i + 3] for i in range(0, len(orf), 3)]
        total += len(codons_in_orf)
        for codon in codons_in_orf:
            codon_alphabet[codon] += 1
    for key in codon_alphabet:
        codon_alphabet[key] = (codon_alphabet[key] / total) * 100
    return codon_alphabet


def find_dicodon_frequency(list):
    standard_table = unambiguous_dna_by_name["Bacterial"]
    letters = standard_table.nucleotide_alphabet
    alphabet_list = []
    # get all dicodon variations
    for l1 in letters:
        for l2 in letters:
            for l3 in letters:
                for l4 in letters:
                    for l5 in letters:
                        for l6 in letters:
                            dicodon = l1 + l2 + l3 + l4 + l5 + l6
                            alphabet_list.append(dicodon)
    dicodon_alphabet = dict.fromkeys(alphabet_list, 0.0)

    total = 0
    for orf in list:
        dicodons_in_orf = [orf[i:i + 6] for i in range(0, len(orf), 6)]
        total += len(dicodons_in_orf)
        for dicodon in dicodons_in_orf:
            if len(dicodon) == 6:
                dicodon_alphabet[dicodon] += 1
    for key in dicodon_alphabet:
        dicodon_alphabet[key] = (dicodon_alphabet[key] / total) * 100
    return dicodon_alphabet

def print_distance_matrix(frequencies):
    ids = []

    for id in frequencies:
        ids.append(id)

    data = frequencies[ids[0]]

    print("*******", end = ' ')
    for i in ids:
        print(i, end = ' ')

    for i in range(0,8):
        print("")
        print(ids[i], end = ' ')
        for j in range(0, 8):
            total_points = 0
            for k in data:
                points = frequencies[ids[i]][k] - frequencies[ids[j]][k]
                if points < 0:
                    points = points * (-1)
                total_points += points
            print("%.2f" % total_points, end = ' ')

if __name__ == '__main__':
    codon_frequencies = {}
    dicodon_frequencies = {}

    for i in files:
        record = SeqIO.read("data/" + i, "fasta")
        triplets = []
        triplets.append([record.seq[i:i + 3] for i in range(0, len(record.seq), 3)])
        triplets.append([record.seq[i:i + 3] for i in range(1, len(record.seq), 3)])
        triplets.append([record.seq[i:i + 3] for i in range(2, len(record.seq), 3)])
        triplets.append([record.seq.reverse_complement()[i:i + 3] for i in range(0, len(record.seq), 3)])
        triplets.append([record.seq.reverse_complement()[i:i + 3] for i in range(1, len(record.seq), 3)])
        triplets.append([record.seq.reverse_complement()[i:i + 3] for i in range(2, len(record.seq), 3)])
        # 1.
        codons = find_orfs(triplets)
        codons = np.concatenate(codons)
        # 2.
        farthest_start_codons = find_farthest_starts(triplets)
        # 3.
        longest_fragments = find_longest_fragments(codons)
        # 4.
        codon_frequencies[record.id] = find_codon_frequency(longest_fragments)
        dicodon_frequencies[record.id] = find_dicodon_frequency(longest_fragments)

    print("Codon frequency matrix:")
    print_distance_matrix(codon_frequencies)
    print("")
    print("Dicodon frequency matrix:")
    print_distance_matrix(dicodon_frequencies)