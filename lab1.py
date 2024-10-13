from Bio import SeqIO
from collections import Counter
import os
import math


# splice one of the 6 sequences into subsequences with start and stop codons
def get_valid_sequences(sequence):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}

    spliced_sequences = []

    i = 0
    while i < len(sequence) - 2:
        if sequence[i:i + 3] == start_codon:
            for j in range(i + 3, len(sequence) - 2, 3):
                codon = sequence[j:j + 3]
                if codon in stop_codons:
                    spliced_sequences.append(sequence[i:j + 3])
                    i = j + 3
                    break
            else:
                break
        else:
            i += 1
    # filter out sequences which are shorter than 100 bp
    valid_sequences = []
    for sequence in spliced_sequences:
        if len(sequence) >= 100:
            valid_sequences.append(sequence)

    return valid_sequences


def translate_sequence(sequences):
    protein_sequences = []
    for seq in sequences:
        protein_sequences.append(seq.translate())
    return protein_sequences


def get_codon_count(sequences):
    combined_codons = ""
    for seq in sequences:
        for s in seq:
            combined_codons += s

    return Counter(combined_codons)


def get_dicodon_count(sequences):
    dicodon_count = Counter()

    for seq in sequences:
        for s in seq:
            for i in range(len(s) - 1):
                dicodon = str(s[i:i + 2])
                dicodon_count[dicodon] += 1

    return dicodon_count


def get_protein_sequences(sequence):
    sequences = []
    protein_sequences = []

    reverse_sequence = sequence.reverse_complement()

    # divide sequence into 6 parts
    for i in range(3):
        sequences.append(sequence[i:])
    for i in range(3):
        sequences.append(reverse_sequence[i:])

    for seq in sequences:
        valid_sequences = get_valid_sequences(seq)
        protein_sequences.append(translate_sequence(valid_sequences))

    return protein_sequences


def get_frequency_from_counter(counter):
    freq = {}
    for val in counter:
        if not val.endswith('*'):
            freq[val] = counter[val] / sum(counter.values())
    return freq


def get_frequency_from_file(filename):
    path = os.path.abspath(os.curdir) + os.path.sep + "data" + os.path.sep + filename
    seq_frequencies = ()
    it = SeqIO.parse(path, "fasta")

    for record in it:
        protein_sequences = get_protein_sequences(record.seq)
        codon_counter = get_codon_count(protein_sequences)
        dicodon_counter = get_dicodon_count(protein_sequences)

        codon_freq = get_frequency_from_counter(codon_counter)
        dicodon_freq = get_frequency_from_counter(dicodon_counter)

        seq_frequencies = (codon_freq, dicodon_freq)

    return seq_frequencies


def calculate_key_distance(freq1, freq2):
    all_keys = set(freq1.keys()).union(set(freq2.keys()))
    distance = 0.0

    for key in all_keys:
        freq1 = freq1.get(key, 0.0)
        freq2 = freq2.get(key, 0.0)

        distance += (freq1 - freq2) ** 2

    return math.sqrt(distance)


def calculate_distance(sequence1_freq, sequence2_freq):
    codon_distance = calculate_key_distance(sequence1_freq[0], sequence2_freq[0])
    dicodon_distance = calculate_key_distance(sequence1_freq[1], sequence2_freq[1])

    return (codon_distance + dicodon_distance) / 2


bacterial1_freq = get_frequency_from_file("bacterial1.fasta")
bacterial2_freq = get_frequency_from_file("bacterial2.fasta")
print(calculate_distance(bacterial1_freq, bacterial2_freq))
# surasti daznius ✓
# pagal daznius apskaiciuoti atstuma pagal random kazkokia funkcija [ Distance(A, B) sqrt(sum(i=1, n)(Ai - Bi)^2) ] ✓
# atstumai yra tarp failu
