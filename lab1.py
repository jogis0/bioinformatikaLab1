from Bio import SeqIO
from collections import Counter

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


def get_codon_frequency(sequences):
    combined_codons = ""
    for seq in sequences:
        for s in seq:
            combined_codons += s

    return Counter(combined_codons)


def get_dicodon_frequency(sequences):
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


file = "C:\\Users\\viliu\\PycharmProjects\\bioinformatikaLab1\\data\\bacterial1.fasta"
it = SeqIO.parse(file, "fasta")

for record in it:
    protein_sequences = get_protein_sequences(record.seq)
    codon_counter = get_codon_frequency(protein_sequences)
    dicodon_counter = get_dicodon_frequency(protein_sequences)

    print("Codon counter\n")
    print(codon_counter)
    print("Dicodon counter\n")
    print(dicodon_counter)
    # surasti daznius
    # pagal daznius apskaiciuoti atstuma pagal random kazkokia funkcija
    # atstumai yra tarp failu
