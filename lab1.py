from Bio import SeqIO


def get_valid_sequences(dna_sequence):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}

    valid_sequences = []

    i = 0
    while i < len(dna_sequence) - 2:
        if dna_sequence[i:i + 3] == start_codon:
            # Look for the next stop codon
            for j in range(i + 3, len(dna_sequence) - 2, 3):
                codon = dna_sequence[j:j + 3]
                if codon in stop_codons:
                    valid_sequences.append(dna_sequence[i:j + 3])
                    i = j + 3  # Move index to after the stop codon
                    break
            else:
                # If no valid stop codon was found, break out of the loop
                break
        else:
            i += 1

    return valid_sequences

file = "C:\\Users\\viliu\PycharmProjects\\bioinformatikaLab1\data\\test.fasta"
it = SeqIO.parse(file, "fasta")

for record in it:
    for sequence in get_valid_sequences(record.seq):
        print(sequence)