from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
import numpy as np

def calculate_pairwise_identity(alignment_file):
    alignment = AlignIO.read(alignment_file, "fasta")
    num_seqs = len(alignment)
    length = alignment.get_alignment_length()
    
    # Initialize a matrix to store identity values
    identity_matrix = np.zeros((num_seqs, num_seqs))

    # Loop over each pair of sequences and calculate identity
    for i in range(num_seqs):
        for j in range(i, num_seqs):
            seq1 = alignment[i].seq
            seq2 = alignment[j].seq
            identity_count = sum(base1 == base2 for base1, base2 in zip(seq1, seq2))
            identity = identity_count / length
            identity_matrix[i][j] = identity
            identity_matrix[j][i] = identity
    
    return identity_matrix, [record.id for record in alignment]

def save_identity_matrix(identity_matrix, seq_ids, output_file):
    with open(output_file, 'w') as f:
        f.write("ID\t" + "\t".join(seq_ids) + "\n")
        for i, row in enumerate(identity_matrix):
            f.write(seq_ids[i] + "\t" + "\t".join(map(str, row)) + "\n")

# Usage
alignment_file = "aligned_sequences.fasta"
identity_matrix, seq_ids = calculate_pairwise_identity(alignment_file)
save_identity_matrix(identity_matrix, seq_ids, "pairwise_identity_matrix.txt")
