# -------------------------------
# Abgabegruppe:
# Personen:
# HU-Accountname:
# -------------------------------
import time
from Bio import SeqIO


def build_suffix_array(dna_sequence):
    """
    Implements the Manber and Myers algorithm to efficiently construct a suffix array for dna_sequence. 
    
    Input:
        dna_sequence = str

    Output:
        result = List[int]
    """
    n = len(dna_sequence)
    k = 1  # Initial comparison length
    sa = list(range(n))  # Initial suffix array: [0, 1, 2, ..., n-1]

    # TODO: add your implementation here

    return sa


def find_all(dna_sequence, suffix_array, dna_pattern):
    """
    Implements the binary search algorithm to efficiently find all occurrences of dna_pattern within dna_sequence using its suffix array. 
    
    Input:
        dna_sequence = str,
        suffix_array = List[int],
        dna_pattern = str

    Output:
        result = List[int]
    """
    result = []

    # TODO: add your implementation here

    return result


def find_protein_coding_genes(dna_sequence, suffix_array):
    """
    Locates the position of protein coding genes within dna_sequence using its suffix array. 
    
    Input:
        dna_sequence = str,
        suffix_array = List[int]

    Output:
        genes = List[str],
        offsets = List[int]
    """
    # Gene Sequences
    genes = ["white (w)", "sex-lethal (sxl)"]
    offsets = []

    # TODO: add your implementation here

    return genes, offsets


def test_find_all(dna_sequence, suffix_array):
    """
    Tests the find_all function against a reference implementation to ensure correctness.

    Input:
        dna_sequence = str
        suffix_array = List[int]
    """
    # Reference function
    def _builtin_find_all(sequence, pattern):
        offsets = []
        index = sequence.find(pattern)
        while index != -1:
            offsets.append(index)
            index = sequence.find(pattern, index + 1)
        return offsets
    
    dna_pattern = "CCGCCCATCGCGATTCATCGACACAGAGGGCAAGGTGCGCAAACCGGAGTACTTCATACCC"
    assert sorted(find_all(dna_sequence, suffix_array, dna_pattern)) == sorted(_builtin_find_all(dna_sequence, dna_pattern))

    dna_pattern = "ACCGCTGACATAAAGTAGATCGCTGTTGTTATGf"
    assert len(find_all(dna_sequence, suffix_array, dna_pattern)) == 0

    dna_pattern = "GTCTC"
    assert sorted(find_all(dna_sequence, suffix_array, dna_pattern)) == sorted(_builtin_find_all(dna_sequence, dna_pattern))


if __name__ == "__main__":
    # Log runtimes for competition
    runtimes = []

    # FASTA file
    sequence_file = 'sequence.fasta'

    # Parse FASTA file
    records = list(SeqIO.parse(sequence_file, 'fasta'))
    assert len(records) == 1

    # Store DNA data
    dna_sequence = str(records[0].seq)

    # Build suffix array
    runtimes.append(time.perf_counter())
    suffix_array = build_suffix_array(dna_sequence)
    runtimes[-1] = time.perf_counter() - runtimes[-1]
    assert len(suffix_array) == len(dna_sequence)

    # Compute sanity checks
    test_find_all(dna_sequence, suffix_array)

    # Find protein coding genes
    runtimes.append(time.perf_counter())
    genes, offsets = find_protein_coding_genes(dna_sequence, suffix_array)
    runtimes[-1] = time.perf_counter() - runtimes[-1]

    for gene, offset in zip(genes, offsets):
        print(f"{gene}: {offset}")

    # Competition: Set flag based on your preference
    participate_in_competition = True

    if participate_in_competition:
        for idx, runtime in enumerate(runtimes):
            print(f"Task {2*idx+1}: {round(runtime, 2)} seconds")
    
