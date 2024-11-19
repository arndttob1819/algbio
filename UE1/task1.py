# -------------------------------
# Abgabegruppe:
# Personen:
# HU-Accountname:
# -------------------------------
import time
from Bio import SeqIO
import matplotlib.pyplot as plt


def find_all(dna_sequence, dna_pattern):
    """
    Implements the Knuth-Morris-Pratt (KMP) algorithm to efficiently find all occurrences of dna_pattern within dna_sequence. 
    
    Input:
        dna_sequence = str,
        dna_pattern = str

    Output:
        result = List[int]
    """
    result = []

    # TODO: add your implementation

    return result


def count_restriction_enzyme_sites(dna_sequnece):
    """
    Counts recognition sites for specified restriction enzymes within dna_sequence.

    Input:
        dna_sequence = str

    Output:
        enzymes = List[str],
        sites = List[str],
        counts = List[int]
    """
    enzymes = ["EcoRI", "BamHI", "HindIII", "NotI", "I-SceI"]
    sites = ["GAATTC", "GGATCC", "AAGCTT", "GCGGCCGC", "TAGGGATAACAGGGTAAT"]
    counts = []

    # TODO: add your implementation

    return enzymes, sites, counts


def compute_NotI_site_distribution(dna_sequence):
    """
    Plots the distribution of NotI restriction enzyme recognition sites within dna_sequence.

    Input:
        dna_sequence = str

    Output:
        ax = matplotlib.axes.Axes
    """
    enzyme = "NotI"
    site = "GCGGCCGC"
    _, ax = plt.subplots()
    
    # TODO: add your implementation
    
    return ax


def find_longest_microsatellite_repeat(dna_sequence):
    """
    Finds the longest consecutive repeat of a specified trinucleotide microsatellite in dna_sequence.

    Input:
        dna_sequence = str

    Output:
        trinucleotide = str,
        longest_repeat = int
    """
    trinucleotide = "AAT"
    longest_repeat = 0
    
    # TODO: add your implementation

    return trinucleotide, longest_repeat


def count_palindromic_sequences(dna_sequence):
    """
    Finds palindromic sequences of a specified length in dna_sequence and plots their counts.

    Input:
        dna_sequence = str

    Output:
        ax = matplotlib.axes.Axes
    """
    length = 4
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    _, ax = plt.subplots()

    # TODO: add your implementation

    return ax


def test_find_all(dna_sequence):
    """
    Tests the find_all function against a reference implementation to ensure correctness.

    Input:
        dna_sequence = str
    """
    # Reference function
    def _builtin_find_all(sequence, pattern):
        offsets = []
        index = sequence.find(pattern)
        while index != -1:
            offsets.append(index)
            index = sequence.find(pattern, index + 1)
        return offsets
    
    dna_pattern = "AGGACTCAGTCTGTCAGATACTTAGGACTCGACATGCATAAAGGAGAA"
    assert sorted(find_all(dna_sequence, dna_pattern)) == sorted(_builtin_find_all(dna_sequence, dna_pattern))

    dna_pattern = "AGGACTCAGTCTGTCAGATACTTAGGACTCGACATGCATAAAGGAGAAf"
    assert len(find_all(dna_sequence, dna_pattern)) == 0

    dna_pattern = "TCCCAGC"
    assert sorted(find_all(dna_sequence, dna_pattern)) == sorted(_builtin_find_all(dna_sequence, dna_pattern))


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

    # Compute sanity checks
    test_find_all(dna_sequence)

    # Restriction enzyme sites
    '''
    EcoRI, GAATTC: 52432
    BamHI, GGATCC: 36449
    HindIII, AAGCTT: 54061
    NotI, GCGGCCGC: 209
    I-SceI, TAGGGATAACAGGGTAAT: 0
    '''
    runtimes.append(time.process_time())
    enzymes, sites, counts = count_restriction_enzyme_sites(dna_sequence)
    runtimes[-1] = time.process_time() - runtimes[-1]

    # print("Restriction enzyme sites")
    for enzyme, site, count in zip(enzymes, sites, counts):
        print(f"{enzyme}, {site}: {count}")

    # NotI distribution
    runtimes.append(time.process_time())
    ax = compute_NotI_site_distribution(dna_sequence)
    runtimes[-1] = time.process_time() - runtimes[-1]

    plt.savefig("NotI_distribution.pdf", bbox_inches="tight")

    # Microsatellite Repeats
    runtimes.append(time.process_time())
    microsatellite, longest_repeat = find_longest_microsatellite_repeat(dna_sequence)
    runtimes[-1] = time.process_time() - runtimes[-1]
    
    # Longest AAT repeat: 19
    print(f"Longest {microsatellite} repeat: {longest_repeat}")

    # Palindromic sequences
    runtimes.append(time.process_time())
    ax = count_palindromic_sequences(dna_sequence)
    runtimes[-1] = time.process_time() - runtimes[-1]

    plt.savefig("palindromic_sequences_counts.pdf", bbox_inches="tight")

    # Competition: Set flag based on your preference
    participate_in_competition = True

    if participate_in_competition:
        for idx, runtime in enumerate(runtimes):
            print(f"Task {idx+2}: {round(runtime, 2)} seconds")