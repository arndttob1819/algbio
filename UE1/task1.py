# -------------------------------
# Abgabegruppe: 19
# Personen: Tobias Arndt, Benjamin Höhnisch, Tim Patzak
# HU-Accountname: arndttob, hoehnisb, patzakti
# -------------------------------
from Bio import SeqIO
import matplotlib.pyplot as plt
import time


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

    # Setup look-up-table (LPS)
    lookUp = [0]
    for i in range(1, len(dna_pattern)):
        # Check if smaller suffix was already found
        if lookUp[i-1] != 0:
            # Check if suffix increased
            if dna_pattern[lookUp[i-1]] == dna_pattern[i]:
                lookUp.append(lookUp[i-1] + 1)

                # LPS'
                lookUp[i-1] = 0
            # Check if suffix is length 1
            elif dna_pattern[0] == dna_pattern[i]:
                lookUp.append(1)
            else:
                lookUp.append(0)
        else:
            # Check if current char is a prefix
            if dna_pattern[0] == dna_pattern[i]:
                lookUp.append(1)
            else:
                lookUp.append(0)

    # KMP
    i = 0
    j = 0
    while i + len(dna_pattern) <= len(dna_sequence):
        if dna_sequence[i] == dna_pattern[j]:
            # If chars match move both pointers
            i += 1
            j += 1
        else:
            if j == 0:
                # If first pattern char doesn't match move only sequence pointer
                i += 1
            else:
                # Jump
                j = lookUp[j-1]

        # Pattern was found
        if j == len(dna_pattern):
            result.append(i-len(dna_pattern))
            j = lookUp[-1]

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

    # Search for patterns individually and count occurrences
    for pattern in sites:
        offsets = find_all(dna_sequence, pattern)
        counts.append(len(offsets))

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
    
    # Find all occurrences of the enzyme
    offsets = find_all(dna_sequence, site)

    # Add data to plot
    plt.hist(offsets, bins=100)
    
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

    offsets = find_all(dna_sequence, trinucleotide)

    # Check if pattern is at least found once
    if len(offsets) > 0:
        longest_repeat = 1
    for i in range(0, len(offsets)):
        j = 1
        # Check if offsets have a distance of exactly j * pattern length
        while (i+j) < len(offsets) and (offsets[i] + len(trinucleotide)*j) == offsets[i+j]:
            j += 1
            if j > longest_repeat:
                longest_repeat = j
        # Skip already found repetitions
        i += j

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

    offsets = []
    palindromes = {}
    for i in range(0, len(dna_sequence)-length):
        isPalindrome = True

        # Compare letters inward to check if they are complement of each other
        for j in range(0, int(length/2)):
            if dna_sequence[i+length-1 - j] not in complement.keys():
                isPalindrome = False
                break

            if complement[dna_sequence[i+length-1 - j]] != dna_sequence[i + j]:
                isPalindrome = False
                break

        if isPalindrome:
            offsets.append(i)

            palindrome = dna_sequence[i:i+length]
            if palindrome in palindromes.keys():
                palindromes[palindrome] += 1
            else:
                palindromes[palindrome] = 1

    plt.bar(palindromes.keys(), palindromes.values())
    plt.xticks(rotation=90)

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