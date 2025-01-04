# -------------------------------
# Abgabegruppe: 19
# Personen: Tobias Arndt, Benjamin Höhnisch, Tim Patzak
# HU-Accountname: arndttob, hoehnisb, patzakti
# -------------------------------
import time
from Bio import SeqIO
from collections import defaultdict


def radix(elements, k):
    # Setup alphabet
    keySet = {'$': 0, 'A': 1, 'C': 2, 'G': 3, 'N': 4, 'T': 5}
    sortedElements = elements

    # Check letters from right to left
    for i in range(k-1, -1, -1):
        # Put elements in bucket according to current letter
        buckets = [[], [], [], [], [], []]
        for element in sortedElements:
            key = keySet[element[i]]
            buckets[key].append(element)

        # Merge buckets as one list
        sortedElements = []
        for bucket in buckets:
            sortedElements += bucket

    return sortedElements

# Sorts content of a bucket recursively
def sortBucket(seq, bucket, k):
    sa = []

    # Create new needed buckets as lists
    buckets = defaultdict(list)

    # Put suffixes in buckets according to current depth
    for i in bucket:
        # Pad key with $ to the right if not long enough
        key = seq[i:i+k].ljust(k, '$')
        buckets[key].append(i)

    # Sort bucket keys with radix sort
    sortedKeys = radix(buckets.keys(), k)

    for key in sortedKeys:
        # If bucket has more than one element sort that bucket
        if len(buckets[key]) > 1:
            sa += sortBucket(seq, buckets[key], k*2)
        # Otherwise add to current final result
        else:
            sa += buckets[key]

    return sa


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

    # Start recursive sorting with whole string in one bucket
    sa = sortBucket(dna_sequence, sa, k)

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

    low = 0
    high = len(suffix_array)
    oldMiddle = -1
    while True:
        # Get current suffix to check
        middle = int(low + ((high-low)/2))
        suffix = dna_sequence[suffix_array[middle]:]

        # End search if it got stuck
        if middle == oldMiddle:
            break

        found = True
        # Compare char by char pattern with suffix
        for i in range(len(dna_pattern)):
            # Is suffix shorter than pattern look right
            if i >= len(suffix):
                found = False
                low = middle + 1
                break
            # Is current pattern char greater than suffix char look right
            elif dna_pattern[i] > suffix[i]:
                found = False
                low = middle + 1
                break
            # Is current pattern char lower than suffix char look left
            elif dna_pattern[i] < suffix[i]:
                found = False
                high = middle
                break

        # Keep track of last check position
        oldMiddle = middle

        if found:
            # Add found position to results
            result.append(suffix_array[middle])

            # Check right neighbours for pattern
            counter = 1
            while found:
                neighbour = dna_sequence[suffix_array[middle+counter]:]
                # Check if pattern is prefix of neighbour
                for i in range(len(dna_pattern)):
                    if i >= len(neighbour) or dna_pattern[i] != neighbour[i]:
                        found = False
                        break

                if found:
                    # If pattern was found again check next right neighbour
                    result.append(suffix_array[middle+counter])
                    counter += 1

            # Check left neighbours for pattern
            found = True
            counter = 1
            while found:
                neighbour = dna_sequence[suffix_array[middle-counter]:]
                # Check if pattern is prefix of neighbour
                for i in range(len(dna_pattern)):
                    if i >= len(neighbour) or dna_pattern[i] != neighbour[i]:
                        found = False
                        break

                if found:
                    # If pattern was found again check next left neighbour
                    result.append(suffix_array[middle-counter])
                    counter += 1

            # End search after all occurrences where found
            break

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

    # Loop over all genes
    for gene in genes:
        # Read gene file and extract pattern
        file = gene + '.fasta'
        records = list(SeqIO.parse(file, 'fasta'))
        genePattern = str(records[0].seq)

        # Search for gene pattern in sequence
        offsets += find_all(dna_sequence, suffix_array, genePattern)

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
    print(sorted(_builtin_find_all(dna_sequence, dna_pattern)))
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
    
