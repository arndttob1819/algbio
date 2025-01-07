# -------------------------------
# Abgabegruppe:
# Personen:
# HU-Accountname:
# -------------------------------
import time
from Bio import SeqIO


PATHS = []


def find_local_alignments(dna_sequence1, dna_sequence2, match=1, mismatch=-1, gap=-1, max_alignments=10_000):
    """
    Implements the Smith-Waterman algorithm to find all local alignments between two dna sequences.

    Input:
        dna_sequence1 = str,
        dna_sequence2 = str,
        match = int,
        mismatch = int,
        gap = int,
        max_alignments = int

    Output:
        alignments = List[Tuple[str,str]],
        max_score = int
    """
    alignments = []

    matrix = [[0 for i in range(len(dna_sequence2)+1)] for j in range(len(dna_sequence1)+1)]
    max_scores = []
    maximum_score = 0

    matrix, max_scores = compute_matrix(dna_sequence1, dna_sequence2)

    for score in max_scores:
        traceback(matrix, dna_sequence1, dna_sequence2, path = [score])

    alignments = calculate_alignments(dna_sequence1, dna_sequence2)

    return alignments, maximum_score

def compute_matrix(dna_sequence1, dna_sequence2, match=1, mismatch=-1, gap=-1):
    matrix = [[0 for i in range(len(dna_sequence2)+1)] for j in range(len(dna_sequence1)+1)]
    max_scores = []
    maximum_score = 0

    for i in range(1, len(dna_sequence1)+1):
        for j in range(1, len(dna_sequence2)+1):
            left_score = matrix[i][j-1] + gap
            up_score = matrix[i-1][j] + gap
            if dna_sequence1[i-1] == dna_sequence2[j-1]:
                diag_score = match +  matrix[i-1][j-1]
            else:
                diag_score = mismatch +  matrix[i-1][j-1]
            maximum = max(0, left_score, up_score, diag_score)
            matrix[i][j] = maximum

            # find max
            if maximum > maximum_score:
                max_scores = []
                max_scores.append([i,j])
                maximum_score = maximum
            elif maximum == maximum_score:
                max_scores.append([i,j])

    return matrix, max_scores

def traceback(matrix, dna_sequence1, dna_sequence2, path):
    if 0 in path[-1]:
        PATHS.append(path)
        return 
    
    i = path[-1][0]
    j = path[-1][1]
    left_score = matrix[i][j-1] - 1
    up_score = matrix[i-1][j] - 1
    if dna_sequence1[i-1] == dna_sequence2[j-1]:
        diag_score = 1 +  matrix[i-1][j-1]
    else:
        diag_score = -1 +  matrix[i-1][j-1]
    
    maximum = max(left_score, up_score, diag_score)

    if left_score == maximum:
        traceback(matrix, dna_sequence1, dna_sequence2, path + [[i, j-1]])
    if up_score == maximum:
        traceback(matrix, dna_sequence1, dna_sequence2, path + [[i-1, j]])
    if diag_score == maximum:
        traceback(matrix, dna_sequence1, dna_sequence2, path + [[i-1, j-1]])

def calculate_alignments(dna_sequence1, dna_sequence2):
    alignments = []
    for path in PATHS:
        string1 = ""
        string2 = ""
        last_i = 0
        last_j = 0
        path = path[::-1]
        for tile in path[1:]:
            # print(tile)
            if last_i == tile[0]:
                string1 += "-"
                string2 += dna_sequence2[tile[1]-1]
            elif last_j == tile[1]:
                string1 += dna_sequence1[tile[0]-1]
                string2 += "-"
            else:
                string1 += dna_sequence1[tile[0]-1]
                string2 += dna_sequence2[tile[1]-1]
            last_i = tile[0]
            last_j = tile[1]
        alignments.append((string1, string2))
    return alignments

def compute_hemoglobin_ranking():
    """
    Computes a ranking (most similar to least similar) of hemoglobin gene conservation between human and animals.

    Output:
        ranking = List[str],
        scores = List[int]
    """
    # FASTA file
    human_file = "human.fasta"
    comparates = ("catfish", "chimpanzee", "cow", "rat")

    ranking = []
    scores = []

    # TODO: add your implementation here
    file = comparates[0] + ".fasta"
    human_records = list(SeqIO.parse(human_file, 'fasta'))
    human_genePattern = str(human_records[0].seq)

    records1 = list(SeqIO.parse(file, 'fasta'))
    genePattern1 = str(records1[0].seq)
    print(find_local_alignments(genePattern1, human_genePattern))

    return ranking, scores


def compute_most_similar_hemoglobin_alignments():
    """
    Computes local alignments of hemoglobin gene between human and most similiar animal (regarding the ranking).

    Output:
        alignments = List[Tuple[str,str]],
        max_score = int
    """
    human_file = "human.fasta"
    comparate_file = "xxx.fasta"

    alignments = []
    max_score = 0

    # TODO: add your implementation here

    return alignments, max_score


def test_find_local_alignments():
    """
    Tests the test_find_local_alignments function against a reference result to ensure correctness.

    Input:
        dna_sequence = str
    """
    seq1 = "GATTACA"
    seq2 = "AAGATCAA"

    alignments, score = find_local_alignments(seq1, seq2, max_alignments=10)

    assert len(alignments) == 5
    assert score == 3

    assert ('GATTACA', 'GATCA-A') in alignments
    assert ('GATTACA', 'GAT--CA') in alignments
    assert ('GATTACA', 'GA-T-CA') in alignments


if __name__ == "__main__":
    # Log runtimes for competition
    runtimes = []

    # Compute sanity checks
    test_find_local_alignments()

    # Hemoglobin ranking with different animals
    runtimes.append(time.perf_counter())
    ranking, scores = compute_hemoglobin_ranking()
    runtimes[-1] = time.perf_counter() - runtimes[-1]

    print(f"Ranking: {ranking}")
    print(f"Scores: {scores}")
    assert scores == [1542, 650, 519, 166]

    # Hemoglobin alignments with most similar animal
    runtimes.append(time.perf_counter())
    alignments, score = compute_most_similar_hemoglobin_alignments()
    runtimes[-1] = time.perf_counter() - runtimes[-1]

    print(f"Found {len(alignments)} optimal alignments with maximal score: {score}")
    assert score == 1542

    print("Examples:")
    for idx, (align1, align2) in enumerate(alignments[:3]):
        print(f"\nAlignment {idx+1}:")
        print("Sequence 1 Alignment: ", align1)
        print("Sequence 2 Alignment: ", align2)

    # Competition: Set flag based on your preference
    participate_in_competition = True

    if participate_in_competition:
        print(f"\nRuntimes:")

        for idx, runtime in enumerate(runtimes):
            print(f"Task {idx+2}: {round(runtime, 2)} seconds")