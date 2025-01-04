# -------------------------------
# Abgabegruppe:
# Personen:
# HU-Accountname:
# -------------------------------
import time
from Bio import SeqIO


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
    max_score = 0

    # TODO: add your implementation here

    return alignments, max_score


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
