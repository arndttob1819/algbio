# -------------------------------
# Abgabegruppe: 19
# Personen: Tobias Arndt, Benjamin Höhnisch, Tim Patzak
# HU-Accountname: arndttob, hoehnisb, patzakti
# -------------------------------
import time
from Bio import SeqIO

maxAlignmentsLeft = -1


def getAlignment(matrix, x, y, seq1, seq2, match, mismatch, gap, currentAlignements=None):
    global maxAlignmentsLeft

    if currentAlignements is None:
        currentAlignements = ([], [])

    alignments = []
    align1, align2 = currentAlignements

    # Follow path until score is 0
    currentScore = matrix[x][y]
    while currentScore > 0:
        # Monitors if recursive call is needed
        foundPath = False
        oldX = x
        oldY = y

        # Check for mismatch
        m = mismatch
        if seq1[y-1] == seq2[x-1]:
            m = match

        # Get scores from potential paths
        leftScore = matrix[x][y-1]+gap
        topScore = matrix[x-1][y]+gap
        digScore = matrix[x-1][y-1]+m

        if currentScore == leftScore:
            # Move left
            align1.append(seq1[y-1])
            align2.append('-')
            y -= 1
            foundPath = True

        if currentScore == topScore:
            # Move up
            if foundPath:
                # Check if max allowed alignments have been found already
                if maxAlignmentsLeft > 0:
                    maxAlignmentsLeft -= 1

                    # Search another path recursively
                    tmp1 = align1[:-1] + ['-']
                    tmp2 = align2[:-1] + [seq2[x-1]]

                    alignments += getAlignment(matrix, oldX-1, oldY, seq1, seq2, match, mismatch, gap, (tmp1, tmp2))
            else:
                align1.append('-')
                align2.append(seq2[x-1])
                x -= 1
                foundPath = True

        if currentScore == digScore:
            # Move diagonally
            if foundPath:
                # Check if max allowed alignments have been found already
                if maxAlignmentsLeft > 0:
                    maxAlignmentsLeft -= 1

                    # Search another path recursively
                    tmp1 = align1[:-1] + [seq1[oldY-1]]
                    tmp2 = align2[:-1] + [seq2[oldX-1]]

                    alignments += getAlignment(matrix, oldX-1, oldY-1, seq1, seq2, match, mismatch, gap, (tmp1, tmp2))
            else:
                align1.append(seq1[y-1])
                align2.append(seq2[x-1])
                x -= 1
                y -= 1

        # Update current score
        currentScore = matrix[x][y]

    # Format alignments
    alignments += [(''.join(align1[::-1]), ''.join(align2[::-1]))]

    return alignments


def getScoreAndAlignmentEnd(seq1, seq2, match=1, mismatch=-1, gap=-1, maxAlignments=10000):
    maxScore = 0
    alignmentEnds = []

    # Setup matrix
    matrix = [[0] * (len(seq1)+1) for i in range(len(seq2)+1)]

    for i in range(1, len(seq2)+1):
        for j in range(1, len(seq1)+1):
            # Check for match
            m = mismatch
            if seq1[j-1] == seq2[i-1]:
                m = match

            # Get scores from older paths
            leftScore = matrix[i][j-1] + gap
            topScore = matrix[i-1][j] + gap
            digScore = matrix[i-1][j-1] + m

            # Determine best score and insert value in current row
            newValue = max(0, leftScore, topScore, digScore)
            matrix[i][j] = newValue

            # Keep track of currently best alignments
            if newValue > maxScore:
                alignmentEnds = [(i, j)]
                maxScore = newValue
            elif newValue == maxScore:
                if len(alignmentEnds) < maxAlignments:
                    alignmentEnds.append((i, j))

    return matrix, maxScore, alignmentEnds


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

    # Compute matrix of patterns and get end points of highest possible alignments
    matrix, max_score, alignmentEnds = getScoreAndAlignmentEnd(dna_sequence1, dna_sequence2, match, mismatch, gap, max_alignments)

    global maxAlignmentsLeft
    maxAlignmentsLeft = max_alignments

    # Get alignments from their end points
    for k in range(len(alignmentEnds)):
        if maxAlignmentsLeft > 0:
            maxAlignmentsLeft -= 1

            i, j = alignmentEnds[k]
            alignments += getAlignment(matrix, i, j, dna_sequence1, dna_sequence2, match, mismatch, gap)

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

    # Get human sequence
    humanRecords = list(SeqIO.parse(human_file, 'fasta'))
    humanPattern = str(humanRecords[0].seq)

    # Iterate over all animals
    values = []
    for animal in comparates:
        # Get animal sequence
        file = animal + '.fasta'
        animalRecords = list(SeqIO.parse(file, 'fasta'))
        animalPattern = str(animalRecords[0].seq)

        # Get max possible score for alignment between human and animal
        _, maxScore, _ = getScoreAndAlignmentEnd(humanPattern, animalPattern)
        values.append((maxScore, animal))

    # Sort calculated scores descending and store them in ranking
    for value in sorted(values, reverse=True):
        score, animal = value
        ranking.append(animal)
        scores.append(score)

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

    # Get scores from ranking
    ranking, scores = compute_hemoglobin_ranking()

    # Get human sequence
    humanRecords = list(SeqIO.parse(human_file, 'fasta'))
    humanPattern = str(humanRecords[0].seq)

    # Get highest scoring animal pattern
    animal = ranking[0]
    animalRecords = list(SeqIO.parse(animal + '.fasta', 'fasta'))
    animalPattern = str(animalRecords[0].seq)

    # Get alignments between human and animal
    alignments, max_score = find_local_alignments(humanPattern, animalPattern)

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
