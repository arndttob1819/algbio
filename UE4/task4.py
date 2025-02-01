# -------------------------------
# Abgabegruppe: 19
# Personen: Tobias Arndt, Benjamin Höhnisch, Tim Patzak
# HU-Accountname: arndttob, hoehnisb, patzakti
# -------------------------------
import time
from Bio import SeqIO
from itertools import product

import re


# Find all occurrences of a substring in a string
def findAll(txt, word):
    return [match.start() for match in re.finditer(word, txt)]


# Get the score between two words of same length
def getScore(word1, word2, sequenceMatrix):
    if len(word1) != len(word2):
        return -1

    score = 0
    for i in range(len(word1)):
        score += sequenceMatrix[(word1[i], word2[i])]

    return score


# Get all acceptable variants of given word
def getVariants(word, threshold, sequenceMatrix):
    alphabet = ['A', 'C', 'G', 'T']
    acceptableVariants = []

    # Get all possible variants of word
    variants = [''.join(comb) for comb in product(alphabet, repeat=len(word))]

    for variant in variants:
        # The word is not a variant of itself
        if word == variant:
            continue

        score = getScore(word, variant, sequenceMatrix)
        # If score is not lower than threshold its acceptable
        if score >= threshold:
            acceptableVariants.append(variant)

    return acceptableVariants


def generate_words(sequence, word_length):
    '''
    Generates all possible subwords of length word_length from sequence.

    Input:
        sequence = str,
        word_length = int

    Output:
        words = List[str]
    '''
    words = [sequence[i:i+word_length] for i in range(len(sequence)-word_length+1)]

    return words


def search_sequence(words, sequence, threshold, scoring_matrix):
    '''
    Searches for matches of words in sequence. Returns tuples of matching words
    with a score above or equal to the threshold and their offsets.

    Input:
        words = List[str],
        sequence = str,
        threshold = int,
        scoring_matrix = Dict[Tuple[str,str],int]

    Output:
        hits = List[Tuple[str, int]]
    '''
    hits = []

    # Count offset from subword
    i = 0
    for word in words:
        for hit in findAll(sequence, word):
            hits.append((word, hit, i))

        # Get all variants of a words that are still acceptable
        variants = getVariants(word, threshold, scoring_matrix)
        for variant in variants:
            for hit in findAll(sequence, variant):
                hits.append((variant, hit, i))

        i += 1
            
    return hits


def extend_hit(query_sequence, db_sequence, hit_position_query, hit_position_db, word_length, scoring_matrix):
    '''
    Extends a hit (matching word) to find its maximal segment pair.

    Input:
        query_sequence = str,
        db_sequence = str,
        hit_position_query = int,
        hit_position_db = int,
        word_length = int,
        scoring_matrix = Dict[Tuple[str,str],int]

    Output:
        alignment = Dict[str, int],
        max_score = int
    '''
    max_score = 0
    current_score = 0

    # Setup current best alignment
    alignment = {
        'query_start': hit_position_query,
        'db_start': hit_position_db,
        'query_end': hit_position_query + word_length-1,
        'db_end': hit_position_db + word_length-1,
        'length': word_length
    }

    threshold = 0

    i = 1
    length = word_length
    # Look left for extensions
    while hit_position_query-i >= 0 and hit_position_db-i >= 0:
        # Calculate new score
        current_score += scoring_matrix[(query_sequence[hit_position_query-i], db_sequence[hit_position_db-i])]
        length += 1

        # Update new best alignment
        if current_score > max_score:
            max_score = current_score
            alignment['query_start'] = hit_position_query - i
            alignment['db_start'] = hit_position_db - i
            alignment['length'] = length
        # Stop looking left if current score is far too low
        elif current_score < max_score-threshold:
            break

        i += 1

    i = 1
    length = alignment['length']
    current_score = max_score

    queryHitEnd = hit_position_query + word_length-1
    dbHitEnd = hit_position_db + word_length-1
    # Look right for extensions
    while queryHitEnd+i < len(query_sequence) and dbHitEnd+i < len(db_sequence):
        # Calculate new score
        current_score += scoring_matrix[(query_sequence[queryHitEnd+i], db_sequence[dbHitEnd+i])]
        length += 1

        # Update new best alignment
        if current_score > max_score:
            max_score = current_score
            alignment['query_end'] = queryHitEnd + i
            alignment['db_end'] = dbHitEnd + i
            alignment['length'] = length
        # Stop looking right if current score is far too low
        elif current_score < max_score-threshold:
            break

        i += 1

    # Add score of subword to final score
    max_score += getScore(query_sequence[hit_position_query:queryHitEnd+1], db_sequence[hit_position_db:dbHitEnd+1], scoring_matrix)

    return alignment, max_score


def blast(query_sequence, db_sequences, scoring_matrix, word_length=3, threshold=1, alphabet='ACGT', top_k=10):
    '''
    Executes the BLAST algorithm. Returns best top_k maximal segment pairs with descending score.

    Input:
        query_sequence = str,
        db_sequences = List[Tuple[str,str]],
        word_length = int,
        threshold = int,
        scoring_matrix = Dict[Tuple[str,str],int]
        alphabet = str,
        top_k = int

    Output:
        msps = List[Tuple[str, Dict[str,int], int]]
    '''
    msps = []

    words = generate_words(query_sequence, word_length)
    currentMinScore = 10000000000
    for seq in db_sequences:
        name, pattern = seq

        hits = search_sequence(words, pattern, threshold, scoring_matrix)

        foundAlignments = set()
        for hit in hits:
            word, posDB, posQuery = hit
            alignment, score = extend_hit(query_sequence, pattern, posQuery, posDB, word_length, scoring_matrix)

            # Check if element should be considered for top k scores
            # Either if not enough elements in list or score is higher than lowest score in list
            if len(foundAlignments) < top_k or score > currentMinScore:
                # Keep track of lowest added score
                if score < currentMinScore:
                    currentMinScore = score

                # Prevent duplicate alignments
                alignment['score'] = score
                alignmentElements = tuple(alignment.items())
                if alignmentElements not in foundAlignments:
                    foundAlignments.add(alignmentElements)

        # Reconstruct alignments from set
        for result in foundAlignments:
            alignment = {}
            score = 0
            for element in result:
                key, value = element
                if key == 'score':
                    score = value
                else:
                    alignment[key] = value

            msps.append((name, alignment, score))

    # Sort alignments with score
    msps.sort(key=lambda x: x[2], reverse=True)

    return msps[:top_k]


def find_blaTEM1_msps():
    '''
    Searches for blaTEM-1 gene in 10 different bacteria sequences using BLAST.

    Output:
        msps = List[Tuple[str, Dict[str,int], int]]
    '''
    blaTEM1_file = 'blaTEM-1.fasta'
    db_sequence_names = (
        'escherichia_coli', 'klebsiella_pneumoniae', 'proteus_mirabilis', 'salmonella_enterica', 'enterobacter_cloacae',
        'haemophilus_influenzae', 'serratia_marcescens', 'pseudomonas_aeruginosa', 'acinetobacter_baumannii', 'citrobacter_freundii'
    )

    msps = []

    blaTEMRecord = list(SeqIO.parse(blaTEM1_file, 'fasta'))
    querySequence = str(blaTEMRecord[0].seq)

    dbSequences = []

    for name in db_sequence_names:
        record = list(SeqIO.parse(name + '.fasta', 'fasta'))
        pattern = str(record[0].seq)
        dbSequences.append((name, pattern))

    msps = blast(querySequence, dbSequences, scoring_matrix)

    return msps


def test_blast(scoring_matrix):
    """
    Tests the BLAST function to ensure correctness.

    Input:
        scoring_matrix = Dict[Tuple[str,str],int]
    """
    query_sequence = "ATGCG"
    db_sequences = [("test", "ATGACGATGCGT")]

    msps = blast(query_sequence, db_sequences, scoring_matrix)
    best_msp, second_best_msp = msps[0], msps[1]

    assert best_msp[1] == {'query_start': 0, 'db_start': 6, 'query_end': 4, 'db_end': 10, 'length': 5}
    assert second_best_msp[1] == {'query_start': 0, 'db_start': 0, 'query_end': 2, 'db_end': 2, 'length': 3}


if __name__ == "__main__":
    # Scoring matrix
    scoring_matrix = {
        ('A', 'A'): 1,
        ('C', 'C'): 1,
        ('G', 'G'): 1,
        ('T', 'T'): 1,
        ('A', 'C'): -1,
        ('A', 'G'): -1,
        ('A', 'T'): -1,
        ('C', 'A'): -1,
        ('C', 'G'): -1,
        ('C', 'T'): -1,
        ('G', 'A'): -1,
        ('G', 'C'): -1,
        ('G', 'T'): -1,
        ('T', 'A'): -1,
        ('T', 'C'): -1,
        ('T', 'G'): -1,
    }

    # Compute sanity checks
    test_blast(scoring_matrix)

    # Compute MSPs for blaTEM-1
    runtime = time.perf_counter()
    msps = find_blaTEM1_msps()
    runtime = time.perf_counter() - runtime

    # Show results
    for seq_name, alignment, score in msps:
        print(f"Sequence: {seq_name}: Alignment: Query Start={alignment['query_start']}, Query End={alignment['query_end']}, "
              f"DB Start={alignment['db_start']}, DB End={alignment['db_end']}, "
              f"Länge={alignment['length']}, Score={score}")

    # Competition: Set flag based on your preference
    participate_in_competition = True

    if participate_in_competition:
        print(f"\nTask 4: {round(runtime, 2)} seconds")
