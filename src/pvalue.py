"""
pvalue.py
permutation test to check if the best match score is statistically significant.
shuffles the mystery sequence N times and realigns against the best match.
p-value = proportion of shuffled scores >= real score.
p=0.0 means none of the random shuffles matched as well — significant result.
"""
import numpy as np
from src.MSA_score import create_aligner


def shuffle_sequence(sequence_string):
    """shuffle dna sequence using numpy """
    bases = np.array(list(sequence_string))
    np.random.shuffle(bases)
    return "".join(bases)


def compute_pvalue(mystery_record, best_match_record, real_score, n_shuffles=50):
    """
    permutation test against the best match.
    shuffles mystery sequence n times, realigns each against best match,
    counts how many shuffled scores beat the real score.

    mystery_record: SeqRecord of the mystery sequence
    best_match_record: SeqRecord of the best matching db sequence
    real_score: alignment score from find_best_match
    n_shuffles: 50 is plenty — result is obvious with a score this high

    returns p-value and the full null distribution
    """
    aligner = create_aligner()
    null_scores = []

    print(f"\nRunning permutation test ({n_shuffles} shuffles)...")

    for i in range(n_shuffles):
        shuffled = shuffle_sequence(str(mystery_record.seq))
        null_score = aligner.score(shuffled, str(best_match_record.seq))
        null_scores.append(null_score)

        if (i + 1) % 10 == 0:
            print(f"  {i + 1}/{n_shuffles} done")

    # p-value = how many null scores are >= real score
    pvalue = sum(1 for s in null_scores if s >= real_score) / n_shuffles

    print(f"\nReal score    : {real_score:.1f}")
    print(f"Null mean     : {sum(null_scores)/len(null_scores):.1f}")
    print(f"Null max      : {max(null_scores):.1f}")
    print(f"P-value       : {pvalue}")

    return pvalue, null_scores