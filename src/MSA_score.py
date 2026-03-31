"""
compute_similarity.py
compute pairwise alignment scores between the mystery sequence and all db sequences.
all sequences are 16,735 bp so no normalisation needed — raw score is fine.
https://biopython.org/docs/1.76/api/Bio.Align.html
"""
from Bio.Align import PairwiseAligner


def create_aligner():
    """biopython presets"""
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2
    aligner.mismatch_score = -4
    aligner.open_gap_score = -4
    aligner.extend_gap_score = -1
    return aligner


def find_best_match(mystery_record, sequence_database):
    """
    align and score mystery sequence against every sequence in the database.
    returns results decending from best match down.
    mystery_record: SeqRecord of the mystery sequence
    sequence_database: dict of {breed_name: SeqRecord}
    list of dicts sorted by score descending
    """
    aligner = create_aligner()
    results = []

    print(f"\nAligning against database sequences")

    for breed_name, db_record in sequence_database.items():
        score = aligner.score(str(mystery_record.seq), str(db_record.seq))

        results.append({
            "breed": breed_name,
            "ID": db_record.id,
            "score": score,
        })

        print(f"  {breed_name:40s}  {score:.1f}")

    results.sort(key=lambda x: x["score"], reverse=True)

    best = results[0]
    print(f"\nBest match: {best['breed']} ({best['ID']})")
    print(f"Score: {best['score']:.1f}")

    return results