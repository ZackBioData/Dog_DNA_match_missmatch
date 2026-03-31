from src.load_fasta import load_sequence_database, load_test_sequence
from src.MSA_score import find_best_match
from src.pvalue import compute_pvalue

database = load_sequence_database("data/dog_breeds.fa")
mystery  = load_test_sequence("data/mystery.fa")

results  = find_best_match(mystery, database)

best_match_record = database[results[0]["ID"]]
pvalue, null_scores = compute_pvalue(mystery, best_match_record, results[0]["score"])