from src.load_fasta import load_sequence_database, load_test_sequence
from src.sequence_alignment import find_best_match
from src.pvalue import compute_pvalue
from src.visualise_results import save_all_results
from src.phylogeny import run_phylogeny

database = load_sequence_database("data/dog_breeds.fa")
mystery  = load_test_sequence("data/mystery.fa")

results  = find_best_match(mystery, database)
best_match_record = database[results[0]["ID"]]
pvalue, null_scores = compute_pvalue(mystery, best_match_record, results[0]["score"])
save_all_results(results, mystery, pvalue, null_scores)
tree = run_phylogeny(database, results[0]["breed"])