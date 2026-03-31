from src.load_fasta import load_sequence_database, load_test_sequence
from src.MSA_score import find_best_match

database = load_sequence_database("data/dog_breeds.fa")
mystery  = load_test_sequence("data/mystery.fa")

results  = find_best_match(mystery, database)

