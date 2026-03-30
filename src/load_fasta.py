from Bio import SeqIO


def load_sequence_database(database_fasta_path):
    """
    Load all sequences from a FASTA file into a dictionary.
    Keys are sequence IDs, values are SeqRecord objects.
    """
    sequence_database = {}

    for sequence_record in SeqIO.parse(database_fasta_path, "fasta"):
        sequence_database[sequence_record.id] = sequence_record

    if not sequence_database:
        raise ValueError(f"No sequences found in database file: {database_fasta_path}")

    print(f"Loaded {len(sequence_database)} sequences from database.")
    return sequence_database


def load_test_sequence(test_fasta_path):
    """
    Load a single test sequence from a FASTA file.
    Raises an error if the file contains more than one sequence.
    """
    sequences = list(SeqIO.parse(test_fasta_path, "fasta"))

    if len(sequences) == 0:
        raise ValueError(f"No sequence found in test file: {test_fasta_path}")

    if len(sequences) > 1:
        raise ValueError(f"Expected 1 test sequence, found {len(sequences)} in: {test_fasta_path}")

    test_sequence = sequences[0]
    print(f"Loaded test sequence: {test_sequence.id} ({len(test_sequence.seq)} bp)")
    return test_sequence