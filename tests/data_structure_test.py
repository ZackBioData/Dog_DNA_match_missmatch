import re
from collections import Counter
from Bio import SeqIO


# inspect data 

def extract_breed(description):
    """pull breed name out of the fasta header tags"""
    match = re.search(r'\[breed=([^\]]+)\]', description)
    return match.group(1) if match else "Unknown"



# database
print("DATABASE: dog_breeds.fa")


db_records = list(SeqIO.parse("data/dog_breeds.fa", "fasta"))

lengths = [len(r.seq) for r in db_records]
breeds = [extract_breed(r.description) for r in db_records]
breed_counts = Counter(breeds)

print(f"total sequences : {len(db_records)}")
print(f"unique breeds   : {len(breed_counts)}")
print(f"seq lengths     : min={min(lengths)} max={max(lengths)} — all same: {min(lengths)==max(lengths)}")
print()
print("breed breakdown:")
for breed, count in sorted(breed_counts.items(), key=lambda x: -x[1]):
    print(f"  {count:3d}x   {breed}")


# mystery sequence


print()
print("MYSTERY: mystery.fa")

mystery_records = list(SeqIO.parse("data/mystery.fa", "fasta"))
mystery = mystery_records[0]

print(f"accession : {mystery.id}")
print(f"breed tag : {extract_breed(mystery.description)}")
print(f"length    : {len(mystery.seq)} bp")
print(f"same length as database? {len(mystery.seq) == min(lengths)}")
