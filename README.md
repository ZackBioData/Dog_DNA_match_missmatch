# Dog DNA Sequence Identifier

> This project will take unknown mitocondrial DNA and use pairwise sequence alignment to find its most familiar relative in a database of dog breeds, compute a p-value to confirm the match isnt just noise, and build a phylogenetic tree showing how all the breeds in the database relate to each other maternally.

---

## workflow

1. Load the FASTA database of reference sequences and the mystery sequence
2. Align and score the mystery sequence against every sequence in the database using biopython PairwiseAligner
3. Rank results and identify the best match
4. Run a permutation test to compute a p-value, shuffles the mystery sequence 50 times and checks if any random shuffle matches as well as the real sequence
5. Build a neighbour-joining phylogenetic tree from a pairwise mismatch distance matrix across all breeds
6. Save everything to results/

---

## Project structure

```
Dog_DNA_match_missmatch/
│
├── data/
│   ├── dog_breeds.fa       # 99 reference sequences across 39 breeds
│   └── mystery.fa          # the sequence we want to identify
│
├── src/
│   ├── main.py             # control panel — runs the full pipeline
│   ├── load_fasta.py       # reads FASTA files into memory, extracts breed names
│   ├── sequence_alignment.py  # pairwise alignment against all 99 sequences
│   ├── pvalue.py           # permutation test
│   ├── visualise_results.py   # bar chart + p-value plot + text summary
│   └── phylogeny.py        # distance matrix + neighbour-joining tree
│
├── tests/
│   └── data_structure_test.py
│
├── results/                # output files land here (gitignored)
├── requirements.txt
└── README.txt              # you are here
```

---

## output

running the pipeline produces 4 files in `results/`:

- `results_TIMESTAMP.txt` — closest match, alignment score, p-value, full ranked list
- `scores_TIMESTAMP.png` — bar chart of all 99 alignment scores (best match in red)
- `pvalue_TIMESTAMP.png` — null distribution vs real score
- `phylogeny_TIMESTAMP.png` — neighbour-joining tree showing maternal lineage relationships

---

## how to run

install dependencies from requirements.txt
```
pip install -r requirements.txt
```
run the pipeline using "python -m src.main" from the root 


---

## note on the phylogenetic tree

The tree is built from mitochondrial DNA so it reflects maternal lineage not overall breed similarity.
two breeds that look completely different can cluster together if they share a common maternal ancestor.
this is why the Irish Soft Coated Wheaten Terrier and English Springer Spaniel sit next to each other —
their maternal lines converge even though they look nothing alike.

---

## dependencies

- `biopython` — sequence parsing, alignment, phylogenetics
- `matplotlib` — plotting
- `numpy` — faster shuffle in permutation test
- `pytest` — tests

---

## references


https://rosalind.info/problems/clus/
https://biopython.org/docs/1.75/api/Bio.Phylo.TreeConstruction.html