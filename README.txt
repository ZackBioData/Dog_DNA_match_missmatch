# Dog DNA Sequence Identifier

> This project will take unknown DNA and use Multiple sequence alignment to show the most its most familiar relative to a collection of dogs in a database, and a metric to compute how distant it is from the other sequences in the database


---

## workflow

1. Loads a FASTA database of reference sequences and a test sequence
2. Align and score the test sequence against every sequence in the database using biopython .align
3. Compute a p-value 
4. Build a visual showing the distance between test sequence and the database sequences
---

## Project structure

```
dna_identifier/
│
├── data/                    
│   ├── dog_breeds.fasta       # Reference sequences 
│   └── mystery.fasta  # The sequence we want to identify
│
├── src/                     # Main folder for all different functions
│   ├── main.py              # control panel runs the full pipeline
│   ├── load_fasta.py        # Reads FASTA files into memory
│   ├── msa_score.py         # Alignment scoring + p-value
│   ├── rank_similarity_scores.py# Ranks all sequences, picks the winner
│   └── write_results.py     # Saves text output + generates plots
│   └── phylogeny.py         # https://www.ebi.ac.uk/jdispatcher/phylogeny              https://biopython.org/docs/1.75/api/Bio.Phylo.TreeConstruction.html
│
│
├── results/                 # Outputs (gitignored)
├── requirements.txt         # Python dependencies
├── DESIGN.md                # Architecture, pseudocode, and design decisions
└── README.md                # You are here
```

---

## Output

Running the pipeline produces four files in `results/`:

- `similarity_results.txt` — text summary of the closest match, score, and p-value
- `alignment_scores.png` — bar chart of all alignment scores (closest match in red)
- `pvalue.png` — null distribution vs real score
- `phylogenetic_tree.png` — tree 

---

## Dependencies

- `biopython` — all logic
- `matplotlib` — plotting

---

## references
https://rosalind.info/problems/clus/

