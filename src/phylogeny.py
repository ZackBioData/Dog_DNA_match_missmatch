"""
phylogeny.py
builds a neighbour-joining phylogenetic tree from one representative per breed.
deduplicates the database first — keeps the first sequence found per breed name
so the tree is readable rather than having 29 Mixed Breed leaves.

uses pairwise mismatch counting for the distance matrix all seqs are 16,735 bp
so position-by-position comparison works fine without full alignment.

--- molecular clock note ---
the obvious next step would be to use mismatch distances to estimate divergence times.

two approaches we looked at:

1. generic mammalian mtdna clock — Brown et al. (1979) PNAS estimated 0.02
   substitutions per bp per million years. easy to apply but Galtier et al. (2009)
   warns this can be off by a factor of 10 for specific lineages so probably not
   worth pretending its accurate.

2. dog-specific calibration — Frantz et al. (2016) estimated the nuclear mutation
   rate in wolves at ~4.5e-9 per bp per generation and placed dog/wolf divergence
   at 25,000-33,000 years ago. with a wolf outgroup sequence we could back-calculate
   the actual mtdna rate for this dataset from a known dated event and use that to
   date all pairwise divergences. we didnt implement this because translating nuclear
   rate to mtdna rate requires additional assumptions and adding a wolf outgroup
   would need more validation. left as a potential extension.
"""
import os
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.ticker import FixedLocator

from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from src.load_fasta import extract_breed


def deduplicate_by_breed(sequence_database):
    """
    keep one sequence per breed.
    returns a new dict keyed by breed name.
    """
    seen_breeds = {}
    for accession, record in sequence_database.items():
        breed = extract_breed(record)
        if breed not in seen_breeds:
            seen_breeds[breed] = record
    print(f"  deduplicated: {len(sequence_database)} sequences → {len(seen_breeds)} unique breeds")
    return seen_breeds


def count_mismatches(seq1_str, seq2_str):
    """position-by-position mismatch count — works since all seqs are same length"""
    return sum(a != b for a, b in zip(seq1_str, seq2_str))


def build_distance_matrix(breed_database):
    """
    pairwise mismatch distance matrix.
    distance = mismatches / seq_length, between 0 (identical) and 1 (completely different).
    lower triangle only since matrix is symmetric.
    """
    records    = list(breed_database.values())
    labels     = list(breed_database.keys())
    seq_length = len(records[0].seq)
    n          = len(records)

    print(f"\nBuilding {n}x{n} distance matrix ({n*(n-1)//2} comparisons)...")

    lower_triangle = []
    for i in range(n):
        row = []
        for j in range(i + 1):
            if i == j:
                row.append(0.0)
            else:
                mismatches = count_mismatches(str(records[i].seq), str(records[j].seq))
                row.append(mismatches / seq_length)
        lower_triangle.append(row)

        if (i + 1) % 10 == 0:
            print(f"  {i + 1}/{n} rows done")

    print("Distance matrix complete.")
    return DistanceMatrix(labels, lower_triangle)


def build_tree(distance_matrix):
    """
    neighbour-joining tree.
    NJ doesnt assume equal rates across branches, more appropriate than UPGMA here.
    midpoint rooted after since NJ produces an unrooted tree by default.
    """
    print("\nBuilding neighbour-joining tree...")
    constructor = DistanceTreeConstructor()
    tree        = constructor.nj(distance_matrix)
    tree.root_at_midpoint()
    print(f"Tree built — {len(tree.get_terminals())} leaf nodes")
    return tree


def plot_tree(tree, best_match_breed, output_dir="results"):
    """
    draw the tree.
    - best match in red, everything else steelblue
    - inner node labels hidden, only show breed names at leaves
    - x axis labelled as % divergence (mismatch distance * 100)
    - x axis reversed so most similar is on the right
    """
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # colour leaves red for best match, steelblue for everything else
    label_colours = {}
    for leaf in tree.get_terminals():
        if best_match_breed.lower() in leaf.name.lower():
            label_colours[leaf.name] = "red"
        else:
            label_colours[leaf.name] = "steelblue"

    # hide inner node labels — only show terminal breed names
    label_func = lambda node: node.name if node.is_terminal() else ""

    fig, ax = plt.subplots(figsize=(14, 18))

    Phylo.draw(
        tree,
        axes=ax,
        label_colors=label_colours,
        label_func=label_func,
        do_show=False
    )

    ax.set_title(
        "Neighbour-joining phylogenetic tree — dog mitochondrial genomes\n"
        "best match to mystery sequence highlighted in red",
        fontsize=11
    )

    # remove both axes — branch lengths are on the x axis but the labels
    # overlap badly, cleaner to just show the tree shape
    ax.set_axis_off()

    red_patch  = mpatches.Patch(color="red",       label=f"best match ({best_match_breed})")
    blue_patch = mpatches.Patch(color="steelblue", label="database sequences")
    ax.legend(handles=[red_patch, blue_patch], loc="lower right", fontsize=9)

    plt.tight_layout()
    path = os.path.join(output_dir, f"phylogeny_{timestamp}.png")
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"Phylogeny plot saved: {path}")
    return path


def run_phylogeny(sequence_database, best_match_breed, output_dir="results"):
    """entry point — called from main.py after find_best_match"""
    print("\nRunning phylogeny...")
    breed_db = deduplicate_by_breed(sequence_database)
    dm       = build_distance_matrix(breed_db)
    tree     = build_tree(dm)
    plot_tree(tree, best_match_breed, output_dir)
    return tree