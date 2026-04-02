"""
visualise_results.py
saves outputs after the pipeline runs:
    1. txt file — best match, score, p-value, full ranked list
    2. bar chart of all 99 scores — best match in red so its obvious
    3. histogram showing the distribution of  null scores vs real score for the p-value
"""
import os
from datetime import datetime
import matplotlib.pyplot as plt


def write_text_results(results, mystery_record, pvalue, output_dir="results"):
    """
    writes the results to a txt file.
    results are just the key numbers and a ranked list.
    """
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_path = os.path.join(output_dir, f"results_{timestamp}.txt")

    with open(output_path, "w") as f:
        f.write(f"mystery sequence : {mystery_record.id}\n")
        f.write(f"length           : {len(mystery_record.seq)} bp\n")
        f.write(f"date             : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        f.write(f"best match       : {results[0]['breed']}\n")
        f.write(f"accession        : {results[0]['ID']}\n")
        f.write(f"score            : {results[0]['score']:.1f}\n")
        f.write(f"p-value          : {pvalue}\n\n")
        f.write("full rankings:\n")
        for i, r in enumerate(results, 1):
            f.write(f"  {i:03d}  {r['breed']:45s}  {r['score']:.1f}\n")

    print(f"results written to: {output_path}")
    return output_path


def plot_scores(results, output_dir="results"):
    """
    bar chart including  the alignment scores of 1 bar per species.
    mixed breeds removed
    best match bar is red, everything else is grey.
    """
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    # keep best score per breed — removes duplicate + Mixed Breed bars
    seen = {}
    for r in results:
        breed = r["breed"]
        if breed not in seen or r["score"] > seen[breed]["score"]:
            seen[breed] = r
    deduped = sorted(seen.values(), key=lambda x: x["score"], reverse=True)
    breeds = [r["breed"] for r in deduped]
    scores = [r["score"] for r in deduped]
    colours = ["red" if i == 0 else "steelblue" for i in range(len(deduped))]

    fig, ax = plt.subplots(figsize=(18, 6))
    ax.bar(range(len(breeds)), scores, color=colours)
    ax.set_xticks(range(len(breeds)))
    ax.set_xticklabels(breeds, rotation=90, fontsize=7)
    ax.set_ylabel("alignment score")
    ax.set_title("mystery sequence vs dog breed database")

    # zoom y axis — without this all bars look the same height
    ax.set_ylim(min(scores) - 50, max(scores) + 50)

    plt.tight_layout()
    path = os.path.join(output_dir, f"scores_{timestamp}.png")
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"score plot saved: {path}")
    return path


def plot_pvalue(real_score, null_scores, output_dir="results"):
    """
    histogram of null distribution vs the real score.
    if the match is good the real score (red line) should be
    way off to the right with all the null scores bunched up on the left.
    """
    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(null_scores, bins=20, color="steelblue", alpha=0.7, label="null scores (shuffled)")
    ax.axvline(x=real_score, color="red", linewidth=2, label=f"real score ({real_score:.1f})")
    ax.set_xlabel("alignment score")
    ax.set_ylabel("frequency")
    ax.set_title("permutation test — real score vs null distribution")
    ax.legend()
    plt.tight_layout()

    path = os.path.join(output_dir, f"pvalue_{timestamp}.png")
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"pvalue plot saved: {path}")
    return path


def save_all_results(results, mystery_record, pvalue, null_scores, output_dir="results"):
    """runs all three — txt, score plot, pvalue plot"""
    write_text_results(results, mystery_record, pvalue, output_dir)
    plot_scores(results, output_dir)
    plot_pvalue(results[0]["score"], null_scores, output_dir)