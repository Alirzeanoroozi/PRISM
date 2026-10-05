"""Plot distribution of num_chains from templates_analysis_sorted_unique.csv."""

import csv
import os
from collections import defaultdict

import matplotlib.pyplot as plt

CSV_PATH = "templates/templates_analysis_sorted_unique.csv"
OUTPUT_PATH = "templates/num_chains_distogram.png"

COLORS = {"homo": "#4C72B0", "hetero": "#DD8452"}


def main():
    base = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(base, CSV_PATH)
    output_path = os.path.join(base, OUTPUT_PATH)

    by_chains = defaultdict(lambda: {"homo": 0, "hetero": 0, "unknown": 0})
    with open(csv_path, newline="") as f:
        for row in csv.DictReader(f):
            n = int(row["num_chains"])
            kind = row.get("homo_hetero", "unknown").strip().lower()
            if kind not in COLORS:
                kind = "unknown"
            by_chains[n][kind] += 1

    xs = sorted(by_chains)
    homo = [by_chains[x]["homo"] for x in xs]
    hetero = [by_chains[x]["hetero"] for x in xs]

    fig, ax = plt.subplots(figsize=(12, 5))
    ax.bar(
        xs,
        homo,
        width=0.85,
        label="homo",
        color=COLORS["homo"],
        edgecolor="white",
        linewidth=0.5,
    )
    ax.bar(
        xs,
        hetero,
        width=0.85,
        bottom=homo,
        label="hetero",
        color=COLORS["hetero"],
        edgecolor="white",
        linewidth=0.5,
    )
    ax.set_xlabel("Number of chains")
    ax.set_ylabel("Template count")
    ax.set_title("Distribution of chain counts (templates)")
    ax.set_xticks(xs)
    ax.legend(title="Template chains")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    n = sum(h + t for h, t in zip(homo, hetero))
    n_unknown = sum(by_chains[x]["unknown"] for x in xs)
    print(f"saved -> {output_path}  (n={n} homo+hetero, unknown={n_unknown})")


if __name__ == "__main__":
    main()
