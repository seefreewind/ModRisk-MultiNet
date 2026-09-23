#!/usr/bin/env python3
"""Validate R2 headline counts from public, frozen derived tables."""

import csv
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def rows(relative):
    with (ROOT / relative).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def check(label, actual, expected):
    assert actual == expected, f"{label}: observed {actual}, expected {expected}"
    print(f"PASS {label}: {actual}")


def main():
    data = rows("results/phase3_classification_repair_R2/PAIR_CLASSIFICATION_REBUILT.tsv")
    check("all pairs", len(data), 2485)
    check("unique pairs", len({r["pair_id"] for r in data}), 2485)
    check("change classes", Counter(r["change_class"] for r in data),
          {"ATTENUATED": 900, "ENHANCED": 119, "NO_SIGNIFICANT_CHANGE": 1466})
    check("residual states", Counter(r["residual_state"] for r in data),
          {"RESIDUAL_SIGNIFICANT": 1100, "NEAR_NULL_CLASSIFIED": 605,
           "RESIDUAL_NON_SIGNIFICANT_OTHER": 780})
    significant = [r for r in data if float(r["difference_fdr"]) < .05]
    check("FDR-significant signed differences", len(significant), 1019)
    check("significant change directions", Counter(r["change_class"] for r in significant),
          {"ATTENUATED": 900, "ENHANCED": 119})
    for cls, expected in (("ATTENUATED", 495), ("ENHANCED", 47), ("NO_SIGNIFICANT_CHANGE", 558)):
        check(f"{cls} with residual FDR significance", sum(r["change_class"] == cls and
              r["residual_significant"].lower() == "true" for r in data), expected)
    check("baseline FDR-significant", sum(float(r["baseline_fdr"]) < .05 for r in data), 1535)
    check("equivalence ±0.05", sum(r["equivalent_005"].lower() == "true" for r in data), 2)
    check("equivalence ±0.10", sum(r["equivalent_010"].lower() == "true" for r in data), 162)
    print("ALL R2 CHECKS PASSED")


if __name__ == "__main__":
    main()
