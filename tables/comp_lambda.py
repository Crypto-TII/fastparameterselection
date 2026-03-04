import csv
import statistics
import sys
import os

# --- Configuration ---
filenames = [
    "lambda_bin_2_10.csv", "lambda_bin_2_11.csv", "lambda_ter_2_10.csv", "lambda_ter_2_15.csv",
    "lambda_ter_2_16.csv", "lambda_ter_2_17.csv", "lambda_hybrid_2_13.csv", "lambda_hybrid_2_15.csv",
    "lambda_bin_2_10_num.csv", "lambda_bin_2_11_num.csv", "lambda_ter_2_10_num.csv",
    "lambda_ter_2_15_num.csv"
]

# filenames = sys.argv[1:]  # You can use command line arguments too.

if not filenames:
    print("Usage: python stats.py <file1.csv> [file2.csv ...]")
    sys.exit(1)

# --- Define possible column groups ---
groups = {
    "usvp": ["usvp", "usvp_s", "usvp num"],
    "bdd": ["bdd", "bdd_s", "bdd num"],
    "hybrid": ["hybrid"]
}

results = {}

# Initialize accumulators for each group
for group_name, cols in groups.items():
    results[group_name] = {
        "col_diffs": {col: [] for col in cols},
        "all_diffs": [],
        # (diff, row, file, value, est)
        "col_worst": {col: (0, None, None, None, None) for col in cols},
        # (diff, row, col, file, est)
        "overall_worst": (0, None, None, None, None),
    }

# --- Process all files ---
for filename in filenames:
    if not os.path.exists(filename):
        print(f"Warning: file not found '{filename}', skipping.")
        continue

    with open(filename, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        headers = reader.fieldnames or []
        rows = list(reader)

    # Detect which groups exist in this file
    available_groups = [
        g for g in groups.keys()
        if (
            g == "usvp" and (
                "est usvp" in headers or "est usvp_s" in headers or "usvp num" in headers)
        ) or (
            g == "bdd" and (
                "est bdd" in headers or "est bdd_s" in headers or "bdd num" in headers)
        ) or (
            g == "hybrid" and ("hybrid" in headers or "est hybrid" in headers)
        )
    ]

    if not available_groups:
        print(f"Skipping '{filename}': no recognized groups found.")
        continue

    for group_name in available_groups:
        est_col = f"est {group_name}"
        if est_col not in headers:
            continue  # Skip if the estimate column isn't present

        cols = groups[group_name]
        for i, row in enumerate(rows, start=1):
            try:
                est_value = float(row[est_col])
            except (ValueError, KeyError):
                continue

            for col in cols:
                if col not in row or row[col] == "":
                    continue
                try:
                    val = float(row[col])
                except ValueError:
                    continue

                diff = abs(val - est_value)
                results[group_name]["col_diffs"][col].append(diff)
                results[group_name]["all_diffs"].append(diff)

                # Track worst per column
                if diff > results[group_name]["col_worst"][col][0]:
                    results[group_name]["col_worst"][col] = (
                        diff, i, filename, val, est_value)

                # Track worst overall
                if diff > results[group_name]["overall_worst"][0]:
                    results[group_name]["overall_worst"] = (
                        diff, i, col, filename, est_value)

# --- Compute stats ---
for group_name, data in results.items():
    if not data["all_diffs"]:
        continue

    col_stats = {}
    for col, diffs in data["col_diffs"].items():
        if not diffs:
            continue
        mean_ = statistics.mean(diffs)
        std_ = statistics.stdev(diffs) if len(diffs) > 1 else 0.0
        worst_ = data["col_worst"][col]
        col_stats[col] = (mean_, std_, worst_)

    overall_mean = statistics.mean(data["all_diffs"])
    overall_std = statistics.stdev(data["all_diffs"]) if len(
        data["all_diffs"]) > 1 else 0.0

    data["col_stats"] = col_stats
    data["overall_mean"] = overall_mean
    data["overall_std"] = overall_std

# --- Print results ---
all_diffs_combined = []

for group_name, res in results.items():
    if not res["all_diffs"]:
        continue

    print(f"\n{group_name.upper()} Differences (across {len(filenames)} file(s)):")
    for col, (mean_, std_, worst_) in res["col_stats"].items():
        print(f"  {col:10s} -> mean: {mean_:.3f}, std: {std_:.3f}, "
              f"worst diff: {worst_[0]:.3f} "
              f"(row {worst_[1]}, file '{worst_[2]}', value={worst_[3]}, est={worst_[4]})")

    print(
        f"  OVERALL     -> mean: {res['overall_mean']:.3f}, std: {res['overall_std']:.3f}")
    w = res["overall_worst"]
    print(f"  OVERALL worst diff: {w[0]:.3f} "
          f"(row {w[1]}, column '{w[2]}', file '{w[3]}', est={w[4]})")

    # collect diffs for overall stats across all groups
    all_diffs_combined.extend(res["all_diffs"])

# --- Global statistics across all groups ---
if all_diffs_combined:
    global_mean = statistics.mean(all_diffs_combined)
    global_std = statistics.stdev(all_diffs_combined) if len(
        all_diffs_combined) > 1 else 0.0
    print("\n=== OVERALL STATS (USVP + BDD + HYBRID combined) ===")
    print(f"Global mean difference: {global_mean:.3f}")
    print(f"Global std deviation:   {global_std:.3f}")
else:
    print("\nNo data available for global statistics.")
