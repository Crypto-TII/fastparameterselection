import csv
import statistics
import sys
import os

# --- Configuration ---
# You can run with:
#   python stats_lambda_all.py file1.csv file2.csv ...
#filenames = sys.argv[1:]

filenames = [
    "n_bin_80.csv", "n_bin_100.csv", "n_bin_110.csv", "n_bin_120.csv",
    "n_bin_128.csv", "n_bin_140.csv", "n_ter_2_10.csv", "n_ter_2_15.csv"
]

if not filenames:
    print("Usage: python stats_lambda_all.py <file1.csv> [file2.csv ...]")
    sys.exit(1)

# --- Define which columns belong to which group ---
groups = {
    "usvp": ["est usvp", "est usvp_s", "est num"],
    "bdd": ["est bdd", "est bdd_s", "est bdd num"],
}

# --- Initialize data containers ---
results = {}
for group_name, cols in groups.items():
    results[group_name] = {
        "col_diffs": {col: [] for col in cols},
        "all_diffs": [],
        "col_worst": {col: (0, None, None, None, None) for col in cols},  # (diff, row, file, value, lambda)
        "overall_worst": (0, None, None, None, None),
    }

# --- Process each file ---
for filename in filenames:
    if not os.path.exists(filename):
        print(f"⚠️  Warning: file not found '{filename}', skipping.")
        continue

    with open(filename, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        headers = reader.fieldnames or []
        if "lambda" not in headers:
            print(f"⚠️  Skipping '{filename}': missing 'lambda' column.")
            continue

        for i, row in enumerate(reader, start=1):
            try:
                lam = float(row["lambda"])
            except ValueError:
                continue

            for group_name, cols in groups.items():
                for col in cols:
                    if col not in headers or row[col] == "":
                        continue
                    try:
                        val = float(row[col])
                    except ValueError:
                        continue

                    diff = abs(val - lam)

                    results[group_name]["col_diffs"][col].append(diff)
                    results[group_name]["all_diffs"].append(diff)

                    # Track worst per column
                    if diff > results[group_name]["col_worst"][col][0]:
                        results[group_name]["col_worst"][col] = (diff, i, filename, val, lam)

                    # Track overall worst for the group
                    if diff > results[group_name]["overall_worst"][0]:
                        results[group_name]["overall_worst"] = (diff, i, col, filename, lam)

# --- Compute and print results ---
all_diffs_combined = []  # store all diffs across all groups

for group_name, data in results.items():
    if not any(data["all_diffs"]):
        continue

    print(f"\n{group_name.upper()} ESTIMATES vs LAMBDA (across {len(filenames)} file(s)):")

    # Per-column stats
    for col, diffs in data["col_diffs"].items():
        if not diffs:
            continue
        mean_ = statistics.mean(diffs)
        std_ = statistics.stdev(diffs) if len(diffs) > 1 else 0.0
        worst_ = data["col_worst"][col]
        print(f"  {col:10s} -> mean: {mean_:.3f}, std: {std_:.3f}, "
              f"worst diff: {worst_[0]:.3f} "
              f"(row {worst_[1]}, file '{worst_[2]}', value={worst_[3]}, lambda={worst_[4]})")

    # Overall stats per group
    overall_mean = statistics.mean(data["all_diffs"])
    overall_std = statistics.stdev(data["all_diffs"]) if len(data["all_diffs"]) > 1 else 0.0
    w = data["overall_worst"]
    print(f"  OVERALL     -> mean: {overall_mean:.3f}, std: {overall_std:.3f}")
    print(f"  OVERALL worst diff: {w[0]:.3f} "
          f"(row {w[1]}, column '{w[2]}', file '{w[3]}', lambda={w[4]})")

    # Add to global stats
    all_diffs_combined.extend(data["all_diffs"])

# --- Global statistics across all groups ---
if all_diffs_combined:
    global_mean = statistics.mean(all_diffs_combined)
    global_std = statistics.stdev(all_diffs_combined) if len(all_diffs_combined) > 1 else 0.0
    print("\n=== OVERALL STATS (USVP + BDD combined) ===")
    print(f"Global mean difference: {global_mean:.3f}")
    print(f"Global std deviation:   {global_std:.3f}")
else:
    print("\nNo data available for global statistics.")

