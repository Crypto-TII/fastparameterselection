import csv
import statistics
import sys
import os

# --- Usage ---
# python compare_lambda_estimates_global_split.py file1.csv file2.csv ...
filenames = ["std_e.csv"]

if not filenames:
    print("Usage: python compare_lambda_estimates_global_split.py <file1.csv> [file2.csv ...]")
    sys.exit(1)

# --- Define groups of columns to compare ---
groups = {
    "usvp": ["est usvp", "* est usvp"],
    "bdd": ["est bdd", "* est bdd"],
}

# --- Initialize results ---
results = {}
for group_name, cols in groups.items():
    results[group_name] = {
        "col_diffs": {col: [] for col in cols},
        "all_diffs": [],
        "col_worst": {col: (0, None, None, None, None) for col in cols},  # (diff, row, file, value, lambda)
        "overall_worst": (0, None, None, None, None),
    }

# --- Process all files ---
for filename in filenames:
    if not os.path.exists(filename):
        print(f"⚠️ Warning: file '{filename}' not found, skipping.")
        continue

    with open(filename, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        headers = reader.fieldnames or []
        if "lambda" not in headers:
            print(f"⚠️ Skipping '{filename}': missing 'lambda' column.")
            continue

        for i, row in enumerate(reader, start=1):
            # Skip empty lambda
            if row["lambda"] == "":
                continue

            try:
                lam = float(row["lambda"])
            except ValueError:
                continue

            # Ignore rows where lambda is 0
            if lam == 0:
                continue

            for group_name, cols in groups.items():
                for col in cols:
                    if col not in headers or row[col] == "":
                        continue
                    try:
                        val = float(row[col])
                    except ValueError:
                        continue

                    # Ignore rows where estimate is 0
                    if val == 0:
                        continue

                    diff = abs(val - lam)

                    results[group_name]["col_diffs"][col].append(diff)
                    results[group_name]["all_diffs"].append(diff)

                    # Track worst per column
                    if diff > results[group_name]["col_worst"][col][0]:
                        results[group_name]["col_worst"][col] = (diff, i, filename, val, lam)

                    # Track overall worst
                    if diff > results[group_name]["overall_worst"][0]:
                        results[group_name]["overall_worst"] = (diff, i, col, filename, lam)

# --- Compute and print stats ---
# Collect for global analysis
global_by_type = {
    "plain": [],   # e.g. est usvp, est bdd
    "star": []     # e.g. * est usvp, * est bdd
}
global_worst_plain = (0, None, None, None, None, None)  # diff, group, row, col, file, lambda
global_worst_star = (0, None, None, None, None, None)

for group_name, data in results.items():
    if not any(data["all_diffs"]):
        continue

    print(f"\n{group_name.upper()} ESTIMATES vs LAMBDA (ignoring 0s, across {len(filenames)} file(s)):")

    for col, diffs in data["col_diffs"].items():
        if not diffs:
            continue
        mean_ = statistics.mean(diffs)
        std_ = statistics.stdev(diffs) if len(diffs) > 1 else 0.0
        worst_ = data["col_worst"][col]
        print(f"  {col:12s} -> mean: {mean_:.3f}, std: {std_:.3f}, "
              f"worst diff: {worst_[0]:.3f} "
              f"(row {worst_[1]}, file '{worst_[2]}', est={worst_[3]}, lambda={worst_[4]})")

        # Add to global accumulator by type
        if col.strip().startswith("*"):
            global_by_type["star"].extend(diffs)
            if worst_[0] > global_worst_star[0]:
                global_worst_star = (worst_[0], group_name, worst_[1], col, worst_[2], worst_[4])
        else:
            global_by_type["plain"].extend(diffs)
            if worst_[0] > global_worst_plain[0]:
                global_worst_plain = (worst_[0], group_name, worst_[1], col, worst_[2], worst_[4])

    # Per-group overall
    overall_mean = statistics.mean(data["all_diffs"])
    overall_std = statistics.stdev(data["all_diffs"]) if len(data["all_diffs"]) > 1 else 0.0
    w = data["overall_worst"]
    print(f"  OVERALL     -> mean: {overall_mean:.3f}, std: {overall_std:.3f}")
    print(f"  OVERALL worst diff: {w[0]:.3f} "
          f"(row {w[1]}, column '{w[2]}', file '{w[3]}', lambda={w[4]})")

# --- Global stats by category ---
print("\n=== GLOBAL STATS ===")

# Plain (non-starred)
if global_by_type["plain"]:
    mean_plain = statistics.mean(global_by_type["plain"])
    std_plain = statistics.stdev(global_by_type["plain"]) if len(global_by_type["plain"]) > 1 else 0.0
    w = global_worst_plain
    print(f"Non-starred (est usvp, est bdd):")
    print(f"  mean: {mean_plain:.3f}, std: {std_plain:.3f}")
    print(f"  worst diff: {w[0]:.3f} (group '{w[1]}', row {w[2]}, col '{w[3]}', file '{w[4]}', lambda={w[5]})")
else:
    print("No non-starred data found.")

# Starred (* est ...)
if global_by_type["star"]:
    mean_star = statistics.mean(global_by_type["star"])
    std_star = statistics.stdev(global_by_type["star"]) if len(global_by_type["star"]) > 1 else 0.0
    w = global_worst_star
    print(f"\nStarred (* est usvp, * est bdd):")
    print(f"  mean: {mean_star:.3f}, std: {std_star:.3f}")
    print(f"  worst diff: {w[0]:.3f} (group '{w[1]}', row {w[2]}, col '{w[3]}', file '{w[4]}', lambda={w[5]})")
else:
    print("\nNo starred data found.")

