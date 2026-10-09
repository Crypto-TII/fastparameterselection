#!/bin/bash
#
# Regenerate the tables of "A Tool for Fast and Secure LWE Parameter Selection:
# the FHE case" from scratch.
#
#   bash regenerate_tables.sh              # every table
#   bash regenerate_tables.sh 2 6 10       # only the listed tables
#
# Results land in tables/ as table_NN_<name>.csv. Where one command produces
# the columns of two paper tables (a uSVP table and its BDD counterpart), the
# file is named after both; tables/README.md maps columns to tables.
#
# Every run is seeded (SEED, default 0) so the output is reproducible: the
# hybrid solver and --param std_e restart from randomised initial points.
#
# Runtime: tables 2-13 and 20 take roughly an hour in total. Tables 18 and 19
# call the Lattice Estimator's primal_hybrid and take several hours.
#
# Environment:
#   PY=python3     interpreter to use (sage --python3 on some macOS setups)
#   SEED=0         seed for the randomised solvers

set -u

PY="${PY:-python3}"
SEED="${SEED:-0}"
OUT="tables"
STD="3.19"

mkdir -p "$OUT"

_csv=""

_start() {          # _start <filename>
    _csv="$OUT/$1"
    : > "$_csv"
    echo "==> $_csv"
}

_run() {            # _run <estimate.py args...>
    echo "    $*"
    if ! $PY -m fastparameterselection.estimate "$@" --seed "$SEED" >/dev/null 2>&1; then
        echo "    !! FAILED: $*" >&2
        return 1
    fi
    if [[ ! -f output.csv ]]; then
        echo "    !! no output.csv for: $*" >&2
        return 1
    fi
    if [[ -s "$_csv" ]]; then
        tail -n +2 output.csv >> "$_csv"
    else
        cat output.csv >> "$_csv"
    fi
}

# --- Tables 2 and 4: security level, binary secret ---------------------------
# Table 2 = columns usvp / usvp_s / est usvp;  Table 4 = bdd / bdd_s / est bdd.
table_2() { table_2_4; }
table_4() { table_2_4; }
table_2_4() {
    _start "table_02_04_lambda_binary.csv"
    _run --param lambda --n 1024 --logq "20;24;25;26;27;28;30;33;37;42" \
         --secret binary --error gaussian --std "$STD" --table -v
    _run --param lambda --n 2048 --logq "37;46;50;53;54;57;62;67;74;84" \
         --secret binary --error gaussian --std "$STD" --table -v
}

# --- Tables 3 and 5: security level, ternary secret --------------------------
table_3() { table_3_5; }
table_5() { table_3_5; }
table_3_5() {
    _start "table_03_05_lambda_ternary.csv"
    _run --param lambda --n 1024 --logq "16;18;19;25;27;28;30;32;43;48" \
         --secret ternary --error gaussian --std "$STD" --table -v
    _run --param lambda --n 32768 --logq "650;760;810;880;930;1000;1050;1200;1400;1450" \
         --secret ternary --error gaussian --std "$STD" --table -v
    _run --param lambda --n 65536 --logq "1776;1229;955" \
         --secret ternary --error gaussian --std "$STD" --table -v
    _run --param lambda --n 131072 --logq "3576;2469;1918" \
         --secret ternary --error gaussian --std "$STD" --table -v
}

# --- Tables 6 and 8: LWE dimension, binary secret ----------------------------
# Table 6 = usvp / est usvp / usvp_s / est usvp_s;  Table 8 = the bdd columns.
table_6() { table_6_8; }
table_8() { table_6_8; }
table_6_8() {
    _start "table_06_08_n_binary.csv"
    _run --param n --lambda 80  --logq "42;58;71;84" --secret binary --error gaussian --std "$STD" --table -v
    _run --param n --lambda 100 --logq "34;46;57;67" --secret binary --error gaussian --std "$STD" --table -v
    _run --param n --lambda 110 --logq "31;42;52;61" --secret binary --error gaussian --std "$STD" --table -v
    _run --param n --lambda 120 --logq "28;39;48;57" --secret binary --error gaussian --std "$STD" --table -v
    _run --param n --lambda 128 --logq "27;37;45;54" --secret binary --error gaussian --std "$STD" --table -v
    _run --param n --lambda 140 --logq "24;34;41;49" --secret binary --error gaussian --std "$STD" --table -v
}

# --- Tables 7 and 9: LWE dimension, ternary secret ---------------------------
table_7() { table_7_9; }
table_9() { table_7_9; }
table_7_9() {
    _start "table_07_09_n_ternary.csv"
    for pair in 80:43 100:34 110:32 120:29 128:27 140:25 \
                80:1325 100:1025 110:1000 120:930 130:880 140:810; do
        _run --param n --lambda "${pair%%:*}" --logq "${pair#*:}" \
             --secret ternary --error gaussian --std "$STD" --table -v
    done
}

# --- Tables 10 and 14: security level from the numerical solver --------------
# Table 10 = usvp num / est usvp;  Table 14 = bdd num / est bdd.
table_10() { table_10_14; }
table_14() { table_10_14; }
table_10_14() {
    _start "table_10_14_lambda_num.csv"
    _run --param lambda --n 1024  --logq "13;18;27;32;42;54" \
         --secret binary --error gaussian --std "$STD" --table -v --num-only
    _run --param lambda --n 2048  --logq "28;32;37;53;64;84" \
         --secret binary --error gaussian --std "$STD" --table -v --num-only
    _run --param lambda --n 1024  --logq "14;19;27;34;43" \
         --secret ternary --error gaussian --std "$STD" --table -v --num-only
    _run --param lambda --n 32768 --logq "475;611;880;1050;1400" \
         --secret ternary --error gaussian --std "$STD" --table -v --num-only
}

# --- Tables 11 and 15: LWE dimension from the numerical solver ---------------
# Table 11 = usvp num / est usvp;  Table 15 = bdd num / est bdd.
table_11() { table_11_15; }
table_15() { table_11_15; }
table_11_15() {
    _start "table_11_15_n_num_ternary.csv"
    for pair in 80:43 100:34 110:32 120:29 128:27 140:25 \
                80:1325 100:1025 110:1000 120:930 130:880 140:810; do
        _run --param n --lambda "${pair%%:*}" --logq "${pair#*:}" \
             --secret ternary --error gaussian --std "$STD" --table -v --num-only
    done
}

# --- Tables 12 and 16: maximum log q -----------------------------------------
# Table 12 = logq usvp / est usvp;  Table 16 = logq bdd / est bdd.
table_12() { table_12_16; }
table_16() { table_12_16; }
table_12_16() {
    # Dimension outer, security level inner: the order the paper prints.
    _start "table_12_16_logq_binary.csv"
    for n in 1024 2048; do
        for l in 100 128 192 256; do
            _run --param logq --lambda "$l" --n "$n" --secret binary --error gaussian --std "$STD" --table -v
        done
    done
    _start "table_12_16_logq_ternary.csv"
    for n in 1024 32768; do
        for l in 100 128 192 256; do
            _run --param logq --lambda "$l" --n "$n" --secret ternary --error gaussian --std "$STD" --table -v
        done
    done
}

# --- Tables 13 and 17: minimum standard deviation of the error ---------------
# Table 13 = the usvp columns;  Table 17 = the bdd columns.
table_13() { table_13_17; }
table_17() { table_13_17; }
table_13_17() {
    _start "table_13_17_std_e_binary.csv"
    for n in 1024 2048; do
        for l in 100 128 192 256; do
            for q in 32 64; do
                _run --param std_e --lambda "$l" --n "$n" --logq "$q" --secret binary --table -v -c
            done
        done
    done
}

# --- Table 18: security level, hybrid attack ---------------------------------
table_18() {
    _start "table_18_lambda_hybrid.csv"
    for t in 8192:200:128 8192:119:128 8192:87:128 8192:210:192 8192:128:192 8192:91:192 \
             32768:850:128 32768:500:128 32768:330:128 32768:850:192 32768:565:192 32768:410:192; do
        IFS=: read -r n q h <<< "$t"
        _run --param lambda --n "$n" --logq "$q" --hw "$h" --secret sparse -v --table
    done
}

# --- Table 19: maximum log q, hybrid attack ----------------------------------
table_19() {
    _start "table_19_logq_hybrid.csv"
    for n in 1024 32768; do
        for l in 100 128 192; do
            for h in 64 128 192; do
                _run --param logq --lambda "$l" --n "$n" --hw "$h" --secret sparse -v --table
            done
        done
    done
}

# --- Table 20: parameters outside the usual range ----------------------------
table_20() {
    _start "table_20_examples.csv"
    for pair in 110:512 128:512 110:3072 128:3072; do
        _run --param logq --lambda "${pair%%:*}" --n "${pair#*:}" \
             --secret ternary --error gaussian --std "$STD" --table -v
    done
}

ALL="2_4 3_5 6_8 7_9 10_14 11_15 12_16 13_17 18 19 20"

if [[ $# -eq 0 ]]; then
    targets="$ALL"
else
    targets="$*"
fi

for t in $targets; do
    fn="table_$t"
    if ! declare -F "$fn" > /dev/null; then
        echo "no such table: $t" >&2
        exit 1
    fi
    "$fn"
done

echo "done."
