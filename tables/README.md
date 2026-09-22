# Regenerated tables

Everything here is produced by `regenerate_tables.sh` at the repository root:

```bash
bash regenerate_tables.sh            # all tables
bash regenerate_tables.sh 2_4 18     # selected tables
```

Runs are seeded (`SEED`, default 0), so repeated invocations give identical
files. Without a seed the hybrid solver and `--param std_e` restart from
randomised initial points and drift by a bit or two between runs.

## Which file holds which table

Several paper tables are two views of one run: a uSVP table and its BDD
counterpart share a command and differ only in the columns quoted. Those are
stored once, under both numbers.

| File | Paper table | Columns for that table |
| --- | --- | --- |
| `table_02_04_lambda_binary.csv` | 2 | `usvp` = Eq. (25), `usvp_s` = Eq. (27), `est usvp` |
| | 4 | `bdd` = Eq. (28), `bdd_s` = Eq. (31), `est bdd` |
| `table_03_05_lambda_ternary.csv` | 3 | `usvp`, `usvp_s`, `est usvp` |
| | 5 | `bdd`, `bdd_s`, `est bdd` |
| `table_06_08_n_binary.csv` | 6 | `usvp` = Eq. (32), `usvp_s` = Eq. (33), `est usvp`, `est usvp_s` |
| | 8 | `bdd` = Eq. (34), `bdd_s` = Eq. (35), `est bdd`, `est bdd_s` |
| `table_07_09_n_ternary.csv` | 7 | `usvp`, `usvp_s`, `est usvp`, `est usvp_s` |
| | 9 | `bdd`, `bdd_s`, `est bdd`, `est bdd_s` |
| `table_10_14_lambda_num.csv` | 10 | `usvp num`, `est usvp` |
| | 14 | `bdd num`, `est bdd` |
| `table_11_15_n_num_ternary.csv` | 11 | `usvp num`, `est usvp` |
| | 15 | `bdd num`, `est bdd` |
| `table_12_16_logq_binary.csv`, `table_12_16_logq_ternary.csv` | 12 | `logq usvp`, `est usvp` |
| | 16 | `logq bdd`, `est bdd` |
| `table_13_17_std_e_binary.csv` | 13 | `log2(std_e) usvp`, `est usvp`, `* log2(std_e) usvp`, `* est usvp`, `est calls usvp` |
| | 17 | `log2(std_e) bdd`, `est bdd`, `* log2(std_e) bdd`, `* est bdd`, `est calls bdd` |
| `table_18_lambda_hybrid.csv` | 18 | `hybrid`, `est hybrid` |
| `table_20_examples.csv` | 20 | `output` |

Columns prefixed `est ` come from the Lattice Estimator with the BDGL16 cost
model, reported as `floor(log2(rop))`. A `*` prefix marks a value after the
`-c` correction loop.

## Status

The BDD `eta` question is settled: the code is correct and Section 5 of the
paper is not. Everything in this directory uses the settled equation.

Section 5's BDD system should read `lambda - (0.292*eta + 16.4) = 0` for its
third equation, matching `T_SVP(eta) = 2^(0.292*eta + 16.4)` in Section 3.2.
Tables 14, 15, 16 and 17 in the paper must be replaced from the files here.

## Known issues these files will show

- `--param logq` reports `output = min(logq usvp, logq bdd)`. A larger modulus
  is weaker, so the binding constraint is the attack tolerating the smaller
  one. Note that `est usvp` and `est bdd` are each measured at their own
  attack's log q rather than at `output`, so read the row as two independent
  results plus a recommendation, not as three views of one value.

- Table 20's four rows do not reproduce the values printed in the paper, and
  two of the printed rows are short of their claimed security level: (110,
  512, 17) measures 102 bits and (110, 3072, 101) measures 101 bits.
