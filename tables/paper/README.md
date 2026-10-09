# Tables under the paper's own file names

The files in the directory above are named after the paper tables they feed,
and each holds one run of the tool. The paper's LaTeX source splits those same
runs by dimension or by security level, so one file above becomes several here.

Regenerate this directory with:

    python3 tables/to_paper_names.py

The columns and their order match the paper's CSVs, and the LaTeX selects
columns by name, so these can be copied straight over the paper's `tables/`
directory.

## Two differences to know about before copying

**The `n_*` files carry one extra column.** The paper's `n_bin_*.csv` and
`n_ter_*.csv` have a single `est num`; these have `est usvp num` and
`est bdd num`. The single column was a bug: one dictionary key held both
numerical estimates, so the uSVP value was overwritten by the BDD one. Tables
11 and 15 then read the same column while reporting different attacks.

An extra column is harmless on its own, because `pgfplotstable` selects by
name. But four table environments in `sections/numerical.tex` ask for the
column that no longer exists, and need editing:

| Line | Currently | Should be | Table |
| --- | --- | --- | --- |
| 181 | `columns={log q,lambda,est num,usvp num}` | `est usvp num` | 11 |
| 201 | `columns={log q,lambda,est num,usvp num}` | `est usvp num` | 11 |
| 517 | `columns={log q,lambda,est num,bdd num}` | `est bdd num` | 15 |
| 537 | `columns={log q,lambda,est num,bdd num}` | `est bdd num` | 15 |

**The secret distribution is labelled differently.** These files say `Binary`
and `Ternary` where the paper's say `Uniform (-1 0)` and `Uniform (-1 1)`. The
tool's own output changed; no table typesets that column, so this only matters
if one starts to.

## Not produced here

* `logq_hybrid_2_10.csv`, `logq_hybrid_2_15.csv` — Table 19. Its data file is
  no longer in the repository, so these cannot be regenerated.
* Table 20 has no CSV in the paper; it is typeset inline in
  `sections/tool.tex`.
