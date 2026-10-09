Limitations
===========

The tool reports a security level, so it is worth being precise about what
that number does and does not cover.

What is not modelled
--------------------

**Dual attacks.** Only primal attacks are considered. The efficient dual
attacks of [GJ21]_ and [MATZOV]_ rely on heuristics shown to fail in
[DP23]_, and the provable variants are weaker than the primal attacks at
these parameters. If a corrected dual analysis lands, these estimates would
need revisiting.

**Quantum attacks.** All costs are classical, in the core-SVP model with the
BDGL16 sieving cost.

**Correctness of decryption.** The tool answers only the security question.
Whether a parameter set supports the circuit you want to evaluate depends on
noise growth in your scheme, which is out of scope.

Where the formulas apply
------------------------

The constants in ``src/const.py`` were fitted against Lattice Estimator data
at :math:`\sigma_e = 3.19`, for binary and ternary secrets, over
:math:`n \in [2^{10}, 2^{15}]`. They hold up beyond that range -- the paper
checks :math:`n = 2^{16}` and :math:`2^{17}` -- but the further you go, the
more you are extrapolating. Outside the fitted setting the tool falls back to
the numerical solvers automatically.

The model was tuned to slightly *under*-estimate security rather than over-
estimate it. Discrepancies of a few bits against the Lattice Estimator are
expected and acceptable.

Known weaknesses
----------------

**The hybrid log q search is not reliable.** ``--param "logq"`` with
``--secret "sparse"`` can return a modulus whose measured security is below
the target. The solver checks each candidate against its own lambda model, but
that model and the Lattice Estimator disagree by anywhere from 2 to 100 bits,
so the check cannot underwrite the answer. Always run this combination with
``-v`` and compare ``est hybrid`` against your target.

**The Lattice Estimator loses precision for large errors.** When
:math:`\sigma_e` is far from `3.19` -- as it is for much of
``--param "std_e"`` -- the ``est`` columns should be read as indicative.

**Numerical solvers can fail to converge.** They restart from randomised
initial points and give up after a bounded number of attempts. A failure is
reported as a fallback value with an estimator verdict of `0`, which is a
signal to rerun with a different ``--seed`` or to use the closed forms.

**Local optima.** The solvers may settle on a solution that is not the best
available. Restarting with a different ``--seed``, or raising ``--nrestart``
for the hybrid, sometimes finds a better one.

Verifying a result
------------------

The tool is a fast way to narrow the search, not a replacement for an
estimator. Once you have candidate parameters, confirm them:

.. code-block:: bash

   fastparams --param "lambda" --n "1024" --logq "27" \
      --secret "binary" --error "gaussian" --std "3.19" --table -v

The ``est`` columns are the Lattice Estimator's verdict under the BDGL16 cost
model. You can also check against an independent tool such as the
`Leaky-LWE Estimator <https://github.com/lducas/leaky-LWE-Estimator>`_.

.. [GJ21] Guo and Johansson, *Faster dual lattice attacks for solving LWE
   with applications to CRYSTALS*, ASIACRYPT 2021.
.. [MATZOV] MATZOV, *Report on the security of LWE: improved dual lattice
   attack*, 2022.
.. [DP23] Ducas and Pulles, *Does the dual-sieve attack on learning with
   errors even work?*, CRYPTO 2023.
