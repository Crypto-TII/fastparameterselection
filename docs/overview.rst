What the tool computes
======================

The problem
-----------

An LWE instance is described by four quantities:

.. list-table::
   :header-rows: 1
   :widths: 12 88

   * - Symbol
     - Meaning
   * - :math:`n`
     - the LWE dimension (in Ring-LWE terms, the polynomial degree)
   * - :math:`q`
     - the ciphertext modulus, always handled here as :math:`\log_2 q`
   * - :math:`\sigma_e`
     - the standard deviation of the error distribution
   * - :math:`\sigma_s`
     - the standard deviation of the secret distribution

together with the security level :math:`\lambda` they achieve, meaning the
best known attack costs about :math:`2^\lambda` operations.

These are not independent. Raising :math:`q` lets a scheme evaluate a deeper
circuit before the noise swamps the plaintext, but it also weakens the
instance, so :math:`n` must grow to compensate. Choosing FHE parameters means
navigating that trade-off.

The tool inverts the relation in whichever direction you need: **fix any three
and it solves for the fourth.** :math:`\sigma_s` is the exception -- it is
dictated by the scheme rather than chosen, so it is always an input.

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Option
     - Solves for
   * - ``--param "lambda"``
     - the security level, given :math:`n`, :math:`q`, :math:`\sigma_e`
   * - ``--param "n"``
     - the dimension, given :math:`\lambda`, :math:`q`, :math:`\sigma_e`
   * - ``--param "logq"``
     - the largest modulus, given :math:`\lambda`, :math:`n`, :math:`\sigma_e`
   * - ``--param "std_e"``
     - the smallest error, given :math:`\lambda`, :math:`n`, :math:`q`

The attacks
-----------

Three lattice attacks are modelled, all of them *primal*:

**uSVP** -- embed the instance in a lattice where the secret is an unusually
short vector, then run BKZ until it surfaces.

**BDD** -- treat the instance as bounded distance decoding: reduce, solve SVP
on a projected sublattice, then lift with Babai's algorithm.

**Hybrid** -- guess part of a sparse secret and solve BDD on what remains.
Only relevant when the secret has low Hamming weight, so it is used for
``--secret "sparse"``.

Dual attacks are deliberately **not** modelled; see :doc:`limitations`.

Since an attacker uses whichever is cheapest, the ``output`` column reports
the weakest link: the *lowest* security level, the *largest* dimension, the
*smallest* modulus, the *largest* error.

Two ways of answering
---------------------

For each attack the tool carries two independent routes to the same relation,
and ``--table`` shows both:

**Closed-form formulas** (columns ``usvp``, ``bdd``, and the simplified
``usvp_s``, ``bdd_s``) come from the complexity analysis in Section 3 of the
paper, with their lower-order terms fitted against the Lattice Estimator in
Section 4. They are instant and numerically stable, but the fit was done at
:math:`\sigma_e = 3.19` for binary and ternary secrets.

**Numerical solvers** (columns ``usvp num``, ``bdd num``) solve the same
systems directly with ``scipy``. They are slower and occasionally fail to
converge, but they carry no fitted constants, so they apply to any parameters.
Outside the fitted setting the tool selects them automatically; ``--num-only``
forces them.

Agreement between the two is a useful signal. Where they diverge, ``-v`` gives
the Lattice Estimator's own verdict as the tie-breaker.
