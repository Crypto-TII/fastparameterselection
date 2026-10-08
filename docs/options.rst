Command line reference
======================

Every invocation names the quantity to solve for with ``--param`` and supplies
the others. Run ``fastparams -h`` for the same list in the
terminal.

Choosing what to solve for
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Option
     - Meaning
   * - ``--param <p>``
     - ``lambda``, ``n``, ``logq``, ``std_e``, or ``est`` to run the Lattice
       Estimator directly and print its own report
   * - ``--n <n>``
     - LWE dimension, e.g. ``1024``
   * - ``--lambda <l>``
     - target security level, e.g. ``128``
   * - ``--logq <q>``
     - :math:`\log_2 q`. Accepts a list and ranges: ``20;24-28;30``
   * - ``--hw <h>``
     - Hamming weight of a sparse secret, e.g. ``128``

The parameter named by ``--param`` is the one being solved for, so it is the
one you leave out.

Distributions
-------------

Options without the ``s-`` prefix configure the **error**; the ``s-`` prefixed
form configures the **secret**. See :doc:`tutorial_non_fhe` for the full
table.

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Option
     - Meaning
   * - ``--secret <d>``
     - ``binary``, ``ternary``, ``sparse``, ``uniform``, ``uniformmod``,
       ``gaussian``, ``binomial``
   * - ``--error <d>``
     - same choices; defaults to ``gaussian``
   * - ``--std <s>`` / ``--s-std <s>``
     - standard deviation for ``gaussian``; defaults to ``3.19``
   * - ``--eta <e>`` / ``--s-eta <e>``
     - width for ``binomial``; defaults to ``1``
   * - ``-a <a>`` ``-b <b>`` / ``--s-a`` ``--s-b``
     - bounds for ``uniform``
   * - ``--q <p>``
     - modulus for ``uniformmod``

Output
------

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Option
     - Meaning
   * - ``--table``
     - show every formula and solver separately rather than one recommendation
   * - ``-v``
     - also run the Lattice Estimator and show its verdict in the ``est``
       columns
   * - ``-c``
     - walk the result towards the target using the Lattice Estimator; see
       :doc:`tutorial_correction`
   * - ``--num-only``
     - use only the numerical solvers, skipping the fitted closed forms

Results are also written to ``output.csv`` in the working directory.

Solver control
--------------

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Option
     - Meaning
   * - ``--seed <n>``
     - base seed for the randomised solvers; defaults to ``0``. See
       :doc:`reproducibility`
   * - ``--nrestart <n>``
     - restart budget for the hybrid optimisation
   * - ``--coreSVP <m>``
     - cost model, ``BDGL`` or ``MATZOV``. Work in progress
   * - ``--mitm``
     - meet-in-the-middle guessing in the hybrid attack. Work in progress
   * - ``--ntru``
     - check for the overstretched NTRU regime. **Currently broken**: it
       reports every parameter set as overstretched

Refitting
---------

These apply only with ``--fit``, which relearns the constants in
``src/const.py`` from ``dataset/``:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Option
     - Meaning
   * - ``--fit``
     - refit rather than estimate
   * - ``--attack <a>``
     - ``usvp`` or ``bdd``
   * - ``--simpl <0|1>``
     - ``0`` for the full formula, ``1`` for the simplified one
