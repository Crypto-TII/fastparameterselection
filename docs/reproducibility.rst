Reproducibility
===============

Seeding
-------

Several solvers restart their optimisation from randomised initial points: the
hybrid attack, and ``--param "std_e"``. Each such call derives its own
generator from a base seed together with the call's own arguments, so a given
set of LWE parameters always produces the same answer, independently of how
many other rows were computed in the same invocation:

.. code-block:: bash

   # these agree
   python3 src/estimate.py --param "std_e" --lambda "128" --n "1024" --logq "32;48;64" --secret "binary"
   python3 src/estimate.py --param "std_e" --lambda "128" --n "1024" --logq "32" --secret "binary"

Pass ``--seed`` to change the base seed. See ``src/rng.py``.

Regenerating the tables of the paper
------------------------------------

``regenerate_tables.sh`` maps each table to the commands that produce it:

.. code-block:: bash

   bash regenerate_tables.sh            # every table
   bash regenerate_tables.sh 2_4 18     # selected tables

Results land in ``tables/`` as ``table_NN_<name>.csv``. ``tables/README.md``
maps columns to tables and records the known issues. Tables 18 and 19 use the
hybrid solver and take several hours; the rest take about an hour in total.

Refitting the constants
-----------------------

The constants in ``src/const.py`` are fitted to the Lattice Estimator data in
``dataset/``. To regenerate all sixteen sets:

.. code-block:: bash

   bash find_all_constants.sh

They reproduce the published values of Equations (25), (27), (28), (31), (32),
(33), (34) and (35) to about seven decimal places.

Tests
-----

.. code-block:: bash

   pip install pytest
   python3 -m pytest                  # everything
   python3 -m pytest -m "not slow"    # skip refitting the constants

The suite pins values printed in the paper, so a change to a formula that
would move a published table fails the build. Tests marked ``needs_sage`` run
the command line and skip themselves when SageMath is absent; the rest cover
the formulas, the numerical solvers and the generators, and need only numpy
and scipy.
