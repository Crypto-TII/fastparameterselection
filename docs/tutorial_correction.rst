Finding optimal values
======================

The tool gives a close approximation to the largest secure modulus and the
smallest secure error. Its output can then be used as the starting point for
an exhaustive search with an LWE estimator; here we use the Lattice Estimator.

Take the error standard deviation for a binary secret at
:math:`\lambda = 192`, `n = 2048`, :math:`\log q = 64`:

.. code-block:: bash

   python3 -m fastparameterselection.estimate --param "std_e" --lambda "192" --n "2048" --logq "64" --secret "binary" --table -v

.. code-block:: text

      secret dist. | lambda | lwe dim. | log q | log2(std_e) usvp  | est usvp | log2(std_e) bdd   | est bdd | output           
      -------------+--------+----------+-------+-------------------+----------+-------------------+---------+------------------
      Binary       | 192    | 2048     | 64    | 28.60             | 192      | 26.02             | 177     | 28.60            

The uSVP solver's value already reaches the target, but the BDD solver's
`26.02` measures at only 177 bits against a target of 192.

**Applying the correction**

Adding ``-c`` walks each value towards the target in steps of 0.1, calling the
Lattice Estimator at each step, and reports both the corrected value and how
many calls it took.

.. code-block:: bash

   python3 -m fastparameterselection.estimate --param "std_e" --lambda "192" --n "2048" --logq "64" --secret "binary" --table -v -c

.. code-block:: text

      secret dist. | lambda | lwe dim. | log q | log2(std_e) usvp  | est usvp | * log2(std_e) usvp | * est usvp | est calls usvp | log2(std_e) bdd   | est bdd | * log2(std_e) bdd | * est bdd | est calls bdd | output           
      -------------+--------+----------+-------+-------------------+----------+--------------------+------------+----------------+-------------------+---------+-------------------+-----------+---------------+------------------
      Binary       | 192    | 2048     | 64    | 28.60             | 192      | 28.50              | 192        | 3              | 26.02             | 177     | 28.92             | 192       | 30            | 28.92            

Columns prefixed ``*`` are the corrected values. The BDD value moves from
`26.02` to `28.92`, which measures at 192, and took 30 estimator calls to
find. The correction is bounded at ``MAX_CORRECTION_CALLS`` in
``src/const.py``.

The same option works for ``--param "logq"``.

.. seealso::

   :doc:`tutorial_logq`, :doc:`tutorial_stde`, :doc:`limitations`.
