Estimating Modulus
==================

Use ``--param "logq"`` to estimate the modulus (`logq`) of the scheme.

The following commands estimate the modulus for a binary secret distribution, error distribution as a discrete Gaussian with standard deviation of `3.19`, target security levels (`lambda`), and LWE dimensions:

.. code-block:: bash

      python3 src/estimate.py --param "logq" --lambda "100" --n "1024" --secret "binary" --error "gaussian" --std "3.19"

**Output Example**

.. code-block:: text

      secret dist.   | lambda | lwe dim. | output
      ---------------+--------+----------+-------
      Uniform (-1 0) | 100    | 1024     | 33    

**Show all estimations from formulas/numerical methods**

When adding the option ``--table``, the output shows the individual results for each of the formulas and numerical methods.

.. code-block:: text

      secret dist.   | lambda | lwe dim. | logq usvp | logq bdd | output
      ---------------+--------+----------+-----------+----------+-------
      Uniform (-1 0) | 100    | 1024     | 33        | 33       | 33    

**Compare results against the Lattice Estimator**

Use option ``-v``, to see the result (columns est 'name of attack / num') of running the Lattice Estimator with the given parameters and compare them with our formulas/numerical methods.

.. code-block:: text

      secret dist.   | lambda | lwe dim. | logq usvp | est usvp | logq bdd | est bdd | output
      ---------------+--------+----------+-----------+----------+----------+---------+-------
      Uniform (-1 0) | 100    | 1024     | 33        | 103      | 33       | 101     | 33  