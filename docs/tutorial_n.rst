Estimating LWE Dimension
========================

Use ``--param "n"`` to estimate the LWE dimension (`n`) of the scheme.

The following commands estimate the LWE dimension for a binary secret distribution, error distribution as a discrete Gaussian with standard deviation of `3.19`, target security levels (`lambda`), and various sizes for the ciphertext modulus:

.. code-block:: bash

      python3 src/estimate.py --param "n" --lambda "128" --logq "27;37;45;54" --secret "binary" --error "gaussian" --std "3.19"

**Output Example**

.. code-block:: text

      secret dist.   | lambda | log q | output | pow 
      ---------------+--------+-------+--------+-----
      Uniform (-1 0) | 128    | 27    | 1063   | 1024
      Uniform (-1 0) | 128    | 37    | 1442   | 1024
      Uniform (-1 0) | 128    | 45    | 1746   | 2048
      Uniform (-1 0) | 128    | 54    | 2090   | 2048


**Show all estimations from formulas/numerical methods**

When adding the option ``--table``, the output shows the individual results for each of the formulas and numerical methods.

.. code-block:: text

      secret dist.   | lambda | log q | usvp | usvp_s | usvp num | bdd  | bdd_s | bdd num | output | pow 
      ---------------+--------+-------+------+--------+----------+------+-------+---------+--------+-----
      Uniform (-1 0) | 128    | 27    | 1045 | 1043   | 1032     | 1062 | 1063  | 1029    | 1063   | 1024
      Uniform (-1 0) | 128    | 37    | 1424 | 1425   | 1416     | 1442 | 1442  | 1411    | 1442   | 1024
      Uniform (-1 0) | 128    | 45    | 1727 | 1730   | 1723     | 1746 | 1746  | 1715    | 1746   | 2048
      Uniform (-1 0) | 128    | 54    | 2069 | 2072   | 2068     | 2088 | 2090  | 2057    | 2090   | 2048

**Compare results against the Lattice Estimator**

Use option ``-v``, to see the result (columns est 'name of attack / num') of running the Lattice Estimator with the given parameters and compare them with our formulas/numerical methods.

.. code-block:: text

      secret dist.   | lambda | log q | usvp | est usvp | usvp_s | est usvp_s | usvp num | est num | bdd  | est bdd | bdd_s | est bdd_s | bdd num | output | pow 
      ---------------+--------+-------+------+----------+--------+------------+----------+---------+------+---------+-------+-----------+---------+--------+-----
      Uniform (-1 0) | 128    | 27    | 1045 | 130      | 1043   | 130        | 1032     | 126     | 1062 | 130     | 1063  | 130       | 1029    | 1063   | 1024
      Uniform (-1 0) | 128    | 37    | 1424 | 129      | 1425   | 129        | 1416     | 126     | 1442 | 129     | 1442  | 129       | 1411    | 1442   | 1024
      Uniform (-1 0) | 128    | 45    | 1727 | 129      | 1730   | 129        | 1723     | 127     | 1746 | 129     | 1746  | 129       | 1715    | 1746   | 2048
      Uniform (-1 0) | 128    | 54    | 2069 | 129      | 2072   | 129        | 2068     | 127     | 2088 | 129     | 2090  | 129       | 2057    | 2090   | 2048

