Estimating LWE Dimension
========================

Use ``--param "n"`` to estimate the LWE dimension `n`.

The following command estimates the dimension for a binary secret, a discrete
Gaussian error of standard deviation `3.19`, a target security level, and
various sizes of the ciphertext modulus:

.. code-block:: bash

   python3 src/estimate.py --param "n" --lambda "128" --logq "27;37;45;54" --secret "binary" --error "gaussian" --std "3.19"

.. code-block:: text

      secret dist. | lambda | log q | output | pow 
      -------------+--------+-------+--------+-----
      Binary       | 128    | 27    | 1063   | 1024
      Binary       | 128    | 37    | 1442   | 1024
      Binary       | 128    | 45    | 1746   | 2048
      Binary       | 128    | 54    | 2090   | 2048

The last column is the closest power of two to the recommended dimension.
Here the reported dimension is the largest across the attacks, since a larger
dimension is more secure.

**Show all estimations from formulas/numerical methods**

.. code-block:: bash

   python3 src/estimate.py --param "n" --lambda "128" --logq "27;37;45;54" --secret "binary" --error "gaussian" --std "3.19" --table

.. code-block:: text

      secret dist. | lambda | log q | usvp | usvp_s | usvp num | bdd  | bdd_s | bdd num | output | pow 
      -------------+--------+-------+------+--------+----------+------+-------+---------+--------+-----
      Binary       | 128    | 27    | 1046 | 1042   | 1032     | 1062 | 1063  | 1049    | 1063   | 1024
      Binary       | 128    | 37    | 1425 | 1424   | 1416     | 1442 | 1441  | 1431    | 1442   | 1024
      Binary       | 128    | 45    | 1728 | 1729   | 1723     | 1746 | 1746  | 1735    | 1746   | 2048
      Binary       | 128    | 54    | 2069 | 2072   | 2068     | 2088 | 2090  | 2077    | 2090   | 2048

**Compare results against the Lattice Estimator**

Each formula's dimension is fed back into the Lattice Estimator, so every
``est`` column is measured at the dimension in the column to its left.

.. code-block:: bash

   python3 src/estimate.py --param "n" --lambda "128" --logq "27;37;45;54" --secret "binary" --error "gaussian" --std "3.19" --table -v

.. code-block:: text

      secret dist. | lambda | log q | usvp | est usvp | usvp_s | est usvp_s | usvp num | est usvp num | bdd  | est bdd | bdd_s | est bdd_s | bdd num | est bdd num | output | pow 
      -------------+--------+-------+------+----------+--------+------------+----------+--------------+------+---------+-------+-----------+---------+-------------+--------+-----
      Binary       | 128    | 27    | 1046 | 130      | 1042   | 130        | 1032     | 128          | 1062 | 130     | 1063  | 130       | 1049    | 128         | 1063   | 1024
      Binary       | 128    | 37    | 1425 | 129      | 1424   | 129        | 1416     | 128          | 1442 | 129     | 1441  | 129       | 1431    | 128         | 1442   | 1024
      Binary       | 128    | 45    | 1728 | 129      | 1729   | 129        | 1723     | 129          | 1746 | 129     | 1746  | 129       | 1735    | 128         | 1746   | 2048
      Binary       | 128    | 54    | 2069 | 129      | 2072   | 129        | 2068     | 129          | 2088 | 129     | 2090  | 129       | 2077    | 128         | 2090   | 2048
