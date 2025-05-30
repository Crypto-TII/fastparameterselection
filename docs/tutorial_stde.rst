Estimating Error Standard Deviation
===================================

Use ``--param "std_e"`` to estimate the error standard deviation (`std_e`) of the scheme.

The following commands estimate the error standard deviation for a binary secret distribution, target security levels (`lambda`), LWE dimensions, and various sizes for the ciphertext modulus:

.. code-block:: bash

      python3 src/estimate.py --param "std_e" --lambda "192" --n "2048" --logq "64" --secret "binary"

**Output Example**

.. code-block:: text

      secret dist.   | lambda | lwe dim. | log q | output           
      ---------------+--------+----------+-------+-------
      Uniform (-1 0) | 192    | 2048     | 64    | 28.60 

**Show all estimations from formulas/numerical methods**

When adding the option ``--table``, the output shows the individual results for each of the formulas and numerical methods.

.. code-block:: text

      secret dist.   | lambda | lwe dim. | log q | log2(std_e) usvp  | log2(std_e) bdd    | output           
      ---------------+--------+----------+-------+-------------------+--------------------+-------
      Uniform (-1 0) | 192    | 2048     | 64    | 28.60             | 25.71              | 28.60   

**Compare results against the Lattice Estimator**

Use option ``-v``, to see the result (columns est 'name of attack / num') of running the Lattice Estimator with the given parameters and compare them with our formulas/numerical methods.

.. code-block:: text

      secret dist.   | lambda | lwe dim. | log q | log2(std_e) usvp   | est usvp | log2(std_e) bdd   | est bdd | output           
      ---------------+--------+----------+-------+--------------------+----------+-------------------+---------+-------
      Uniform (-1 0) | 192    | 2048     | 64    | 25.71              | 177      | 28.60             | 190     | 28.60 