Finding optimal values 
======================

Our tool offers a close approximation to the optimal values for the ciphertext modulus and the standard deviation of the error distribution. However,  we can use the output offered by the tool as a starting point for an exhaustive search of the optimal parameters using a LWE estimator.
In our tool, we use the Lattice Estimator. 

The following command estimate the error standard deviation for a binary secret distribution, target security levels (`lambda`), LWE dimensions, and various sizes for the ciphertext modulus:

.. code-block:: bash

      python3 src/estimate.py --param "std_e" --lambda "192" --n "2048" --logq "64" --secret "binary"

.. code-block:: text

      secret dist.   | lambda | lwe dim. | log q | log2(std_e) usvp   | est usvp | log2(std_e) bdd   | est bdd | output           
      ---------------+--------+----------+-------+--------------------+----------+-------------------+---------+------------------
      Uniform (-1 0) | 192    | 2048     | 64    | 25.71              | 177      | 28.60             | 190     | 28.60  

Observe that given the output of the numerical methods, we obtain a security level of `177` for the USVP attack and `190` for the BDD attack. We can get closer to the optimal values by using the Lattice Estimator. 

**Appling the correction**

By adding the option ``--correct`` or ``-c``, we can find the optimal values for the ciphertext modulus and the standard deviation of the error distribution. The output will show the corrected values.

.. code-block:: bash

      python3 src/estimate.py --param "std_e" --lambda "192" --n "2048" --logq "64" --secret "binary" --table -v -c

.. code-block:: text

      secret dist.   | lambda | lwe dim. | log q | log2(std_e) usvp   | est usvp | log2(std_e) bdd   | est bdd | output           
      ---------------+--------+----------+-------+--------------------+----------+-------------------+---------+------------------
      Uniform (-1 0) | 192    | 2048     | 64    | 26.21              | 195      | 29.10             | 194     | 29.10  