Estimating Security Level
==============================

Use ``--param "lambda"`` to estimate the security level of the LWE scheme.

The following commands estimate the security level for a binary secret distribution, error distribution a discrete Gaussian with standard deviation of `3.19`, lwe dimension `1024`, and various sizes for the ciphertext modulus:

.. code-block:: bash

   python3 src/estimate.py --param "lambda" --n "1024" --logq "20;24-28;30;33;37;42" --secret "binary" --error "gaussian" --std "3.19"

.. code-block:: text

      secret dist.   | lwe dim. | log q | output
      ---------------+----------+-------+-------
      Uniform (-1 0) | 1024     | 20    | 173   
      Uniform (-1 0) | 1024     | 24    | 142   
      Uniform (-1 0) | 1024     | 25    | 136   
      Uniform (-1 0) | 1024     | 26    | 130   
      Uniform (-1 0) | 1024     | 27    | 125   
      Uniform (-1 0) | 1024     | 28    | 120   
      Uniform (-1 0) | 1024     | 30    | 112   
      Uniform (-1 0) | 1024     | 33    | 101   
      Uniform (-1 0) | 1024     | 37    | 90    
      Uniform (-1 0) | 1024     | 42    | 79 

**Show all estimations from formulas/numerical methods**

When adding the option ``--table``, the output shows the individual results for each of the formulas and numerical methods.

.. code-block:: text

      secret dist.   | lwe dim. | log q | usvp | usvp_s | usvp num | bdd | bdd_s | bdd num | output
      ---------------+----------+-------+------+--------+----------+-----+-------+---------+-------
      Uniform (-1 0) | 1024     | 20    | 178  | 172    | 174      | 166 | 173   | 175     | 173   
      Uniform (-1 0) | 1024     | 24    | 145  | 142    | 144      | 138 | 142   | 144     | 142   
      Uniform (-1 0) | 1024     | 25    | 139  | 136    | 137      | 132 | 136   | 138     | 136   
      Uniform (-1 0) | 1024     | 26    | 133  | 130    | 132      | 127 | 130   | 132     | 130   
      Uniform (-1 0) | 1024     | 27    | 128  | 125    | 126      | 122 | 125   | 127     | 125   
      Uniform (-1 0) | 1024     | 28    | 123  | 120    | 122      | 117 | 120   | 122     | 120   
      Uniform (-1 0) | 1024     | 30    | 114  | 112    | 113      | 109 | 112   | 113     | 112   
      Uniform (-1 0) | 1024     | 33    | 103  | 101    | 102      | 99  | 101   | 102     | 101   
      Uniform (-1 0) | 1024     | 37    | 92   | 90     | 90       | 88  | 90    | 90      | 90    
      Uniform (-1 0) | 1024     | 42    | 81   | 80     | 78       | 77  | 79    | 79      | 79 

**Compare results against the Lattice Estimator**

Use option ``-v``, to see the result (columns est 'name of attack / num') of running the Lattice Estimator with the given parameters and compare them with our formulas/numerial methods.

.. code-block:: text

      secret dist.   | lwe dim. | log q | usvp | usvp_s | usvp num | est usvp | bdd | bdd_s | bdd num | est bdd | output
      ---------------+----------+-------+------+--------+----------+----------+-----+-------+---------+---------+-------
      Uniform (-1 0) | 1024     | 20    | 178  | 172    | 174      | 174      | 166 | 173   | 175     | 172     | 173   
      Uniform (-1 0) | 1024     | 24    | 145  | 142    | 144      | 144      | 138 | 142   | 144     | 142     | 142   
      Uniform (-1 0) | 1024     | 25    | 139  | 136    | 137      | 138      | 132 | 136   | 138     | 136     | 136   
      Uniform (-1 0) | 1024     | 26    | 133  | 130    | 132      | 132      | 127 | 130   | 132     | 130     | 130   
      Uniform (-1 0) | 1024     | 27    | 128  | 125    | 126      | 127      | 122 | 125   | 127     | 125     | 125   
      Uniform (-1 0) | 1024     | 28    | 123  | 120    | 122      | 122      | 117 | 120   | 122     | 120     | 120   
      Uniform (-1 0) | 1024     | 30    | 114  | 112    | 113      | 114      | 109 | 112   | 113     | 112     | 112   
      Uniform (-1 0) | 1024     | 33    | 103  | 101    | 102      | 103      | 99  | 101   | 102     | 101     | 101   
      Uniform (-1 0) | 1024     | 37    | 92   | 90     | 90       | 91       | 88  | 90    | 90      | 90      | 90    
      Uniform (-1 0) | 1024     | 42    | 81   | 80     | 78       | 80       | 77  | 79    | 79      | 79      | 79    
