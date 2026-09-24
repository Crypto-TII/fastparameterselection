Estimating Security Level
==============================

Use ``--param "lambda"`` to estimate the security level of the LWE scheme.

The following command estimates the security level for a binary secret, a
discrete Gaussian error of standard deviation `3.19`, LWE dimension `1024`,
and various sizes of the ciphertext modulus:

.. code-block:: bash

   python3 src/estimate.py --param "lambda" --n "1024" --logq "20;24-28;30;33;37;42" --secret "binary" --error "gaussian" --std "3.19"

.. code-block:: text

      secret dist. | lwe dim. | log q | output
      -------------+----------+-------+-------
      Binary       | 1024     | 20    | 173   
      Binary       | 1024     | 24    | 142   
      Binary       | 1024     | 25    | 136   
      Binary       | 1024     | 26    | 130   
      Binary       | 1024     | 27    | 125   
      Binary       | 1024     | 28    | 121   
      Binary       | 1024     | 30    | 112   
      Binary       | 1024     | 33    | 101   
      Binary       | 1024     | 37    | 90    
      Binary       | 1024     | 42    | 79    

The reported level is the lowest across the attacks considered.

**Show all estimations from formulas/numerical methods**

Adding ``--table`` shows the individual result of each formula and numerical
method: ``usvp`` and ``bdd`` are the fitted closed forms, Eq. (25) and
Eq. (28); ``usvp_s`` and ``bdd_s`` their simplified versions, Eq. (27) and
Eq. (31); ``usvp num`` and ``bdd num`` the numerical solvers of Section 5.

.. code-block:: bash

   python3 src/estimate.py --param "lambda" --n "1024" --logq "20;24-28;30;33;37;42" --secret "binary" --error "gaussian" --std "3.19" --table

.. code-block:: text

      secret dist. | lwe dim. | log q | usvp | usvp_s | usvp num | bdd | bdd_s | bdd num | output
      -------------+----------+-------+------+--------+----------+-----+-------+---------+-------
      Binary       | 1024     | 20    | 178  | 177    | 174      | 171 | 173   | 171     | 173   
      Binary       | 1024     | 24    | 145  | 145    | 144      | 142 | 142   | 141     | 142   
      Binary       | 1024     | 25    | 139  | 139    | 137      | 136 | 136   | 135     | 136   
      Binary       | 1024     | 26    | 133  | 133    | 132      | 130 | 130   | 129     | 130   
      Binary       | 1024     | 27    | 128  | 128    | 126      | 125 | 125   | 124     | 125   
      Binary       | 1024     | 28    | 123  | 123    | 122      | 121 | 120   | 119     | 121   
      Binary       | 1024     | 30    | 114  | 114    | 113      | 112 | 112   | 111     | 112   
      Binary       | 1024     | 33    | 103  | 103    | 102      | 101 | 101   | 100     | 101   
      Binary       | 1024     | 37    | 92   | 92     | 90       | 90  | 90    | 88      | 90    
      Binary       | 1024     | 42    | 81   | 81     | 78       | 79  | 79    | 77      | 79    

**Compare results against the Lattice Estimator**

Add ``-v`` to run the Lattice Estimator on the same parameters. The ``est``
columns are its verdict, as ``floor(log2(rop))`` under the BDGL16 cost model.

.. code-block:: bash

   python3 src/estimate.py --param "lambda" --n "1024" --logq "20;24-28;30;33;37;42" --secret "binary" --error "gaussian" --std "3.19" --table -v

.. code-block:: text

      secret dist. | lwe dim. | log q | usvp | usvp_s | usvp num | est usvp | bdd | bdd_s | bdd num | est bdd | output
      -------------+----------+-------+------+--------+----------+----------+-----+-------+---------+---------+-------
      Binary       | 1024     | 20    | 178  | 177    | 174      | 174      | 171 | 173   | 171     | 171     | 173   
      Binary       | 1024     | 24    | 145  | 145    | 144      | 144      | 142 | 142   | 141     | 141     | 142   
      Binary       | 1024     | 25    | 139  | 139    | 137      | 138      | 136 | 136   | 135     | 135     | 136   
      Binary       | 1024     | 26    | 133  | 133    | 132      | 132      | 130 | 130   | 129     | 130     | 130   
      Binary       | 1024     | 27    | 128  | 128    | 126      | 127      | 125 | 125   | 124     | 125     | 125   
      Binary       | 1024     | 28    | 123  | 123    | 122      | 122      | 121 | 120   | 119     | 120     | 121   
      Binary       | 1024     | 30    | 114  | 114    | 113      | 114      | 112 | 112   | 111     | 112     | 112   
      Binary       | 1024     | 33    | 103  | 103    | 102      | 103      | 101 | 101   | 100     | 101     | 101   
      Binary       | 1024     | 37    | 92   | 92     | 90       | 91       | 90  | 90    | 88      | 90      | 90    
      Binary       | 1024     | 42    | 81   | 81     | 78       | 80       | 79  | 79    | 77      | 79      | 79    
