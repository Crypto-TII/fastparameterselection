Estimating Error Standard Deviation
===================================

Use ``--param "std_e"`` to estimate the smallest error standard deviation that
still reaches a target security level. Values are reported as
:math:`\log_2 \sigma_e`.

.. code-block:: bash

   python3 -m fastparameterselection.estimate --param "std_e" --lambda "192" --n "2048" --logq "64" --secret "binary"

.. code-block:: text

      secret dist. | lambda | lwe dim. | log q | output           
      -------------+--------+----------+-------+------------------
      Binary       | 192    | 2048     | 64    | 28.60            

A larger error is more secure, so the reported value is the largest across the
attacks.

**Show all estimations from formulas/numerical methods**

.. code-block:: bash

   python3 -m fastparameterselection.estimate --param "std_e" --lambda "192" --n "2048" --logq "64" --secret "binary" --table

.. code-block:: text

      secret dist. | lambda | lwe dim. | log q | log2(std_e) usvp  | log2(std_e) bdd   | output           
      -------------+--------+----------+-------+-------------------+-------------------+------------------
      Binary       | 192    | 2048     | 64    | 28.60             | 26.02             | 28.60            

**Compare results against the Lattice Estimator**

.. code-block:: bash

   python3 -m fastparameterselection.estimate --param "std_e" --lambda "192" --n "2048" --logq "64" --secret "binary" --table -v

.. code-block:: text

      secret dist. | lambda | lwe dim. | log q | log2(std_e) usvp  | est usvp | log2(std_e) bdd   | est bdd | output           
      -------------+--------+----------+-------+-------------------+----------+-------------------+---------+------------------
      Binary       | 192    | 2048     | 64    | 28.60             | 192      | 26.02             | 177     | 28.60            

Here the BDD solver's value of `26.02` measures at 177 bits against a target of
192, so it needs the correction described in :doc:`tutorial_correction`. Note
also that the Lattice Estimator loses precision when :math:`\sigma_e` is far
from `3.19`, so treat these ``est`` columns as indicative.

.. seealso::

   :doc:`overview`, :doc:`tutorial_lambda`, :doc:`tutorial_n`,
   :doc:`tutorial_logq`, :doc:`tutorial_correction`, :doc:`limitations`.
