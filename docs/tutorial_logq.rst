Estimating Modulus
==================

Use ``--param "logq"`` to estimate the largest ciphertext modulus that still
reaches a target security level.

.. code-block:: bash

   fastparams --param "logq" --lambda "100" --n "1024" --secret "binary" --error "gaussian" --std "3.19"

.. code-block:: text

      secret dist. | lambda | lwe dim. | output
      -------------+--------+----------+-------
      Binary       | 100    | 1024     | 33    

A larger modulus is weaker, so the binding constraint is the attack tolerating
the *smallest* modulus, and that is what the ``output`` column reports.

**Show all estimations from formulas/numerical methods**

.. code-block:: bash

   fastparams --param "logq" --lambda "100" --n "1024" --secret "binary" --error "gaussian" --std "3.19" --table

.. code-block:: text

      secret dist. | lambda | lwe dim. | logq usvp | logq bdd | output
      -------------+--------+----------+-----------+----------+-------
      Binary       | 100    | 1024     | 33        | 33       | 33    

**Compare results against the Lattice Estimator**

.. code-block:: bash

   fastparams --param "logq" --lambda "100" --n "1024" --secret "binary" --error "gaussian" --std "3.19" --table -v

.. code-block:: text

      secret dist. | lambda | lwe dim. | logq usvp | est usvp | logq bdd | est bdd | output
      -------------+--------+----------+-----------+----------+----------+---------+-------
      Binary       | 100    | 1024     | 33        | 103      | 33       | 101     | 33    

``est usvp`` and ``est bdd`` are each measured at their own attack's modulus,
not at ``output``. Read the row as two independent results plus a
recommendation, rather than as three views of one value.

.. seealso::

   :doc:`overview`, :doc:`tutorial_lambda`, :doc:`tutorial_n`,
   :doc:`tutorial_stde`, :doc:`tutorial_correction`, :doc:`limitations`.
