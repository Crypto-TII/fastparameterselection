General LWE Parameters
======================

The formulas are fine-tuned for FHE settings, but the tool accepts arbitrary
LWE parameters :math:`(\lambda, n, \sigma_e, \sigma_s, q)`.

Both the secret and the error can be drawn from any of the distributions
below. Options without the ``s-`` prefix set the **error**; the ``s-`` prefixed
form sets the **secret**.

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Distribution
     - Secret
     - Error
   * - Uniform binary :math:`\mathcal{U}_2`
     - ``--secret "binary"``
     - ``--error "binary"``
   * - Uniform ternary :math:`\mathcal{U}_3`
     - ``--secret "ternary"``
     - ``--error "ternary"``
   * - Uniform mod p :math:`\mathcal{U}_p`
     - ``--secret "uniformmod" --q "5"``
     - ``--error "uniformmod" --q "5"``
   * - Uniform :math:`\mathcal{U}_{[a,b]}`
     - ``--secret "uniform" --s-a "-1" --s-b "1"``
     - ``--error "uniform" -a "-1" -b "1"``
   * - Sparse ternary :math:`\mathcal{HWT}(h)`
     - ``--secret "sparse" --hw "128"``
     - not applicable
   * - Discrete Gaussian :math:`\mathcal{DG}(0, \sigma^2)`
     - ``--secret "gaussian" --s-std "3.19"``
     - ``--error "gaussian" --std "3.19"``
   * - Centered binomial :math:`\psi_\eta`
     - ``--secret "binomial" --s-eta "3"``
     - ``--error "binomial" --eta "3"``

Outside the binary and ternary FHE settings the tool switches to
``--num-only`` automatically, since the fitted constants of Section 4 were
learned at :math:`\sigma_e = 3.19` for those two secret distributions only.

To evaluate Kyber-512, whose secret and error are both centered binomial with
:math:`\eta = 3`, and :math:`\log q = 12`:

.. code-block:: bash

   fastparams --param "lambda" --n "512" --logq "12" \
      --secret "binomial" --s-eta "3" --error "binomial" --eta "3" --num-only

.. code-block:: text

      secret dist. | lwe dim. | log q | output
      -------------+----------+-------+-------
      CB           | 512      | 12    | 138   

.. note::

   ``--eta`` and ``--s-eta`` must be given explicitly. They default to 1, and
   ``--std`` is ignored for a binomial distribution, so omitting them silently
   estimates a different scheme.

.. seealso::

   :doc:`overview`, :doc:`options`, :doc:`limitations`.
