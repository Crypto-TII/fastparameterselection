General LWE Parameters
======================

While our formulas are fine-tuned for FHE settings, our tool allows users to input arbitrary LWE parameters :math:`(\lambda, n, \sigma_e, \sigma_s, q)`.

In particular, :math:`\sigma_s` and :math:`\sigma_e` can be chosen from the following distributions:

- **Uniform binary distribution** :math:`\mathcal{U}_2`:
  Use the command ``--secret "binary"``.

- **Uniform ternary distribution** :math:`\mathcal{U}_3`:
  Use the command ``--secret "ternary"``.

- **Uniform modulus distribution** :math:`\mathcal{U}_p`:
  Use the command ``--secret "uniformmod" --p "5"``.

- **Uniform distribution** :math:`\mathcal{U}_{[a,b]}`:
  Use the command ``--secret "uniform" --a "-1" --b "1"``.

- **Sparse ternary distribution** :math:`\mathcal{HWT}(h)`:
  Use the command ``--secret "sparse" --hw "128"``, where ``hw`` is the Hamming weight :math:`h`.

- **Discrete Gaussian distribution** :math:`\mathcal{D}_{\mathbb{Z}, \sigma^2}`:
  Use the command ``--secret "gaussian" --sig "3.19"``.

- **Centered binomial distribution** :math:`\psi_\eta`:
  Use the command ``--secret "binomial" --eta "21"``.

In this case, we recommend using the option ``--num-only``, which disables the fine-tuned formulas and runs only our numerical solver.

To evaluate the security of the Kyber scheme with :math:`n = 512` and :math:`\log q = 11`, use the following command:

.. code-block:: bash

   python3 estimate.py --param "lambda" --n "512" --logq "11" \
      --secret "binomial" --error "binomial" --std "1.22" --num-only

.. code-block:: text

      secret dist.     | lwe dim. | log q | usvp num | est usvp | bdd num | est bdd | output
      -----------------+----------+-------+----------+----------+---------+---------+-------
      CenteredBinomial | 512      | 11    | 136      | 137      | 136     | 131     | 136   