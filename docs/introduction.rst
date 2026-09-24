Introduction
============

A tool for choosing secure parameters for LWE-based schemes, quickly.

Selecting parameters for Fully Homomorphic Encryption means balancing a
ciphertext modulus large enough to evaluate your circuit against a lattice
problem hard enough to be secure. The usual route is to run an LWE estimator
over candidate parameters until something fits, which is slow, or to pick from
a published table, which is rigid.

This tool takes a third route. Starting from the complexity of the best known
lattice attacks, it derives closed-form relations between the LWE parameters,
so any one of them can be computed directly from the others. See
:doc:`overview` for what that means in practice, and the paper for the
derivations: `A Tool for Fast and Secure LWE Parameter Selection: the FHE case
<https://eprint.iacr.org/2024/1895>`_.

.. code-block:: bash

   python3 src/estimate.py --param "lambda" --n "1024" --logq "20;35;40" --secret "binary" --std "3.19"

.. code-block:: text

   secret dist. | lwe dim. | log q | output
   -------------+----------+-------+-------
   Binary       | 1024     | 20    | 173
   Binary       | 1024     | 35    | 95
   Binary       | 1024     | 40    | 83

Read :doc:`limitations` before relying on a result.

Installation
------------



Common errors
-------------

Some MacOS users may encounter an error when running the tool using `python3 src/estimate.py`. This is due to the fact that the tool requires SageMath to run. To resolve this issue, you can run the tool using SageMath directly:

Use with Docker
---------------



Where to go next
----------------

- :doc:`overview` -- what the tool computes, and how the columns relate
- :doc:`options` -- every command line option
- :doc:`limitations` -- what the numbers do and do not cover
- :doc:`reproducibility` -- seeding, regenerating the paper's tables, tests

Citing
------

.. code-block:: bibtex

   @misc{cryptoeprint:2024/1895,
         author = {Beatrice Biasioli and Elena Kirshanova and Chiara Marcolla and Sergi Rovira},
         title = {A Tool for Fast and Secure {LWE} Parameter Selection: the {FHE} case},
         howpublished = {Cryptology {ePrint} Archive, Paper 2024/1895},
         year = {2024},
         url = {https://eprint.iacr.org/2024/1895}
   }

The paper associated with our tool is a follow-up and extension of the following paper presented at Africacrypt 2024. The pre-prints of both papers are available at:

- `Cryptology ePrint Archive, Report 2024/1001 <https://eprint.iacr.org/2024/1001>`_
- `Cryptology ePrint Archive, Report 2024/1895 <https://eprint.iacr.org/2024/1895>`_
