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

   python3 -m fastparameterselection.estimate --param "lambda" --n "1024" --logq "20;35;40" --secret "binary" --std "3.19"

.. code-block:: text

   secret dist. | lwe dim. | log q | output
   -------------+----------+-------+-------
   Binary       | 1024     | 20    | 173
   Binary       | 1024     | 35    | 95
   Binary       | 1024     | 40    | 83

Read :doc:`limitations` before relying on a result.

Installation
------------

SageMath is required and **cannot be installed by pip**: ``nd.py`` imports
``sage.all``, and the Lattice Estimator used by ``-v`` is a Sage library.
Install Sage first, by whichever route suits your system:

.. code-block:: bash

   sudo apt install sagemath          # Debian/Ubuntu
   conda install -c conda-forge sage  # conda

Then install this package into that interpreter:

.. code-block:: bash

   pip install git+https://github.com/Crypto-TII/fastparameterselection

which provides a ``fastparams`` command:

.. code-block:: bash

   fastparams --param "lambda" --n "1024" --logq "27" --secret "binary" --std "3.19"

From a checkout, without installing, the module form works too:

.. code-block:: bash

   python3 -m fastparameterselection.estimate --param "lambda" --n "1024" --logq "27" --secret "binary" --std "3.19"

Building the package needs ``setuptools >= 61`` for PEP 621 support; older
versions silently produce an empty ``UNKNOWN`` wheel.

The remaining Python dependencies (numpy, scipy, lmfit) are installed
automatically. The Lattice Estimator is vendored in the repository, so ``-v``
works from a checkout; after a pip install it is used only if an ``estimator``
package is importable.

On some MacOS setups ``python3`` cannot see Sage even when it is installed. Run
the tool through Sage itself in that case:

.. code-block:: bash

   sage --python3 -m fastparameterselection.estimate --param "lambda" --n "1024" --logq "27" --secret "binary"

Use with Docker
---------------

You can build and run the repository with Docker using the following command:

.. code-block:: bash

   docker-compose -f ./docker/docker-compose.yaml up --build

To only run the container, use:

.. code-block:: bash

   docker-compose -f ./docker/docker-compose.yaml up

Currently, it runs `estimate.py` to obtain the parameter ``lambda``, given ``n = 1024``, ``logq = 35``, binary secret distribution, and standard deviation of the error distribution ``3.19``. To run the estimation with your parameters, you can modify the command line in `docker-compose.yaml` as follows:

Find an estimation of the security level:

.. code-block:: yaml

   command: [ "sage", "--python3", "-m", "fastparameterselection.estimate", "--param", "lambda",  "--n", "1024", "--logq", "20-30;35;40-60", "--secret", "binary", "--error", "gaussian", "--std", "3.19"]

Find an estimation of the security level and verify it against the Lattice Estimator:

.. code-block:: yaml

   command: [ "sage", "--python3", "-m", "fastparameterselection.estimate", "--param", "lambda",  "--n", "1024", "--logq", "20-30;35;40-60", "--secret", "binary", "--error", "gaussian", "--std", "3.19", "-v" ]

Find an estimation of the LWE dimension:

.. code-block:: yaml

   command: ["sage", "--python3", "-m", "fastparameterselection.estimate", "--param", "n",  "--lambda", "80", "--logq", "20", "--secret", "binary", "--error", "gaussian", "--std", "3.19"]

Find an estimation of the size of the modulus q:

.. code-block:: yaml

   command: ["sage", "--python3", "-m", "fastparameterselection.estimate", "--param", "logq",  "--lambda", "80", "--n", "1024", "--secret", "binary", "--error", "gaussian", "--std", "3.19"]

Find an estimation of the standard deviation of the error distribution:

.. code-block:: yaml

   command: ["sage", "--python3", "-m", "fastparameterselection.estimate", "--param", "std_e",  "--lambda", "80", "--n", "1024", "--logq", "20", "--secret", "binary"]

**Note**: for a standard deviation of the error much larger than ``3.19``, the Lattice Estimator becomes imprecise, so the ``est bdd`` column of ``--param std_e`` should be treated as indicative only.

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
