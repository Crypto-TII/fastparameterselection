A Tool for Fast and Secure LWE Parameter Selection
===================================================

Select secure parameters for LWE-based schemes without running a lattice
estimator over every candidate. The security level, the LWE dimension, the
ciphertext modulus and the error standard deviation are mutually constrained:
fix any three and this tool solves for the fourth.

It is built on an analysis of the uSVP, BDD and hybrid attacks, from which
closed-form relations between the parameters are derived. The details are in
the paper: `A Tool for Fast and Secure LWE Parameter Selection: the FHE case
<https://eprint.iacr.org/2024/1895>`_.

Full documentation, including what the tool does **not** cover, is at
`readthedocs <https://fastparameterselection.readthedocs.io/en/latest/>`_.

Installation
------------

SageMath is required and **cannot be installed by pip**: ``nd.py`` imports
``sage.all``, and the Lattice Estimator used by ``-v`` is a Sage library.
Install Sage first, by whichever route suits your system:

.. code-block:: bash

   sudo apt install sagemath          # Debian/Ubuntu
   conda install -c conda-forge sage  # conda

Then install this package into that interpreter. Install **from a checkout, in
editable mode**, so that ``-v`` keeps working:

.. code-block:: bash

   git clone https://github.com/Crypto-TII/fastparameterselection
   cd fastparameterselection
   pip install -e .

This provides a ``fastparams`` command that runs from any directory:

.. code-block:: bash

   fastparams --param "lambda" --n "1024" --logq "27" --secret "binary" --std "3.19"

Editable mode matters for ``-v``. The Lattice Estimator is vendored in the
repository and is located relative to the package, so an editable install finds
it, while a copying install (``pip install .`` or
``pip install git+https://...``) does not. With a copying install ``-v`` prints
``Failed to import lattice_estimator`` and the verification columns are skipped;
everything else still works.

Without installing at all, the module form works from a checkout:

.. code-block:: bash

   python3 -m fastparameterselection.estimate --param "lambda" --n "1024" --logq "27" --secret "binary" --std "3.19"

Build tooling
~~~~~~~~~~~~~

Installing needs ``setuptools >= 61`` for PEP 621 metadata, and ``>= 64`` for
the editable install above. Older versions fail in two confusing ways: an
editable install reports a missing ``build_editable`` hook, and a regular
install silently produces an empty wheel named ``UNKNOWN-0.0.0`` and exits
successfully.

On distributions that ship an old pip (Ubuntu 22.04 ships pip 22.0.2 and
setuptools 59.6), the system setuptools leaks into pip's isolated build
environment, so the build fails this way even when a newer setuptools is
available. Upgrade both first:

.. code-block:: bash

   pip install --user -U pip setuptools

The remaining Python dependencies (numpy, scipy, lmfit) are installed
automatically.

On some MacOS setups ``python3`` cannot see Sage even when it is installed. Run
the tool through Sage itself in that case:

.. code-block:: bash

   sage --fastparams --param "lambda" --n "1024" --logq "27" --secret "binary"

Basic usage
-----------

Estimate the security level of a parameter set:

.. code-block:: bash

   fastparams --param "lambda" --n "1024" --logq "20;35;40" --secret "binary" --std "3.19"

.. code-block:: text

   secret dist. | lwe dim. | log q | output
   -------------+----------+-------+-------
   Binary       | 1024     | 20    | 173   
   Binary       | 1024     | 35    | 95    
   Binary       | 1024     | 40    | 83    

Estimate the LWE dimension needed for a target security level:

.. code-block:: bash

   fastparams --param "n" --lambda "80" --logq "20-23" --secret "binary" --std "3.19"

.. code-block:: text

   secret dist. | lambda | log q | output | pow
   -------------+--------+-------+--------+----
   Binary       | 80     | 20    | 513    | 512
   Binary       | 80     | 21    | 537    | 512
   Binary       | 80     | 22    | 561    | 512
   Binary       | 80     | 23    | 586    | 512

Estimate the largest usable modulus:

.. code-block:: bash

   fastparams --param "logq" --lambda "80" --n "1024" --secret "binary" --error "gaussian" --std "3.19"

.. code-block:: text

   secret dist. | lambda | lwe dim. | output
   -------------+--------+----------+-------
   Binary       | 80     | 1024     | 41    

Estimate the smallest usable error:

.. code-block:: bash

   fastparams --param "std_e" --lambda "192" --n "2048" --logq "64" --secret "binary"

.. code-block:: text

   secret dist. | lambda | lwe dim. | log q | output           
   -------------+--------+----------+-------+------------------
   Binary       | 192    | 2048     | 64    | 28.60            

Add ``--table`` to see each formula and solver separately, and ``-v`` to check
the result against the Lattice Estimator.

Before you rely on a result
---------------------------

The tool answers the security question only, and not every part of it. It does
not model dual or quantum attacks, says nothing about whether a parameter set
supports your circuit, and its fitted constants were learned at a standard
deviation of 3.19 for binary and ternary secrets. The hybrid ``--param logq``
search can return a modulus below the target, and ``--ntru`` is known to be
broken.

The full list is on the `limitations
<https://fastparameterselection.readthedocs.io/en/latest/limitations.html>`_
page. Confirm candidate parameters with ``-v`` before deploying them.

Reproducibility
---------------

Several numerical solvers (the hybrid attack, and ``--param std_e``) restart from randomised initial points. Each such call derives its own generator from a base seed and the call's own parameters, so a given set of LWE parameters always produces the same answer -- independently of how many other rows were computed in the same invocation. Pass ``--seed`` to change the base seed; see ``fastparameterselection/rng.py``.

The hybrid solvers are sensitive to their restart budget. At ``--param logq --lambda 192 --n 1024 --hw 64 --secret sparse`` the default budget can settle on a log q that the Lattice Estimator scores far below the target, while ``--nrestart 1000`` finds the right one. Run the hybrid with ``-v`` and check the reported ``est hybrid`` against your target.

The fitted constants in ``fastparameterselection/const.py`` can be regenerated from ``dataset/`` with:

.. code-block:: bash

   bash find_all_constants.sh

The tables of the paper can be regenerated with:

.. code-block:: bash

   bash regenerate_tables.sh

Results land in ``tables/``; ``tables/README.md`` maps columns to tables.

Tests
-----

.. code-block:: bash

   pip install pytest
   python3 -m pytest              # everything
   python3 -m pytest -m "not slow"  # skip refitting the constants

The suite pins values printed in the paper, so a change to a formula that
would move a published table fails the build. Tests marked ``needs_sage`` run
the command line and skip themselves when SageMath is not installed; the rest
cover the formulas, the numerical solvers and the per-call generators and need
only numpy and scipy.

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

ToDo list
---------

* Include meet-in-the-middle for the hybrid attack. Challenge: derive a compact formula for the admissibility probability. Current status: equations for the optimisation including the mitm speed-up for enumeration are in place, untested.
* Improve the runtime of ``numerical_lambda_hybrid()``. Computing ``probability_enum()`` and ``ss_enum()`` on the log scale is both faster and removes their dependency on Sage; a prototype agreed with the current code to 1e-10 and ran up to 750x faster.
* Make the hybrid ``--param logq`` search reject a modulus its own lambda model does not vouch for, instead of returning it with a warning.
* Fix ``--ntru``: ``check_overstreched()`` reports every parameter set as overstretched.

Bugs
----

Please report bugs through the `GitHub issue tracker <https://github.com/Crypto-TII/fastparameterselection/issues>`_.

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
