A Tool for Fast and Secure LWE Parameter Selection
===================================================

We offer a tool to select secure parameters for LWE-based applications in a fast and flexible way. The tool can provide you with any of the following parameters: security level, size of the ciphertext modulus, LWE dimension, and standard deviation of the error distribution.

Our tool is constructed by studying the uSVP, BDD and Hybrid attacks against LWE. From this study, we derive formulas that describe each of the aforementioned parameters as a function of the others. You can find all the details in this paper: `A Tool for Fast and Secure LWE Parameter Selection: the FHE case <https://eprint.iacr.org/2024/1895>`_.

Basic Usage
----------------------------

We present the basic usage of the tool below. For more advanced usage, please refer to the `readthedocs <https://fastparameterselection.readthedocs.io/en/latest/>`_ section.

Find an estimation of the security level by running:

.. code-block:: bash

   python3 -m fastparameterselection.estimate --param "lambda" --n "1024" --logq "20;35;40" --secret "binary" --std "3.19"

.. code-block:: text

   secret dist.   | lwe dim. | log q | output
   ---------------+----------+-------+-------
   Uniform (-1 0) | 1024     | 20    | 173   
   Uniform (-1 0) | 1024     | 35    | 95    
   Uniform (-1 0) | 1024     | 40    | 83    

Find an estimation of the LWE dimension required to obtain a given security level:

.. code-block:: bash

   python3 -m fastparameterselection.estimate --param "n" --lambda "80" --logq "20-23" --secret "binary" --std "3.19"

.. code-block:: text

   secret dist.   | lambda | log q | output | pow
   ---------------+--------+-------+--------+----
   Uniform (-1 0) | 80     | 20    | 514    | 512
   Uniform (-1 0) | 80     | 21    | 538    | 512
   Uniform (-1 0) | 80     | 22    | 562    | 512
   Uniform (-1 0) | 80     | 23    | 586    | 512

Find an estimation of the size of the modulus q:

.. code-block:: bash

   python3 -m fastparameterselection.estimate --param "logq" --lambda "80" --n "1024" --secret "binary" --error "3.19"

.. code-block:: text

   secret dist.   | lambda | lwe dim. | output
   ---------------+--------+----------+-------
   Uniform (-1 0) | 80     | 1024     | 42    


Find an estimation of the standard deviation of the error distribution:

.. code-block:: bash

   python3 -m fastparameterselection.estimate --param "std_e" --lambda "192" --n "2048" --logq "64" --secret "binary"

.. code-block:: text

   secret dist.   | lambda | lwe dim. | log q | output           
   ---------------+--------+----------+-------+-------
   Uniform (-1 0) | 192    | 2048     | 64    | 28.60    


Common errors
----------------------------
Some MacOS users may encounter an error when running the tool using `python3 -m fastparameterselection.estimate`. This is due to the fact that the tool requires SageMath to run. To resolve this issue, you can run the tool using SageMath directly:
`sage-python3 -m fastparameterselection.estimate`



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

Common errors
-------------

Some MacOS users may encounter an error running the tool with ``python3``.
This is because the tool requires SageMath. Run it through Sage instead:

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

   command: [ "sage", "--python3", "-m", "fastparameterselection.estimate", "--param", "lambda",  "--n", "1024", "--logq", "20-30\\;35\\;40-60", "--secret", "binary", "--error", "gaussian", "--std", "3.19"]

Find an estimation of the security level and verify it against the Lattice Estimator:

.. code-block:: yaml

   command: [ "sage", "--python3", "-m", "fastparameterselection.estimate", "--param", "lambda",  "--n", "1024", "--logq", "20-30\\;35\\;40-60", "--secret", "binary", "--error", "gaussian", "--std", "3.19", "-v" ]

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

Reproducibility
---------------

Several numerical solvers (the hybrid attack, and ``--param std_e``) restart from randomised initial points. Each such call derives its own generator from a base seed and the call's own parameters, so a given set of LWE parameters always produces the same answer -- independently of how many other rows were computed in the same invocation. Pass ``--seed`` to change the base seed; see ``src/rng.py``.

The hybrid solvers are sensitive to their restart budget. At ``--param logq --lambda 192 --n 1024 --hw 64 --secret sparse`` the default budget can settle on a log q that the Lattice Estimator scores far below the target, while ``--nrestart 1000`` finds the right one. Run the hybrid with ``-v`` and check the reported ``est hybrid`` against your target.

The fitted constants in ``src/const.py`` can be regenerated from ``dataset/`` with:

.. code-block:: bash

   bash find_all_constants.sh

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

ToDo list
---------

* Resolve the BDD ``eta`` equation (see the OPEN ITEM note at the top of ``src/numerical_solver.py``).
* Include meet-in-the-middle for the hybrid attack. Challenge: derive compact formula for admissibility probability. Current status: added equations for optimization that include mitm speed-up for enumeration. Not tested.
* Improve runtime of numerical_lambda_hybrid(). Compute probability_enum() and babai_prob() on the log scale directly. 

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
