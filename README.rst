A Tool for Fast and Secure LWE Parameter Selection
===================================================

We offer a tool to select secure parameters for LWE-based applications in a fast and flexible way. The tool can provide you with any of the following parameters: security level, size of the ciphertext modulus, LWE dimension, and standard deviation of the error distribution.

Our tool is constructed by studying the uSVP, BDD and Hybrid attacks against LWE. From this study, we derive formulas that describe each of the aforementioned parameters as a function of the others. You can find all the details in this paper: `A Tool for Fast and Secure LWE Parameter Selection: the FHE case <https://eprint.iacr.org/2024/1895>`_.

Basic Usage
----------------------------

We present the basic usage of the tool below. For more advanced usage, please refer to the `readthedocs <https://fastparameterselection.readthedocs.io/en/latest/>`_ section.

Find an estimation of the security level by running:

.. code-block:: bash

   python3 src/estimate.py --param "lambda" --n "1024" --logq "20;35;40" --secret "binary" --std "3.19"

.. code-block:: text

   secret dist.   | lwe dim. | log q | output
   ---------------+----------+-------+-------
   Uniform (-1 0) | 1024     | 20    | 173   
   Uniform (-1 0) | 1024     | 35    | 95    
   Uniform (-1 0) | 1024     | 40    | 83    

Find an estimation of the LWE dimension required to obtain a given security level:

.. code-block:: bash

   python3 src/estimate.py --param "n" --lambda "80" --logq "20-23" --secret "binary" --std "3.19"

.. code-block:: text

   secret dist.   | lambda | log q | output | pow
   ---------------+--------+-------+--------+----
   Uniform (-1 0) | 80     | 20    | 514    | 512
   Uniform (-1 0) | 80     | 21    | 538    | 512
   Uniform (-1 0) | 80     | 22    | 562    | 512
   Uniform (-1 0) | 80     | 23    | 586    | 512

Find an estimation of the size of the modulus q:

.. code-block:: bash

   python3 src/estimate.py --param "logq" --lambda "80" --n "1024" --secret "binary" --error "3.19"

.. code-block:: text

   secret dist.   | lambda | lwe dim. | output
   ---------------+--------+----------+-------
   Uniform (-1 0) | 80     | 1024     | 42    


Find an estimation of the standard deviation of the error distribution:

.. code-block:: bash

   python3 src/estimate.py --param "std_e" --lambda "192" --n "2048" --logq "64" --secret "binary"

.. code-block:: text

   secret dist.   | lambda | lwe dim. | log q | output           
   ---------------+--------+----------+-------+-------
   Uniform (-1 0) | 192    | 2048     | 64    | 28.60    


Common errors
----------------------------
Some MacOS users may encounter an error when running the tool using `python3 src/estimate.py`. This is due to the fact that the tool requires SageMath to run. To resolve this issue, you can run the tool using SageMath directly:
`sage-python3 src/estimate.py`



Dependencies
------------

We have added the functionality to compare the output of our formulas against the `Lattice Estimator <https://github.com/malb/lattice-estimator>`_. We also use it for some of the more advanced features of our tool. Please use the custom version of the Estimator provided in this repository to avoid unexpected errors. Our code also depends on SageMath, you can follow this `link for an installation guide <https://doc.sagemath.org/html/en/installation/index.html>`_. 


The following Python libraries are required:
- `Numpy <https://numpy.org/>`_
- `SciPy <https://scipy.org/>`_

You can install the dependencies by running:

.. code-block:: bash

   pip install -r requirements.txt

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

   command: [ "sage", "--python3", "src/estimate.py", "--param", "lambda",  "--n", "1024", "--logq", "20-30\\;35\\;40-60", "--secret", "binary", "--error", "3.19"]

Find an estimation of the security level and verify it against the Lattice Estimator:

.. code-block:: yaml

   command: [ "sage", "--python3", "src/estimate.py", "--param", "lambda",  "--n", "1024", "--logq", "20-30\\;35\\;40-60", "--secret", "binary", "--error", "3.19", "--verify", "1" ]

Find an estimation of the LWE dimension:

.. code-block:: yaml

   command: ["sage", "--python3", "src/estimate.py", "--param", "n",  "--lambda", "80", "--logq", "20", "--secret", "binary", "--error", "3.19"]

Find an estimation of the size of the modulus q:

.. code-block:: yaml

   command: ["sage", "--python3", "src/estimate.py", "--param", "logq",  "--lambda", "80", "--n", "1024", "--secret", "binary", "--error", "3.19"]

Find an estimation of the standard deviation of the error distribution:

.. code-block:: yaml

   command: ["sage", "--python3", "src/estimate.py", "--param", "std_e",  "--lambda", "80", "--n", "1024", "--logq", "20", "--secret", "binary"]

Find an estimation of the security level, given the example parameters in `example_lambda_binary.csv`:

.. code-block:: yaml

   command: ["sage", "--python3", "src/estimate.py", "--param", "lambda", "--file", "./examples/example_lambda_binary.csv", "--verify", "1"]

Find an estimation of the LWE dimension, given the example parameters in `example_n_ternary.csv`:

.. code-block:: yaml

   command: ["sage", "--python3", "src/estimate.py", "--param", "n", "--file", "./examples/example_n_ternary.csv"]

Find an estimation of the size of the modulus q, given the example parameters in `example_logq_binary.csv`:

.. code-block:: yaml

   command: ["sage", "--python3", "src/estimate.py", "--param", "logq", "--file", "./examples/example_logq_binary.csv"]

Find an estimation of the standard deviation of the error distribution, given the example parameters in `example_error_binary.csv`:

.. code-block:: yaml

   command: ["sage", "--python3", "src/estimate.py", "--param", "error", "--file", "./examples/example_error_binary.csv"]

**Note**: In the Docker version, we applied a change to the `Lattice Estimator <https://github.com/malb/lattice-estimator>`_ to address an imprecision in the case where the standard deviation of the error distribution is much larger than ``3.19``.

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
