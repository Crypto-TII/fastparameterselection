A Tool for Fast and Secure LWE Parameter Selection
===================================================

We offer a tool to select secure parameters for LWE-based applications in a fast and flexible way. The tool can provide you with any of the following parameters: security level, size of the ciphertext modulus, LWE dimension, and standard deviation of the error distribution.

Our tool is constructed by studying the uSVP, BDD and Hybrid attacks against LWE. From this study, we derive formulas that describe each of the aforementioned parameters as a function of the others. You can find all the details in this paper: `A Tool for Fast and Secure LWE Parameter Selection: the FHE case <https://eprint.iacr.org/2024/1895>`_.

Options
----------------------------

Our tool supports the following command-line options:

General Options
---------------

- **-h, --help**:
  Display the help message and exit.

- **--param PARAM**:
  Specify the parameter to estimate. Options include:
   - **"lambda"**: Security level.
   - **"n"**: LWE dimension.
   - **"logq"**: Modulus size.
   - **"std_e"**: Standard deviation of the error distribution.

  **--logq LOGQ**:
  Specify the modulus size (e.g., "32;64" or "20-30").

- **--lambda LAMBDA**:
  Specify the target security level (e.g., "80", "128").

- **--n N**:
  Specify the LWE dimension (e.g., "1024", "2048").

Secret and Error Distribution Options
----------------------------

- **--secret SECRET and --error ERROR** 
  Specify the secret distribution. Options include:

   - **"binary"**: Uniform binary distribution.
   - **"ternary"**: Uniform ternary distribution.
   - **"uniformmod"**: Uniform modulus distribution (requires **--p**).
   - **"uniform"**: Uniform distribution (requires **--a** and **--b**).
   - **"sparse"**: Sparse ternary distribution (requires **--hw**).
   - **"gaussian"**: Discrete Gaussian distribution (requires **--s-std**).
   - **"binomial"**: Centered binomial distribution (requires **--s-eta**).


- **--s-std S_STD**:
  Standard deviation for the secret distribution (used with "gaussian").

- **--s-a S_A**:
  Lower bound for the secret distribution (used with "uniform").

- **--s-b S_B**:
  Upper bound for the secret distribution (used with "uniform").

- **--hw HW**:
  Hamming weight for the secret distribution (used with "sparse").

- **--s-eta S_ETA**:
  Parameter for the centered binomial distribution (used with "binomial").

- **--std STD**:
  Standard deviation for the error distribution (used with "gaussian").

- **--a A**:
  Lower bound for the error distribution (used with "uniform").

- **--b B**:
  Upper bound for the error distribution (used with "uniform").

- **--eta ETA**:
  Parameter for the centered binomial distribution (used with "binomial").

Additional Options
-------------------

- **--table**:
  Output results in a table format.

- **--num-only**:
  Output only numerical results (no headers).

- **-v**:
  Compare results against the Lattice Estimator.

- **--ntru**:
  Use NTRU-specific parameter estimation.

- **-c, --correction**:
  Apply correction factors to the estimation.

Table Description
----------------------------

The output table may contain any of the following columns:

- **secret dist.**: The distribution of the secret (can be binary, ternary, or sparse).
- **lwe dim.**: The Learning With Errors (LWE) dimension.
- **lambda**: The security level.
- **log q**: The size of the modulus q in bits.
- **lwe est**: The output of running the Lattice Estimator using the output of our formulas and the rest of the LWE parameters.
- **usvp**: Output of the formula which estimates the cost of the (unique) SVP attack.
- **usvp_s**: Output of the simplified formula (removing dependency on beta) which estimates the cost of the (unique) SVP attack.
- **bdd**: Output of the formula which estimates the cost of the BDD attack.
- **bdd_s**: Output of the simplified formula (removing dependency on beta) which estimates the cost of the BDD attack.
- **logq usvp**: Output of the numerical approximation of log q for the (unique) SVP attack.
- **logq bdd**: Output of the numerical approximation of log q for the BDD attack.
- **usvp num**: Output of the numerical approximation of the (unique) SVP attack.
- **bdd num**: Output of the numerical approximation of the BDD attack.
- **log2(std_e) usvp**: Output of the numerical approximation of the (log2) standard deviation of the error for the (unique) SVP attack.
- **log2(std_e) bdd**: Output of the numerical approximation of the (log2) standard deviation of the error for the BDD attack.
- **bdd 3.19**: The result of running the Lattice Estimator with standard deviation of the error 3.19 and primal_bdd.
- **usvp 3.19**: The result of running the Lattice Estimator with standard deviation of the error 3.19 and primal_usvp.
- **est usvp**: Output of the Lattice Estimator for the (unique) SVP attack.
- **est bdd**: Output of the Lattice Estimator for the BDD attack.
- **est usvp_s**: Output of the Lattice Estimator using the result from the simplified formula for the (unique) SVP attack.
- **est bdd_s**: Output of the Lattice Estimator using the result from the simplified formula for the BDD attack.
- **output**: Recommended value to be used considering all the outputs of the formulas and numerical methods.
- **pow**: Closest power of 2 to the LWE dimension recommended in Output.
- **hw**: Hamming weight of the secret.
- **hybrid**: Output of the numerical approximation for lambda of the hybrid attack.
- **logq hybrid**: Output of the numerical approximation for logq of the hybrid attack.
- **est hybrid**: Output of the Lattice Estimator for the hybrid attack.
- **est**: Output of the Lattice Estimator.


Basic Usage
----------------------------

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

We have added the functionality to compare the output of our formulas against the `Lattice Estimator <https://github.com/malb/lattice-estimator>`_. Please download the Estimator if you want to use such functionality.

**Note**: At present, the Estimator is also needed to run one of the formulas. This will be fixed shortly.

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