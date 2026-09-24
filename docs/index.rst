Fast Parameter Selection
========================

Choose secure LWE parameters for FHE without running a lattice estimator over
every candidate. Fix any three of the security level, the dimension, the
modulus and the error, and the tool solves for the fourth.

.. code-block:: bash

   python3 -m fastparameterselection.estimate --param "n" --lambda "128" --logq "27" --secret "binary" --std "3.19"

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   introduction
   overview
   options
   limitations

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorial_lambda
   tutorial_n
   tutorial_logq
   tutorial_stde
   tutorial_correction
   tutorial_non_fhe

.. toctree::
   :maxdepth: 2
   :caption: Reference

   reproducibility
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
