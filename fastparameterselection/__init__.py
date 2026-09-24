"""Fast and secure LWE parameter selection.

Command line::

    fastparams --param "lambda" --n 1024 --logq 27 --secret binary --std 3.19

From Python, the closed-form formulas and numerical solvers are importable
directly::

    from fastparameterselection import const
    from fastparameterselection.formulas import model_lambda_usvp

    model_lambda_usvp(1024, 27, 0.5, 3.19, const.LAMBDA_USVP_BIN)

See https://fastparameterselection.readthedocs.io for the full documentation.
"""

__version__ = "1.0.0"

__all__ = ["const", "formulas", "numerical_solver", "rng", "__version__"]
