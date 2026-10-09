"""Test configuration.

Tests import the package from the checkout, so the repository root goes on
sys.path rather than the package directory itself.
"""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
PKG = ROOT / "fastparameterselection"

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


@pytest.fixture(scope="session")
def repo_root():
    """Absolute path to the repository root."""
    return ROOT


@pytest.fixture(scope="session")
def estimate_py():
    """Argument list that invokes the command line entry point.

    The package uses relative imports, so the module form is required; running
    the file by path would fail on the first ``from .const import ...``.
    """
    return ["-m", "fastparameterselection.estimate"]


def _sage_available():
    """nd.py imports sage.all, so anything running the CLI needs SageMath."""
    try:
        import sage.all  # noqa: F401
    except Exception:
        return False
    return True


SAGE_AVAILABLE = _sage_available()


def pytest_collection_modifyitems(config, items):
    if SAGE_AVAILABLE:
        return
    skip = pytest.mark.skip(reason="SageMath not available (nd.py imports sage.all)")
    for item in items:
        if "needs_sage" in item.keywords:
            item.add_marker(skip)
