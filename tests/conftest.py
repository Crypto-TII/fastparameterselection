"""Test configuration.

The tool is a set of scripts that import each other by bare module name, so
``src/`` has to be importable directly rather than as a package.
"""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


@pytest.fixture(scope="session")
def repo_root():
    """Absolute path to the repository root."""
    return ROOT


@pytest.fixture(scope="session")
def estimate_py():
    """Absolute path to the command line entry point."""
    return SRC / "estimate.py"


def _sage_available():
    """src/nd.py imports sage.all, so anything running the CLI needs SageMath."""
    try:
        import sage.all  # noqa: F401
    except Exception:
        return False
    return True


SAGE_AVAILABLE = _sage_available()


def pytest_collection_modifyitems(config, items):
    if SAGE_AVAILABLE:
        return
    skip = pytest.mark.skip(reason="SageMath not available (src/nd.py imports sage.all)")
    for item in items:
        if "needs_sage" in item.keywords:
            item.add_marker(skip)
