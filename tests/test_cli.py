"""End to end checks of the command line, including the examples in Section 6.

These run ``src/estimate.py`` as a subprocess, so they cover option parsing,
rounding and formatting as well as the formulas.
"""

import re
import subprocess
import sys

import pytest

pytestmark = pytest.mark.needs_sage


def run(estimate_py, *args, cwd=None):
    """Invoke the tool and return its stdout, asserting it exited cleanly."""
    out = subprocess.run(
        [sys.executable, str(estimate_py), *args],
        capture_output=True, text=True, timeout=1800, cwd=cwd,
    )
    assert out.returncode == 0, f"exited {out.returncode}\n{out.stdout}\n{out.stderr}"
    return out.stdout


def outputs(stdout):
    """Last column of each result row, as strings."""
    rows = [l for l in stdout.splitlines() if re.match(r"^\S+\s+\|", l)]
    return [r.rsplit("|", 1)[-1].strip() for r in rows]


# --------------------------------------------------------------------------
# Section 6, "How to use our results in practice"
# --------------------------------------------------------------------------

def test_section6_security_level(estimate_py):
    out = run(estimate_py, "--param", "lambda", "--n", "1024",
              "--logq", "20;35;40", "--secret", "binary", "--std", "3.19")
    assert outputs(out) == ["173", "95", "83"]


def test_section6_lwe_dimension(estimate_py):
    out = run(estimate_py, "--param", "n", "--lambda", "128",
              "--logq", "27;37;45;54", "--secret", "binary",
              "--error", "gaussian", "--std", "3.19")
    # Last column is the closest power of two; the dimension is the one before.
    rows = [l for l in out.splitlines() if l.startswith("Binary")]
    dims = [l.split("|")[-2].strip() for l in rows]
    assert dims == ["1063", "1442", "1746", "2090"]
    assert outputs(out) == ["1024", "1024", "2048", "2048"]


def test_section6_modulus(estimate_py):
    """Section 6 prints 881 here; the paper needs updating to 872.

    881 came from taking the maximum of the two attacks' moduli. The binding
    constraint is the minimum, and the uSVP solver gives 872 for these
    parameters.
    """
    out = run(estimate_py, "--param", "logq", "--lambda", "128", "--n", "32768",
              "--secret", "ternary", "--error", "gaussian", "--std", "3.19")
    assert outputs(out) == ["872"]


def test_section6_error_stddev(estimate_py):
    out = run(estimate_py, "--param", "std_e", "--lambda", "192", "--n", "2048",
              "--logq", "64", "--secret", "binary")
    assert outputs(out) == ["28.60"]


# --------------------------------------------------------------------------
# Behaviour of the command line itself
# --------------------------------------------------------------------------

def test_runs_from_any_directory(estimate_py, tmp_path):
    """The vendored Lattice Estimator must resolve relative to src/, not cwd."""
    out = run(estimate_py, "--param", "lambda", "--n", "1024", "--logq", "20",
              "--secret", "binary", "--std", "3.19", cwd=tmp_path)
    assert outputs(out) == ["173"]


def test_logq_range_syntax(estimate_py):
    out = run(estimate_py, "--param", "lambda", "--n", "1024",
              "--logq", "20;24-26", "--secret", "binary", "--std", "3.19")
    # The reported level is the lowest across the attacks, so these are the
    # `output` column of Tables 2 and 4 for log q 20, 24, 25 and 26.
    assert outputs(out) == ["173", "142", "136", "130"]


@pytest.mark.parametrize("args", [
    ("--param", "lambda", "--n", "1024", "--logq", "20", "--secret", "nosuchdist"),
    ("--param", "lambda", "--n", "1024", "--logq", "not-a-number", "--secret", "binary"),
    ("--param", "lambda", "--n", "1024", "--logq", "20", "--nosuchflag", "1"),
])
def test_bad_input_exits_nonzero(estimate_py, args):
    """Failures must be visible to a shell script, not reported as success."""
    out = subprocess.run([sys.executable, str(estimate_py), *args],
                         capture_output=True, text=True, timeout=600)
    assert out.returncode != 0


def test_binomial_secret_honours_eta(estimate_py):
    """--s-eta and --eta were unreachable: getopt did not accept them."""
    low = run(estimate_py, "--param", "lambda", "--n", "512", "--logq", "12",
              "--secret", "binomial", "--s-eta", "2",
              "--error", "binomial", "--eta", "2")
    high = run(estimate_py, "--param", "lambda", "--n", "512", "--logq", "12",
               "--secret", "binomial", "--s-eta", "3",
               "--error", "binomial", "--eta", "3")
    assert outputs(low) != outputs(high)
