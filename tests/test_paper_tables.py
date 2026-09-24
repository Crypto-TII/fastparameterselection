"""Values printed in the paper, pinned so they cannot drift unnoticed.

Every expectation here is a cell of a table in "A Tool for Fast and Secure LWE
Parameter Selection: the FHE case", regenerated on this branch. None of them
needs the Lattice Estimator, so the whole file runs in a couple of seconds.

If one of these fails, either the paper and the code have parted company or a
formula changed on purpose. Both are worth stopping for.
"""

import math

import pytest

from fastparameterselection import const
from fastparameterselection.formulas import (
    model_lambda_bdd,
    model_lambda_bdd_s,
    model_lambda_usvp,
    model_lambda_usvp_s,
    model_n_bdd_rev1,
    model_n_bdd_s,
    model_n_usvp,
    model_n_usvp_s,
)
from fastparameterselection.numerical_solver import (
    numerical_lambda_bdd_rev1,
    numerical_lambda_usvp,
    numerical_logq_bdd,
    numerical_logq_usvp,
    numerical_n_usvp,
)

# Standard deviations of the secret distributions, as src/nd.py defines them.
STD_BINARY = 0.5
STD_TERNARY = math.sqrt(2.0 / 3.0)
STD_E = 3.19


def _lam(value):
    """Round a security level the way param_calls.estimate_usvp_bdd does."""
    return int(round(value))


def _dim(value):
    """Round an LWE dimension the way param_calls.process_n_param does."""
    return int(math.ceil(value))


# --------------------------------------------------------------------------
# Tables 2 and 4: security level, binary secret, Eq. (25), (27), (28), (31)
# --------------------------------------------------------------------------

# (n, log q, Eq25, Eq27, Eq28, Eq31)
LAMBDA_BINARY = [
    (1024, 20, 178, 177, 171, 173),
    (1024, 27, 128, 128, 125, 125),
    (1024, 42, 81, 81, 79, 79),
    (2048, 37, 193, 193, 190, 190),
    (2048, 54, 128, 128, 127, 126),
    (2048, 84, 82, 81, 80, 80),
]


@pytest.mark.parametrize("n,logq,eq25,eq27,eq28,eq31", LAMBDA_BINARY)
def test_lambda_binary(n, logq, eq25, eq27, eq28, eq31):
    assert _lam(model_lambda_usvp(n, logq, STD_BINARY, STD_E, const.LAMBDA_USVP_BIN)) == eq25
    assert _lam(model_lambda_usvp_s(n, logq, const.LAMBDA_USVP_S_BIN)) == eq27
    assert _lam(model_lambda_bdd(n, logq, STD_BINARY, STD_E, const.LAMBDA_BDD_BIN)[0].real) == eq28
    assert _lam(model_lambda_bdd_s(n, logq, const.LAMBDA_BDD_S_BIN)) == eq31


# --------------------------------------------------------------------------
# Tables 3 and 5: security level, ternary secret
# --------------------------------------------------------------------------

LAMBDA_TERNARY = [
    (1024, 16, 250, 231, 228, 232),
    (1024, 48, 71, 69, 70, 70),
    (32768, 880, 128, 128, 128, 128),
    (32768, 1450, 77, 77, 77, 76),
    (65536, 1776, 128, 128, 128, 127),
    (131072, 1918, 259, 255, 256, 254),
]


@pytest.mark.parametrize("n,logq,eq25,eq27,eq28,eq31", LAMBDA_TERNARY)
def test_lambda_ternary(n, logq, eq25, eq27, eq28, eq31):
    assert _lam(model_lambda_usvp(n, logq, STD_TERNARY, STD_E, const.LAMBDA_USVP_TER)) == eq25
    assert _lam(model_lambda_usvp_s(n, logq, const.LAMBDA_USVP_S_TER)) == eq27
    assert _lam(model_lambda_bdd(n, logq, STD_TERNARY, STD_E, const.LAMBDA_BDD_TER)[0].real) == eq28
    assert _lam(model_lambda_bdd_s(n, logq, const.LAMBDA_BDD_S_TER)) == eq31


# --------------------------------------------------------------------------
# Tables 6 and 8: LWE dimension, binary secret, Eq. (32), (33), (34), (35)
# --------------------------------------------------------------------------

# (lambda, log q, Eq32, Eq33, Eq34, Eq35)
N_BINARY = [
    (80, 42, 1036, 1037, 1047, 1050),
    (80, 84, 2075, 2052, 2090, 2088),
    (128, 27, 1046, 1042, 1062, 1063),
    (128, 54, 2069, 2072, 2088, 2090),
    (140, 49, 2042, 2043, 2059, 2061),
]


@pytest.mark.parametrize("lam,logq,eq32,eq33,eq34,eq35", N_BINARY)
def test_n_binary(lam, logq, eq32, eq33, eq34, eq35):
    assert _dim(model_n_usvp(lam, logq, STD_BINARY, STD_E, const.N_USVP_BIN)) == eq32
    assert _dim(model_n_usvp_s(lam, logq, const.N_USVP_S_BIN)) == eq33
    assert _dim(model_n_bdd_rev1(lam, logq, STD_BINARY, STD_E, const.N_BDD_BIN)) == eq34
    assert _dim(model_n_bdd_s(lam, logq, const.N_BDD_S_BIN)) == eq35


# --------------------------------------------------------------------------
# Tables 7 and 9: LWE dimension, ternary secret
# --------------------------------------------------------------------------

N_TERNARY = [
    (80, 43, 1054, 1061, 1060, 1065),
    (128, 27, 1050, 1045, 1067, 1064),
    (80, 1325, 32413, 30343, 32333, 32542),
    (140, 810, 33207, 32710, 33137, 33874),
]


@pytest.mark.parametrize("lam,logq,eq32,eq33,eq34,eq35", N_TERNARY)
def test_n_ternary(lam, logq, eq32, eq33, eq34, eq35):
    assert _dim(model_n_usvp(lam, logq, STD_TERNARY, STD_E, const.N_USVP_TER)) == eq32
    assert _dim(model_n_usvp_s(lam, logq, const.N_USVP_S_TER)) == eq33
    assert _dim(model_n_bdd_rev1(lam, logq, STD_TERNARY, STD_E, const.N_BDD_TER)) == eq34
    assert _dim(model_n_bdd_s(lam, logq, const.N_BDD_S_TER)) == eq35


# --------------------------------------------------------------------------
# Tables 10 and 14: security level from the numerical solvers
# --------------------------------------------------------------------------

# (n, log q, secret stddev, uSVP num, BDD num)
LAMBDA_NUM = [
    (1024, 13, STD_BINARY, 268, 262),
    (1024, 54, STD_BINARY, 60, 59),
    (2048, 37, STD_BINARY, 192, 190),
    (1024, 14, STD_TERNARY, 263, 258),
    (32768, 475, STD_TERNARY, 255, 255),
]


@pytest.mark.parametrize("n,logq,std_s,usvp,bdd", LAMBDA_NUM)
def test_lambda_numerical(n, logq, std_s, usvp, bdd):
    assert math.floor(numerical_lambda_usvp(n, logq, std_s, STD_E)) == usvp
    assert math.floor(numerical_lambda_bdd_rev1(n, logq, std_s, STD_E)) == bdd


# --------------------------------------------------------------------------
# Table 11: LWE dimension from the uSVP numerical solver
# --------------------------------------------------------------------------

N_NUM_TERNARY = [
    (80, 43, 1047),
    (128, 27, 1008),
    (110, 1000, 32681),
    (140, 810, 32999),
]


@pytest.mark.parametrize("lam,logq,expected", N_NUM_TERNARY)
def test_n_numerical_ternary(lam, logq, expected):
    assert _dim(numerical_n_usvp(lam, logq, STD_TERNARY, STD_E)) == expected


# --------------------------------------------------------------------------
# Tables 12 and 16: maximum log q
# --------------------------------------------------------------------------

# (lambda, n, secret stddev, uSVP log q, BDD log q)
LOGQ = [
    (100, 1024, STD_BINARY, 33, 33),
    (256, 2048, STD_BINARY, 28, 28),
    (128, 1024, STD_TERNARY, 27, 27),
    (128, 32768, STD_TERNARY, 872, 880),
    (256, 32768, STD_TERNARY, 474, 475),
]


@pytest.mark.parametrize("lam,n,std_s,usvp,bdd", LOGQ)
def test_logq(lam, n, std_s, usvp, bdd):
    assert math.floor(numerical_logq_usvp(lam, n, std_s, STD_E)) == usvp
    assert math.floor(numerical_logq_bdd(lam, n, std_s, STD_E)) == bdd


def test_logq_is_reported_conservatively():
    """The binding attack is the one tolerating the smallest modulus.

    Reporting the maximum instead returned a modulus the other attack already
    breaks: at lambda=110, n=512, ternary it gave 16, which the Lattice
    Estimator scores at 109 bits.
    """
    usvp = math.floor(numerical_logq_usvp(110, 512, STD_TERNARY, STD_E))
    bdd = math.floor(numerical_logq_bdd(110, 512, STD_TERNARY, STD_E))
    assert min(usvp, bdd) == 15
