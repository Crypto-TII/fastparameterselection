"""Deterministic per-call random number generators."""


import hashlib
import random

DEFAULT_SEED = 0

_base_seed = DEFAULT_SEED


def set_base_seed(seed):
    """Set the process-wide base seed that every per-call generator derives from."""
    global _base_seed
    _base_seed = int(seed)


def get_base_seed():
    """Return the current base seed."""
    return _base_seed


def call_rng(name, *args):
    """Return a :class:`random.Random` for one solver call.

    :param name: identifier of the call site, so two solvers given the same
        parameters do not share a stream.
    :param args: the call's own parameters; anything with a stable ``repr``.
    :return: a generator seeded from the base seed, ``name`` and ``args``.
    """
    key = repr((_base_seed, name) + args).encode("utf-8")
    digest = hashlib.blake2b(key, digest_size=16).digest()
    return random.Random(int.from_bytes(digest, "big"))
