from numpy import pi, exp, log, log2, sqrt, divide
from scipy.optimize import fsolve, minimize

import numpy as np
import warnings

const = 2 * pi * exp(1)
ln2 = log(2)
e = exp(1)


def _delta(beta):
    """
    Calculate the delta value for a given beta.

    :param beta: Beta value.
    :return: Delta value.
    """
    return (beta / (2 * pi * e) * (pi * beta) ** (1 / beta)) ** (1 / (2 * (beta - 1)))


def numerical_lambda_bdd(n, logq, std_s, std_e):
    """
    Estimate the lambda value for the BDD model using numerical methods.

    :param n: Dimension.
    :param logq: Logarithm of the modulus.
    :param std_s: Standard deviation of the secret.
    :param std_e: Standard deviation of the error.
    :return: Lambda value.
    """
    lnq = logq * ln2
    zeta = max(1, round(std_e / std_s))

    # Initial guess for beta
    beta_initial_guess = n / 4

    def nom(beta): return 2 * n * lnq * log(beta / const)
    def denom(beta): return log(beta / const) + 2 * lnq - 2 * log(std_e) - log(const) - \
        2 * (lnq - log(zeta)) * sqrt(n * log(beta / const) / (2 * lnq * beta))

    def eq6(beta): return beta - nom(beta) / (denom(beta) ** 2)

    beta_solution = fsolve(eq6, beta_initial_guess, full_output=False)

    d_optimal = sqrt(
        2 * n * lnq * beta_solution[0] / log(beta_solution[0] / const))

    # Compute lambda
    l_solution = 0.292 * beta_solution[0] + log2(8 * d_optimal) + 16.4

    return l_solution


def numerical_lambda_usvp(n, logq, std_s, std_e):
    """
    Estimate the lambda value for the USVP model using numerical methods.

    :param n: Dimension.
    :param logq: Logarithm of the modulus.
    :param std_s: Standard deviation of the secret.
    :param std_e: Standard deviation of the error.
    :return: Lambda value.
    """
    lnq = logq * ln2
    zeta = max(1, round(std_e / std_s))

    # Initial guess for beta
    beta_initial_guess = n / 4

    def nom(beta): return 2 * n * (lnq - log(zeta)) * log(beta / const)
    def denom(beta): return lnq + log(sqrt(beta) / (const * std_e))
    def eq11(beta): return beta - nom(beta) / (denom(beta) ** 2)

    beta_solution = fsolve(eq11, beta_initial_guess, full_output=False)

    # Compute d (as substitute in eq12)
    d_optimal = sqrt(2 * n * (lnq - log(zeta)) *
                     beta_solution[0] / log(beta_solution[0] / const))

    # Compute lambda
    l_solution = 0.292 * beta_solution[0] + log2(8 * d_optimal) + 16.4

    return l_solution


def numerical_n_bdd(l, logq, std_s, std_e):
    """
    Estimate the n value for the BDD model using numerical methods.

    :param l: Security parameter.
    :param logq: Logarithm of the modulus.
    :param std_s: Standard deviation of the secret.
    :param std_e: Standard deviation of the error.
    :return: n value.
    """
    eta_initial_guess = (l - 16.4) / 0.292
    def eta_eq(eta): return l - (0.292 * eta + log2(8 * eta) + 16.4)
    eta_solution = fsolve(eta_eq, eta_initial_guess, full_output=False)
    eta = eta_solution[0]

    lnq = logq * ln2
    def d_optimal(n, beta): return sqrt(n * lnq / log(_delta(beta)))
    def eq8(n, beta): return l - (0.292 * beta +
                                  log2(8 * d_optimal(n, beta)) + 16.4)

    ln_std_e = log(std_e)
    zeta = max(0, ln_std_e - log(std_s))
    def d(n, beta): return max(d_optimal(n, beta), n)

    def eq_for_n_and_beta(n, beta): return eta - d(n, beta) + 1 / log(_delta(
        beta)) * (lnq - ln_std_e - 0.5 * log(const) - n / d(n, beta) * (lnq - zeta))

    def system_bdd_eta_n(x):  # x[0] = n, x[1] = beta
        f1 = eq_for_n_and_beta(x[0], x[1])
        f2 = eq8(x[0], x[1])
        return f1, f2

    n_initial_guess = 100
    solutions_n_and_beta = fsolve(
        system_bdd_eta_n, [n_initial_guess, eta], full_output=True)
    n_solution = solutions_n_and_beta[0][0]
    if not solutions_n_and_beta[2] == 1:
        print(solutions_n_and_beta[3])
    return n_solution


def numerical_n_usvp(l, logq, std_s, std_e):
    """
    Estimate the n value for the USVP model using numerical methods.

    :param l: Security parameter.
    :param logq: Logarithm of the modulus.
    :param std_s: Standard deviation of the secret.
    :param std_e: Standard deviation of the error.
    :return: n value.
    """
    lnq = logq * ln2
    zeta = max(1, round(std_e / std_s))

    # Initial guesses for n and beta
    n_initial_guess = 100
    beta_initial_guess = (l - 16.4) / 0.292

    def nom(n, beta): return 2 * n * (lnq - log(zeta)) * log(beta / const)
    def denom(beta): return lnq + log(sqrt(beta) / (const * std_e))
    def eq11(n, beta): return beta - nom(n, beta) / (denom(beta) ** 2)

    def d_optimal(n, beta):
        if (beta < const):
            print("Error: could not find optimal d, maybe lambda is too small")
            exit(0)
        else:
            return sqrt(2 * n * lnq * beta / log(beta / const))

    def eq12(n, beta): return l - (0.292 * beta +
                                   log2(8 * d_optimal(n, beta)) + 16.4)

    def system_usvp_l(x):  # x[0] = n, x[1] = beta
        f1 = eq11(x[0], x[1])
        f2 = eq12(x[0], x[1])
        return f1, f2

    n_solution, beta_solution = fsolve(
        system_usvp_l, [n_initial_guess, beta_initial_guess], full_output=False)
    return n_solution


def numerical_logq_bdd(l, n, std_s, std_e):
    """
    Estimate the logq value for the BDD model using numerical methods.

    :param l: Security parameter.
    :param n: Dimension.
    :param std_s: Standard deviation of the secret.
    :param std_e: Standard deviation of the error.
    :return: logq value.
    """
    
    beta = (l - 16.4) / 0.292

    if (beta < const):
        print("Error: will not start optimization, most likely provided lambda is too small")
        exit(0)

    lnq_initial_guess = 2*n*log(beta/const) / beta #comes from a simplfied the formula below: (remove log(beta / const), log(zeta) and - 2 * log(std_e) - log(const) from denom)
    zeta = std_e / std_s
    
    def nom(lnq): return 2 * n * lnq * log(beta / const)
    def denom(lnq): return log(beta / const) + 2 * lnq - 2 * log(std_e) - log(const) - \
        2 * (lnq - log(zeta)) * sqrt(n * log(beta / const) / (2 * lnq * beta))

    def eq(lnq): return beta - nom(lnq) / (denom(lnq) ** 2)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lnq_solution, info, ier, msg = fsolve(
            eq, lnq_initial_guess, full_output=True)
        if ier != 1:
            print(f"Warning: numerical solver for bdd did not converge: {msg}")
        return divide(lnq_solution[0], ln2)


def numerical_logq_usvp(l, n, std_s, std_e):
    """
    Estimate the logq value for the USVP model using numerical methods.

    :param l: Security parameter.
    :param n: Dimension.
    :param std_s: Standard deviation of the secret.
    :param std_e: Standard deviation of the error.
    :return: logq value.
    """
    
    beta = (l - 16.4) / 0.292
    if (beta < const):
        print("Error: will not start optimization, most likely provided lambda is too small")
        exit(0)

    zeta = std_e / std_s

    lnq_initial_guess = ( 2*n*log(beta/log(beta/const)) - beta/log(beta/const) + sqrt(n**2*log(beta/log(beta/const))**2 - n*log(beta/log(beta/const))*beta/log(beta/const)) )/(2*beta/log(beta/const))
    assert(lnq_initial_guess>1)
    
    def nom(lnq): return 2 * n * (lnq - log(zeta)) * log(beta / const)
    def denom(lnq): return lnq + log(sqrt(beta) / (const * std_e))
    def eq(lnq): return beta - nom(lnq) / (denom(lnq) ** 2)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        lnq_solution, info, ier, msg = fsolve(
            eq, lnq_initial_guess, full_output=True)
        if ier != 1:
            print(
                f"Warning: numerical solver for usvp did not converge: {msg}")
        return divide(lnq_solution[0], ln2)


def numerical_std_e_bdd(l, n, logq, std_s):
    """
    Estimate the std_e value for the BDD model using numerical methods.

    :param l: Security parameter.
    :param n: Dimension.
    :param logq: Logarithm of the modulus.
    :param std_s: Standard deviation of the secret.
    :return: std_e value.
    """
    eta_initial_guess = (l - 16.4) / 0.292
    def eta_eq(eta): return l - (0.292 * eta + log2(8 * eta) + 16.4)
    eta_solution = fsolve(eta_eq, eta_initial_guess, full_output=False)

    lnq = logq * ln2
    def d_optimal(beta): return sqrt(n * lnq / log(_delta(beta)))
    def eq8(beta): return l - (0.292 * beta + log2(8 * d_optimal(beta)) + 16.4)

    beta_solution = fsolve(eq8, eta_solution, full_output=False)
    d = max(d_optimal(beta_solution[0]), n)
    eta = eta_solution[0]
    beta = beta_solution[0]

    std_e_initial_guess = 5.502177429822036
    def zeta(ln_std_e): return max(0, ln_std_e - log(std_s))
    def eq_for_sigmae(ln_std_e): return eta - d + 1 / log(_delta(beta)) * \
        (lnq - ln_std_e - 0.5 * log(const) - n / d * (lnq - zeta(ln_std_e)))

    std_e_solution = fsolve(
        eq_for_sigmae, std_e_initial_guess, full_output=True)

    if not std_e_solution[2] == 1:
        print(std_e_solution[3])

    return exp(std_e_solution[0][0])


def numerical_std_e_usvp(l, n, logq, std_s):
    """
    Estimate the std_e value for the USVP model using numerical methods.

    :param l: Security parameter.
    :param n: Dimension.
    :param logq: Logarithm of the modulus.
    :param std_s: Standard deviation of the secret.
    :return: std_e value.
    """
    lnq = logq * ln2

    # Initial guesses
    beta_initial_guess = (l - 16.4) / 0.292
    std_e_initial_guess = 3.19

    def zeta(std_e): return max(1, np.round(std_e / std_s))

    def nom(std_e, beta): return 2 * n * \
        (lnq - log(zeta(std_e))) * log(beta / const)

    def denom(std_e, beta):
        value_inside_log = sqrt(beta) / (const * std_e)
        return lnq + log(value_inside_log)

    def eq11(x):
        std_e, beta = x
        return beta - nom(std_e, beta) / (denom(std_e, beta) ** 2)

    def d_optimal(std_e, beta):
        if (beta < const):
            print("Error: could not find optimal d, maybe lambda is too small")
            exit(0)
        else:
            return sqrt(2 * n * (lnq - log(zeta(std_e))) * beta / log(beta / const))

    def eq12(x):
        std_e, beta = x
        return l - (0.292 * beta + log2(8 * d_optimal(std_e, beta)) + 16.4)

    def objective(x):
        return eq11(x) ** 2 + eq12(x) ** 2

    # Bounds to ensure std_e is positive
    bounds = [(1e-3, None), (2, None)]

    initial_guess = [std_e_initial_guess, beta_initial_guess]

    std_e_solution = minimize(objective, initial_guess, bounds=bounds)

    if not std_e_solution.success:
        print(f"Optimization did not converge: {std_e_solution.message}")

    return std_e_solution.x[0]
