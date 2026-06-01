from numpy import pi, exp, log, log2, sqrt
from scipy.optimize import fsolve

import numpy as np
import random

const = 2 * pi * exp(1)
ln2 = log(2)
e = exp(1)


def _delta(beta):
    """
    Calculate the delta value for a given beta.

    :param beta: Beta value.
    :return: Delta value.
    """
    beta = max(beta, 100)  # Ensure beta is at positive
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

    def denom(beta):
        # print(n * log(beta / const) / (2 * lnq * beta))
        return log(beta / const) + 2 * lnq - 2 * log(std_e) - log(const) - 2 * (lnq - log(zeta)) * sqrt(n * log(beta / const) / (2 * lnq * beta))

    def eq6(beta): return beta - nom(beta) / (denom(beta) ** 2)

    beta_solution = fsolve(eq6, beta_initial_guess, full_output=False)

    d_optimal = sqrt(
        2 * n * lnq * beta_solution[0] / log(beta_solution[0] / const))

    # Compute lambda
    l_solution = 0.292 * beta_solution[0] + log2(8 * d_optimal) + 16.4

    return l_solution


def numerical_lambda_bdd_rev1(n, logq, std_s, std_e):

    # print("in numerical_lambda_bdd_rev1!")
    lnq = logq * ln2
    zeta = max(1, round(std_e / std_s))

    # Computing approximate beta for the initial guess
    if n/(lnq) >= 1:
        beta_approx = 2*n/lnq*(log(n/(lnq)))
    else:
        beta_approx = n / 6

    # Initial guess for beta
    beta_initial_guess = beta_approx

    def d_opt(beta):
        if beta < 2:
            return 10*5
        return sqrt(n*(lnq - log(zeta))*2*beta/log(beta/const))

    def eta(beta):
        if beta < 2:
            return 10*5
        return beta + (log2(d_opt(beta)))/0.292


    def main_eq(beta):
        """
        eta = d - 1/(ln(delta_beta))( ln(q/(sigma_1 * sqrt(const))) - n*ln(q/zeta)/d )
        """
        return eta(beta) - d_opt(beta) + 2*beta/(log(beta/const))*(lnq - log(std_e) - 0.5*log(const) - n/d_opt(beta)*(lnq-log(zeta)))

    beta_solution = fsolve(main_eq, beta_initial_guess, full_output=False)


    l_solution = 0.292 * beta_solution[0] + \
        log2(8 * d_opt(beta_solution[0])) + 16.4

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
    def eta_eq(eta): return l - (0.292 * eta + 16.4)
    eta_solution = fsolve(eta_eq, eta_initial_guess, full_output=False)
    eta = eta_solution[0]

    lnq = logq * ln2

    def d_optimal(n, beta): return sqrt(
        n * lnq / log(_delta(beta)))  # TODO missing log(zeta) here? Same issue in the Estimator
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
            print(
                f"Error in {numerical_n_usvp.__name__}: could not find optimal d, maybe lambda is too small"
            )
            return 1
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

    eta_initial_guess = (l - 16.4) / 0.292
    def eta_eq(eta): return l - (0.292 * eta + 16.4)
    eta_solution = fsolve(eta_eq, eta_initial_guess, full_output=False)
    eta = eta_solution[0]

    def d_optimal(lnq, beta):
        lnq = max(lnq, 1)  # Ensure lnq is at positive
        return sqrt(n * lnq / log(_delta(beta)))
    def eq8(lnq, beta): return l - (0.292 * beta +
                                    log2(8 * d_optimal(lnq, beta)) + 16.4)

    ln_std_e = log(std_e)
    zeta = max(0, ln_std_e - log(std_s))
    def d(lnq, beta): return max(d_optimal(lnq, beta), n)

    def eq_for_lnq_and_beta(lnq, beta): return eta - d(lnq, beta) + 1 / log(_delta(
        beta)) * (lnq - ln_std_e - 0.5 * log(const) - n / d(lnq, beta) * (lnq - zeta))

    def system_bdd_eta_n(x):  # x[0] = n, x[1] = beta
        f1 = eq_for_lnq_and_beta(x[0], x[1])
        f2 = eq8(x[0], x[1])
        return f1, f2

    lnq_initial_guess = 121
    solutions_lnq_and_beta = fsolve(
        system_bdd_eta_n, [lnq_initial_guess, eta], full_output=True)
    lnq_solution = solutions_lnq_and_beta[0][0]
    if not solutions_lnq_and_beta[2] == 1:
        print(solutions_lnq_and_beta[3])
    return lnq_solution/ln2


def numerical_std_e_usvp(l, n, logq, std_s):
    """
    Estimate the n value for the USVP model using numerical methods.

    :param l: Security parameter.
    :param logq: Logarithm of the modulus.
    :param std_s: Standard deviation of the secret.
    :param std_e: Standard deviation of the error.
    :return: n value.
    """
    lnq = logq * ln2
    # zeta = max(1, round(std_e / std_s))

    # Initial guesses for n and beta
    std_e_initial_guess = 3.19
    beta_initial_guess = (l - 16.4) / 0.292

    def zeta(std_e):
        # print("std_e: ", std_e, "std_s: ", std_s)
        return max(1, np.round(std_e / std_s))

    def nom(std_e, beta): return 2 * n * \
        (lnq - log(zeta(std_e))) * log(beta / const)

    def denom(std_e, beta): return lnq + log(sqrt(beta) / (const * std_e))

    def eq11(std_e, beta): return beta - \
        nom(std_e, beta) / (denom(std_e, beta) ** 2)

    def d_optimal(beta):
        if (beta < const):
            print(
                f"Error in {numerical_n_usvp.__name__}: could not find optimal d, maybe lambda is too small"
            )
            return 1
        else:
            # return sqrt(2 * n * (lnq - log(zeta(std_e))) * beta / log(beta / const))
            return sqrt(2 * n * lnq * beta / log(beta / const))

    def eq12(beta): return l - (0.292 * beta +
                                log2(8 * d_optimal(beta)) + 16.4)

    def system_usvp_l(x):  # x[0] = n, x[1] = beta
        f1 = eq11(x[0], x[1])
        f2 = eq12(x[1])
        return f1, f2

    # std_e_solution, beta_solution = fsolve(
    #     system_usvp_l, [std_e_initial_guess, beta_initial_guess], full_output=False)
    # return std_e_solution, 1

    # Initial guesses for fsolve
    initial_guess = [std_e_initial_guess, beta_initial_guess]

    # Solve using fsolve
    guess_range = 50
    guess_range_beta = beta_initial_guess + 200
    max_attempts = 1000

    ier = 0
    msg = ""
    solution = 0

    for i in range(max_attempts):
        # Random guess

        if i != 0:
            initial_guess[0] = np.random.uniform(-guess_range,
                                                 guess_range, size=1)
            initial_guess[1] = np.random.uniform(-guess_range_beta,
                                                 guess_range_beta, size=1)

        print(f"Checking {i+1} with initial guess {initial_guess}")
        try:
            solution, info, ier, msg = fsolve(
                system_usvp_l, initial_guess, full_output=True, xtol=1e-3, maxfev=1000)
        except Exception as e:
            continue

        if ier == 1:
            break

    if ier != 1:  # Check if fsolve converged
        print(
            f"Warning: fsolve did not converge in {numerical_std_e_usvp.__name__}: {msg}")
        print("Solution found (if any):", solution)
        return std_e_initial_guess, False  # Default value if fsolve fails

    std_e_solution, beta_solution = solution
    print(
        f"numerical_std_e_usvp Converged: std_e={std_e_solution}, beta={beta_solution}")

    return std_e_solution, True


def numerical_logq_usvp(l, n, std_s, std_e):
    """
    Estimate the n value for the USVP model using numerical methods.

    :param l: Security parameter.
    :param logq: Logarithm of the modulus.
    :param std_s: Standard deviation of the secret.
    :param std_e: Standard deviation of the error.
    :return: n value.
    """
    zeta = max(1, round(float(std_e / std_s)))

    # Initial guesses for n and beta
    lnq_initial_guess = 20
    beta_initial_guess = (l - 16.4) / 0.292

    def nom(lnq, beta): return 2 * n * (lnq - log(zeta)) * \
        log(beta / const)  # TODO: check if we want to keep log(zeta) here

    def denom(lnq, beta): return lnq + log(sqrt(beta) / (const * std_e))
    def eq11(lnq, beta): return beta - nom(lnq, beta) / (denom(lnq, beta) ** 2)

    def d_optimal(lnq, beta):
        if (beta < const):
            print(
                f"Error in {numerical_n_usvp.__name__}: could not find optimal d, maybe lambda is too small"
            )
            return 1
        else:
            if lnq < 0:  # Ensure lnq is at positive
                print(f"lnq is negative: {lnq}")
                return 1
            return sqrt(2 * n * lnq * beta / log(beta / const))

    def eq12(lnq, beta): return l - (0.292 * beta +
                                     log2(8 * d_optimal(lnq, beta)) + 16.4)

    def system_usvp_l(x):  # x[0] = n, x[1] = beta
        f1 = eq11(x[0], x[1])
        f2 = eq12(x[0], x[1])
        return f1, f2

    lnq_solution, beta_solution = fsolve(
        system_usvp_l, [lnq_initial_guess, beta_initial_guess], full_output=False)
    return lnq_solution/ln2


def numerical_std_e_bdd(l, n, logq, std_s):
    """
    Estimate the std_e value for the BDD model using numerical methods, with retry logic for convergence.

    :param l: Security parameter.
    :param n: Dimension.
    :param logq: Logarithm of the modulus.
    :param std_s: Standard deviation of the secret.
    :param retry_range: Tuple specifying the range for random starting points.
    :param max_retries: Maximum number of retries.
    :return: std_e value.
    """
    retry_range = [-20, 60]
    std_e_initial_guess = 1

    max_retries = 100
    eta_initial_guess = (l - 16.4) / 0.292

    def eta_eq(eta):
        return l - (0.292 * eta + 16.4)

    eta_solution = fsolve(eta_eq, eta_initial_guess, full_output=False)

    lnq = logq * ln2

    def d_optimal(beta):
        return sqrt(n * lnq / log(_delta(beta)))

    def eq8(beta):
        return l - (0.292 * beta + log2(8 * d_optimal(beta)) + 16.4)

    beta_solution = fsolve(eq8, eta_solution, full_output=False)
    d = max(d_optimal(beta_solution[0]), n)
    eta = eta_solution[0]
    beta = beta_solution[0]

    def zeta(ln_std_e):
        return max(0, ln_std_e - log(std_s))

    def eq_for_sigmae(ln_std_e):
        return eta - d + 1 / log(_delta(beta)) * \
            (lnq - ln_std_e - 0.5 * log(const) - n / d * (lnq - zeta(ln_std_e)))

    # Retry logic
    for attempt in range(max_retries):

        # Generate a random starting point within the specified range
        if attempt != 0:
            std_e_initial_guess = random.uniform(*retry_range)

        print(
            f"Attempt {attempt + 1} of {max_retries}... initial guess {std_e_initial_guess}")
        # Solve using fsolve
        std_e_solution, info, ier, msg = fsolve(
            eq_for_sigmae, std_e_initial_guess, full_output=True)

        std_e_solution = std_e_solution[0]

        if ier == 1:  # Check if fsolve converged
            print(
                f"numerical_std_e_bdd Converged on attempt {attempt + 1} with std_e={std_e_solution, exp(std_e_solution)}")
            print("solution: ", eq_for_sigmae(std_e_solution))
            return exp(std_e_solution), bool(ier)

        print(f"Attempt {attempt + 1} failed to converge: {msg}")

    print("All attempts failed to converge.")
    return None
