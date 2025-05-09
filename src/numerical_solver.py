import traceback
from numpy import pi, exp, log, log2, sqrt, divide
from scipy.optimize import fsolve, minimize

import numpy as np
import warnings
import random

import inspect
import matplotlib.pyplot as plt

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
            print(
                f"Error in {numerical_n_usvp.__name__}: could not find optimal d, maybe lambda is too small"
            )
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
        print(
            "Error: will not start optimization, most likely provided lambda is too small")
        exit(0)

    # comes from a simplfied the formula below: (remove log(beta / const), log(zeta) and - 2 * log(std_e) - log(const) from denom)
    lnq_initial_guess = 2*n*log(beta/const) / beta
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
        print(
            "Error: will not start optimization, most likely provided lambda is too small")
        exit(0)

    zeta = std_e / std_s

    lnq_initial_guess = (2*n*log(beta/log(beta/const)) - beta/log(beta/const) + sqrt(n**2*log(beta/log(
        beta/const))**2 - n*log(beta/log(beta/const))*beta/log(beta/const)))/(2*beta/log(beta/const))
    assert (lnq_initial_guess > 1)

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


def numerical_std_e_bdd_minimize(l, n, logq, std_s):
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
    retry_range = [1.0, 10.0]
    max_retries = 5
    eta_initial_guess = (l - 16.4) / 0.292

    def eta_eq(eta):
        return l - (0.292 * eta + log2(8 * eta) + 16.4)

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

    def objective(ln_std_e):
        return eq_for_sigmae(ln_std_e) ** 2

    # Retry logic
    for attempt in range(max_retries):
        print(f"Attempt {attempt + 1} of {max_retries}...")

        # Generate a random starting point within the specified range
        std_e_initial_guess = random.uniform(*retry_range)

        # Minimize the objective function
        # bounds = [(1e-3, None)]  # Ensure std_e is positive

        std_e_solution = minimize(
            objective, [std_e_initial_guess])

        if std_e_solution.success:
            print(
                f"numerical_std_e_bdd_minimize Converged on attempt {attempt + 1} with std_e={exp(std_e_solution.x[0])}")
            return exp(std_e_solution.x[0]), std_e_solution.success

        print(
            f"Attempt {attempt + 1} failed to converge: {std_e_solution.message}")

    print("All attempts failed to converge.")
    return None


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

    max_retries = 100
    eta_initial_guess = (l - 16.4) / 0.292

    def eta_eq(eta):
        return l - (0.292 * eta + log2(8 * eta) + 16.4)

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
        std_e_initial_guess = random.uniform(*retry_range)

        print(
            f"Attempt {attempt + 1} of {max_retries}... initial guess {std_e_initial_guess}")
        # Solve using fsolve
        std_e_solution, info, ier, msg = fsolve(
            eq_for_sigmae, std_e_initial_guess, full_output=True)

        std_e_solution = std_e_solution[0]

        if ier == 1:  # Check if fsolve converged
            print(
                f"numerical_std_e_bdd Converged on attempt {attempt + 1} with std_e={exp(std_e_solution)}")
            return exp(std_e_solution), bool(ier)

        print(f"Attempt {attempt + 1} failed to converge: {msg}")

    print("All attempts failed to converge.")
    return None


def numerical_std_e_usvp_minimize(l, n, logq, std_s):
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

    def nom(std_e, beta):
        # Ensure zeta(std_e) is positive
        zeta_value = zeta(std_e)
        if zeta_value <= 0:
            print(
                f"Warning: zeta(std_e) is non-positive (zeta={zeta_value}). Setting zeta to a small positive value.")
            zeta_value = 1e-10  # Set to a small positive value

        # Ensure beta / const is positive
        beta_ratio = beta / const
        if beta_ratio <= 0:
            print(
                f"Warning: beta / const is non-positive (beta_ratio={beta_ratio}). Setting beta_ratio to a small positive value.")
            beta_ratio = 1e-10  # Set to a small positive value

        return 2 * n * (lnq - log(zeta_value)) * log(beta_ratio)

    def denom(std_e, beta):
        # Ensure beta is positive
        if beta <= 0:
            print(
                f"Warning: beta is non-positive (beta={beta}). Setting beta to a small positive value.")
            beta = 1e-10  # Set beta to a small positive value

        value_inside_log = sqrt(beta) / (const * std_e)

        # Ensure value_inside_log is positive
        if value_inside_log <= 0:
            print(
                f"Warning: value_inside_log is non-positive (value_inside_log={value_inside_log}). Setting it to a small positive value.")
            value_inside_log = 1e-10  # Set to a small positive value

        return lnq + log(value_inside_log)

    def eq11(x):
        std_e, beta = x
        return beta - nom(std_e, beta) / (denom(std_e, beta) ** 2)

    def d_optimal(std_e, beta):
        # This check is not need if we add bounds to minimize
        if (beta < const):
            print(
                f"Error in {numerical_std_e_usvp.__name__}: could not find optimal d, maybe lambda is too small"
            )
            exit(0)
        else:
            return sqrt(2 * n * (lnq - log(zeta(std_e))) * beta / log(beta / const))

    def eq12(x):
        std_e, beta = x
        return l - (0.292 * beta + log2(8 * d_optimal(std_e, beta)) + 16.4)

    def objective(x):
        return eq11(x) ** 2 + eq12(x) ** 2

    # Bounds to ensure std_e is positive and beta is bigger than const
    # bounds = [(1e-3, const), (None, None)]

    initial_guess = [std_e_initial_guess, beta_initial_guess]

    std_e_solution = minimize(objective, initial_guess)

    if not std_e_solution.success:
        sol = std_e_solution.x
        residuals = objective(sol)
        print(
            f"Optimization did not converge in {numerical_std_e_usvp.__name__}: {std_e_solution.message}")
        print("Solution found:", sol)
        print("Residuals at solution:", residuals)

        # Default value if optimization fails
        return std_e_initial_guess, std_e_solution.success

    return exp(std_e_solution.x[0]), std_e_solution.success


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

    max_retries = 100

    # Initial guesses
    beta_initial_guess = (l - 16.4) / 0.292
    std_e_initial_guess = 3.19

    retry_range = [2**-20, 2**0.34]
    retry_range_beta = [const, beta_initial_guess + 100]

    def zeta(std_e):
        return max(1, np.round(std_e / std_s))

    def nom(std_e, beta):
        zeta_value = zeta(std_e)
        if zeta_value <= 0:
            print(
                f"Warning: zeta(std_e) is non-positive (zeta={zeta_value}). Setting zeta to a small positive value.")
            zeta_value = 1e-10  # Set to a small positive value

        beta_ratio = beta / const
        if beta_ratio <= 0:
            print(
                f"Warning: beta / const is non-positive (beta_ratio={beta_ratio}). Setting beta_ratio to a small positive value.")
            beta_ratio = 1e-10  # Set to a small positive value

        return 2 * n * (lnq - log(zeta_value)) * log(beta_ratio)

    def denom(std_e, beta):
        if beta <= 0:
            print(
                f"Warning: beta is non-positive (beta={beta}). Setting beta to a small positive value.")
            beta = 1e-10  # Set beta to a small positive value

        value_inside_log = sqrt(beta) / (const * std_e)
        if value_inside_log <= 0:
            print(
                f"Warning: value_inside_log is non-positive (value_inside_log={value_inside_log}). Setting it to a small positive value.")
            value_inside_log = 1e-10  # Set to a small positive value

        return lnq + log(value_inside_log)

    def eq11(x):
        std_e, beta = x
        nom_value = nom(std_e, beta)
        denom_value = denom(std_e, beta)
        if nom_value is None or denom_value is None:
            return None  # Invalid value
        return beta - nom_value / (denom_value ** 2)

    def d_optimal(std_e, beta):
        zeta_value = zeta(std_e)
        beta_ratio = beta / const
        if zeta_value <= 0:
            print(
                f"Warning: zeta(std_e) is non-positive (zeta={zeta_value}). Setting zeta to a small positive value.")
            zeta_value = 1e-10  # Set to a small positive value

        if beta_ratio <= 0:
            print(
                f"Warning: beta / const is non-positive (beta_ratio={beta_ratio}). Setting beta_ratio to a small positive value.")
            beta_ratio = 1e-10  # Set to a small positive value

        value_inside_sqrt = 2 * n * \
            (lnq - log(zeta_value)) * beta / log(beta_ratio)
        if value_inside_sqrt <= 0:
            print(
                f"Warning: value_inside_sqrt is non-positive (value_inside_sqrt={value_inside_sqrt}). Setting it to a small positive value.")
            value_inside_sqrt = 1e-10  # Set to a small positive value

        return sqrt(value_inside_sqrt)

    def eq12(x):
        std_e, beta = x
        return l - (0.292 * beta + log2(8 * d_optimal(std_e, beta)) + 16.4)

    def system_of_equations(x):
        """
        Combine eq11 and eq12 into a system of equations for fsolve.
        """
        eq11_value = eq11(x)
        eq12_value = eq12(x)
        if eq11_value is None or eq12_value is None:
            return None  # Invalid system
        return [eq11_value, eq12_value]

    # Retry logic
    valid_attempts = 0
    for attempt in range(max_retries):
        try:
            # Generate a random starting point within the specified range
            initial_guess = [std_e_initial_guess, beta_initial_guess]
            if eq11(std_e_solution) < 0:
                initial_guess[1] += step_beta

                # Validate the system before calling fsolve
            if system_of_equations(initial_guess) is None:
                print(f"Skipping invalid initial guess: {initial_guess}")
                continue  # Skip invalid guesses without counting toward retries

            print(
                f"Attempt {valid_attempts + 1} of {max_retries}... initial guess {initial_guess}")

            # Solve using fsolve
            std_e_solution, info, ier, msg = fsolve(
                system_of_equations, initial_guess, full_output=True
            )

            if ier == 1:  # Check if fsolve converged
                std_e_solution = std_e_solution[0]
                print(
                    f"numerical_std_e_usvp Converged on attempt {valid_attempts + 1} with std_e={std_e_solution}")
                if std_e_solution < 1e-3:
                    print(f"Warning: std_e is too small ({std_e_solution})")
                else:
                    return exp(std_e_solution), bool(ier)

            print(f"Attempt {valid_attempts + 1} failed to converge: {msg}")
            print("eq11:", eq11(std_e_solution))
            print("eq12:", eq12(std_e_solution))
            exit(0)
            valid_attempts += 1  # Increment only for valid attempts

        except Exception as e:
            print(f"Error during attempt {valid_attempts + 1}: {e}")
            traceback.print_exc()  # Print the full traceback for debugging

    print("All attempts failed to converge.")
    return None


def plot_system_of_equations(l, n, logq, std_s, min_valid_points=500, solution_threshold=1e-3):
    """
    Plot the system of equations defined in numerical_std_e_usvp and highlight solutions.

    :param l: Security parameter.
    :param n: Dimension.
    :param logq: Logarithm of the modulus.
    :param std_s: Standard deviation of the secret.
    :param min_valid_points: Minimum number of valid points required to plot.
    :param solution_threshold: Threshold for identifying solutions (close to zero).
    """
    lnq = logq * ln2

    def zeta(std_e):
        return max(1, round(std_e / std_s))

    def nom(std_e, beta):
        zeta_value = zeta(std_e)
        beta_ratio = beta / const
        if zeta_value <= 0 or beta_ratio <= 0:
            return None  # Invalid value
        return 2 * n * (lnq - np.log(zeta_value)) * np.log(beta_ratio)

    def denom(std_e, beta):
        value_inside_log = np.sqrt(beta) / (const * std_e)
        if value_inside_log <= 0:
            return None  # Invalid value
        return lnq + np.log(value_inside_log)

    def eq11(std_e, beta):
        nom_value = nom(std_e, beta)
        denom_value = denom(std_e, beta)
        if nom_value is None or denom_value is None:
            return None  # Invalid value
        return beta - nom_value / (denom_value ** 2)

    def d_optimal(std_e, beta):
        zeta_value = zeta(std_e)
        beta_ratio = beta / const
        if zeta_value <= 0 or beta_ratio <= 0:
            return None  # Invalid value
        value_inside_sqrt = 2 * n * \
            (lnq - np.log(zeta_value)) * beta / np.log(beta_ratio)
        if value_inside_sqrt <= 0:
            return None  # Invalid value
        return np.sqrt(value_inside_sqrt)

    def eq12(std_e, beta):
        d_value = d_optimal(std_e, beta)
        if d_value is None:
            return None  # Invalid value
        return l - (0.292 * beta + np.log2(8 * d_value) + 16.4)

    def system_value(std_e, beta):
        """
        Combine eq11 and eq12 into a single value.
        The closer this value is to zero, the closer the point is to solving the system.
        """
        eq11_value = eq11(std_e, beta)
        eq12_value = eq12(std_e, beta)
        if eq11_value is None or eq12_value is None:
            return None  # Invalid value
        return np.sqrt(eq11_value**2 + eq12_value**2)  # Combine the equations

    # Generate points until we have enough valid data
    valid_points = []
    while len(valid_points) < min_valid_points:
        std_e_values = np.random.uniform(
            1e-3, 10, 1000)  # Avoid zero for std_e
        beta_values = np.random.uniform(const, 100, 1000)
        for std_e, beta in zip(std_e_values, beta_values):
            system_val = system_value(std_e, beta)
            if system_val is not None:
                valid_points.append((std_e, beta, system_val))
                if len(valid_points) >= min_valid_points:
                    break

    # Extract valid data for plotting
    std_e_plot = [p[0] for p in valid_points]
    beta_plot = [p[1] for p in valid_points]
    system_plot = [p[2] for p in valid_points]

    # Identify solution points
    solution_points = [(p[0], p[1])
                       for p in valid_points if p[2] < solution_threshold]
    solution_std_e = [p[0] for p in solution_points]
    solution_beta = [p[1] for p in solution_points]

    # Create a scatter plot for the system
    plt.figure(figsize=(10, 8))
    scatter = plt.scatter(std_e_plot, beta_plot,
                          c=system_plot, cmap='viridis', alpha=0.7)
    plt.colorbar(scatter, label="System Value (sqrt(eq11^2 + eq12^2))")

    # Highlight solution points
    plt.scatter(solution_std_e, solution_beta, color='red',
                label='Solutions', edgecolor='black', s=50)

    # Add labels and title
    plt.xlabel("std_e")
    plt.ylabel("beta")
    plt.title("System of Equations (eq11 and eq12)")
    plt.legend()
    plt.grid(True)
    plt.show()
