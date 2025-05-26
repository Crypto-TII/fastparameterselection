import math
import csv
import sys
import getopt
from formulas import check_overstreched
from const import (
    LAMBDA_USVP_BIN, LAMBDA_USVP_TER, LAMBDA_USVP_S_BIN,
    LAMBDA_USVP_S_TER, LAMBDA_BDD_BIN, LAMBDA_BDD_TER, LAMBDA_BDD_S_BIN,
    LAMBDA_BDD_S_TER, N_USVP_BIN, N_USVP_TER, N_USVP_S_BIN, N_USVP_S_TER,
    N_BDD_BIN, N_BDD_TER, N_BDD_S_BIN, N_BDD_S_TER, POINTS, VALUES
)

from nd import NoiseDistribution as ND

sys.path.append('./latticeestimator')


def check_estimator_installed():
    try:
        global LWE, RC
        from estimator import LWE, RC
        return True
    except ImportError:
        print("Warning: Failed to import lattice_estimator, some options will not work")
        return False


def check_ntru(output_dict):
    beta_ = check_overstreched(output_dict)
    if beta_ > 0 and output_dict['lambda'] > 0 and (output_dict['lambda'] - 0.292 * beta_) > 20:
        print("Warning: the ntru parameters are in the overstretched regime")


def print_warnings(verify, estimator_installed):
    print("\n")
    if verify and not estimator_installed:
        print("Warning: Verification not possible, Lattice Estimator not installed")
    print("\n")


def handle_errors(std_e, logq, lwe_d, l, param):
    errors = check_parameters(std_e, logq, lwe_d, l, param)
    if errors:
        for error in errors:
            print(error)
        return True
    return False


def set_functions_params(secret):
    """
    Set the lambda functions based on the secret distribution.

    :param secret: Secret distribution (binary or ternary).
    :return: Tuple of lambda functions.
    """
    if secret == "binary":
        return LAMBDA_USVP_BIN, LAMBDA_USVP_S_BIN, LAMBDA_BDD_BIN, LAMBDA_BDD_S_BIN, N_USVP_BIN, N_USVP_S_BIN, N_BDD_BIN, N_BDD_S_BIN
    else:
        return LAMBDA_USVP_TER, LAMBDA_USVP_S_TER, LAMBDA_BDD_TER, LAMBDA_BDD_S_TER, N_USVP_TER, N_USVP_S_TER, N_BDD_TER, N_BDD_S_TER


def parse_options(argv):
    """
    Parse command-line options.

    :param argv: List of command-line arguments.
    :return: List of options and arguments.
    """
    try:
        opts, args = getopt.getopt(argv, "a,b,h,v,c", [
                                   "secret=", "error=", "param=", "n=", "lambda=", "logq=", "file=", "hw=",  "std=", "eta=", "ntru", "table", "num-only"])
    except Exception as e:
        print(e)
        helper()
    if len(opts) == 0:
        helper()

    return opts


def check_parameters(std_e, logq, lwe_d, l, param):
    """
    Check the validity of the parameters.

    :param std_e: Standard deviation of the error.
    :param logq: List of log q values.
    :param lwe_d: LWE dimension.
    :param l: Security parameter.
    :param param: Parameter to be checked.
    :return: List of error messages.
    """
    errors = []

    if param != 'std_e':
        if std_e == 0:
            errors.append(
                "Error: std_e = 0, the LWE problem can be solved in polynomial time, impossible to reach the desired security level")
        elif std_e < 0:
            errors.append("Error: std_e must be bigger than 0")

    if param != 'logq':
        for q in logq:
            if q <= 0:
                errors.append("Error: logq must be bigger than 0")

    if param != 'n':
        if lwe_d <= 0:
            errors.append("Error: LWE dimension must be bigger than 0")

    if param != 'lambda':
        if l <= 0 and param != 'est':
            errors.append("Error: lambda must be bigger than 0")

    return errors


def set_distribution(secret, params):
    """
    Set the secret/error distribution and its standard deviation.
    """

    if secret == 'binary':
        secret_dist = ND.UniformMod(2)
    elif secret == 'ternary':
        secret_dist = ND.UniformMod(3)
    elif secret == 'sparse':
        try:
            secret_dist = ND.SparseTernary(
                n=params['n'], p=params['hw']/2, m=params['hw']/2)
        except:
            print("Error: hamming weight --hw is required for sparse secret")
            sys.exit()
    elif secret == 'uniformmod':
        secret_dist = ND.UniformMod(params['q'])
    elif secret == 'uniform':
        try:
            secret_dist = ND.Uniform(params['a'], params['b'])
        except:
            print("Error: interval bounds --a and --b are required for uniform secret")
            sys.exit()
    elif secret == 'gaussian':
        secret_dist = ND.DiscreteGaussian(params['std'])
        # try:
        #     secret_dist = ND.DiscreteGaussian(params['std'])
        # except:
        #     print("Error: standard deviation --std is required for gaussian secret")
        #     sys.exit()
    elif secret == 'binomial':
        try:
            secret_dist = ND.CenteredBinomial(params['eta'])
        except:
            print("Error: parameter --eta is required for binomial secret")
            sys.exit()
    else:
        print("Secret distribution not supported")
        print("Available options: binary, ternary, sparse, uniformmod, uniform, gaussian, binomial")
        print("dist: ", secret)
        sys.exit()
    return secret_dist


def get_secret_value(opts):

    for opt, arg in opts:
        if opt == '--secret':
            return arg
    return None


def handle_options(opts):
    """
    Handle the command-line options.

    :param opts: List of options and arguments.
    :return: Tuple of output dictionary and various parameters.
    """
    output_dict = {}
    verify = 0
    ntru_flag = False
    lwe_d = 0
    hw = 0
    logq = 0
    secret_dist_tag = "binary"  # Default value for the secret distribution
    error_dist_tag = "gaussian"  # Default value for the error distribution
    # Default value for the standard deviation of the secret
    params = {
        'n': lwe_d,       # LWE dimension
        'hw': hw,         # Hamming weight
        'std': 3.19,     # Standard deviation for Gaussian
        'a': 0,           # Lower bound for uniform distribution
        'b': 1,           # Upper bound for uniform distribution
        'eta': 1,         # Parameter for binomial distribution
        'q': 2            # Modulus for uniformmod distribution
    }
    secret_dist = set_distribution(secret_dist_tag, params)
    error_dist = set_distribution(error_dist_tag, params)
    l = 0
    table = False
    num_only = False
    correction = False

    for opt, arg in opts:
        if opt == '--help' or opt == '-h':
            helper()
        elif opt == '--hw':
            hw = int(arg)
            params['hw'] = hw
        elif opt == '--param':
            param = arg
        elif opt == '--n':
            try:
                lwe_d = int(arg)
            except:
                print("Error: Invalid LWE dimension format")
                sys.exit()
            output_dict['n'] = lwe_d
            params['n'] = lwe_d
        elif opt == '--lambda':
            l = int(arg)
            output_dict['lambda'] = l
        elif opt == '--logq':
            logq = parse_logq(arg)
            output_dict['logq'] = logq
        elif opt == '-v':
            verify = 1
        elif opt == '--ntru':
            ntru_flag = True
        elif opt == '--table':
            table = True
        elif opt == '--num-only':
            num_only = True
        elif opt == '-c':
            correction = True
        elif opt == '--std':
            params['std'] = float(arg)
        elif opt == '-a':
            params['a'] = float(arg)
        elif opt == '-b':
            params['b'] = float(arg)
        elif opt == '--eta':
            params['eta'] = float(arg)
        elif opt == '--secret':
            secret_dist_tag = str(arg)
        elif opt == '--error':
            error_dist_tag = str(arg)
        else:
            helper()

    secret_dist = set_distribution(secret_dist_tag, params)
    error_dist = set_distribution(error_dist_tag, params)

    if secret_dist_tag != 'binary' and secret_dist_tag != 'ternary':
        num_only = True
    if error_dist_tag != 'gaussian':
        num_only = True

    return output_dict, l, secret_dist, error_dist, param, lwe_d, logq, verify, ntru_flag, table, hw, num_only, correction, error_dist_tag


def export_to_csv(data, output_file):
    """
    Export data to a CSV file.

    :param data: List of dictionaries containing table rows.
    :param output_file: Path to the output CSV file.
    """
    if not data:
        print("No data to export.")
        return

    # Get the headers from the keys of the first dictionary
    headers = data[0].keys()

    with open(output_file, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=headers)
        writer.writeheader()
        writer.writerows(data)

    print(f"Data exported to {output_file}")


def closest_power_of_2(n):
    """
    Find the closest power of 2 to a given number.

    :param n: Input number.
    :return: Closest power of 2.
    """
    if n <= 0:
        raise ValueError("Input must be a positive number.")

    # Calculate the power of 2 just below and above the number
    lower_pow = 2 ** math.floor(math.log2(n))
    upper_pow = 2 ** math.ceil(math.log2(n))

    # Determine which is closer
    if abs(n - lower_pow) < abs(n - upper_pow):
        return lower_pow
    else:
        return upper_pow


def print_table(headers, rows):
    """
    Print a table with headers and rows.

    :param headers: List of headers.
    :param rows: List of rows.
    """
    # Calculate the maximum width for each column
    col_widths = [max(len(str(cell)) for cell in col)
                  for col in zip(headers, *rows)]

    # Create a format string for each row
    row_format = " | ".join(["{:<" + str(width) + "}" for width in col_widths])

    # Print the header
    print(row_format.format(*headers))

    # Print the separator
    print("-+-".join(['-' * width for width in col_widths]))

    # Print the rows
    for row in rows:
        formatted_row = [f"{value:.2f}" if isinstance(
            value, float) else str(value) for value in row]
        print(row_format.format(*formatted_row))


def parse_logq(logq_str):
    """
    Parse the logq string into a list of logq values.

    :param logq_str: Logq string.
    :return: List of logq values.
    """
    logq = []
    parts = logq_str.split(';')
    for part in parts:
        if '-' in part:
            try:
                start, end = map(int, part.split('-'))
            except:
                print("Error: Invalid logq format")
                sys.exit()
            logq.extend(range(start, end + 1))
        else:
            logq.append(int(part))
    return logq


def helper():
    """
    Print the helper message and exit.
    """
    print("Usage: python3 src/estimate.py [OPTIONS]")
    print("\nOptions:")
    print("  --param <param>         Parameter to estimate (lambda, n, logq, std_e, est)")
    print("  --n <n>                 LWE dimension (e.g., 1024)")
    print("  --lambda <lambda>       Security parameter (e.g., 80)")
    print("  --logq <logq>           Log q values (e.g., 20;24-28;30;33;37;42)")
    print("  --secret <secret>       Secret distribution (binary, ternary, sparse, uniformmod, uniform, gaussian, binomial)")
    print("  --error <error>         Error distribution (binary, ternary, sparse, uniformmod, uniform, gaussian, binomial)")
    print("  --hw <hw>               Hamming weight (for sparse secrets) (e.g., 64)")
    print("  --std <std>             Standard deviation for Gaussian distribution (e.g., 3.19)")
    print("  --eta <eta>             Parameter for binomial distribution (e.g., 1)")
    print("  -a <a>                  Lower bound for uniform distribution (e.g., 0)")
    print("  -b <b>                  Upper bound for uniform distribution (e.g., 1)")
    print("  --verify                Verify results against the Lattice Estimator")
    print("  --table                 Output results from all the formulas")
    print("  --ntru                  Check NTRU parameters")
    print("  --num-only              Output only numerical results")
    print("  -c                      Apply correction logic")
    print("  -h, --help              Show this help message and exit")
    print("\nExamples:")
    print('  # Example 1: Estimate lambda with binary secret')
    print('  python3 src/estimate.py --param "lambda" --n "1024" --logq "20;24-28;30;33;37;42" --secret "binary" --error "gaussian" --std "3.19"')
    print('\n  # Example 2: Estimate n with sparse secret')
    print('  python3 src/estimate.py --param "n" --lambda "80" --logq "20-30" --secret "sparse" --hw "64"')
    print('\n  # Example 3: Estimate logq with ternary secret')
    print('  python3 src/estimate.py --param "logq" --lambda "80" --n "1024" --secret "ternary" --error "gaussian" --std "3.19"')
    print('\n  # Example 4: Estimate std_e')
    print('  python3 src/estimate.py --param "std_e" --lambda "80" --n "1024" --logq "20" --secret "binary" --error "gaussian"')
    print('\n  # Example 5: Check NTRU parameters')
    print('  python3 src/estimate.py --param "lambda" --n "1024" --logq "40" --hw "64" --secret "sparse" --ntru')
    print('\n  # Example 6: Output results in table format with verification using the Lattice Estimator')
    print('  python3 src/estimate.py --param "lambda" --n "1024" --logq "20" --secret "binary" --error "gaussian" --std "3.19" -v --table')
    print('\n  # Example 7: Apply correction logic for logq')
    print('  python3 src/estimate.py --param "logq" --lambda "80" --n "1024" --secret "binary" --error "gaussian" --std "3.19" -c')
    sys.exit()


def create_explanation_dict(headers):
    """
    Create a dictionary of explanations for the headers.

    :param headers: List of headers.
    :return: Dictionary of explanations.
    """
    explanations = {
        "secret dist.": "The distribution of the secret (can be binary, ternary or sparse)",
        "lwe dim.": "The Learning With Errors (LWE) dimension",
        "lambda": "The security level",
        "log q": "The size of the modulus q in bits",
        "lwe est": "The output of running the Lattice Estimator using the output of our formulas and the rest of the LWE parameters",
        "usvp": "Output of the formula which estimates the cost of the (unique) SVP attack",
        "usvp_s": "Output of the simplified formula (removing dependency on beta) which estimates the cost of the (unique) SVP attack",
        "bdd": "Output of the formula which estimates the cost of the BDD attack",
        "bdd_s": "Output of the simplified formula (removing dependency on beta) which estimates the cost of the BDD attack",
        "logq usvp": "Output of the numerical approximation of log q for the (unique) SVP attack",
        "logq bdd": "Output of the numerical approximation of log q for the BDD attack",
        "usvp num": "Output of the numerical approximation of the (unique) SVP attack",
        "bdd num": "Output of the numerical approximation of the BDD attack",
        "log2(std_e) usvp": "Output of the numerical approximation of the (log2) standard deviation of the error for the (unique) SVP attack",
        "log2(std_e) bdd": "Output of the numerical approximation of the (log2) standard deviation of the error for the BDD attack",
        "bdd 3.19": "The result of running the Lattice Estimator with standard deviation of the error 3.19 and primal_bdd",
        "usvp 3.19": "The result of running the Lattice Estimator with standard deviation of the error 3.19 and primal_usvp",
        "diff": "The difference between the output of the previous column and the output of the Lattice Estimator",
        "est usvp": "Output of the Lattice Estimator for the (unique) SVP attack",
        "est bdd": "Output of the Lattice Estimator for the BDD attack",
        "est usvp_s": "Output of the Lattice Estimator using the result from the simplified formula for the (unique) SVP attack",
        "est bdd_s": "Output of the Lattice Estimator using the result from the simplified formula for the BDD attack",
        "output": "Recommended value to be used considering all the outputs of the formulas and numerical methods",
        "pow": "Closest power of 2 to the LWE dimension recommended in Output",
        "hw": "Hamming weight of the secret",
        "hybrid": "Output of the numerical approximation for lambda of the hybrid attack",
        "logq hybrid": "Output of the numerical approximation for logq of the hybrid attack",
        "est hybrid": "Output of the Lattice Estimator for the hybrid attack",
        "est": "Output of the Lattice Estimator"
    }

    # Create a dictionary using the headers and explanations
    explanation_dict = {}
    for header in headers:
        # Add the explanation if it exists in the explanations dictionary, otherwise use a default message
        explanation_dict[header] = explanations.get(
            header, "No explanation available for this header.")

    return explanation_dict


def helper_headers(header):
    """
    Print the headers and their explanations.

    :param header: List of headers.
    """
    explanation_dict = create_explanation_dict(header)

    max_length = max(len(header) for header in explanation_dict.keys())
    max_length_exp = max(len(explanation)
                         for explanation in explanation_dict.values())

    # Print each header and its explanation with proper formatting
    for header, explanation in explanation_dict.items():
        print(f"{header:<{max_length}}: {explanation}")

    print("." * max_length_exp)
    print('\n')


def get_parameters(lwe_d, lnq, secret_dist, error_dist, est_usvp_numerical, est_bdd_numerical, error_dist_tag, param):

    lwe_parameters_usvp = None
    lwe_parameters_bdd = None

    if param == 'logq':
        lwe_parameters_usvp = LWE.Parameters(
            lwe_d, 2 ** est_usvp_numerical, secret_dist, error_dist)
        lwe_parameters_bdd = LWE.Parameters(
            lwe_d, 2 ** est_bdd_numerical, secret_dist, error_dist)
    if param == 'std_e':
        # TODO: numerics for std only works for gaussian distribution
        error_dist_usvp = set_distribution(
            error_dist_tag, {'std': 2**est_usvp_numerical})
        lwe_parameters_usvp = LWE.Parameters(
            lwe_d, 2**lnq, secret_dist, error_dist_usvp)
        error_dist_bdd = set_distribution(
            error_dist_tag, {'std': 2**est_bdd_numerical})

        print("Error dist bdd", error_dist_bdd.stddev)

        lwe_parameters_bdd = LWE.Parameters(
            lwe_d, 2**lnq, secret_dist, error_dist_bdd)

    return lwe_parameters_usvp, lwe_parameters_bdd


def correction_logic(l, lwe_d, lnq, lwe_usvp, lwe_bdd, secret_dist, error_dist, est_usvp_numerical, est_bdd_numerical, param, num_calls, error_dist_tag=None):

    print("lwe usvp", lwe_usvp, "lwe bdd", lwe_bdd, "l", l)

    if lwe_usvp >= l:
        while (lwe_usvp >= l):
            corrected_lwe_usvp = lwe_usvp
            corrected_logq_usvp = est_usvp_numerical
            print("Applying correction lwe usvp > l", "logq ",
                  corrected_logq_usvp, "std e ", error_dist.stddev, " est", corrected_lwe_usvp)

            if param == 'logq':
                est_usvp_numerical += 1
            elif param == 'std_e':
                est_usvp_numerical -= 0.5

            num_calls += 1
            lwe_parameters_usvp, _ = get_parameters(
                lwe_d, lnq, secret_dist, error_dist, est_usvp_numerical, est_bdd_numerical, error_dist_tag, param)

            lwe_usvp = math.floor(math.log2(LWE.primal_usvp(
                lwe_parameters_usvp, red_cost_model=RC.BDGL16)["rop"]))
    else:
        while (lwe_usvp < l):
            print("lwe usvp", lwe_usvp, "lwe bdd", lwe_bdd, "l", l)
            if param == 'logq':
                est_usvp_numerical -= 1
            elif param == 'std_e':
                est_usvp_numerical += 0.5
            num_calls += 1
            lwe_parameters_usvp, _ = get_parameters(
                lwe_d, lnq, secret_dist, error_dist, est_usvp_numerical, est_bdd_numerical, error_dist_tag, param)
            lwe_usvp = math.floor(math.log2(LWE.primal_usvp(
                lwe_parameters_usvp, red_cost_model=RC.BDGL16)["rop"]))
            corrected_lwe_usvp = lwe_usvp
            corrected_logq_usvp = est_usvp_numerical
            print("Applying correction lwe usvp < l", "logq ",
                  corrected_logq_usvp, " est", corrected_lwe_usvp)

    if lwe_bdd >= l:
        while (lwe_bdd >= l):
            corrected_lwe_bdd = lwe_bdd
            corrected_logq_bdd = est_bdd_numerical
            print("Applying correction lwe bdd > l", "logq ",
                  corrected_logq_bdd, " est", corrected_lwe_bdd)

            if param == 'logq':
                est_bdd_numerical += 1
            elif param == 'std_e':
                est_bdd_numerical -= 0.5

            num_calls += 1
            _, lwe_parameters_bdd = get_parameters(
                lwe_d, lnq, secret_dist, error_dist, est_usvp_numerical, est_bdd_numerical, error_dist_tag, param)
            try:
                lwe_bdd = math.floor(math.log2(LWE.primal_bdd(
                    lwe_parameters_bdd, red_cost_model=RC.BDGL16)["rop"]))
            except Exception:
                continue
    else:
        while (lwe_bdd < l):

            if param == 'logq':
                est_bdd_numerical -= 1
            elif param == 'std_e':
                est_bdd_numerical += 0.5

            num_calls += 1
            _, lwe_parameters_bdd = get_parameters(
                lwe_d, lnq, secret_dist, error_dist, est_usvp_numerical, est_bdd_numerical, error_dist_tag, param)
            lwe_bdd = math.floor(math.log2(LWE.primal_bdd(
                lwe_parameters_bdd, red_cost_model=RC.BDGL16)["rop"]))
            corrected_lwe_bdd = lwe_bdd
            corrected_logq_bdd = est_bdd_numerical
            print("Applying correction lwe bdd < l", "logq ",
                  corrected_logq_bdd, " est", corrected_lwe_bdd)

    print("Number of calls to the estimator: ", num_calls)

    # TODO: verify that returning max makes sense in all cases
    return max(corrected_logq_bdd, corrected_logq_usvp), corrected_logq_bdd, corrected_logq_usvp, corrected_lwe_bdd, corrected_lwe_usvp
