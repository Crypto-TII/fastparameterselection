import math
import csv
import sys
import getopt
from formulas import check_overstreched
from const import (
    LAMBDA_USVP_BIN, LAMBDA_USVP_TER, LAMBDA_USVP_S_BIN,
    LAMBDA_USVP_S_TER, LAMBDA_BDD_BIN, LAMBDA_BDD_TER, LAMBDA_BDD_S_BIN,
    LAMBDA_BDD_S_TER, N_USVP_BIN, N_USVP_TER, N_USVP_S_BIN, N_USVP_S_TER,
    N_BDD_BIN, N_BDD_TER, N_BDD_S_BIN, N_BDD_S_TER
)

sys.path.append('./latticeestimator')


def check_estimator_installed():
    try:
        global LWE, ND
        from estimator import LWE, ND
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
        opts, args = getopt.getopt(argv, "h,v", [
                                   "secret=", "error=", "param=", "n=", "lambda=", "logq=", "file=", "hw=", "ntru", "table", "num-only"])
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
        if l <= 0:
            errors.append("Error: lambda must be bigger than 0")

    return errors


def set_secret(secret, output_dict):
    """
    Set the secret distribution and its standard deviation.

    :param secret: Secret distribution (binary or ternary).
    :param output_dict: Dictionary to store the output values.
    :return: Tuple of standard deviation of the secret and secret modulus.
    """
    std_s = 0
    secret_q = 0
    if secret == 'binary':
        std_s = UniformModStd(2)
        secret_q = 2
        output_dict['std_s'] = 0.5
    elif secret == 'ternary':
        std_s = UniformModStd(3)
        secret_q = 3
        output_dict['std_s'] = math.sqrt(2./3)
    else:
        if secret != 'sparse':
            print("Secret distribution not supported")
            sys.exit()
    return std_s, secret_q


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
    logq = 0
    std_e = 3.19  # Default value for the standard deviation of the error
    secret = "binary"  # Default value for the secret distribution
    # Default value for the standard deviation of the secret
    std_s, secret_q = set_secret(secret, output_dict)
    l = 0
    hw = 0
    table = False
    num_only = False
    for opt, arg in opts:
        if opt == '--help' or opt == '-h':
            helper()
        elif opt == '--hw':
            hw = int(arg)
        elif opt == '--secret':
            secret = arg
            std_s, secret_q = set_secret(secret, output_dict)
        elif opt == '--error':
            std_e = float(arg)
            output_dict['std_e'] = std_e
        elif opt == '--param':
            param = arg
        elif opt == '--n':
            try:
                lwe_d = int(arg)
            except:
                print("Error: Invalid LWE dimension format")
                sys.exit()
            output_dict['n'] = lwe_d
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
        else:
            helper()
    return output_dict, l, secret, param, lwe_d, logq, verify, ntru_flag, std_s, std_e, secret_q, table, hw, num_only

# Exctracted from the Lattice Estimator


def UniformModStd(q):
    """
    Calculate the standard deviation of a uniform distribution modulo q.

    :param q: Modulus.
    :return: Standard deviation.
    """
    a = -(q // 2)
    b = -a - 1 if q % 2 == 0 else -a

    if b < a:
        raise ValueError(
            f"upper limit must be larger than lower limit but got: {b} < {a}")
    m = b - a + 1
    mean = (a + b) / float(2)
    stddev = math.sqrt((m**2 - 1) / float(12))

    return stddev


def load_all_from_csv(file_path):
    """
    Load all entries from a CSV file.

    :param file_path: Path to the CSV file.
    :return: List of entries.
    """
    with open(file_path, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        entries = [row for row in reader]
    return entries


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
    print(
        "  --secret <secret>       Secret distribution (binary, ternary, sparse) [optional, default: binary]")
    print(
        "  --error <error>         Standard deviation of the error (e.g., 3.19) [optional, default: 3.19]")
    print("  --hw <hw>               Hamming weight (for sparse secrets) (e.g., 64)")
    print("  --verify                Verify results against the Lattice Estimator")
    print("  --table                 Output results from all the formulas")
    print("  --ntru                  Check NTRU parameters")
    print("  -h, --help              Show this help message and exit")
    print("\nExamples:")
    print('  python3 src/estimate.py --param "lambda" --n "1024" --logq "20;24-28;30;33;37;42" --secret "binary" --error "3.19"')
    print('  python3 src/estimate.py --param "n" --lambda "80" --logq "20-30" --secret "binary" --error "3.19"')
    print('  python3 src/estimate.py --param "logq" --lambda "80" --n "1024" --secret "binary" --error "3.19"')
    print('  python3 src/estimate.py --param "std_e" --lambda "80" --n "1024" --logq "20" --secret "binary"')
    print('  python3 src/estimate.py --param "lambda" --n "1024" --logq "40" --hw "64" --secret "sparse"')
    print('  python3 src/estimate.py --param "logq" --lambda "80" --n "1024" --hw "64" --secret "sparse"')
    print('\nWith verification and table output:')
    print('  python3 src/estimate.py --param "lambda" --n "1024" --logq "20" --secret "binary" --error "3.19" -v --table')
    print('  python3 src/estimate.py --param "n" --lambda "80" --logq "20" --secret "binary" --error "3.19" -v --table')
    print('  python3 src/estimate.py --param "logq" --lambda "80" --n "1024" --secret "binary" --error "3.19" -v --table')
    print('  python3 src/estimate.py --param "std_e" --lambda "80" --n "1024" --logq "20" --secret "binary" -v --table')
    print('  python3 src/estimate.py --param "lambda" --n "1024" --logq "40" --hw "64" --secret "sparse" -v --table')
    print('  python3 src/estimate.py --param "logq" --lambda "80" --n "1024" --hw "64" --secret "sparse" -v --table')
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
        "bdd": "Output of the formula which estimates the cost of the (unique) BDD attack",
        "bdd_s": "Output of the simplified formula (removing dependency on beta) which estimates the cost of the (unique) BDD attack",
        "logq usvp": "Output of the numerical approximation of log q for the (unique) SVP attack",
        "logq bdd": "Output of the numerical approximation of log q for the (unique) BDD attack",
        "usvp num": "Output of the numerical approximation of the (unique) SVP attack",
        "bdd num": "Output of the numerical approximation of the (unique) BDD attack",
        "log2(std_e) usvp": "Output of the numerical approximation of the (log2) standard deviation of the error for the (unique) SVP attack",
        "log2(std_e) bdd": "Output of the numerical approximation of the (log2) standard deviation of the error for the (unique) BDD attack",
        "bdd 3.19": "The result of running the Lattice Estimator with standard deviation of the error 3.19 and primal_bdd",
        "usvp 3.19": "The result of running the Lattice Estimator with standard deviation of the error 3.19 and primal_usvp",
        "diff": "The difference between the output of the previous column and the output of the Lattice Estimator",
        "est usvp": "Output of the Lattice Estimator for the (unique) SVP attack",
        "est bdd": "Output of the Lattice Estimator for the (unique) BDD attack",
        "est usvp_s": "Output of the Lattice Estimator using the result from the simplified formula for the (unique) SVP attack",
        "est bdd_s": "Output of the Lattice Estimator using the result from the simplified formula for the (unique) BDD attack",
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
