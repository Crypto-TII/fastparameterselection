import sys
from aux_functions import (
    check_estimator_installed, set_functions_params, parse_options, handle_options,
    print_table, helper_headers, handle_errors,
    update_headers, print_warnings, check_ntru
)
from param_calls import process_parameters, process_subdata


def main(argv):

    opts = parse_options(argv)
    output_dict, l, secret, param, lwe_d, logq, verify, ntru_flag, std_s, std_e, secret_q, table, hw = handle_options(
        opts)

    estimator_installed = check_estimator_installed()

    if handle_errors(std_e, logq, lwe_d, l, param):
        return

    if not estimator_installed and verify:
        print("Lattice Estimator not installed, can't run verification")
        return

    lambda_usvp, lambda_usvp_s, lambda_bdd, lambda_bdd_s, n_usvp, n_usvp_s, n_bdd, n_bdd_s = set_functions_params(
        secret)

    model_values = {
        'lambda_usvp': lambda_usvp,
        'lambda_usvp_s': lambda_usvp_s,
        'lambda_bdd': lambda_bdd,
        'lambda_bdd_s': lambda_bdd_s,
        'n_usvp': n_usvp,
        'n_usvp_s': n_usvp_s,
        'n_bdd': n_bdd,
        'n_bdd_s': n_bdd_s
    }

    params = {
        'param': param,
        'logq': logq,
        'l': l,
        'lwe_d': lwe_d,
        'std_s': std_s,
        'std_e': std_e,
        'model_values': model_values,
        'verify': verify,
        'estimator_installed': estimator_installed,
        'secret': secret,
        'secret_q': secret_q,
        'hw': hw,
        'output_dict': output_dict
    }

    headers, data = process_parameters(params)

    if param in ['n', 'logq', 'std_e', 'lambda']:
        if not table:
            subdata = process_subdata(
                param, data, verify, estimator_installed, secret)
            headers = update_headers(param, verify, estimator_installed)
            helper_headers(headers)
            print_table(headers, subdata)
        else:
            helper_headers(headers)
            print_table(headers, data)

    if ntru_flag:
        check_ntru(output_dict)

    print_warnings(verify, estimator_installed)

    if param != "est":
        return output_dict[param]


if __name__ == "__main__":
    main(sys.argv[1:])
