import sys
from aux_functions import (
    check_estimator_installed, set_functions_params, parse_options, handle_options,
    print_table, helper_headers, handle_errors,
    print_warnings, check_ntru, export_to_csv, get_secret_value
)
from param_calls import process_parameters
from fit_formula import find_constants


def main(argv):

    opts = parse_options(argv)

    if any(opt == "--fit" for opt, _ in opts):
        find_constants(opts)
        return

    output_dict, l, secret_dist, error_dist, param, lwe_d, logq, verify, ntru_flag, table, hw, num_only, correction, error_dist_tag = handle_options(
        opts)

    if handle_errors(error_dist.stddev, logq, lwe_d, l, param):
        return

    estimator_installed = check_estimator_installed()
    if not estimator_installed and verify:
        print("Lattice Estimator not installed, running without verification.")
        verify = False

    lambda_usvp, lambda_usvp_s, lambda_bdd, lambda_bdd_s, n_usvp, n_usvp_s, n_bdd, n_bdd_s = set_functions_params(
        get_secret_value(opts))

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
        'model_values': model_values,
        'verify': verify,
        'estimator_installed': estimator_installed,
        'correction': correction,
        'secret_dist': secret_dist,
        'error_dist': error_dist,
        'error_tag': error_dist_tag,
        'hw': hw,
        'output_dict': output_dict,
        'num_only': num_only,
    }

    if ntru_flag:
        check_ntru(params)

    data = process_parameters(params, table)

    export_to_csv(data, "output.csv")

    if param in ['n', 'logq', 'std_e', 'lambda']:
        headers = list(data[0].keys()) if data else []
        data_values = [list(d.values()) for d in data]
        helper_headers(headers)
        print_table(headers, data_values)
    else:
        print("Parameter " + param + " not valid")
        return

    print_warnings(verify, estimator_installed)

    if param != "est":
        return output_dict[param]


if __name__ == "__main__":
    main(sys.argv[1:])
