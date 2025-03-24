import warnings
from estimator import LWE, RC, ND
import math
from numpy import log2

from formulas import (
    model_lambda_usvp, model_lambda_usvp_s, model_lambda_bdd, model_lambda_bdd_s,
    model_n_usvp, model_n_usvp_s, model_n_bdd, model_n_bdd_s,
)
from numerical_solver import numerical_n_usvp, numerical_n_bdd, numerical_logq_usvp, numerical_logq_bdd, numerical_std_e_usvp, numerical_std_e_bdd
from numerical_hybrid import numerical_lambda_hybrid_v2, numerical_logq_hybrid
from aux_functions import closest_power_of_2


from const import (
    HEADERS_N_VERIFY, HEADERS_N_NO_VERIFY, HEADERS_LOGQ_VERIFY, HEADERS_LOGQ_NO_VERIFY, HEADERS_LOGQ_HYBRID_VERIFY, HEADERS_LOGQ_HYBRID_NO_VERIFY,
    HEADERS_STD_E_VERIFY, HEADERS_STD_E_NO_VERIFY, HEADERS_LAMBDA_VERIFY, HEADERS_LAMBDA_NO_VERIFY, HEADERS_LAMBDA_ALT_VERIFY, HEADERS_LAMBDA_ALT_NO_VERIFY,
    HEADERS_LAMBDA_HYBRID_VERIFY, HEADERS_LAMBDA_HYBRID_NO_VERIFY, SECRET_DIST, LAMBDA, LOG_Q, USVP, LWE_USVP, USVP_S, LWE_USVP_S, USVP_NUM, LWE_NUM, BDD, LWE_BDD, BDD_S, LWE_BDD_S, BDD_NUM, OUTPUT, POW, LWE_DIM, LOGQ_USVP, LOGQ_BDD, HW, HYBRID, LWE_HYBRID, STD_E_USVP, STD_E_BDD, NUM_USVP, NUM_BDD, EST
)

import sys
sys.path.append('./latticeestimator')

warnings.filterwarnings('error')


def run_verification(lq, secret, est_usvp, est_bdd, est_usvp_pow, est_bdd_pow):
    """
    Run verification using the Lattice Estimator.

    :param lq: Logarithm of the modulus.
    :param secret: Secret distribution.
    :param est_usvp: Estimated usvp value.
    :param est_bdd: Estimated bdd value.
    :param est_usvp_pow: Estimated usvp power of 2 value.
    :param est_bdd_pow: Estimated bdd power of 2 value.
    :return: Tuple of verification results.
    """
    lwe_parameters_usvp = []
    lwe_parameters_bdd = []
    if secret == 'binary':
        lwe_parameters_usvp = LWE.Parameters(
            est_usvp, 2 ** lq, ND.UniformMod(2), ND.DiscreteGaussian(3.19))
        lwe_parameters_bdd = LWE.Parameters(
            est_bdd, 2 ** lq, ND.UniformMod(2), ND.DiscreteGaussian(3.19))
        lwe_parameters_usvp_pow = LWE.Parameters(
            est_usvp_pow, 2 ** lq, ND.UniformMod(2), ND.DiscreteGaussian(3.19))
        lwe_parameters_bdd_pow = LWE.Parameters(
            est_bdd_pow, 2 ** lq, ND.UniformMod(2), ND.DiscreteGaussian(3.19))
    else:
        lwe_parameters_usvp = LWE.Parameters(
            est_usvp, 2 ** lq, ND.UniformMod(3), ND.DiscreteGaussian(3.19))
        lwe_parameters_bdd = LWE.Parameters(
            est_bdd, 2 ** lq, ND.UniformMod(3), ND.DiscreteGaussian(3.19))
        lwe_parameters_usvp_pow = LWE.Parameters(
            est_usvp_pow, 2 ** lq, ND.UniformMod(3), ND.DiscreteGaussian(3.19))
        lwe_parameters_bdd_pow = LWE.Parameters(
            est_bdd_pow, 2 ** lq, ND.UniformMod(3), ND.DiscreteGaussian(3.19))

    lwe_usvp = math.floor(math.log2(LWE.primal_usvp(
        lwe_parameters_usvp, red_cost_model=RC.BDGL16)["rop"]))
    lwe_bdd = math.floor(math.log2(LWE.primal_bdd(
        lwe_parameters_bdd, red_cost_model=RC.BDGL16)["rop"]))
    lwe_usvp_pow = math.floor(math.log2(LWE.primal_usvp(
        lwe_parameters_usvp_pow, red_cost_model=RC.BDGL16)["rop"]))
    lwe_bdd_pow = math.floor(math.log2(LWE.primal_bdd(
        lwe_parameters_bdd_pow, red_cost_model=RC.BDGL16)["rop"]))

    return lwe_usvp, lwe_bdd, lwe_usvp_pow, lwe_bdd_pow


def process_parameters(params, table):
    param = params['param']
    logq = params['logq']
    l = params['l']
    lwe_d = params['lwe_d']
    std_s = params['std_s']
    std_e = params['std_e']
    model_values = params['model_values']
    verify = params['verify']
    estimator_installed = params['estimator_installed']
    secret = params['secret']
    secret_q = params['secret_q']
    hw = params['hw']
    output_dict = params['output_dict']

    if param == 'n':
        headers, data = process_n(logq, l, std_s, std_e, model_values['n_usvp'], model_values['n_usvp_s'],
                                  model_values['n_bdd'], model_values['n_bdd_s'], verify, estimator_installed, secret, secret_q, output_dict)
    elif param == 'logq':
        headers, data = process_logq(
            l, lwe_d, std_s, std_e, verify, estimator_installed, secret, secret_q, hw, output_dict)
    elif param == 'std_e':
        headers, data = process_std_e(
            logq, l, lwe_d, std_s, verify, estimator_installed, secret, secret_q, output_dict)
    elif param == 'lambda':
        data = process_lambda(logq, lwe_d, std_s, std_e, model_values['lambda_usvp'], model_values['lambda_usvp_s'], model_values[
            'lambda_bdd'], model_values['lambda_bdd_s'], verify, estimator_installed, secret, secret_q, hw, table, output_dict)
    elif param == "est":
        headers, data = process_est(logq, lwe_d, std_e, secret_q)
    else:
        headers, data = helper(), []
    return data


def process_n(logq, l, std_s, std_e, n_usvp, n_usvp_s, n_bdd, n_bdd_s, verify, estimator_installed, secret, secret_q, output_dict):
    headers = HEADERS_N_VERIFY if (
        verify and estimator_installed) else HEADERS_N_NO_VERIFY
    data = process_n_param(logq, l, std_s, std_e, n_usvp, n_usvp_s, n_bdd,
                           n_bdd_s, verify, estimator_installed, secret, secret_q, output_dict)
    return headers, data


def process_logq(l, lwe_d, std_s, std_e, verify, estimator_installed, secret, secret_q, hw, output_dict):
    if secret != 'sparse':
        headers = HEADERS_LOGQ_VERIFY if (
            verify and estimator_installed) else HEADERS_LOGQ_NO_VERIFY
        data = process_logq_param(
            l, lwe_d, std_s, std_e, verify, estimator_installed, secret, secret_q, output_dict)
    else:
        headers = HEADERS_LOGQ_HYBRID_VERIFY if verify and estimator_installed else HEADERS_LOGQ_HYBRID_NO_VERIFY
        data = process_logq_param_hybrid(
            l, lwe_d, std_s, std_e, verify, estimator_installed, secret, secret_q, hw, output_dict)
    return headers, data


def process_std_e(logq, l, lwe_d, std_s, verify, estimator_installed, secret, secret_q, output_dict):
    headers = HEADERS_STD_E_VERIFY if (
        verify and estimator_installed) else HEADERS_STD_E_NO_VERIFY
    data = process_std_e_param(
        logq, l, lwe_d, std_s, verify, estimator_installed, secret, secret_q, output_dict)
    return headers, data


def process_lambda(logq, lwe_d, std_s, std_e, lambda_usvp, lambda_usvp_s, lambda_bdd, lambda_bdd_s, verify, estimator_installed, secret, secret_q, h, table, output_dict):
    if secret != 'sparse':
        data = process_lambda_param(logq, lwe_d, std_s, std_e, lambda_usvp, lambda_usvp_s,
                                    lambda_bdd, lambda_bdd_s, verify, estimator_installed, secret, secret_q, table, output_dict)
    else:
        data = process_lambda_param_hybrid(
            logq, lwe_d, std_s, h, secret, verify, estimator_installed, output_dict)

    return data


def process_est(logq, lwe_d, std_e, secret_q):
    print("Running Lattice Estimator")
    for lq in logq:
        parameters = LWE.Parameters(
            lwe_d, 2 ** lq, ND.UniformMod(secret_q), ND.DiscreteGaussian(std_e))
        # LWE.estimate(parameters, red_cost_model=RC.BDGL16)
        print(LWE.arora_gb(parameters))
        # print(LWE.primal_usvp(parameters))
    return [], []


def process_n_param(logq, l, std_s, std_e, n_usvp, n_usvp_s, n_bdd, n_bdd_s, verify, estimator_installed, secret, secret_q, output_dict):
    """
    Process the parameter 'n' and estimate its value using various models and numerical solvers.

    :param file_path: Path to the input file.
    :param logq: List of log q values.
    :param l: Security parameter.
    :param std_s: Standard deviation of the secret.
    :param std_e: Standard deviation of the error.
    :param n_usvp: Parameter for the usvp model.
    :param n_usvp_s: Parameter for the usvp_s model.
    :param n_bdd: Parameter for the bdd model.
    :param n_bdd_s: Parameter for the bdd_s model.
    :param verify: Boolean flag to indicate if verification is needed.
    :param estimator_installed: Boolean flag to indicate if the estimator is installed.
    :param secret: Secret distribution.
    :param secret_q: Secret modulus.
    :param output_dict: Dictionary to store the output values.

    :return: List of data points with estimated values for 'n'.
    """
    data = []
    if len(logq) > 1:
        output_dict['n'] = []
    for lq in logq:
        # Run the formulas for usvp and bdd together with their simplified versions and the numerical solver. Since we are interested in n we round up the result.
        est_usvp = int(math.ceil(model_n_usvp(l, lq, std_s, std_e, n_usvp)))
        est_usvp_s = int(math.ceil(model_n_usvp_s(l, lq, n_usvp_s)))
        est_bdd = int(math.ceil(model_n_bdd(l, lq, std_s, std_e, n_bdd)))
        est_bdd_s = int(math.ceil(model_n_bdd_s(l, lq, std_s, std_e, n_bdd_s)))
        est_usvp_numerical = int(
            math.ceil(numerical_n_usvp(l, lq, std_s, std_e)))
        est_bdd_numerical = int(
            math.ceil(numerical_n_bdd(l, lq, std_s, std_e)))

        # Store the minimum value of n provided from all the formulas
        optimal_value = max(est_usvp, est_usvp_s, est_bdd,
                            est_bdd_s, est_usvp_numerical, est_bdd_numerical)
        if len(logq) > 1:
            output_dict['n'].append(optimal_value)
        else:
            output_dict['n'] = optimal_value

        if verify and estimator_installed:
            lwe_parameters_usvp = LWE.Parameters(
                est_usvp, 2 ** lq, ND.UniformMod(secret_q), ND.DiscreteGaussian(std_e))
            lwe_parameters_bdd = LWE.Parameters(
                est_bdd, 2 ** lq, ND.UniformMod(secret_q), ND.DiscreteGaussian(std_e))
            lwe_parameters_usvp_s = LWE.Parameters(
                est_usvp_s, 2 ** lq, ND.UniformMod(secret_q), ND.DiscreteGaussian(std_e))
            lwe_parameters_bdd_s = LWE.Parameters(
                est_bdd_s, 2 ** lq, ND.UniformMod(secret_q), ND.DiscreteGaussian(std_e))
            lwe_parameters_usvp_num = LWE.Parameters(
                est_usvp_numerical, 2 ** lq, ND.UniformMod(secret_q), ND.DiscreteGaussian(std_e))
            lwe_parameters_bdd_num = LWE.Parameters(
                est_bdd_numerical, 2 ** lq, ND.UniformMod(secret_q), ND.DiscreteGaussian(std_e))

            lwe_usvp = math.floor(math.log2(LWE.primal_usvp(
                lwe_parameters_usvp, red_cost_model=RC.BDGL16)["rop"]))
            lwe_bdd = math.floor(math.log2(LWE.primal_bdd(
                lwe_parameters_bdd, red_cost_model=RC.BDGL16)["rop"]))
            lwe_usvp_s = math.floor(math.log2(LWE.primal_usvp(
                lwe_parameters_usvp_s, red_cost_model=RC.BDGL16)["rop"]))
            lwe_bdd_s = math.floor(math.log2(LWE.primal_bdd(
                lwe_parameters_bdd_s, red_cost_model=RC.BDGL16)["rop"]))
            lwe_usvp_numerical = math.floor(math.log2(LWE.primal_usvp(
                lwe_parameters_usvp_num, red_cost_model=RC.BDGL16)["rop"]))
            lwe_bdd_numerical = math.floor(math.log2(LWE.primal_bdd(
                lwe_parameters_bdd_num, red_cost_model=RC.BDGL16)["rop"]))

            data_point = {
                SECRET_DIST: secret, LAMBDA: l, LOG_Q: lq, USVP: est_usvp, LWE_USVP: lwe_usvp, USVP_S: est_usvp_s, LWE_USVP_S: lwe_usvp_s,
                USVP_NUM: est_usvp_numerical, LWE_NUM: lwe_usvp_numerical, BDD: est_bdd, LWE_BDD: lwe_bdd, BDD_S: est_bdd_s, LWE_BDD_S: lwe_bdd_s,
                BDD_NUM: est_bdd_numerical, LWE_NUM: lwe_bdd_numerical, OUTPUT: optimal_value, POW: closest_power_of_2(optimal_value)
            }
        else:
            data_point = {
                SECRET_DIST: secret, LAMBDA: l, LOG_Q: lq, USVP: est_usvp, USVP_S: est_usvp_s, USVP_NUM: est_usvp_numerical,
                BDD: est_bdd, BDD_S: est_bdd_s, BDD_NUM: est_bdd_numerical, OUTPUT: optimal_value, POW: closest_power_of_2(optimal_value)
            }
        data.append(data_point)

    return data


def process_logq_param(l, lwe_d, std_s, std_e, verify, estimator_installed, secret, secret_q, output_dict):
    """
    Process the parameter 'logq' and estimate its value using various models and numerical solvers.

    :param l: Security parameter.
    :param lwe_d: LWE dimension.
    :param std_s: Standard deviation of the secret.
    :param std_e: Standard deviation of the error.
    :param verify: Boolean flag to indicate if verification is needed.
    :param estimator_installed: Boolean flag to indicate if the estimator is installed.
    :param secret: Secret distribution.
    :param secret_q: Secret modulus.
    :param output_dict: Dictionary to store the output values.

    :return: List of data points with estimated values for 'logq'.
    """
    data = []
    est_usvp_numerical = int(math.floor(
        numerical_logq_usvp(l, lwe_d, std_s, std_e)))
    est_bdd_numerical = int(math.floor(
        numerical_logq_bdd(l, lwe_d, std_s, std_e)))

    optimal_value = max(est_usvp_numerical, est_bdd_numerical)
    output_dict['logq'] = optimal_value

    if verify and estimator_installed:
        lwe_parameters_bdd = LWE.Parameters(
            lwe_d, 2 ** est_bdd_numerical, ND.UniformMod(secret_q), ND.DiscreteGaussian(std_e))
        lwe_parameters_usvp = LWE.Parameters(
            lwe_d, 2 ** est_usvp_numerical, ND.UniformMod(secret_q), ND.DiscreteGaussian(std_e))
        lwe_bdd = math.floor(math.log2(LWE.primal_bdd(
            lwe_parameters_bdd, red_cost_model=RC.BDGL16)["rop"]))
        lwe_usvp = math.floor(math.log2(LWE.primal_usvp(
            lwe_parameters_usvp, red_cost_model=RC.BDGL16)["rop"]))
        data_point = {
            SECRET_DIST: secret, LAMBDA: l, LWE_DIM: lwe_d, LOG_Q: est_usvp_numerical,
            LOGQ_BDD: est_bdd_numerical, LWE_USVP: lwe_usvp, LWE_BDD: lwe_bdd, OUTPUT: optimal_value
        }
    else:
        data_point = {
            SECRET_DIST: secret, LAMBDA: l, LWE_DIM: lwe_d, LOG_Q: est_usvp_numerical,
            LOGQ_BDD: est_bdd_numerical, OUTPUT: optimal_value
        }
    data.append(data_point)

    print(data)  # Print the result
    return data


def process_logq_param_hybrid(l, lwe_d, std_s, std_e, verify, estimator_installed, secret, secret_q, h, output_dict):
    """
    Process the parameter 'logq' and estimate its value using various models and numerical solvers.

    :param l: Security parameter.
    :param lwe_d: LWE dimension.
    :param std_s: Standard deviation of the secret.
    :param std_e: Standard deviation of the error.
    :param verify: Boolean flag to indicate if verification is needed.
    :param estimator_installed: Boolean flag to indicate if the estimator is installed.
    :param secret: Secret distribution.
    :param secret_q: Secret modulus.
    :param h: Hamming weight.
    :param output_dict: Dictionary to store the output values.

    :return: List of data points with estimated values for 'logq'.
    """
    data = []
    est_hybrid = numerical_logq_hybrid(lwe_d, l, h)

    if verify and estimator_installed:
        FHEParam = LWE.Parameters(
            n=lwe_d,
            q=2**est_hybrid,
            Xs=ND.SparseTernary(lwe_d, p=h/2, m=h/2),
            Xe=ND.DiscreteGaussian(stddev=3.19)
        )
        primal_hybrid_cost = math.floor(math.log2(LWE.primal_hybrid(
            FHEParam, red_cost_model=RC.BDGL16, mitm=False)["rop"]))
        data_point = {
            SECRET_DIST: secret, LWE_DIM: lwe_d, LAMBDA: l, HW: h, HYBRID: est_hybrid, LWE_HYBRID: primal_hybrid_cost
        }
    else:
        data_point = {
            SECRET_DIST: secret, LWE_DIM: lwe_d, LAMBDA: l, HW: h, HYBRID: est_hybrid
        }

    optimal_value = est_hybrid
    output_dict['logq'] = optimal_value

    data.append(data_point)

    return data


def process_std_e_param(logq, l, lwe_d, std_s, verify, estimator_installed, secret, secret_q, output_dict):
    """
    Process the parameter 'std_e' and estimate its value using various models and numerical solvers.

    :param logq: List of log q values.
    :param l: Security parameter.
    :param lwe_d: LWE dimension.
    :param std_s: Standard deviation of the secret.
    :param verify: Boolean flag to indicate if verification is needed.
    :param estimator_installed: Boolean flag to indicate if the estimator is installed.
    :param secret: Secret distribution.
    :param secret_q: Secret modulus.
    :param output_dict: Dictionary to store the output values.

    :return: List of data points with estimated values for 'std_e'.
    """
    data = []
    for lq in logq:
        est_usvp_numerical = log2(numerical_std_e_usvp(l, lwe_d, lq, std_s))
        est_bdd_numerical = log2(numerical_std_e_bdd(l, lwe_d, lq, std_s))

        optimal_value = max(est_usvp_numerical, est_bdd_numerical)
        output_dict['std_e'] = optimal_value

        if verify and estimator_installed:
            lwe_parameters_bdd = LWE.Parameters(
                lwe_d, 2 ** lq, ND.UniformMod(secret_q), ND.DiscreteGaussian(2**est_bdd_numerical))
            lwe_parameters_usvp = LWE.Parameters(
                lwe_d, 2 ** lq, ND.UniformMod(secret_q), ND.DiscreteGaussian(2**est_usvp_numerical))
            lwe_bdd = math.floor(math.log2(LWE.primal_bdd(
                lwe_parameters_bdd, red_cost_model=RC.BDGL16)["rop"]))
            lwe_usvp = math.floor(math.log2(LWE.primal_bdd(
                lwe_parameters_usvp, red_cost_model=RC.BDGL16)["rop"]))
            data_point = {
                SECRET_DIST: secret, LAMBDA: l, LWE_DIM: lwe_d, LOG_Q: lq, STD_E_USVP: est_usvp_numerical,
                STD_E_BDD: est_bdd_numerical, LWE_USVP: lwe_usvp, LWE_BDD: lwe_bdd, OUTPUT: optimal_value
            }
        else:
            data_point = {
                SECRET_DIST: secret, LAMBDA: l, LWE_DIM: lwe_d, LOG_Q: lq,
                STD_E_USVP: est_usvp_numerical, STD_E_BDD: est_bdd_numerical, OUTPUT: optimal_value
            }
        data.append(data_point)

    print(data)  # Print the result
    return data


def process_lambda_param(logq, lwe_d, std_s, std_e, lambda_usvp, lambda_usvp_s, lambda_bdd, lambda_bdd_s, verify, estimator_installed, secret, secret_q, table, output_dict):
    data = []
    if len(logq) > 1:
        output_dict['lambda'] = []
    for lq in logq:
        est_usvp, est_bdd = estimate_usvp_bdd(
            lwe_d, lq, std_s, std_e, lambda_usvp, lambda_bdd, secret_q)
        if abs(std_e - 3.19) < 1e-9:
            est_usvp_s, est_bdd_s = estimate_usvp_s_bdd_s(
                lwe_d, lq, lambda_usvp_s, lambda_bdd_s)
            optimal_value = max(est_usvp, est_usvp_s, est_bdd, est_bdd_s)
            if verify and estimator_installed:
                lwe_parameters_bdd = LWE.Parameters(
                    lwe_d, 2 ** lq, ND.UniformMod(secret_q), ND.DiscreteGaussian(std_e))
                lwe_parameters_usvp = LWE.Parameters(
                    lwe_d, 2 ** lq, ND.UniformMod(secret_q), ND.DiscreteGaussian(std_e))
                lwe_bdd = math.floor(math.log2(LWE.primal_bdd(
                    lwe_parameters_bdd, red_cost_model=RC.BDGL16)["rop"]))
                lwe_usvp = math.floor(math.log2(LWE.primal_bdd(
                    lwe_parameters_usvp, red_cost_model=RC.BDGL16)["rop"]))
                if table:
                    data_point = {
                        SECRET_DIST: secret, LWE_DIM: lwe_d, LOG_Q: lq, USVP: est_usvp, LWE_USVP: lwe_usvp,
                        USVP_S: est_usvp_s, BDD: est_bdd, LWE_BDD: lwe_bdd, BDD_S: est_bdd_s, OUTPUT: optimal_value
                    }
                else:
                    data_point = {
                        SECRET_DIST: secret, LWE_DIM: lwe_d, LOG_Q: lq, OUTPUT: optimal_value, LWE_USVP: lwe_usvp
                    }
            else:
                if table:
                    data_point = {
                        SECRET_DIST: secret, LWE_DIM: lwe_d, LOG_Q: lq, USVP: est_usvp,
                        USVP_S: est_usvp_s, BDD: est_bdd, BDD_S: est_bdd_s, OUTPUT: optimal_value
                    }
                else:
                    data_point = {
                        SECRET_DIST: secret, LWE_DIM: lwe_d, LOG_Q: lq, OUTPUT: optimal_value
                    }
        else:
            optimal_value = max(est_usvp, est_bdd)
            if verify and estimator_installed:
                lwe_parameters_bdd = LWE.Parameters(
                    lwe_d, 2 ** lq, ND.UniformMod(secret_q), ND.DiscreteGaussian(std_e))
                lwe_parameters_usvp = LWE.Parameters(
                    lwe_d, 2 ** lq, ND.UniformMod(secret_q), ND.DiscreteGaussian(std_e))
                lwe_bdd = math.floor(math.log2(LWE.primal_bdd(
                    lwe_parameters_bdd, red_cost_model=RC.BDGL16)["rop"]))
                lwe_usvp = math.floor(math.log2(LWE.primal_bdd(
                    lwe_parameters_usvp, red_cost_model=RC.BDGL16)["rop"]))
                data_point = {
                    SECRET_DIST: secret, LWE_DIM: lwe_d, LOG_Q: lq, LWE_USVP: est_usvp, LWE_USVP_S: lwe_usvp,
                    LWE_BDD: est_bdd, LWE_BDD_S: lwe_bdd, OUTPUT: optimal_value
                }
            else:
                data_point = {
                    SECRET_DIST: secret, LWE_DIM: lwe_d, LOG_Q: lq, LWE_USVP: est_usvp,
                    LWE_BDD: est_bdd, OUTPUT: optimal_value
                }

        if len(logq) > 1:
            output_dict['lambda'].append(optimal_value)
        else:
            output_dict['lambda'] = optimal_value

        data.append(data_point)
    return data


def process_lambda_param_hybrid(logq, lwe_d, std_s, h, secret, verify, estimator_installed, output_dict):
    """
    Process the parameter 'lambda' and estimate its value using various models and numerical solvers.

    :param logq: List of log q values.
    :param lwe_d: LWE dimension.
    :param std_s: Standard deviation of the secret.
    :param h: Hamming weight.
    :param secret: Secret distribution.
    :param verify: Boolean flag to indicate if verification is needed.
    :param estimator_installed: Boolean flag to indicate if the estimator is installed.
    :param output_dict: Dictionary to store the output values.

    :return: List of data points with estimated values for 'lambda'.
    """
    data = []
    if len(logq) > 1:
        output_dict['lambda'] = []
    for lq in logq:
        est_hybrid = math.floor(numerical_lambda_hybrid_v2(lwe_d, lq, 3.19, h))

        if verify and estimator_installed:
            FHEParam = LWE.Parameters(
                n=lwe_d,
                q=2**lq,
                Xs=ND.SparseTernary(lwe_d, h/2, h/2),
                Xe=ND.DiscreteGaussian(stddev=3.19)
            )
            primal_hybrid_cost = math.floor(math.log2(LWE.primal_hybrid(
                FHEParam, red_cost_model=RC.BDGL16, mitm=False)["rop"]))
            data_point = {
                SECRET_DIST: secret, LWE_DIM: lwe_d, LOG_Q: lq, HW: h, HYBRID: est_hybrid, LWE_HYBRID: primal_hybrid_cost
            }
        else:
            data_point = {
                SECRET_DIST: secret, LWE_DIM: lwe_d, LOG_Q: lq, HW: h, HYBRID: est_hybrid
            }

        if len(logq) > 1:
            output_dict['lambda'].append(est_hybrid)
        else:
            output_dict['lambda'] = est_hybrid

        data.append(data_point)

    return data


def estimate_usvp_bdd(lwe_d, lq, std_s, std_e, lambda_usvp, lambda_bdd, secret_q):
    est_usvp = int(round(model_lambda_usvp(
        lwe_d, lq, std_s, std_e, lambda_usvp)))
    est_bdd = 0
    try:
        est_bdd = int(round(model_lambda_bdd(lwe_d, lq, std_s,
                      std_e, secret_q, lambda_bdd)[0].real))
    except Exception as e:
        print(e)
    return est_usvp, est_bdd


def estimate_usvp_s_bdd_s(lwe_d, lq, lambda_usvp_s, lambda_bdd_s):
    est_usvp_s = int(round(model_lambda_usvp_s(lwe_d, lq, lambda_usvp_s)))
    est_bdd_s = int(round(model_lambda_bdd_s(lwe_d, lq, lambda_bdd_s)))
    return est_usvp_s, est_bdd_s
