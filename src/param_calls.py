import warnings
from estimator import LWE, RC, ND
import math
from numpy import log2
import traceback

from formulas import (
    model_lambda_usvp, model_lambda_usvp_s, model_lambda_bdd, model_lambda_bdd_s,
    model_n_usvp, model_n_usvp_s, model_n_bdd, model_n_bdd_s,
)
from numerical_solver import numerical_n_usvp, numerical_n_bdd, numerical_logq_usvp, numerical_logq_bdd, numerical_std_e_usvp, numerical_std_e_bdd, numerical_lambda_bdd, numerical_lambda_usvp
from numerical_hybrid import numerical_lambda_hybrid_v2, numerical_logq_hybrid
from aux_functions import closest_power_of_2, helper, set_distribution

from const import (
    SECRET_DIST, LAMBDA, LOG_Q, USVP, LWE_USVP, LWE_USVP_F, USVP_S, LWE_USVP_S, USVP_NUM, LWE_NUM, BDD, LWE_BDD, BDD_S, LWE_BDD_S, BDD_NUM, OUTPUT, POW, LWE_DIM, LOGQ_BDD, LOGQ_USVP, LOGQ_USVP_F, HW, HYBRID, LOGQ_HYBRID, LWE_HYBRID, STD_E_USVP, STD_E_BDD, EST
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
    secret_dist = params['secret_dist']
    error_dist = params['error_dist']
    model_values = params['model_values']
    verify = params['verify']
    estimator_installed = params['estimator_installed']
    hw = params['hw']
    output_dict = params['output_dict']
    num_only = params['num_only']
    correction = params['correction']
    error_dist_tag = params['error_tag']

    if param == 'n':
        data = process_n(logq, l, error_dist, model_values['n_usvp'], model_values['n_usvp_s'],
                         model_values['n_bdd'], model_values['n_bdd_s'], verify, estimator_installed, secret_dist, table, num_only, output_dict)
    elif param == 'logq':
        data = process_logq(
            l, lwe_d, error_dist, verify, estimator_installed, correction, secret_dist, hw, output_dict)
    elif param == 'std_e':
        data = process_std_e(
            logq, l, lwe_d, verify, estimator_installed, secret_dist, error_dist, correction, error_dist_tag, table, output_dict)
    elif param == 'lambda':
        data = process_lambda(logq, lwe_d, error_dist, model_values['lambda_usvp'], model_values['lambda_usvp_s'], model_values[
            'lambda_bdd'], model_values['lambda_bdd_s'], verify, estimator_installed, secret_dist, hw, table, num_only, output_dict)
    elif param == "est":
        data = process_est(logq, lwe_d, error_dist, secret_dist)
    else:
        data = helper(), []
    return data


def process_n(logq, l, std_e, n_usvp, n_usvp_s, n_bdd, n_bdd_s, verify, estimator_installed, secret_dist, table, num_only, output_dict):

    data = process_n_param(logq, l, secret_dist, std_e, n_usvp, n_usvp_s, n_bdd,
                           n_bdd_s, verify, estimator_installed, table, num_only, output_dict)
    return data


def process_logq(l, lwe_d, std_e, verify, estimator_installed, correction, secret_dist, hw, output_dict):

    secret = secret_dist.tag

    if secret != 'SparseTernary':
        data = process_logq_param(
            l, lwe_d, std_e, verify, estimator_installed, correction, secret_dist, output_dict)
    else:
        data = process_logq_param_hybrid(
            l, lwe_d, verify, estimator_installed, secret_dist, hw, output_dict)
    return data


def process_std_e(logq, l, lwe_d, verify, estimator_installed, secret_dist, error_dist, correction, error_dist_tag, table, output_dict):
    data = process_std_e_param(
        logq, l, lwe_d, verify, estimator_installed, secret_dist, error_dist, correction, error_dist_tag, table, output_dict)
    return data


def process_lambda(logq, lwe_d, std_e, lambda_usvp, lambda_usvp_s, lambda_bdd, lambda_bdd_s, verify, estimator_installed, secret_dist, h, table, num_only, output_dict):
    secret = secret_dist.tag
    print("Secret: ", secret)
    if secret != 'SparseTernary':
        data = process_lambda_param(logq, lwe_d, std_e, lambda_usvp, lambda_usvp_s,
                                    lambda_bdd, lambda_bdd_s, verify, estimator_installed, secret_dist, table, num_only, output_dict)
    else:
        data = process_lambda_param_hybrid(
            logq, lwe_d, h, secret_dist, verify, estimator_installed, output_dict)

    return data


def process_est(logq, lwe_d, std_e, secret_q):
    print("Running Lattice Estimator")
    for lq in logq:
        parameters = LWE.Parameters(
            lwe_d, 2 ** lq, secret_dist, ND.DiscreteGaussian(std_e))
        # LWE.estimate(parameters, red_cost_model=RC.BDGL16)
        print(LWE.arora_gb(parameters))
        # print(LWE.primal_usvp(parameters))
    return [], []


def process_n_param(logq, l, secret_dist, error_dist, n_usvp, n_usvp_s, n_bdd, n_bdd_s, verify, estimator_installed, table, num_only, output_dict):
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

    std_s = float(secret_dist.stddev)
    std_e = float(error_dist.stddev)
    secret = secret_dist.tag

    if len(logq) > 1:
        output_dict['n'] = []
    for lq in logq:
        # Run the formulas for usvp and bdd together with their simplified versions and the numerical solver. Since we are interested in n we round up the result.
        est_bdd, est_usvp, est_bdd_s, est_usvp_s = 0, 0, 0, 0
        if not num_only:
            est_usvp = int(
                math.ceil(model_n_usvp(l, lq, std_s, std_e, n_usvp)))
            est_usvp_s = int(math.ceil(model_n_usvp_s(l, lq, n_usvp_s)))
            est_bdd = int(math.ceil(model_n_bdd(l, lq, std_s, std_e, n_bdd)))
            est_bdd_s = int(
                math.ceil(model_n_bdd_s(l, lq, std_s, std_e, n_bdd_s)))
            # Store the minimum value of n provided from all the formulas

        est_usvp_numerical = int(
            math.ceil(numerical_n_usvp(l, lq, std_s, std_e)))
        est_bdd_numerical = int(
            math.ceil(numerical_n_bdd(l, lq, std_s, std_e)))

        return_value = max(est_usvp, est_usvp_s, est_bdd,
                           est_bdd_s, est_usvp_numerical, est_bdd_numerical)

        if len(logq) > 1:
            output_dict['n'].append(return_value)
        else:
            output_dict['n'] = return_value

        if verify and estimator_installed:
            lwe_bdd, lwe_usvp, lwe_bdd_s, lwe_usvp_s = 0, 0, 0, 0
            if not num_only:
                lwe_parameters_usvp = LWE.Parameters(
                    est_usvp, 2 ** lq, secret_dist, ND.DiscreteGaussian(std_e))
                lwe_parameters_bdd = LWE.Parameters(
                    est_bdd, 2 ** lq, secret_dist, ND.DiscreteGaussian(std_e))
                lwe_parameters_usvp_s = LWE.Parameters(
                    est_usvp_s, 2 ** lq, secret_dist, ND.DiscreteGaussian(std_e))
                lwe_parameters_bdd_s = LWE.Parameters(
                    est_bdd_s, 2 ** lq, secret_dist, ND.DiscreteGaussian(std_e))
                lwe_usvp = math.floor(math.log2(LWE.primal_usvp(
                    lwe_parameters_usvp, red_cost_model=RC.BDGL16)["rop"]))
                lwe_bdd = math.floor(math.log2(LWE.primal_bdd(
                    lwe_parameters_bdd, red_cost_model=RC.BDGL16)["rop"]))
                lwe_usvp_s = math.floor(math.log2(LWE.primal_usvp(
                    lwe_parameters_usvp_s, red_cost_model=RC.BDGL16)["rop"]))
                lwe_bdd_s = math.floor(math.log2(LWE.primal_bdd(
                    lwe_parameters_bdd_s, red_cost_model=RC.BDGL16)["rop"]))

            lwe_parameters_usvp_num = LWE.Parameters(
                est_usvp_numerical, 2 ** lq, secret_dist, ND.DiscreteGaussian(std_e))
            lwe_parameters_bdd_num = LWE.Parameters(
                est_bdd_numerical, 2 ** lq, secret_dist, ND.DiscreteGaussian(std_e))

            lwe_usvp_numerical = math.floor(math.log2(LWE.primal_usvp(
                lwe_parameters_usvp_num, red_cost_model=RC.BDGL16)["rop"]))
            lwe_bdd_numerical = math.floor(math.log2(LWE.primal_bdd(
                lwe_parameters_bdd_num, red_cost_model=RC.BDGL16)["rop"]))

            estimates = {
                est_usvp: lwe_usvp, est_usvp_s: lwe_usvp_s,
                est_usvp_numerical: lwe_usvp_numerical, est_bdd: lwe_bdd,
                est_bdd_s: lwe_bdd_s, est_bdd_numerical: lwe_bdd_numerical
            }

            if table:
                if not num_only:
                    data_point = {
                        SECRET_DIST: secret, LAMBDA: l, LOG_Q: lq, USVP: est_usvp, LWE_USVP: lwe_usvp, USVP_S: est_usvp_s, LWE_USVP_S: lwe_usvp_s,
                        USVP_NUM: est_usvp_numerical, LWE_NUM: lwe_usvp_numerical, BDD: est_bdd, LWE_BDD: lwe_bdd, BDD_S: est_bdd_s, LWE_BDD_S: lwe_bdd_s,
                        BDD_NUM: est_bdd_numerical, LWE_NUM: lwe_bdd_numerical, OUTPUT: return_value, POW: closest_power_of_2(return_value)
                    }
                else:
                    data_point = {
                        SECRET_DIST: secret, LAMBDA: l, LOG_Q: lq, USVP_NUM: est_usvp_numerical, LWE_USVP: lwe_usvp_numerical,
                        BDD_NUM: est_bdd_numerical, LWE_BDD: lwe_bdd_numerical, OUTPUT: return_value, POW: closest_power_of_2(return_value)
                    }
            else:
                data_point = {
                    SECRET_DIST: secret, LAMBDA: l, LOG_Q: lq, OUTPUT: return_value, EST: estimates[return_value], POW: closest_power_of_2(return_value)
                }
        else:
            if table:
                if not num_only:
                    data_point = {
                        SECRET_DIST: secret, LAMBDA: l, LOG_Q: lq, USVP: est_usvp, USVP_S: est_usvp_s, USVP_NUM: est_usvp_numerical,
                        BDD: est_bdd, BDD_S: est_bdd_s, BDD_NUM: est_bdd_numerical, OUTPUT: return_value, POW: closest_power_of_2(return_value)
                    }
                else:
                    data_point = {
                        SECRET_DIST: secret, LAMBDA: l, LOG_Q: lq, USVP_NUM: est_usvp_numerical, BDD_NUM: est_bdd_numerical, OUTPUT: return_value, POW: closest_power_of_2(return_value)
                    }
            else:
                data_point = {
                    SECRET_DIST: secret, LAMBDA: l, LOG_Q: lq, OUTPUT: return_value, POW: closest_power_of_2(return_value)
                }
        data.append(data_point)

    return data


def get_parameters(lwe_d, lnq, secret_dist, error_dist, est_usvp_numerical, est_bdd_numerical, error_dist_tag, param):

    lwe_parameters_usvp = None
    lwe_parameters_bdd = None

    if param == 'logq':
        lwe_parameters_usvp = LWE.Parameters(
            lwe_d, 2 ** est_usvp_numerical, secret_dist, error_dist)
        lwe_parameters_bdd = LWE.Parameters(
            lwe_d, 2 ** est_bdd_numerical, secret_dist, error_dist)
    if param == 'std_e':
        # TODO: numerics for std only works for gaussian distribution right?
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

    return max(corrected_logq_bdd, corrected_logq_usvp), corrected_logq_bdd, corrected_logq_usvp, corrected_lwe_bdd, corrected_lwe_usvp


def process_logq_param(l, lwe_d, error_dist, verify, estimator_installed, correction, secret_dist, output_dict):
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
    std_s = secret_dist.stddev
    secret = secret_dist.tag
    std_e = error_dist.stddev

    est_usvp_numerical = int(math.floor(
        numerical_logq_usvp(l, lwe_d, std_s, std_e)))
    est_bdd_numerical = int(math.floor(
        numerical_logq_bdd(l, lwe_d, std_s, std_e)))

    return_value = max(est_usvp_numerical, est_bdd_numerical)
    output_dict['logq'] = return_value

    if verify and estimator_installed:
        lwe_parameters_bdd = LWE.Parameters(
            lwe_d, 2 ** est_bdd_numerical, secret_dist, error_dist)
        lwe_parameters_usvp = LWE.Parameters(
            lwe_d, 2 ** est_usvp_numerical, secret_dist, error_dist)
        lwe_bdd = math.floor(math.log2(LWE.primal_bdd(
            lwe_parameters_bdd, red_cost_model=RC.BDGL16)["rop"]))
        lwe_usvp = math.floor(math.log2(LWE.primal_usvp(
            lwe_parameters_usvp, red_cost_model=RC.BDGL16)["rop"]))

        corrected_logq_bdd, corrected_logq_usvp, corrected_lwe_bdd, corrected_lwe_usvp = est_usvp_numerical, est_bdd_numerical, lwe_bdd, lwe_usvp

        num_calls = 2

        if correction:
            return_value, corrected_logq_bdd, corrected_logq_usvp, corrected_lwe_bdd, corrected_lwe_usvp = correction_logic(
                l, lwe_d, None, lwe_usvp, lwe_bdd, secret_dist, error_dist, est_usvp_numerical, est_bdd_numerical, 'logq', num_calls)

        data_point = {
            SECRET_DIST: secret, LAMBDA: l, LWE_DIM: lwe_d, LOGQ_USVP: corrected_logq_usvp,
            LWE_USVP: corrected_lwe_usvp, LOGQ_BDD: corrected_logq_bdd, LWE_BDD: corrected_lwe_bdd, OUTPUT: return_value
        }
    else:
        data_point = {
            SECRET_DIST: secret, LAMBDA: l, LWE_DIM: lwe_d, LOGQ_USVP: est_usvp_numerical,
            LOGQ_BDD: est_bdd_numerical, OUTPUT: return_value
        }
    data.append(data_point)

    print(data)  # Print the result
    return data


def process_logq_param_hybrid(l, lwe_d, verify, estimator_installed, secret_dist, h, output_dict):
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
    secret = secret_dist.tag
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
            SECRET_DIST: secret, LWE_DIM: lwe_d, LAMBDA: l, HW: h, LOGQ_HYBRID: est_hybrid, LWE_HYBRID: primal_hybrid_cost
        }
    else:
        data_point = {
            SECRET_DIST: secret, LWE_DIM: lwe_d, LAMBDA: l, HW: h, LOGQ_HYBRID: est_hybrid
        }

    return_value = est_hybrid
    output_dict['logq'] = return_value

    data.append(data_point)

    return data


def process_std_e_param(logq, l, lwe_d, verify, estimator_installed, secret_dist, error_dist, correction, error_dist_tag, table, output_dict):
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
    std_s = secret_dist.stddev
    secret = secret_dist.tag

    for lq in logq:

        lwe_bdd_minimize, lwe_usvp_minimize, lwe_bdd, lwe_usvp = 0, 0, 0, 0
        est_usvp_numerical, est_bdd_numerical, est_bdd_numerical_minimize, est_usvp_numerical_minimize = 1, 1, 1, 1

        try:
            est_usvp_numerical, est_usvp_numerical_status = numerical_std_e_usvp(
                l, lwe_d, lq, std_s)
            est_usvp_numerical = log2(est_usvp_numerical)
        except Exception as e:
            print(f"Error in numerical_std_e_usvp: {e}")
            traceback.print_exc()  # Print the full traceback, including the line number
            est_usvp_numerical, est_usvp_numerical_status = 0, False

        try:
            est_bdd_numerical, est_bdd_numerical_status = numerical_std_e_bdd(
                l, lwe_d, lq, std_s)
            est_bdd_numerical = log2(est_bdd_numerical)
        except Exception as e:
            print(f"Error in numerical_std_e_bdd: {e}")
            traceback.print_exc()  # Print the full traceback, including the line number
            est_bdd_numerical, est_bdd_numerical_status = 0, False

        # if est_bdd_numerical_status:
        #     lwe_parameters_bdd = LWE.Parameters(
        #         lwe_d, 2 ** lq, secret_dist, ND.DiscreteGaussian(2**est_bdd_numerical))
        #     lwe_bdd = math.floor(math.log2(LWE.primal_bdd(
        #         lwe_parameters_bdd, red_cost_model=RC.BDGL16)["rop"]))

        # if est_usvp_numerical_status:
        #     lwe_parameters_usvp = LWE.Parameters(
        #         lwe_d, 2 ** lq, secret_dist, ND.DiscreteGaussian(2**est_usvp_numerical))
        #     lwe_usvp = math.floor(math.log2(LWE.primal_usvp(
        #         lwe_parameters_usvp, red_cost_model=RC.BDGL16)["rop"]))

        return_value = max(est_usvp_numerical, est_bdd_numerical)

        # if correction:

        #     if lwe_usvp >= l:
        #         while (lwe_usvp >= l):
        #             corrected_lwe_usvp = lwe_usvp
        #             corrected_logq_usvp = est_usvp_numerical
        #             print("Applying correction lwe usvp > l", "logq ",
        #                   corrected_logq_usvp, " est", corrected_lwe_usvp)
        #             est_usvp_numerical += 1
        #             num_calls += 1
        #             lwe_parameters_usvp = LWE.Parameters(
        #                 lwe_d, 2 ** est_usvp_numerical, secret_dist, ND.DiscreteGaussian(std_e))
        #             lwe_usvp = math.floor(math.log2(LWE.primal_usvp(
        #                 lwe_parameters_usvp, red_cost_model=RC.BDGL16)["rop"]))
        #     else:
        #         while (lwe_usvp < l):
        #             est_usvp_numerical -= 1
        #             num_calls += 1
        #             lwe_parameters_usvp = LWE.Parameters(
        #                 lwe_d, 2 ** est_usvp_numerical, secret_dist, ND.DiscreteGaussian(std_e))
        #             lwe_usvp = math.floor(math.log2(LWE.primal_usvp(
        #                 lwe_parameters_usvp, red_cost_model=RC.BDGL16)["rop"]))
        #             corrected_lwe_usvp = lwe_usvp
        #             corrected_logq_usvp = est_usvp_numerical
        #             print("Applying correction lwe usvp < l", "logq ",
        #                   corrected_logq_usvp, " est", corrected_lwe_usvp)

        #     if lwe_bdd >= l:
        #         while (lwe_bdd >= l):
        #             corrected_lwe_bdd = lwe_bdd
        #             corrected_logq_bdd = est_bdd_numerical
        #             print("Applying correction lwe bdd > l", "logq ",
        #                   corrected_logq_bdd, " est", corrected_lwe_bdd)
        #             est_bdd_numerical += 1
        #             num_calls += 1
        #             lwe_parameters_bdd = LWE.Parameters(
        #                 lwe_d, 2 ** est_bdd_numerical, secret_dist, ND.DiscreteGaussian(std_e))
        #             lwe_bdd = math.floor(math.log2(LWE.primal_usvp(
        #                 lwe_parameters_bdd, red_cost_model=RC.BDGL16)["rop"]))
        #     else:
        #         while (lwe_bdd < l):
        #             est_bdd_numerical -= 1
        #             num_calls += 1
        #             lwe_parameters_bdd = LWE.Parameters(
        #                 lwe_d, 2 ** est_bdd_numerical, secret_dist, ND.DiscreteGaussian(std_e))
        #             lwe_bdd = math.floor(math.log2(LWE.primal_usvp(
        #                 lwe_parameters_bdd, red_cost_model=RC.BDGL16)["rop"]))
        #             corrected_lwe_bdd = lwe_bdd
        #             corrected_logq_bdd = est_bdd_numerical
        #             print("Applying correction lwe bdd < l", "logq ",
        #                   corrected_logq_bdd, " est", corrected_lwe_bdd)

        #     print("Number of calls to the estimator: ", num_calls)
        #     return_value = max(corrected_logq_bdd, corrected_logq_usvp)

        output_dict['std_e'] = return_value

        if verify and estimator_installed:
            if est_bdd_numerical_status:
                lwe_parameters_bdd = LWE.Parameters(
                    lwe_d, 2 ** lq, secret_dist, set_distribution(error_dist_tag, {'std': 2**est_bdd_numerical}))
                lwe_bdd = math.floor(math.log2(LWE.primal_bdd(
                    lwe_parameters_bdd, red_cost_model=RC.BDGL16)["rop"]))

            if est_usvp_numerical_status:
                lwe_parameters_usvp = LWE.Parameters(
                    lwe_d, 2 ** lq, secret_dist, set_distribution(error_dist_tag, {'std': 2**est_bdd_numerical}))
                lwe_usvp = math.floor(math.log2(LWE.primal_usvp(
                    lwe_parameters_usvp, red_cost_model=RC.BDGL16)["rop"]))

            estimates = {
                est_usvp_numerical: lwe_usvp, est_bdd_numerical: lwe_bdd
            }

            corrected_std_e_bdd, corrected_std_e_usvp, corrected_lwe_bdd, corrected_lwe_usvp = est_usvp_numerical, est_bdd_numerical, lwe_bdd, lwe_usvp

            num_calls = 2

            if correction:
                return_value, corrected_std_e_usvp, corrected_std_e_bdd, corrected_lwe_bdd, corrected_lwe_usvp = correction_logic(
                    l, lwe_d, lq, lwe_usvp, lwe_bdd, secret_dist, error_dist, est_usvp_numerical, est_bdd_numerical, 'std_e', num_calls, error_dist_tag)

            if table:
                data_point = {
                    SECRET_DIST: secret, LAMBDA: l, LWE_DIM: lwe_d, LOG_Q: lq,
                    STD_E_USVP: corrected_std_e_usvp, LWE_USVP: corrected_lwe_usvp, STD_E_BDD: corrected_std_e_bdd, LWE_BDD: corrected_lwe_bdd, OUTPUT: return_value
                }
            else:
                data_point = {
                    SECRET_DIST: secret, LAMBDA: l, LWE_DIM: lwe_d, LOG_Q: lq, OUTPUT: return_value, EST: estimates[return_value]
                }
        else:
            if table:
                data_point = {
                    SECRET_DIST: secret, LAMBDA: l, LWE_DIM: lwe_d, LOG_Q: lq,
                    STD_E_USVP: est_usvp_numerical, STD_E_BDD: est_bdd_numerical, OUTPUT: return_value
                }
            else:
                data_point = {
                    SECRET_DIST: secret, LAMBDA: l, LWE_DIM: lwe_d, LOG_Q: lq, OUTPUT: return_value
                }
        data.append(data_point)

    return data


def process_lambda_param(logq, lwe_d, std_e, lambda_usvp, lambda_usvp_s, lambda_bdd, lambda_bdd_s, verify, estimator_installed, secret_dist, table, num_only, output_dict):
    data = []
    if len(logq) > 1:
        output_dict['lambda'] = []

    for lq in logq:
        data_point = process_lambda_for_lq(lq, lwe_d, std_e, lambda_usvp, lambda_usvp_s,
                                           lambda_bdd, lambda_bdd_s, verify, estimator_installed, secret_dist, table, num_only)
        return_value = data_point[OUTPUT]
        if len(logq) > 1:
            output_dict['lambda'].append(return_value)
        else:
            output_dict['lambda'] = return_value
        data.append(data_point)

    return data


def process_lambda_for_lq(lq, lwe_d, error_dist, lambda_usvp, lambda_usvp_s, lambda_bdd, lambda_bdd_s, verify, estimator_installed, secret_dist, table, num_only):
    if not num_only:
        est_usvp, est_bdd = estimate_usvp_bdd(
            lwe_d, lq, error_dist, lambda_usvp, lambda_bdd, secret_dist)

    std_s = float(secret_dist.stddev)
    std_e = float(error_dist.stddev)

    est_num_bdd = math.floor(
        numerical_lambda_bdd(lwe_d, lq, std_s, std_e))
    est_num_usvp = math.floor(
        numerical_lambda_usvp(lwe_d, lq, std_s, std_e))

    if abs(std_e - 3.19) < 1e-9:
        if not num_only:
            est_usvp_s, est_bdd_s = estimate_usvp_s_bdd_s(
                lwe_d, lq, lambda_usvp_s, lambda_bdd_s)
            return_value = min(max(est_usvp, est_usvp_s),
                               max(est_bdd, est_bdd_s))
        else:
            return_value = min(est_num_bdd, est_num_usvp)
            est_usvp_s = None
            est_bdd_s = None
            est_bdd = None
            est_usvp = None
        data_point = create_data_point(lq, lwe_d, std_e, secret_dist, est_usvp,
                                       est_usvp_s, est_bdd, est_bdd_s, est_num_bdd, est_num_usvp, return_value, verify, estimator_installed, table, num_only)
    else:
        if not num_only:
            return_value = min(est_usvp, est_bdd)
        else:
            return_value = min(est_num_bdd, est_num_usvp)
            est_usvp_s = None
            est_bdd_s = None
            est_bdd = None
            est_usvp = None
        data_point = create_data_point(lq, lwe_d, std_e, secret_dist, est_usvp,
                                       None, est_bdd, None, est_num_bdd, est_num_usvp, return_value, verify, estimator_installed, table, num_only)

    return data_point


def create_data_point(lq, lwe_d, std_e, secret_dist, est_usvp, est_usvp_s, est_bdd, est_bdd_s, est_num_bdd, est_num_usvp, return_value, verify, estimator_installed, table, num_only):

    secret = secret_dist.tag

    if verify and estimator_installed:
        lwe_parameters = LWE.Parameters(
            lwe_d, 2 ** lq, secret_dist, ND.DiscreteGaussian(std_e))
        lwe_bdd = math.floor(math.log2(LWE.primal_bdd(
            lwe_parameters, red_cost_model=RC.BDGL16)["rop"]))
        lwe_usvp = math.floor(math.log2(LWE.primal_usvp(
            lwe_parameters, red_cost_model=RC.BDGL16)["rop"]))

        if not num_only:
            estimates = {
                est_usvp: lwe_usvp, est_bdd: lwe_bdd, est_usvp_s: lwe_usvp, est_bdd_s: lwe_bdd, est_num_bdd: lwe_bdd, est_num_usvp: lwe_usvp
            }
        else:
            estimates = {
                est_num_bdd: lwe_bdd, est_num_usvp: lwe_usvp
            }

        if table:
            if not num_only:
                data_point = {
                    SECRET_DIST: secret, LWE_DIM: lwe_d, LOG_Q: lq, USVP: est_usvp, LWE_USVP: lwe_usvp,
                    USVP_S: est_usvp_s, USVP_NUM: est_num_usvp, BDD: est_bdd, LWE_BDD: lwe_bdd, BDD_S: est_bdd_s, BDD_NUM: est_num_bdd, OUTPUT: return_value
                }
            else:
                data_point = {
                    SECRET_DIST: secret, LWE_DIM: lwe_d, LOG_Q: lq, USVP_NUM: est_num_usvp, LWE_USVP: lwe_usvp, BDD_NUM: est_num_bdd, LWE_BDD: lwe_bdd,
                    OUTPUT: return_value
                }
        else:
            data_point = {
                SECRET_DIST: secret, LWE_DIM: lwe_d, LOG_Q: lq, OUTPUT: return_value, EST: estimates[return_value]
            }
    else:
        if table:
            if not num_only:
                data_point = {
                    SECRET_DIST: secret, LWE_DIM: lwe_d, LOG_Q: lq, USVP: est_usvp,
                    USVP_S: est_usvp_s, USVP_NUM: est_num_usvp, BDD: est_bdd, BDD_S: est_bdd_s, BDD_NUM: est_num_bdd, OUTPUT: return_value
                }
            else:
                data_point = {
                    SECRET_DIST: secret, LWE_DIM: lwe_d, LOG_Q: lq, USVP_NUM: est_num_usvp, BDD_NUM: est_num_bdd, OUTPUT: return_value
                }
        else:
            data_point = {
                SECRET_DIST: secret, LWE_DIM: lwe_d, LOG_Q: lq, OUTPUT: return_value
            }

    return data_point


def process_lambda_param_hybrid(logq, lwe_d, h, secret_dist, verify, estimator_installed, output_dict):
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

    secret = secret_dist.tag

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


def estimate_usvp_bdd(lwe_d, lq, error_dist, lambda_usvp, lambda_bdd, secret_dist):

    std_s = float(secret_dist.stddev)
    std_e = float(error_dist.stddev)

    est_usvp = int(round(model_lambda_usvp(
        lwe_d, lq, std_s, std_e, lambda_usvp)))
    est_bdd = 0
    try:
        est_bdd = int(round(model_lambda_bdd(
            lwe_d, lq, std_s, std_e, lambda_bdd)[0].real))
    except Exception as e:
        print(e)
    return est_usvp, est_bdd


def estimate_usvp_s_bdd_s(lwe_d, lq, lambda_usvp_s, lambda_bdd_s):
    est_usvp_s = int(round(model_lambda_usvp_s(lwe_d, lq, lambda_usvp_s)))
    est_bdd_s = int(round(model_lambda_bdd_s(lwe_d, lq, lambda_bdd_s)))
    return est_usvp_s, est_bdd_s
