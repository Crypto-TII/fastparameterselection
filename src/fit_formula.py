from nd import NoiseDistribution as ND
import lmfit
import numpy as np
from aux_functions import helper_fit
import pickle

import sys

from formulas import (
    model_n_bdd_rev1, model_n_bdd_s, model_n_usvp, model_n_usvp_s, model_lambda_usvp, model_lambda_bdd, model_lambda_bdd_s, model_lambda_usvp_s
)

#!/usr/bin/env python3


secret = "binary"
param = None
lwe_d = None
logQ = 64
name_file = None
verify = 0
lwe_parameters = []
std_s = 0.5
std_e = 3.19
secret_q = 2
ntru_flag = False
simpl = 0
ymin = 1000
ymax = 2050
xmin = 20
xmax = 65
degrees = []
attack = 'usvp'

headers = []
data = []

params_list = [('P0', 0.07, '-inf'), ('P1', 0.34, '-inf'), ('P2', 1, '-inf'),
               ('P3', 1, '-inf'), ('P4', 1, '-inf')]

security_levels = [145, 140, 135, 128, 125,
                   120, 115, 110, 105, 100, 95, 90, 85, 80]

degrees = [1024, 2048]

verbose = 0

# AUX FUNCTIONS FOR FITTING


def degree_fit(params, logq, levels, e_std, s_std, std_s_num):
    new_params = [params[p[0]] for p in params_list]

    func_map = {
        ('n', '0', 'bdd'): model_n_bdd_rev1,
        ('n', '0', 'usvp'): model_n_usvp,
        ('n', '1', 'bdd'): model_n_bdd_s,
        ('n', '1', 'usvp'): model_n_usvp_s,
        ('lambda', '0', 'usvp'): model_lambda_usvp,
        ('lambda', '0', 'bdd'): model_lambda_bdd,
        ('lambda', '1', 'usvp'): model_lambda_usvp_s,
        ('lambda', '1', 'bdd'): model_lambda_bdd_s,
    }

    residues = np.array([])

    simpl_key = str(simpl)

    for i, level in enumerate(levels):
        key = (param, simpl_key, attack)
        func = func_map.get(key)

        if func is None:
            model = None
        else:
            x = security_levels[i] if param == 'n' else degrees[i]

            if simpl_key == '1':
                model = func(x, logq, new_params)
            else:
                model = func(x, logq, s_std, e_std, new_params)

        if level is not None and model is not None:
            residues = np.concatenate((residues, level - model))

    return residues


def fit_formula(points_est, e_std, s_std, std_s_num, params):

    args = {'logq': logQ}

    args['levels'] = points_est
    args['e_std'] = e_std
    args['s_std'] = s_std
    args['std_s_num'] = std_s_num

    if (verbose):
        print(points_est)

    fit = []

    fit = lmfit.minimize(degree_fit, params, nan_policy='omit', kws=args)

    fit_results = []

    for p in params_list:
        res = fit.params[p[0]].value
        fit_results.append(res)

    return fit_results


def filter_points(d_estimates):

    d_estimates_filtered = []

    if verbose:
        print("d_estimates unfiltered: ", d_estimates)

    for d in d_estimates[1:]:

        d = sorted(d, key=lambda tup: (tup[0], -tup[1]))

        filtered_list = [tup for i, tup in enumerate(d) if i == 0 or (
            tup[0] != d[i-1][0] or tup[1] > d[i-1][1])]
        filtered_list2 = [tup for i, tup in enumerate(
            filtered_list) if i == 0 or tup[1] != filtered_list[i-1][1]]
        d_estimates_filtered.append(filtered_list2)

    return d_estimates_filtered


def _parse_options(options):
    std_e = 3.19
    std_s = 0.5
    secret = "binary"
    secret_q = 2
    name_file = None
    param_val = None
    simpl_val = 0
    attack_val = None

    def _handle_secret(arg):
        nonlocal secret, std_s, secret_q
        secret = arg
        if secret == 'binary':
            std_s = float(ND.UniformMod(2).stddev)
            secret_q = 2
        elif secret == 'ternary':
            std_s = float(ND.UniformMod(3).stddev)
            secret_q = 3
        else:
            print("Secret distribution not supported")
            sys.exit()

    handlers = {
        '--attack': ('attack_val', lambda a: a),
        '--error': ('std_e', lambda a: float(a)),
        '--param': ('param_val', lambda a: a),
        '--file': ('name_file', lambda a: a),
        '--simpl': ('simpl_val', lambda a: a),
    }
    results = {}

    for opt, arg in options:
        if opt == '-h':
            helper_fit()
            continue
        if opt == '--secret':
            _handle_secret(arg)
            continue
        if opt == '--fit':
            continue
        if opt in handlers:
            key, fn = handlers[opt]
            results[key] = fn(arg)
            continue
        helper_fit()

    attack_val = results.get('attack_val', attack_val)
    std_e = results.get('std_e', std_e)
    param_val = results.get('param_val', param_val)
    name_file = results.get('name_file', name_file)
    simpl_val = results.get('simpl_val', simpl_val)

    return std_e, std_s, secret, secret_q, name_file, param_val, simpl_val, attack_val


def _configure_ranges(secret, attack, param_value):
    def _logq_large():
        return list(range(10, 200)) + list(range(200, 500, 10)) + list(range(500, 1000, 10)) + list(range(1000, 1600, 50))

    dispatch = {
        ('n', 'binary'): {
            'degrees': list(range(2**10, 2**12, 32)),
            'logQ': list(range(20, 65)),
            'ymin': 1000, 'ymax': 2050, 'xmin': 20, 'xmax': 65,
            'bdd': ('dataset/data_binary_bdd.pkl', [('P0', 0.07, '-inf'), ('P1', 0.34, '-inf'), ('P2', 1, '-inf'), ('P3', 1, '-inf'), ('P4', 1, '-inf')]),
            'usvp': ('dataset/data_binary_usvp.pkl', [('P0', 0.07, '-inf'), ('P1', 0.34, 0), ('P2', 1, 0), ('P3', 1, '-inf'), ('P4', 1, '-inf')]),
        },
        ('n', 'ternary'): {
            'degrees': list(range(2**10, 2**15, 2**2)),
            'logQ': _logq_large(),
            'ymin': 1024, 'ymax': 2048, 'xmin': 10, 'xmax': 90,
            'bdd': ('dataset/data_ternary_bdd.pkl', [('P0', 0.07, '-inf'), ('P1', 0.34, '-inf'), ('P2', 1, '-inf'), ('P3', 1, '-inf'), ('P4', 1, '-inf')]),
            'usvp': ('dataset/data_ternary_usvp.pkl', [('P0', 0.07, '-inf'), ('P1', 0.34, 0), ('P2', 1, 0.931202), ('P3', 1, '-inf'), ('P4', 1, '-inf')]),
        },
        ('lambda', 'binary'): {
            'degrees': [2**10, 2**11],
            'logQ': list(range(20, 65)),
            'ymin': 80, 'ymax': 256, 'xmin': 10, 'xmax': 90,
            'bdd': ('dataset/data_binary_bdd.pkl', [('P0', 0.07, 0), ('P1', 0.34, '-inf'), ('P2', 1, '-inf'), ('P3', 1, '-inf'), ('P4', 1, '-inf')]),
            'usvp': ('dataset/data_binary_usvp.pkl', [('P0', 0.07, '-inf'), ('P1', 0.34, '-inf'), ('P2', 1, '-inf'), ('P3', 1, '-inf'), ('P4', 1, '-inf')]),
        },
        ('lambda', 'ternary'): {
            'degrees': [2**10, 2**11, 2032, 2**13, 9216, 10240, 11264, 12609, 13633,
                        14657, 15681, 16384, 17408, 19456, 22528, 24194, 25218, 28290, 32386, 32768],
            'logQ': _logq_large(),
            'ymin': 80, 'ymax': 256, 'xmin': 20, 'xmax': 1400,
            'bdd': ('dataset/data_ternary.pkl', [('P0', 0.07, '-inf'), ('P1', 0.34, '-inf'), ('P2', 1, '-inf'), ('P3', 1, '-inf'), ('P4', 1, '-inf')]),
            'usvp': ('dataset/data_ternary.pkl', [('P0', 0.07, '-inf'), ('P1', 0.34, '-inf'), ('P2', 1, '-inf'), ('P3', 1, '-inf'), ('P4', 1, '-inf')]),
        },
    }
    if param_value is None:
        cfg = None
    else:
        key = (param_value, 'binary' if secret == 'binary' else 'ternary')
        cfg = dispatch.get(key)

    if cfg is None:
        degrees_local = globals().get('degrees', degrees)
        logQ_local = globals().get('logQ', logQ)
        ymin_local = globals().get('ymin', ymin)
        ymax_local = globals().get('ymax', ymax)
        xmin_local = globals().get('xmin', xmin)
        xmax_local = globals().get('xmax', xmax)
        name_file = None
        params_list_local = None
        return name_file, params_list_local, degrees_local, logQ_local, ymin_local, ymax_local, xmin_local, xmax_local

    attack_entry = cfg.get(attack)
    if attack_entry is None:
        print("Attack " + str(attack) + " not considered")
        exit(1)

    name_file, params_list_local = attack_entry
    degrees_local = cfg['degrees']
    logQ_local = cfg['logQ']
    ymin_local = cfg['ymin']
    ymax_local = cfg['ymax']
    xmin_local = cfg['xmin']
    xmax_local = cfg['xmax']

    return name_file, params_list_local, degrees_local, logQ_local, ymin_local, ymax_local, xmin_local, xmax_local


def _load_and_filter(name_file, secret, verbose):
    with open(name_file, 'rb') as file:
        d_estimates = pickle.load(file)
    if verbose:
        print("file name", name_file)
        print("d_estimates", d_estimates)
    if secret == 'tfhe':
        d_estimates_filtered = [d_estimates]
        degrees_local = []
        for d in d_estimates:
            degrees_local.append(d[0])
        return d_estimates_filtered, degrees_local
    else:
        d_estimates_filtered = filter_points(d_estimates)
        return d_estimates_filtered, None


def _build_params(params_list_src):
    params = lmfit.Parameters()
    for p in params_list_src:
        if p[2] != '-inf':
            params.add(p[0], value=p[1], min=p[2])
        else:
            params.add(p[0], value=p[1])
    return params


def _build_point_grids(d_estimates_filtered, param_value, attack_value, degrees_val, logQ_val, verbose_flag):
    def _init_grids(rows, cols):
        row_template = [np.nan] * cols
        return [row_template.copy() for _ in range(rows)], [row_template.copy() for _ in range(rows)], [row_template.copy() for _ in range(rows)]

    def _populate_grid_n():
        rows = len(security_levels)
        pe, pa, ps = _init_grids(rows, len(logQ_val))
        for i in range(rows):
            tuples = d_estimates_filtered[i]
            for tup in tuples:
                if verbose_flag:
                    print(tup)
                atk = tup[3]
                match_attack = (attack_value == 'bdd' and atk == 1) or (
                    attack_value == 'usvp' and atk == 0)
                if not match_attack:
                    continue
                try:
                    col = logQ_val.index(tup[1])
                except ValueError:
                    continue
                pe[i][col] = tup[0]
                pa[i][col] = atk
        return pe, pa, ps

    def _populate_grid_lambda():
        rows = len(degrees_val)
        pe, pa, ps = _init_grids(rows, len(logQ_val))
        deg_index = {d: i for i, d in enumerate(degrees_val)}
        for filtered in d_estimates_filtered:
            for tup in filtered:
                if verbose_flag:
                    print("tup: ", tup)
                d = tup[0]
                i = deg_index.get(d)
                if i is None:
                    continue
                atk = tup[3]
                match_attack = (attack_value == 'usvp' and atk == 0) or (
                    attack_value == 'bdd' and atk == 1)
                if not match_attack:
                    continue
                try:
                    col = logQ_val.index(tup[1])
                except ValueError:
                    continue
                pe[i][col] = tup[2]
                pa[i][col] = atk
        return pe, pa, ps

    if param_value == 'n':
        points_est, points_atk, points_secret_dist = _populate_grid_n()
    elif param_value == 'lambda':
        points_est, points_atk, points_secret_dist = _populate_grid_lambda()
    else:
        points_est, points_atk, points_secret_dist = [], [], []

    if verbose_flag:
        print("points_est", points_est)

    return points_est, points_atk, points_secret_dist

# MAIN FUNCTION TO FIND CONSTANTS


def find_constants(opts):

    # estimator_installed = check_estimator_installed()
    # if not estimator_installed:
    #     print("Lattice Estimator not installed, can't find constants")
    #     exit(0)

    std_e, std_s, secret, secret_q, name_file, param_val, simpl_val, attack_val = _parse_options(
        opts)

    if attack_val is not None:
        globals()['attack'] = attack_val
    if param_val is not None:
        globals()['param'] = param_val
    if simpl_val is not None:
        globals()['simpl'] = simpl_val

    param_in_use = globals().get('param', param)
    attack_in_use = globals().get('attack', attack)

    name_file_conf, params_list_local, degrees_local, logQ_local, ymin_local, ymax_local, xmin_local, xmax_local = _configure_ranges(
        secret, attack_in_use, param_in_use)

    globals().update({'degrees': degrees_local, 'logQ': logQ_local,
                      'ymin': ymin_local, 'ymax': ymax_local, 'xmin': xmin_local, 'xmax': xmax_local})

    if name_file is not None:
        name_file_conf = name_file

    d_estimates_filtered, tfhe_degrees = _load_and_filter(
        name_file_conf, secret, verbose)
    if tfhe_degrees:
        globals()['degrees'] = tfhe_degrees
        degrees_local = tfhe_degrees

    params_src = params_list_local if params_list_local is not None else params_list
    params = _build_params(params_src)

    points_est, points_atk, points_secret_dist = _build_point_grids(
        d_estimates_filtered, param_in_use, attack_in_use, degrees_local, logQ_local, verbose)

    # perform fitting
    results = fit_formula(points_est, std_e, std_s, secret_q, params)

    print("\nFormula: ", opts)
    print("Params: ", results, "\n")
