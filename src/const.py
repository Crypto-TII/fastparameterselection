import numpy as np

SECRET_DIST = "secret dist."
LAMBDA = "lambda"
LOG_Q = "log q"
USVP = "usvp"
LWE_USVP = "est usvp"
LWE_USVP_C = "* est usvp"
LWE_USVP_F = "est usvp formula"
LWE_BDD_F = "est bdd formula"
USVP_S = "usvp_s"
LWE_USVP_S = "est usvp_s"
USVP_NUM = "usvp num"
LWE_NUM = "est num"
BDD = "bdd"
LWE_BDD = "est bdd"
LWE_BDD_C = "* est bdd"
BDD_S = "bdd_s"
LWE_BDD_S = "est bdd_s"
BDD_NUM = "bdd num"
OUTPUT = "output"
POW = "pow"
LWE_DIM = "lwe dim."
LOGQ_USVP = "logq usvp"
LOGQ_USVP_F = "logq usvp formula"
LOGQ_BDD_F = "logq bdd formula"
LOGQ_BDD = "logq bdd"
HW = "hw"
LOGQ_HYBRID = "logq hybrid"
HYBRID = "hybrid"
LWE_HYBRID = "est hybrid"
STD_E_USVP = "log2(std_e) usvp"
STD_E_BDD = "log2(std_e) bdd"
STD_E_USVP_C = "* log2(std_e) usvp"
STD_E_BDD_C = "* log2(std_e) bdd"
NUM_USVP = "num usvp"
NUM_BDD = "num bdd"
EST = "est"
NUM_CALLS_BDD = "est calls bdd"
NUM_CALLS_USVP = "est calls usvp"

# Lambda parameters for USVP and BDD models
LAMBDA_USVP_BIN = [0.317747, 2.071129, 1.849214]
LAMBDA_USVP_TER = [0.296208, 0.800603, 12.09086]
LAMBDA_USVP_S_BIN = [0.445309, 1.486982, 0.950115, 11.21416]
LAMBDA_USVP_S_TER = [0.833542, 0.154947, 1.469823, 18.09877]

LAMBDA_BDD_BIN = [0.26497, 3.25511, -13.69437]
LAMBDA_BDD_TER = [0.28891, 0.87868, 19.1069]
LAMBDA_BDD_S_BIN = [0.424578, 2.122152, 1.959558, 1.155390]
LAMBDA_BDD_S_TER = [0.606897, 0.476667, 0.667667, 15.20932]

# LWE dimension parameters for USVP and BDD models
N_USVP_BIN = [1.02575, 0.17241, 34.84910]
N_USVP_TER = [1.05153, 0.52652, 43.20997]
N_USVP_S_BIN = [-1.142080, 0.231197, 1.106616, -0.233138]
N_USVP_S_TER = [-1.073049, 0.278319, 0.931202, 0.792882]

N_BDD_BIN = [1.154587, -46.18551, -4.457340, 0.809972]
N_BDD_TER = [1.417954, -48.44275, -2.871196, 1.884925]
N_BDD_S_BIN = [0.463730, -1.634159, 5.236220, 1.818256]
N_BDD_S_TER = [2.755987, -10.41781, 0.869780, 0.318689]

# Experimentally found values to offer a correction to the output of logq numerical


# def init_correction_values():
#     x_vals = np.array([80, 100, 128, 192, 256])
#     y_vals = np.array([2**10, 2**11, 2**15])
#     z = np.array([4, 3, 2, 1, 0.5])
#     t = np.array([190, 100, 80, 30, 20])

#     # Build the input grid and function values
#     points = []
#     values = []
#     for y in y_vals:
#         for xi, zi, ti in zip(x_vals, z, t):
#             val = zi + ti if y == 2**15 else zi
#             points.append([y, xi])
#             values.append(val)

#     points = np.array(points)
#     values = np.array(values)
#     return points, values


# POINTS, VALUES = init_correction_values()
