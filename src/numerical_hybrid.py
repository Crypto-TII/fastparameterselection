import random
from sage.all import RR, binomial, RealDistribution, prod
from numpy import pi, exp, log, log2, sqrt, ceil
from scipy.optimize import fsolve

import warnings


# EXACT EQUATIONS:
# eq1a = lambda ng_, beta_, d_: l - (0.292*beta_+16.4+3+log2(d_)) + probability_enum(n, h, ng_, wg)-1 # log2(binomial(n-h,ng_-wg)*binomial(h,wg)/binomial(n, ng_) + binomial(n-h,ng_-wg+1)*binomial(h,wg-1)/binomial(n, ng_)) #approx_binom(n-h, ng_-wg)+approx_binom(h, wg)-approx_binom(n, ng_) - 2   #log2(binomial(n-h,ng_-wg)*binomial(h,wg)/binomial(n, ng_) + binomial(n-h,ng_-wg+1)*binomial(h,wg-1)/binomial(n, ng_)) #(approx_binom(ng_,wg)+wg)+log2(d_) - (0.292*beta_+16.4+3)+2  #0.5*(ng_*entropy(wg/ng_)+wg) #- 0.5*log2(8*n*wg/ng_*(1-wg/ng_)) #0.5*(log2(math.comb(ng_, wg),2)+wg)
# eq1b = lambda ng_, beta_, d_:  l - ss_enum(ng_,wg)-2*log2(d_) + probability_enum(n, h, ng_, wg)-12 #approx_binom(n-h, ng_-wg)+approx_binom(h, wg)-approx_binom(n, ng_)
# eq2 = lambda ng_, beta_, d_, logq: d_ - ceil(sqrt(n * int(logq) * ln2 / log(_delta(beta_)))) + ng_ - 1 #adding ceil over sqrt() makes the eq. precise
# eq3 = lambda ng_, beta_, d_, logq: (-d_+1)*log2(_delta(beta_))+ ((d_-n+ng_-1)*int(logq)+ (n-ng_)*log2(xi))/d_ - log2(2*sigma_e*sigma_e) #success probability of Babai=1, i.e. ||b_d*|| = sigma_e \approx 4

# [logq, h, beta, ng, d, wg, lambda]
all_params14 = [[250.0, 64, 311, 7315, 17222, 8, 152.10], [300.0, 64, 243, 7105, 17690, 6, 135.44],
                [400.0, 64, 197, 5496, 21313, 5, 113.43], [
                    450.0, 64, 222, 3181, 26328, 6, 106.69],
                [500.0, 64, 160, 4346, 23821, 4, 97.83],  [
                    550.0, 64, 175, 2416, 27915, 5, 91.67],

                [700, 64, 116, 2253, 28240, 4, 76.15782878506808], [
                    750, 64, 106, 1920, 28924, 4, 72.00766336007604],
                [800, 64, 108, 798, 31207, 4, 68.46461135817178], [
    850, 64, 104, 164, 32515, 5, 64.77966727705204],
    [900, 64, 93, 106, 32623, 5, 61.58602469289659], [
    950, 64, 91, 0, 33456, 0, 61.00200074471519],
    [1000, 64, 91, 0, 34325, 0, 61.038996729309524],


    [250.0, 128, 417, 5299, 21819, 11, 190.64], [
    300.0, 128, 386, 3750, 25173, 11, 165.33],
    [350.0, 128, 363, 2236, 28353, 11, 146.87], [
    400.0, 128, 305, 1984, 28852, 9, 130.17],
    [450.0, 128, 262, 1714, 29399, 8, 117.08], [
    500.0, 128, 231, 1339, 30159, 7, 106.55],
    [550.0, 128, 202, 1153, 30526, 6, 97.89],
    [700, 128, 148, 265, 32324, 7, 77.77841686845393], [
    750, 128, 133, 125, 32616, 7, 73.27134882222275],
    [800, 128, 120, 0, 32890, 0, 69.44536145158476], [
    850, 128, 105, 88, 32670, 6, 65.10664522401244],
    [900, 128, 94, 43, 32769, 7, 61.99789614522741], [
    950, 128, 91, 0, 33456, 0, 61.00200074471519],
    [1000, 128, 91, 0, 34325, 0, 61.038996729309524],

    [250.0, 192, 532, 3231, 26322, 16, 208.80], [
    300.0, 192, 469, 1920, 29035, 15, 178.02],
    [350.0, 192, 388, 1578, 29718, 12, 153.90], [
    400.0, 192, 330, 1217, 30443, 11, 135.22],
    [450.0, 192, 280, 1073, 30722, 9, 120.75], [
    500.0, 192, 255, 391, 32118, 10, 108.88],
    [550.0, 192, 219, 406, 32078, 8, 99.33],
    [700, 192, 153, 0, 32900, 0, 79.08179996359026], [
    750, 192, 133, 125, 32616, 7, 73.27134882222275],
    [800, 192, 119, 85, 32732, 7, 69.21360645769066], [
    850, 192, 106, 59, 32776, 7, 65.44942196868422],
    [900, 192, 94, 62, 32750, 6, 61.925730501613074], [
    950, 192, 91, 0, 33456, 0, 61.00200074471519],
    [1000, 192, 91, 0, 34325, 0, 61.038996729309524],

    [250.0, 256, 417, 5299, 21819, 11, 190.64], [
    300.0, 256, 386, 3750, 25173, 11, 165.33],
    [350.0, 256, 363, 2236, 28353, 11, 146.87], [
    400.0, 256, 305, 1984, 28852, 9, 130.17],
    [450.0, 256, 262, 1714, 29399, 8, 117.08], [
    500.0, 256, 231, 1339, 30159, 7, 106.55],
    [550.0, 256, 202, 1153, 30526, 6, 97.89],
    [700, 256, 150, 151, 32563, 8, 78.35999644908343], [
    750, 256, 133, 125, 32616, 7, 73.27134882222275],
    [800, 256, 119, 59, 32758, 8, 69.25413802826124], [
    850, 256, 107, 0, 32913, 0, 65.65037080108321],
    [900, 256, 94, 62, 32750, 6, 61.925730501613074], [
    950, 256, 91, 0, 33456, 0, 61.00200074471519],
    [1000, 256, 91, 0, 34325, 0, 61.038996729309524],

    [700, 512, 152, 58, 32780, 11, 78.95141727596588], [
    750, 512, 134, 63, 32745, 9, 73.64342664307411],
    [800, 512, 119, 59, 32758, 8, 69.25413802826124], [
    850, 512, 106, 59, 32776, 7, 65.44942196868422],
    [900, 512, 95, 0, 32894, 0, 62.145546900238294], [950, 512, 91, 10, 33446, 10, 61.0016100106636], [1000, 512, 91, 10, 34315, 10, 61.038617969222045]]

all_params15 = [[700, 64, 186, 14350, 34934, 5, 126], [750, 64, 228, 10330, 43999, 6, 121],
                [800, 64, 218, 9493, 45832, 6, 116], [
                    850, 64, 208, 8745, 47454, 6, 111],
                [900, 64, 208, 7327, 50501, 6, 107], [
                    950, 64, 182, 8014, 49024, 5, 103],
                [1000, 64, 200, 5162, 55061, 6, 100],

                [700, 128, 366, 4166, 57182, 11, 149], [
                    750, 128, 339, 3677, 58191, 10, 140],
                [800, 128, 284, 5189, 55043, 8, 132], [
    850, 128, 289, 3112, 59331, 9, 125],
    [900, 128, 264, 3175, 59203, 8, 119], [
    950, 128, 255, 2207, 61164, 8, 113],
    [1000, 128, 229, 2734, 60086, 7, 108],

    [700, 192, 390, 2897, 59804, 12, 155], [
    750, 192, 357, 2642, 60322, 11, 145],
    [800, 192, 328, 2416, 60773, 10, 138], [
    850, 192, 311, 1652, 62330, 10, 129],
    [900, 192, 283, 1806, 62005, 9, 122], [
    950, 192, 272, 917, 63796, 10, 116],
    [1000, 192, 253, 803, 64049, 9, 110],

    [700, 256, 403, 2224, 61191, 12, 158], [
    750, 256, 371, 1842, 61955, 11, 148],
    [800, 256, 335, 1984, 61657, 10, 139], [
    850, 256, 322, 924, 63807, 11, 131],
    [900, 256, 297, 812, 64024, 10, 123], [
    950, 256, 277, 563, 64537, 10, 116],
    [1000, 256, 258, 390, 64874, 10, 111],

    [700, 512, 436, 589, 64588, 17, 163], [
    750, 512, 396, 473, 64771, 16, 151],
    [800, 512, 362, 392, 64946, 15, 141], [
    850, 512, 331, 361, 64973, 14, 132],
    [900, 512, 306, 241, 65241, 14, 125], [
    950, 512, 285, 0, 65712, -1, 119],
    [1000, 512, 261, 204, 65304, 12, 112]
]

all_params16 = [[700, 64, 502, 29972, 66900, 12, 204.3086856519731], [750, 64, 496, 27694, 72150, 12, 199.37220382253997],
                [800, 64, 299, 37512, 49128, 7, 174.33098274157615], [
                    850, 64, 292, 36075, 52534, 7, 168.79487539660414],
                [900, 64, 251, 37402, 49384, 6, 163.78559399692077], [
                    950, 64, 303, 31484, 63348, 7, 160.90719271637388],

                [1000, 64, 247, 34312, 56697, 6, 154.64053954098628], [
                    700, 128, 539, 28087, 71276, 13, 249.88276051600195],
                [750, 128, 501, 27433, 72768, 12, 238.6364372746093], [
                    800, 128, 566, 21107, 86998, 14, 231.7410914426918],
                [850, 128, 492, 22833, 83154, 12, 219.74323550817536], [
                    900, 128, 491, 20366, 88615, 12, 211.576871641367],
                [950, 128, 458, 20137, 89112, 11, 203.43093137911046], [
                    1000, 128, 475, 16489, 97048, 12, 196.5746228989299],

                [700, 192, 693, 20518, 88335, 17, 286.296735151414], [
                    750, 192, 644, 19810, 89879, 16, 271.35354500033776],
                [800, 192, 674, 15094, 100094, 17, 258.8319896335122], [
                    850, 192, 598, 16406, 97262, 15, 245.46310719147453],
                [900, 192, 588, 14138, 102114, 15, 233.9091332247279], [
                    950, 192, 575, 12149, 106327, 15, 223.40413201505788],
                [1000, 192, 541, 11737, 107182, 14, 213.62315678892801],

                [700, 256, 840, 13621, 103275, 22, 308.4712511484923], [
                    750, 256, 762, 13788, 102902, 20, 290.042959469996],
                [800, 256, 718, 12699, 105197, 19, 273.9726320045856], [
                    850, 256, 679, 11646, 107409, 18, 259.33443985905603],
                [900, 256, 707, 6747, 117592, 20, 248.70589123795554], [
                    950, 256, 664, 6265, 118573, 19, 235.78921886091717],
                [1000, 256, 623, 5967, 119170, 18, 223.9543039378368],

                [700, 512, 1030, 5050, 121177, 32, 342.8540670705324], [
                    750, 512, 947, 4672, 121897, 29, 318.9133585374426],
                [800, 512, 877, 4260, 122739, 27, 297.85611374146856], [
                    850, 512, 813, 3997, 123255, 25, 279.268050680993],
                [900, 512, 760, 3539, 124163, 24, 262.81633359803743], [
                    950, 512, 707, 3498, 124248, 22, 248.03356251829015],
                [1000, 512, 664, 3162, 124919, 21, 234.78705890453406]
                ]

all_params17 = [
    [2500.0, 64, 218, 58398, 137200, 4, 138.81], [
        2600.0, 64, 219, 55194, 144564, 4, 135.71],
    [2700.0, 64, 221, 51762, 152383, 4, 132.85], [
        2800.0, 64, 211, 51160, 153744, 4, 129.66],
    [2900.0, 64, 252, 38341, 182253, 5, 128.32], [
        3000.0, 64, 216, 44115, 169535, 4, 124.48],
    [3100.0, 64, 216, 41197, 175985, 4, 122.03], [
        3200.0, 64, 219, 37476, 184136, 4, 119.93],
    [3300.0, 64, 202, 39284, 180182, 4, 117.13], [
        3400.0, 64, 217, 32180, 195597, 4, 115.49],

    [800, 128, 699, 81027, 84065, 15, 350.09968802677224], [
        850, 128, 696, 77483, 92421, 15, 339.0132890998411],
    [900, 128, 649, 77111, 93297, 14, 327.8948900547392], [
        950, 128, 609, 76624, 94446, 13, 318.0851735688354],
    [1000, 128, 688, 67530, 115975, 15, 311.64359617713654],
    [2500, 128, 392, 23673, 213728, 8, 167.18], [
        2600, 128, 378, 22134, 216956, 8, 161.46],
    [2700, 128, 393, 14848, 232083, 9, 158.49], [
        2800, 128, 346, 20623, 220113, 7, 152.01],
    [2900, 128, 359, 13770, 234320, 8, 148.17], [
        3000, 128, 347, 12482, 236947, 8, 144.15],
    [3100, 128, 303, 19210, 223051, 6, 139.6], [
        3200, 128, 319, 11549, 238853, 7, 135.73],

    [700, 192, 882, 78641, 89721, 19, 456.5103580338542], [
        750, 192, 1036, 66339, 118846, 23, 441.2132564899414],
    [800, 192, 830, 73100, 102840, 18, 421.8560209553864], [
        850, 192, 1027, 57837, 138655, 23, 412.354268400881],
    [900, 192, 862, 63281, 125989, 19, 392.2665644085177], [
        950, 192, 1009, 50151, 156190, 23, 387.4267214325475],
    [1000, 192, 854, 56030, 142783, 19, 367.54015963955896],

    [700, 256, 1203, 63227, 126183, 27, 518.5415474215844],
    [750, 256, 1040, 66150, 119306, 23, 496.875156416131], [
        800, 256, 1072, 59937, 133821, 24, 475.4616704949971],
    [850, 256, 1064, 55805, 143348, 24, 456.56651201034595], [
        900, 256, 1016, 54080, 147284, 23, 439.07310472202175],
    [950, 256, 1048, 47774, 161558, 24, 422.94916480509175], [
        1000, 256, 1034, 44243, 169433, 24, 407.94532844240103],

    [700, 512, 1848, 35889, 187954, 46, 662.928570615052], [
        750, 512, 1756, 33078, 194012, 44, 625.5100119944597],
    [800, 512, 1667, 30700, 199099, 42, 591.5489890228017], [
        850, 512, 1831, 16256, 229510, 51, 579.4601853420198],
    [900, 512, 1504, 26837, 207301, 38, 532.2872140012998], [
        950, 512, 1498, 21375, 218804, 39, 506.1205897513737],
    [1000, 512, 1394, 21859, 217769, 36, 481.9016547943525]]


const = 2 * pi * exp(1)
ln2 = log(2)
e = exp(1)

epsilon = 1e-10
# def safe_log2(x):
#     return log2(max(epsilon, x))


def _delta(beta):
    small = (
        (2, 1.02190),
        (5, 1.01862),
        (10, 1.01616),
        (15, 1.01485),
        (20, 1.01420),
        (25, 1.01342),
        (28, 1.01331),
        (40, 1.01295),
    )
    return (beta / (2 * pi * e) * (pi * beta) ** (1 / beta)) ** (1 / (2 * (beta - 1)))


def log2_delta(beta):
    """

    :param beta: bkz block size
    :return: log2(_delta(beta))
    """
    return ((1 / (2 * (beta - 1))) * log2(beta / (2 * pi * e)) + (1 / (2 * (beta - 1))) * (1 / beta) * log2(pi * beta))


def entropy(x):
    """
    :param x:
    :return: Binary entropy function of x
    """
    if x == 0 or x==1: return 0
    return -x*log2(x) - (1-x)*log2(1-x)


def approx_binom(n, k):
    """

    :param n: integer
    :param k: integer
    :return: apprixmation to binom(n,k) from Wiki
    """
    if n<=0 or k<=0 or n<=k: return 0
    return n*entropy(k/n)+0.5*log2(n/(8*k*(n-k)))


def approx_startpoint_bdd(n, logq, h):
    """
    This function is used to compute the starting point for numerical_lambda_hybrid_v2()
    :param n: LWE dimension
    :param logq: LWE modulus
    :param h: Ternary secret weight
    :return: tuple (beta, d)
    """
    sigma_s = sqrt(h/n)
    sigma_e = 3.19

    lnq = logq * ln2
    zeta = round(sigma_e/sigma_s)

    # very crude approximation for beta
    beta_initial_guess = n / 4

    def nom(beta): return 2 * n * lnq * log(beta/const)

    def denom(beta): return log(beta/const) + 2 * lnq - 2 * log(sigma_e) - \
        log(const) - 2 * (lnq - log(zeta)) * \
        sqrt(n * log(beta/const) / (2 * lnq * beta))

    def eq6(beta): return beta - nom(beta) / (denom(beta)**2)

    beta_solution = fsolve(eq6, beta_initial_guess, full_output=False)

    d_optimal = sqrt(
        2 * n * lnq * beta_solution[0] / log(beta_solution[0] / const))
    return [beta_solution[0], d_optimal]


def probability_enum(n, h, ng, w):
    """
    This function used to compute exact complexity of hybrid
    :param n: LWE dimension
    :param h: Ternary secret weight
    :param ng: Number of guessed coordinates
    :param w: Weight of guessed secret
    :return: Probability that the secret has weight w on ng coordinates
    """
    prob = 0
    ng = int(ng)
    if ng<0: return  -1
    for i in range(0, w+1):
        prob += RR(binomial(n-h, ng-i)*binomial(h, i)/binomial(n, ng))
    if prob == 0:
        return -1
    return log2(prob)


def ss_enum(ng, w):
    """
    This function used to compute exact complexity of hybrid
    :param ng: Number of guessed coordinates
    :param w: Weight of guessed secret
    :return: Number of vectros of dim ng of weight up to (including) w
    """
    ss = 0
    for i in range(0, w+1):
        ss += RR(binomial(ng, i) * 2**i)
    return log2(ss)


def babai_prob(n, logq, sigma_e, h, beta, d):
    """
    The function is taken from the LatticeEstimator (GSA + Babai probability)
    :param n: LWE dimension
    :param logq: LWE modulus
    :param sigma_e: LWE error st. dev
    :param h:  Ternary secret weight
    :param beta: BKZ parameter
    :param d: optimal lattice dimension
    :return: probability of Babai algorithm
    """
    sigma_s = sqrt(h/n)
    xi = sigma_e / sigma_s
    log_vol = RR(logq * (d - n - 1) + log2(xi) * n)
    r_log = [(d - 1 - 2 * i) * RR(log2_delta(beta)) + log_vol / d for i in range(d)]
    r = [2 ** (2 * r_) for r_ in r_log]

    norm = sqrt(d) * sigma_e
    denom = float(2 * norm) ** 2
    T = RealDistribution("beta", ((len(r) - 1) / 2, 1.0 / 2))
    probs = [1 - T.cum_distribution_function(1 - r_ / denom) for r_ in r]
    P = prod(probs)
    if P <= 0: return -500 #float('-inf')
    else: return log2(P)


# d = m+n-ng
def exact_runtime(n, logq, sigma_e, h, beta, d, ng, w):
    """
    The function is needed to filter solutions
    :param n: LWE dimension
    :param logq: LWE modulus
    :param sigma_e: LWE error st. dev
    :param h: Ternary secret weight
    :param beta: BKZ parameter
    :param d: optimal lattice dimension
    :param ng: Number of guessed coordinates
    :param w: Weight of guessed secret
    :return: tuple (runtime of bkz, runtime of enumeration, babai probability) on log2-scale
    """
    ng = int(ng)
    prob = probability_enum(n, h, ng, w)
    t_bkz  = 0.292*beta + log2(8*d) + 16.4 - prob + 2
    t_enum = ss_enum(ng,w)+2*log2(d) - prob + 1
    babai = babai_prob(n-ng, logq, sigma_e, h, beta, d)
    return t_bkz, t_enum, babai


def numerical_lambda_hybrid_v2(n, logq, sigma_e, h):
    """
    Caller's function for sec. level of hybrid attack
    :param n: LWE dimension
    :param logq: LWE modulus
    :param sigma_e: LWE error st. dev
    :param h: Ternary secret weight
    :return: security level for hybrid attack. Return float('inf') if numerical solver fails
    """

    sigma_s = sqrt(h / n)
    xi = sigma_e / sigma_s

    rt_min = float('inf')
    sol_tolerance = 100 #max solution tolerance

    #Brute-forcing for w
    for wg in range(2, 52):
        def eq1(ng_, beta_, d_): return (approx_binom(ng_, wg) + wg) + \
                                        log2(d_) - (0.292 * beta_ + 16.4 + 3) + 2

        def eq2(ng_, beta_, d_): return d_ - \
                                        ceil(sqrt(n * logq  / log2_delta(beta_))) + ng_ - 1

        def eq3(ng_, beta_, d_): return (-d_ + 1) * log2(_delta(beta_)) + ((d_ - n +
                                                                            ng_ - 1) * logq + (n - ng_) * log2(xi)) / d_ - 2*log2(2*sigma_e)

        def system(x):
            f1 = eq1(x[0], x[1], x[2])
            f2 = eq2(x[0], x[1], x[2])
            f3 = eq3(x[0], x[1], x[2])
            return f1, f2, f3

        # Finding initial point using bdd
        initial_bdd = approx_startpoint_bdd(n, logq, h)
        beta_start = initial_bdd[0]
        d_start = sqrt(2 * n * logq * ln2 * beta_start / log(beta_start / const))
        initial_guess = [n / 4, beta_start, d_start - n / 4]

        bound_trials = 250
        while bound_trials > 0:
            bound_trials -= 1
            # print("Initial guess test:", initial_guess)
            try:
                with warnings.catch_warnings(record=True) as w:
                    warnings.simplefilter("always")
                    res = fsolve(system, initial_guess,
                                 maxfev=2**21, full_output=False)
                    # if len(w) > 0 and issubclass(w[-1].category, RuntimeWarning):
                    # print(f"Warning in fsolve: {w[-1].message}")
                # print("Result from fsolve:", res)
            except Exception as e:
                print(f"Error in fsolve: {e}")
                continue

            if len(res) < 3:
                print("Error: Result from fsolve has fewer than 3 elements")
                continue

            # TODO: compute sol_tolerance only if a solution is found
            sol_tolerance = abs(eq1(res[0], res[1], res[2])) + abs(
                eq2(res[0], res[1], res[2])) + abs(eq3(res[0], res[1], res[2]))

            if sol_tolerance < 2:
                break

            beta_start = max(40, initial_guess[1] + random.randint(-1500, 30))
            d_start = sqrt(2 * n * logq * ln2 * beta_start / log(beta_start / const))
            initial_guess = [n / 4, beta_start, d_start - n / 4]

        prob = probability_enum(n, h, res[0], wg)

        if len(res) >= 3 and prob != -1:

            rt = 0.292 * res[1] + log2(8 * res[2]) + \
                 16.4 - probability_enum(n, h, res[0], wg)
            beta = int(res[1])
            d = int(res[2])
            ng = ceil(sqrt(n * logq / log2_delta(beta))) - d + 1
            bkz_exact, enum_exact, babai = exact_runtime(n, logq, sigma_e, h, beta, d, ng, wg)
            rt_exact = bkz_exact - babai
            #print(wg, rt,  rt_exact, sol_tolerance)
            if rt_exact < rt_min and sol_tolerance < 7:
                #print("min!:", res, rt, wg, sol_tolerance )
                rt_min = rt

    return rt_min


def numerical_logq_starting_point(n, l, sigma_e, h):
    """
    Starting point for logq computations adapted from numerical_logq_bdd()
    :param n: LWE dimension
    :param l: Target sec. level
    :param sigma_e: LWE error st. dev
    :param h: Ternary secret weight
    :return: tuple (logq, beta).
    """
    eta_initial_guess = (l - 16.4) / 0.292
    def eta_eq(eta): return l - (0.292 * eta + log2(8 * eta) + 16.4)
    eta_solution = fsolve(eta_eq, eta_initial_guess, full_output=False)
    eta = eta_solution[0]

    std_s = sqrt(h / n)

    def d_optimal(logq, beta):
        if logq<=0 or beta<=0: return 0
        return sqrt(n * logq/ log2_delta(beta))

    def d(logq, beta): return max(d_optimal(logq, beta), 2*n)
    def eq8(logq, beta): return l - (0.292 * beta +
                                    log2(8 * d(logq, beta)) + 16.4)

    log_std_e = log2(sigma_e)
    zeta = max(0, log_std_e - log2(std_s))

    def eq_for_lnq_and_beta(logq, beta): return eta - d(logq, beta) + 1 / log2_delta(beta) * (logq - log_std_e - 0.5 * log2(const) - n / d(logq, beta) * (logq - zeta))

    def system_bdd_eta_n(x):  # x[0] = n, x[1] = beta
        f1 = eq_for_lnq_and_beta(x[0], x[1])
        f2 = eq8(x[0], x[1])
        return f1, f2

    logq_initial_guess = 121
    solutions_lnq_and_beta = fsolve(
        system_bdd_eta_n, [logq_initial_guess, eta_initial_guess], full_output=True)
    logq_solution = solutions_lnq_and_beta[0][0]
    beta_solustion = solutions_lnq_and_beta[0][1]
    if not solutions_lnq_and_beta[2] == 1:
        return 700, 80
    return logq_solution, beta_solustion

def numerical_logq_hybrid_runoptimize(n, l, sigma_e, h, initial_guess):
    """
    Internal function. Brute-forces over wg-weight
    :param n: LWE dimension
    :param l: Target sec. level
    :param sigma_e: LWE error st.dev
    :param h:  Ternary secret weight
    :param initial_guess: tuple (logq, beta)
    :return: list of tuples [l, n, h, logq, beta, d, ng, wg]
    """
    sigma_s = sqrt(h/n)
    xi = sigma_e / sigma_s

    beta_initial_guess = initial_guess[1] #(l - log2(8*2*n) - 16.4) / 0.292
    logq_initial_guess = initial_guess[0]
    d_initial_guess = ceil(
        sqrt(n * logq_initial_guess / log2_delta(beta_initial_guess)))+ 1 - n/16

    sols = []

    for wg in range(2, 34):
        def eq1a(ng_, beta_, d_): return l - (0.292*beta_+16.4+3+log2(d_)) + approx_binom(n-h,
                                                                                          ng_-wg)+approx_binom(h, wg) - approx_binom(n, ng_) - 2

        def eq1b(ng_, d_): return l - approx_binom(ng_, wg) - wg - 2*log2(d_) + approx_binom(
            n-h, ng_-wg)+approx_binom(h, wg) - approx_binom(n, ng_) - 1

        def eq2(ng_, beta_, d_, logq): return d_ - \
            ceil(sqrt(n * logq / log2_delta(beta_))) + ng_ - 1

        def eq3(ng_, beta_, d_, logq): return (-d_+1)*log2_delta(beta_) + ((d_-n+ng_-1)
                                                                          * logq + (n-ng_)*log2(xi))/d_  - 2*log2(sigma_e) #- log2(2*sigma_e) - 0.5*log2(d_)

        def system(x):
            f1a = eq1a(x[0], x[1], x[2])
            f1b = eq1b(x[0],  x[2])
            f2 = eq2(x[0], x[1], x[2], x[3])
            f3 = eq3(x[0], x[1], x[2], x[3])
            return f1a, f1b, f2, f3

        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                res = fsolve(system, [n/16, beta_initial_guess, d_initial_guess,
                             logq_initial_guess], maxfev=2**21, full_output=False)
                # if len(w) > 0 and issubclass(w[-1].category, RuntimeWarning):
                # print(f"Warning in fsolve: {w[-1].message}")
            # print("Result from fsolve:", res)
        except Exception as e:
            print(f"Error in fsolve: {e}")
            continue

        sol_tolerance = abs(eq1a(res[0], res[1], res[2]))+abs(eq1b(res[0], res[2]))+abs(
            eq2(res[0], res[1], res[2], res[3]))+abs(eq3(res[0], res[1], res[2], res[3]))

        if sol_tolerance > 45: continue


        beta = int(res[1])
        d = int(res[2])
        logq = int(res[3])
        ng = ceil(sqrt(n * logq / log2_delta(beta))) - d + 1
        if ng<=0: continue
        bkz_exact, enum_exact, babai = exact_runtime(n, logq, sigma_e, h, beta, d, ng, wg)
        rt = max(bkz_exact, enum_exact) - babai

        eq1a_tolerance_exact = abs(eq1a(ng, beta, d))

        sol_tolerance_exact = abs(eq1a(ng, beta, d))+abs(eq1b(ng, d))+abs(
            eq2(ng, beta, d, logq))+abs(eq3(ng, beta, d, logq))

        #print(wg, rt, logq, " : ", bkz_exact, babai, ";", sol_tolerance, sol_tolerance_exact, eq1a_tolerance_exact, abs(eq1b(ng, d)),
        #     abs(eq3(ng, beta, d, logq)))

        if (abs(sol_tolerance) < 35 and eq1a_tolerance_exact < 7 and abs(l-rt) < 15):
            sols.append([l, n, h, logq, beta, d, ng, wg])
            #print(l, n, h, logq, beta, d, ng, wg )
    return sols


#
# caller's function: runs numerical_logq_hybrid_runoptimize with different initial guesses
# checks obtained candidates for logq by calling check_candidates_logq()
#
def numerical_logq_hybrid(n, l, h):
    """
    :param n: LWE dimension
    :param l: Target sec. level
    :param h: Ternary secret weight
    :return: logq best scored by check_candidates_logq_exact()
    """

    #TODO: make 3.19 as input parameter
    initial_guess = numerical_logq_starting_point(n, l, 3.19, h)
    res = numerical_logq_hybrid_runoptimize(n, l, 3.19, h, initial_guess)

    # experimentally verified bounds to shift initial points
    if n <= 2**11: B = 100
    elif n <= 2**13: B = 200
    else: B = 300

    # if no candidates found, try again with new initial guess up to bound_trials  times
    bound_trials = 200
    while len(res) == 0 and bound_trials > 0:
        bound_trials -= 1
        shift = random.randint(-B, B)
        init_guess2 = [max(12, initial_guess[0]+shift), max(40, initial_guess[1]-shift)]
        res = numerical_logq_hybrid_runoptimize(n, l, 3.19, h, init_guess2)
        print("re-start numerical_logq_hybrid_runoptimize #", 200-bound_trials, "...")
        if len(res) > 0:
            break

    if len(res) == 0:
        print("numerical_logq_hybrid_runoptimize couldn't find a solution after",
              200-bound_trials, "trial for this parameters set")
        return 450 #float('inf')


    best_logq = check_candidates_logq_exact(res)

    return best_logq

def check_candidates_logq_exact(candidate_list):
    """
    receives a list of logqs and their corresponding [ng, beta, d] found by numerical_logq_hybrid_runoptimize
    outputs logq that gives l_ closes to the target security level l according to exact_runtime()
    :param candidate_list:
    :return: smallest(?) logq scored by exact_runtime()
    """
    best_diff = 70
    best_logq = float("inf")
    found = False

    #sols.append([l, n, h, logq, beta, d, ng, wg])
    for el in candidate_list:
        l = el[0] #target sec. level
        n = el[1]
        h = el[2]
        logq = el[3]
        beta = el[4]
        d = el[5]
        ng = el[6]
        w = el[7]
        tbz, tenum, babai =  exact_runtime(n, logq, 3.19, h, beta, d, ng, w)
        rt = max(tbz,tenum) - babai
        if abs( rt - l) < best_diff+5 and logq<best_logq: # for some reason, it's better to output the smallest found logq. I don't know why
            best_diff = abs( max(tbz,tenum) - babai  - l)
            best_logq = logq
            #print(logq, l, max(tbz,tenum)-babai, best_diff  )
            found = True
            print('found logq:', logq)

    if not found:
        print('No solution found setting best_diff to', best_diff, 'returning fixed value')
        return candidate_list[0][3]

    return round(best_logq)





