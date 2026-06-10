from latticeestimator.estimator import *
from latticeestimator.estimator.nd import NoiseDistribution
from latticeestimator.estimator.lwe_parameters import LWEParameters
from math import log2


def run_estimator_hybrid():
    nlogs = [14, 15, 16, 17]
    hs = [128, 256, 512]
    logqs = [700,750,800,850,900,950,1000]
    res = []
    for nlog in nlogs:
        for h in hs:
            for logq in logqs:
                FHEParam = LWEParameters(
                    n =2**nlog,
                    q= 2**logq,
                    Xs=ND.SparseTernary(p = h/2, m=h/2, n=2**nlog),
                    Xe=ND.DiscreteGaussian(stddev=3.19)
                )
                try:
                    primal_hybrid_cost = LWE.primal_hybrid(FHEParam, red_cost_model=RC.BDGL16, mitm=False)
                    #print(primal_hybrid_cost)
                    ll = log2(primal_hybrid_cost['rop'])
                    beta = primal_hybrid_cost['beta']
                    eta  = primal_hybrid_cost['eta']
                    zeta = primal_hybrid_cost['zeta']
                    d = primal_hybrid_cost['d']
                    svp_cost = log2(float(primal_hybrid_cost["svp"]))
                    red_cost = log2(float(primal_hybrid_cost["red"]))
                    prob =  primal_hybrid_cost["prob"]
                    #[logq, h, beta, ng, d, wg, lambda]
                    res.append([nlog, h, beta, zeta, d, round(ll)])
                    print(f"n:{nlog} logq:{logq } h:{h} eta: {eta} beta:{beta} lambda:{ll} zeta:{zeta} d:{d} red:{red_cost} svp:{svp_cost} prob:{prob} ")
                    print("------------------------------")
                except Exception as e:
                     print('error in the estimator:', e)
                     print("------------------------------")
    print(res)

def run_estimator():
    nlogs = [14, 15, 16]
    logqs = [700, 750, 800, 850, 900, 950, 1000]
    for nlog in nlogs:
        for logq in logqs:
            FHEParam = LWEParameters(
                n =2**nlog,
                q= 2**logq,
                Xs=ND.UniformMod(2),
                Xe=ND.DiscreteGaussian(stddev=3.19),
            )
            primal_uSVP_cost = LWE.primal_usvp(FHEParam, red_cost_model=RC.BDGL16)
            print(primal_uSVP_cost)
            primal_hybrid_cost_bdd = LWE.bdd(FHEParam, red_cost_model=RC.BDGL16)
            print(primal_hybrid_cost_bdd)


if __name__ == "__main__":
    #run_estimator()
    #run_estimator_hybrid()
    FHEParam = LWEParameters(
                    n =2**16,
                    q= 2**1200,
                    #Xs=ND.SparseTernary(p = 192/2, m=192/2, n=2**15),
                    Xs=ND.Ternary,
                    Xe=ND.DiscreteGaussian(stddev=3.19)
                )
    primal_bdd_cost = LWE.primal_usvp(FHEParam)
    print(primal_bdd_cost)
