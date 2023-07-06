from arguments.setup import (
    G,
    G2,
    GT,
    generate_CRS,
    random_fp,
    random_fp_seeded,
    vanishing_poly,
    omega_base,
    polynomialsEvalRep,
    eval_poly
)
from utils.ssbls12 import Poly, Fp
from arguments.pcs.kzg import KZG
import random

# closest power of 2 for length of list of points
n = 8

# In multiset argument, prover will try to convience verifier that two list of points are equal
# To do that, we will be using grand product check, xs1 * xs2 * ... = ys1 * ys2 * ...
# With random from verifier, (xs1 + γ) * (xs2 + γ) * ... = (ys1 + γ) * (ys2 + γ) * ...
#
# 1 = (xs1 + γ) * (xs2 + γ) * ... / (ys1 + γ) * (ys2 + γ) * ...
#

xs = [Fp(i) * Fp(i) for i in range(1, n + 1)]
print(xs)
ys = xs.copy()
random.shuffle(ys)
print(ys)


def accumulator_factor(i, xs, ys, gamma):
    res = Fp(1)
    for j in range(i + 1):
        res *= (xs[j] + gamma) / (ys[j] + gamma)
    return res


def prove():
    # polys represented with n points
    omega = omega_base ** (2**32 // n)
    omega2 = omega_base ** (2**32 // (n * 2))
    ROOTS = [omega**i for i in range(n)]
    ROOTS2 = [omega**i for i in range(n * 2)]
    PolyEvalRep = polynomialsEvalRep(Fp, omega, n)
    PolyEvalRep2 = polynomialsEvalRep(Fp, omega2, n * 2)

    # get vanishing polynomial
    ZH = vanishing_poly(n)
    ZH_ext = PolyEvalRep2.from_coeffs(ZH)
    ZH_coeffs = ZH_ext.to_coeffs()

    # get gamma from verifier after commiting to xsp and ysp
    gamma = random_fp()
    gamma_poly = PolyEvalRep.from_coeffs(Poly([gamma]))
    L_1 = PolyEvalRep(ROOTS, [Fp(1)] + [Fp(0) for i in range(len(ROOTS) - 1)])

    # create polynomial from points
    xsp = PolyEvalRep(ROOTS, xs)
    ysp = PolyEvalRep(ROOTS, ys)

    accumulator_poly_eval = [Fp(1)]
    accumulator_poly_eval += [
        accumulator_factor(i, xs, ys, gamma) for i in range(n - 1)
    ]

    # create accumulator polynomial using accumulator evaluations
    accumulator_poly = PolyEvalRep(ROOTS, accumulator_poly_eval)
    accumulator_poly_shift_evals = eval_poly(accumulator_poly, ROOTS, ROOTS[1])
    accumulator_poly_shift = PolyEvalRep(ROOTS, accumulator_poly_shift_evals)

    # (f(x) + gamma) * Z(x)
    fz = (xsp + gamma_poly) * accumulator_poly
    # (g(x) + gamma) * Z(xw)
    qz = (ysp + gamma_poly) * accumulator_poly_shift

    print("fz ==>", fz.to_coeffs().coefficients)
    print("qz ==>", qz.to_coeffs().coefficients)

    # quotient polynomial
    # ((f(x) + gamma) * Z(x) - (g(x) + gamma) * Z(xw)) / ZH(x)
    # tz = (fz - qz)
    # print("tz.to_coeffs() ==>", tz.to_coeffs().coefficients)


    # print(accumulator_poly.to_coeffs().coefficients)
    # print(accumulator_poly_shift.to_coeffs().coefficients)
    # print("ZH_coeffs", ZH_coeffs)
    # print("fz ==>", fz.to_coeffs().coefficients)
    # print("qz ==>", qz.to_coeffs().coefficients)
    # print("fz ==>", PolyEvalRep.divideWithCoset(fz.to_coeffs(), ZH_coeffs).coefficients)
    # print("ZH_coeffs ==>", ZH_coeffs)
    # tz = PolyEvalRep.divideWithCoset(tz.to_coeffs(), ZH_coeffs)

    # commit xsp, ysp, accumulator_poly
    # kzg = KZG(CRS)
    # commitment_xsp = kzg.commit(xsp)
    # commitment_ysp = kzg.commit(ysp)
    # commitment_accumulator_poly = kzg.commit(accumulator_poly)


if __name__ == "__main__":
    # setup with degree 100 (all polynomial must be of degree < 100)
    _, CRS = generate_CRS(100)
    prove()
