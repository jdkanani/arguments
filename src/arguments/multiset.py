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
    eval_poly,
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

    # While multiplying polynomials, to maintain their degree we need to perform product in larger domain
    # since we will be performing two n polynomial multiplication, we need to use 2n domain
    PolyEvalRep2 = polynomialsEvalRep(Fp, omega2, n * 2)

    # get vanishing polynomial
    ZH = vanishing_poly(n)
    ZH_ext = PolyEvalRep2.from_coeffs(ZH)
    ZH_coeffs = ZH_ext.to_coeffs()

    # get gamma from verifier after commiting to xsp and ysp
    gamma = random_fp()
    gamma_poly = PolyEvalRep2.from_coeffs(Poly([gamma]))
    L_1 = PolyEvalRep(ROOTS, [Fp(1)] + [Fp(0) for i in range(len(ROOTS) - 1)])
    L_1_ext = PolyEvalRep2.from_coeffs(L_1.to_coeffs())
    ONE = PolyEvalRep2.from_coeffs(Poly([Fp(1)]))

    # create polynomial from points
    xsp = PolyEvalRep(ROOTS, xs)
    ysp = PolyEvalRep(ROOTS, ys)

    accumulator_poly_eval = [Fp(1)]
    accumulator_poly_eval += [
        accumulator_factor(i, xs, ys, gamma) for i in range(n - 1)
    ]

    # create accumulator polynomial using accumulator evaluations
    accumulator_poly = PolyEvalRep(ROOTS, accumulator_poly_eval)

    # shift accumulator polynomial to 2n points
    accumulator_poly = PolyEvalRep2.from_coeffs(accumulator_poly.to_coeffs())
    accumulator_poly_shift_evals = eval_poly(accumulator_poly, ROOTS2, ROOTS[1])
    accumulator_poly_shift = PolyEvalRep2(ROOTS2, accumulator_poly_shift_evals)

    # extend accumulator polynomial to 2n points
    xsp_ext = PolyEvalRep2.from_coeffs(xsp.to_coeffs())
    ysp_ext = PolyEvalRep2.from_coeffs(ysp.to_coeffs())

    # This is where we need larger domain for product as we are multiplying f and Z (accumulator_poly)
    # (f(x) + gamma) * Z(x)
    fz = (xsp_ext + gamma_poly) * accumulator_poly
    # (g(x) + gamma) * Z(xw)
    qz = (ysp_ext + gamma_poly) * accumulator_poly_shift

    # There are two constraints
    # 1. f(x) * Z(x) = g(x) * Z(xw) 
    #       equivalent to 
    #   f(x) * Z(x) - g(x) * Z(xw) = 0
    # 2. L1(x) * (Z(x) - 1) = 0
    # 
    # We will create linear combination of these two constraints

    # create random for linear combination of constraints
    alpha = random_fp_seeded("alpha")

    # tz = α * ((f(x) + gamma) * Z(x) - (g(x) + gamma) * Z(xw)) + α^2 (L1(x) * (Z(x) - 1)))
    tz = ((fz - qz) * alpha) + (L_1_ext * (accumulator_poly - ONE) * alpha**2)

    # Now, we want to check if tz is zero at all points in roots of unity (all of omegas),
    # it must be divisible by vanishing polynomial
    # q = tz / ZH
    q = PolyEvalRep2.divideWithCoset(tz.to_coeffs(), ZH_coeffs)

    # commit xsp, ysp, accumulator_poly
    # kzg = KZG(CRS)
    # commitment_xsp = kzg.commit(xsp)
    # commitment_ysp = kzg.commit(ysp)
    # commitment_accumulator_poly = kzg.commit(accumulator_poly)


if __name__ == "__main__":
    # setup with degree 100 (all polynomial must be of degree < 100)
    _, CRS = generate_CRS(100)
    prove()
