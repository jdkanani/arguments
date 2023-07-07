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
ys = xs.copy()
random.shuffle(ys)


def accumulator_factor(i, xs, ys, gamma):
    res = Fp(1)
    for j in range(i + 1):
        res *= (xs[j] + gamma) / (ys[j] + gamma)
    return res


def prove(kzg):
    # proof object
    proof = {}

    # polys represented with n points
    omega = omega_base ** (2**32 // n)
    omega3 = omega_base ** (2**32 // (n * 8))
    ROOTS = [omega**i for i in range(n)]
    ROOTS3 = [omega3**i for i in range(n * 8)]
    PolyEvalRep = polynomialsEvalRep(Fp, omega, n)

    # While multiplying polynomials, to maintain their degree we need to perform product in larger domain
    # since we will be performing three n polynomial multiplication, we need to use 3n domain
    PolyEvalRep3 = polynomialsEvalRep(Fp, omega3, n * 8)

    # get vanishing polynomial
    ZH = vanishing_poly(n)
    ZH_ext = PolyEvalRep3.from_coeffs(ZH)
    ZH_coeffs = ZH_ext.to_coeffs()

    # get gamma (γ) from verifier after commiting to xsp and ysp
    gamma = random_fp_seeded("gamma")
    gamma_poly = PolyEvalRep3.from_coeffs(Poly([gamma]))
    L_1 = PolyEvalRep(ROOTS, [Fp(1)] + [Fp(0) for i in range(len(ROOTS) - 1)])
    L_1_ext = PolyEvalRep3.from_coeffs(L_1.to_coeffs())
    ONE = PolyEvalRep3.from_coeffs(Poly([Fp(1)]))

    # create polynomial from points
    xsp = PolyEvalRep(ROOTS, xs)
    ysp = PolyEvalRep(ROOTS, ys)

    accumulator_poly_eval = [Fp(1)]
    accumulator_poly_eval += [
        accumulator_factor(i, xs, ys, gamma) for i in range(n - 1)
    ]

    # create accumulator polynomial using accumulator evaluations
    accumulator_poly = PolyEvalRep(ROOTS, accumulator_poly_eval)
    # evaluate accumulator polynomial to 3n points
    accumulator_poly = PolyEvalRep3.from_coeffs(accumulator_poly.to_coeffs())
    # shift accumulator polynomial by ω for Z(xω)
    accumulator_poly_shift_evals = eval_poly(accumulator_poly, ROOTS3, ROOTS[1])
    accumulator_poly_shift = PolyEvalRep3(ROOTS3, accumulator_poly_shift_evals)

    # extend accumulator polynomial to 3n points
    xsp_ext = PolyEvalRep3.from_coeffs(xsp.to_coeffs())
    ysp_ext = PolyEvalRep3.from_coeffs(ysp.to_coeffs())

    # commit xsp (f), ysp (g), accumulator_poly (z)
    proof["f"] = {"commitment": kzg.commit(xsp.to_coeffs())}
    proof["g"] = {"commitment": kzg.commit(ysp.to_coeffs())}
    proof["z"] = {"commitment": kzg.commit(accumulator_poly.to_coeffs())}

    # This is where we need larger domain for product as we are multiplying f and Z (accumulator_poly)
    # (f(x) + gamma) * Z(x)
    fz = (xsp_ext + gamma_poly) * accumulator_poly
    # (g(x) + gamma) * Z(xω)
    gz = (ysp_ext + gamma_poly) * accumulator_poly_shift

    # There are two constraints
    # 1. f(x) * Z(x) = g(x) * Z(xω)
    #       equivalent to
    #   f(x) * Z(x) - g(x) * Z(xω) = 0
    # 2. L1(x) * (Z(x) - 1) = 0
    #
    # We will create linear combination of these two constraints

    # create random for linear combination of constraints
    alpha = random_fp_seeded("alpha")

    # tz = α * ((f(x) + gamma) * Z(x) - (g(x) + gamma) * Z(xω)) + α^2 (L1(x) * (Z(x) - 1)))
    tz = ((fz - gz) * alpha) + (L_1_ext * (accumulator_poly - ONE) * alpha**2)

    # Now, we want to check if tz is zero at all points in roots of unity (all of ω/Ω),
    # it must be divisible by vanishing polynomial (this is similar to Round 3 of Plonk)
    # t = tz / ZH
    t = PolyEvalRep3.divideWithCoset(tz.to_coeffs(), ZH_coeffs)
    t = PolyEvalRep3.from_coeffs(t)

    # commit quotient polynomial
    proof["t"] = {"commitment": kzg.commit(t.to_coeffs())}

    # create zeta (ζ) from verifier
    zeta = random_fp_seeded("zeta")

    # # =============================================================================
    # #
    # # Verify locally before creating commitments for verifier everything at zeta (ζ)
    # #

    # # verify fz, gz, L1, tz at zeta
    # assert fz.to_coeffs()(zeta) == (
    #     xsp.to_coeffs()(zeta) + gamma
    # ) * accumulator_poly.to_coeffs()(zeta), "fz != (f(x) + gamma) * Z(x)"

    # accumulator_shift_zeta = eval_poly(accumulator_poly, [zeta * ROOTS[1]])[0]
    # assert gz.to_coeffs()(zeta) == (
    #     ysp.to_coeffs()(zeta) + gamma
    # ) * accumulator_shift_zeta , "gz != (g(x) + gamma) * Z(xω)"

    # assert (fz.to_coeffs()(zeta) - gz.to_coeffs()(zeta)) * alpha + (
    #     L_1_ext.to_coeffs()(zeta)
    #     * (accumulator_poly.to_coeffs()(zeta) - ONE.to_coeffs()(zeta))
    # ) * alpha**2 == tz.to_coeffs()(
    #     zeta
    # ), "tz != α * ((f(x) + gamma) * Z(x) - (g(x) + gamma) * Z(xω)) + α^2 (L1(x) * (Z(x) - 1)))"

    # # verify tz and t at zeta
    # assert tz.to_coeffs()(zeta) == t.to_coeffs()(zeta) * ZH_coeffs(zeta), "tz != t * ZH"
    # #
    # # =============================================================================

    # Create proof for f, g, accumulator_poly, t at zeta
    #
    # This is similar to Round 4 and 5 of Plonk.
    # We are doing unoptimized version of those round and have multiple openings for each polynomial.
    # This can be optimized by using batch opening with linear combination of polynomials using new challenge Nu (ν).

    # Note that kzg.open(poly, a) returns proof and value of polynomial (b) at point a

    # open f, g, t, accumulator_poly and store proof and value at zeta
    proof["f"]["proof"], proof["f"]["zeta_value"] = kzg.open(xsp.to_coeffs(), zeta)
    proof["g"]["proof"], proof["g"]["zeta_value"] = kzg.open(ysp.to_coeffs(), zeta)
    proof["t"]["proof"], proof["t"]["zeta_value"] = kzg.open(t.to_coeffs(), zeta)
    proof["z"]["proof"], proof["z"]["zeta_value"] = kzg.open(
        accumulator_poly.to_coeffs(), zeta
    )
    proof["z"]["shift_proof"], proof["z"]["shift_zeta_value"] = kzg.open(
        accumulator_poly.to_coeffs(), zeta * ROOTS[1]
    )

    # verify if -
    # t * ZH == α * ((f(x) + gamma) * Z(x) - (g(x) + gamma) * Z(xω)) + α^2 (L1(x) * (Z(x) - 1)))
    # at zeta (ζ)
    L1Constrainst = alpha**2 * (
        # L1(x) * (Z(x) - 1)
        L_1.to_coeffs()(zeta)
        * (proof["z"]["zeta_value"] - Fp(1))
    )
    TransitionConstraint = alpha * (
        # (f(x) + gamma) * Z(x)
        (proof["f"]["zeta_value"] + gamma) * proof["z"]["zeta_value"]
        # (g(x) + gamma) * Z(xω)
        - (proof["g"]["zeta_value"] + gamma) * proof["z"]["shift_zeta_value"]
    )
    # t * ZH == α * ((f(x) + gamma) * Z(x) - (g(x) + gamma) * Z(xω)) + α^2 (L1(x) * (Z(x) - 1)))
    assert proof["t"]["zeta_value"] * ZH(zeta) == TransitionConstraint + L1Constrainst


if __name__ == "__main__":
    # setup with degree 100 (all polynomial must be of degree < 100)
    _, CRS = generate_CRS(100)

    # create KZG object
    kzg = KZG(CRS)

    print("xs", xs)
    print("ys", ys)

    # create proof
    proof = prove(kzg)
