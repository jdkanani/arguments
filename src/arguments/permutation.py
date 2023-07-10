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

# In permutation argument, prover will try to convience verifier that two list of points are equal
# with given permutation function (σ).
# σ : [n] → [n]
# g = σ(f)


# f = [33, 55, 11, 77, 88, 66, 44, 22]
# g = [11, 22, 33, 44, 55, 66, 77, 88]
#
# i = [0, 1, 2, 3, 4, 5, 6, 7]
# σ = lambda i: [2, 7, 0, 6, 1, 5, 3, 4][i]
#
# g = [f[σ(i)] for i in range(n)]
#
#
# To do that, we will be using grand product check with index,
# (xs1 + i1) * (xs2 + i2) * ... = (ys1 + σ(i1)) * (ys2 + σ(i2)) * ...
#
# With random from verifier,
# (xs1 + β * i1 + γ) * (xs2 + β * i2 + γ) * ... = (ys1 + β * σ(i1) + γ) * (ys2 + β * σ(i2) + γ) * ...
#
# 1 = (xs1 + β * i1 + γ) * (xs2 + β * i2 + γ) * ... / * (ys1 + β * σ(i1) + γ) * (ys2 + β * σ(i2) + γ) * ...
#


xs = [Fp(i) * Fp(i) for i in range(1, n + 1)]
ys = xs.copy()
random.shuffle(ys)

# σ(i) = j
# σ = lambda i: sigma[i]
sigma = [xs.index(y) for y in ys]
#
# In production and plonk, σ-domain and σ-id are calculated using rotating given wires with same values
# def permute_indices(wires):
#     size = len(wires)
#     permutation = [i for i in range(size)]
#     for i in range(size):
#         for j in range(i, size):
#             if wires[i] == wires[j]:
#                 permutation[i], permutation[j] = permutation[j], permutation[i]
#                 break
#     return permutation
# permutation = permute_indices(wires)
# Sid = range(0)
# Sdomain = [Sid[i] for i in permutation]
#


def accumulator_factor(i, xs, ys, beta, id_domain, perm_domain, gamma):
    res = Fp(1)
    for j in range(i + 1):
        # res *= (xs1 + β * i1 + γ) / (ys1 + β * σ(i1) + γ)
        numerator = xs[j] + beta * id_domain[j] + gamma
        denominator = ys[j] + beta * perm_domain[j] + gamma
        res *= numerator / denominator
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

    L_1 = PolyEvalRep(ROOTS, [Fp(1)] + [Fp(0) for i in range(len(ROOTS) - 1)])
    L_1_ext = PolyEvalRep3.from_coeffs(L_1.to_coeffs())
    ONE = PolyEvalRep3.from_coeffs(Poly([Fp(1)]))

    # create id domain and permutation domains
    # id domain = [0, 1, 2, ..., n-1]
    id_domain = [Fp(i) for i in range(n)]
    # permutation domain = [σ(0), σ(1), σ(2), ..., σ(n-1)]
    permutation_domain = [Fp(sigma[i]) for i in range(n)]

    # create Sid and Ssigma polynomials
    Sid = PolyEvalRep(ROOTS, id_domain)
    Ssigma = PolyEvalRep(ROOTS, permutation_domain)

    # get commitment Sid and Ssigma polynomials
    # ideally, this is preprocessed and stored with prover and verifier
    proof["Sid"] = {"commitment": kzg.commit(Sid.to_coeffs())}
    proof["Ssigma"] = {"commitment": kzg.commit(Ssigma.to_coeffs())}

    # create polynomial from points
    xsp = PolyEvalRep(ROOTS, xs)
    ysp = PolyEvalRep(ROOTS, ys)

    # commit xsp and ysp as f and g
    proof["f"] = {"commitment": kzg.commit(xsp.to_coeffs())}
    proof["g"] = {"commitment": kzg.commit(ysp.to_coeffs())}

    # get beta (β) and gamma (γ) from verifier after commiting to xsp and ysp
    transcript = [proof["f"]["commitment"], proof["g"]["commitment"]]
    beta = random_fp_seeded(str(transcript) + "0")
    gamma = random_fp_seeded(str(transcript) + "1")
    # beta_poly = β
    beta_poly = PolyEvalRep3.from_coeffs(Poly([beta]))
    # gamma_poly = γ
    gamma_poly = PolyEvalRep3.from_coeffs(Poly([gamma]))

    accumulator_poly_eval = [Fp(1)]
    accumulator_poly_eval += [
        accumulator_factor(i, xs, ys, beta, id_domain, permutation_domain, gamma)
        for i in range(n - 1)
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

    # extend Sid and Ssigma to 3n points
    Sid_ext = PolyEvalRep3.from_coeffs(Sid.to_coeffs())
    Ssigma_ext = PolyEvalRep3.from_coeffs(Ssigma.to_coeffs())

    # commit accumulator_poly (z)
    proof["z"] = {"commitment": kzg.commit(accumulator_poly.to_coeffs())}

    # This is where we need larger domain for product as we are multiplying f and Z (accumulator_poly)
    # (f(x) + β * Sid(x) + γ) * Z(x)
    fz = (xsp_ext + beta_poly * Sid_ext + gamma_poly) * accumulator_poly
    # (g(x) + β * Sσ(x) + γ) * Z(xω)
    gz = (ysp_ext + beta_poly * Ssigma_ext + gamma_poly) * accumulator_poly_shift

    # There are two constraints
    # 1. f(x) * Z(x) = g(x) * Z(xω)
    #       equivalent to
    #   f(x) * Z(x) - g(x) * Z(xω) = 0
    # 2. L1(x) * (Z(x) - 1) = 0
    #
    # We will create linear combination of these two constraints

    # get random α from verifier using transcript
    transcript += [proof["z"]["commitment"]]
    alpha = random_fp_seeded(str(transcript))

    # create random for linear combination of constraints
    # tz = α * ((f(x) + β * Sid(x) + γ) * Z(x) - (g(x) + β * Sσ(x) + γ) * Z(xω)) + α^2 (L1(x) * (Z(x) - 1)))
    tz = ((fz - gz) * alpha) + (L_1_ext * (accumulator_poly - ONE) * alpha**2)

    # Now, we want to check if tz is zero at all points in roots of unity (all of ω/Ω),
    # it must be divisible by vanishing polynomial (this is similar to Round 3 of Plonk)
    # t = tz / ZH
    t = PolyEvalRep3.divideWithCoset(tz.to_coeffs(), ZH_coeffs)
    t = PolyEvalRep3.from_coeffs(t)

    # commit quotient polynomial
    proof["t"] = {"commitment": kzg.commit(t.to_coeffs())}

    # create zeta (ζ) from verifier using transcript
    transcript += [proof["t"]["commitment"]]
    zeta = random_fp_seeded(str(transcript))

    # # =============================================================================
    # #
    # # Verify locally before creating commitments for verifier everything at zeta (ζ)
    # #

    # # verify fz, gz, L1, tz at zeta
    # assert fz.to_coeffs()(zeta) == (
    #     xsp.to_coeffs()(zeta) + beta * Sid.to_coeffs()(zeta) + gamma
    # ) * accumulator_poly.to_coeffs()(zeta), "fz != (f(x) + β * Sid(x) + γ) * Z(x)"

    # accumulator_shift_zeta = eval_poly(accumulator_poly, [zeta * ROOTS[1]])[0]
    # assert (
    #     gz.to_coeffs()(zeta)
    #     == (ysp.to_coeffs()(zeta) + beta * Ssigma.to_coeffs()(zeta) + gamma)
    #     * accumulator_shift_zeta
    # ), "gz != (g(x) + β * Sσ(x) + γ) * Z(xω)"

    # assert (fz.to_coeffs()(zeta) - gz.to_coeffs()(zeta)) * alpha + (
    #     L_1_ext.to_coeffs()(zeta)
    #     * (accumulator_poly.to_coeffs()(zeta) - ONE.to_coeffs()(zeta))
    # ) * alpha**2 == tz.to_coeffs()(
    #     zeta
    # ), "tz != α * ((f(x) + β * Sid(x) + γ) * Z(x) - (g(x) + β * Sσ(x) + γ) * Z(xω)) + α^2 (L1(x) * (Z(x) - 1)))"

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

    # open sid, ssigma at zeta
    proof["Sid"]["proof"], proof["Sid"]["zeta_value"] = kzg.open(Sid.to_coeffs(), zeta)
    proof["Ssigma"]["proof"], proof["Ssigma"]["zeta_value"] = kzg.open(
        Ssigma.to_coeffs(), zeta
    )

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

    print("Proof created")
    return proof


def verify(proof, kzg):
    # polys represented with n points
    omega = omega_base ** (2**32 // n)
    ROOTS = [omega**i for i in range(n)]
    PolyEvalRep = polynomialsEvalRep(Fp, omega, n)

    # create id domain and permutation domains
    # id domain = [0, 1, 2, ..., n-1]
    id_domain = [Fp(i) for i in range(n)]
    # permutation domain = [σ(0), σ(1), σ(2), ..., σ(n-1)]
    permutation_domain = [Fp(sigma[i]) for i in range(n)]

    # create Sid and Ssigma polynomials
    Sid = PolyEvalRep(ROOTS, id_domain)
    Ssigma = PolyEvalRep(ROOTS, permutation_domain)

    # get commitment Sid and Ssigma polynomials
    sid_commitment = kzg.commit(Sid.to_coeffs())
    ssigma_commitment = kzg.commit(Ssigma.to_coeffs())

    # get random γ from verifier using transcript
    transcript = [proof["f"]["commitment"], proof["g"]["commitment"]]
    beta = random_fp_seeded(str(transcript) + "0")
    gamma = random_fp_seeded(str(transcript) + "1")
    # get random α from verifier using transcript
    transcript += [proof["z"]["commitment"]]
    alpha = random_fp_seeded(str(transcript))
    # get random ζ from verifier using transcript
    transcript += [proof["t"]["commitment"]]
    zeta = random_fp_seeded(str(transcript))

    # verify all opening proofs at zeta
    print("Verifying all opening proofs at zeta")
    assert kzg.verify(
        sid_commitment,  # use preprocessed commitment
        proof["Sid"]["proof"],
        zeta,
        proof["Sid"]["zeta_value"],
    ), "Sid opening proof is invalid"
    assert kzg.verify(
        ssigma_commitment,  # use preprocessed commitment
        proof["Ssigma"]["proof"],
        zeta,
        proof["Ssigma"]["zeta_value"],
    ), "Ssigma opening proof is invalid"
    assert kzg.verify(
        proof["f"]["commitment"], proof["f"]["proof"], zeta, proof["f"]["zeta_value"]
    ), "f opening proof is invalid"
    assert kzg.verify(
        proof["g"]["commitment"], proof["g"]["proof"], zeta, proof["g"]["zeta_value"]
    ), "g opening proof is invalid"
    assert kzg.verify(
        proof["t"]["commitment"], proof["t"]["proof"], zeta, proof["t"]["zeta_value"]
    ), "t opening proof is invalid"
    assert kzg.verify(
        proof["z"]["commitment"], proof["z"]["proof"], zeta, proof["z"]["zeta_value"]
    ), "z opening proof is invalid"
    assert kzg.verify(
        proof["z"]["commitment"],
        proof["z"]["shift_proof"],
        zeta * ROOTS[1],
        proof["z"]["shift_zeta_value"],
    ), "z shift opening proof is invalid"

    # get vanishing polynomial
    ZH = vanishing_poly(n)

    # evaluate vanishing polynomial at zeta
    vanishing_poly_eval = ZH(zeta)

    # evaluate L1 lagrange polynomial L1(x) = (x^n - 1) / (n * (x - 1)) at zeta
    L_1_zeta = (zeta**n - Fp(1)) / (n * (zeta - Fp(1)))

    # verify if
    # t * ZH == α * ((f(x) + β * Sid(x) + γ) * Z(x) - (g(x) + β * Sσ(x) + γ) * Z(xω)) + α^2 (L1(x) * (Z(x) - 1)))
    # at zeta (ζ)
    L1Constrainst = alpha**2 * (
        # L1(x) * (Z(x) - 1)
        L_1_zeta
        * (proof["z"]["zeta_value"] - Fp(1))
    )

    TransitionConstraint = alpha * (
        # (f(x) + β * Sid(x) + γ) * Z(x)
        (proof["f"]["zeta_value"] + beta * proof["Sid"]["zeta_value"] + gamma)
        * proof["z"]["zeta_value"]
        # (g(x) + β * Sσ(x) + γ) * Z(xω)
        - (proof["g"]["zeta_value"] + beta * proof["Ssigma"]["zeta_value"] + gamma)
        * proof["z"]["shift_zeta_value"]
    )

    print("Verifying all constraints at zeta")
    # t * ZH == α * ((f(x) + β * Sid(x) + γ) * Z(x) - (g(x) + β * Sσ(x) + γ) * Z(xω)) + α^2 (L1(x) * (Z(x) - 1)))
    return (
        proof["t"]["zeta_value"] * vanishing_poly_eval
        == TransitionConstraint + L1Constrainst
    )


if __name__ == "__main__":
    # setup with degree 100 (all polynomial must be of degree < 100)
    _, CRS = generate_CRS(100)

    # create KZG object
    kzg = KZG(CRS)

    print("xs", xs)
    print("ys", ys)

    # create proof
    proof = prove(kzg)

    # verify proof
    assert verify(proof, kzg), "Verify proof: FAIL"
    print("Verify proof: OK")
