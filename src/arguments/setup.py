import random

from utils.ssbls12 import Fp, Poly, Group
from finitefield.polynomial_evalrep import get_omega, polynomialsEvalRep

# e(G,G) I am using a Type1 Bilinear group for simplicity
G = Group.G
G2 = G
GT = Group.GT

# | # Choosing roots of unity
# | The BLS12-381 is chosen in part because it's FFT friendly. To use radix-2
# | FFT, we need to find m^th roots of unity, where m is a power of two, and
# | m is the degree bound of the polynomial we want to represent.
# |
# | In the BLS12-381, we can find primitive n^th roots of unity, for any
# | power of two n up to n <= 2^**32.
# | This follows because for the ssbls12-381 exponent field Fp, we have
# |    2^32 divides (p - 1).

omega_base = get_omega(Fp, 2**32, seed=0)


def random_fp():
    return Fp(random.randint(0, Fp.p - 1))


def random_fp_seeded(seeded):
    random.seed(seeded)
    return Fp(random.randint(0, Fp.p - 1))


# Generate CRS using tau and returns them
def generate_CRS(n):
    tau = random_fp()
    CRS = [G * (tau**i) for i in range(n + 3)]
    return tau, CRS


def vanishing_poly(n):  # degree = 8, coeffs = [1, 0, 0, 0, 0, 0, 0, 0, -1] = x^8 - 1
    # For the special case of evaluating at all n powers of omega,
    # the vanishing poly has a special form.
    #  t(X) = (X-omega^0)(X-omega^1)....(X-omega^(n-1)) = X^n - 1
    return Poly([Fp(-1)] + [Fp(0)] * (n - 1) + [Fp(1)])


def eval_poly(poly, domain, shift=Fp(1)):
    poly_coeff = poly.to_coeffs()
    eval = []
    for j in range(len(domain)):
        eval += [
            sum(
                [
                    (domain[j] * shift) ** i * poly_coeff.coefficients[i]
                    for i in range(poly_coeff.degree() + 1)
                ]
            )
        ]
    return eval


# Evaluate a polynomial in exponent
def evaluate_in_exponent(powers_of_tau, poly):
    # powers_of_tau:
    #    [G*0, G*tau, ...., G*(tau**m)]
    # poly:
    #    degree-m bound polynomial in coefficient form
    # print('P.degree:', poly.degree())
    # print('taus:', len(powers_of_tau))
    assert poly.degree() + 1 < len(powers_of_tau)
    return sum(
        [powers_of_tau[i] * poly.coefficients[i] for i in range(poly.degree() + 1)],
        G * 0,
    )


############################
# PLONK related setup
############################


def setup_for_plonk(gates_matrix, permutation, L, p_i):
    print("Starting Setup Phase...")
    (m, n) = gates_matrix.shape

    assert n & n - 1 == 0, "n must be a power of 2"
    omega = omega_base ** (2**32 // n)
    ROOTS = [omega**i for i in range(n)]

    PolyEvalRep = polynomialsEvalRep(Fp, omega, n)

    # Generate polynomials from columns of gates_matrix
    q_L = PolyEvalRep(ROOTS, [Fp(i) for i in gates_matrix[0]])
    q_R = PolyEvalRep(ROOTS, [Fp(i) for i in gates_matrix[1]])
    q_M = PolyEvalRep(ROOTS, [Fp(i) for i in gates_matrix[2]])
    q_O = PolyEvalRep(ROOTS, [Fp(i) for i in gates_matrix[3]])
    q_C = PolyEvalRep(ROOTS, [Fp(i) for i in gates_matrix[4]])
    Qs = [q_L, q_R, q_M, q_O, q_C]

    # The public input poly vanishes everywhere except for the position of the
    # public input gate where it evaluates to -(public_input)
    public_input = [Fp(0) for i in range(len(ROOTS))]
    for i in L:
        public_input[i] = Fp(-p_i)
    p_i_poly = PolyEvalRep(ROOTS, public_input)

    # We generate domains on which we can evaluate the witness polynomials
    k = random_fp()
    id_domain_a = ROOTS
    id_domain_b = [k * root for root in ROOTS]
    id_domain_c = [k**2 * root for root in ROOTS]
    id_domain = id_domain_a + id_domain_b + id_domain_c

    # We permute the positions of the domain generated above
    perm_domain = [id_domain[i - 1] for i in permutation]
    perm_domain_a = perm_domain[:n]
    perm_domain_b = perm_domain[n : 2 * n]
    perm_domain_c = perm_domain[2 * n : 3 * n]

    # Generate polynomials that return the permuted index when evaluated on the
    # domain
    S_sigma_1 = PolyEvalRep(ROOTS, perm_domain_a)
    S_sigma_2 = PolyEvalRep(ROOTS, perm_domain_b)
    S_sigma_3 = PolyEvalRep(ROOTS, perm_domain_c)
    Ss = [S_sigma_1, S_sigma_2, S_sigma_3]

    perm_precomp = [id_domain, perm_domain, k, Ss]

    # We perform the trusted setup
    tau, CRS = generate_CRS(n)

    # We take some work off the shoulders of the verifier
    print("Starting Verifier Preprocessing...")
    q_exp = [evaluate_in_exponent(CRS, q.to_coeffs()) for q in Qs]
    s_exp = [evaluate_in_exponent(CRS, s.to_coeffs()) for s in Ss]
    x_exp = G2 * tau
    verifier_preprocessing = [q_exp, s_exp, x_exp]
    print("Setup Phase Finished!")

    return CRS, Qs, p_i_poly, perm_precomp, verifier_preprocessing
