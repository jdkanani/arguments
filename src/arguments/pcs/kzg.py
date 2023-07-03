from arguments.setup import G, G2, GT, generate_CRS, random_fp, evaluate_in_exponent
from utils.ssbls12 import Poly, Fp

# f = main polynomial
# q = quotient polynomial
#
# CRS = [g, g^tau, g^(tau^2), ...]
# [s] = g^tau = CRS[1]
# commitment = f(CRS)
#
# f(a) = b
# q = (f(x) - b) / (x - a))
# proof = q(CRS)
# [x] = g^a
# [y] = g^b
#
# Verify the pairing equation
#
# e([commitment - y], [1]) == e([proof],  [s - x])
#    equivalent to
# e([commitment - y]^(-1), [1]) * e([proof],  [s - x]) == 1_T
#


# quotien polynomial
# q = (f(x) - b) / (x - a))
def quotient(p, a, b):
    return (p - Poly([b])) / (Poly([0, 1]) - Poly([a]))


# setup CRS
def setup(n):
    _, CRS = generate_CRS(n)
    return CRS


# commit polynomial
def commit(poly, crs):
    return evaluate_in_exponent(crs, poly)


# open polynomial at a
def open(poly, crs, a):
    b = poly(a)
    q = quotient(poly, a, b)
    proof = commit(q, crs)
    return proof, b


def verify(commitment, proof, a, b, crs):
    # [x] = g^a
    x = G * a
    # [y] = g^b
    y = G * b
    # [s - x] = smx
    smx = crs[1] + -x
    # [commitment - y] = cmy
    cmy = commitment + -y
    # check pairing equation between smx and cmy and proof
    return cmy.pair(G * 1) == proof.pair(smx)


if __name__ == "__main__":
    # setup with degree 100 (all polynomial must be of degree < 100)
    CRS = setup(100)

    # create main polynomial and commit
    f = Poly([1, 2, 3, 4, 5])
    commitment = commit(f, CRS)

    a = random_fp()
    # show example how to use commit, open, and verify
    proof, b = open(f, CRS, a)
    print("Verify: ", verify(commitment, proof, a, b, CRS))
