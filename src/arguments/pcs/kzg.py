from arguments.setup import G, generate_CRS, random_fp, evaluate_in_exponent
from utils.ssbls12 import Poly, Fp
from functools import reduce


# quotient polynomial
# q = (f(x) - b) / (x - a))
def quotient(p, a, b):
    return (p - Poly([b])) / (Poly([0, 1]) - Poly([a]))


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


class KZG:
    """
    Usage
    -----
    kzg = KZG.with_degree(n)
    commitment = kzg.commit(f)
    proof, b = kzg.open(f, a)
    assert verify(commitment, proof, a, b)
    """

    def __init__(self, CRS):
        self.CRS = CRS

    @classmethod
    def with_degree(cls, n):
        _, CRS = generate_CRS(n)
        return cls(CRS)

    def commit(self, poly):
        return evaluate_in_exponent(self.CRS, poly)

    def open(self, poly, a):
        b = poly(a)
        q = quotient(poly, a, b)
        proof = self.commit(q)
        return proof, b

    def verify(self, commitment, proof, a, b):
        # [x] = g^a
        x = G * a
        # [y] = g^b
        y = G * b
        # [s - x] = smx
        smx = self.CRS[1] + -x
        # [commitment - y] = cmy
        cmy = commitment + -y
        # check pairing equation between smx and cmy and proof
        return cmy.pair(G * 1) == proof.pair(smx)


# f = main polynomial
# q = quotient polynomial
# I = interpolation polynomial using [a1, a2, ...] and [b1, b2, ...] with lagrange interpolation
# Z = zero polynomial using [a1, a2, ...] = (x - a1)(x - a2)...
#
# CRS = [g, g^tau, g^(tau^2), ...]
# [s] = g^tau = CRS[1]
# commitment = f(CRS)
#
# f(a1) = b1, f(a2) = b2, ...
# I(x) = langrange_interpolation([a1, a2, ...], [b1, b2, ...])
# Z(x) = (x - a1)(x - a2)...
# q = (f(x) - I(x)) / Z(x)
#
# proof = q(CRS)
# Prover will send proof over to verifier
#
# [I(x)] = I(CRS)
# [Z(x)] = Z(CRS)
#
# Verify the pairing equation
#
# e([commitment - I(CRS)], [1]) == e([proof],  [Z(CRS)])
#


class BatchedKZG:
    """
    Usage
    -----
    batchedKzg = BatchedKZG.with_degree(n)
    commitment = batchedKzg.commit(f)
    proof, [b1, b2, ...] = batchedKzg.open(f, [a1, a2, ...])
    assert verify(commitment, proof, [a1, a2, ...], [b1, b1, ...])
    """

    def __init__(self, CRS):
        self.CRS = CRS

    @classmethod
    def with_degree(cls, n):
        _, CRS = generate_CRS(n)
        return cls(CRS)

    def commit(self, poly):
        return evaluate_in_exponent(self.CRS, poly)

    def open(self, poly, ax):
        # evaluate polynomial f(x) at [ax]
        bx = [poly(a) for a in ax]
        # create interpolation polynomial I(X) using [ax] and [bx]
        ip = Poly.interpolate(ax, bx)
        # create zero polynomial Z(X) using [ax]
        zp = reduce(lambda x, y: x * y, [(Poly([0, 1]) - Poly([a])) for a in ax])
        # quotient polynomial q = (f(x) - I(x)) / Z(x)
        q = (poly - ip) / zp
        # create proof
        proof = self.commit(q)
        return proof, bx

    def verify(self, commitment, proof, ax, bx):
        # create interpolation polynomial I(X) using [ax] and [bx]
        ip = Poly.interpolate(ax, bx)
        # create zero polynomial Z(X) using [ax]
        zp = reduce(lambda x, y: x * y, [(Poly([0, 1]) - Poly([a])) for a in ax])
        # [I(x)] = I(CRS)
        ipCRS = evaluate_in_exponent(self.CRS, ip)
        # [Z(x)] = Z(CRS)
        zpCRS = evaluate_in_exponent(self.CRS, zp)

        # [commitment - I(CRS)] = cmi
        cmi = commitment + -ipCRS
        # check pairing equation between cmi and proof and zpCRS
        return cmi.pair(G * 1) == proof.pair(zpCRS)


if __name__ == "__main__":

    def verify_single_proof():
        # setup with degree 100 (all polynomial must be of degree < 100)
        kzg = KZG.with_degree(100)

        # create main polynomial and commit
        f = Poly([1, 2, 3, 4, 5])

        # random a
        a = random_fp()

        # commit and open
        commitment = kzg.commit(f)
        proof, b = kzg.open(f, a)
        print("Verify: ", kzg.verify(commitment, proof, a, b))

    def verify_batched_proof():
        # create main polynomial and commit
        f = Poly([1, 2, 3, 4, 5])
        # random ax
        ax = [random_fp() for _ in range(10)]

        # commit and open with batched kzg
        batchedKzg = BatchedKZG.with_degree(100)
        commitment = batchedKzg.commit(f)
        proof, bx = batchedKzg.open(f, ax)
        print("Batched verify: ", batchedKzg.verify(commitment, proof, ax, bx))

    verify_single_proof()
    verify_batched_proof()
