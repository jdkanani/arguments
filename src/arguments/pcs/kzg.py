from arguments.setup import G, G2, GT, generate_CRS, random_fp, evaluate_in_exponent
from utils.ssbls12 import Poly, Fp

#
# s = g * tau
# c = commitment = f(s)
# f(a) = b
# x = g * a
# y = g * b
#
# Verify the pairing equation
#
# e([commitment - y], [1]) = e([proof],  [s - x])
#    equivalent to
# e([commitment - y]^(-1), [1]) * e([proof],  [s - x]) = 1_T
#


# polynomial
def f():
    # f(2) = 57
    return Poly([1, 2, 3, 4, 5])


# quotien polynomial
# q = (f(x) - b) / (x - a))
def q(a, b):
    return (f() - Poly([b])) / (Poly([0, 1]) - Poly([a]))


if __name__ == "__main__":
    # main polynomial
    cp = f()

    # quotient polynomial at a = 2, b = 57; f(a) = b
    a = Fp(2)
    b = cp(a)
    qp = q(a, b)

    # create CRS/tau
    tau, CRS = generate_CRS(cp.degree())
    # [commitment] = f(CRS)
    commitment = evaluate_in_exponent(CRS, cp)
    # [proof] = q(CRS)
    proof = evaluate_in_exponent(CRS, qp)
    # s = G * tau
    s = CRS[1]

    # [s - x] = smx
    smx = s + -G * a

    # [commitment - y] = cmy
    cmy = commitment + -G * b

    # check pairing equation between smx and cmy and proof
    assert cmy.pair(G * 1) == proof.pair(smx)
