def isPowerOfTwo(n):
    # bit-arithmetic trick
    return n & (n - 1) == 0


def nearestPowerOfTwo(n):
    if isPowerOfTwo(n):
        return n
    return 2 ** n.bit_length()


if __name__ == "__main__":
    res = isPowerOfTwo(8)
    print(res)
    res = nearestPowerOfTwo(12)
    print(res)
