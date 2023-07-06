def is_power_of_two(n):
    # bit-arithmetic trick
    return n & (n - 1) == 0


def nearest_power_of_two(n):
    if is_power_of_two(n):
        return n
    return 2 ** n.bit_length()


if __name__ == "__main__":
    res = is_power_of_two(8)
    print(res)
    res = nearest_power_of_two(12)
    print(res)
