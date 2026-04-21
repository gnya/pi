import os
from decimal import Decimal

from .time import timeit

PI_PATH = f"{os.path.dirname(__file__)}\\pi.txt"


def int_to_str(x: int, n: int = 0) -> str:
    if n == 0:
        n = int(x.bit_length() * 0.3)
    if n < 100:
        return str(x)
    n >>= 1
    l = 10**n
    a, b = divmod(x, l)
    upper = int_to_str(a, n)
    lower = int_to_str(b, n).rjust(n, "0")
    return upper + lower


@timeit
def check_pi(proximate: Decimal | str | int) -> int:
    with open(PI_PATH) as f:
        prox_pi = str(proximate).replace(".", "")
        real_pi = f.read().replace(".", "")

        for i, (prox_pi_digit, real_pi_digit) in enumerate(zip(prox_pi, real_pi)):
            if prox_pi_digit != real_pi_digit:
                return i - 1

        return min(len(prox_pi) - 1, len(real_pi) - 1)
