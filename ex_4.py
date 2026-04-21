# 円周率を1000000桁計算してみる（ガウス＝ルジャンドルのアルゴリズム）
from decimal import Decimal, localcontext

from utils import check_pi, timeit


# この関数はvalueの逆数平方根と平方根を計算します
@timeit
def calc_rsqrt(value: Decimal, last_rsqrt: Decimal) -> tuple[Decimal, Decimal]:
    x = Decimal(0)
    new_x = last_rsqrt

    while x != new_x:
        x = new_x
        new_x = (3 - value * x * x) * x / 2

    return new_x, new_x * value


@timeit
def calc_pi(digits: int, max_loop: int) -> Decimal | None:
    with localcontext() as context:
        context.prec = digits + 5

        a = Decimal(1)
        inv_b = Decimal(2).sqrt()
        b = 1 / inv_b
        t = Decimal(0.25)
        p = Decimal(1)

        for n in range(max_loop):
            new_a = (a + b) / 2
            # bの値は収束していくのでニュートン法を使ったほうが早く計算できる
            inv_b, b = calc_rsqrt(a * b, inv_b)
            t -= p * (a - new_a) ** 2
            p *= 2

            a = new_a

            pi = (a + b) ** 2 / (4 * t)
            match_digits = check_pi(pi)

            print(f"n: {n + 1}, digits: {match_digits}")

            if match_digits >= digits:
                return pi

    return None


# 100万桁の計算に60秒ほどかかる
calc_pi(1000000, 10000)
