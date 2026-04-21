# 円周率を10000桁計算してみる（マチンの公式）
import math
from decimal import Decimal, localcontext

from utils import check_pi, timeit


@timeit
def calc_pi(digits: int, max_loop: int) -> Decimal | None:
    with localcontext() as context:
        context.prec = digits + 5

        # 最大ループ回数が必ず100の倍数+1になるようにする
        max_loop = math.ceil(max_loop / 100) * 100 + 1

        A = Decimal(5)
        B = Decimal(239)
        x = Decimal(0)

        for n in range(max_loop):
            m = 2 * Decimal(n) + 1
            delta_x = (4 / A**m - 1 / B**m) / m

            if n % 2 == 0:
                x += delta_x
            else:
                x -= delta_x

            # check_piは毎回pi.txtを読むので100回に1回だけチェックする
            if n % 100 == 0:
                pi = 4 * x
                match_digits = check_pi(pi)

                print(f"n: {n + 1}, digits: {match_digits}")

                if match_digits >= digits:
                    return pi

    return None


# 1万桁の計算に30秒ほどかかる
calc_pi(10000, 10000)
