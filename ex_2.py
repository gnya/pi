# 円周率を1000桁計算してみる（アルキメデスの手法）
# 浮動小数点型で表せる桁数を超えるのでdecimalを利用する
import math
from decimal import Decimal, localcontext

from utils import check_pi, timeit


@timeit
def calc_pi(digits: int, max_loop: int) -> Decimal | None:
    with localcontext() as context:
        context.prec = digits + 5

        # 最大ループ回数が必ず100の倍数+1になるようにする
        max_loop = math.ceil(max_loop / 100) * 100 + 1

        # xは正2^(n+1)角形の1辺の長さの2乗
        x = Decimal(2.0)

        for n in range(max_loop):
            # 元の式では 1 - sqrt(a) が含まれる形になる
            # nが大きくなるとsqrt(a)が1に近付いて誤差が大きくなるので式変形した
            a = Decimal(1.0) - x / Decimal(4.0)
            x /= Decimal(2.0) * (Decimal(1.0) + a.sqrt())

            # check_piは毎回pi.txtを読むので100回に1回だけチェックする
            if n % 100 == 0:
                pi = (Decimal(2.0) ** Decimal(n + 2)) * x.sqrt()

                print(f"n: {n + 1}, digits: {check_pi(pi)}")

                if check_pi(pi) >= digits:
                    return pi

    return None


calc_pi(1700, 10000)
