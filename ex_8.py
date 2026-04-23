# Prime Factorization Optimizationに挑戦してみる
# まず素数と仲良くなる
import math

from utils import timeit


# この関数はエラトステネスの篩を使ってnまでの素数を計算します
@timeit
def calc_primes(n: int) -> list[int]:
    if n < 2:
        raise RuntimeError("nは2以上の値である必要があります")

    # 奇数のみのテーブルを作成する
    size = (n + 1) // 2
    table = [True] * size

    table[0] = False  # 1は素数ではない
    p = 3

    while p <= math.isqrt(n):
        pp = p * p // 2
        # 奇数テーブル上でp移動することは2p移動することに相当する
        # 奇数p移動するときそのインデックスは常に偶数になるので無視できる
        table[pp:size:p] = [False] * ((size - pp - 1) // p + 1)

        while True:
            p += 2

            if p > n or table[p // 2]:
                break

    result = [2]

    for i, b in enumerate(table):
        if b:
            result.append(2 * i + 1)

    return result


if __name__ == "__main__":
    # 素数の数(検算用)：
    #  10**2まで25個
    #  10**4まで1229個
    #  10**8まで5761455個

    # 10**9までの素数を求めるのに30秒ほどかかります
    print(f"primes: {len(calc_primes(10**9))}")
