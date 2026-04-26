# Prime Factorization Optimizationに挑戦してみる
# 素因数分解をやってみる
import math
from collections import defaultdict

import numpy as np
import numpy.typing as npt
from gmpy2 import mpz

from utils import timeit

LPF = npt.NDArray[np.unsignedinteger]
Factors = defaultdict[mpz, int]


# 与えられた辞書からFactorsを初期化します
def factors(f: dict[int, int] = {}) -> Factors:
    return defaultdict(int, {(mpz(k), v) for k, v in f.items()})


# エラトステネスの篩を使ってnまでの最小素因数のテーブルを計算します
@timeit
def calc_lpf_table(n: int, dtype: type[np.unsignedinteger] = np.uint32) -> LPF:
    if n < 2:
        raise RuntimeError("nは2以上の値である必要があります")
    elif n > np.iinfo(dtype).max:
        raise RuntimeError("nがdtypeで表せる範囲を超えています")

    # 奇数のみのテーブルを作成する
    table = np.arange(1, n + 1, 2, dtype)
    size = len(table)

    p = 3

    while p <= math.isqrt(n):
        # 奇数テーブル上でp移動することは2p移動することに相当する
        # 奇数p移動するときそのインデックスは常に偶数になるので無視できる
        pp = p * p // 2
        table[pp:size:p] = p
        last_p = p

        while True:
            p += 2

            if p > n or last_p < table[p // 2]:
                break

    return table


# nの素因数分解を計算します
def factorize(n: mpz, table: LPF) -> Factors:
    f = factors()

    if n < 0:
        # -1を符号を扱うために使います
        f[mpz(-1)] += 1
        n = -n

    while n & 1 == 0:
        f[mpz(2)] += 1
        n >>= 1

    if n > len(table) * 2:
        # 最小素因数のテーブルより大きな値だった場合の処理
        # TODO 必要ならフォールバックの処理を書く
        f[n] += 1
        print("n > len(table):", n)
    else:
        while n > 1:
            p = mpz(table[n // 2])
            f[p] += 1
            n //= p

    return f


# 素因数の辞書から元の値を復元します
def to_int(f: Factors) -> mpz:
    result = mpz(1)

    for p, n in f.items():
        result *= p ** mpz(n)

    return result


# 素因数の辞書で表された値同士の積を求めます
def mul(a: Factors, b: Factors) -> Factors:
    c = a.copy()

    for p, n in b.items():
        c[p] += n

    if mpz(-1) in c:
        # 符号を整理する
        c[mpz(-1)] %= 2

        if c[mpz(-1)] == 0:
            del c[mpz(-1)]

    return c


# 素因数の辞書で表された値同士の和を求めます
def add(a: Factors, b: Factors, table: LPF):
    c = factors()
    c_a = a.copy()
    c_b = b.copy()

    # 共通部分を括りだす
    for p in a.keys() & b.keys():
        c[p] += min(a[p], b[p])
        c_a[p] -= c[p]

        if c_a[p] == 0:
            del c_a[p]

        c_b[p] -= c[p]

        if c_b[p] == 0:
            del c_b[p]

    c_c = to_int(c_a) + to_int(c_b)

    return mul(c, factorize(c_c, table))


if __name__ == "__main__":
    # 素数の数(検算用)：
    #  10**2まで25個
    #  10**4まで1229個
    #  10**8まで5761455個

    # 最小素因数のテーブルを計算する
    table = calc_lpf_table(10**9)

    """
    # テーブルをキャッシュして使用するなら以下のコードを参考にする
    # NOTE 例えば10^10までのテーブルの容量は20GB近くになる
    import os
    import pickle

    path = f"{os.path.dirname(__file__)}\\lpf_table.pkl"

    with open(path, "wb") as f:
        pickle.dump(table, f, protocol=pickle.HIGHEST_PROTOCOL)

    with open(path, "rb") as f:
       table = pickle.load(f)
    """

    n = 5761455
    print(f"factorize({n}): {factorize(n, table)}")
