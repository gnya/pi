# 円周率を10^8桁計算してみる（チュドノフスキー・アルゴリズム）
# 参考: https://bellard.org/pi/pi2700e9/pipcrecord.pdf
import math
from multiprocessing import Pool
from typing import Callable

import gmpy2
from gmpy2 import mpfr, mpz

from utils import split_range, timeit

PQT = tuple[mpz, mpz, mpz]


class PiCalculator:
    BITS_PER_DIGIT = 3.32192809488736234787
    DIGITS_PER_ITER = 14.1816474627254776555

    CONST_A = mpz(545140134)
    CONST_B = mpz(13591409)
    CONST_C = mpz(640320**3 / 24)
    CONST_D = mpz(640320**3 / 144)

    def __init__(self, depth: int, n_jobs: int = 8, n_merge_jobs: int = 2) -> None:
        if n_jobs < 1:
            raise RuntimeError("n_jobsには1以上の値を設定してください")

        min_depth = 0 if n_jobs == 1 else (n_jobs - 2).bit_length()

        if depth < min_depth:
            raise RuntimeError(f"depthには{min_depth}以上の値を設定してください")
        elif n_merge_jobs < 1:
            raise RuntimeError("n_merge_jobsには1以上の値を設定してください")
        elif n_merge_jobs > n_jobs:
            raise RuntimeError("n_merge_jobsにはn_jobs以下の値を設定してください")

        self.max_loop = 1 + (1 << depth)
        self.n_jobs = n_jobs
        self.n_merge_jobs = n_merge_jobs

        digits = self.max_loop * self.DIGITS_PER_ITER
        self.precision = int(math.ceil(digits * self.BITS_PER_DIGIT))

        self._pqt_list: list[PQT] = []

    # この関数はアルゴリズムのn項目のPQTを返します
    def pqt_term(self, n: mpz) -> PQT:
        a = self.CONST_A * n + self.CONST_B
        p = (2 * n - 1) * (6 * n - 1) * (6 * n - 5)
        q = self.CONST_C * n**3

        if n % 2 == 0:
            t = a * p
        else:
            t = -a * p

        return (p, q, t)  # type: ignore

    # この関数はn要素目の計算済みPQTを返します
    def pqt_list(self, n: mpz) -> PQT:
        return self._pqt_list[int(n)]

    # この関数はアルゴリズムのPQTの合成を行います
    def merge_pqt(self, a: PQT, b: PQT) -> PQT:
        return (
            a[0] * b[0],
            a[1] * b[1],
            a[2] * b[1] + b[2] * a[0],
        )

    # この関数は2分法でアルゴリズムのPQTを計算します
    # NOTE 範囲はstart項目からend-1項目までになります
    def calc_pqt(self, start: mpz, end: mpz, term: Callable[mpz, PQT]) -> PQT:
        if end <= 0 or start > end - 1:
            raise RuntimeError("endとstartの値が正しくありません:", start, end)
        elif start == end - 1:
            return term(start)
        else:
            middle = start + (end - start) // 2

            return self.merge_pqt(
                self.calc_pqt(start, middle, term),
                self.calc_pqt(middle, end, term),
            )

    # この関数はアルゴリズムのPQTを計算します
    @timeit
    def calc_all_pqt(self) -> PQT:
        with Pool(self.n_jobs) as pool:
            args = split_range(1, self.max_loop, self.n_jobs)
            args = [(mpz(s), mpz(e), self.pqt_term) for s, e in args]
            self._pqt_list = pool.starmap(self.calc_pqt, args)

        with Pool(self.n_merge_jobs) as pool:
            args = split_range(0, self.n_jobs, self.n_merge_jobs)
            args = [(mpz(s), mpz(e), self.pqt_list) for s, e in args]
            self._pqt_list = pool.starmap(self.calc_pqt, args)

        return self.calc_pqt(0, len(self._pqt_list), self.pqt_list)

    # この関数は円周率を計算します
    @timeit
    def calc_pi(self) -> mpfr:
        with gmpy2.context(precision=self.precision):
            _, num, den = self.calc_all_pqt()
            pi = num / ((den + num * self.CONST_B) * gmpy2.rec_sqrt(self.CONST_D))

        return pi


if __name__ == "__main__":
    import os

    k = 23
    calc = PiCalculator(k)

    print("digits:", int(2**k * calc.DIGITS_PER_ITER))

    # 1億桁の計算に40秒ほどかかる
    pi = calc.calc_pi()
    str_pi = str(pi)

    path = f"{os.path.dirname(__file__)}\\pi_{k}.txt"

    with open(path, "w") as f:
        # 計算した円周率をテキストに書き出す
        f.write(str_pi)

    # 最後から20桁を表示する
    # TODO 500万桁を超える円周率をチェックできるようにする
    print(f"pi {len(str_pi) - 1 - 20}-{len(str_pi) - 1 - 1}:", str_pi[-20:])
