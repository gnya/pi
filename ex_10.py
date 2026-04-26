# 円周率を10^8桁計算してみる（チュドノフスキー・アルゴリズム）
# 参考: https://bellard.org/pi/pi2700e9/pipcrecord.pdf
from multiprocessing import Pool
from typing import Callable

import gmpy2
from gmpy2 import mpfr, mpz

from utils import check_pi, split_range, timeit

PQT = tuple[mpz, mpz, mpz]


class PiCalculator:
    DIGITS_PER_ITER = 14.1816474627254776555

    CONST_A = 545140134
    CONST_B = 13591409
    CONST_C = int(640320**3 / 24)
    CONST_D = int(640320**3 / 16)

    def __init__(self, digits: int, n_jobs: int = 8, n_merge_jobs: int = 2) -> None:
        if n_jobs < 1:
            raise RuntimeError("n_jobsには1以上の値を設定してください")

        if n_merge_jobs < 1:
            raise RuntimeError("n_merge_jobsには1以上の値を設定してください")
        elif n_merge_jobs > n_jobs:
            raise RuntimeError("n_merge_jobsにはn_jobs以下の値を設定してください")

        self.digits = digits
        self.n_jobs = n_jobs
        self.n_merge_jobs = n_merge_jobs

        self.max_loop = int(digits / self.DIGITS_PER_ITER) + 5
        self.precision = self.digits * 5

        self._pqt_list = []

    # この関数はアルゴリズムのn項目のPQTを返します
    def pqt_term(self, n: int) -> PQT:
        a = mpz((self.CONST_A * n + self.CONST_B))
        p = mpz((2 * n - 1) * (6 * n - 1) * (6 * n - 5))
        q = mpz(self.CONST_C * n**3)

        if n % 2 == 0:
            t = a * p
        else:
            t = -a * p

        return (p, q, t)

    # この関数はn要素目の計算済みPQTを返します
    def pqt_list(self, n: int) -> PQT:
        return self._pqt_list[n]

    # この関数はアルゴリズムのPQTの合成を行います
    def merge_pqt(self, pqt_l: PQT, pqt_r: PQT) -> PQT:
        return (
            pqt_l[0] * pqt_r[0],
            pqt_l[1] * pqt_r[1],
            pqt_l[2] * pqt_r[1] + pqt_r[2] * pqt_l[0],
        )

    # この関数は2分法でアルゴリズムのPQTを計算します
    # NOTE 範囲はstart項目からend-1項目までになります
    def calc_pqt(self, start: int, end: int, term: Callable[int, PQT]) -> PQT:
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
            args = [(s, e, self.pqt_term) for s, e in args]
            self._pqt_list = pool.starmap(self.calc_pqt, args)

        with Pool(self.n_merge_jobs) as pool:
            args = split_range(0, self.n_jobs, self.n_merge_jobs)
            args = [(s, e, self.pqt_list) for s, e in args]
            self._pqt_list = pool.starmap(self.calc_pqt, args)

        return self.calc_pqt(0, len(self._pqt_list), self.pqt_list)

    # この関数は円周率を計算します
    @timeit
    def calc_pi(self) -> mpfr:
        _, num, t = self.calc_all_pqt()
        den = mpz(3) * (t + num * mpz(self.CONST_B))

        with gmpy2.context(precision=self.precision):
            # 逆数平方根のほうが速い
            pi = num / (den * gmpy2.rec_sqrt(mpz(self.CONST_D)))

        return pi


if __name__ == "__main__":
    import ex_7

    n = 10**7
    calc = PiCalculator(n)

    """
    timer = time_start("pqt_term (ex_10)")
    for i in range(n):
        calc.pqt_term(i + 1)
    timer.stop()

    timer = time_start("pqt_term (ex_7)")
    for i in range(n):
        ex_7.pqt_term(i + 1)
    timer.stop()
    """

    # 古い方の実装
    ex_7.calc_pi_old(n, 8)

    # 1億桁の計算に40秒ほどかかる
    pi = calc.calc_pi()

    # 今の方法だと500万桁を超える円周率をチェックできない…
    print(f"digits: {check_pi(pi)}")
