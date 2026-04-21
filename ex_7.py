# 円周率を100000000桁計算してみる（チュドノフスキー・アルゴリズム）
# Python標準では限界がきたのでgmpy2を使っています
from multiprocessing import Pool

import gmpy2
from gmpy2 import mpfr, mpz

from utils import check_pi, timeit

DIGITS_PER_ITER = 14.1816474627254776555
CONST_A = 545140134
CONST_B = 13591409
CONST_C = int(320160**3)


# この関数はチュドノフスキー・アルゴリズムの分母分子の合成を行います
def comp_pqt(
    p_l: mpz, q_l: mpz, t_l: mpz, p_r: mpz, q_r: mpz, t_r: mpz
) -> tuple[mpz, mpz, mpz]:
    return p_l * p_r, q_l * q_r, t_l * q_r + t_r * p_l


# この関数は2分法でチュドノフスキー・アルゴリズムの分母分子を計算します
def calc_pqt(start: int, end: int) -> tuple[mpz, mpz, mpz]:
    if end < start:
        raise RuntimeError("endはstart以上の値である必要があります")
    elif end > start:
        middle = int((end - start) * 0.5) + start

        return comp_pqt(
            *calc_pqt(start, middle),
            *calc_pqt(middle + 1, end),
        )

    # end == startであることに注意
    if start == 0:
        return mpz(1), mpz(1), mpz(CONST_B)
    else:
        n = start
        p = -(6 * n - 1) * (6 * n - 3) * (6 * n - 5)
        q = CONST_C * n**3
        t = p * (CONST_A * n + CONST_B)

        return mpz(p), mpz(q), mpz(t)


# この関数はチュドノフスキー・アルゴリズムの分母分子を並列で計算します
def calc_pqt_multiprocess(max_loop: int, processes: int) -> tuple[mpz, mpz, mpz]:
    if processes < 1:
        raise RuntimeError("processesには1以上の値を設定してください")

    with Pool(processes) as pool:
        results = pool.starmap(
            calc_pqt,
            [
                (
                    (max_loop * n) // processes,
                    (max_loop * (n + 1)) // processes - 1,
                )
                for n in range(processes)
            ],
        )

    p, q, t = results[0]

    for result_p, result_q, result_t in results[1:]:
        p, q, t = comp_pqt(p, q, t, result_p, result_q, result_t)

    return p, q, t


@timeit
def calc_pi(digits: int, processes: int) -> mpfr:
    max_loop = int(digits / DIGITS_PER_ITER) + 5
    _, q, t = calc_pqt_multiprocess(max_loop, processes)

    with gmpy2.context(precision=digits * 5):
        sqrt_d = gmpy2.sqrt(mpz(2 * CONST_C)).real
        pi = (sqrt_d * q) / (mpz(6) * t)

    return pi


if __name__ == "__main__":
    # 1億桁の計算に60秒ほどかかる
    pi = calc_pi(100000000, 8)

    # 今の方法だと500万桁を超える円周率をチェックできない…
    match_digits = check_pi(pi)

    print(f"digits: {match_digits}")
