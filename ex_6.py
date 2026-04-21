# 円周率を10000000桁計算してみる（チュドノフスキー・アルゴリズム）
# Python標準では限界がきたのでgmpy2を使っています
import gmpy2
from gmpy2 import mpfr, mpz

from utils import check_pi, timeit

DIGITS_PER_ITER = 14.1816474627254776555
CONST_A = 545140134
CONST_B = 13591409
CONST_C = int(320160**3)
CONST_D = 2 * CONST_C


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


@timeit
def calc_pi(digits: int) -> mpfr:
    _, q, t = calc_pqt(0, int(digits / DIGITS_PER_ITER) + 5)

    with gmpy2.context(precision=digits * 5):
        sqrt_d = gmpy2.sqrt(mpz(CONST_D)).real
        pi = (sqrt_d * q) / (mpz(6) * t)

    return pi


# 1000万桁の計算に6秒ほどかかる
pi = calc_pi(10000000)

# 今の方法だと500万桁を超える円周率をチェックできない…
match_digits = check_pi(pi)

print(f"digits: {match_digits}")
