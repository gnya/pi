# 円周率を1000000桁計算してみる（チュドノフスキー・アルゴリズム）
import math

from utils import check_pi, timeit

DIGITS_PER_ITER = 14.1816474627254776555
CONST_A = 545140134
CONST_B = 13591409
CONST_C = int(320160**3)
CONST_D = 2 * CONST_C


# この関数は分割統治法でチュドノフスキー・アルゴリズムの分母分子を計算します
def calc_pqt(start: int, end: int) -> tuple[int, int, int]:
    if end < start:
        raise RuntimeError("endはstart以上の値である必要があります")
    elif end > start:
        middle = int((end - start) / 2) + start
        p_l, q_l, t_l = calc_pqt(start, middle)
        p_r, q_r, t_r = calc_pqt(middle + 1, end)

        return p_l * p_r, q_l * q_r, t_l * q_r + t_r * p_l

    # end == startであることに注意
    if start == 0:
        return 1, 1, CONST_B
    else:
        n = start
        p = -(6 * n - 1) * (6 * n - 3) * (6 * n - 5)
        q = CONST_C * n**3
        t = p * (CONST_A * n + CONST_B)

        return p, q, t


@timeit
def calc_pi(digits: int) -> int:
    _, q, t = calc_pqt(0, int(digits / DIGITS_PER_ITER) + 5)

    # 小数の計算は遅いのですべてintで計算する
    sqrt_d = math.isqrt(CONST_D * (10 ** (2 * digits)))

    return (sqrt_d * q) // (t * 6)


# 100万桁の計算に40秒ほどかかる
pi = calc_pi(1000000)

# intをstrに変換するのに時間がかかるのでcheck_piはDecimalを使った場合より遅い
match_digits = check_pi(pi)

print(f"digits: {match_digits}")
