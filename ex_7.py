# 円周率を10^8桁計算してみる（チュドノフスキー・アルゴリズム）
from multiprocessing import Pool

import gmpy2
from gmpy2 import mpfr, mpz

from utils import check_pi, split_range, timeit

PQT = tuple[mpz, mpz, mpz]

DIGITS_PER_ITER = 14.1816474627254776555
CONST_A = 545140134
CONST_B = 13591409
CONST_C = int(320160**3)
CONST_D = 2 * CONST_C


# この関数はアルゴリズムの分母分子の合成を行います
def merge_pqt(pqt_l: PQT, pqt_r: PQT) -> PQT:
    return (
        pqt_l[0] * pqt_r[0],
        pqt_l[1] * pqt_r[1],
        pqt_l[2] * pqt_r[1] + pqt_r[2] * pqt_l[0],
    )


# この関数はアルゴリズムのn項目の分母分子を返します
def pqt_term(n: int) -> PQT:
    if n == 0:
        return (mpz(1), mpz(1), mpz(CONST_B))
    else:
        p = -(6 * n - 1) * (6 * n - 3) * (6 * n - 5)
        q = CONST_C * n**3
        t = p * (CONST_A * n + CONST_B)

        return (mpz(p), mpz(q), mpz(t))


# この関数は2分法でアルゴリズムの分母分子を計算します
def calc_pqt_bs(start: int, end: int) -> PQT:
    if end < start:
        raise RuntimeError("endはstart以上の値である必要があります")
    elif end > start:
        middle = start + (end - start) // 2

        return merge_pqt(
            calc_pqt_bs(start, middle),
            calc_pqt_bs(middle + 1, end),
        )
    else:
        return pqt_term(start)


# この関数は2分法でアルゴリズムの分母分子を合成します
def merge_pqt_list(start: int, end: int, pqt_list: list[PQT]) -> PQT:
    if end < start:
        raise RuntimeError("endはstart以上の値である必要があります")
    elif end > start:
        middle = start + (end - start) // 2

        return merge_pqt(
            merge_pqt_list(start, middle, pqt_list),
            merge_pqt_list(middle + 1, end, pqt_list),
        )
    else:
        return pqt_list[start]


# この関数はアルゴリズムの分母分子を並列で計算します
# NOTE n_merge_jobsの値は特に理由がない限り2に設定するのが良さそう
def calc_pqt(max_loop: int, n_jobs: int, n_merge_jobs: int) -> PQT:
    if n_jobs < 2:
        raise RuntimeError("processesには2以上の値を設定してください")

    with Pool(n_jobs) as pool:
        args = split_range(max_loop, n_jobs)
        results = pool.starmap(calc_pqt_bs, args)

    # 並列化を入れ子にするよりもこちらのほうが高速でした
    with Pool(n_merge_jobs) as pool:
        args = split_range(n_jobs, n_merge_jobs)
        args = [(s, e, results) for s, e in args]
        results = pool.starmap(merge_pqt_list, args)

    p, q, t = merge_pqt_list(0, len(results) - 1, results)

    return p, q, t


@timeit
def calc_pi(digits: int, n_jobs: int) -> mpfr:
    max_loop = int(digits / DIGITS_PER_ITER) + 5
    _, q, t = calc_pqt(max_loop, n_jobs, 2)

    with gmpy2.context(precision=digits * 5):
        # 逆数平方根のほうが速い
        pi = q / (mpz(6) * t * gmpy2.rec_sqrt(mpz(CONST_D)))

    return pi


if __name__ == "__main__":
    # 1億桁の計算に40秒ほどかかる
    pi = calc_pi(100000000, 8)

    # 今の方法だと500万桁を超える円周率をチェックできない…
    print(f"digits: {check_pi(pi)}")
