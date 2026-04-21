# 円周率をまず10桁計算してみる（アルキメデスの手法）
import math

# xは正2^(n+1)角形の1辺の長さの2乗
x = 2.0

for n in range(30):
    a = 1.0 - x / 4.0
    # 元の式では 1 - sqrt(a) が含まれる形になる
    # nが大きくなるとsqrt(a)が1に近付いて誤差が大きくなるので式変形した
    x /= 2.0 * (1.0 + math.sqrt(a))

    pi = math.pow(2.0, n + 2) * math.sqrt(x)
    error = round(math.log10(abs(math.pi - pi)))

    print(f"n: {n + 1}, pi: {pi}, error (log10): {error}")
