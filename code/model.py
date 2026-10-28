import numpy as np, math
OEIS = {2: 11, 3: 17, 4: 23, 5: 71, 6: 71, 7: 599, 8: 599, 9: 3011, 10: 27901, 11: 52163, 12: 778699, 13: 2374649, 14: 10428007, 15: 34614961}
N = 400_000_000
s = np.ones(N//2, dtype=bool); s[0]=False  # odd numbers 2i+1
for i in range(1, int(N**0.5)//2+1):
    if s[i]:
        p = 2*i+1; s[p*p//2::p] = False
primes = np.concatenate(([2], 2*np.nonzero(s)[0]+1)).astype(np.int64)
cum = np.cumsum((primes-1).astype(np.float64))
def tail(k):  # P(Poisson(1) >= k+1)
    return 1 - sum(math.exp(-1)/math.factorial(j) for j in range(k+1))
def tailx(k):
    return sum(math.exp(-1)/math.factorial(j) for j in range(k+1, k+40))
print("k  actual a(k)    model median   ratio")
for k in range(2, 17):
    t = tailx(k)
    target = math.log(2)/t
    idx = np.searchsorted(cum, target)
    pred = primes[idx] if idx < len(primes) else float('nan')
    a = OEIS.get(k)
    print(f"{k:2d} {a if a else '?':>12} {pred:>14,} {a/pred if a else float('nan'):6.2f}")
    # 90% interval: E in [-ln(0.95), -ln(0.05)]
    lo = primes[np.searchsorted(cum, -math.log(0.95)/t)]; hi_i = np.searchsorted(cum, -math.log(0.05)/t)
    hi = primes[hi_i] if hi_i < len(primes) else None
    print(f"     90% range {lo:,} .. {hi if hi is None else f'{hi:,}'}")
