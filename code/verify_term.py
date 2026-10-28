#!/usr/bin/env python3
"""Independent check of one prime for OEIS A060427 / Erdos problem 1056.
Written separately from the scan programs (f1056c.c, f1056p.c) and the reference (ref1056.c).

Usage: python3 verify_term.py p [k]
  Computes m! mod p for m = 1..p-1, finds the most repeated value, lists every m where it occurs,
  checks that each block (m_i, m_{i+1}] multiplies to 1 mod p (the interval form of the problem),
  and runs Wilson's theorem as a sanity check of the loop.
  If k is given, also reports whether p gives k intervals (k+1 equal factorials).
Memory: about 8*p bytes (1.8 GB at p = 2.2e8). Time: a few minutes at p ~ 2e8."""
import sys, time
import numpy as np


def is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17):          # deterministic for n < 3.4e14
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


p = int(sys.argv[1]); k = int(sys.argv[2]) if len(sys.argv) > 2 else None
assert is_prime(p), f"{p} is not prime"
t = time.time()
f = np.empty(p - 1, dtype=np.int64)          # f[m-1] = m! mod p
x = 1
for m in range(1, p):
    x = x * m % p
    f[m - 1] = x
vals, counts = np.unique(f, return_counts=True)
mx = int(counts.max()); top = [int(v) for v in vals[counts == mx]]
print(f"p = {p}: most repeated value(s) of m! mod p: {top}, each {mx} times ({time.time() - t:.0f} s)")
v = top[0]
ms = [int(i) + 1 for i in np.nonzero(f == v)[0]]
print(f"m with m! = {v} (mod p): {ms}")
ok = True
for a, b in zip(ms, ms[1:]):
    prod = 1
    for n in range(a + 1, b + 1):
        prod = prod * n % p
    ok &= (prod == 1)
print(f"each of the {len(ms) - 1} blocks (m_i, m_(i+1)] multiplies to 1 mod p: {ok}")
print("Wilson check (p-1)! = -1 mod p:", int(f[-1]) == p - 1)
if k is not None:
    print(f"gives {k} intervals ({k + 1} equal factorials): {mx >= k + 1}")
