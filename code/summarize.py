#!/usr/bin/env python3
"""Records a(k) = least prime with mx >= k+1, from all chunk files (and the coverage they prove)."""
import glob
OEIS = {2: 11, 3: 17, 4: 23, 5: 71, 6: 71, 7: 599, 8: 599, 9: 3011, 10: 27901, 11: 52163, 12: 778699, 13: 2374649, 14: 10428007, 15: 34614961}
rows = sorted((int(l.split()[0]), int(l.split()[1])) for l in open('progress.txt') if l[0] != '#')
cov = int(open("START").read())
for a, b in rows:
    assert a == cov, f"gap at {cov}"; cov = b
best = {}
for f in sorted(glob.glob('chunks/c_*.txt'), key=lambda s: int(s.split('_')[-1][:-4])):
    for l in open(f):
        if l[0] == '#': continue
        p, mx = map(int, l.split()[:2])
        for k in range(2, mx):
            best.setdefault(k, p)
best = {k: min(v, best[k]) for k, v in best.items()}
print(f"covered every prime in [{int(open('START').read()):,}, {cov:,})")
for k in sorted(best):
    o = OEIS.get(k); tag = "matches OEIS" if o == best[k] else ("NEW" if o is None else f"MISMATCH (OEIS {o})")
    print(f"a({k}) = {best[k]:,}  {tag}")
