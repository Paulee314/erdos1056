#!/usr/bin/env python3
"""Certificate for the #1056 scan: one line per finished block in certify.txt.
For each block listed in progress.txt (not yet certified):
  - sha256 of chunks/c_<P0>.txt, and its '# done' line must match progress.txt (range, primes, max_mx)
  - spot check: re-run a window of ~600 numbers inside the block with ref1056 (plain %, no shortcuts,
    written separately) and compare the lines byte for byte
  - close calls: re-run every prime in the block with 15 or more equal factorials with ref1056
Line: P0 P1 primes max_mx sha256 window_lo window_hi window_primes close_calls result
Usage: certify.py [dir] [ref1056 binary]"""
import hashlib, os, random, subprocess, sys
D = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__)))
REF = os.path.abspath(sys.argv[2] if len(sys.argv) > 2 else os.path.join(D, 'ref1056'))
OUT = os.path.join(D, 'certify.txt')
done = set()
if os.path.exists(OUT):
    done = {l.split()[0] for l in open(OUT) if l.strip() and not l.startswith('#')}
else:
    open(OUT, 'w').write("# P0 P1 primes max_mx sha256 window_lo window_hi window_primes close_calls result\n")
def ref(lo, hi):
    r = subprocess.run([REF, str(lo), str(hi)], capture_output=True, text=True, check=True)
    return [l for l in r.stdout.splitlines() if l.strip()]
rows = [l.split() for l in open(os.path.join(D, 'progress.txt')) if not l.startswith('#')]
for p0, p1, npr, mx, *_ in rows:
    if p0 in done: continue
    path = os.path.join(D, 'chunks', f'c_{p0}.txt')
    data = open(path, 'rb').read()
    sha = hashlib.sha256(data).hexdigest()
    lines = data.decode().rstrip('\n').split('\n')
    ok = lines[-1] == f"# done [{p0}, {p1}) primes={npr} max_mx={mx}"
    body = [l for l in lines if not l.startswith('#')]
    ok = ok and len(body) == int(npr)
    lo = random.Random(int(p0)).randrange(int(p0), max(int(p0) + 1, int(p1) - 600))
    hi = min(lo + 600, int(p1))
    mine = [l for l in body if lo <= int(l.split()[0]) < hi]
    theirs = ref(lo, hi)
    ok = ok and mine == theirs
    close = [l for l in body if int(l.split()[1]) >= 15]
    for l in close:
        p = int(l.split()[0])
        ok = ok and ref(p, p + 1) == [l]
    with open(OUT, 'a') as f:
        f.write(f"{p0} {p1} {npr} {mx} {sha} {lo} {hi} {len(mine)} {len(close)} {'PASS' if ok else 'FAIL'}\n")
    print(p0, 'PASS' if ok else 'FAIL', len(mine), len(close), flush=True)
