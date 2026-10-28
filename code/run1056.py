#!/usr/bin/env python3
"""Drive the #1056 scan in chunks of about W steps (sum of p over the chunk's primes), resuming from progress.txt.
Usage: run1056.py [minutes_limit] [threads] [binary] [start_P0] [steps_per_chunk]
  cloud: run1056.py 105            (2 threads, f1056c, resumes from progress.txt)
  Mac:   run1056.py 100000 28 ./f1056p <cloud frontier> 3e12
Stops before the time limit, or as soon as a prime with mx >= TARGET_MX (17 for a(16)) is found."""
import math, os, subprocess, sys, time, hashlib, datetime
HERE = os.path.dirname(os.path.abspath(__file__))
A = sys.argv[1:]
THREADS = A[1] if len(A) > 1 else "2"
BIN = os.path.abspath(A[2]) if len(A) > 2 else f"{HERE}/f1056c"
START = int(float(A[3])) if len(A) > 3 else 5
W = float(A[4]) if len(A) > 4 else 4e11   # steps per chunk (about 17 min at 2.6 ns/step on 2 cloud threads)
PROG, CH = f"{HERE}/progress.txt", f"{HERE}/chunks"
TARGET_MX = int(os.environ.get("TARGET_MX", "17"))   # a(16) needs 17 equal factorial values
limit = float(sys.argv[1]) * 60 if len(sys.argv) > 1 else 100 * 60
os.makedirs(CH, exist_ok=True)
if not os.path.exists(PROG):
    with open(PROG, "w") as f:
        f.write(f"# {os.path.basename(BIN)} sha256 {hashlib.sha256(open(BIN, 'rb').read()).hexdigest()} threads {THREADS}\n# P0 P1 nprimes max_mx seconds finished_utc\n")
t_start, rate = time.time(), None
while True:
    rows = [l.split() for l in open(PROG) if not l.startswith("#")]
    P0 = int(rows[-1][1]) if rows else START
    if any(int(r[3]) >= TARGET_MX for r in rows):
        print("target reached; stopping"); break
    P1 = int(math.sqrt(P0 * P0 + 2 * W * math.log(max(P0, 100)))) + 1
    est = rate * W if rate else 0
    if time.time() - t_start + est * 1.15 > limit:
        print(f"time limit; next chunk would start at {P0}"); break
    out = f"{CH}/c_{P0}.txt"
    t0 = time.time()
    subprocess.run([BIN, str(P0), str(P1), THREADS, out], check=True)
    dt = time.time() - t0
    last = open(out).read().rstrip().split("\n")[-1]
    assert last.startswith(f"# done [{P0}, {P1})"), last
    nprimes = int(last.split("primes=")[1].split()[0]); mxx = int(last.split("max_mx=")[1])
    steps = sum(int(l.split()[0]) - 1 for l in open(out) if not l.startswith("#"))
    rate = dt / steps
    with open(PROG, "a") as f:
        f.write(f"{P0} {P1} {nprimes} {mxx} {dt:.0f} {datetime.datetime.now(datetime.timezone.utc):%Y-%m-%dT%H:%M:%SZ}\n")
    print(f"[{P0}, {P1}) primes={nprimes} max_mx={mxx} {dt:.0f}s {rate*1e9:.2f} ns/step", flush=True)
