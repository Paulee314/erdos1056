# Erdős problem 1056: new terms a(15) and a(16) of OEIS A060427

**Problem** ([erdosproblems.com/1056](https://www.erdosproblems.com/1056), Guy's *Unsolved Problems in Number
Theory* A15): Let k ≥ 2. Does there exist a prime p and consecutive intervals I₁, …, I_k such that
∏_{n ∈ I_i} n ≡ 1 (mod p) for all 1 ≤ i ≤ k?

The problem page notes that it is open and cannot be settled by a finite computation. This repository adds data,
not a proof: the smallest such prime for k = 15 and k = 16, from an exhaustive search with independent checks.

[OEIS A060427](https://oeis.org/A060427) lists the least such prime for each k. Its terms run to
a(14) = 10,428,007. k adjacent intervals with product 1 mod p are the same thing as k + 1 equal values among
m! mod p (m = 1, …, p − 1): if m₀! ≡ m₁! ≡ … ≡ m_k! with m₀ < m₁ < … < m_k, the blocks (m_{i−1}, m_i] each multiply
to 1. So a(k) is the least prime for which some value occurs k + 1 times among 1!, 2!, …, (p − 1)! mod p. This
convention reproduces every listed term a(2), …, a(14).

## Results

| k | a(k) | repeated value of m! mod p | found |
|---|---|---|---|
| 15 | **34,614,961** | 33,919,633 (16 times) | 26 Oct 2026, cloud search |
| 16 | **169,916,011** | 140,320,727 (17 times) | 28 Oct 2026, Mac Studio search |

a(15): the value 33,919,633 occurs at m = 45566, 291444, 1119644, 3438999, 3614660, 4448849, 5767103, 11242336,
11684541, 12800791, 18784810, 19897411, 19969551, 23762668, 33505794, 34394461, so the 15 blocks between them each
multiply to 1 mod 34,614,961. No other value reaches 16 for this prime.

a(16): the value 140,320,727 occurs at m = 1309627, 18412281, 44736402, 44831451, 46789287, 57243054, 59449242,
86843187, 108076990, 108853377, 111309733, 121283195, 129898421, 144445796, 145598172, 151727983, 157387059, so the
16 blocks between them each multiply to 1 mod 169,916,011. No other value reaches 17 for this prime.

**Coverage.** Every prime in [5, 34,677,737) was scanned in the cloud (25–26 Oct 2026) and every prime in
[34,677,737, 170,201,504) on a Mac Studio M3 Ultra (27–28 Oct 2026). The block ranges chain with no gaps
(`summarize.py` asserts this), and no prime below a(15) has 16 equal values, and none below a(16) has 17.

## How the search works

For each prime p the scan computes every m! mod p and finds the largest number of m sharing one value. That is
about p multiplications per prime, so the whole search is roughly Σ p ≈ p²/(2 ln p) steps. The fast programs use
Montgomery multiplication with four primes interleaved per thread, and count repeats with a two-pass bucket sort
that stays in cache. The output has one line per prime: `p mx nmx [m ...]`, where mx is the largest repeat count,
nmx how many values reach it, and the m's are those with m! ≡ −1 (mod p) below p − 2 (kept as data for problems
1072 and 1074).

## Checks

1. **Independent reference.** `ref1056.c` is a separate, plain program (ordinary `%`, a direct count array, no
   shortcuts) with the same output format.
2. **Every finished block is certified** (`certify.py`): the SHA-256 of its output, its closing line against the
   progress log, a randomly placed window of about 600 numbers re-run with `ref1056` and compared byte for byte,
   and every prime in the block with 15 or more equal values re-run with `ref1056`. Cloud: 90 of 90 blocks pass.
   Mac: 250 of 250 pass.
3. **The two builds agree.** The cloud program (`f1056c.c`, x86) and the portable Mac program (`f1056p.c`) give
   byte-identical output on [5, 300000), [10⁷, 10⁷ + 400) and the whole block containing a(15).
4. **Each new term checked a third way.** `verify_term.py` (Python and NumPy, written separately) recomputes all
   m! mod p for the single prime, lists the m's, checks that each block multiplies to 1 mod p, and runs Wilson's
   theorem as a sanity check.
5. **Known terms.** The scan reproduces OEIS a(2) through a(14) exactly.

## A heuristic for where the terms fall

If the values m! mod p behaved like p − 1 independent random residues, the number of residues hit k + 1 or more
times would be close to Poisson, giving a per-prime chance of about (p − 1)·P(Poisson(1) ≥ k + 1). `model.py`
turns this into a predicted median and 90% range for each a(k). Of the 14 terms a(2), …, a(15), 13 fall inside
their 90% range; a(8) = 599 falls just below it. It is a heuristic only, not a proof of anything.

## Contents

| Path | What |
|---|---|
| `code/f1056c.c` | Cloud scan (x86). Build: `cc -O3 -march=native -o f1056c f1056c.c -lpthread -lm` |
| `code/f1056p.c` | Portable scan used on the Mac. Build: `cc -O3 -o f1056p f1056p.c -lpthread -lm` |
| `code/ref1056.c` | Independent plain reference. Build: `cc -O2 -o ref1056 ref1056.c` |
| `code/run1056.py` | Runs the scan in blocks and appends to `progress.txt`; stops when the target is found. |
| `code/certify.py` | The per-block certificate described above. |
| `code/summarize.py` | Checks that the blocks chain with no gaps and prints the least prime for each k. |
| `code/verify_term.py` | Separate single-prime check of a term. |
| `code/model.py` | The random-factorials heuristic. |
| `results/a15/`, `results/a16/` | The two new terms, with every repeated m and the checks run on them. |
| `results/cloud/`, `results/mac/` | Progress logs (block ranges, prime counts, timings), certificates, and every prime with 15 or more equal values. |
| `results/validation/` | Outputs used to confirm the two builds agree. |
| `data/` | Every per-prime output line, compressed (about 90 MB). |

The full per-prime outputs (one line for each of the 9,513,810 primes in [5, 170,201,504)) are in `data/`:
`cloud_chunks.tgz` (blocks below 34,677,737) and `mac_chunks_1.tgz`, `mac_chunks_2.tgz` (Mac blocks starting below and
above 100,000,000). `results/*/certify.txt` lists the SHA-256 of every block file.

**Re-check a term:** `cc -O2 -o ref1056 code/ref1056.c && ./ref1056 169916011 169916012` prints
`169916011 17 1 84958005` in about 2 seconds; `python3 code/verify_term.py 169916011 16` lists the 17 m's and checks
the blocks (about a minute; 2 GB of memory). The same works for a(15) with `34614961 34614962` and `34614961 15`.

## Credits

The search programs, checks and this write-up were written by Claude (Anthropic) at the direction of Paul Lee, who
ran the Mac search. Related: problem page [erdosproblems.com/1056](https://www.erdosproblems.com/1056);
Makowski (1983) gave the k = 3 example mod 17; Erdős observed the k = 2 example mod 11.

License: Apache 2.0.
