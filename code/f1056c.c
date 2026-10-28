// Erdos #1056 / OEIS A060427 scan.
// For each prime p in [P0, P1): the largest number of m in [1, p-1] sharing one value of m! mod p ("mx").
// a(k) = least prime with mx >= k+1 (k adjacent intervals with product 1 mod p <=> k+1 equal factorials).
// Also records every m in [1, p-2) with m! = -1 (mod p) (data for #1072/#1074; (p-1)! = -1 always).
// Method: K primes advance in lockstep (independent Montgomery chains, for instruction-level parallelism).
// Pass 1 writes each m! (in Montgomery form: a bijection of residues, so equal counts are unchanged) into
// buckets by its high bits; pass 2 counts each bucket in a cache-resident table.
// Usage: f1056 P0 P1 threads outfile
// outfile lines: "p mx nmx [m ...]"  (nmx = how many residues reach mx; m = all m < p-2 with m! = -1)
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <immintrin.h>
#define BATCH 2048
#define SLOT 32   /* u16 per staging slot = one 64-byte line */
typedef uint64_t u64; typedef uint32_t u32; typedef uint16_t u16; typedef uint8_t u8;
#define K 4
#define BSH 16
#define MAXNEG 16

static u32 *primes; static long nprimes; static long next_group = 0;
static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
static FILE *out;
typedef struct { u32 p, mx, nmx, nneg; u32 neg[MAXNEG]; } res_t;
static res_t *results;

static inline u32 montmul(u32 a, u32 b, u32 p, u32 pinv) {   // a*b*2^-32 mod p
    u64 t = (u64)a * b; u32 q = (u32)t * pinv; u64 u = (t + (u64)q * p) >> 32;
    return (u32)(u >= p ? u - p : u);
}
static u32 neg_inv32(u32 p) { u32 x = p; for (int i = 0; i < 5; i++) x *= 2 - p * x; return -x; }

typedef struct { u16 *buf; u16 **ptr; u16 **end; u32 nb, cap; size_t bufsz; u16 *stage; u8 *sfill; u32 stagenb; } lane_t;

static void lane_setup(lane_t *L, u32 p) {
    u32 nb = (p >> BSH) + 1;
    u32 cap = (u32)((1u << BSH) * ((double)(p - 1) / p) + 12 * sqrt((double)(1u << BSH)) + 64);
    cap = (cap + SLOT - 1) / SLOT * SLOT;                 // whole 64-byte lines per bucket
    size_t need = (size_t)nb * cap;
    if (need > L->bufsz) { free(L->buf); L->buf = aligned_alloc(64, need * sizeof(u16)); L->bufsz = need; if (!L->buf) { perror("malloc"); exit(1); } }
    if (nb > L->nb) { free(L->ptr); free(L->end); L->ptr = malloc(nb * sizeof(u16*)); L->end = malloc(nb * sizeof(u16*)); }
    if (nb > L->stagenb) { free(L->stage); free(L->sfill); L->stage = aligned_alloc(64, (size_t)nb * SLOT * sizeof(u16)); L->sfill = malloc(nb); L->stagenb = nb; }
    L->nb = nb; L->cap = cap; memset(L->sfill, 0, nb);
    for (u32 b = 0; b < nb; b++) { L->ptr[b] = L->buf + (size_t)b * cap; L->end[b] = L->ptr[b] + cap; }
}

static __thread int overflow_flag;
static inline void put(lane_t *L, u32 v) {            // stage in a 64-byte slot; full slots go out with a non-temporal store
    u32 b = v >> BSH; u32 f = L->sfill[b]; u16 *s = L->stage + (size_t)b * SLOT;
    s[f] = (u16)v;
    if (++f == SLOT) {
        u16 *q = L->ptr[b];
        if (__builtin_expect(q + SLOT > L->end[b], 0)) { overflow_flag = 1; f = 0; }
        else { _mm512_stream_si512((__m512i *)q, _mm512_load_si512((const void *)s)); L->ptr[b] = q + SLOT; f = 0; }
    }
    L->sfill[b] = (u8)f;
}
static void flush_lane(lane_t *L) {
    for (u32 b = 0; b < L->nb; b++) {
        u32 f = L->sfill[b]; u16 *s = L->stage + (size_t)b * SLOT, *q = L->ptr[b];
        if (q + f > L->end[b]) { overflow_flag = 1; continue; }
        memcpy(q, s, f * sizeof(u16)); L->ptr[b] = q + f; L->sfill[b] = 0;
    }
    _mm_sfence();
}

// pass 2: count each bucket in a 64K-entry table (L1/L2), then take the max with AVX-512 and clear with memset.
static void count_lane(lane_t *L, u8 *cnt, u32 *mx_out, u32 *nmx_out) {
    u32 mx = 0, nmx = 0; const u32 W = 1u << BSH;
    for (u32 b = 0; b < L->nb; b++) {
        u16 *s = L->buf + (size_t)b * L->cap, *e = L->ptr[b];
        for (u16 *q = s; q < e; q++) cnt[*q]++;
        __m512i vm = _mm512_setzero_si512();
        for (u32 i = 0; i < W; i += 64) vm = _mm512_max_epu8(vm, _mm512_load_si512((const void *)(cnt + i)));
        u8 tmp[64] __attribute__((aligned(64))); _mm512_store_si512((void *)tmp, vm);
        u32 bmx = 0; for (int i = 0; i < 64; i++) if (tmp[i] > bmx) bmx = tmp[i];
        if (bmx >= mx && bmx > 0) {
            u32 n = 0; for (u32 i = 0; i < W; i++) n += (cnt[i] == bmx);
            if (bmx > mx) { mx = bmx; nmx = n; } else nmx += n;
        }
        memset(cnt, 0, W);
    }
    *mx_out = mx; *nmx_out = nmx;
}

static void process_group(long g0, long gn, lane_t *lanes, u8 *cnt) {
    u32 p[K], pinv[K], F[K], M[K], R1[K], NEG[K]; int nl = (int)gn;
    for (int j = 0; j < nl; j++) {
        p[j] = primes[g0 + j]; pinv[j] = neg_inv32(p[j]);
        R1[j] = (u32)((1ull << 32) % p[j]);              // Montgomery form of 1
        F[j] = R1[j]; M[j] = 0;                           // F = Mont(0!) ; M = Mont(0), incremented to Mont(1) first
        NEG[j] = p[j] - R1[j];                            // Montgomery form of -1
        lane_setup(&lanes[j], p[j]);
        results[g0 + j].p = p[j]; results[g0 + j].nneg = 0;
    }
    u32 pmin = p[0]; for (int j = 1; j < nl; j++) if (p[j] < pmin) pmin = p[j];
    overflow_flag = 0;
    u32 m = 1;
    if (nl == K) {
        static __thread u32 vb[K][BATCH];
        while (m + BATCH <= pmin - 1) {                    // full batches inside m = 1 .. pmin-2
            for (int i = 0; i < BATCH; i++, m++) {
                for (int j = 0; j < K; j++) {
                    u32 t = M[j] + R1[j]; M[j] = t >= p[j] ? t - p[j] : t;
                    F[j] = montmul(F[j], M[j], p[j], pinv[j]); vb[j][i] = F[j];
                    if (__builtin_expect(F[j] == NEG[j], 0)) { res_t *r = &results[g0 + j]; if (r->nneg < MAXNEG) r->neg[r->nneg] = m; r->nneg++; }
                }
            }
            for (int j = 0; j < K; j++) { lane_t *Lj = &lanes[j]; for (int i = 0; i < BATCH; i++) put(Lj, vb[j][i]); }
        }
    }
    for (int j = 0; j < nl; j++) {                         // finish each lane individually up to m = p-1
        for (u32 mm = m; mm < p[j]; mm++) {
            u32 t = M[j] + R1[j]; M[j] = t >= p[j] ? t - p[j] : t;
            F[j] = montmul(F[j], M[j], p[j], pinv[j]);
            put(&lanes[j], F[j]);
            if (F[j] == NEG[j] && mm < p[j] - 1) { res_t *r = &results[g0 + j]; if (r->nneg < MAXNEG) r->neg[r->nneg] = mm; r->nneg++; }
        }
        if (F[j] != NEG[j]) { fprintf(stderr, "Wilson check failed at p=%u\n", p[j]); exit(2); }
        flush_lane(&lanes[j]);
    }
    if (overflow_flag) { fprintf(stderr, "bucket overflow in group at p=%u\n", p[0]); exit(3); }
    for (int j = 0; j < nl; j++) { u32 mx, nmx; count_lane(&lanes[j], cnt, &mx, &nmx); results[g0 + j].mx = mx; results[g0 + j].nmx = nmx; }
}

static void *worker(void *arg) {
    (void)arg;
    lane_t lanes[K]; memset(lanes, 0, sizeof lanes);
    u8 *cnt = aligned_alloc(64, 1u << BSH); memset(cnt, 0, 1u << BSH);
    for (;;) {
        pthread_mutex_lock(&mu); long g = next_group; next_group += K; pthread_mutex_unlock(&mu);
        if (g >= nprimes) break;
        long gn = nprimes - g < K ? nprimes - g : K;
        process_group(g, gn, lanes, cnt);
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 5) { fprintf(stderr, "usage: f1056 P0 P1 threads outfile\n"); return 1; }
    u64 P0 = strtoull(argv[1], 0, 10), P1 = strtoull(argv[2], 0, 10); int T = atoi(argv[3]);
    if (P0 < 5) P0 = 5;
    if (P1 >= (1ull << 31)) { fprintf(stderr, "P1 too large\n"); return 1; }
    // segmented sieve of [P0, P1)
    u64 L = P1 - P0; u8 *comp = calloc(L, 1); u64 r = (u64)sqrt((double)P1) + 2;
    u8 *small = calloc(r + 1, 1);
    for (u64 i = 2; i <= r; i++) if (!small[i]) {
        for (u64 j = i * i; j <= r; j += i) small[j] = 1;
        u64 s = (P0 + i - 1) / i * i; if (s < i * i) s = i * i;
        for (u64 j = s; j < P1; j += i) comp[j - P0] = 1;
    }
    nprimes = 0; for (u64 i = 0; i < L; i++) if (!comp[i]) nprimes++;
    primes = malloc(nprimes * sizeof(u32)); long k = 0;
    for (u64 i = 0; i < L; i++) if (!comp[i]) primes[k++] = (u32)(P0 + i);
    free(comp); free(small);
    results = calloc(nprimes, sizeof(res_t));
    pthread_t th[64]; for (int t = 0; t < T; t++) pthread_create(&th[t], 0, worker, 0);
    for (int t = 0; t < T; t++) pthread_join(th[t], 0);
    out = fopen(argv[4], "w"); if (!out) { perror("open"); return 1; }
    u32 best = 0;
    for (long i = 0; i < nprimes; i++) {
        res_t *q = &results[i];
        fprintf(out, "%u %u %u", q->p, q->mx, q->nmx);
        for (u32 j = 0; j < q->nneg && j < MAXNEG; j++) fprintf(out, " %u", q->neg[j]);
        if (q->nneg > MAXNEG) fprintf(out, " +%u", q->nneg - MAXNEG);
        fprintf(out, "\n");
        if (q->mx > best) best = q->mx;
    }
    fprintf(out, "# done [%llu, %llu) primes=%ld max_mx=%u\n", (unsigned long long)P0, (unsigned long long)P1, nprimes, best);
    fclose(out);
    return 0;
}
