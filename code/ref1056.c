// Reference (independent, no shortcuts): for each prime p in [P0,P1): plain m! mod p with %, direct count array.
// Same output format as f1056: "p mx nmx [m ...]" with every m < p-1 having m! = p-1 (mod p).
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
int isprime(uint64_t n){ if(n<2) return 0; for(uint64_t d=2; d*d<=n; d++) if(n%d==0) return 0; return 1; }
int main(int argc, char **argv){
  uint64_t P0=strtoull(argv[1],0,10), P1=strtoull(argv[2],0,10);
  uint32_t *cnt = calloc(P1, 4);
  for(uint64_t p=(P0<5?5:P0); p<P1; p++){ if(!isprime(p)) continue;
    uint64_t f=1; uint32_t neg[1000]; int nn=0;
    for(uint64_t m=1; m<p; m++){ f=f*m%p; cnt[f]++; if(f==p-1 && m<p-1 && nn<1000) neg[nn++]=m; }
    uint32_t mx=0, nmx=0; for(uint64_t v=0; v<p; v++){ if(cnt[v]>mx){mx=cnt[v]; nmx=1;} else if(cnt[v]==mx) nmx++; cnt[v]=0; }
    printf("%llu %u %u", (unsigned long long)p, mx, nmx); for(int i=0;i<nn && i<16;i++) printf(" %u", neg[i]); if(nn>16) printf(" +%d", nn-16); printf("\n");
  }
  return 0;
}
