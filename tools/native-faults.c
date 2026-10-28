#include "decom_h.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
void *__real_malloc(size_t);
void __real_free(void *);
/* Test-only allocator state, never part of the package library. */
static int attempt,fail_at,live,enabled;
static void *owned[4096];
void *__wrap_malloc(size_t n) {
    if (enabled && ++attempt==fail_at) return NULL;
    void *p=__real_malloc(n);
    if(p && enabled) { assert(live<4096); owned[live++]=p; }
    return p;
}
void __wrap_free(void *p) {
    for(int i=0;i<live;++i) if(owned[i]==p) { owned[i]=owned[--live]; break; }
    __real_free(p);
}
int main(void) {
    igraph_set_error_handler(igraph_error_handler_ignore);
    for(int mode=0;mode<2;++mode) {
    int max_attempt=0;
    for(int f=0;f<=max_attempt;++f) {
        igraph_t g,t; igraph_vector_ptr_t a,s,o;
        assert(!igraph_ring(&g,8,0,0,mode)); assert(!igraph_empty(&t,0,0));
        assert(!igraph_vector_ptr_init(&a,0)); assert(!igraph_vector_ptr_init(&s,0)); assert(!igraph_vector_ptr_init(&o,0));
        attempt=0;fail_at=f;enabled=1;
        igraph_error_t ret=decompose_atoms(&g,&a,&s,&t,&o);
        enabled=0;
        if(f==0) { assert(!ret); max_attempt=attempt; }
        else assert(ret==IGRAPH_ENOMEM);
        free_vertex_list(&a);free_vertex_list(&s);free_vertex_list(&o);igraph_destroy(&t);igraph_destroy(&g);
        if(live) { fprintf(stderr,"failure position %d: live=%d attempts=%d result=%d\n", f, live, attempt, (int)ret); return 1; }
    }
    printf("native allocation failures (circular=%d): %d checked failure positions, no owned allocations leaked\n",mode,max_attempt);
    }
    return 0;
}
