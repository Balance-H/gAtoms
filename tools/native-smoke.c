#include "decom_h.h"
#include <assert.h>
#include <stdio.h>
/* Builds the same undirected path (circular == 0) or cycle (circular != 0) on n
 * vertices that igraph_ring() would, but from igraph_empty() plus
 * igraph_add_edge(). Both are part of the object set the package actually
 * compiles (igraph_ring is not), so this harness links against exactly the
 * same igraph subset as the package and cannot silently depend on extra code. */
static igraph_error_t test_ring(igraph_t *g, igraph_integer_t n, igraph_bool_t circular) {
    igraph_integer_t i;
    igraph_error_t e = igraph_empty(g, n, IGRAPH_UNDIRECTED);
    if (e != IGRAPH_SUCCESS) return e;
    for (i = 0; i + 1 < n; ++i) {
        e = igraph_add_edge(g, i, i + 1);
        if (e != IGRAPH_SUCCESS) return e;
    }
    if (circular) e = igraph_add_edge(g, n - 1, 0);
    return e;
}
int main(void) {
    igraph_error_handler_t *old=igraph_set_error_handler(igraph_error_handler_ignore);
    for(int i=0;i<200;++i) {
        igraph_t g,t; igraph_vector_ptr_t a,s,o; igraph_vector_int_t seed,h;
        assert(!test_ring(&g,8,i%2)); assert(!igraph_empty(&t,0,0));
        assert(!igraph_vector_ptr_init(&a,0)); assert(!igraph_vector_ptr_init(&s,0)); assert(!igraph_vector_ptr_init(&o,0));
        assert(!igraph_vector_int_init_int(&seed,2,0,2)); assert(!igraph_vector_int_init(&h,20));
        assert(!decompose_atoms(&g,&a,&s,&t,&o));
        assert(igraph_vector_ptr_size(&o)==2*igraph_vector_ptr_size(&a)-1);
        assert(!get_minimal_collapsible(&g,&seed,&h));
        VECTOR(seed)[0]=-1; assert(get_minimal_collapsible(&g,&seed,&h)==IGRAPH_EINVVID);
        assert(build_prime_ordering(&a,&g,&o)==IGRAPH_EINVAL);
        igraph_vector_int_destroy(&seed); igraph_vector_int_destroy(&h);
        free_vertex_list(&a); free_vertex_list(&s); free_vertex_list(&o); igraph_destroy(&t); igraph_destroy(&g);
    }
    igraph_set_error_handler(old); puts("native smoke: 200 repetitions passed"); return 0;
}
