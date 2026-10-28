#include "decom_h.h"
#include <assert.h>
#include <stdio.h>
int main(void) {
    igraph_error_handler_t *old=igraph_set_error_handler(igraph_error_handler_ignore);
    for(int i=0;i<200;++i) {
        igraph_t g,t; igraph_vector_ptr_t a,s,o; igraph_vector_int_t seed,h;
        assert(!igraph_ring(&g,8,0,0,i%2)); assert(!igraph_empty(&t,0,0));
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
