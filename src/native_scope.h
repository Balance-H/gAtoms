#ifndef GATOMS_NATIVE_SCOPE_H
#define GATOMS_NATIVE_SCOPE_H
#include "decom_h.h"
#include <stdlib.h>

/* Per-call ownership stack, never global. Each scope has a bounded number of
 * local objects; loops use/reset child scopes. Registration failure destroys
 * the object immediately. Outputs are deep copies owned by the caller. */
typedef void (*ga_destructor)(void *);
typedef struct { void *p[32]; ga_destructor destroy[32]; int n; } ga_scope;
static void ga_clear(ga_scope *s) {
    while (s->n) { --s->n; s->destroy[s->n](s->p[s->n]); }
}
static igraph_error_t ga_own(ga_scope *s, void *p, ga_destructor d) {
    if (s->n == 32) { d(p); return IGRAPH_EOVERFLOW; }
    s->p[s->n] = p; s->destroy[s->n++] = d; return IGRAPH_SUCCESS;
}
static void ga_vec_free(void *p) { igraph_vector_int_destroy(p); free(p); }
static void ga_graph_free(void *p) { igraph_destroy(p); free(p); }
static void ga_list_free(void *p) { free_vertex_list(p); free(p); }
static igraph_error_t ga_vec(ga_scope *s, igraph_vector_int_t **v, igraph_integer_t n) {
    igraph_error_t e;
    *v = malloc(sizeof(**v));
    if (!*v) return IGRAPH_ENOMEM;
    e = igraph_vector_int_init(*v, n);
    if (e) { free(*v); return e; }
    return ga_own(s, *v, ga_vec_free);
}
static igraph_error_t ga_list(ga_scope *s, igraph_vector_ptr_t **v) {
    igraph_error_t e;
    *v = malloc(sizeof(**v));
    if (!*v) return IGRAPH_ENOMEM;
    e = igraph_vector_ptr_init(*v, 0);
    if (e) { free(*v); return e; }
    return ga_own(s, *v, ga_list_free);
}
static igraph_error_t ga_subgraph(ga_scope *s, const igraph_t *g,
    const igraph_vector_int_t *nodes, igraph_t **sub,
    igraph_vector_int_t *map, igraph_vector_int_t *invmap) {
    igraph_error_t e;
    *sub = malloc(sizeof(**sub));
    if (!*sub) return IGRAPH_ENOMEM;
    /* igraph_vss_vector is a non-owning selector; no allocation to destroy. */
    e = igraph_induced_subgraph_map(g, *sub, igraph_vss_vector(nodes),
                                  IGRAPH_SUBGRAPH_AUTO, map, invmap);
    if (e) { free(*sub); return e; }
    return ga_own(s, *sub, ga_graph_free);
}
#define TRY(expr) do { ret = (expr); if (ret != IGRAPH_SUCCESS) goto cleanup; } while (0)
#define VEC(s, v, n) igraph_vector_int_t *v; TRY(ga_vec(&(s), &(v), (n)))
#define LIST(s, v) igraph_vector_ptr_t *v; TRY(ga_list(&(s), &(v)))
#endif
