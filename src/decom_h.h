#ifndef DECOM_H_H
#define DECOM_H_H
#include <igraph/igraph.h>
/* gAtoms is built against the igraph C library bundled in src/igraph, whose
 * vendored copy is the 1.0.1 release. The native core is written against the
 * 1.0 series API and its index conventions (see the subgraph-mapping note in
 * decom_h.c), so the guard below accepts the 1.0.x series from 1.0.1 upwards:
 * a later 1.0.z patch release keeps the same API and mapping semantics, while a
 * new major or minor series does not. 1.0.1 is the version the test suite has
 * actually been run against. */
#if IGRAPH_VERSION_MAJOR != 1 || IGRAPH_VERSION_MINOR != 0 || IGRAPH_VERSION_PATCH < 1
#error "gAtoms requires the igraph C 1.0.x API, version 1.0.1 or newer"
#endif
/* Zero-based IDs. Initialized outputs; returned list vectors are deep copies
 * owned by caller. Partial outputs on error are valid and must be freed. */
void free_vertex_list(igraph_vector_ptr_t *);
igraph_error_t components_forbidden(const igraph_t *, igraph_vector_ptr_t *, igraph_vector_ptr_t *, const igraph_vector_int_t *);
igraph_error_t close_separator(const igraph_t *, igraph_integer_t, const igraph_vector_int_t *, igraph_vector_int_t *);
igraph_error_t get_minimal_collapsible(const igraph_t *, const igraph_vector_int_t *, igraph_vector_int_t *);
/* Connected simple undirected graph; outputs and undirected tree must be empty. */
igraph_error_t decompose_atoms(const igraph_t *, igraph_vector_ptr_t *, igraph_vector_ptr_t *, igraph_t *, igraph_vector_ptr_t *);
igraph_error_t build_prime_ordering(const igraph_vector_ptr_t *, const igraph_t *, igraph_vector_ptr_t *);
#endif
