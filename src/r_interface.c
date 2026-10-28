#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <limits.h>
#include "decom_h.h"

/* R errors (including allocation failures while converting results) unwind
 * through cleanup. No persistent/global state, no duplicate fallback backend. */
typedef struct {
    SEXP edges,seeds,atom_input;
    int n,mode;
    igraph_t graph,tree;
    igraph_vector_int_t ev,sv,h;
    igraph_vector_ptr_t atoms,separators,ordering;
    int graph_ok,tree_ok,ev_ok,sv_ok,h_ok,a_ok,s_ok,o_ok,handler_ok;
    igraph_error_handler_t *old_handler;
} call_data;
#define RCHECK(expr) do { igraph_error_t rc=(expr); if(rc) Rf_error("gAtoms native %s: %s (graph must be simple/undirected; atom tree must be connected with k vertices and k-1 edges)", #expr,igraph_strerror(rc)); } while(0)
static int whole(SEXP x,R_xlen_t i,const char *name) {
    double d;
    if(TYPEOF(x)==INTSXP) { int z=INTEGER(x)[i]; if(z!=NA_INTEGER) return z; }
    else if(TYPEOF(x)==REALSXP) {
        d=REAL(x)[i]; if(R_FINITE(d) && d>=0 && d<=INT_MAX && d==(int)d) return (int)d;
    }
    Rf_error("%s must contain finite whole-number IDs",name); return 0;
}
static void validate_edges(SEXP e,int n) {
    SEXP dim=Rf_getAttrib(e,R_DimSymbol);
    if(!Rf_isMatrix(e) || XLENGTH(dim)!=2 || INTEGER(dim)[1]!=2 ||
       (TYPEOF(e)!=INTSXP && TYPEOF(e)!=REALSXP)) Rf_error("edge_mat must be a numeric matrix with two columns");
    for(R_xlen_t i=0;i<XLENGTH(e);++i) {
        int z=whole(e,i,"edge_mat"); if(z<0 || z>=n) Rf_error("edge endpoints must lie in [0, n_vertices - 1]");
    }
}
static void cleanup(void *data,Rboolean jump) {
    call_data *d=data; (void)jump;
    if(d->o_ok) free_vertex_list(&d->ordering);
    if(d->s_ok) free_vertex_list(&d->separators);
    if(d->a_ok) free_vertex_list(&d->atoms);
    if(d->h_ok) igraph_vector_int_destroy(&d->h);
    if(d->sv_ok) igraph_vector_int_destroy(&d->sv);
    if(d->ev_ok) igraph_vector_int_destroy(&d->ev);
    if(d->tree_ok) igraph_destroy(&d->tree);
    if(d->graph_ok) igraph_destroy(&d->graph);
    if(d->handler_ok) igraph_set_error_handler(d->old_handler);
}
static SEXP to_vector(const igraph_vector_int_t *v) {
    SEXP out=PROTECT(Rf_allocVector(INTSXP,igraph_vector_int_size(v)));
    for(igraph_integer_t i=0;i<igraph_vector_int_size(v);++i) INTEGER(out)[i]=(int)VECTOR(*v)[i]+1;
    UNPROTECT(1); return out;
}
static SEXP to_list(const igraph_vector_ptr_t *v) {
    SEXP out=PROTECT(Rf_allocVector(VECSXP,igraph_vector_ptr_size(v)));
    for(igraph_integer_t i=0;i<igraph_vector_ptr_size(v);++i) SET_VECTOR_ELT(out,i,to_vector(VECTOR(*v)[i]));
    UNPROTECT(1); return out;
}
static SEXP run(void *data) {
    call_data *d=data;
    d->old_handler=igraph_set_error_handler(igraph_error_handler_ignore); d->handler_ok=1;
    RCHECK(igraph_vector_int_init(&d->ev,XLENGTH(d->edges))); d->ev_ok=1;
    R_xlen_t m=XLENGTH(d->edges)/2;
    for(R_xlen_t i=0;i<m;++i) {
        VECTOR(d->ev)[2*i]=whole(d->edges,i,"edges"); VECTOR(d->ev)[2*i+1]=whole(d->edges,i+m,"edges");
    }
    RCHECK(igraph_create(&d->graph,&d->ev,d->n,0)); d->graph_ok=1;
    RCHECK(igraph_vector_int_init(&d->h,0)); d->h_ok=1;
    if(d->mode==1) {
        RCHECK(igraph_vector_int_init(&d->sv,XLENGTH(d->seeds))); d->sv_ok=1;
        for(R_xlen_t i=0;i<XLENGTH(d->seeds);++i) VECTOR(d->sv)[i]=whole(d->seeds,i,"r_nodes");
        RCHECK(get_minimal_collapsible(&d->graph,&d->sv,&d->h));
        return to_vector(&d->h);
    }
    RCHECK(igraph_vector_ptr_init(&d->atoms,0)); d->a_ok=1;
    RCHECK(igraph_vector_ptr_init(&d->separators,0)); d->s_ok=1;
    RCHECK(igraph_vector_ptr_init(&d->ordering,0)); d->o_ok=1;
    if(d->mode==2) {
        for(R_xlen_t i=0;i<XLENGTH(d->atom_input);++i) {
            SEXP a=VECTOR_ELT(d->atom_input,i);
            igraph_vector_int_t *v=malloc(sizeof(*v)); if(!v) Rf_error("memory allocation failed");
            igraph_error_t rc=igraph_vector_int_init(v,XLENGTH(a));
            if(rc) { free(v); RCHECK(rc); }
            for(R_xlen_t j=0;j<XLENGTH(a);++j) VECTOR(*v)[j]=whole(a,j,"atoms");
            rc=igraph_vector_ptr_push_back(&d->atoms,v);
            if(rc) { igraph_vector_int_destroy(v); free(v); RCHECK(rc); }
        }
        RCHECK(build_prime_ordering(&d->atoms,&d->graph,&d->ordering));
        return to_list(&d->ordering);
    }
    RCHECK(igraph_empty(&d->tree,0,0)); d->tree_ok=1;
    RCHECK(decompose_atoms(&d->graph,&d->atoms,&d->separators,&d->tree,&d->ordering));
    SEXP out=PROTECT(Rf_allocVector(VECSXP,4)),names=PROTECT(Rf_allocVector(STRSXP,4));
    const char *labels[]={"atoms","separators","tree","prime_ordering"};
    for(int i=0;i<4;++i) SET_STRING_ELT(names,i,Rf_mkChar(labels[i]));
    Rf_setAttrib(out,R_NamesSymbol,names);
    SET_VECTOR_ELT(out,0,to_list(&d->atoms)); SET_VECTOR_ELT(out,1,to_list(&d->separators));
    RCHECK(igraph_get_edgelist(&d->tree,&d->h,0)); SET_VECTOR_ELT(out,2,to_vector(&d->h));
    SET_VECTOR_ELT(out,3,to_list(&d->ordering)); UNPROTECT(2); return out;
}
static SEXP invoke(call_data *d) {
    SEXP token=PROTECT(R_MakeUnwindCont());
    SEXP out=R_UnwindProtect(run,d,cleanup,d,token);
    UNPROTECT(1); return out;
}
static int inputs(SEXP edges,SEXP directed,SEXP nv) {
    if(XLENGTH(nv)!=1) Rf_error("n_vertices must be scalar");
    int n=whole(nv,0,"n_vertices"); if(n<0) Rf_error("n_vertices must be nonnegative");
    if(XLENGTH(directed)!=1 || Rf_asLogical(directed)!=FALSE) Rf_error("directed graphs are unsupported");
    validate_edges(edges,n); return n;
}
static SEXP c_decompose(SEXP edges,SEXP directed,SEXP nv) {
    call_data d={0}; d.n=inputs(edges,directed,nv); d.edges=edges; return invoke(&d);
}
static SEXP c_collapsible(SEXP edges,SEXP directed,SEXP nv,SEXP seeds) {
    call_data d={0}; d.n=inputs(edges,directed,nv); d.edges=edges; d.seeds=seeds; d.mode=1;
    if(!XLENGTH(seeds)) Rf_error("r_nodes must not be empty");
    for(R_xlen_t i=0;i<XLENGTH(seeds);++i) {
        int z=whole(seeds,i,"r_nodes"); if(z<0 || z>=d.n) Rf_error("r_nodes outside graph");
    }
    return invoke(&d);
}
/* Internal test seam exercises invalid trees and unsorted atom sets directly. */
static SEXP c_ordering(SEXP atoms,SEXP edges) {
    call_data d={0}; if(TYPEOF(atoms)!=VECSXP || XLENGTH(atoms)>INT_MAX) Rf_error("atoms must be a list");
    d.n=(int)XLENGTH(atoms); d.edges=edges; d.atom_input=atoms; d.mode=2; validate_edges(edges,d.n);
    for(R_xlen_t i=0;i<XLENGTH(atoms);++i) for(R_xlen_t j=0;j<XLENGTH(VECTOR_ELT(atoms,i));++j)
        if(whole(VECTOR_ELT(atoms,i),j,"atoms")<0 || whole(VECTOR_ELT(atoms,i),j,"atoms")==INT_MAX)
            Rf_error("invalid atom vertex (must fit one-based R integer output)");
    return invoke(&d);
}
static const R_CallMethodDef calls[]={
    {"_gAtoms_decompose_atoms_impl",(DL_FUNC)&c_decompose,3},
    {"_gAtoms_get_minimal_collapsible_impl",(DL_FUNC)&c_collapsible,4},
    {"_gAtoms_build_prime_ordering_impl",(DL_FUNC)&c_ordering,2},
    {NULL,NULL,0}
};
void attribute_visible R_init_gAtoms(DllInfo *dll) {
    R_registerRoutines(dll,NULL,calls,NULL,NULL); R_useDynamicSymbols(dll,FALSE); R_forceSymbols(dll,TRUE);
}
