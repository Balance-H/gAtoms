#include "native_scope.h"
/* The three igraph calls below take arguments that were introduced in the
 * 1.0 series. gAtoms targets the bundled 1.0.1 library, so the 1.0.1 forms
 * are used directly:
 *
 *   igraph_neighbors()  takes loops + multiple.  IGRAPH_NO_LOOPS with
 *       IGRAPH_MULTIPLE reproduces the historical behaviour: self-loops are
 *       dropped and parallel edges are kept, so a repeated neighbour still
 *       appears once per edge.
 *   igraph_degree()     takes an igraph_loops_t.  IGRAPH_NO_LOOPS ignores
 *       self-loops (the enum value is 0, the previous boolean meaning).
 *   igraph_is_simple()  takes a directed flag.  IGRAPH_DIRECTED judges
 *       simplicity with edge directions considered, which is exactly what the
 *       0.10 implementation did (it always walked IGRAPH_OUT and therefore did
 *       not treat a mutual edge pair as a multi-edge). */
/* Original RDA/CMSA: MCS selects a closed neighbourhood, CMSA expands via
 * nonadjacent boundary pairs and two close separators, RDA distributes
 * (separator,parent) pairs. Serial execution, local ownership, checked errors.
 * Caller installs a returning igraph error handler (the R wrapper does). */
void free_vertex_list(igraph_vector_ptr_t *v) {
    for(igraph_integer_t i=0;i<igraph_vector_ptr_size(v);++i) ga_vec_free(VECTOR(*v)[i]);
    igraph_vector_ptr_destroy(v);
}
static igraph_error_t unique_vec(igraph_vector_int_t *v) {
    igraph_integer_t n=igraph_vector_int_size(v),w=0;
    igraph_vector_int_sort(v);
    for(igraph_integer_t i=0;i<n;++i)
        if(!w || VECTOR(*v)[i]!=VECTOR(*v)[w-1]) VECTOR(*v)[w++]=VECTOR(*v)[i];
    return igraph_vector_int_resize(v,w);
}
static igraph_bool_t contains(const igraph_vector_int_t *v,igraph_integer_t x) {
    for(igraph_integer_t i=0;i<igraph_vector_int_size(v);++i) if(VECTOR(*v)[i]==x) return 1;
    return 0;
}
static igraph_bool_t subset(const igraph_vector_int_t *s,const igraph_vector_int_t *u) {
    for(igraph_integer_t i=0;i<igraph_vector_int_size(s);++i) if(!contains(u,VECTOR(*s)[i])) return 0;
    return 1;
}
/* ---------------------------------------------------------------------------
 * igraph 1.0.1 subgraph mapping.
 *
 * igraph_induced_subgraph_map() returns two maps for the subgraph on `vids`:
 *
 *   forward map (map):    indexed by the ORIGINAL graph's vertex id; the value
 *                         is the 0-based index of that vertex inside the
 *                         subgraph, or -1 when the vertex was not kept.  Its
 *                         length is igraph_vcount(graph).
 *   inverse map (invmap): indexed by the subgraph's 0-based vertex index; the
 *                         value is the ORIGINAL graph's 0-based vertex id.
 *                         Its length is igraph_vcount(subgraph).
 *
 * gAtoms addresses vertices 0-based internally and only converts to 1-based
 * IDs at the R boundary, so both maps are used exactly as igraph defines them,
 * with no offset applied on either side.  Applying a shift here would silently
 * mis-attribute vertices and break the decomposition.
 *
 * Verified against the bundled 1.0.1 library with a direct C probe
 * (graph v=4, keep {0,1,2} -> subgraph v=3):
 *     map {0,1,2,-1}   invmap {0,1,2}
 * --------------------------------------------------------------------------- */

/* Forward map: original-graph vertex id -> 0-based subgraph index, or -1. */
static igraph_integer_t ga_map_index(const igraph_vector_int_t *map,igraph_integer_t v) {
    return VECTOR(*map)[v];
}
/* Inverse map: 0-based subgraph index -> original-graph 0-based vertex id. */
static igraph_integer_t ga_inv_vertex(const igraph_vector_int_t *inv,igraph_integer_t local) {
    return VECTOR(*inv)[local];
}
static igraph_error_t valid_ids(const igraph_t *g,const igraph_vector_int_t *v) {
    if(v) for(igraph_integer_t i=0;i<igraph_vector_int_size(v);++i)
        if(VECTOR(*v)[i]<0 || VECTOR(*v)[i]>=igraph_vcount(g)) return IGRAPH_EINVVID;
    return IGRAPH_SUCCESS;
}
static igraph_error_t copy_to_list(igraph_vector_ptr_t *out,const igraph_vector_int_t *v) {
    igraph_vector_int_t *p=malloc(sizeof(*p)); if(!p) return IGRAPH_ENOMEM;
    igraph_error_t e=igraph_vector_int_init_copy(p,v);
    if(e) { free(p); return e; }
    e=igraph_vector_ptr_push_back(out,p); if(e) ga_vec_free(p); return e;
}
static igraph_error_t map_vertices(igraph_vector_int_t *out,const igraph_vector_int_t *v,const igraph_vector_int_t *inv) {
    igraph_error_t e=igraph_vector_int_resize(out,igraph_vector_int_size(v)); if(e) return e;
    for(igraph_integer_t i=0;i<igraph_vector_int_size(v);++i) {
        igraph_integer_t x=VECTOR(*v)[i];
        if(x<0 || x>=igraph_vector_int_size(inv)) return IGRAPH_EINVVID;
        VECTOR(*out)[i]=ga_inv_vertex(inv,x);
    }
    return unique_vec(out);
}
igraph_error_t components_forbidden(const igraph_t *g,igraph_vector_ptr_t *comps,
    igraph_vector_ptr_t *bounds,const igraph_vector_int_t *forbidden) {
    ga_scope sc={0}; igraph_error_t ret=IGRAPH_SUCCESS; igraph_integer_t n=igraph_vcount(g);
    TRY(valid_ids(g,forbidden));
    VEC(sc,visited,n); VEC(sc,marked,n); VEC(sc,queue,n);
    VEC(sc,nei,0); VEC(sc,comp,0); VEC(sc,bound,0);
    if(forbidden) for(igraph_integer_t i=0;i<igraph_vector_int_size(forbidden);++i) VECTOR(*visited)[VECTOR(*forbidden)[i]]=2;
    for(igraph_integer_t v=0;v<n;++v) {
        if(VECTOR(*visited)[v]) continue;
        igraph_vector_int_clear(comp); igraph_vector_int_clear(bound);
        igraph_integer_t head=0,tail=0; VECTOR(*queue)[tail++]=v; VECTOR(*visited)[v]=1;
        while(head<tail) {
            igraph_integer_t cur=VECTOR(*queue)[head++];
            TRY(igraph_vector_int_push_back(comp,cur)); TRY(igraph_neighbors(g,nei,cur,IGRAPH_ALL,IGRAPH_NO_LOOPS,IGRAPH_MULTIPLE));
            for(igraph_integer_t i=0;i<igraph_vector_int_size(nei);++i) {
                igraph_integer_t w=VECTOR(*nei)[i];
                if(VECTOR(*visited)[w]==2) {
                    if(!VECTOR(*marked)[w]) { TRY(igraph_vector_int_push_back(bound,w)); VECTOR(*marked)[w]=1; }
                } else if(!VECTOR(*visited)[w]) { VECTOR(*visited)[w]=1; VECTOR(*queue)[tail++]=w; }
            }
        }
        TRY(copy_to_list(comps,comp)); TRY(copy_to_list(bounds,bound));
        for(igraph_integer_t i=0;i<igraph_vector_int_size(bound);++i) VECTOR(*marked)[VECTOR(*bound)[i]]=0;
    }
cleanup: ga_clear(&sc); return ret;
}
igraph_error_t close_separator(const igraph_t *g,igraph_integer_t vertex,
    const igraph_vector_int_t *forbidden,igraph_vector_int_t *bound) {
    ga_scope sc={0}; igraph_error_t ret=IGRAPH_SUCCESS; igraph_integer_t n=igraph_vcount(g);
    if(vertex<0 || vertex>=n) return IGRAPH_EINVVID;
    TRY(valid_ids(g,forbidden)); if(forbidden && contains(forbidden,vertex)) return IGRAPH_EINVAL;
    VEC(sc,visited,n); VEC(sc,queue,n); VEC(sc,nei,0);
    if(forbidden) for(igraph_integer_t i=0;i<igraph_vector_int_size(forbidden);++i) VECTOR(*visited)[VECTOR(*forbidden)[i]]=2;
    igraph_vector_int_clear(bound);
    igraph_integer_t head=0,tail=0; VECTOR(*queue)[tail++]=vertex; VECTOR(*visited)[vertex]=1;
    while(head<tail) {
        TRY(igraph_neighbors(g,nei,VECTOR(*queue)[head++],IGRAPH_ALL,IGRAPH_NO_LOOPS,IGRAPH_MULTIPLE));
        for(igraph_integer_t i=0;i<igraph_vector_int_size(nei);++i) {
            igraph_integer_t w=VECTOR(*nei)[i];
            if(VECTOR(*visited)[w]==2) { TRY(igraph_vector_int_push_back(bound,w)); VECTOR(*visited)[w]=1; }
            else if(!VECTOR(*visited)[w]) { VECTOR(*visited)[w]=1; VECTOR(*queue)[tail++]=w; }
        }
    }
cleanup: ga_clear(&sc); return ret;
}
/* Child scope frees each CMSA induced graph immediately. */
static igraph_error_t expand_pair(const igraph_t *g,const igraph_vector_int_t *comp,
    igraph_integer_t hj,igraph_integer_t hk,igraph_vector_int_t *added) {
    ga_scope sc={0}; igraph_error_t ret=IGRAPH_SUCCESS; igraph_t *sub;
    VEC(sc,nodes,0); VEC(sc,map,0); VEC(sc,inv,0); VEC(sc,na,0); VEC(sc,nb,0);
    VEC(sc,sa,0); VEC(sc,sb,0); VEC(sc,global,0);
    TRY(igraph_vector_int_append(nodes,comp)); TRY(igraph_vector_int_push_back(nodes,hj)); TRY(igraph_vector_int_push_back(nodes,hk));
    TRY(ga_subgraph(&sc,g,nodes,&sub,map,inv));
    igraph_integer_t a=ga_map_index(map,hj),b=ga_map_index(map,hk);
    if(a<0 || b<0) { ret=IGRAPH_EINVVID; goto cleanup; }
    TRY(igraph_neighbors(sub,na,a,IGRAPH_ALL,IGRAPH_NO_LOOPS,IGRAPH_MULTIPLE)); TRY(igraph_neighbors(sub,nb,b,IGRAPH_ALL,IGRAPH_NO_LOOPS,IGRAPH_MULTIPLE));
    TRY(close_separator(sub,b,na,sa)); TRY(close_separator(sub,a,nb,sb));
    TRY(map_vertices(global,sa,inv)); TRY(igraph_vector_int_append(added,global));
    TRY(map_vertices(global,sb,inv)); TRY(igraph_vector_int_append(added,global));
cleanup: ga_clear(&sc); return ret;
}
igraph_error_t get_minimal_collapsible(const igraph_t *g,const igraph_vector_int_t *seeds,igraph_vector_int_t *out) {
    ga_scope sc={0},round={0}; igraph_error_t ret=IGRAPH_SUCCESS; igraph_bool_t simple;
    if(igraph_is_directed(g)) return IGRAPH_EINVAL;
    TRY(igraph_is_simple(g,&simple,IGRAPH_DIRECTED)); if(!simple) return IGRAPH_EINVAL;
    TRY(valid_ids(g,seeds)); VEC(sc,h,0); TRY(igraph_vector_int_append(h,seeds)); TRY(unique_vec(h));
    for(;;) {
        LIST(round,comps); LIST(round,bounds); VEC(round,added,0);
        TRY(components_forbidden(g,comps,bounds,h));
        for(igraph_integer_t c=0;c<igraph_vector_ptr_size(comps);++c) {
            igraph_vector_int_t *bound=VECTOR(*bounds)[c]; igraph_bool_t found=0; igraph_vector_int_sort(bound);
            for(igraph_integer_t j=0;j<igraph_vector_int_size(bound) && !found;++j)
                for(igraph_integer_t k=j+1;k<igraph_vector_int_size(bound) && !found;++k) {
                    igraph_bool_t adjacent;
                    TRY(igraph_are_adjacent(g,VECTOR(*bound)[j],VECTOR(*bound)[k],&adjacent));
                    if(!adjacent) { TRY(expand_pair(g,VECTOR(*comps)[c],VECTOR(*bound)[j],VECTOR(*bound)[k],added)); found=1; }
                }
        }
        igraph_integer_t old=igraph_vector_int_size(h);
        TRY(igraph_vector_int_append(h,added)); TRY(unique_vec(h)); ga_clear(&round);
        if(igraph_vector_int_size(h)==old) break;
    }
    igraph_vector_int_clear(out); TRY(igraph_vector_int_append(out,h));
cleanup: ga_clear(&round); ga_clear(&sc); return ret;
}
/* Parallel vectors replace pair_t allocations, preserving RDA pair semantics. */
typedef struct { igraph_vector_int_t vertices,parents; igraph_vector_ptr_t separators; } task_t;
static void task_free(task_t *t) {
    igraph_vector_int_destroy(&t->vertices); igraph_vector_int_destroy(&t->parents); free_vertex_list(&t->separators); free(t);
}
static igraph_error_t task_new(task_t **out,const igraph_vector_int_t *nodes) {
    task_t *t=malloc(sizeof(*t)); igraph_error_t e; if(!t) return IGRAPH_ENOMEM;
    e=igraph_vector_int_init_copy(&t->vertices,nodes); if(e) { free(t); return e; }
    e=igraph_vector_ptr_init(&t->separators,0);
    if(e) { igraph_vector_int_destroy(&t->vertices); free(t); return e; }
    e=igraph_vector_int_init(&t->parents,0);
    if(e) { free_vertex_list(&t->separators); igraph_vector_int_destroy(&t->vertices); free(t); return e; }
    *out=t; return IGRAPH_SUCCESS;
}
static igraph_error_t task_pair(task_t *t,const igraph_vector_int_t *s,igraph_integer_t p) {
    igraph_error_t e=copy_to_list(&t->separators,s); if(e) return e;
    return igraph_vector_int_push_back(&t->parents,p);
}
static igraph_error_t process_task(const igraph_t *g,const igraph_vector_int_t *index,
    task_t *task,igraph_vector_ptr_t *pending,igraph_vector_ptr_t *atoms,igraph_vector_ptr_t *separators,igraph_t *tree) {
    ga_scope sc={0}; igraph_error_t ret=IGRAPH_SUCCESS; igraph_t *sub; task_t *child=NULL;
    VEC(sc,map,0); VEC(sc,inv,0); VEC(sc,nv,0); VEC(sc,h,0); VEC(sc,hglobal,0);
    VEC(sc,cglobal,0); VEC(sc,anew,0); VEC(sc,aglobal,0); LIST(sc,comps); LIST(sc,bounds);
    TRY(ga_subgraph(&sc,g,&task->vertices,&sub,map,inv));
    igraph_integer_t v=VECTOR(task->vertices)[0];
    for(igraph_integer_t i=1;i<igraph_vector_int_size(&task->vertices);++i) {
        igraph_integer_t w=VECTOR(task->vertices)[i]; if(VECTOR(*index)[w]<VECTOR(*index)[v]) v=w;
    }
    igraph_integer_t local=ga_map_index(map,v); if(local<0) { ret=IGRAPH_EINVVID; goto cleanup; }
    TRY(igraph_neighbors(sub,nv,local,IGRAPH_ALL,IGRAPH_NO_LOOPS,IGRAPH_MULTIPLE)); TRY(igraph_vector_int_push_back(nv,local)); TRY(unique_vec(nv));
    igraph_bool_t clique; TRY(igraph_is_clique(sub,igraph_vss_vector(nv),0,&clique));
    if(clique) { TRY(igraph_vector_int_append(h,nv)); } else { TRY(get_minimal_collapsible(sub,nv,h)); }
    TRY(map_vertices(hglobal,h,inv)); TRY(copy_to_list(atoms,hglobal));
    igraph_integer_t u=igraph_vector_ptr_size(atoms)-1; TRY(igraph_add_vertices(tree,1,NULL));
    for(igraph_integer_t i=0;i<igraph_vector_ptr_size(&task->separators);++i)
        if(subset(VECTOR(task->separators)[i],hglobal)) {
            TRY(igraph_add_edge(tree,u,VECTOR(task->parents)[i])); VECTOR(task->parents)[i]=-1;
        }
    TRY(components_forbidden(sub,comps,bounds,h));
    for(igraph_integer_t c=0;c<igraph_vector_ptr_size(comps);++c) {
        TRY(map_vertices(cglobal,VECTOR(*bounds)[c],inv)); TRY(copy_to_list(separators,cglobal));
        igraph_vector_int_clear(anew); TRY(igraph_vector_int_append(anew,VECTOR(*comps)[c]));
        TRY(igraph_vector_int_append(anew,VECTOR(*bounds)[c])); TRY(map_vertices(aglobal,anew,inv));
        TRY(task_new(&child,aglobal)); TRY(task_pair(child,cglobal,u));
        for(igraph_integer_t i=0;i<igraph_vector_ptr_size(&task->separators);++i)
            if(VECTOR(task->parents)[i]>=0 && subset(VECTOR(task->separators)[i],aglobal))
                TRY(task_pair(child,VECTOR(task->separators)[i],VECTOR(task->parents)[i]));
        TRY(igraph_vector_ptr_push_back(pending,child)); child=NULL;
    }
cleanup: if(child) task_free(child); ga_clear(&sc); return ret;
}
igraph_error_t decompose_atoms(const igraph_t *g,igraph_vector_ptr_t *atoms,
    igraph_vector_ptr_t *separators,igraph_t *tree,igraph_vector_ptr_t *ordering) {
    ga_scope sc={0}; igraph_error_t ret=IGRAPH_SUCCESS; igraph_vector_ptr_t pending;
    igraph_bool_t pending_ok=0,simple,connected; task_t *task=NULL;
    if(igraph_is_directed(g) || igraph_is_directed(tree) || igraph_vcount(tree) ||
       igraph_vector_ptr_size(atoms) || igraph_vector_ptr_size(separators) || igraph_vector_ptr_size(ordering)) return IGRAPH_EINVAL;
    TRY(igraph_is_simple(g,&simple,IGRAPH_DIRECTED)); if(!simple) return IGRAPH_EINVAL;
    if(!igraph_vcount(g)) return build_prime_ordering(atoms,tree,ordering);
    TRY(igraph_is_connected(g,&connected,IGRAPH_WEAK)); if(!connected) return IGRAPH_EINVAL;
    VEC(sc,alpha,0); VEC(sc,order,0); VEC(sc,index,igraph_vcount(g)); VEC(sc,full,igraph_vcount(g));
    TRY(igraph_maximum_cardinality_search(g,alpha,order));
    for(igraph_integer_t i=0;i<igraph_vcount(g);++i) { VECTOR(*index)[VECTOR(*order)[i]]=i; VECTOR(*full)[i]=i; }
    TRY(igraph_vector_ptr_init(&pending,0)); pending_ok=1;
    TRY(task_new(&task,full)); TRY(igraph_vector_ptr_push_back(&pending,task)); task=NULL;
    while(!igraph_vector_ptr_empty(&pending)) {
        task=igraph_vector_ptr_pop_back(&pending);
        TRY(process_task(g,index,task,&pending,atoms,separators,tree)); task_free(task); task=NULL;
    }
    TRY(build_prime_ordering(atoms,tree,ordering));
cleanup:
    if(task) task_free(task);
    if(pending_ok) {
        for(igraph_integer_t i=0;i<igraph_vector_ptr_size(&pending);++i) task_free(VECTOR(pending)[i]);
        igraph_vector_ptr_destroy(&pending);
    }
    ga_clear(&sc); return ret;
}
/* Leaf peeling is separate from RDA. Stable original atom IDs, logical removal,
 * sorted intersections and deep copies, followed by reversal of ALL entries. */
igraph_error_t build_prime_ordering(const igraph_vector_ptr_t *atoms,const igraph_t *tree,igraph_vector_ptr_t *ordering) {
    ga_scope sc={0}; igraph_error_t ret=IGRAPH_SUCCESS; igraph_bool_t connected,simple;
    igraph_integer_t k=igraph_vector_ptr_size(atoms);
    if(igraph_vector_ptr_size(ordering) || igraph_is_directed(tree) || igraph_vcount(tree)!=k || igraph_ecount(tree)!=(k?k-1:0)) return IGRAPH_EINVAL;
    if(!k) return IGRAPH_SUCCESS;
    TRY(igraph_is_connected(tree,&connected,IGRAPH_WEAK)); TRY(igraph_is_simple(tree,&simple,IGRAPH_DIRECTED));
    if(!connected || !simple) return IGRAPH_EINVAL;
    VEC(sc,active,k); VEC(sc,degree,k); VEC(sc,nei,0); VEC(sc,sep,0);
    TRY(igraph_degree(tree,degree,igraph_vss_all(),IGRAPH_ALL,IGRAPH_NO_LOOPS)); igraph_vector_int_fill(active,1);
    for(igraph_integer_t remaining=k;remaining>1;--remaining) {
        igraph_integer_t leaf=-1,parent=-1;
        for(igraph_integer_t i=0;i<k;++i) if(VECTOR(*active)[i] && VECTOR(*degree)[i]==1) { leaf=i; break; }
        if(leaf<0) { ret=IGRAPH_EINVAL; goto cleanup; }
        TRY(igraph_neighbors(tree,nei,leaf,IGRAPH_ALL,IGRAPH_NO_LOOPS,IGRAPH_MULTIPLE));
        for(igraph_integer_t i=0;i<igraph_vector_int_size(nei);++i) if(VECTOR(*active)[VECTOR(*nei)[i]]) { parent=VECTOR(*nei)[i]; break; }
        if(parent<0) { ret=IGRAPH_EINVAL; goto cleanup; }
        igraph_vector_int_t *l=VECTOR(*atoms)[leaf],*q=VECTOR(*atoms)[parent]; igraph_vector_int_clear(sep);
        for(igraph_integer_t i=0;i<igraph_vector_int_size(l);++i) if(contains(q,VECTOR(*l)[i])) TRY(igraph_vector_int_push_back(sep,VECTOR(*l)[i]));
        TRY(unique_vec(sep)); TRY(copy_to_list(ordering,l)); TRY(copy_to_list(ordering,sep));
        VECTOR(*active)[leaf]=0; --VECTOR(*degree)[parent];
    }
    for(igraph_integer_t i=0;i<k;++i) if(VECTOR(*active)[i]) TRY(copy_to_list(ordering,VECTOR(*atoms)[i]));
    igraph_integer_t len=igraph_vector_ptr_size(ordering);
    for(igraph_integer_t i=0;i<len/2;++i) {
        void *tmp=VECTOR(*ordering)[i]; VECTOR(*ordering)[i]=VECTOR(*ordering)[len-1-i]; VECTOR(*ordering)[len-1-i]=tmp;
    }
cleanup: ga_clear(&sc); return ret;
}
