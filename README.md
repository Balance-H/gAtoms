# gAtoms

gAtoms decomposes simple undirected graphs into atoms using the existing C +
igraph RDA/CMSA backend. The structural API validates inputs, preserves vertex
labels and returns component trees and deterministic perfect orderings.

On top of that shared backend the package provides three statistical modules,
each of which accepts the same `atom_decomp` object through a `decomp`
argument: fixed-graph Gaussian concentration models (`fit_gaussian_graph()`),
graphical log-linear models (`fit_loglinear_graph()`), and decomposition-aware
G-Wishart log normalizing constants (`gwishart_log_norm()`).

```r
library(gAtoms)
g <- igraph::make_ring(6, circular = FALSE)
dec <- atom_decompose(g)
print(dec)
atoms(dec)
separators(dec)
prime_ordering(dec)
stats::residuals(dec)
summary(dec)
plot(dec)
validate(dec, expensive = TRUE)
decomp_cost(dec)
```

## Inputs and indexing

Inputs may be igraph objects, base binary adjacency matrices or
`Matrix::dgCMatrix` objects. Sparse inputs remain sparse. Directed graphs,
loops, multiple edges, nonbinary/nonsymmetric matrices, NA structure, invalid
dimensions and duplicate/missing vertex names are rejected rather than repaired.
Matrix row and column labels must agree when both are present.

Public vertex sets preserve original names; unnamed inputs use one-based integer
positions. `vertex_map` records `internal_id` (zero-based C ID), `index` (R
position), and `label` (original name or the position converted to text).
The explicit `_index` fields always use one-based R positions.

```r
sparse <- as_dgCMatrix(dec)
g_again <- as_igraph(sparse)
```

## Ordering and object fields

For a connected graph with k atoms, `prime_ordering` is the vertex-set list
**P1, S2, P2, ..., Sk, Pk**, of length 2k-1. After RDA finishes the atom tree,
the C helper repeatedly peels the smallest-index active leaf, records its atom
and its intersection with its active neighbour, then reverses the entire
sequence. No tree vertex is physically deleted or renumbered.

`atom_decompose()` returns class `atom_decomp`. `atoms` are ordered by that
sequence; `separators` has one entry per nonroot atom, mapped by
`separator_atoms`. `residuals`, `parent`, `children` and `complete` align with
ordered atoms. Root parents are NA. `tree` is an edge matrix using those atom
IDs; `trees` holds component igraph trees. `components`, `native_atom_ids`,
`diagnostics`, `call` and `version` retain provenance and mappings.

The expert `decompose_atoms()` retains `atoms`, raw RDA `separators`, `tree`
and earlier summary fields, and adds `prime_ordering`. Raw separators do not
necessarily follow tree-edge order. The high-level object's primary separators
are final parent-child intersections; raw sets live under
`diagnostics$raw_separators_index`.

Disconnected graphs have one `component_orderings[[i]]` and tree per graph
component, with no artificial connecting edges or empty separators. Their
single `prime_ordering` field is NULL; the accessor `prime_ordering(dec)` returns
the component list. An isolated vertex has one atom and ordering length one;
the empty graph has empty lists.

`validate()` returns TRUE or a structured list of error messages. It checks
coverage, mappings, forest topology, clique separators, perfect-ordering
intersections, running intersection, parent/children, residuals and complete
flags. Expensive mode adds graph-edge coverage and atom maximality checks.

## Build and development

Install R and a C compiler. That is all: the igraph C library is **bundled** with
this package (`src/igraph`, version 1.0.1) and linked statically into
`gAtoms`. No external igraph installation, shared library, `pkg-config` file,
vcpkg tree or environment variable is needed, on any platform. OpenMP is not
used by gAtoms.

Install with:

```r
install.packages(".", repos = NULL, type = "source")
```

The bundled copy is configured once, offline, with igraph's own CMake build
(`BUILD_SHARED_LIBS=OFF`, `IGRAPH_GLPK_SUPPORT=OFF`, `IGRAPH_OPENMP_SUPPORT=OFF`,
vendored ARPACK/BLAS/LAPACK/GMP/PLFIT/QHULL/Infomap, no libxml2). CMake is a
development-time tool only: it is not used when installing gAtoms. The whole
igraph source tree ships, but only the 189 translation units the package actually
needs are compiled; that selection was derived from the link-time dependency
closure of the gAtoms shared object, not guessed. The object list is recorded in
`src/Makevars.igraph`, which `src/Makevars` and `src/Makevars.win` both follow,
so all three platforms share one build description. See `inst/NOTICE` for the
bundled components, their licenses, the exact configuration and how the subset
was derived.

Because gAtoms links the GPL-2.0-or-later igraph library into its own shared
object, the package as a whole is distributed under **GPL-2.0-or-later**
(`inst/licenses/` collects the license texts of igraph and of the third-party
components igraph itself vendors).

`tools/check.ps1` and `build.bat` are development helpers for local builds and
checks.

CI targets Linux, Windows and macOS. Unit tests include optional comparisons
against `gRbase::mpd()` on 300 seeded graphs; gRbase is not required by the core.
See CONTRIBUTING.md and replication/benchmark.R for development and benchmarking.

The older `fit_ggm_localized()` / `atom_ggm_fit` interface has been removed and
is no longer part of the package; the current Gaussian API is
`gaussian_suffstats()`, `fit_gaussian_graph()` and `precision()`.

## Gaussian concentration models

Fit a fixed graph using `fit_gaussian_graph(X, G)`, optionally supplying
`decomp = atom_decompose(G)`. Complete atoms use closed form; incomplete atoms
use `engine = "gRim"` (the default) or `engine = "ggm"`, and the requested engine
is used as given. The new API uses SSD/n,
requires complete finite observations, and applies no automatic ridge.
Use `precision(fit)`, `vcov(fit)`, `logLik(fit)`, `AIC(fit)` and `BIC(fit)`.
`timing$decomposition` measures only `atom_decompose(graph)` and is exactly zero
when a `decomp` is supplied.
See `vignette("gaussian-graphs", package = "gAtoms")` for a worked example.

## Graphical log-linear models

Use `fit_loglinear_graph(dat, G)` for complete categorical data with column names
matching graph vertex names. `margin_counts()` builds reusable atom/separator
margins directly from rows, preserving declared factor levels and sampling zeros.
Complete atoms use closed form; incomplete atoms use the requested engine,
`engine = "gRim"` (default, `gRim::effloglin()`) or `engine = "loglin"`
(`stats::loglin()`), with no automatic fallback between engines.
`timing$decomposition` measures only `atom_decompose(graph)` and is exactly zero
when a `decomp` is supplied.
`fitted(fit)` returns factors; explicit `fitted(fit, type="full")` is guarded by
a cell limit and returns probabilities. See `vignette("loglinear-graphs",
package="gAtoms")` and `replication/loglinear-benchmark.R`.

## G-Wishart normalizing constants

`gwishart_log_norm(G, b, D)` computes `log I_G(b, D)` for the G-Wishart
distribution with `p(K | G, b, D) propto |K|^((b-2)/2) exp(-tr(K D)/2)`, the same
parameterization used by `BDgraph::gnorm()`. `b` must be an integer of at least
3, and `D` a symmetric positive definite scale matrix ordered like the graph
vertices (named `D` is realigned by name). Complete prime components and clique
separators use the exact complete-graph closed form, so fully decomposable
graphs are exact and need no BDgraph installation; an incomplete prime
component is the only term passed to the optional `BDgraph::gnorm()`, and the
returned object then reports `approximate = TRUE` with the BDgraph call count in
`diagnostics`. See `vignette("gwishart-normalizing-constants",
package="gAtoms")`.
