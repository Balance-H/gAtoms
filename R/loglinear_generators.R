.loglinear_generators <- function(graph) {
  cliques <- lapply(igraph::max_cliques(graph),function(v) sort(as.integer(v)))
  keys <- vapply(cliques,.set_key,"")
  cliques <- cliques[!duplicated(keys)][order(keys[!duplicated(keys)],method="radix")]
  if (!setequal(unlist(cliques),seq_len(igraph::vcount(graph))) ||
      !all(vapply(cliques,function(s) .is_complete(graph,s),logical(1))))
    .loglinear_failure("invalid_generator_margins","Invalid local graphical generators.")
  cliques
}
