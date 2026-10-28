.loglinear_control <- function(control) {
  defaults <- list(tol=1e-8, contract_tol=1e-6, maxit=1000L,
                   local_table_max_cells=1e6, full_table_max_cells=1e6)
  if (!is.list(control) || (length(control) && (is.null(names(control)) ||
      anyNA(names(control)) || anyDuplicated(names(control)) || any(!names(control)%in%names(defaults)))))
    stop("control must be a named list of recognized options.",call.=FALSE)
  defaults[names(control)] <- control
  for (v in defaults) if (!is.numeric(v) || length(v)!=1L || !is.finite(v) || v<=0)
    stop("Control values must be positive finite scalars.",call.=FALSE)
  if (defaults$maxit!=trunc(defaults$maxit) || defaults$maxit>.Machine$integer.max)
    stop("maxit must be a positive integer.",call.=FALSE)
  defaults
}

.fit_loglinear_loglin <- function(observed, generators, control) {
  # loglin also allocates a 2^q parameter bookkeeping vector, even for one-level variables.
  if (2^length(dim(observed))>control$local_table_max_cells)
    .loglinear_failure("solver_failed","loglin parameter bookkeeping exceeds the local cell limit.")
  warnings <- character()
  z <- withCallingHandlers(stats::loglin(observed,margin=generators,fit=TRUE,
    eps=control$tol,iter=as.integer(control$maxit),print=FALSE),warning=function(w) {
      warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")
    })
  if (any(grepl("converg",warnings,ignore.case=TRUE)))
    .loglinear_failure("nonconverged",paste(warnings,collapse="; "))
  list(fitted_counts=z$fit,iterations=NA_integer_,message=paste(warnings,collapse="; "))
}

.fit_loglinear_grim <- function(observed, generators, control) {
  if (!requireNamespace("gRim",quietly=TRUE) || !exists("effloglin",asNamespace("gRim"),inherits=FALSE))
    .loglinear_failure("solver_failed","gRim::effloglin is unavailable.")
  local <- observed
  vn <- paste0("V",seq_along(dim(local))); names(dimnames(local)) <- vn
  margin <- lapply(generators,function(s) vn[s])
  sufficient <- lapply(generators,function(s) base::margin.table(local,s))
  warnings <- character()
  # Installed effloglin emits progress even with print=FALSE; do not expose it.
  invisible(utils::capture.output(z <- withCallingHandlers(
    gRim::effloglin(sufficient,margin=margin,fit=TRUE,eps=control$tol,
      iter=as.integer(control$maxit),print=FALSE),warning=function(w) {
        warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")
      })))
  a <- z$fit
  if (is.null(dim(a)) || !setequal(names(dimnames(a)),vn))
    .loglinear_failure("solver_failed","gRim returned invalid fitted dimensions.")
  a <- aperm(a,match(vn,names(dimnames(a))))
  names(dimnames(a)) <- names(dimnames(observed))
  list(fitted_counts=a,iterations=NA_integer_,message=paste(warnings,collapse="; "))
}

.validate_loglinear_fit <- function(fitted, observed, generators, n, tol) {
  if (!is.numeric(fitted) || !identical(dim(fitted),dim(observed)) ||
      !identical(dimnames(fitted),dimnames(observed)) || any(!is.finite(fitted)) || any(fitted<0))
    .loglinear_failure("invalid_counts","Local engine returned invalid counts or level dimensions.")
  if (abs(sum(fitted)-n)>tol*max(1,n))
    .loglinear_failure("invalid_counts","Local fitted counts do not sum to n.")
  for (s in generators) {
    error <- max(abs(base::margin.table(fitted,s)-base::margin.table(observed,s)))
    if (error>tol*max(1,n))
      .loglinear_failure("invalid_generator_margins","Local fitted generator margins fail the IPF equations.")
  }
  if (any(fitted[observed>0]<=0))
    .loglinear_failure("solver_failed","A positive observed count has zero fitted probability.")
  prob <- fitted/n
  if (any(!is.finite(prob)) || any(prob<0) || abs(sum(prob)-1)>tol)
    .loglinear_failure("invalid_counts","Invalid local fitted probabilities.")
  invisible(TRUE)
}

# Build an atom_loglinear_error directly so that provenance fields can be
# attached before the condition is signalled.
.engine_failure <- function(status, message, engine = NA_character_) {
  structure(class=c("atom_loglinear_error","error","condition"),
    list(message=message,call=NULL,status=status,atom_id=NA_integer_,
      vertices=character(),engine=engine))
}

.fit_loglinear_atom <- function(observed, graph, complete, engine, control,
                                atom_id=1L, vertices=names(dimnames(observed))) {
  start <- proc.time()[[3L]]; used <- engine
  attempted <- character(); engine_errors <- character(); first_status <- NULL
  tryCatch({
    n <- sum(observed)
    generators <- if (complete) list(seq_along(dim(observed))) else .loglinear_generators(graph)
    if (complete) {
      used <- "closed_form"; method <- "saturated"
      z <- list(fitted_counts=observed,iterations=0L,message="")
      attempted <- character(); fell_back <- FALSE
      reason <- NA_character_; fallback_message <- ""
    } else {
      run <- function(which) {
        z <- if (which=="gRim") .fit_loglinear_grim(observed,generators,control) else
          .fit_loglinear_loglin(observed,generators,control)
        .validate_loglinear_fit(z$fitted_counts,observed,generators,n,control$contract_tol)
        z
      }
      # The engine selected by the caller is used as given; an explicit engine
      # never falls back to the other one.
      attempted <- engine
      z <- tryCatch(run(engine),error=identity)
      if (inherits(z,"error")) {
        first_status <- z$status; engine_errors <- conditionMessage(z)
        stop(z)
      }
      used <- engine; fell_back <- FALSE
      reason <- NA_character_; fallback_message <- ""
      method <- if (used=="gRim") "effloglin" else "IPF"
    }
    .validate_loglinear_fit(z$fitted_counts,observed,generators,n,control$contract_tol)
    list(atom_id=atom_id,vertices=vertices,complete=complete,engine=used,method=method,
      observed=observed,fitted_counts=z$fitted_counts,fitted_prob=z$fitted_counts/n,
      generators=lapply(generators,function(s) vertices[s]),iterations=z$iterations,
      converged=TRUE,status=if (complete) "closed_form" else "converged",
      requested_engine=engine,attempted_engines=attempted,fallback=fell_back,
      fallback_reason=reason,fallback_message=fallback_message,
      message=z$message,elapsed=unname(proc.time()[[3L]]-start))
  },error=function(e) {
    if (!inherits(e,"atom_loglinear_error"))
      e <- .engine_failure("solver_failed",conditionMessage(e),engine)
    if (!length(attempted)) attempted <- if (is.na(used)) engine else used
    if (!length(engine_errors)) engine_errors <- conditionMessage(e)
    names(engine_errors) <- attempted
    e$atom_id <- atom_id; e$vertices <- vertices
    # The requested engine is the only engine that was executed.
    e$engine <- used
    e$requested_engine <- engine; e$attempted_engines <- attempted
    e$engine_errors <- engine_errors
    e$message <- engine_errors[[1L]]
    stop(e)
  })
}
