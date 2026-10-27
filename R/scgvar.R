# Subgroup chain graphical VAR.
scgvar = function(
  data,
  vars,
  idvar,
  beepvar = NULL,
  dayvar = NULL,
  search = TRUE,
  scale = TRUE,
  maxit.in = 100,
  centerWithin = TRUE,
  gamma = 0.5,
  verbose = TRUE,
  redo = FALSE,
  subjectNetworks = TRUE,
  lambda_min_kappa_fixed = 0.001,
  lambda_min_beta_fixed = 0.001,
  lambda_min_kappa = 0.05,
  lambda_min_beta = lambda_min_kappa,
  lambda_min_glasso = 0.01,
  dir = NULL,
  btwnNetworks = FALSE,
  mlgvclust = TRUE,
  indplot = FALSE,
  elim.groups = FALSE,
  g.ID = NULL,
  nclust = NULL,
  returnResults = TRUE,
  ...
){
  Call = match.call()
  if (missing(vars) || missing(idvar)){
    stop("Both 'vars' and 'idvar' must be supplied.", call. = FALSE)
  }
  flags = list(search = search, scale = scale, centerWithin = centerWithin,
               verbose = verbose, redo = redo, btwnNetworks = btwnNetworks,
               mlgvclust = mlgvclust, indplot = indplot,
               elim.groups = elim.groups, returnResults = returnResults)
  for (name in names(flags)){
    value = flags[[name]]
    if (!is.logical(value) || length(value) != 1L || is.na(value)){
      stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if (redo || !mlgvclust || !is.null(g.ID)){
    stop("The experimental 'redo', 'mlgvclust = FALSE', and 'g.ID' options are not supported.",
         call. = FALSE)
  }
  if (!is.numeric(maxit.in) || length(maxit.in) != 1L ||
      !is.finite(maxit.in) || maxit.in < 1 || maxit.in != floor(maxit.in)){
    stop("'maxit.in' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(gamma) || length(gamma) != 1L || !is.finite(gamma) || gamma < 0){
    stop("'gamma' must be a finite nonnegative number.", call. = FALSE)
  }
  ratios = list(lambda_min_kappa_fixed = lambda_min_kappa_fixed,
                lambda_min_beta_fixed = lambda_min_beta_fixed,
                lambda_min_kappa = lambda_min_kappa,
                lambda_min_beta = lambda_min_beta,
                lambda_min_glasso = lambda_min_glasso)
  for (name in names(ratios)){
    value = ratios[[name]]
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value <= 0 || value >= 1){
      stop("'", name, "' must be in (0, 1).", call. = FALSE)
    }
  }
  dots = list(...)
  allowed = setdiff(names(formals(graphicalVAR)),
                    c("data", "vars", "idvar", "beepvar", "dayvar", "scale",
                      "centerWithin", "gamma", "verbose", "maxit.in",
                      "lambda_min_kappa", "lambda_min_beta"))
  if (length(dots) && (is.null(names(dots)) || anyNA(names(dots)) ||
      any(!nzchar(names(dots))) || anyDuplicated(names(dots)) ||
      any(!names(dots) %in% allowed))){
    stop("Arguments in '...' must have unique, exact graphicalVAR argument names.",
         call. = FALSE)
  }
  if (!is.null(dots$lags) && !identical(as.numeric(dots$lags), 1)){
    stop("SCGVAR currently supports only 'lags = 1'.", call. = FALSE)
  }
  dots$lags = 1
  if (!returnResults && is.null(dir)){
    stop("Use 'returnResults = TRUE' or supply 'dir' to retain the results.", call. = FALSE)
  }
  output.dir = NULL
  if (!is.null(dir)){
    if (!is.character(dir) || length(dir) != 1L || is.na(dir) || !nzchar(dir)){
      stop("'dir' must be NULL or a nonempty directory path.", call. = FALSE)
    }
    output.dir = file.path(path.expand(dir), "graphical VAR")
    if (file.exists(output.dir)){
      stop("The output directory already exists; choose a new 'dir': ", output.dir,
           call. = FALSE)
    }
  } else if (indplot){
    stop("'indplot = TRUE' requires an output 'dir'.", call. = FALSE)
  }
  input = .scgvarData(data, vars, idvar, beepvar, dayvar)
  .scgvarArguments(dots, length(vars))
  data = input$data
  beepvar = input$beepvar
  dayvar = input$dayvar
  IDs = unique(data[[idvar]])
  if (!isTRUE(subjectNetworks)){
    if (is.logical(subjectNetworks) || !is.atomic(subjectNetworks) ||
        anyNA(subjectNetworks) || anyDuplicated(subjectNetworks) ||
        length(subjectNetworks) < 2L || any(!subjectNetworks %in% IDs)){
      stop("'subjectNetworks' must be TRUE or at least two distinct observed subject IDs.",
           call. = FALSE)
    }
    IDs = IDs[match(subjectNetworks, IDs)]
  }
  if (length(IDs) < 2L){
    stop("At least two subjects are needed for subgroup clustering.", call. = FALSE)
  }
  if (!is.null(nclust) && (!is.numeric(nclust) || length(nclust) != 1L ||
      !is.finite(nclust) || nclust != floor(nclust) ||
      nclust < 1 || nclust > length(IDs))){
    stop("'nclust' must be NULL or an integer from 1 to the number of selected subjects.",
         call. = FALSE)
  }
  settings = list(vars = vars, idvar = idvar, beepvar = beepvar, dayvar = dayvar,
                  scale = scale, centerWithin = centerWithin, gamma = gamma,
                  maxit.in = maxit.in, verbose = verbose, dots = dots,
                  lambda_min_kappa_fixed = lambda_min_kappa_fixed,
                  lambda_min_beta_fixed = lambda_min_beta_fixed,
                  lambda_min_kappa = lambda_min_kappa,
                  lambda_min_beta = lambda_min_beta,
                  lambda_min_glasso = lambda_min_glasso,
                  btwnNetworks = btwnNetworks)
  if (verbose) message("Estimating pooled and subject-specific networks")
  group.output = .scgvarFit(data, IDs, settings, "Pooled sample")
  similarity = .scgvarSimilarity(group.output$subjectPCC, group.output$subjectPDC,
                                 group.output$fixedPCC, group.output$fixedPDC,
                                 elim.groups)
  dimnames(similarity) = list(as.character(IDs), as.character(IDs))
  if (verbose) message("Detecting subgroups")
  clusters = .scgvarCluster(similarity, search, nclust)
  membership = data.frame(IDs = IDs, membership = unname(clusters$membership))
  groups = sort(unique(membership$membership))
  subgroup.outputs = vector("list", length(groups))
  names(subgroup.outputs) = as.character(groups)
  for (i in seq_along(groups)){
    label = paste("Subgroup", groups[i])
    if (verbose) message("Estimating ", tolower(label), " networks")
    selected = IDs[membership$membership == groups[i]]
    subdata = data[data[[idvar]] %in% selected, , drop = FALSE]
    subgroup.outputs[[i]] = .scgvarFit(subdata, selected[FALSE], settings, label)
  }
  Results = list(group.output = group.output, similarity = similarity,
                 sim.mat = clusters$matrix, membership = membership,
                 subgroup.outputs = subgroup.outputs,
                 clustering = clusters$diagnostics, output.dir = output.dir,
                 call = Call)
  class(Results) = "scgvar"
  if (!is.null(output.dir)){
    .scgvarWrite(Results, vars, indplot)
    Results$output.dir = normalizePath(output.dir, winslash = "/", mustWork = TRUE)
  }
  if (returnResults) return(invisible(Results))
  invisible(NULL)
}

# Validate before fitting. Missing beeps are retained as gaps by tsData().
.scgvarData = function(data, vars, idvar, beepvar, dayvar){
  if (!is.data.frame(data) && !is.matrix(data)){
    stop("'data' must be a matrix or data frame.", call. = FALSE)
  }
  data = as.data.frame(data)
  if (!nrow(data) || is.null(names(data)) || anyNA(names(data)) ||
      any(!nzchar(names(data))) || anyDuplicated(names(data))){
    stop("'data' must have rows and unique, nonempty column names.", call. = FALSE)
  }
  if (!is.character(vars) || length(vars) < 2L || anyNA(vars) ||
      any(!nzchar(vars)) || anyDuplicated(vars)){
    stop("'vars' must contain at least two distinct column names.", call. = FALSE)
  }
  columns = list(idvar = idvar, beepvar = beepvar, dayvar = dayvar)
  for (name in names(columns)){
    value = columns[[name]]
    if (is.null(value) && name != "idvar") next
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)){
      stop("'", name, "' must be a single column name.", call. = FALSE)
    }
  }
  selected = c(vars, idvar, beepvar, dayvar)
  if (anyDuplicated(selected)){
    stop("Measurement, ID, day, and beep columns must be distinct.", call. = FALSE)
  }
  absent = setdiff(selected, names(data))
  if (length(absent)){
    stop("Columns missing from 'data': ", paste(absent, collapse = ", "), call. = FALSE)
  }
  data = data[, selected, drop = FALSE]
  if (anyNA(data)) stop("Missing values are not supported in selected columns.", call. = FALSE)
  if (any(!vapply(data[vars], function(x) is.numeric(x) && all(is.finite(x)), logical(1)))){
    stop("Columns in 'vars' must be finite numeric measurements.", call. = FALSE)
  }
  for (name in c(idvar, dayvar)){
    x = data[[name]]
    if (!is.atomic(x) || is.matrix(x) ||
        (is.numeric(x) && any(!is.finite(x)))){
      stop("ID and day columns must contain finite atomic values.", call. = FALSE)
    }
  }
  # Explicit unique names avoid tsData's reserved ID/DAY/BEEP defaults.
  if (is.null(dayvar)){
    dayvar = utils::tail(make.unique(c(names(data), ".scgvar_day")), 1L)
    data[[dayvar]] = 1L
  }
  if (is.null(beepvar)){
    beepvar = utils::tail(make.unique(c(names(data), ".scgvar_beep")), 1L)
    data[[beepvar]] = integer(nrow(data))
    for (id in unique(data[[idvar]])){
      for (day in unique(data[[dayvar]][data[[idvar]] == id])){
        rows = which(data[[idvar]] == id & data[[dayvar]] == day)
        data[[beepvar]][rows] = seq_along(rows)
      }
    }
  }
  beeps = data[[beepvar]]
  if (!is.numeric(beeps) || any(!is.finite(beeps)) ||
      any(beeps != floor(beeps)) || any(beeps < 0)){
    stop("The beep column must contain nonnegative finite integers.", call. = FALSE)
  }
  if (anyDuplicated(data[, c(idvar, dayvar, beepvar), drop = FALSE])){
    stop("Duplicate subject/day/beep combinations are not supported.", call. = FALSE)
  }
  for (id in unique(data[[idvar]])){
    person = data[data[[idvar]] == id, , drop = FALSE]
    variances = vapply(person[vars], stats::var, numeric(1))
    if (any(!is.finite(variances) | variances <= 0)){
      stop("Subject '", id, "' has a constant or insufficiently observed variable.", call. = FALSE)
    }
    pairs = 0L
    for (day in unique(person[[dayvar]])){
      times = sort(person[[beepvar]][person[[dayvar]] == day])
      pairs = pairs + sum(diff(times) == 1)
    }
    if (pairs < 3L){
      stop("Subject '", id, "' needs at least three usable lag-1 pairs within days.",
           call. = FALSE)
    }
  }
  list(data = data, beepvar = beepvar, dayvar = dayvar)
}

# Prevent invalid controls from entering the numerical solver (including its
# interactive debugging branch for malformed beta regularization matrices).
.scgvarArguments = function(dots, p){
  for (name in intersect(c("nLambda", "maxit.out"), names(dots))){
    x = dots[[name]]
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 || x != floor(x)){
      stop("'", name, "' must be a positive integer.", call. = FALSE)
    }
  }
  for (name in intersect(c("lambda_beta", "lambda_kappa"), names(dots))){
    x = dots[[name]]
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || any(x < 0) ||
        anyDuplicated(x)){
      stop("'", name, "' must contain distinct finite nonnegative tuning parameters.",
           call. = FALSE)
    }
  }
  for (name in intersect(c("regularize_mat_beta", "regularize_mat_kappa"), names(dots))){
    x = dots[[name]]
    cols = if (name == "regularize_mat_beta") c(p, p + 1L) else p
    if (!is.matrix(x) || (!is.numeric(x) && !is.logical(x)) || anyNA(x) ||
        any(!x %in% c(0, 1)) || nrow(x) != p || !ncol(x) %in% cols){
      stop("'", name, "' must be a 0/1 or logical matrix with ", p,
           " rows and ", paste(cols, collapse = " or "), " columns.", call. = FALSE)
    }
    if (name == "regularize_mat_kappa" && !isTRUE(all.equal(x, t(x), check.attributes = FALSE))){
      stop("'regularize_mat_kappa' must be symmetric.", call. = FALSE)
    }
  }
  for (name in intersect(c("deleteMissings", "penalize.diagonal"), names(dots))){
    x = dots[[name]]
    if (!is.logical(x) || length(x) != 1L || is.na(x)){
      stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
    }
  }
  if ("ebic_tol" %in% names(dots)){
    x = dots$ebic_tol
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0){
      stop("'ebic_tol' must be a finite positive number.", call. = FALSE)
    }
  }
  invisible(NULL)
}

# Pooled and individual estimation with optional between-person networks.
.scgvarFit = function(data, IDs, settings, context){
  s = settings
  prepped = tsData(data, vars = s$vars, idvar = s$idvar, beepvar = s$beepvar,
                   dayvar = s$dayvar, scale = s$scale, centerWithin = s$centerWithin)
  estimate = function(input, kappa, beta, label, extra = list()){
    args = c(list(data = input, lambda_min_kappa = kappa,
                  lambda_min_beta = beta, gamma = s$gamma,
                  maxit.in = s$maxit.in, verbose = FALSE), extra, s$dots)
    fit = tryCatch(do.call(graphicalVAR, args), error = function(e){
      stop(label, ": ", conditionMessage(e), call. = FALSE)
    })
    matrices = c("PCC", "PDC", "beta", "kappa")
    valid = vapply(matrices, function(name){
      is.matrix(fit[[name]]) && all(is.finite(fit[[name]]))
    }, logical(1))
    if (!all(valid) || length(fit$EBIC) != 1L || !is.finite(fit$EBIC)){
      stop(label, ": estimation returned non-finite network estimates.", call. = FALSE)
    }
    fit
  }
  fixed = estimate(prepped, s$lambda_min_kappa_fixed, s$lambda_min_beta_fixed, context)
  between = NULL
  if (s$btwnNetworks){
    means = as.matrix(prepped$data_means[, s$vars, drop = FALSE])
    covariance = stats::cov(means)
    if (nrow(means) <= ncol(means) || any(!is.finite(covariance)) ||
        inherits(try(chol(covariance), silent = TRUE), "try-error")){
      stop(context, ": between-person estimation requires more subjects than variables and a positive-definite covariance of subject means; use 'btwnNetworks = FALSE'.",
           call. = FALSE)
    }
    between = tryCatch(qgraph::EBICglasso(covariance, nrow(means), s$gamma,
                                         returnAllResults = TRUE,
                                         lambda.min.ratio = s$lambda_min_glasso),
                       error = function(e) stop(context, ": ", conditionMessage(e), call. = FALSE))
  }
  individual = vector("list", length(IDs))
  names(individual) = as.character(IDs)
  for (i in seq_along(IDs)){
    individual[[i]] = estimate(
      prepped$data[prepped$data[[s$idvar]] == IDs[i], , drop = FALSE],
      s$lambda_min_kappa, s$lambda_min_beta, paste0("Subject '", IDs[i], "'"),
      list(vars = s$vars, idvar = s$idvar, beepvar = s$beepvar, dayvar = s$dayvar,
           scale = s$scale, centerWithin = s$centerWithin))
  }
  Results = list(fixedPCC = fixed$PCC, fixedPDC = fixed$PDC,
                 fixedbeta = fixed$beta, fixedkappa = fixed$kappa,
                 fixedResults = fixed, betweenNet = between$optnet,
                 betweenResults = between,
                 ids = if (length(IDs)) IDs else unique(data[[s$idvar]]),
                 sample.ids = unique(data[[s$idvar]]), subject.ids = IDs,
                 subjectPCC = lapply(individual, "[[", "PCC"),
                 subjectPDC = lapply(individual, "[[", "PDC"),
                 subjectbeta = lapply(individual, "[[", "beta"),
                 subjectkappa = lapply(individual, "[[", "kappa"),
                 subjecResults = individual, subjectResults = individual)
  class(Results) = "mlGraphicalVAR"
  Results
}

.scgvarSimilarity = function(PCC, PDC, fixedPCC, fixedPDC, elim.groups = FALSE){
  contemporary = lapply(PCC, function(x){
    if (elim.groups) x[abs(fixedPCC) > 0] = 0
    sign(x) * (abs(x) > 0)
  })
  temporal = lapply(PDC, function(x){
    if (elim.groups) x[abs(fixedPDC) > 0] = 0
    sign(x) * (abs(x) > 0.01)
  })
  features = do.call(rbind, lapply(seq_along(PCC), function(i){
    c(as.vector(contemporary[[i]]), as.vector(temporal[[i]]))
  }))
  similarity = tcrossprod(1 * (features > 0)) + tcrossprod(1 * (features < 0))
  diag(similarity) = 0
  similarity
}

.scgvarCluster = function(similarity, search, nclust){
  values = similarity[upper.tri(similarity)]
  if (!any(values > 0)){
    stop("The similarity matrix has no positive between-subject edges; clustering is not possible.",
         call. = FALSE)
  }
  partition = function(m){
    if (!any(m > 0)) return(NULL)
    graph = igraph::graph_from_adjacency_matrix(m, mode = "undirected",
                                               weighted = TRUE, diag = FALSE)
    communities = igraph::cluster_walktrap(graph, steps = 4)
    if (!is.null(nclust) && igraph::components(graph)$no > nclust) return(NULL)
    membership = if (is.null(nclust)) igraph::membership(communities) else
      igraph::cut_at(communities, no = nclust)
    membership = match(membership, unique(membership))
    conductance = vapply(unique(membership), function(group){
      inside = membership == group
      volume = sum(m[inside, , drop = FALSE])
      if (volume == 0) return(1)
      sum(m[inside, !inside, drop = FALSE]) / volume
    }, numeric(1))
    list(membership = membership, score = 1 - mean(conductance),
         modularity = igraph::modularity(graph, membership),
         matrix = m)
  }
  candidates = data.frame(threshold = numeric(), score = numeric())
  threshold = NA_real_
  fallback = FALSE
  if (search){
    # Enumerate the distinct graphs in the quantile threshold interval.
    bounds = stats::quantile(values, c(0.01, 0.99), names = FALSE, type = 7)
    cutoffs = sort(unique(c(bounds[1], values[values >= bounds[1] & values <= bounds[2]])))
    best = NULL
    for (cutoff in cutoffs){
      m = similarity
      m[m <= cutoff] = 0
      candidate = partition(m)
      score = if (is.null(candidate)) NA_real_ else candidate$score
      candidates = rbind(candidates, data.frame(threshold = cutoff, score = score))
      if (!is.null(candidate) && (is.null(best) || score > best$score + 1e-12)){
        best = candidate
        threshold = cutoff
      }
    }
    if (is.null(best)){
      best = partition(similarity)
      fallback = TRUE
      if (!is.null(best)){
        warning("No feasible nonempty thresholded graph; using the original similarity graph.",
                call. = FALSE)
      }
    }
  } else {
    m = similarity
    positive = m > 0
    m[positive] = m[positive] - min(m[positive]) + .Machine$double.eps
    best = partition(m)
  }
  if (is.null(best)){
    stop("'nclust' cannot be attained: it is smaller than the number of disconnected components.",
         call. = FALSE)
  }
  list(matrix = best$matrix, membership = best$membership,
       diagnostics = list(search = search, threshold = threshold,
                           score = best$score, modularity = best$modularity,
                           fallback = fallback, candidates = candidates,
                           nclust = length(unique(best$membership))))
}

.scgvarWrite = function(Results, vars, indplot){
  destination = Results$output.dir
  parent = dirname(destination)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  root = tempfile(".scgvar-", tmpdir = parent)
  if (file.exists(destination) || !dir.create(root, showWarnings = FALSE)){
    stop("Could not create a new output directory: ", destination, call. = FALSE)
  }
  # Stage a complete export before publishing it; failures never leave a
  # partially populated graphical VAR folder or overwrite previous results.
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  writeMatrix = function(x, path){
    utils::write.table(x, file = path, sep = ",", row.names = FALSE, col.names = FALSE)
  }
  plotNetwork = function(pcc, pdc, path){
    grDevices::pdf(path)
    device = grDevices::dev.cur()
    on.exit(grDevices::dev.off(device), add = TRUE)
    qgraph::qgraphMixed(pcc, pdc, parallelEdge = TRUE, theme = "gimme",
                        ltyUndirected = 1, ltyDirected = 2, labels = vars,
                        layout = "circle")
    invisible(NULL)
  }
  writeFixed = function(fit, path, filename){
    if (!dir.create(path, recursive = TRUE, showWarnings = FALSE)){
      stop("Could not create output directory: ", path, call. = FALSE)
    }
    for (name in c("fixedPCC", "fixedPDC", "fixedbeta", "fixedkappa")){
      writeMatrix(fit[[name]], file.path(path, paste0(name, ".csv")))
    }
    if (!is.null(fit$betweenNet)) writeMatrix(fit$betweenNet, file.path(path, "BTW Net.csv"))
    plotNetwork(fit$fixedPCC, fit$fixedPDC, file.path(path, filename))
  }
  writeFixed(Results$group.output, file.path(root, "Group Solution"), "Group Plot.pdf")
  individual = file.path(root, "Individual Solutions")
  dir.create(individual, showWarnings = FALSE)
  fits = Results$group.output$subjectResults
  for (i in seq_along(fits)){
    for (name in c("PCC", "PDC", "beta", "kappa")){
      writeMatrix(fits[[i]][[name]], file.path(individual, paste0("Subj", i, name, ".csv")))
    }
    if (indplot){
      plotNetwork(fits[[i]]$PCC, fits[[i]]$PDC, file.path(individual, paste0("Subj", i, "Plot.pdf")))
    }
  }
  utils::write.csv(data.frame(file.index = seq_along(fits), IDs = Results$membership$IDs),
                   file.path(individual, "Subject Index.csv"), row.names = FALSE)
  subgroups = file.path(root, "Subgroup Solutions")
  dir.create(subgroups, showWarnings = FALSE)
  writeMatrix(Results$membership, file.path(subgroups, "Group Membership.csv"))
  for (i in seq_along(Results$subgroup.outputs)){
    writeFixed(Results$subgroup.outputs[[i]], file.path(subgroups, paste("Subgroup", i)),
               "Subgroup Plot.pdf")
  }
  utils::write.csv(Results$similarity, file.path(root, "sim.mat.test.csv"))
  utils::write.csv(Results$sim.mat, file.path(root, "sim.mat.csv"))
  utils::write.csv(Results$clustering$candidates, file.path(root, "threshold-search.csv"),
                   row.names = FALSE)
  utils::write.csv(as.data.frame(Results$clustering[c("search", "threshold", "score",
                                                    "modularity", "fallback", "nclust")]),
                   file.path(root, "clustering.csv"), row.names = FALSE)
  graph = igraph::graph_from_adjacency_matrix(Results$sim.mat, mode = "undirected",
                                             weighted = TRUE, diag = FALSE)
  plotClusters = function(){
    grDevices::pdf(file.path(root, "Subgroups Plot.pdf"))
    device = grDevices::dev.cur()
    on.exit(grDevices::dev.off(device), add = TRUE)
    graphics::plot(graph, layout = igraph::layout_in_circle(graph),
                    vertex.color = Results$membership$membership,
                    vertex.label = as.character(Results$membership$IDs))
    invisible(NULL)
  }
  plotClusters()
  if (file.exists(destination) || !file.rename(root, destination)){
    stop("Could not finalize output directory: ", destination, call. = FALSE)
  }
  invisible(NULL)
}
