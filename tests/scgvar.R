# Regression and integration tests for the SCGVAR extension.
# Run with R CMD check, or source this file after installing the package.
library(graphicalVAR)

local({
  expectError = function(expr, pattern = NULL) {
    result = tryCatch(force(expr), error = identity)
    stopifnot(inherits(result, "error"))
    if (!is.null(pattern)) {
      stopifnot(grepl(pattern, conditionMessage(result), ignore.case = TRUE))
    }
    invisible(result)
  }
  expectEqual = function(actual, expected, tolerance = 1e-8) {
    stopifnot(isTRUE(all.equal(actual, expected, tolerance = tolerance)))
  }
  similarity = getFromNamespace(".scgvarSimilarity", "graphicalVAR")
  cluster = getFromNamespace(".scgvarCluster", "graphicalVAR")
  writeResults = getFromNamespace(".scgvarWrite", "graphicalVAR")

  # Independently counted edges: PCC uses both triangles; PDC includes
  # autoregressions and excludes edges whose absolute value is exactly 0.01.
  pccA = pccB = pccC = matrix(0, 3, 3)
  pccA[1, 2] = pccA[2, 1] = 0.4
  pccA[1, 3] = pccA[3, 1] = -0.3
  pccB[1, 2] = pccB[2, 1] = 0.2
  pccB[2, 3] = pccB[3, 2] = 0.6
  pccC[1, 2] = pccC[2, 1] = -0.2
  pccC[1, 3] = pccC[3, 1] = -0.2
  pdcA = pdcB = pdcC = diag(0.2, 3)
  pdcA[1, 2] = pdcB[1, 2] = 0.4
  pdcC[1, 2] = -0.4
  pdcA[2, 3] = pdcB[2, 3] = pdcC[2, 3] = 0.01
  fixedPCC = matrix(0, 3, 3)
  fixedPCC[1, 2] = fixedPCC[2, 1] = 0.2
  fixedPDC = diag(0.2, 3)
  rawExpected = matrix(c(0, 6, 5, 6, 0, 3, 5, 3, 0), 3, 3)
  removedExpected = matrix(c(0, 1, 2, 1, 0, 0, 2, 0, 0), 3, 3)
  expectEqual(unname(similarity(list(pccA, pccB, pccC),
                               list(pdcA, pdcB, pdcC),
                               fixedPCC, fixedPDC)), rawExpected)
  expectEqual(unname(similarity(list(pccA, pccB, pccC),
                               list(pdcA, pdcB, pdcC),
                               fixedPCC, fixedPDC, elim.groups = TRUE)),
              removedExpected)

  # Community recovery, singleton handling, and an unidentifiable graph.
  block = matrix(1, 6, 6)
  block[1:3, 1:3] = block[4:6, 4:6] = 8
  diag(block) = 0
  blockResult = cluster(block, search = TRUE, nclust = 2)
  blockMembership = blockResult$membership
  stopifnot(length(unique(blockMembership[1:3])) == 1L,
            length(unique(blockMembership[4:6])) == 1L,
            blockMembership[1] != blockMembership[4])
  stopifnot(length(unique(cluster(block, search = TRUE,
                                  nclust = NULL)$membership)) == 2L)
  disconnected = block
  disconnected[1:3, 4:6] = disconnected[4:6, 1:3] = 0
  stopifnot(length(unique(cluster(disconnected, search = FALSE,
                                  nclust = 2)$membership)) == 2L)
  expectError(cluster(disconnected, search = FALSE, nclust = 1))
  expectError(cluster(matrix(0, 3, 3), search = TRUE, nclust = NULL))
  isolated = matrix(0, 7, 7)
  isolated[1:6, 1:6] = disconnected
  stopifnot(length(unique(cluster(isolated, search = FALSE,
                                  nclust = 3)$membership)) == 3L)
  equalWeights = matrix(1, 3, 3)
  diag(equalWeights) = 0
  warningState = new.env(parent = emptyenv())
  warningState$text = character()
  fallback = withCallingHandlers(cluster(equalWeights, search = TRUE, nclust = 1),
    warning = function(w) {
      warningState$text = c(warningState$text, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  stopifnot(fallback$diagnostics$fallback, length(warningState$text) == 1L,
            length(unique(fallback$membership)) == 1L)

  # A reproducible two-group VAR: triangular beta ensures stationarity, and
  # the groups have opposite signs on two cross-lagged effects.
  set.seed(20261026)
  vars = c("x", "y", "z")
  ids = c("S20", "S3", "S91", "S7", "S40", "S2", "S82", "S6")
  trueGroup = rep(1:2, each = 4)
  beta = diag(0.15, 3)
  beta[2, 1] = beta[3, 2] = 0.65
  kappa = diag(3)
  kappa[1, 2] = kappa[2, 1] = -0.25
  records = lapply(seq_along(ids), function(i) {
    subjectBeta = beta
    if (trueGroup[i] == 2L) {
      subjectBeta[2, 1] = subjectBeta[3, 2] = -0.65
    }
    values = as.data.frame(graphicalVARsim(180, subjectBeta, kappa))
    names(values) = vars
    values$id = ids[i]
    values$day = 1L
    values$beep = seq_len(nrow(values))
    values
  })
  dat = do.call(rbind, records)
  baseArgs = list(data = dat, vars = vars, idvar = "id", dayvar = "day",
                  beepvar = "beep", search = FALSE, verbose = FALSE,
                  lambda_beta = 0.08, lambda_kappa = 0.10, gamma = 0,
                  likelihood = "penalized", ebic_tol = 2e-4, nclust = 2)
  run = function(...) {
    changes = list(...)
    args = baseArgs
    for (name in names(changes)) args[name] = changes[name]
    do.call(scgvar, args)
  }

  # Invalid data/settings must fail before creating output or fitting models.
  expectError(run(vars = "x"))
  expectError(run(vars = c("x", "x")))
  expectError(run(vars = c("x", "absent")))
  expectError(run(idvar = "absent"))
  expectError(run(search = NA))
  expectError(run(nclust = 1.5))
  expectError(run(nclust = 9))
  expectError(run(subjectNetworks = FALSE))
  expectError(run(subjectNetworks = "S20"))
  expectError(run(subjectNetworks = c("S20", "absent")))
  expectError(run(subjectNetworks = c("S20", "S20")))
  expectError(run(lags = 2))
  expectError(run(lags = c(1, 2)))
  expectError(run(not_an_argument = 1))
  expectError(run(redo = TRUE))
  expectError(run(mlgvclust = FALSE))
  expectError(run(g.ID = "id"))
  expectError(run(maxit.in = 0))
  expectError(run(gamma = Inf))
  expectError(run(lambda_min_kappa_fixed = 0))
  expectError(run(lambda_min_kappa = 1))
  expectError(run(indplot = TRUE))
  expectError(run(returnResults = FALSE))
  expectError(run(nLambda = 0))
  expectError(run(maxit.out = 1.5))
  expectError(run(lambda_beta = NA_real_))
  expectError(run(lambda_beta = c(0.1, 0.1)))
  expectError(run(lambda_kappa = -0.1))
  expectError(run(ebic_tol = 0))
  expectError(run(penalize.diagonal = NA))
  expectError(run(regularize_mat_beta = matrix(TRUE, 2, 2)))
  expectError(run(regularize_mat_kappa = matrix(TRUE, 2, 2)))
  asymmetricMask = matrix(TRUE, 3, 3)
  asymmetricMask[1, 2] = FALSE
  expectError(run(regularize_mat_kappa = asymmetricMask))
  bad = dat
  bad$x[1] = NA_real_
  expectError(run(data = bad))
  bad$x[1] = Inf
  expectError(run(data = bad))
  bad = dat
  bad$x = as.character(bad$x)
  expectError(run(data = bad))
  bad = dat
  bad$id[1] = NA_character_
  expectError(run(data = bad))
  bad = dat
  bad$beep[2] = bad$beep[1]
  expectError(run(data = bad))
  bad = dat
  bad$beep[1] = 1.5
  expectError(run(data = bad))
  bad = dat
  bad$beep = bad$beep * 2L
  expectError(run(data = bad), "lag-1 pairs")
  bad = dat
  bad$x[bad$id == ids[1]] = 1
  expectError(run(data = bad))
  expectError(run(data = dat[dat$id == ids[1], , drop = FALSE]))

  # Calling the function without dir has no file, working-directory, or
  # graphics-device effects. Use an explicit test device so par() never
  # implicitly opens Rplots.pdf.
  deviceFile = tempfile(fileext = ".pdf")
  grDevices::pdf(deviceFile)
  testDevice = grDevices::dev.cur()
  on.exit(if (testDevice %in% grDevices::dev.list()) {
    grDevices::dev.off(testDevice)
  }, add = TRUE)
  oldDirectory = getwd()
  oldDevices = grDevices::dev.list()
  oldPar = graphics::par(no.readonly = TRUE)
  oldFiles = list.files(all.files = TRUE)
  resultWithVisibility = withVisible(run())
  result = resultWithVisibility$value
  stopifnot(!resultWithVisibility$visible, inherits(result, "scgvar"),
            identical(getwd(), oldDirectory),
            identical(grDevices::dev.list(), oldDevices),
            identical(list.files(all.files = TRUE), oldFiles),
            identical(graphics::par(no.readonly = TRUE), oldPar),
            is.null(result$output.dir), is.data.frame(result$membership),
            identical(result$membership$IDs, ids),
            identical(result$group.output$sample.ids, ids),
            identical(result$group.output$subject.ids, ids),
            identical(names(result$group.output$subjectPCC), ids),
            identical(result$group.output$subjectResults,
                      result$group.output$subjecResults))
  inferred = result$membership$membership
  stopifnot(length(unique(inferred[trueGroup == 1])) == 1L,
            length(unique(inferred[trueGroup == 2])) == 1L,
            inferred[1] != inferred[5],
            length(result$subgroup.outputs) == 2L)
  stopifnot(all(diag(result$similarity) == 0),
            isTRUE(all.equal(result$similarity, t(result$similarity))),
            all(result$sim.mat >= 0))

  # Fixed and individual estimators agree numerically with direct calls to
  # the current public graphicalVAR API, including modern optional arguments.
  direct = graphicalVAR(dat, vars = vars, idvar = "id", dayvar = "day",
                        beepvar = "beep", verbose = FALSE, gamma = 0,
                        lambda_beta = 0.08, lambda_kappa = 0.10,
                        likelihood = "penalized", ebic_tol = 2e-4)
  for (field in c("beta", "kappa", "PCC", "PDC", "EBIC", "N")) {
    expectEqual(result$group.output$fixedResults[[field]], direct[[field]])
  }
  for (i in seq_along(ids)) {
    directSubject = graphicalVAR(dat[dat$id == ids[i], , drop = FALSE],
                                 vars = vars, idvar = "id", dayvar = "day",
                                 beepvar = "beep", verbose = FALSE, gamma = 0,
                                 lambda_beta = 0.08, lambda_kappa = 0.10,
                                 likelihood = "penalized", ebic_tol = 2e-4)
    for (field in c("beta", "kappa", "PCC", "PDC", "EBIC", "N")) {
      expectEqual(result$group.output$subjectResults[[i]][[field]],
                  directSubject[[field]])
    }
  }
  for (group in seq_along(result$subgroup.outputs)) {
    selected = result$membership$IDs[inferred == group]
    stopifnot(identical(result$subgroup.outputs[[group]]$ids, selected),
              identical(result$subgroup.outputs[[group]]$sample.ids, selected),
              length(result$subgroup.outputs[[group]]$subject.ids) == 0L)
    directGroup = graphicalVAR(dat[dat$id %in% selected, , drop = FALSE],
                               vars = vars, idvar = "id", dayvar = "day",
                               beepvar = "beep", verbose = FALSE, gamma = 0,
                               lambda_beta = 0.08, lambda_kappa = 0.10,
                               likelihood = "penalized", ebic_tol = 2e-4)
    for (field in c("beta", "kappa", "PCC", "PDC", "EBIC", "N")) {
      expectEqual(result$subgroup.outputs[[group]]$fixedResults[[field]],
                  directGroup[[field]])
    }
  }

  # Current regularization masks (including the implicit intercept column)
  # reach the public estimator unchanged.
  betaMask = matrix(TRUE, 3, 3)
  diag(betaMask) = FALSE
  maskResult = run(regularize_mat_beta = betaMask)
  directMask = graphicalVAR(dat, vars = vars, idvar = "id", dayvar = "day",
                            beepvar = "beep", verbose = FALSE, gamma = 0,
                            lambda_beta = 0.08, lambda_kappa = 0.10,
                            regularize_mat_beta = betaMask,
                            likelihood = "penalized", ebic_tol = 2e-4)
  expectEqual(maskResult$group.output$fixedbeta, directMask$beta)

  unscaled = run(scale = FALSE, centerWithin = FALSE)
  directUnscaled = graphicalVAR(dat, vars = vars, idvar = "id", dayvar = "day",
                                beepvar = "beep", verbose = FALSE, gamma = 0,
                                scale = FALSE, centerWithin = FALSE,
                                lambda_beta = 0.08, lambda_kappa = 0.10,
                                likelihood = "penalized", ebic_tol = 2e-4)
  expectEqual(unscaled$group.output$fixedbeta, directUnscaled$beta)
  expectEqual(unscaled$group.output$fixedkappa, directUnscaled$kappa)

  # Separate fixed and individual tuning ranges also apply to subgroup fits.
  gridArgs = baseArgs
  gridArgs$lambda_beta = gridArgs$lambda_kappa = NULL
  gridArgs$nLambda = 3
  gridArgs$lambda_min_beta_fixed = 0.2
  gridArgs$lambda_min_kappa_fixed = 0.3
  gridArgs$lambda_min_beta = 0.4
  gridArgs$lambda_min_kappa = 0.5
  gridResult = do.call(scgvar, gridArgs)
  checkRange = function(fit, betaRatio, kappaRatio) {
    stopifnot(nrow(fit$path) == 9L)
    expectEqual(min(fit$path$beta) / max(fit$path$beta), betaRatio)
    expectEqual(min(fit$path$kappa) / max(fit$path$kappa), kappaRatio)
  }
  checkRange(gridResult$group.output$fixedResults, 0.2, 0.3)
  for (fit in gridResult$group.output$subjectResults) checkRange(fit, 0.4, 0.5)
  for (fit in gridResult$subgroup.outputs) checkRange(fit$fixedResults, 0.2, 0.3)

  # A single singular person must be reported by ID, never silently converted
  # into an empty network or an exported partial analysis.
  singular = dat
  firstRows = singular$id == ids[1]
  singular$y[firstRows] = singular$x[firstRows]
  failedExport = tempfile("scgvar-failed-")
  expectError(run(data = singular, lambda_beta = 0, lambda_kappa = 0,
                  dir = failedExport), "Subject.*S20")
  stopifnot(!dir.exists(file.path(failedExport, "graphical VAR")))

  # Between-person estimation is optional and uses current EBICglasso.
  # Tiny subgroup samples can produce EBICglasso's dense-network warning;
  # that warning is expected in this numerical compatibility fixture.
  betweenResult = suppressWarnings(run(btwnNetworks = TRUE))
  means = as.matrix(betweenResult$group.output$fixedResults$data$data_means[, vars])
  directBetween = suppressWarnings(qgraph::EBICglasso(
    stats::cov(means), nrow(means), gamma = 0, returnAllResults = TRUE,
    lambda.min.ratio = 0.01))
  expectEqual(betweenResult$group.output$betweenNet, directBetween$optnet)
  stopifnot(is.null(result$group.output$betweenNet),
            all(vapply(betweenResult$subgroup.outputs,
                        function(x) is.matrix(x$betweenNet), logical(1))))

  # Subsetting must map actual IDs, not their positions or sorted row order.
  selectedIDs = ids[c(8, 2, 5, 3)]
  selected = run(subjectNetworks = selectedIDs)
  stopifnot(identical(selected$membership$IDs, selectedIDs),
            identical(names(selected$group.output$subjectPCC), selectedIDs),
            identical(selected$group.output$sample.ids, ids),
            identical(selected$group.output$subject.ids, selectedIDs))
  expectEqual(selected$group.output$fixedPCC, result$group.output$fixedPCC)
  for (id in selectedIDs) {
    expectEqual(selected$group.output$subjectPCC[[id]],
                result$group.output$subjectPCC[[id]])
  }
  numericData = dat
  numericData$id = c(20L, 3L, 91L, 7L, 40L, 2L, 82L, 6L)[match(dat$id, ids)]
  numericResult = run(data = numericData, subjectNetworks = c(82L, 3L),
                      nclust = 1)
  stopifnot(identical(numericResult$membership$IDs, c(82L, 3L)))
  factorData = dat
  factorData$id = factor(factorData$id, levels = rev(ids))
  factorResult = run(data = factorData, subjectNetworks = ids[c(8, 2)], nclust = 1)
  stopifnot(is.factor(factorResult$membership$IDs),
            identical(as.character(factorResult$membership$IDs), ids[c(8, 2)]))
  suppressWarnings(expectError(run(subjectNetworks = ids[c(8, 2)],
                                    btwnNetworks = TRUE, nclust = 1),
                                "Subgroup.*between-person"))

  # Omitted timing uses each subject's row order. Explicit timing sorts rows
  # without creating lagged links across days or absent beeps.
  implicitTime = run(beepvar = NULL, dayvar = NULL)
  expectEqual(implicitTime$group.output$fixedPCC, result$group.output$fixedPCC)
  expectEqual(implicitTime$group.output$fixedPDC, result$group.output$fixedPDC)
  # Old tsData default column names must not overwrite measured variables.
  reservedNames = dat
  names(reservedNames)[match(vars, names(reservedNames))] = c("ID", "DAY", "BEEP")
  reservedResult = run(data = reservedNames, vars = c("ID", "DAY", "BEEP"),
                       beepvar = NULL, dayvar = NULL)
  expectEqual(unname(reservedResult$group.output$fixedPCC),
              unname(result$group.output$fixedPCC))
  timed = dat
  timed$day = ifelse(timed$beep <= 90, 1L, 2L)
  timed$beep = (timed$beep - 1L) %% 90L + 1L
  timed = timed[timed$beep != 45L, , drop = FALSE]
  timed = timed[rev(seq_len(nrow(timed))), , drop = FALSE]
  timeResult = run(data = timed)
  # Per subject/day: 89 observations, two contiguous sequences -> 87 pairs.
  stopifnot(timeResult$group.output$fixedResults$N == 8L * 2L * 87L)

  # Explicit export completes only after estimation and leaves caller state
  # intact. A second call must refuse to overwrite its existing output.
  exportParent = tempfile("scgvar-test-")
  exportResult = run(dir = exportParent, indplot = TRUE)
  stopifnot(dir.exists(exportResult$output.dir),
            file.exists(file.path(exportResult$output.dir, "Group Solution",
                                  "fixedPCC.csv")),
            file.exists(file.path(exportResult$output.dir, "Subgroup Solutions",
                                  "Group Membership.csv")),
            identical(getwd(), oldDirectory),
            identical(grDevices::dev.list(), oldDevices),
            identical(graphics::par(no.readonly = TRUE), oldPar))
  exportFiles = list.files(exportResult$output.dir, recursive = TRUE)
  individualPlots = file.path(exportResult$output.dir, "Individual Solutions",
                               paste0("Subj", seq_along(ids), "Plot.pdf"))
  stopifnot(any(grepl("[.]pdf$", exportFiles, ignore.case = TRUE)),
            all(file.exists(individualPlots)),
            all(file.info(individualPlots)$size > 0))
  expectError(run(dir = exportParent))
  # Inject an export error after opening a plotting device. The staging
  # directory and device must be cleaned, without publishing partial files.
  invalidExport = result
  invalidParent = tempfile("scgvar-export-failure-")
  invalidExport$output.dir = file.path(invalidParent, "graphical VAR")
  invalidExport$group.output$fixedPCC = matrix("not numeric", 3, 3)
  suppressWarnings(expectError(writeResults(invalidExport, vars, indplot = FALSE)))
  remaining = setdiff(list.files(invalidParent, all.files = TRUE), c(".", ".."))
  stopifnot(!dir.exists(invalidExport$output.dir), length(remaining) == 0L,
            identical(grDevices::dev.list(), oldDevices),
            identical(graphics::par(no.readonly = TRUE), oldPar),
            identical(getwd(), oldDirectory))
  silentResult = run(returnResults = FALSE, dir = tempfile("scgvar-silent-"))
  stopifnot(is.null(silentResult))
  cat("SCGVAR integration and regression tests passed.\n")
})
