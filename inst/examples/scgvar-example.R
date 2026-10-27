# Reproducible SCGVAR recovery case with no forced subgroup count.
# Run with source() for no file output, or:
# Rscript scgvar-example.R path/to/new-run
library(graphicalVAR)

set.seed(20261026)
vars = paste0("V", seq_len(4))
ids = c("S08", "S02", "S11", "S05", "S01", "S12",
        "S04", "S10", "S03", "S09", "S06", "S07")
truth = data.frame(IDs = ids, true.group = rep(c("positive", "negative"), each = 6),
                   stringsAsFactors = FALSE)
nTime = 320L
burnin = 100L

# beta rows are outcomes; columns are lagged predictors (no intercept here).
beta.positive = diag(0.35, length(vars))
beta.positive[2, 1] = 0.8
beta.positive[4, 3] = 0.8
beta.negative = beta.positive
beta.negative[2, 1] = -0.8
beta.negative[4, 3] = -0.8
kappa = diag(length(vars))
stopifnot(max(Mod(eigen(beta.positive)$values)) < 1,
          max(Mod(eigen(beta.negative)$values)) < 1)

simulateSubject = function(i){
  beta = if (truth$true.group[i] == "positive") beta.positive else beta.negative
  observations = matrix(0, nTime + burnin, length(vars))
  for (time in seq.int(2L, nrow(observations))){
    observations[time, ] = as.numeric(beta %*% observations[time - 1L, ]) +
      stats::rnorm(length(vars))
  }
  observations = observations[seq.int(burnin + 1L, nrow(observations)), , drop = FALSE]
  colnames(observations) = vars
  data.frame(id = ids[i], day = 1L, beep = seq_len(nTime), observations,
             check.names = FALSE)
}
data = do.call(rbind, lapply(seq_along(ids), simulateSubject))
row.names(data) = NULL

# Only direct Rscript invocation consumes the optional output-path argument.
arguments = commandArgs(trailingOnly = TRUE)
scriptMode = sys.nframe() == 0L && any(grepl("^--file=", commandArgs()))
outputDir = if (scriptMode && length(arguments)) arguments[1L] else NULL
if (scriptMode && length(arguments) > 1L){
  stop("Supply at most one argument: a fresh output directory.", call. = FALSE)
}

# Truth labels are separate and are not supplied to the estimator.
# nLambda = 6 is a short integration-test grid, not a tuning recommendation.
elapsed = system.time({
  fit = scgvar(data, vars = vars, idvar = "id", beepvar = "beep", dayvar = "day",
                search = TRUE, nclust = NULL, gamma = 0.5, nLambda = 6,
                lambda_min_kappa = 0.03, lambda_min_beta = 0.03,
                verbose = FALSE, dir = outputDir)
})
comparison = fit$membership
comparison$true.group = truth$true.group[match(comparison$IDs, truth$IDs)]
stopifnot(!anyNA(comparison), identical(comparison$IDs, ids),
          all(is.finite(fit$similarity)), all(is.finite(fit$sim.mat)),
          all(diag(fit$sim.mat) == 0))
adjustedRand = igraph::compare(comparison$membership,
                               as.integer(factor(comparison$true.group)),
                               method = "adjusted.rand")
cat("Detected subgroups:", fit$clustering$nclust, "\n")
cat("Adjusted Rand index:", format(adjustedRand, digits = 6), "\n")
cat("Threshold fallback:", fit$clustering$fallback, "\n")
cat("Elapsed seconds:", elapsed[["elapsed"]], "\n\n")
print(comparison, row.names = FALSE)
cat("\nIndependent truth versus recovered membership:\n")
print(table(Truth = comparison$true.group, Cluster = comparison$membership))

coefficients = do.call(rbind, lapply(seq_along(fit$subgroup.outputs), function(i){
  nativeFit = fit$subgroup.outputs[[i]]$fixedResults
  stopifnot(inherits(nativeFit, "graphicalVAR"),
            all(is.finite(nativeFit$beta)), all(is.finite(nativeFit$kappa)))
  data.frame(cluster = i,
             V1.to.V2.beta = nativeFit$beta["V2", "V1_lag1"],
             V3.to.V4.beta = nativeFit$beta["V4", "V3_lag1"],
             V1.to.V2.PDC = nativeFit$PDC["V1", "V2"],
             V3.to.V4.PDC = nativeFit$PDC["V3", "V4"])
}))
cat("\nSubgroup coefficients and directed partial correlations:\n")
print(coefficients, row.names = FALSE)
cat("\nRecovery is measured without a forced cluster count.\n")
cat("This generated case is an implementation check, not a general accuracy guarantee.\n")

if (!is.null(outputDir)){
  saveRDS(list(data = data, truth = truth, beta.positive = beta.positive,
                beta.negative = beta.negative, kappa = kappa),
           file.path(outputDir, "simulation-input.rds"))
  saveRDS(fit, file.path(outputDir, "scgvar-fit.rds"))
  utils::write.csv(comparison, file.path(outputDir, "membership-comparison.csv"),
                   row.names = FALSE)
  utils::write.csv(coefficients, file.path(outputDir, "subgroup-coefficients.csv"),
                   row.names = FALSE)
  summaryText = c(paste("Detected subgroups:", fit$clustering$nclust),
                  paste("Adjusted Rand index:", adjustedRand),
                  paste("Threshold fallback:", fit$clustering$fallback),
                  paste("Elapsed seconds:", elapsed[["elapsed"]]),
                  "No cluster count was supplied to estimation.",
                  capture.output(sessionInfo()))
  writeLines(summaryText, file.path(outputDir, "example-results.txt"))
}
