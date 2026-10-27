\name{scgvar}
\alias{scgvar}
\title{
Subgroup chain graphical VAR estimation
}
\description{
Estimates pooled and individual graphical VAR networks, clusters subjects by
signed network similarity, and refits a pooled graphical VAR in each subgroup.
Estimation uses the current \code{\link{graphicalVAR}} implementation and its
LASSO and EBIC controls. Files are written only when \code{dir} is supplied.
}
\usage{
scgvar(data, vars, idvar, beepvar = NULL, dayvar = NULL,
       search = TRUE, scale = TRUE, maxit.in = 100,
       centerWithin = TRUE, gamma = 0.5, verbose = TRUE,
       redo = FALSE, subjectNetworks = TRUE,
       lambda_min_kappa_fixed = 0.001, lambda_min_beta_fixed = 0.001,
       lambda_min_kappa = 0.05, lambda_min_beta = lambda_min_kappa,
       lambda_min_glasso = 0.01, dir = NULL, btwnNetworks = FALSE,
       mlgvclust = TRUE, indplot = FALSE, elim.groups = FALSE,
       g.ID = NULL, nclust = NULL, returnResults = TRUE, ...)
}
\arguments{
  \item{data}{Matrix or data frame of repeated measures, with observations in
  rows. Selected columns must not contain missing values. Measurement columns
  must be finite and numeric. Every subject must have nonconstant measurements
  and at least three usable lag-1 pairs within days. This minimum is an input
  check, not a guarantee of reliable estimation.}
  \item{vars}{Character vector of at least two distinct measurement column
  names, different from the ID and time columns.}
  \item{idvar}{Name of the subject ID column. IDs must be nonmissing atomic
  values. At least two subjects must be selected for clustering.}
  \item{beepvar}{Name of a column of nonnegative integer assessment numbers
  within each subject and day, or \code{NULL}. Gaps prevent lag pairs across
  missing assessments. Duplicate subject/day/beep combinations are rejected.
  If \code{NULL}, the observed row order within each subject and day defines
  consecutive assessments.}
  \item{dayvar}{Name of the day column, or \code{NULL} to treat each subject's
  series as one day. Lag pairs never cross day boundaries. Supply a day column
  only when multiple observations within a day are available.}
  \item{search}{Logical. If \code{TRUE}, enumerate feasible thresholded
  similarity graphs and select the highest conductance-based score. If
  \code{FALSE}, subtract the smallest positive similarity from every positive
  entry and add machine epsilon. See Details.}
  \item{scale}{Logical. Standardize measurements during estimation, following
  the preparation used by \code{\link{mlGraphicalVAR}}.}
  \item{maxit.in}{Positive integer giving the maximum inner iterations for
  every \code{\link{graphicalVAR}} fit.}
  \item{centerWithin}{Logical. Center measurements within subjects before
  estimating pooled fixed effects.}
  \item{gamma}{Finite nonnegative EBIC hyperparameter. Use zero for BIC.}
  \item{verbose}{Logical. Print progress messages for the SCGVAR stages.
  Diagnostic warnings and estimator messages can still be emitted.}
  \item{redo}{Legacy experimental argument. Only \code{FALSE} is supported.}
  \item{subjectNetworks}{\code{TRUE} to cluster all subjects, or a vector of
  at least two distinct observed subject IDs. The supplied order is preserved.
  Pooled fixed effects still use all supplied data; individual fits, clustering,
  and subgroup refits use only the selected subjects. \code{FALSE} is invalid.}
  \item{lambda_min_kappa_fixed}{Minimum-to-maximum tuning parameter ratio
  for kappa in pooled and subgroup fixed-effect fits, in \eqn{(0, 1)}.}
  \item{lambda_min_beta_fixed}{Minimum-to-maximum tuning parameter ratio
  for beta in pooled and subgroup fixed-effect fits, in \eqn{(0, 1)}.}
  \item{lambda_min_kappa}{Minimum-to-maximum tuning parameter ratio for
  kappa in individual fits, in \eqn{(0, 1)}.}
  \item{lambda_min_beta}{Minimum-to-maximum tuning parameter ratio for
  beta in individual fits, in \eqn{(0, 1)}.}
  \item{lambda_min_glasso}{Minimum-to-maximum glasso tuning parameter ratio
  for optional between-subject networks, in \eqn{(0, 1)}.}
  \item{dir}{\code{NULL} for no file output, or a parent directory under which
  a new directory named \file{graphical VAR} will be created after estimation.
  An existing \file{graphical VAR} directory is never overwritten. Supplying
  \code{dir} exports CSV files and pooled, subgroup, and membership PDF plots.}
  \item{btwnNetworks}{Logical. Estimate between-subject networks using
  \code{\link[qgraph]{EBICglasso}} for the full sample and each subgroup.
  Each such fit requires more subjects than measurement variables and a
  positive-definite covariance matrix of subject means; otherwise an error
  identifies the affected fit. Defaults to \code{FALSE}.}
  \item{mlgvclust}{Legacy experimental argument. Only \code{TRUE} is supported.}
  \item{indplot}{Logical. Also export an individual-network PDF for each
  selected subject. Requires a non-\code{NULL} \code{dir}.}
  \item{elim.groups}{Logical. Before computing similarities, set individual
  PCC and PDC entries to zero wherever the corresponding pooled fixed network
  has a nonzero entry. This removes shared fixed-edge positions regardless of
  sign and corrects the reversed mask in the archived implementation.}
  \item{g.ID}{Legacy experimental argument. Only \code{NULL} is supported.}
  \item{nclust}{\code{NULL} for the walktrap-selected number of communities,
  or an integer from one to the number of selected subjects specifying a
  dendrogram cut. A cut cannot have fewer clusters than the graph has connected
  components. Infeasible candidate graphs are excluded from the search.}
  \item{returnResults}{Logical. Invisibly return the structured result when
  \code{TRUE}. If \code{FALSE}, \code{dir} is required and the function
  invisibly returns \code{NULL} after export.}
  \item{\dots}{Additional uniquely and exactly named arguments to
  \code{\link{graphicalVAR}}, such as \code{nLambda}, \code{lambda_beta},
  \code{lambda_kappa}, \code{regularize_mat_beta},
  \code{regularize_mat_kappa}, \code{maxit.out}, and \code{ebic_tol}.
  These apply to pooled, individual, and subgroup fits. \code{nLambda} must
  be one positive integer. Unknown, unnamed,
  duplicate, and abbreviated argument names are rejected. Only
  \code{lags = 1} is supported. Explicit lambda vectors must contain distinct,
  finite, nonnegative values; they override generated grids and their minimum
  ratios as in \code{graphicalVAR}.}
}
\details{
The VAR(1) model is
\deqn{y_t = a + B y_{t-1} + \epsilon_t, \qquad
       \epsilon_t \sim N(0, K^{-1}).}
The pooled and individual estimation steps follow the data preparation used by
\code{\link{mlGraphicalVAR}}. Individual data are taken from the pooled
prepared data and prepared again for subject-specific estimation. Subgroup
models are refitted from the original rows belonging to each subgroup.

Let \eqn{C_i} and \eqn{D_i} denote subject \eqn{i}'s PCC and PDC matrices,
after optional fixed-edge removal. Subject similarity is
\deqn{S_{ij} = \sum_{a,b} I(|C_{i,ab}| > 0, |C_{j,ab}| > 0,
                          \mathrm{sign}(C_{i,ab}) = \mathrm{sign}(C_{j,ab}))
             + \sum_{a,b} I(|D_{i,ab}| > 0.01, |D_{j,ab}| > 0.01,
                          \mathrm{sign}(D_{i,ab}) = \mathrm{sign}(D_{j,ab})).}
The diagonal of \eqn{S} is zero. The full PCC matrices count each undirected
off-diagonal edge twice. The full PDC matrices count each directed edge once,
including autoregressive diagonal entries. Shared zeros do not count.

For \code{search = TRUE}, let \eqn{q_{.01}} and \eqn{q_{.99}} be the type-7
empirical quantiles of the upper triangle of \eqn{S}. Every distinct graph
obtainable by setting entries less than or equal to a threshold in
\eqn{[q_{.01}, q_{.99}]} to zero is evaluated. Weighted walktrap uses four
steps. For a partition into \eqn{k} communities, the score is
\deqn{1 - \frac{1}{k}\sum_{g=1}^{k}
  \frac{\sum_{i \in g, j \notin g} S^*_{ij}}
       {\sum_{i \in g, j} S^*_{ij}},}
where \eqn{S^*} is the candidate weighted adjacency matrix. A community with
zero total degree is assigned an escape fraction of one. Scores within
\eqn{10^{-12}} are tied and the lower threshold wins. If \code{nclust} is
given, both the score and reported modularity use that cut. The search
enumerates discrete graphs instead of using the archive's continuous optimizer;
it does not guarantee the archive's selected partition.

Empty or infeasible graphs are excluded. If no thresholded graph is feasible,
the original similarity graph is used with a warning when it permits the
requested cut; \code{clustering$fallback} records this choice. Otherwise the
function errors. No edges are added to force a requested cluster count. An
all-zero similarity matrix is an error even when \code{nclust} is supplied.
The conductance-based criterion can legitimately prefer a single cluster.

For \code{search = FALSE}, zero weights remain zero and positive weights become
\eqn{S_{ij} - \min(S_{uv}: S_{uv}>0) + \epsilon}, where \eqn{\epsilon} is
\code{.Machine$double.eps}. Walktrap is applied to this weighted graph.

Pooled, individual, and subgroup estimates are native \code{graphicalVAR}
objects accessible through the returned fit components. Their beta matrices
have outcome variables in rows and an intercept column named \code{"1"}
followed by columns named \code{paste0(vars, "_lag1")}. PDC matrices use the
network convention: rows are origins and columns are destinations. Thus
\code{beta["V2", "V1_lag1"]} and \code{PDC["V1", "V2"]} describe the same
direction. Use the existing \code{plot} method on a \code{fixedResults} or
\code{subjectResults} object.

No graphics device or output file is opened by default. Explicit export uses
separate PDF devices and closes them after plotting, preserving the caller's
device and plotting state. Estimation errors occur before export begins.
Files are first written to a temporary sibling directory, which is renamed to
the final output directory only after every export succeeds. On an export
error, temporary files are removed and opened PDF devices are closed.
}
\value{
Invisibly returns a list of class \code{"scgvar"}, unless
\code{returnResults = FALSE}. Its components are:
\item{group.output}{Pooled and selected individual estimates in an
\code{mlGraphicalVAR}-compatible list. Includes \code{fixedResults},
\code{fixedPCC}, \code{fixedPDC}, \code{fixedbeta}, \code{fixedkappa},
\code{betweenNet}, \code{betweenResults}, \code{ids}, individual matrix lists,
and \code{subjectResults}. \code{sample.ids} identifies subjects contributing
to pooled estimation, and \code{subject.ids} identifies fitted individuals.
The historical spelling \code{subjecResults} is
also available. Individual result lists are named by subject ID.}
\item{similarity}{Raw signed-edge similarity matrix, with subject ID dimnames.}
\item{sim.mat}{The weighted adjacency matrix actually used for clustering.}
\item{membership}{Data frame with \code{IDs} and integer \code{membership}
columns. Subject ID types are preserved. Community labels are consecutive
integers in first-occurrence order.}
\item{subgroup.outputs}{Named list of subgroup pooled fits, with the same
fixed-effect components as \code{group.output}. Individual result lists are
empty because subjects are not refitted individually within each subgroup.
Here \code{ids} and \code{sample.ids} contain subgroup subjects, and
\code{subject.ids} is empty.}
\item{clustering}{Diagnostics: \code{search}, selected \code{threshold},
\code{score}, \code{modularity}, \code{fallback}, candidate threshold/score
table \code{candidates}, and detected \code{nclust}. The threshold is
\code{NA} when no threshold was selected; infeasible candidates have
\code{NA} scores.}
\item{output.dir}{\code{NULL}, or the normalized export directory.}
\item{call}{The matched function call.}
}
\note{
The installed \file{doc/scgvar-guide.txt} describes migration from the local
SCGVAR archive, exports, and interpretation. The installed
\file{examples/scgvar-example.R} contains a reproducible two-group simulation
with independent truth labels and adjusted Rand index reporting. A successful
simulation is an implementation check, not evidence of accuracy for every
sample size or data-generating process.
}
\author{
Jonathan J. Park and Sacha Epskamp
}
\seealso{
\code{\link{graphicalVAR}}, \code{\link{mlGraphicalVAR}},
\code{\link{graphicalVARsim}}, \code{\link[igraph]{cluster_walktrap}}
}
\examples{
\dontrun{
# Run the reproducible recovery example without exporting files:
source(system.file("examples", "scgvar-example.R", package = "graphicalVAR"))

# Inspect and plot a native pooled graphicalVAR fit:
fit$membership
plot(fit$group.output$fixedResults)

# Select lag coefficients, excluding the intercept:
fit$subgroup.outputs[[1]]$fixedResults$beta[, paste0(vars, "_lag1")]
}
}
