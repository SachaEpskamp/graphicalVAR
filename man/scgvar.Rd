\name{scgvar}
\alias{scgvar}
\title{
Subgroup chain graphical VAR estimation
}
\description{
Estimates pooled and individual graphical VAR networks, clusters subjects by
signed network similarity, and refits a pooled graphical VAR in each subgroup.
Estimation uses \code{\link{graphicalVAR}} and its
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
  and at least three usable lag-1 pairs within days.}
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
  series as one day. Lag pairs never cross day boundaries.}
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
  \item{verbose}{Logical. Print progress messages.}
  \item{redo}{Reserved argument; must be \code{FALSE}.}
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
  \item{mlgvclust}{Reserved argument; must be \code{TRUE}.}
  \item{indplot}{Logical. Also export an individual-network PDF for each
  selected subject. Requires a non-\code{NULL} \code{dir}.}
  \item{elim.groups}{Logical. Before computing similarities, set individual
  PCC and PDC entries to zero wherever the corresponding pooled fixed network
  has a nonzero entry, regardless of sign.}
  \item{g.ID}{Reserved argument; must be \code{NULL}.}
  \item{nclust}{\code{NULL} for the walktrap-selected number of communities,
  or an integer from one to the number of selected subjects specifying a
  dendrogram cut. A cut cannot have fewer clusters than the graph has connected
  components. Infeasible candidate graphs are excluded from the search.}
  \item{returnResults}{Logical. Invisibly return the structured result when
  \code{TRUE}. If \code{FALSE}, \code{dir} is required and the function
  invisibly returns \code{NULL} after export.}
  \item{\dots}{Additional named arguments to
  \code{\link{graphicalVAR}}, such as \code{nLambda}, \code{lambda_beta},
  \code{lambda_kappa}, \code{regularize_mat_beta},
  \code{regularize_mat_kappa}, \code{maxit.out}, and \code{ebic_tol}.
  These apply to pooled, individual, and subgroup fits. \code{nLambda} must
  be one positive integer. Explicit lambda vectors must contain distinct,
  finite, nonnegative values; they override generated grids and their minimum
  ratios as in \code{graphicalVAR}.}
}
\details{
The VAR(1) model is
\deqn{y_t = a + B y_{t-1} + \epsilon_t, \qquad
       \epsilon_t \sim N(0, K^{-1}).}
Pooled and individual estimation follows the data preparation used by
\code{\link{mlGraphicalVAR}}. Individual data are prepared again for each
subject-specific fit. Subgroup models are refitted from their original rows.

Similarity counts same-sign PCC entries that are nonzero in both subjects,
plus same-sign PDC entries with absolute values greater than 0.01 in both.
The full PCC matrix counts each contemporaneous edge twice. PDC entries count
once and include autoregressive effects. Shared zeros do not contribute, and
the similarity diagonal is zero.

With \code{search = TRUE}, distinct graphs are formed by zeroing similarities
at or below thresholds between the 0.01 and 0.99 type-7 quantiles of the upper
triangle. Weighted walktrap uses four steps. For each community, its external
edge weight is divided by its total degree; a zero-degree community receives
a value of one. The selected graph maximizes one minus the mean of these
values. Scores within \eqn{10^{-12}} are tied, favoring the lower threshold.
When supplied, \code{nclust} determines the partition used for scoring.

Empty graphs and infeasible cuts are excluded. If none of the thresholded
graphs is feasible, the original graph is used with a warning when feasible.
An all-zero similarity matrix or an unattainable cut produces an error.
With \code{search = FALSE}, the minimum positive similarity is subtracted
from positive entries and \code{.Machine$double.eps} is added; zeros stay zero.

Returned model fits are \code{graphicalVAR} objects. Beta matrices have outcome
rows, an intercept column named \code{"1"}, and lagged predictor columns named
\code{paste0(vars, "_lag1")}. PDC matrices have origin rows and destination
columns. The existing \code{plot} and \code{summary} methods apply to these fits.
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
\author{
Jonathan J. Park and Sacha Epskamp
}
\references{
Park, J. J., Chow, S. M., Epskamp, S., & Molenaar, P. C. (2024).
Subgrouping with chain graphical VAR models.
\emph{Multivariate Behavioral Research}, \bold{59}(3), 543--565.
}
\seealso{
\code{\link{graphicalVAR}}, \code{\link{mlGraphicalVAR}},
\code{\link{graphicalVARsim}}, \code{\link[igraph]{cluster_walktrap}}
}
