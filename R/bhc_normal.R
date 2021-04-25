#' Bayesian hierarchical clustering for normal variates
#'
#' @param x A numeric matrix with rows representing observations
#' to cluster.
#' @param alpha Dirichlet process hyperparameter.
#' @param nu Normal-Inverse-Wishart (NIV) prior degrees of freedom.
#' @param kappa NIV prior number of observations.
#' @param mu NIV prior mean.
#' @param lambda NIV prior sum of squares matrix.
#' @details Data observations are assumed to follow a multivariate normal
#' distribution. A normal-inverse-wishart distribution is used as the prior
#' distribution for the mean and covariance.
bhc.normal = function(
    x
    , alpha=1
    , nu=ncol(x)
    , kappa=1
    , mu=colMeans(x, na.rm=TRUE)
    , lambda=var(x, na.rm=TRUE)*ncol(x))
{
    stopifnot(is.matrix(x) && is.numeric(x))
    stopifnot(is.numeric(alpha) && alpha[1L] > 0)
    stopifnot(is.numeric(nu) && nu[1L] >= ncol(x))
    stopifnot(is.numeric(kappa) && kappa[1L] > 0)
    stopifnot(is.numeric(mu) && length(mu) == ncol(x))
    stopifnot(is.matrix(lambda)
        && is.numeric(lambda)
        && isSymmetric(lambda)
        && ncol(lambda) == ncol(x))

    storage.mode(x) = "double"
    storage.mode(alpha) = "double"
    storage.mode(nu) = "double"
    storage.mode(kappa) = "double"
    storage.mode(mu) = "double"
    storage.mode(lambda) = "double"

    obj = .Call(bhc_normal
        , t(x)
        , alpha[1L]
        , nu[1L]
        , kappa[1L]
        , mu
        , lambda)

    obj$labels = rownames(x)

    return (obj)
}
