#' Bayesian hierarchical clustering for poisson counts
#'
#' @param x An integer matrix with rows representing observations
#' to cluster.
#' @param sz The number of sampled units for each row in \code{x}.
#' @param shape Shape of gamma prior.
#' @param rate Rate of gamma prior.
#' @param alpha Dirichlet process hyperparameter.
#' @details Data observations are assumed to follow a Poisson
#' distribution. A Gamma distribution is used as the prior
#' distribution for the Poisson rates.
#' @examples
#' data(snakediet)
#' x = xtabs(snakediet$count ~ snakediet$species + snakediet$food)
#' res = bhc.poisson(x, rowSums(x))
bhc.poisson = function(
    x
    , sz
    , shape=1
    , rate=1
    , alpha=1)
{
    stopifnot(is.matrix(x) && is.numeric(x))
    stopifnot(is.numeric(shape) && shape[1L] > 0)
    stopifnot(is.numeric(rate) && rate[1L] > 0)
    stopifnot(is.numeric(alpha) && alpha[1L] > 0)

    storage.mode(x) = "integer"
    storage.mode(sz) = "integer"
    storage.mode(shape) = "double"
    storage.mode(rate) = "double"
    storage.mode(alpha) = "double"

    obj = .Call(bhc_poisson
        , t(x)
        , sz
        , shape
        , rate
        , alpha)

    obj$labels = rownames(x)

    return (obj)
}
