#' Bayesian hierarchical clustering for multinomial variates
#'
#' @param x An integer matrix with rows representing observations
#' to cluster.
#' @param alpha Dirichlet process hyperparameter.
#' @param beta Dirichlet distribution hyperparameter.
#' @details Data observations are assumed to follow a multinomial
#' distribution. A Dirichlet distribution is used as the prior
#' distribution for the multinomial probabilities.
#' @examples
#' data(snakediet)
#' x = xtabs(snakediet$count ~ snakediet$species + snakediet$food)
#' res = bhc.multinomial(x)
bhc.multinomial = function(
    x
    , alpha=1
    , beta=rep(1, ncol(x)))
{
    stopifnot(is.matrix(x) && is.numeric(x))
    stopifnot(is.numeric(alpha) && alpha[1L] > 0)
    stopifnot(is.numeric(beta) && length(beta) == ncol(x))

    storage.mode(x) = "integer"
    storage.mode(alpha) = "double"
    storage.mode(beta) = "double"

    obj = .Call(bhc_multinomial
        , t(x)
        , alpha[1L]
        , beta)

    obj$labels = rownames(x)

    return (obj)
}
