#ifndef BHC_H
#define BHC_H

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>

struct hyperparam;

struct datum {
    struct datum *next;
};

struct cluster {
    // cluster index
    int k;
    // number of obs in cluster
    unsigned int n;

    // log Dirichlet hyperparameter
    double log_alpha;

    // expected number of clusters
    // coeffient
    double log_dk;

    // prior probability of data being
    // generated from a single cluster
    // at this level
    double log_pk;

    // log posterior of data being
    // generated from single cluster
    double log_r;

    // log probability that data are
    // from single cluster
    double log_h;

    // log marginal probability of data
    double log_d;

    // linked list of data points belonging to cluster
    struct datum *datum_first;
    struct datum *datum_last;

    // left and right descendants in cluster hierarchy
    struct cluster *lf;
    struct cluster *rt;
    // ancestral cluster
    struct cluster *parent;

    // maintain all clusters in doubly-linked list where
    // left-to-right reading reveals merge order
    struct cluster *next;
    struct cluster *prev;
};

// function to merge two clusters
typedef void (*merge_fn)(int k, struct cluster *a, struct cluster *b, struct hyperparam *h, struct cluster **c);
// marginal log likelihood function that all points from one cluster
typedef double (*ml_fn)(struct cluster *c, struct hyperparam *h);
// cluster destructor function
typedef void (*destroy_fn)(struct cluster *c);

void bhc(unsigned int n, struct cluster **head, struct cluster **tail,
    struct cluster **fl, struct hyperparam *h, merge_fn merge, ml_fn ml,
    destroy_fn destroy);

void plt_order(int index, int n, int *merge, int *order, int *xpos);

#endif
