#include "bhc.h"

// hyperparameters for the Normal-inverse Wishart prior
// and the Dirichlet process prior
struct hyperparam {
    // dimension of the NIV distribution
    unsigned int p;
    // prior degrees of freedom
    double nu0;
    // prior number of observations
    double kappa0;
    // prior mean vector
    double *mu0;
    // prior sum of squares matrix
    double *lambda0;
    // log normalization constant
    double logz0;
    // Dirichlet process hyperparameter
    double alpha;
};


struct datum_normal {
    struct datum base;
    unsigned int p;
    double *x;
};


struct cluster_normal {
    struct cluster base;

    // mean of all data in cluster
    double *m;

    // sum of squares matrix of all data in cluster
    double *S;
};


static void hyperparam_free(struct hyperparam *h) {
    if (h) {
        free(h->mu0);
        free(h->lambda0);
        free(h);
    }
}


// modified from function det_ge_real, lines 1260-1316 in
// src/modules/lapack/Lapack.c in R source tree
static double logdet(int n, double *Ain) {
    if (n == 1)
        return log(Ain[0]);
    int i;
    int info;
    int jpvt[n];
    double dii;
    double modulus = 0.0;
    double A[n*n];
    memset(jpvt, 0, n * sizeof(int));
    memcpy(A, Ain, n*n*sizeof(double));
    F77_CALL(dgetrf)(&n, &n, A, &n, jpvt, &info);
    if (info != 0)
        error("Lapack dgetrf() did not return success");
    for (i = 0; i < n; i++) {
        dii = A[i + i * n]; /* i-th diagonal element */
        modulus += log(dii < 0 ? -dii : dii);
    }
    return modulus;
}


static double multigammaln(double x, int d) {
    int i;
    double val = 0.25 * d * (d-1) * log(M_PI);
    for (i = 1; i <= d; ++i) {
        val += lgammafn( x + (1-i)/2 );
    }
    return val;
}

// compute log normalization constant for Normal-inverse Wishart
static double compute_log_z(
    double p, double nu, double kappa, double *mu, double *lambda)
{
    return M_LN2 * (nu*p/2)
        + (p/2) * log(M_2PI/kappa)
        + multigammaln(nu/2, p)
        - (nu/2) * logdet(p, lambda);
}


static struct hyperparam *hyperparam_alloc(
    unsigned int p, double alpha, double nu,
    double kappa, double *mu, double *lambda)
{
    struct hyperparam *h = malloc(sizeof(*h));
    h->p = p;
    h->nu0 = nu;
    h->kappa0 = kappa;
    h->mu0 = malloc(p * sizeof(double));
    h->lambda0 = malloc(p * p * sizeof(double));
    memcpy(h->mu0, mu, p * sizeof(double));
    memcpy(h->lambda0, lambda, p * p * sizeof(double));
    h->alpha = alpha;
    h->logz0 = compute_log_z((double)p, nu, kappa, mu, lambda);
    return h;
}


// compute log marginal probability that all data in c are
// from single Normal-inverse Wishart distribution
static double compute_log_h(struct cluster *cl, struct hyperparam *h)
{
    unsigned int i;
    unsigned int j;
    unsigned int p = h->p;
    double nu;
    double kappa;
    double logz;
    double mu[p];
    double dt[p];
    double lambda[p*p];
    double back[p*p];

    struct cluster_normal *c = (struct cluster_normal *)cl;

    nu = h->nu0 + cl->n;
    kappa = h->kappa0 + cl->n;

    for (j = 0; j < p; ++j) {
        mu[j] = (h->kappa0 * h->mu0[j] + cl->n * c->m[j]) / kappa;
        dt[j] = c->m[j] - h->mu0[j];
    }

    for (i = 0; i < p; ++i) {
        for (j = 0; j <= i; ++j) {
            back[i + j * p] = dt[i] * dt[j];
            back[j + i * p] = back[i + j * p];
        }
    }

    for (j = 0; j < p*p; ++j) {
        lambda[j] = h->lambda0[j] + c->S[j] +
            (h->kappa0 * cl->n/kappa) * back[j];
    }

    logz = compute_log_z((double)p, nu, kappa, mu, lambda);

    return logz - h->logz0 - M_LN_2PI * (cl->n * p / 2);
}


static void datum_alloc(unsigned int p, double *x, struct datum **d) {
    struct datum_normal *datum = malloc(sizeof(*datum));
    memset(datum, 0, sizeof(*datum));
    datum->p = p;
    datum->x = malloc(p * sizeof(double));
    memcpy(datum->x, x, p * sizeof(double));
    *d = (struct datum *)datum;
}


static void datum_free(struct datum *d) {
    struct datum_normal *datum = (struct datum_normal *)d;
    free(datum->x);
    free(datum);
}


static void cluster_free(struct cluster *cl) {
    struct cluster_normal *c = (struct cluster_normal *)cl;
    free(c->m);
    free(c->S);
    free(c);
    c = NULL;
}


static void cluster_alloc(
    int k, double *x, struct hyperparam *h, struct cluster **cl)
{
    struct cluster_normal *c = malloc(sizeof(*c));
    memset(c, 0, sizeof(*c));
    c->base.k = k;
    c->base.log_alpha = log(h->alpha);
    c->m = calloc(h->p, sizeof(double));
    c->S = calloc(h->p*h->p, sizeof(double));

    if (x) {
        c->base.n = 1;
        c->base.log_pk = 0;
        c->base.log_dk = c->base.log_alpha;
        datum_alloc(h->p, x, &(c->base.datum_first));
        c->base.datum_last = c->base.datum_first;
        memcpy(c->m, x, h->p * sizeof(double));
        c->base.log_h = compute_log_h(&c->base, h);
        c->base.log_d = c->base.log_h;
    }
    *cl = (struct cluster *)c;
}


static void cluster_merge(int k, struct cluster *lf, struct cluster *rt,
    struct hyperparam *h, struct cluster **cl)
{
    cluster_alloc(k, 0, h, cl);

    struct datum *datum;
    struct datum_normal *d;
    struct cluster_normal *c = (struct cluster_normal *)(*cl);

    (*cl)->datum_first = lf->datum_first;
    (*cl)->datum_last = rt->datum_last;

    (*cl)->lf = lf;
    (*cl)->rt = rt;
    (*cl)->n = lf->n + rt->n;

    unsigned int i;
    unsigned int j;
    unsigned int n = 0;

    unsigned int p = h->p;
    double x;
    double y;
    double m[p];

    /* compute cluster mean and sum of squares matrix */

    for (datum = (*cl)->datum_first; datum != 0; datum = datum->next)
    {
        d = (struct datum_normal *)datum;
        ++n;
        if (n == 1) {
            for (i = 0; i < p; ++i)
                c->m[i] = ISNAN(d->x[i]) ? 0 : d->x[i];
        } else {
            for (i = 0; i < p; ++i) {
                x = ISNAN(d->x[i]) ? c->m[i] : d->x[i];
                m[i] = c->m[i];
                c->m[i] += (x - m[i]) / n;
            }
            for (i = 0; i < p; ++i) {
                x = ISNAN(d->x[i]) ? c->m[i] : d->x[i];
                for (j = 0; j <= i; ++j) {
                    y = ISNAN(d->x[j]) ? c->m[j] : d->x[j];
                    c->S[i+j*p] += (x - m[i]) * (y - c->m[j]);
                    c->S[j+i*p] = c->S[i+j*p];
                }
            }
        }
    }

    /*

    // here's a way to compute cluster means and sum of squares using
    // only the mean and sos from the two merged clusters so that we
    // avoid iterating over all the data again and again. i don't know
    // how accurate this is, though, so unsure if we should use it.

    for (i = 0; i < p; ++i)
        c->m[i] = (lf->n * lf->m[i] + rt->n * rt->m[i]) / (lf->n + rt->n);

    for (i = 0; i < p; ++i) {
        for (j = 0; j <= i; ++j) {
            c->S[i + j * p] += lf->S[i + j * p];
            c->S[i + j * p] += rt->S[i + j * p];
            c->S[i + j * p] += lf->n * (lf->m[i] - c->m[i]) * (lf->m[j] - c->m[j]);
            c->S[i + j * p] += rt->n * (rt->m[i] - c->m[i]) * (rt->m[j] - c->m[j]);
            c->S[j + i * p] = c->S[i + j * p];
        }
    }

    */

}


SEXP bhc_normal(SEXP x, SEXP alpha, SEXP nu, SEXP kappa, SEXP mu, SEXP lambda)
{
    int i;
    unsigned int p = INTEGER(getAttrib(x, R_DimSymbol))[0];
    unsigned int n = INTEGER(getAttrib(x, R_DimSymbol))[1];

    struct cluster *a = 0;
    struct cluster *c = 0;
    struct cluster *head = 0;
    struct cluster *tail = 0;
    struct cluster *fl = 0;

    struct hyperparam *h = hyperparam_alloc(
        p
        , REAL(alpha)[0]
        , REAL(nu)[0]
        , REAL(kappa)[0]
        , REAL(mu)
        , REAL(lambda));

    // initialize the terminal clusters
    cluster_alloc(-1, REAL(x), h, &head);
    tail = head;
    for (i = 1; i < n; ++i) {
        cluster_alloc(-(i+1), REAL(x) + i*p, h, &c);
        c->prev = tail;
        tail->next = c;
        tail = c;
    }

    bhc(n, &head, &tail, &fl, h, &cluster_merge, &compute_log_h, &cluster_free);

    // tail is now the root of the hierarchy

    SEXP result = PROTECT(allocVector(VECSXP, 4));

    SET_VECTOR_ELT(result, 0, allocMatrix(INTSXP, n-1, 2));
    SET_VECTOR_ELT(result, 1, allocVector(REALSXP, n-1));
    SET_VECTOR_ELT(result, 2, allocVector(INTSXP, n));

    a = fl;

    while (a) {
        if (a->lf && a->rt) {
            INTEGER(VECTOR_ELT(result, 0))[(a->k-1) + 0 * (n-1)] = a->lf->k;
            INTEGER(VECTOR_ELT(result, 0))[(a->k-1) + 1 * (n-1)] = a->rt->k;
            REAL(VECTOR_ELT(result, 1))[a->k-1] = a->log_r;
        }
        a = a->prev;
    }

    INTEGER(VECTOR_ELT(result, 0))[(tail->k-1) + 0 * (n-1)] = tail->lf->k;
    INTEGER(VECTOR_ELT(result, 0))[(tail->k-1) + 1 * (n-1)] = tail->rt->k;
    REAL(VECTOR_ELT(result, 1))[tail->k-1] = tail->log_r;

    int xpos = 0;
    plt_order(tail->k, n-1,
        INTEGER(VECTOR_ELT(result, 0)), INTEGER(VECTOR_ELT(result, 2)), &xpos);

    SET_VECTOR_ELT(result, 3, ScalarReal(tail->log_d));

    SEXP names = PROTECT(allocVector(STRSXP, 4));
    SET_STRING_ELT(names, 0, mkChar("merge"));
    SET_STRING_ELT(names, 1, mkChar("height"));
    SET_STRING_ELT(names, 2, mkChar("order"));
    SET_STRING_ELT(names, 3, mkChar("lnL"));

    setAttrib(result, R_NamesSymbol, names);
    setAttrib(result, R_ClassSymbol, mkString("hclust"));

    cleanup:
        hyperparam_free(h);
        while (fl) {
            c = fl;
            fl = fl->prev;
            if (!c->lf && !c->rt) {
                datum_free(c->datum_first);
            }
            cluster_free(c);
        }
        cluster_free(tail);

    UNPROTECT(2);
    return result;
}
