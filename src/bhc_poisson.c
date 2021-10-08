#include "bhc.h"


// hyperparameters for the Gamma prior on the Poisson distribution
// and the Dirichlet process prior
struct hyperparam {
    // number of categories
    unsigned int p;

    // shape of Gamma
    double shape;

    // rate of Gamma
    double rate;

    // Dirichlet process hyperparameter
    double alpha;
};


static void hyperparam_free(struct hyperparam *h) {
    if (h)
        free(h);
}


static struct hyperparam *hyperparam_alloc(
    unsigned int p, double shape, double rate, double alpha)
{
    struct hyperparam *h = malloc(sizeof(*h));
    memset(h, 0, sizeof(*h));
    h->p = p;
    h->shape = shape;
    h->rate = rate;
    h->alpha = alpha;
    return h;
}

struct datum_poisson {
    struct datum base;
    unsigned int p;
    int size;
    int *x;
};


static void datum_alloc(unsigned int p, int sz, int *x, struct datum **d) {
    struct datum_poisson *datum = malloc(sizeof(*datum));
    memset(datum, 0, sizeof(*datum));
    datum->p = p;
    datum->size = sz;
    datum->x = malloc(p * sizeof(int));
    memcpy(datum->x, x, p * sizeof(int));
    *d = (struct datum *)datum;
}


static void datum_free(struct datum *d) {
    struct datum_poisson *datum = (struct datum_poisson *)d;
    free(datum->x);
    free(datum);
}

struct cluster_poisson {
    struct cluster base;

    // number of observations
    int size;

    // aggregated counts in each category
    int *count;
};

static void cluster_free(struct cluster *cl) {
    struct cluster_poisson *c = (struct cluster_poisson *)cl;
    free(c->count);
    free(c);
    c = NULL;
}

// compute log marginal probability that all data in c are
// from single Poisson-gamma distribution
static double compute_log_h(struct cluster *cl, struct hyperparam *h)
{
    int j;
    double loglk = 0;

    struct cluster_poisson *c = (struct cluster_poisson *)cl;

    for (j = 0; j < h->p; ++j)
    {
        loglk += 
            h->shape * log(h->rate) - 
            (h->shape + c->count[j]) * log(h->rate + c->size) +
            lgammafn(h->shape + c->count[j]) -
            lgammafn(h->shape);
    }

    return loglk;
}

static void cluster_alloc(
    int k, int sz, int *x, struct hyperparam *h, struct cluster **cl)
{
    struct cluster_poisson *c = malloc(sizeof(*c));
    memset(c, 0, sizeof(*c));
    c->base.k = k;
    c->base.log_alpha = log(h->alpha);
    c->count = calloc(h->p, sizeof(int));

    if (x) {
        c->base.n = 1;
        c->base.log_pk = 0;
        c->base.log_dk = c->base.log_alpha;
        datum_alloc(h->p, sz, x, &(c->base.datum_first));
        c->base.datum_last = c->base.datum_first;
        c->size = ((struct datum_poisson *)(c->base.datum_first))->size;
        memcpy(c->count, x, h->p * sizeof(int));
        c->base.log_h = compute_log_h(&c->base, h);
        c->base.log_d = c->base.log_h;
    }
    *cl = (struct cluster *)c;
}

static void cluster_merge(int k, struct cluster *lf, struct cluster *rt,
    struct hyperparam *h, struct cluster **cl)
{
    cluster_alloc(k, 0, 0, h, cl);

    unsigned int i;
    unsigned int p;
    struct cluster_poisson *a = (struct cluster_poisson *)lf;
    struct cluster_poisson *b = (struct cluster_poisson *)rt;
    struct cluster_poisson *c = (struct cluster_poisson *)(*cl);

    (*cl)->datum_first = lf->datum_first;
    (*cl)->datum_last = rt->datum_last;

    (*cl)->lf = lf;
    (*cl)->rt = rt;
    (*cl)->n = lf->n + rt->n;

    p = h->p;

    c->size = a->size + b->size;

    for (i = 0; i < p; ++i)
        c->count[i] = a->count[i] + b->count[i];
}


SEXP bhc_poisson(SEXP x, SEXP sz, SEXP shape, SEXP rate, SEXP alpha)
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
        , REAL(shape)[0]
        , REAL(rate)[0]
        , REAL(alpha)[0]);

    // initialize the terminal clusters
    cluster_alloc(-1, INTEGER(sz)[0], INTEGER(x), h, &head);
    tail = head;
    for (i = 1; i < n; ++i) {
        cluster_alloc(-(i+1), INTEGER(sz)[i], INTEGER(x) + i*p, h, &c);
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
