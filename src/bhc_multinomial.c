#include "bhc.h"


// hyperparameters for the Dirichlet prior on the multinomial distribution
// and the Dirichlet process prior
struct hyperparam {
    // dimension of the Dirichlet distribution
    unsigned int p;

    // sum of over *beta
    double betasum;

    // log gamma of betasum
    double lgamma_betasum;

    // Dirichlet hyperparameter
    double *beta;

    // log gamma of each *beta
    double *lgamma_beta;

    // Dirichlet process hyperparameter
    double alpha;
};


static void hyperparam_free(struct hyperparam *h) {
    if (h) {
        free(h->beta);
        free(h->lgamma_beta);
        free(h);
    }
}


static struct hyperparam *hyperparam_alloc(
    unsigned int p, double alpha, double *beta)
{
    struct hyperparam *h = malloc(sizeof(*h));
    memset(h, 0, sizeof(*h));
    h->p = p;
    h->beta = malloc(p * sizeof(double));
    h->lgamma_beta = malloc(p * sizeof(double));
    memcpy(h->beta, beta, p * sizeof(double));
    for (unsigned int i = 0; i < p; ++i) {
        h->betasum += h->beta[i];
        h->lgamma_beta[i] = lgammafn(h->beta[i]);
    }
    h->lgamma_betasum = lgammafn(h->betasum);
    h->alpha = alpha;
    return h;
}

struct datum_multinom {
    struct datum base;
    unsigned int p;
    int size;
    int *x;
};


static void datum_alloc(unsigned int p, int *x, struct datum **d) {
    struct datum_multinom *datum = malloc(sizeof(*datum));
    memset(datum, 0, sizeof(*datum));
    datum->p = p;
    datum->x = malloc(p * sizeof(int));
    memcpy(datum->x, x, p * sizeof(int));
    for (unsigned int i = 0; i < p; ++i)
        datum->size += x[i];
    *d = (struct datum *)datum;
}


static void datum_free(struct datum *d) {
    struct datum_multinom *datum = (struct datum_multinom *)d;
    free(datum->x);
    free(datum);
}

struct cluster_multinom {
    struct cluster base;

    // sum of aggregated counts
    int size;

    // aggregated counts in each category
    int *count;
};

static void cluster_free(struct cluster *cl) {
    struct cluster_multinom *c = (struct cluster_multinom *)cl;
    free(c->count);
    free(c);
    c = NULL;
}

// compute log marginal probability that all data in c are
// from single Dirichlet-multinomial distribution
static double compute_log_h(struct cluster *cl, struct hyperparam *h)
{
    int j;
    double loglk = 0;

    struct cluster_multinom *c = (struct cluster_multinom *)cl;

    for (j = 0; j < h->p; ++j)
    {
        loglk += lgammafn(c->count[j] + h->beta[j]) - h->lgamma_beta[j];
    }

    loglk += h->lgamma_betasum - lgammafn(c->size + h->betasum);

    return loglk;
}

static void cluster_alloc(
    int k, int *x, struct hyperparam *h, struct cluster **cl)
{
    struct cluster_multinom *c = malloc(sizeof(*c));
    memset(c, 0, sizeof(*c));
    c->base.k = k;
    c->base.log_alpha = log(h->alpha);
    c->count = calloc(h->p, sizeof(int));

    if (x) {
        c->base.n = 1;
        c->base.log_pk = 0;
        c->base.log_dk = c->base.log_alpha;
        datum_alloc(h->p, x, &(c->base.datum_first));
        c->base.datum_last = c->base.datum_first;
        c->size = ((struct datum_multinom *)(c->base.datum_first))->size;
        memcpy(c->count, x, h->p * sizeof(int));
        c->base.log_h = compute_log_h(&c->base, h);
        c->base.log_d = c->base.log_h;
    }
    *cl = (struct cluster *)c;
}

static void cluster_merge(int k, struct cluster *lf, struct cluster *rt,
    struct hyperparam *h, struct cluster **cl)
{
    cluster_alloc(k, 0, h, cl);

    unsigned int i;
    unsigned int p;
    struct cluster_multinom *a = (struct cluster_multinom *)lf;
    struct cluster_multinom *b = (struct cluster_multinom *)rt;
    struct cluster_multinom *c = (struct cluster_multinom *)(*cl);

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


SEXP bhc_multinomial(SEXP x, SEXP alpha, SEXP beta)
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
        , REAL(beta));

    // initialize the terminal clusters
    cluster_alloc(-1, INTEGER(x), h, &head);
    tail = head;
    for (i = 1; i < n; ++i) {
        cluster_alloc(-(i+1), INTEGER(x) + i*p, h, &c);
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
