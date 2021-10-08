#include "bhc.h"

static struct cluster *merge_start(int k, struct cluster *a, struct cluster *b,
    struct hyperparam *h, merge_fn merge, ml_fn ml)
{
    double tmp;

    struct cluster *c = 0;
    struct cluster *lf = a->n <= b->n ? a : b;
    struct cluster *rt = lf == a ? b : a;

    // temporarily link left and right clusters
    lf->datum_last->next = rt->datum_first;

    // merge function malloc's c
    merge(k, lf, rt, h, &c);

    tmp = c->log_alpha + lgammafn(c->n);

    c->log_dk = logspace_add(tmp, a->log_dk + b->log_dk);

    //c->log_pk = tmp - c->log_dk;

    c->log_pk = -log1p(exp(a->log_dk + b->log_dk - tmp));

    if (c->log_pk == 0)
        error("loss of precision in probability calculation");

    c->log_h = ml(c, h);

    c->log_d = logspace_add(
        c->log_pk + c->log_h,
        log(-expm1(c->log_pk)) + a->log_d + b->log_d
    );

    c->log_r = c->log_pk + c->log_h - c->log_d;

    // undo link
    lf->datum_last->next = 0;

    return c;
}


static void merge_finish(struct cluster *c) {
    c->lf->parent = c;
    c->rt->parent = c;
    c->lf->datum_last->next = c->rt->datum_first;
}


static double cluster_height(struct cluster *c) {
    double h = 0;
    while (c) {
        h -= c->log_r;
        c = c->parent;
    }
    return h;
}


void plt_order(int index, int n, int *merge, int *order, int *xpos) {

    if (index > 0) {
        int lf = merge[(index-1) + 0 * n];
        int rt = merge[(index-1) + 1 * n];
        plt_order(lf, n, merge, order, xpos);
        plt_order(rt, n, merge, order, xpos);
    } else {
        order[(*xpos)++] = -index;
    }
    return;
}


/* On input there are n clusters (corresponding to n data points)
** maintained in a doubly-linked list with *hd as the head and *tl
** as the tail. This routine agglomerates clusters until only a
** single cluster remains. It uses the Bayesian Hierarchical
** Clustering algorithm described here:
**
** https://www2.stat.duke.edu/~kheller/bhcnew.pdf
**
** The merge, ml, and destroy pointers refer to functions that
** are specific to a data model. See bhc_normal.c and bhc_multinomial.c
** for examples related to the Normal and Multinomial distibutions,
** respectively.
*/
void bhc(
    unsigned int n,
    struct cluster **hd,
    struct cluster **tl,
    struct cluster **fl,
    struct hyperparam *h,
    merge_fn merge,
    ml_fn ml,
    destroy_fn destroy)
{
    unsigned int i;

    struct cluster *a;
    struct cluster *b;
    struct cluster *c;
    struct cluster *tmp;

    struct cluster *head = *hd;
    struct cluster *tail = *tl;

    for (i = 0; i < n-1; ++i) {
        c = NULL;
        for (a = tail; a != head; a = a->prev) {
            for (b = a->prev; b != 0; b = b->prev) {
                if (!c) {
                    c = merge_start(i+1, a, b, h, merge, ml);
                    continue;
                }
                tmp = merge_start(i+1, a, b, h, merge, ml);
                if (tmp->log_r > c->log_r) {
                    destroy(c);
                    c = tmp;
                    continue;
                }
                destroy(tmp);
            }
        }

        merge_finish(c);

        // remove a and b from the remaining clusters
        a = c->lf;
        b = c->rt;

        if (a == tail) {
            tail = a->prev;
            if (b == tail)
                tail = b->prev;
        }
        if (b == tail) {
            tail = b->prev;
            if (a == tail)
                tail = a->prev;
        }
        if (a == head) {
            head = a->next;
            if (b == head)
                head = b->next;
        }
        if (b == head) {
            head = b->next;
            if (a == head)
                head = a->next;
        }

        if (a->prev)
            a->prev->next = a->next;
        if (a->next)
            a->next->prev = a->prev;
        if (b->prev)
            b->prev->next = b->next;
        if (b->next)
            b->next->prev = b->prev;

        a->next = 0;
        a->prev = 0;
        b->next = 0;
        b->prev = 0;

        // record a and b for later free'ing
        a->prev = *fl;
        *fl = a;
        b->prev = *fl;
        *fl = b;

        // and add c to the remaining clusters
        c->prev = tail;
        if (tail)
            tail->next = c;
        tail = c;

    }
    // tail is now the root of the hierarchy
    *tl = tail;
    *hd = head;
}
