#ifndef BHC_INIT_H
#define BHC_INIT_H

#include <R.h>
#include <Rinternals.h>

SEXP bhc_normal(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP bhc_multinomial(SEXP, SEXP, SEXP);

#endif
