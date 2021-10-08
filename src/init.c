#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

#include "init.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}


static const R_CallMethodDef CallEntries[] = {
    CALLDEF(bhc_normal, 6),
    CALLDEF(bhc_multinomial, 3),
    CALLDEF(bhc_poisson, 5),
    {NULL, NULL, 0}
};


void attribute_visible R_init_bhc(DllInfo *info)
{
    R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
}
