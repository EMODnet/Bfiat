#include<R_ext/Rdynload.h>
#ifndef R_R_H
#  include <R.h>
#endif

void F77_NAME(logistictrawl)(int*, int*, int*, double*, double*,  
      double*, double*, double*, double*, double*);
     
R_FortranMethodDef fortranMethods[] = {
 {"logistictrawl", (DL_FUNC) &F77_SUB(logistictrawl), 8},
 {NULL, NULL, 0}
};
void R_init_bfiat(DllInfo *info) {
  R_registerRoutines(info, NULL, NULL, fortranMethods, NULL);
  R_useDynamicSymbols(info, FALSE); // disable dynamic searching  
}
