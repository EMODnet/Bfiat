#include<R_ext/Rdynload.h>
#ifndef R_R_H
#  include <R.h>
#endif

void F77_NAME(logistictrawl)(int*, int*, int*, double*, double*,  
      double*, double*, double*, double*, double*);
void F77_NAME(logistictrawl2)(int*, int*, double*, double*,  
              double*, double*, double*, double*, double*, double*);
void F77_NAME(steadydensity)(int*, double*, double*,  
              double*, double*, double*, double*, double*,
              int*, double*, double*, double*, double*);
void F77_NAME(eventdensity)(int*, int*, double*, double*, 
              double*, double*, double*, int*, double*);

void F77_NAME(eventdensity2)(int*,  double*, double*, 
              double*, double*, double*, int*, double*);
void F77_NAME(logistic)(int*, int*, double*, double*,  
              double*, double*, double*, double*, double*);

R_FortranMethodDef fortranMethods[] = {
 {"logistictrawl",  (DL_FUNC) &F77_SUB(logistictrawl),   8},
 {"logistictrawl2", (DL_FUNC) &F77_SUB(logistictrawl2),  8},
 {"steadydensity",  (DL_FUNC) &F77_SUB(steadydensity),  13},
 {"eventdensity",   (DL_FUNC) &F77_SUB(eventdensity ),   9},
 {"eventdensity2",  (DL_FUNC) &F77_SUB(eventdensity2),   8},
 {"logistic",       (DL_FUNC) &F77_SUB(logistic),        9},
 {NULL, NULL, 0}
};
void R_init_bfiat(DllInfo *info) {
  R_registerRoutines(info, NULL, NULL, fortranMethods, NULL);
  R_useDynamicSymbols(info, FALSE); // disable dynamic searching  
}
