#include<R_ext/Rdynload.h>
#ifndef R_R_H
#  include <R.h>
#endif

void F77_NAME(perturb_times)(int*, int*, int*, double*, double*,  
      double*, double*, double*, double*, double*);
void F77_NAME(perturb_times2)(int*, int*, double*, double*, double*, 
              double*, double*, double*, double*, double*, double*);

void F77_NAME(perturb_steady)(int*, double*, double*,  
              double*, double*, double*, double*, double*,
              int*, double*, double*, double*, double*);

void F77_NAME(perturb_event)(int*, int*, double*, double*, 
              double*, double*, double*, int*, 
              double*, double*, double*);
void F77_NAME(perturb_event2)(int*,  double*, double*, 
              double*, double*, double*, int*, double*);

void F77_NAME(logistic_time)(int*, int*, double*, double*, double*, 
              double*, double*, double*, double*, double*);

void F77_NAME(metier_time)(int*, int*, int*, double*, double*, double*,  
              double*, double*, double*, double*, double*, double*);
void F77_NAME(metier_event)(int*, int*, int*, double*, double*,   
              double*, double*, double*, int*, 
              double*, double*, double*, double*);

R_FortranMethodDef fortranMethods[] = {
 {"perturb_times",  (DL_FUNC) &F77_SUB(perturb_times),   8},
 {"perturb_times2", (DL_FUNC) &F77_SUB(perturb_times2),  9},
 {"perturb_steady", (DL_FUNC) &F77_SUB(perturb_steady), 13},
 {"perturb_event",  (DL_FUNC) &F77_SUB(perturb_event ), 11},
 {"perturb_event2", (DL_FUNC) &F77_SUB(perturb_event2),  8},
 {"logistic_time",  (DL_FUNC) &F77_SUB(logistic_time),  10},
 {"metier_time",    (DL_FUNC) &F77_SUB(metier_time),    12},
 {"metier_event",   (DL_FUNC) &F77_SUB(metier_event),   13},
 {NULL, NULL, 0}
};
void R_init_bfiat(DllInfo *info) {
  R_registerRoutines(info, NULL, NULL, fortranMethods, NULL);
  R_useDynamicSymbols(info, FALSE); // disable dynamic searching  
}
