#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
  Check these declarations against the C/Fortran source code.
  */

/* .Call calls */

extern SEXP PAFit_normalized_constant(SEXP , SEXP , SEXP , SEXP , 
                                      SEXP , SEXP); 
extern SEXP PAFit_normalized_constant_alpha(SEXP , SEXP , SEXP , SEXP , 
                                     SEXP , SEXP , SEXP , SEXP );

extern SEXP PAFit_get_stats(SEXP , SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP , SEXP , 
                     SEXP , SEXP ); 

extern SEXP PAFit_update_f(SEXP , SEXP , SEXP , SEXP , SEXP , 
                    SEXP , SEXP , SEXP , SEXP , 
                    SEXP ); 

extern SEXP PAFit_update_offset(SEXP , SEXP , 
                         SEXP , SEXP , SEXP , 
                         SEXP , SEXP ); 

extern SEXP PAFit_update_f_alpha(SEXP , SEXP , SEXP , 
                          SEXP , SEXP , SEXP , 
                          SEXP , SEXP , SEXP , 
                          SEXP , SEXP );

extern SEXP PAFit_update_f_alpha_new(SEXP , SEXP , SEXP , 
                              SEXP , SEXP , SEXP , 
                              SEXP , SEXP , SEXP , 
                              SEXP, SEXP , SEXP ); 

extern SEXP PAFit_update_alpha_fast(SEXP , SEXP , SEXP , 
                             SEXP , SEXP , SEXP , SEXP , 
                             SEXP , SEXP , SEXP , SEXP ); 

extern SEXP PAFit_var_alpha(SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP , SEXP ); 

extern SEXP PAFit_coeff_theta(SEXP , SEXP , SEXP , 
                       SEXP , SEXP ); 

extern SEXP PAFit_coeff_var(SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP ); 


extern SEXP PAFit_cal_var_f(SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP , SEXP , 
                     SEXP , SEXP );

extern SEXP PAFit_cal_var_f_new(SEXP , SEXP , SEXP , 
                         SEXP , SEXP , SEXP , SEXP , 
                         SEXP, SEXP , SEXP );
extern SEXP PAFit_update_f_new(SEXP , SEXP , 
                               SEXP , SEXP , 
                               SEXP , SEXP , 
                               SEXP , SEXP , 
                               SEXP , SEXP , 
                               SEXP );

  
static const R_CallMethodDef CallEntries[] = {
  {"PAFit_cal_var_f_new", (DL_FUNC) &PAFit_cal_var_f_new,10},
  {"PAFit_cal_var_f", (DL_FUNC) &PAFit_cal_var_f, 9},
  {"PAFit_coeff_var", (DL_FUNC) &PAFit_coeff_var, 6},
  {"PAFit_coeff_theta", (DL_FUNC) &PAFit_coeff_theta, 5},
  {"PAFit_var_alpha", (DL_FUNC) &PAFit_var_alpha, 11},
  {"PAFit_update_alpha_fast", (DL_FUNC) &PAFit_update_alpha_fast, 11},
  {"PAFit_update_f_alpha_new", (DL_FUNC) &PAFit_update_f_alpha_new, 12},
  {"PAFit_update_f_alpha", (DL_FUNC) &PAFit_update_f_alpha, 11},
  {"PAFit_update_offset", (DL_FUNC) &PAFit_update_offset, 7},
  {"PAFit_update_f", (DL_FUNC) &PAFit_update_f, 10},
  {"PAFit_get_stats", (DL_FUNC) &PAFit_get_stats, 23},
  {"PAFit_normalized_constant_alpha", (DL_FUNC) &PAFit_normalized_constant_alpha, 8},
  {"PAFit_normalized_constant", (DL_FUNC) &PAFit_normalized_constant, 6},
  {"PAFit_update_f_new", (DL_FUNC) &PAFit_update_f_new, 11},
  {NULL, NULL, 0}
};

void R_init_PAFit(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}