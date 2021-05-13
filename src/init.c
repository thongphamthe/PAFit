#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>


extern SEXP _PAFit_normalized_constant(SEXP , SEXP , SEXP , SEXP , 
                                      SEXP , SEXP); 
extern SEXP _PAFit_normalized_constant_alpha(SEXP , SEXP , SEXP , SEXP , 
                                     SEXP , SEXP , SEXP , SEXP );

extern SEXP _PAFit_get_stats(SEXP , SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP , SEXP , 
                     SEXP , SEXP ); 


extern SEXP _PAFit_update_f(SEXP , SEXP , SEXP , SEXP , SEXP , 
                    SEXP , SEXP , SEXP , SEXP , 
                    SEXP ); 

extern SEXP _PAFit_update_offset(SEXP , SEXP , 
                         SEXP , SEXP , SEXP , 
                         SEXP , SEXP ); 

extern SEXP _PAFit_update_f_alpha(SEXP , SEXP , SEXP , 
                          SEXP , SEXP , SEXP , 
                          SEXP , SEXP , SEXP , 
                          SEXP , SEXP );

extern SEXP _PAFit_update_f_alpha_new(SEXP , SEXP , SEXP , 
                              SEXP , SEXP , SEXP , 
                              SEXP , SEXP , SEXP , 
                              SEXP, SEXP , SEXP ); 

extern SEXP _PAFit_update_alpha_fast(SEXP , SEXP , SEXP , 
                             SEXP , SEXP , SEXP , SEXP , 
                             SEXP , SEXP , SEXP , SEXP ); 

extern SEXP _PAFit_var_alpha(SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP , SEXP ); 

extern SEXP _PAFit_coeff_theta(SEXP , SEXP , SEXP , 
                       SEXP , SEXP ); 

extern SEXP _PAFit_coeff_var(SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP ); 


extern SEXP _PAFit_cal_var_f(SEXP , SEXP , SEXP , 
                     SEXP , SEXP , SEXP , SEXP , 
                     SEXP , SEXP );

extern SEXP _PAFit_cal_var_f_new(SEXP , SEXP , SEXP , 
                         SEXP , SEXP , SEXP , SEXP , 
                         SEXP, SEXP , SEXP );
extern SEXP _PAFit_update_f_new(SEXP , SEXP , 
                               SEXP , SEXP , 
                               SEXP , SEXP , 
                               SEXP , SEXP , 
                               SEXP , SEXP , 
                               SEXP );
extern SEXP _PAFit_generate_net_C_with_count_multi_corrected(SEXP, SEXP, SEXP , SEXP , SEXP, SEXP , SEXP , SEXP , SEXP,
                                                               SEXP , SEXP , SEXP, SEXP,SEXP,SEXP); 

  
static const R_CallMethodDef CallEntries[] = {
  {"_PAFit_cal_var_f_new", (DL_FUNC) &_PAFit_cal_var_f_new,10},
  {"_PAFit_cal_var_f", (DL_FUNC) &_PAFit_cal_var_f, 9},
  {"_PAFit_coeff_var", (DL_FUNC) &_PAFit_coeff_var, 6},
  {"_PAFit_coeff_theta", (DL_FUNC) &_PAFit_coeff_theta, 5},
  {"_PAFit_var_alpha", (DL_FUNC) &_PAFit_var_alpha, 11},
  {"_PAFit_update_alpha_fast", (DL_FUNC) &_PAFit_update_alpha_fast, 11},
  {"_PAFit_update_f_alpha_new", (DL_FUNC) &_PAFit_update_f_alpha_new, 12},
  {"_PAFit_update_f_alpha", (DL_FUNC) &_PAFit_update_f_alpha, 11},
  {"_PAFit_update_offset", (DL_FUNC) &_PAFit_update_offset, 7},
  {"_PAFit_update_f", (DL_FUNC) &_PAFit_update_f, 10},
  {"_PAFit_get_stats", (DL_FUNC) &_PAFit_get_stats, 23},
  {"_PAFit_normalized_constant_alpha", (DL_FUNC) &_PAFit_normalized_constant_alpha, 8},
  {"_PAFit_normalized_constant", (DL_FUNC) &_PAFit_normalized_constant, 6},
  {"_PAFit_update_f_new", (DL_FUNC) &_PAFit_update_f_new, 11},
  {"_PAFit_generate_net_C_with_count_multi_corrected", (DL_FUNC) &_PAFit_generate_net_C_with_count_multi_corrected,15},
  {NULL, NULL, 0}
};

void R_init_PAFit(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}