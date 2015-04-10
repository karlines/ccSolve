#include <R.h>
#include <Rdefines.h>
#include "ccSolve.h"

/* based on minpack.lm package file nls_lm.c created by Katharina Mullen
 * adapted and slightly simplified for use with compiled code by 
 * Karline Soetaert */

/* global declarations */
typedef struct nlsopt_struct {
  double ftol;
  double ptol;
  double gtol;
  double epsfcn;
  double *diag;
  double factor;
  double *lower;
  double *upper;
  int nprint;
  int maxiter;
  int niter;
  double rsstrace[1024];
} nlsopt_struct, *nlsOptStruct;

/* globals */

nlsOptStruct OS;
double *rpar, *data;
int *ipar, nrowdat, ncoldat; 

/* use with ccnls 
 * func(npar, ndat,  x, f, rpar, ipar);
 * jac(n, ndat, ldfjac, x, df, rpar, ipar )
*/
typedef void C_fnls_type(int*, int*, double*, double*, double *, int*);
typedef void C_jnls_type(int*, int*, int *, double*, double*, double *, int*);

C_fnls_type *fnls = NULL;
C_jnls_type *jnls = NULL;

/* call when jacobian function is not given by user */

void cc_fcn_lmdif(int *m, int *n, double *par, double *fvec, int *iflag)
{
    int i;
    double sumf;

    /* Rprintf("fcn-lmdif calling...\n"); */
    for (i = 0; i < *n; i++) {
      if (!R_FINITE(par[i]))
	       error("non-finite value supplied by lmdif!");
	    if (par[i] > OS->upper[i])
        par[i] = OS->upper[i];
      if (par[i] < OS->lower[i])
        par[i] = OS->lower[i];
    }  
    if(*iflag == 0){ 
      if(OS->nprint > 0) {
	      Rprintf("It. %4d, RSS = %10g, Par. =", OS->niter, 
		    OS->rsstrace[OS->niter]);
         for (i = 0; i < *n; i++)
          Rprintf(" % 10g", par[i]);
          Rprintf("\n");
      }
      OS->niter++;
    }
    else if (*iflag == 1 || *iflag == 2) {
     fnls( n, m, par, fvec, rpar, ipar );
	   sumf = 0;
	   for (i = 0; i < *m; i++) 
	     sumf += (fvec[i]*fvec[i]);
	   OS->rsstrace[OS->niter] = sumf;
	
    }
    if(OS->niter == OS->maxiter)
      *iflag = -1;
    
}

/* call when jacobian function is given by user */

void cc_fcn_lmder(int *m, int *n, double *par, double *fvec, double *fjac, 
	       int *ldfjac, int *iflag)
{
  int i; 
  double sumf;
  
  /* Rprintf("fcn-lmder calling...\n"); */
  for (i = 0; i < *n; i++) {
      if (!R_FINITE(par[i]))
	       error("non-finite value supplied by lmder!");
	    if (par[i] > OS->upper[i])
        par[i] = OS->upper[i];
      if (par[i] < OS->lower[i])
        par[i] = OS->lower[i];
    }
  
  if(*iflag == 0) {
    if(OS->nprint > 0) {
      Rprintf("It. %4d, RSS = %10g, Par. =", OS->niter, 
	      OS->rsstrace[OS->niter]);
      for (i = 0; i < *n; i++)
   	   Rprintf(" % 10g", par[i]);
       Rprintf("\n");
      
    }
    OS->niter++;
  }
  else if (*iflag == 1) {
    fnls( n, m, par, fvec, rpar, ipar );

    sumf = 0;
    for (i = 0; i < *m; i++) {
      sumf += fvec[i]*fvec[i];
    }
    OS->rsstrace[OS->niter] = sumf;
  }
  else if (*iflag == 2) { /* jacobian */
    jnls( n, m, ldfjac, par, fjac, rpar, ipar );
  }
  
  if(OS->niter == OS->maxiter)
    *iflag = -1;
}


/* initialiser function called by user */
void Initstdat (int *N, double *pdata)
{
  int i;

  if ((*N) > ncoldat)
    {
      warning("Number of data passed to solver, %i; number in DLL, %i\n",
      ncoldat, *N);
      PROBLEM "Confusion over the columns of data"
      ERROR;
    } 
  else
    { 
      for (i = 0; i < nrowdat; i++) pdata[i] = data[i + (*N-1)*nrowdat];
    }      
}


/* Main function */

SEXP cc_nls_lm(SEXP Par, SEXP fn, SEXP jac, SEXP initfunc, 
  SEXP control, SEXP Lower, SEXP Upper, SEXP Ndat, SEXP Ncol, SEXP Dat, 
  SEXP Rpar, SEXP Ipar)
{
    int     i, j;
    int     n, m, ldfjac;
    int     info, nfev, njev, np;

    double  *par, *fvec, *fjac, *qtf, *wa1, *wa2, *wa3, *wa4,
      *perm, *perm_t, *r, *r2, *r2_x_perm_t, *hess;
    int     *ipvt;

//    SEXP    eval_test;
    SEXP    sexp_par, sexp_fvec, sexp_diag, sexp_hess, sexp_info, sexp_niter, 
      sexp_message, sexp_rsstrace;
    SEXP    out, out_names;
    char    lmfun_name[8], message[256];
    int     maxfev;
    int     mode;
    int     npr = 1;
    C_init_dat_type *initializer;

    PROTECT_INDEX ipx;

    OS = (nlsOptStruct) R_alloc(1, sizeof(nlsopt_struct));

    n = LENGTH(Par);
    par = (double *) R_alloc(n, sizeof(double));
    for (i = 0; i < n; i++)
      par[i] = REAL(Par)[i];

    OS->upper = (double *) R_alloc(n, sizeof(double));
    for (i = 0; i < n; i++)
      OS->upper[i] = REAL(Upper)[i];

    OS->lower = (double *) R_alloc(n, sizeof(double));
    for (i = 0; i < n; i++)
      OS->lower[i] = REAL(Lower)[i];

    np = LENGTH(Ipar);
    ipar = (int *) R_alloc(np, sizeof(int));
    for (i = 0; i < np; i++)
      ipar[i] = INTEGER(Ipar)[i];

    np = LENGTH(Rpar);
    rpar = (double *) R_alloc(np, sizeof(double));
    for (i = 0; i < np; i++)
      rpar[i] = REAL(Rpar)[i];

    if (!inherits(fn, "NativeSymbol")) 
       error("'fn' is not a compiled function");
    fnls = (C_fnls_type *)  R_ExternalPtrAddr(fn);  
    if (!inherits(initfunc, "NativeSymbol")) 
      error("'initfunc' is not a compiled function");

    m = INTEGER(Ndat)[0];
    nrowdat = m;
    ncoldat = INTEGER(Ncol)[0];
    np = LENGTH (Dat);  
    data = (double *) R_alloc(np, sizeof(double));
    for (i = 0; i < np; i++)
      data[i] = REAL(Dat)[i];    
    

/* initialiser function for data*/
    initializer = (C_init_dat_type *) R_ExternalPtrAddr(initfunc);
    initializer(&Initstdat, &m);

    ldfjac = m;    

    fvec        = real_vector(m);
    fjac        = real_vector(ldfjac * n);
    qtf         = real_vector(n);
    wa1         = real_vector(n);
    wa2         = real_vector(n);
    wa3         = real_vector(n);
    wa4         = real_vector(m);
    ipvt        =  int_vector(n);
    perm        = real_vector(n * n);
    perm_t      = real_vector(n * n);
    r           = real_vector(n * n);
    r2          = real_vector(n * n);
    r2_x_perm_t = real_vector(n * n);
    hess        = real_vector(n * n);

    OS->ftol    = NUMERIC_VALUE(getListElement(control, "ftol"));
    OS->ptol    = NUMERIC_VALUE(getListElement(control, "ptol"));
    OS->gtol    = NUMERIC_VALUE(getListElement(control, "gtol"));
    OS->epsfcn  = NUMERIC_VALUE(getListElement(control, "epsfcn"));
    OS->factor  = NUMERIC_VALUE(getListElement(control, "factor"));
    OS->diag    = real_vector(n);

    PROTECT_WITH_INDEX(sexp_diag = getListElement(control, "diag"), &ipx);
    switch (TYPEOF(sexp_diag)) {
    case REALSXP:
        if (length(sexp_diag) == n) {
            REPROTECT(sexp_diag = duplicate(sexp_diag), ipx);
            for (i = 0; i < n; i++)
                OS->diag[i] = NUMERIC_POINTER(sexp_diag)[i];
            mode = 2;
        }
        else {
            REPROTECT(sexp_diag = NEW_NUMERIC(n), ipx);
            mode = 1;
        }
        break;
    case VECSXP:
        if (length(sexp_diag) == n) {
            REPROTECT(sexp_diag = duplicate(sexp_diag), ipx);
            for (i = 0; i < n; i++) {
                SET_VECTOR_ELT(sexp_diag, i, AS_NUMERIC(VECTOR_ELT(sexp_diag, i)));
                OS->diag[i] = NUMERIC_VALUE(VECTOR_ELT(sexp_diag, i));
            }
            mode = 2;
        }
        else {
            REPROTECT(sexp_diag = NEW_LIST(n), ipx);
            for (i = 0; i < n; i++)
                SET_VECTOR_ELT(sexp_diag, i, NEW_NUMERIC(1));
            mode = 1;
        }
        break;
    default:
        error("`diag' that you provided is non-list and non-numeric!");
    }
    maxfev = INTEGER_VALUE(getListElement(control, "maxfev"));
    OS->maxiter = INTEGER_VALUE(getListElement(control, "maxiter"));
    if(OS->maxiter > 1024) {
      OS->maxiter = 1024;
      warning("resetting `maxiter' to 1024!");
    }
    OS->nprint = INTEGER_VALUE(getListElement(control, "nprint"));
    if(OS->nprint > 0)
      OS->nprint = 1;

    OS->niter = 0;

    if (!inherits(jac, "NativeSymbol")) {
        F77_CALL(lmdif)(&cc_fcn_lmdif, &m, &n, par, fvec,
                        &OS->ftol, &OS->ptol, &OS->gtol,
                        &maxfev, &OS->epsfcn, OS->diag, &mode,
                        &OS->factor, &npr, &info, &nfev, fjac, &ldfjac,
                         ipvt, qtf, wa1, wa2, wa3, wa4);
        strcpy(lmfun_name, "lmdif");
    }
    else {
        jnls = (C_jnls_type *)  R_ExternalPtrAddr(fn);  
        F77_CALL(lmder)(&cc_fcn_lmder, &m, &n, par, fvec,
                         fjac, &ldfjac,
                        &OS->ftol, &OS->ptol, &OS->gtol,
                        &maxfev, OS->diag, &mode,
                        &OS->factor, &npr, &info, &nfev, &njev,
                         ipvt, qtf, wa1, wa2, wa3, wa4);
        strcpy(lmfun_name, "lmder");
   }   
/*========================================================================*/
    
    fcn_message(message, info, maxfev, OS->maxiter);
    if (info < 1 || 9 < info)
      warning("%s: info = %d. %s\n\n", lmfun_name, info, message);
    
    PROTECT(sexp_hess = NEW_NUMERIC(n*n));
    for (j = 0; j < n; j++)
        for (i = 0; i < n; i++) {
            perm[j*n + i] = (i + 1 == ipvt[j]) ? 1 : 0;
            r   [j*n + i] = (i <= j) ? fjac[j*ldfjac + i] : 0;
        }

    /*  perm %*% t(r) %*% r %*% t(perm)  *
     *    |       |___r2__|         |    *
     *    |           |_r2_x_perm_t_|    *
     *    |_______hess_______|           */

    transpose(perm, n, n, perm_t);
    crossprod(r,    n, n, r,           n, n, r2);
    matprod  (r2,   n, n, perm_t,      n, n, r2_x_perm_t);
    matprod  (perm, n, n, r2_x_perm_t, n, n, hess);

    for (i = 0; i < n*n; i++)
        NUMERIC_POINTER(sexp_hess)[i] = hess[i];

    
    PROTECT(sexp_par = NEW_NUMERIC(n));
    for (i = 0; i < n; i++)
        NUMERIC_POINTER(sexp_par)[i] = par[i];

    PROTECT(sexp_fvec = NEW_NUMERIC(m));
    for (i = 0; i < m; i++)
        NUMERIC_POINTER(sexp_fvec)[i] = fvec[i];

    PROTECT(sexp_rsstrace = NEW_NUMERIC(OS->niter));
    for (i = 0; i < OS->niter; i++)
      NUMERIC_POINTER(sexp_rsstrace)[i] = OS->rsstrace[i];

    PROTECT(sexp_info = NEW_INTEGER(1));
    INTEGER_POINTER(sexp_info)[0] = info;
    
    PROTECT(sexp_niter = NEW_INTEGER(1));
    INTEGER_POINTER(sexp_niter)[0] = OS->niter-1;

    PROTECT(sexp_message = NEW_STRING(1));
    SET_STRING_ELT(sexp_message, 0, mkChar(message));
    if (IS_NUMERIC(sexp_diag)) {
        for (i = 0; i < n; i++)
            NUMERIC_POINTER(sexp_diag)[i] = OS->diag[i];
    }
    else {
      for (i = 0; i < n; i++)
	     NUMERIC_POINTER(VECTOR_ELT(sexp_diag, i))[0] = OS->diag[i];
    }
    
    
    PROTECT(out = NEW_LIST(8));
    SET_VECTOR_ELT(out, 0, sexp_par);
    SET_VECTOR_ELT(out, 1, sexp_hess);
    SET_VECTOR_ELT(out, 2, sexp_fvec);
    SET_VECTOR_ELT(out, 3, sexp_info);
    SET_VECTOR_ELT(out, 4, sexp_message);
    SET_VECTOR_ELT(out, 5, sexp_diag);
    SET_VECTOR_ELT(out, 6, sexp_niter); 
    SET_VECTOR_ELT(out, 7, sexp_rsstrace);
        
    PROTECT(out_names = NEW_STRING(8));
    SET_STRING_ELT(out_names, 0, mkChar("par"));
    SET_STRING_ELT(out_names, 1, mkChar("hessian"));
    SET_STRING_ELT(out_names, 2, mkChar("fvec"));
    SET_STRING_ELT(out_names, 3, mkChar("info"));
    SET_STRING_ELT(out_names, 4, mkChar("message"));
    SET_STRING_ELT(out_names, 5, mkChar("diag"));
    SET_STRING_ELT(out_names, 6, mkChar("niter"));
    SET_STRING_ELT(out_names, 7, mkChar("rsstrace"));
        
    SET_NAMES(out, out_names);

    UNPROTECT(10);

    i = -99;
    initializer(&Initstdat, &i); /* free allocated memory (F95) */

    return out;
}


/* Function that returns a function evaluation */
SEXP cc_nls_func(SEXP Par, SEXP fn, SEXP initfunc, 
  SEXP Ndat, SEXP Ncol, SEXP Dat, SEXP Rpar, SEXP Ipar)
{
    int     i, n, m, np;

    double  *par, *fvec;
    SEXP    sexp_fvec;
    C_init_dat_type *initializer;

    n = LENGTH(Par);
    par = (double *) R_alloc(n, sizeof(double));
    for (i = 0; i < n; i++)
      par[i] = REAL(Par)[i];

    np = LENGTH(Ipar);
    ipar = (int *) R_alloc(np, sizeof(int));
    for (i = 0; i < np; i++)
      ipar[i] = INTEGER(Ipar)[i];

    np = LENGTH(Rpar);
    rpar = (double *) R_alloc(np, sizeof(double));
    for (i = 0; i < np; i++)
      rpar[i] = REAL(Rpar)[i];

    if (!inherits(fn, "NativeSymbol")) 
       error("'fn' is not a compiled function");
    fnls = (C_fnls_type *)  R_ExternalPtrAddr(fn);  
    if (!inherits(initfunc, "NativeSymbol")) 
      error("'initfunc' is not a compiled function");

    m = INTEGER(Ndat)[0];
    nrowdat = m;
    ncoldat = INTEGER(Ncol)[0];
    np = LENGTH (Dat);  
    data = (double *) R_alloc(np, sizeof(double));
    for (i = 0; i < np; i++)
      data[i] = REAL(Dat)[i];    
    
/* initialiser function for data*/
    initializer = (C_init_dat_type *) R_ExternalPtrAddr(initfunc);
    initializer(&Initstdat, &m);

/* function call */
    fvec        = real_vector(m);
    fnls(&n, &m, par, fvec, rpar, ipar );

    PROTECT(sexp_fvec = NEW_NUMERIC(m));
    for (i = 0; i < m; i++)
        NUMERIC_POINTER(sexp_fvec)[i] = fvec[i];

    i = -99;
    initializer(&Initstdat, &i); /* free allocated memory (F95) */

    UNPROTECT(1);

    return sexp_fvec;
}
