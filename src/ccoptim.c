/*
 This is an adaptation of the file "optim.c", optimize.c and "zeroin.c" that is 
 part of the R-core software.

 It has been adapted to work with problems defined in compiled code
 by Karline Soetaert.
  
 .. removed underscore in error(_(
 .. changed .External2 calls to .Call

* original header of optim.c (and same for the other files):
 *
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1999-2013  The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */


#include <R.h>
#include <Rdefines.h>
#include "ccSolve.h"

#define EPSILON DBL_EPSILON

/* global variables */
C_func_type *fcall = NULL;
C_func_type2 *fcall2 = NULL;
C_jac_type *gcall = NULL;

double* parcopy;
double* grads;
int hasgn; 
double *rpar, *data;
int *ipar, nrowdat, ncoldat; 

SEXP getListElement(SEXP list, char *str)
{
    SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
    int i;

    for (i = 0; i < length(list); i++)
	if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	    elmt = VECTOR_ELT(list, i);
	    break;
	}
    return elmt;
}

static double * vect(int n)
{
    return (double *)R_alloc(n, sizeof(double));
}


typedef struct opt_struct
{
    SEXP R_fcall;    /* function */
    SEXP R_gcall;    /* gradient */
    SEXP R_env;      /* where to evaluate the calls */
    double* ndeps;   /* tolerances for numerical derivatives */
    double fnscale;  /* scaling for objective */
    double* parscale;/* scaling for parameters */
    int usebounds;
    double* lower, *upper;
    SEXP names;	     /* names for par */
} opt_struct, *OptStruct;

/* compiled code version of original fminfn */
static double fminfn_cc(int n, double *p, void *ex)
{
    int i;
    double val;

    OptStruct OS = (OptStruct) ex;

    for (i = 0; i < n; i++) {
  	  if (!R_FINITE(p[i])) error("non-finite value supplied by optim");
	    parcopy[i] = p[i] * (OS->parscale[i]);
    }

    fcall(&n, parcopy, &val, rpar, ipar);    


    val = val/(OS->fnscale);
    return val;
}

/* compiled code version of fmingr */
static void fmingr_cc(int n, double *p, double *df, void *ex)
{
  int i;
  double val1, val2, eps, epsused, tmp;
  OptStruct OS = (OptStruct) ex;

  if (hasgn == 1) { /* analytical derivatives */

	  for (i = 0; i < n; i++) {
	    if (!R_FINITE(p[i])) error("non-finite value supplied by optim");
	    parcopy[i] = p[i] * (OS->parscale[i]);
	  }

    gcall(&n, parcopy, grads, rpar, ipar);    

	  for (i = 0; i < n; i++)
	    df[i] = grads[i] * (OS->parscale[i])/(OS->fnscale);

  } else { /* numerical derivatives */
	  for (i = 0; i < n; i++) {
	    if (!R_FINITE(p[i])) error("non-finite value supplied by optim");
	    parcopy[i] = p[i] * (OS->parscale[i]);
	  }

  
  	if(OS->usebounds == 0) {
	    for (i = 0; i < n; i++) {
  		eps = OS->ndeps[i];

		  parcopy[i] = (p[i] + eps) * (OS->parscale[i]);
	    fcall(&n, parcopy, &val1, rpar, ipar);    
   		val1 = val1/(OS->fnscale);

		  parcopy[i] = (p[i] - eps) * (OS->parscale[i]);
	    fcall(&n, parcopy, &val2, rpar, ipar);    
   		val2 = val2/(OS->fnscale);

  		df[i] = (val1 - val2)/(2 * eps);
		  if(!R_FINITE(df[i]))
		    error(("non-finite finite-difference value [%d]"), i+1);
		  parcopy[i] = p[i] * (OS->parscale[i]);
	    }
	  } else { /* usebounds */
	    for (i = 0; i < n; i++) {
		    epsused = eps = OS->ndeps[i];
		    tmp = p[i] + eps;
		    if (tmp > OS->upper[i]) {
		      tmp = OS->upper[i];
		      epsused = tmp - p[i];
		    }
		   parcopy[i] = tmp * (OS->parscale[i]);
  	   fcall(&n, parcopy, &val1, rpar, ipar);    
   		 val1 = val1/(OS->fnscale);

  		 tmp = p[i] - eps;
		   if (tmp < OS->lower[i]) {
		    tmp = OS->lower[i];
		    eps = p[i] - tmp;
	   	 }
		   parcopy[i] = tmp * (OS->parscale[i]);
  	   fcall(&n, parcopy, &val2, rpar, ipar);    
   		 val2 = val2/(OS->fnscale);


   		 df[i] = (val1 - val2)/(epsused + eps);
		   if(!R_FINITE(df[i]))
		    error(("non-finite finite-difference value [%d]"), i+1);
		  parcopy[i] = p[i] * (OS->parscale[i]);
	    }
	  }
  }
}


/* par fn gr method options upper, lower*/                              
SEXP call_optim(SEXP par, SEXP fn, SEXP gr, SEXP initfunc, SEXP method, SEXP options, 
                SEXP slower, SEXP supper, SEXP Ndat, SEXP Ncol, SEXP Dat, SEXP Rpar, SEXP Ipar)
{
    SEXP tmp;
    SEXP res, value, counts, conv;
    int i, np, npar=0, *mask, trace, maxit, fncount = 0, grcount = 0, nREPORT, tmax;
    int ifail = 0;
    double *dpar, *opar, val = 0.0, abstol, reltol, temp;
    const char *tn;
    C_init_dat_type *initializer = NULL;

    OptStruct OS;

    OS = (OptStruct) R_alloc(1, sizeof(opt_struct));
    OS->usebounds = 0;

    OS->names = getAttrib(par, R_NamesSymbol);

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
    fcall = (C_func_type *) R_ExternalPtrAddr(fn);  

    if (! isNull(initfunc)) { /* data are passed with initfunc to inialise them*/
      if (!inherits(initfunc, "NativeSymbol")) 
        error("'initfunc' is not a compiled function");
      nrowdat = INTEGER(Ndat)[0];
      ncoldat = INTEGER(Ncol)[0];
      np = LENGTH (Dat);  
      data = (double *) R_alloc(np, sizeof(double));
      for (i = 0; i < np; i++)
       data[i] = REAL(Dat)[i];    
    

/* initialiser function for data*/
      initializer = (C_init_dat_type *) R_ExternalPtrAddr(initfunc);
      initializer(&Initstdat, &nrowdat);
    }
    if (!isString(method)|| LENGTH(method) != 1)
  	  error("invalid '%s' argument", "method");
    tn = CHAR(STRING_ELT(method, 0));

    npar = LENGTH(par);
    dpar = vect(npar);
    opar = vect(npar);
    
    /* karline : global pars */
    parcopy = vect(npar);
    grads = vect(npar);
    
    trace = asInteger(getListElement(options, "trace"));
    OS->fnscale = asReal(getListElement(options, "fnscale"));
    tmp = getListElement(options, "parscale");
    if (LENGTH(tmp) != npar)
  	  error("'parscale' is of the wrong length");
    
    PROTECT(tmp = coerceVector(tmp, REALSXP));
    OS->parscale = vect(npar);
    for (i = 0; i < npar; i++) OS->parscale[i] = REAL(tmp)[i];
    UNPROTECT(1);
    
    for (i = 0; i < npar; i++)
	    dpar[i] = REAL(par)[i] / (OS->parscale[i]);
    PROTECT(res = allocVector(VECSXP, 5));
    SEXP names;
    PROTECT(names = allocVector(STRSXP, 5));
    SET_STRING_ELT(names, 0, mkChar("par"));
    SET_STRING_ELT(names, 1, mkChar("value"));
    SET_STRING_ELT(names, 2, mkChar("counts"));
    SET_STRING_ELT(names, 3, mkChar("convergence"));
    SET_STRING_ELT(names, 4, mkChar("message"));
    setAttrib(res, R_NamesSymbol, names);
    UNPROTECT(1);
    
    PROTECT(value = allocVector(REALSXP, 1));
    PROTECT(counts = allocVector(INTSXP, 2));
    SEXP countnames;
    PROTECT(countnames = allocVector(STRSXP, 2));
    SET_STRING_ELT(countnames, 0, mkChar("function"));
    SET_STRING_ELT(countnames, 1, mkChar("gradient"));
    setAttrib(counts, R_NamesSymbol, countnames);
    UNPROTECT(1);
    PROTECT(conv = allocVector(INTSXP, 1));
    abstol = asReal(getListElement(options, "abstol"));
    reltol = asReal(getListElement(options, "reltol"));
    maxit = asInteger(getListElement(options, "maxit"));
    if (maxit == NA_INTEGER) error("'maxit' is not an integer");

    if (strcmp(tn, "Nelder-Mead") != 0) {
    
     if (!isNull(gr)) {
       hasgn = 1; 
       if (!inherits(gr, "NativeSymbol")) 
        error("'gr' is not a compiled function");
       gcall = (C_jac_type *) R_ExternalPtrAddr(gr);  
	   } else {
       hasgn = 0;
 	   }
   }
    if (strcmp(tn, "Nelder-Mead") == 0) {
	double alpha, beta, gamm;

	alpha = asReal(getListElement(options, "alpha"));
	beta = asReal(getListElement(options, "beta"));
	gamm = asReal(getListElement(options, "gamma"));
	nmmin(npar, dpar, opar, &val, fminfn_cc, &ifail, abstol, reltol,
	      (void *)OS, alpha, beta, gamm, trace, &fncount, maxit);
	for (i = 0; i < npar; i++)
	    REAL(par)[i] = opar[i] * (OS->parscale[i]);
	grcount = NA_INTEGER;

    }
    else if (strcmp(tn, "SANN") == 0) {
	tmax = asInteger(getListElement(options, "tmax"));
	temp = asReal(getListElement(options, "temp"));
	if (trace) trace = asInteger(getListElement(options, "REPORT"));
	if (tmax == NA_INTEGER || tmax < 1) // PR#15194
	    error("'tmax' is not a positive integer");

  PROTECT(OS->R_gcall = R_NilValue); /* to generate new values... */

	samin (npar, dpar, &val, fminfn_cc, maxit, tmax, temp, trace, (void *)OS);
  UNPROTECT(1);
  for (i = 0; i < npar; i++)
	    REAL(par)[i] = dpar[i] * (OS->parscale[i]);
	fncount = npar > 0 ? maxit : 1;
	grcount = NA_INTEGER;

    } else if (strcmp(tn, "BFGS") == 0) {
	SEXP ndeps;

	nREPORT = asInteger(getListElement(options, "REPORT"));
	if (isNull(gr)) {
	    ndeps = getListElement(options, "ndeps");
	    if (LENGTH(ndeps) != npar)
		error("'ndeps' is of the wrong length");
	    OS->ndeps = vect(npar);
	    PROTECT(ndeps = coerceVector(ndeps, REALSXP));
	    for (i = 0; i < npar; i++) OS->ndeps[i] = REAL(ndeps)[i];
	    UNPROTECT(1);
	}
	mask = (int *) R_alloc(npar, sizeof(int));
	for (i = 0; i < npar; i++) mask[i] = 1;
	vmmin(npar, dpar, &val, fminfn_cc, fmingr_cc, maxit, trace, mask, abstol,
	      reltol, nREPORT, (void *)OS, &fncount, &grcount, &ifail);
	for (i = 0; i < npar; i++)
	    REAL(par)[i] = dpar[i] * (OS->parscale[i]);
    } else if (strcmp(tn, "CG") == 0) {
	int type;
	SEXP ndeps;

	type = asInteger(getListElement(options, "type"));
	if (isNull(gr)) {
	    ndeps = getListElement(options, "ndeps");
	    if (LENGTH(ndeps) != npar)
		error("'ndeps' is of the wrong length");
	    OS->ndeps = vect(npar);
	    PROTECT(ndeps = coerceVector(ndeps, REALSXP));
	    for (i = 0; i < npar; i++) OS->ndeps[i] = REAL(ndeps)[i];
	    UNPROTECT(1);
	}
	cgmin(npar, dpar, opar, &val, fminfn_cc, fmingr_cc, &ifail, abstol,
	      reltol, (void *)OS, type, trace, &fncount, &grcount, maxit);
	for (i = 0; i < npar; i++)
	    REAL(par)[i] = opar[i] * (OS->parscale[i]);

    } else if (strcmp(tn, "L-BFGS-B") == 0) {
	SEXP ndeps, smsg;
	double *lower = vect(npar), *upper = vect(npar);
	int lmm, *nbd = (int *) R_alloc(npar, sizeof(int));
	double factr, pgtol;
	char msg[60];

	nREPORT = asInteger(getListElement(options, "REPORT"));
	factr = asReal(getListElement(options, "factr"));
	pgtol = asReal(getListElement(options, "pgtol"));
	lmm = asInteger(getListElement(options, "lmm"));
	if (isNull(gr)) {
	    ndeps = getListElement(options, "ndeps");
	    if (LENGTH(ndeps) != npar)
		error("'ndeps' is of the wrong length");
	    OS->ndeps = vect(npar);
	    PROTECT(ndeps = coerceVector(ndeps, REALSXP));
	    for (i = 0; i < npar; i++) OS->ndeps[i] = REAL(ndeps)[i];
	    UNPROTECT(1);
	}

	for (i = 0; i < npar; i++) {
	    lower[i] = REAL(slower)[i] / (OS->parscale[i]);
	    upper[i] = REAL(supper)[i] / (OS->parscale[i]);
	    if (!R_FINITE(lower[i])) {
		if (!R_FINITE(upper[i])) nbd[i] = 0; else nbd[i] = 3;
	    } else {
		if (!R_FINITE(upper[i])) nbd[i] = 1; else nbd[i] = 2;
	    }
	}
	OS->usebounds = 1;
	OS->lower = lower;
	OS->upper = upper;
	lbfgsb(npar, lmm, dpar, lower, upper, nbd, &val, fminfn_cc, fmingr_cc,
	       &ifail, (void *)OS, factr, pgtol, &fncount, &grcount,
	       maxit, msg, trace, nREPORT);
	for (i = 0; i < npar; i++)
	    REAL(par)[i] = dpar[i] * (OS->parscale[i]);
	PROTECT(smsg = mkString(msg));
	SET_VECTOR_ELT(res, 4, smsg);
	UNPROTECT(1);
    } else
	error("unknown 'method'");

    if (! isNull(initfunc)) {   /* free allocated memory */
      np = -99;  
      initializer(&Initstdat, &np);
    }

    if(!isNull(OS->names)) setAttrib(par, R_NamesSymbol, OS->names);
    REAL(value)[0] = val * (OS->fnscale);
    SET_VECTOR_ELT(res, 0, par); SET_VECTOR_ELT(res, 1, value);
    INTEGER(counts)[0] = fncount; INTEGER(counts)[1] = grcount;
    SET_VECTOR_ELT(res, 2, counts);
    INTEGER(conv)[0] = ifail;
    SET_VECTOR_ELT(res, 3, conv);
    UNPROTECT(4);
    return res;
}



/* function call */                              
SEXP call_optim_func(SEXP par, SEXP fn, SEXP initfunc, 
    SEXP Ndat, SEXP Ncol, SEXP Dat, SEXP Rpar, SEXP Ipar)
{
    SEXP value;
    int i, np, npar=0;
    double *dpar, val = 0.0;
    C_init_dat_type *initializer = NULL;

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
    fcall = (C_func_type *) R_ExternalPtrAddr(fn);  

    if (! isNull(initfunc)) {   /* data are passed with initfunc to inialise them*/
      if (!inherits(initfunc, "NativeSymbol")) 
        error("'initfunc' is not a compiled function");
      nrowdat = INTEGER(Ndat)[0];
      ncoldat = INTEGER(Ncol)[0];
      np = LENGTH (Dat);  
      data = (double *) R_alloc(np, sizeof(double));
      for (i = 0; i < np; i++)
       data[i] = REAL(Dat)[i];    
    
/* initialiser function for data*/
      initializer = (C_init_dat_type *) R_ExternalPtrAddr(initfunc);
      initializer(&Initstdat, &nrowdat);
    }

    npar = LENGTH(par);
    dpar = vect(npar);
    for (i = 0; i < npar; i++)
	    dpar[i] = REAL(par)[i] ;

    fcall(&npar, dpar, &val, rpar, ipar);    

    if (! isNull(initfunc)) {   /* free allocated memory */
      np = -99;  
      initializer(&Initstdat, &np);
    }
    PROTECT(value = NEW_NUMERIC(1));

    REAL(value)[0] = val;
    UNPROTECT(1);
    return value;
}


/* par fn gr options */ 
SEXP call_optimhess(SEXP par, SEXP fn, SEXP gr, SEXP options, SEXP Rpar, SEXP Ipar)
{
    SEXP tmp, ndeps, ans;
    OptStruct OS;
    int npar, np, i , j;
    double *dpar, *df1, *df2, eps;

    OS = (OptStruct) R_alloc(1, sizeof(opt_struct));
    OS->usebounds = 0;

    npar = LENGTH(par);
    OS->names = getAttrib(par, R_NamesSymbol);

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
    fcall = (C_func_type *) R_ExternalPtrAddr(fn);  

    OS->fnscale = asReal(getListElement(options, "fnscale"));
    tmp = getListElement(options, "parscale");
    
    if (LENGTH(tmp) != npar)
 	    error("'parscale' is of the wrong length");

    PROTECT(tmp = coerceVector(tmp, REALSXP));
    OS->parscale = vect(npar);
    for (i = 0; i < npar; i++) OS->parscale[i] = REAL(tmp)[i];
    UNPROTECT(1);

    /* karline : global pars */
    parcopy = vect(npar);
    grads = vect(npar);

    if (!isNull(gr)) {
      hasgn = 1; 
      if (!inherits(gr, "NativeSymbol")) 
        error("'gr' is not a compiled function");
      gcall = (C_jac_type *) R_ExternalPtrAddr(gr);  
	  } else {
      hasgn = 0;
 	  }

    ndeps = getListElement(options, "ndeps");
    if (LENGTH(ndeps) != npar) error("'ndeps' is of the wrong length");
    OS->ndeps = vect(npar);
    PROTECT(ndeps = coerceVector(ndeps, REALSXP));
    for (i = 0; i < npar; i++) OS->ndeps[i] = REAL(ndeps)[i];
    UNPROTECT(1);
    PROTECT(ans = allocMatrix(REALSXP, npar, npar));
    dpar = vect(npar);
    for (i = 0; i < npar; i++)
	   dpar[i] = REAL(par)[i] / (OS->parscale[i]);
    df1 = vect(npar);
    df2 = vect(npar);
    for (i = 0; i < npar; i++) {
	   eps = OS->ndeps[i]/(OS->parscale[i]);
     dpar[i] = dpar[i] + eps;
	   fmingr_cc(npar, dpar, df1, (void *)OS);
	   dpar[i] = dpar[i] - 2 * eps;
	   fmingr_cc(npar, dpar, df2, (void *)OS);
	   for (j = 0; j < npar; j++)
	    REAL(ans)[i * npar + j] = (OS->fnscale) * (df1[j] - df2[j])/
		(2 * eps * (OS->parscale[i]) * (OS->parscale[j]));
	   dpar[i] = dpar[i] + eps;
    }
    // now symmetrize
    for (i = 0; i < npar; i++) 
	   for (j = 0; j < i; j++) {
	    double tmp =
		0.5 * (REAL(ans)[i * npar + j] + REAL(ans)[j * npar + i]);
	    REAL(ans)[i * npar + j] = REAL(ans)[j * npar + i] = tmp;
	}
    SEXP nm = getAttrib(par, R_NamesSymbol);
    if(!isNull(nm)) {
	SEXP dm;
	PROTECT(dm = allocVector(VECSXP, 2));
	SET_VECTOR_ELT(dm, 0, duplicate(nm));
	SET_VECTOR_ELT(dm, 1, duplicate(nm));
	setAttrib(ans, R_DimNamesSymbol, dm);
	UNPROTECT(1);
    }
    UNPROTECT(1);
    return ans;
}
                     
                     
                     
                     
/*
 This is an adaptation of parts of the file "optimize.c" of the R
 software.
 It has been adapted to work with compiled code solutions 
 by Karline Soetaert.
*/


struct callinfo {
  SEXP R_fcall;
  SEXP R_env;
} ;

/*static SEXP R_fcall1;
  static SEXP R_env1; */

int maximize;

/* compiled code version of fcn1 for use with optimize NOTE: somewhat simpler */
static double cc_fcn1(double x, struct callinfo *info)
{
    double val;

	  if (!R_FINITE(x)) error("non-finite value supplied by optimize");

    fcall2(&x, &val, rpar, ipar);    
	  if (!R_FINITE(val)) {
      warning("NA/inf supplied by optimize, replaced by maximum positive value");
      val = DBL_MAX;
    }
    if (maximize) val = - val;
    return val;
}

/* compiled code version of fcn1 for use with uniroot NOTE: somewhat simpler */
static double cc_fcn2(double x, struct callinfo *info)
{
    double val;

	  if (!R_FINITE(x)) error("non-finite value supplied by uniroot");

    fcall2(&x, &val, rpar, ipar);    
	  if (!R_FINITE(val)) {
      error("NA/inf supplied by root - stopped");
      return 0;
    }
    if (maximize) val = - val;
    return val;
}
                     
                     
/* unfortunately this function has to be copied */                     
static
double Brent_fmin(double ax, double bx, double (*f)(double, void *),
		  void *info, double tol)
{
    /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;

    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

/*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;

    d = 0.;/* -Wall */
    e = 0.;
    fx = (*f)(x, info);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;

/*  main loop starts here ----------------------------------- */

    for(;;) {
	xm = (a + b) * .5;
	tol1 = eps * fabs(x) + tol3;
	t2 = tol1 * 2.;

	/* check stopping criterion */

	if (fabs(x - xm) <= t2 - (b - a) * .5) break;
	p = 0.;
	q = 0.;
	r = 0.;
	if (fabs(e) > tol1) { /* fit parabola */

	    r = (x - w) * (fx - fv);
	    q = (x - v) * (fx - fw);
	    p = (x - v) * q - (x - w) * r;
	    q = (q - r) * 2.;
	    if (q > 0.) p = -p; else q = -q;
	    r = e;
	    e = d;
	}

	if (fabs(p) >= fabs(q * .5 * r) ||
	    p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */

	    if (x < xm) e = b - x; else e = a - x;
	    d = c * e;
	}
	else { /* a parabolic-interpolation step */

	    d = p / q;
	    u = x + d;

	    /* f must not be evaluated too close to ax or bx */

	    if (u - a < t2 || b - u < t2) {
		d = tol1;
		if (x >= xm) d = -d;
	    }
	}

	/* f must not be evaluated too close to x */

	if (fabs(d) >= tol1)
	    u = x + d;
	else if (d > 0.)
	    u = x + tol1;
	else
	    u = x - tol1;

	fu = (*f)(u, info);

	/*  update  a, b, v, w, and x */

	if (fu <= fx) {
	    if (u < x) b = x; else a = x;
	    v = w;    w = x;   x = u;
	    fv = fw; fw = fx; fx = fu;
	} else {
	    if (u < x) a = u; else b = u;
	    if (fu <= fw || w == x) {
		v = w; fv = fw;
		w = u; fw = fu;
	    } else if (fu <= fv || v == x || v == w) {
		v = u; fv = fu;
	    }
	}
    }
    /* end of main loop */

    return x;
}

/* fmin(f, xmin, xmax tol, maximum) */
//cc_do_fmin", f, as.double(lower), as.double(upper), as.double(tol), as.integer(maximum))
SEXP cc_do_fmin(SEXP f, SEXP Lower, SEXP Upper, SEXP Tol, SEXP maximum, SEXP Rpar, SEXP Ipar)
{
    double xmin, xmax, tol;
    SEXP res;
    double val, x;
    int i, np;
    
    struct callinfo info;

    /* the function to be minimized */
    if (!inherits(f, "NativeSymbol")) 
       error("'f' is not a compiled function");
    fcall2 = (C_func_type2 *) R_ExternalPtrAddr(f);  

    np = LENGTH(Ipar);
    ipar = (int *) R_alloc(np, sizeof(int));
    for (i = 0; i < np; i++)
      ipar[i] = INTEGER(Ipar)[i];

    np = LENGTH(Rpar);
    rpar = (double *) R_alloc(np, sizeof(double));
    for (i = 0; i < np; i++)
      rpar[i] = REAL(Rpar)[i];

    /* xmin, xmax, tol */
    xmin = REAL(Lower)[0];
    if (!R_FINITE(xmin))
  	  error("invalid '%s' value", "xmin");
    xmax = REAL(Upper)[0];
    if (!R_FINITE(xmax))
	     error("invalid '%s' value", "xmax");
    if (xmin >= xmax)
	    error("'xmin' not less than 'xmax'");
    tol = REAL(Tol)[0];
    if (!R_FINITE(tol) || tol <= 0.0)
      error("invalid '%s' value", "tol");

    PROTECT(res = allocVector(REALSXP, 2));  /* karline: made 2 long; second = f-value*/
    x = Brent_fmin(xmin, xmax,
			      (double (*)(double, void*)) cc_fcn1, &info, tol);
    REAL(res)[0] = x;   
    fcall2(&x, &val, rpar, ipar);
    REAL(res)[1] = val;   
    UNPROTECT(1);    
			      
    return res;
}



// One Dimensional Root Finding --  just wrapper code for
// Brent's "zeroin"
// ---------------




/* zeroin2(f, ax, bx, f.ax, f.bx, tol, maxiter) */
SEXP cc_zeroin(SEXP f, SEXP ax, SEXP bx, SEXP Tol, SEXP maxiter, SEXP Rpar, SEXP Ipar)
{
    double f_ax, f_bx, val;
    double xmin, xmax, tol;
    int i, iter, np;
    SEXP res;
    struct callinfo info;

    /* the function to be minimized */
    if (!inherits(f, "NativeSymbol")) 
       error("'f' is not a compiled function");
    fcall2 = (C_func_type2 *) R_ExternalPtrAddr(f);  
           /* xmin, xmax */
    xmin = REAL(ax)[0];
    if (!R_FINITE(xmin)) error("invalid '%s' value", "xmin");
    xmax = REAL(bx)[0];
    if (!R_FINITE(xmax)) error("invalid '%s' value", "xmax");
    if (xmin >= xmax) error("'xmin' not less than 'xmax'");

    np = LENGTH(Ipar);
    ipar = (int *) R_alloc(np, sizeof(int));
    for (i = 0; i < np; i++)
      ipar[i] = INTEGER(Ipar)[i];

    np = LENGTH(Rpar);
    rpar = (double *) R_alloc(np, sizeof(double));
    for (i = 0; i < np; i++)
      rpar[i] = REAL(Rpar)[i];

    /* tol */
    tol = REAL(Tol)[0];
    if (!R_FINITE(tol) || tol <= 0.0) error("invalid '%s' value", "tol");
    /* maxiter */
    iter = INTEGER(maxiter)[0];
    if (iter <= 0) error("'maxiter' must be positive");
//    error("till here");

    /* f(ax) = f(xmin) *; f(bx) = f(xmax) */
    fcall2(&xmin, &f_ax, rpar, ipar);    
    fcall2(&xmax, &f_bx, rpar, ipar);    
    if (ISNA(f_ax)) error("NA value for '%s' is not allowed", "f.lower");
    if (ISNA(f_bx)) error("NA value for '%s' is not allowed", "f.upper");

    PROTECT(res = allocVector(REALSXP, 4));  /* one added for f()*/
	  val = R_zeroin2(xmin, xmax, f_ax, f_bx, (double (*)(double, void*)) cc_fcn2,
		 (void *) &info, &tol, &iter);
    REAL(res)[0] = val;
    REAL(res)[1] = (double)iter;
    REAL(res)[2] = tol;
    fcall2(&val, &f_ax, rpar, ipar);    
    REAL(res)[3] = f_ax;

    UNPROTECT(1);
    return res;
}




/* copy from R file zeroin.c */
/* R_zeroin2() is faster for "expensive" f(), in those typical cases where
 *             f(ax) and f(bx) are available anyway : */

double R_zeroin2(			/* An estimate of the root */
    double ax,				/* Left border | of the range	*/
    double bx,				/* Right border| the root is seeked*/
    double fa, double fb,		/* f(a), f(b) */
    double (*f)(double x, void *info),	/* Function under investigation	*/
    void *info,				/* Add'l info passed on to f	*/
    double *Tol,			/* Acceptable tolerance		*/
    int *Maxit)				/* Max # of iterations */
{
    double a,b,c, fc;			/* Abscissae, descr. see above,  f(c) */
    double tol;
    int maxit;

    a = ax;  b = bx;
    c = a;   fc = fa;
    maxit = *Maxit + 1; tol = * Tol;

    /* First test if we have found a root at an endpoint */
    if(fa == 0.0) {
	*Tol = 0.0;
	*Maxit = 0;
	return a;
    }
    if(fb ==  0.0) {
	*Tol = 0.0;
	*Maxit = 0;
	return b;
    }

    while(maxit--)		/* Main iteration loop	*/
    {
	double prev_step = b-a;		/* Distance from the last but one
					   to the last approximation	*/
	double tol_act;			/* Actual tolerance		*/
	double p;			/* Interpolation step is calcu- */
	double q;			/* lated in the form p/q; divi-
					 * sion operations is delayed
					 * until the last moment	*/
	double new_step;		/* Step at this iteration	*/

	if( fabs(fc) < fabs(fb) )
	{				/* Swap data for b to be the	*/
	    a = b;  b = c;  c = a;	/* best approximation		*/
	    fa=fb;  fb=fc;  fc=fa;
	}
	tol_act = 2*EPSILON*fabs(b) + tol/2;
	new_step = (c-b)/2;

	if( fabs(new_step) <= tol_act || fb == (double)0 )
	{
	    *Maxit -= maxit;
	    *Tol = fabs(c-b);
	    return b;			/* Acceptable approx. is found	*/
	}

	/* Decide if the interpolation can be tried	*/
	if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
	    && fabs(fa) > fabs(fb) ) {	/* and was in true direction,
					 * Interpolation may be tried	*/
	    register double t1,cb,t2;
	    cb = c-b;
	    if( a==c ) {		/* If we have only two distinct	*/
					/* points linear interpolation	*/
		t1 = fb/fa;		/* can only be applied		*/
		p = cb*t1;
		q = 1.0 - t1;
	    }
	    else {			/* Quadric inverse interpolation*/

		q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
		p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
		q = (q-1.0) * (t1-1.0) * (t2-1.0);
	    }
	    if( p>(double)0 )		/* p was calculated with the */
		q = -q;			/* opposite sign; make p positive */
	    else			/* and assign possible minus to	*/
		p = -p;			/* q				*/

	    if( p < (0.75*cb*q-fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
		&& p < fabs(prev_step*q/2) )	/* and isn't too large	*/
		new_step = p/q;			/* it is accepted
						 * If p/q is too large then the
						 * bisection procedure can
						 * reduce [b,c] range to more
						 * extent */
	}

	if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
	    if( new_step > (double)0 )	/* than tolerance		*/
		new_step = tol_act;
	    else
		new_step = -tol_act;
	}
	a = b;	fa = fb;			/* Save the previous approx. */
	b += new_step;	fb = (*f)(b, info);	/* Do step to a new approxim. */
	if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
	    /* Adjust c for it to have a sign opposite to that of b */
	    c = a;  fc = fa;
	}

    }
    /* failed! */
    *Tol = fabs(c-b);
    *Maxit = -1;
    return b;
}

