/*
 This is an adaptation of the file "integrate.c" that is 
 part of the R-core software.

 It has been adapted to work with problems defined in compiled code
 by Karline Soetaert.
  
 .. removed underscore in error(_(
 .. changed .External2 calls to .Call
 
 * original header of integrate.c
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 2001-2012  the R Core Team
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

/* global variables 
typedef void C_func_type(int*, double*, double*, double *, int*); */
C_func_type *fint = NULL;
double *rpar;
int *ipar; 


typedef struct int_struct
{
    SEXP f;    /* function */
    SEXP env;  /* where to evaluate the calls */
} int_struct, *IntStruct;


/* This is *the* ``integr_fn f'' : */
static void Rintfn(double *x, int n, void *ex)
{
    int i;
    double *f;

    f = (double *) R_alloc(n, sizeof(double));
    for (i = 0; i < n; i++)  f[i] = 0.;      

    fint(&n, x, f, rpar, ipar);

    for(i = 0; i < n; i++) {
      x[i] = f[i];
    	if(!R_FINITE(x[i])) error("non-finite function value");
    }
    return;
}


/* cc_call_dqags( ff, ,as.double(lower), as.double(upper),
	  		as.double(abs.tol), as.double(rel.tol), limit = limit, 
        as.double(rpar), as.integer(ipar)) */
SEXP cc_call_dqags(SEXP ff, SEXP slower, SEXP supper, SEXP Atol, SEXP Rtol,  SEXP Limit, SEXP Rpar, SEXP Ipar)
{
    int_struct is;
    SEXP ans, ansnames;
    double lower, upper, epsabs, epsrel, result, abserr, *work;
    int i, np, neval, ier, limit, lenw, last, *iwork;

    if (!inherits(ff, "NativeSymbol")) 
       error("'f' is not a compiled function");
    fint = (C_func_type *) R_ExternalPtrAddr(ff);  

    lower = REAL(slower)[0];
    upper = REAL(supper)[0];
    epsabs = REAL(Atol)[0];
    epsrel = REAL(Rtol)[0];
    limit = INTEGER(Limit)[0];
    lenw = 4 * limit;
    iwork = (int *) R_alloc((size_t) limit, sizeof(int));
    work = (double *) R_alloc((size_t) lenw, sizeof(double));

    np = LENGTH(Ipar);
    ipar = (int *) R_alloc(np, sizeof(int));
    for (i = 0; i < np; i++)
      ipar[i] = INTEGER(Ipar)[i];

    np = LENGTH(Rpar);
    rpar = (double *) R_alloc(np, sizeof(double));
    for (i = 0; i < np; i++)
      rpar[i] = REAL(Rpar)[i];

    is.f = NULL;
    is.env = NULL;
    Rdqags(Rintfn, (void*)&is,
	   &lower, &upper, &epsabs, &epsrel, &result,
	   &abserr, &neval, &ier, &limit, &lenw, &last, iwork, work);

    PROTECT(ans = allocVector(VECSXP, 4));
    PROTECT(ansnames = allocVector(STRSXP, 4));
    SET_STRING_ELT(ansnames, 0, mkChar("value"));
    SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, 1));
    REAL(VECTOR_ELT(ans, 0))[0] = result;
    SET_STRING_ELT(ansnames, 1, mkChar("abs.error"));
    SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, 1));
    REAL(VECTOR_ELT(ans, 1))[0] = abserr;
    SET_STRING_ELT(ansnames, 2, mkChar("subdivisions"));
    SET_VECTOR_ELT(ans, 2, allocVector(INTSXP, 1));
    INTEGER(VECTOR_ELT(ans, 2))[0] = last;
    SET_STRING_ELT(ansnames, 3, mkChar("ierr"));
    SET_VECTOR_ELT(ans, 3, allocVector(INTSXP, 1));
    INTEGER(VECTOR_ELT(ans, 3))[0] = ier;
    setAttrib(ans, R_NamesSymbol, ansnames);
    UNPROTECT(2);
    return ans;
}

/* 	cc_call_dqagi(ff, as.double(bound), as.integer(inf),
			as.double(abs.tol), as.double(rel.tol), imit = limit, 
        as.double(rpar), as.integer(ipar))
*/
SEXP cc_call_dqagi(SEXP ff, SEXP Bound, SEXP Inf, SEXP Atol, SEXP Rtol,  SEXP Limit, SEXP Rpar, SEXP Ipar)
{
    int_struct is;
    SEXP ans, ansnames;
    double bound, epsabs, epsrel, result, abserr, *work;
    int i, np, inf, neval, ier, limit, lenw, last, *iwork;

    if (!inherits(ff, "NativeSymbol")) 
       error("'f' is not a compiled function");
    fint = (C_func_type *) R_ExternalPtrAddr(ff);  

    bound = REAL(Bound)[0];
    inf = INTEGER(Inf)[0];
    epsabs = REAL(Atol)[0];
    epsrel = REAL(Rtol)[0];
    limit = INTEGER(Limit)[0];    
    lenw = 4 * limit;
    iwork = (int *) R_alloc((size_t) limit, sizeof(int));
    work = (double *) R_alloc((size_t) lenw, sizeof(double));

    np = LENGTH(Ipar);
    ipar = (int *) R_alloc(np, sizeof(int));
    for (i = 0; i < np; i++)
      ipar[i] = INTEGER(Ipar)[i];

    np = LENGTH(Rpar);
    rpar = (double *) R_alloc(np, sizeof(double));
    for (i = 0; i < np; i++)
      rpar[i] = REAL(Rpar)[i];

    is.f = NULL;
    is.env = NULL;

    Rdqagi(Rintfn, (void*)&is, &bound,&inf,&epsabs,&epsrel,&result,
	   &abserr,&neval,&ier,&limit,&lenw,&last,iwork,work);

    PROTECT(ans = allocVector(VECSXP, 4));
    PROTECT(ansnames = allocVector(STRSXP, 4));
    SET_STRING_ELT(ansnames, 0, mkChar("value"));
    SET_VECTOR_ELT(ans, 0, allocVector(REALSXP, 1));
    REAL(VECTOR_ELT(ans, 0))[0] = result;
    SET_STRING_ELT(ansnames, 1, mkChar("abs.error"));
    SET_VECTOR_ELT(ans, 1, allocVector(REALSXP, 1));
    REAL(VECTOR_ELT(ans, 1))[0] = abserr;
    SET_STRING_ELT(ansnames, 2, mkChar("subdivisions"));
    SET_VECTOR_ELT(ans, 2, allocVector(INTSXP, 1));
    INTEGER(VECTOR_ELT(ans, 2))[0] = last;
    SET_STRING_ELT(ansnames, 3, mkChar("ierr"));
    SET_VECTOR_ELT(ans, 3, allocVector(INTSXP, 1));
    INTEGER(VECTOR_ELT(ans, 3))[0] = ier;
    setAttrib(ans, R_NamesSymbol, ansnames);
    UNPROTECT(2);
    return ans;
}
