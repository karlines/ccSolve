/*
 This is an adaptation of parts of the file "optimize.c" of the R-core 
 software.
 It has been adapted to work with compiled code solutions 
 bu Karline Soetaert.

 .. removed underscore in error(_(
 .. changed .External2 calls to .Call
*/

#include <R.h>
#include <Rdefines.h>
#include "ccSolve.h"

/* One Dimensional Minimization --- just wrapper for
 * Brent's "fmin" --> ../appl/fmin.c */

struct callinfo {
  SEXP R_fcall;
  SEXP R_env;
} ;

/*static SEXP R_fcall1;
  static SEXP R_env1; */

int maximize;


/* compiled code version of fcn1 for use with optimize NOTE: made simpler */
static double cc_fcn1(double x, struct callinfo *info)
{
    int i;
    double val;
    int n = 1;

	  if (!R_FINITE(x)) error("non-finite value supplied by optimize");

    fcall(&n, &x, &val, rpar, ipar);    
	  if (!R_FINITE(val)) {
      error("non-finite value supplied by optimize");
      return 0;
    }
    if (maximize) val = - val;
    return val;
}

/* fmin(f, xmin, xmax tol, maximum) */
//cc_do_fmin", f, as.double(lower), as.double(upper), as.double(tol), as.integer(maximum))
SEXP cc_do_fmin(SEXP f, SEXP Lower, SEXP Upper, SEXP Tol, SEXP maximum, SEXP Rpar, SEXP Ipar)
{
    double xmin, xmax, tol;
    SEXP res;
    double val, x;
    int n = 1;
    
    struct callinfo info;

    /* the function to be minimized */
    if (!inherits(f, "NativeSymbol")) 
       error("'f' is not a compiled function");
    fcall = (C_func_type *) R_ExternalPtrAddr(f);  

    np = LENGTH(Ipar);
    ipar = (int *) R_alloc(np, sizeof(int));
    for (i = 0; i < np; i++)
      ipar[i] = INTEGER(Ipar)[i];

    np = LENGTH(Rpar);
    rpar = (double *) R_alloc(np, sizeof(double));
    for (i = 0; i < np; i++)
      rpar[i] = REAL(Rpar)[i];

    /* xmin, xmax */
    xmin = REAL(Lower);
    if (!R_FINITE(xmin))
  	  error(_("invalid '%s' value"), "xmin");
    xmax = REAL(Upper);
    if (!R_FINITE(xmax))
	     error(_("invalid '%s' value"), "xmax");
    if (xmin >= xmax)
	    error(_("'xmin' not less than 'xmax'"));

    /* tol */

    tol = REAL(Tol);
    if (!R_FINITE(tol) || tol <= 0.0)
	     error(_("invalid '%s' value"), "tol");

    PROTECT(res = allocVector(REALSXP, 2));  /* karline: made 2 long; second = f-value*/

    x = Brent_fmin(xmin, xmax,
			      (double (*)(double, void*)) cc_fcn1, &info, tol);
    REAL(res)[0] = x;   
    fcall(&n, &x, &val, rpar, ipar);
    REAL(res)[1] = val;   
    UNPROTECT(1);    
			      
    return res;
}



// One Dimensional Root Finding --  just wrapper code for
// Brent's "zeroin"
// ---------------

extern double 
R_zeroin2(double ax, double bx, double fa, double fb, 
	  double (*f)(double x, void *info), void *info, 
	  double *Tol, int *Maxit);


static double fcn2(double x, struct callinfo *info)
{
    SEXP s, sx;
    PROTECT(sx = ScalarReal(x));
    SETCADR(info->R_fcall, sx);
    s = eval(info->R_fcall, info->R_env);
    UNPROTECT(1);
    switch(TYPEOF(s)) {
    case INTSXP:
	if (length(s) != 1) goto badvalue;
	if (INTEGER(s)[0] == NA_INTEGER) {
	    warning(_("NA replaced by maximum positive value"));
	    return	DBL_MAX;
	}
	else return INTEGER(s)[0];
	break;
    case REALSXP:
	if (length(s) != 1) goto badvalue;
	if (!R_FINITE(REAL(s)[0])) {
	    if(REAL(s)[0] == R_NegInf) { // keep sign for root finding !
		warning(_("-Inf replaced by maximally negative value"));
		return -DBL_MAX;
	    } else {
		warning(_("NA/Inf replaced by maximum positive value"));
		return DBL_MAX;
	    }
	}
	else return REAL(s)[0];
	break;
    default:
	goto badvalue;
    }
 badvalue:
    error(_("invalid function value in 'zeroin'"));
    return 0;/* for -Wall */

}

/* zeroin2(f, ax, bx, f.ax, f.bx, tol, maxiter) */
SEXP zeroin2(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    double f_ax, f_bx;
    double xmin, xmax, tol;
    int iter;
    SEXP v, res;
    struct callinfo info;

    args = CDR(args);
    PrintDefaults();

    /* the function to be minimized */
    v = CAR(args);
    if (!isFunction(v)) error(_("attempt to minimize non-function"));
    args = CDR(args);

    /* xmin */
    xmin = asReal(CAR(args));
    if (!R_FINITE(xmin)) error(_("invalid '%s' value"), "xmin");
    args = CDR(args);

    /* xmax */
    xmax = asReal(CAR(args));
    if (!R_FINITE(xmax)) error(_("invalid '%s' value"), "xmax");
    if (xmin >= xmax) error(_("'xmin' not less than 'xmax'"));
    args = CDR(args);

    /* f(ax) = f(xmin) */
    f_ax = asReal(CAR(args));
    if (ISNA(f_ax)) error(_("NA value for '%s' is not allowed"), "f.lower");
    args = CDR(args);

    /* f(bx) = f(xmax) */
    f_bx = asReal(CAR(args));
    if (ISNA(f_bx)) error(_("NA value for '%s' is not allowed"), "f.upper");
    args = CDR(args);

    /* tol */
    tol = asReal(CAR(args));
    if (!R_FINITE(tol) || tol <= 0.0) error(_("invalid '%s' value"), "tol");
    args = CDR(args);

    /* maxiter */
    iter = asInteger(CAR(args));
    if (iter <= 0) error(_("'maxiter' must be positive"));

    info.R_env = rho;
    PROTECT(info.R_fcall = lang2(v, R_NilValue)); /* the info used in fcn2() */
    PROTECT(res = allocVector(REALSXP, 3));
    REAL(res)[0] =
	R_zeroin2(xmin, xmax, f_ax, f_bx, (double (*)(double, void*)) fcn2,
		 (void *) &info, &tol, &iter);
    REAL(res)[1] = (double)iter;
    REAL(res)[2] = tol;
    UNPROTECT(2);
    return res;
}



/* General Nonlinear Optimization */

#define FT_SIZE 5		/* default size of table to store computed
				   function values */

typedef struct {
  double   fval;
  double  *x;
  double  *grad;
  double  *hess;
} ftable;

typedef struct {
  SEXP R_fcall;	      /* unevaluated call to R function */
  SEXP R_env;	      /* where to evaluate the calls */
  int have_gradient;
  int have_hessian;
/*  int n;	      -* length of the parameter (x) vector */
  int FT_size;	      /* size of table to store computed
			 function values */
  int FT_last;	      /* Newest entry in the table */
  ftable *Ftable;
} function_info;

/* Initialize the storage in the table of computed function values */

static void FT_init(int n, int FT_size, function_info *state)
{
    int i, j;
    int have_gradient, have_hessian;
    ftable *Ftable;

    have_gradient = state->have_gradient;
    have_hessian = state->have_hessian;

    Ftable = (ftable *)R_alloc(FT_size, sizeof(ftable));

    for (i = 0; i < FT_size; i++) {
	Ftable[i].x = (double *)R_alloc(n, sizeof(double));
				/* initialize to unlikely parameter values */
	for (j = 0; j < n; j++) {
	    Ftable[i].x[j] = DBL_MAX;
	}
	if (have_gradient) {
	    Ftable[i].grad = (double *)R_alloc(n, sizeof(double));
	    if (have_hessian) {
		Ftable[i].hess = (double *)R_alloc(n * n, sizeof(double));
	    }
	}
    }
    state->Ftable = Ftable;
    state->FT_size = FT_size;
    state->FT_last = -1;
}

/* Store an entry in the table of computed function values */

static void FT_store(int n, const double f, const double *x, const double *grad,
		     const double *hess, function_info *state)
{
    int ind;

    ind = (++(state->FT_last)) % (state->FT_size);
    state->Ftable[ind].fval = f;
    Memcpy(state->Ftable[ind].x, x, n);
    if (grad) {
	Memcpy(state->Ftable[ind].grad, grad, n);
	if (hess) {
	    Memcpy(state->Ftable[ind].hess, hess, n * n);
	}
    }
}

/* Check for stored values in the table of computed function values.
   Returns the index in the table or -1 for failure */

static int FT_lookup(int n, const double *x, function_info *state)
{
    double *ftx;
    int i, j, ind, matched;
    int FT_size, FT_last;
    ftable *Ftable;

    FT_last = state->FT_last;
    FT_size = state->FT_size;
    Ftable = state->Ftable;

    for (i = 0; i < FT_size; i++) {
	ind = (FT_last - i) % FT_size;
				/* why can't they define modulus correctly */
	if (ind < 0) ind += FT_size;
	ftx = Ftable[ind].x;
	if (ftx) {
	    matched = 1;
	    for (j = 0; j < n; j++) {
		if (x[j] != ftx[j]) {
		    matched = 0;
		    break;
		}
	    }
	    if (matched) return ind;
	}
    }
    return -1;
}

/* This how the optimizer sees them */

static void fcn(int n, const double x[], double *f, function_info
		*state)
{
    SEXP s, R_fcall;
    ftable *Ftable;
    double *g = (double *) 0, *h = (double *) 0;
    int i;

    R_fcall = state->R_fcall;
    Ftable = state->Ftable;
    if ((i = FT_lookup(n, x, state)) >= 0) {
	*f = Ftable[i].fval;
	return;
    }
				/* calculate for a new value of x */
    s = CADR(R_fcall);
    for (i = 0; i < n; i++) {
	if (!R_FINITE(x[i])) error(_("non-finite value supplied by 'nlm'"));
	REAL(s)[i] = x[i];
    }
    s = PROTECT(eval(state->R_fcall, state->R_env));
    switch(TYPEOF(s)) {
    case INTSXP:
	if (length(s) != 1) goto badvalue;
	if (INTEGER(s)[0] == NA_INTEGER) {
	    warning(_("NA replaced by maximum positive value"));
	    *f = DBL_MAX;
	}
	else *f = INTEGER(s)[0];
	break;
    case REALSXP:
	if (length(s) != 1) goto badvalue;
	if (!R_FINITE(REAL(s)[0])) {
	    warning(_("NA/Inf replaced by maximum positive value"));
	    *f = DBL_MAX;
	}
	else *f = REAL(s)[0];
	break;
    default:
	goto badvalue;
    }
    if (state->have_gradient) {
	g = REAL(PROTECT(coerceVector(getAttrib(s, install("gradient")), REALSXP)));
	if (state->have_hessian) {
	    h = REAL(PROTECT(coerceVector(getAttrib(s, install("hessian")), REALSXP)));
	}
    }
    FT_store(n, *f, x, g, h, state);
    UNPROTECT(1 + state->have_gradient + state->have_hessian);
    return;

 badvalue:
    error(_("invalid function value in 'nlm' optimizer"));
}


static void Cd1fcn(int n, const double x[], double *g, function_info *state)
{
    int ind;

    if ((ind = FT_lookup(n, x, state)) < 0) {	/* shouldn't happen */
	fcn(n, x, g, state);
	if ((ind = FT_lookup(n, x, state)) < 0) {
	    error(_("function value caching for optimization is seriously confused"));
	}
    }
    Memcpy(g, state->Ftable[ind].grad, n);
}


static void Cd2fcn(int nr, int n, const double x[], double *h,
		   function_info *state)
{
    int j, ind;

    if ((ind = FT_lookup(n, x, state)) < 0) {	/* shouldn't happen */
	fcn(n, x, h, state);
	if ((ind = FT_lookup(n, x, state)) < 0) {
	    error(_("function value caching for optimization is seriously confused"));
	}
    }
    for (j = 0; j < n; j++) {  /* fill in lower triangle only */
	Memcpy( h + j*(n + 1), state->Ftable[ind].hess + j*(n + 1), n - j);
    }
}


static double *fixparam(SEXP p, int *n)
{
    double *x;
    int i;

    if (!isNumeric(p))
	error(_("numeric parameter expected"));

    if (*n) {
	if (LENGTH(p) != *n)
	    error(_("conflicting parameter lengths"));
    }
    else {
	if (LENGTH(p) <= 0)
	    error(_("invalid parameter length"));
	*n = LENGTH(p);
    }

    x = (double*)R_alloc(*n, sizeof(double));
    switch(TYPEOF(p)) {
    case LGLSXP:
    case INTSXP:
	for (i = 0; i < *n; i++) {
	    if (INTEGER(p)[i] == NA_INTEGER)
		error(_("missing value in parameter"));
	    x[i] = INTEGER(p)[i];
	}
	break;
    case REALSXP:
	for (i = 0; i < *n; i++) {
	    if (!R_FINITE(REAL(p)[i]))
		error(_("missing value in parameter"));
	    x[i] = REAL(p)[i];
	}
	break;
    default:
	error(_("invalid parameter type"));
    }
    return x;
}

	/* Fatal errors - we don't deliver an answer */

static void opterror(int nerr)
{
    switch(nerr) {
    case -1:
	error(_("non-positive number of parameters in nlm"));
    case -2:
	error(_("nlm is inefficient for 1-d problems"));
    case -3:
	error(_("invalid gradient tolerance in nlm"));
    case -4:
	error(_("invalid iteration limit in nlm"));
    case -5:
	error(_("minimization function has no good digits in nlm"));
    case -6:
	error(_("no analytic gradient to check in nlm!"));
    case -7:
	error(_("no analytic Hessian to check in nlm!"));
    case -21:
	error(_("probable coding error in analytic gradient"));
    case -22:
	error(_("probable coding error in analytic Hessian"));
    default:
	error(_("*** unknown error message (msg = %d) in nlm()\n*** should not happen!"), nerr);
    }
}


	/* Warnings - we return a value, but print a warning */

static void optcode(int code)
{
    switch(code) {
    case 1:
	Rprintf(_("Relative gradient close to zero.\n"));
	Rprintf(_("Current iterate is probably solution.\n"));
	break;
    case 2:
	Rprintf(_("Successive iterates within tolerance.\n"));
	Rprintf(_("Current iterate is probably solution.\n"));
	break;
    case 3:
	Rprintf(_("Last global step failed to locate a point lower than x.\n"));
	Rprintf(_("Either x is an approximate local minimum of the function,\n\
the function is too non-linear for this algorithm,\n\
or steptol is too large.\n"));
	break;
    case 4:
	Rprintf(_("Iteration limit exceeded.  Algorithm failed.\n"));
	break;
    case 5:
	Rprintf(_("Maximum step size exceeded 5 consecutive times.\n\
Either the function is unbounded below,\n\
becomes asymptotic to a finite value\n\
from above in some direction,\n"\
"or stepmx is too small.\n"));
	break;
    }
    Rprintf("\n");
}

/* NOTE: The actual Dennis-Schnabel algorithm `optif9' is in ../appl/uncmin.c */

SEXP nlm(SEXP call, SEXP op, SEXP args, SEXP rho)
{
    SEXP value, names, v, R_gradientSymbol, R_hessianSymbol;

    double *x, *typsiz, fscale, gradtl, stepmx,
	steptol, *xpls, *gpls, fpls, *a, *wrk, dlt;

    int code, i, j, k, itnlim, method, iexp, omsg, msg,
	n, ndigit, iagflg, iahflg, want_hessian, itncnt;


/* .Internal(
 *	nlm(function(x) f(x, ...), p, hessian, typsize, fscale,
 *	    msg, ndigit, gradtol, stepmax, steptol, iterlim)
 */
    function_info *state;

    args = CDR(args);
    PrintDefaults();

    state = (function_info *) R_alloc(1, sizeof(function_info));

    /* the function to be minimized */

    v = CAR(args);
    if (!isFunction(v))
	error(_("attempt to minimize non-function"));
    PROTECT(state->R_fcall = lang2(v, R_NilValue));
    args = CDR(args);

    /* `p' : inital parameter value */

    n = 0;
    x = fixparam(CAR(args), &n);
    args = CDR(args);

    /* `hessian' : H. required? */

    want_hessian = asLogical(CAR(args));
    if (want_hessian == NA_LOGICAL) want_hessian = 0;
    args = CDR(args);

    /* `typsize' : typical size of parameter elements */

    typsiz = fixparam(CAR(args), &n);
    args = CDR(args);

    /* `fscale' : expected function size */

    fscale = asReal(CAR(args));
    if (ISNA(fscale)) error(_("invalid NA value in parameter"));
    args = CDR(args);

    /* `msg' (bit pattern) */
    omsg = msg = asInteger(CAR(args));
    if (msg == NA_INTEGER) error(_("invalid NA value in parameter"));
    args = CDR(args);

    ndigit = asInteger(CAR(args));
    if (ndigit == NA_INTEGER) error(_("invalid NA value in parameter"));
    args = CDR(args);

    gradtl = asReal(CAR(args));
    if (ISNA(gradtl)) error(_("invalid NA value in parameter"));
    args = CDR(args);

    stepmx = asReal(CAR(args));
    if (ISNA(stepmx)) error(_("invalid NA value in parameter"));
    args = CDR(args);

    steptol = asReal(CAR(args));
    if (ISNA(steptol)) error(_("invalid NA value in parameter"));
    args = CDR(args);

    /* `iterlim' (def. 100) */
    itnlim = asInteger(CAR(args));
    if (itnlim == NA_INTEGER) error(_("invalid NA value in parameter"));

    state->R_env = rho;

    /* force one evaluation to check for the gradient and hessian */
    iagflg = 0;			/* No analytic gradient */
    iahflg = 0;			/* No analytic hessian */
    state->have_gradient = 0;
    state->have_hessian = 0;
    R_gradientSymbol = install("gradient");
    R_hessianSymbol = install("hessian");

    /* This vector is shared with all subsequent calls */
    v = allocVector(REALSXP, n);
    for (i = 0; i < n; i++) REAL(v)[i] = x[i];
    SETCADR(state->R_fcall, v);
    SET_NAMED(v, 2); // in case the functions try to alter it
    value = eval(state->R_fcall, state->R_env);

    v = getAttrib(value, R_gradientSymbol);
    if (v != R_NilValue) {
	if (LENGTH(v) == n && (isReal(v) || isInteger(v))) {
	    iagflg = 1;
	    state->have_gradient = 1;
	    v = getAttrib(value, R_hessianSymbol);

	    if (v != R_NilValue) {
		if (LENGTH(v) == (n * n) && (isReal(v) || isInteger(v))) {
		    iahflg = 1;
		    state->have_hessian = 1;
		} else {
		    warning(_("hessian supplied is of the wrong length or mode, so ignored"));
		}
	    }
	} else {
	    warning(_("gradient supplied is of the wrong length or mode, so ignored"));
	}
    }
    if (((msg/4) % 2) && !iahflg) { /* skip check of analytic Hessian */
      msg -= 4;
    }
    if (((msg/2) % 2) && !iagflg) { /* skip check of analytic gradient */
      msg -= 2;
    }
    FT_init(n, FT_SIZE, state);
    /* Plug in the call to the optimizer here */

    method = 1;	/* Line Search */
    iexp = iahflg ? 0 : 1; /* Function calls are expensive */
    dlt = 1.0;

    xpls = (double*)R_alloc(n, sizeof(double));
    gpls = (double*)R_alloc(n, sizeof(double));
    a = (double*)R_alloc(n*n, sizeof(double));
    wrk = (double*)R_alloc(8*n, sizeof(double));

    /*
     *	 Dennis + Schnabel Minimizer
     *
     *	  SUBROUTINE OPTIF9(NR,N,X,FCN,D1FCN,D2FCN,TYPSIZ,FSCALE,
     *	 +	   METHOD,IEXP,MSG,NDIGIT,ITNLIM,IAGFLG,IAHFLG,IPR,
     *	 +	   DLT,GRADTL,STEPMX,STEPTOL,
     *	 +	   XPLS,FPLS,GPLS,ITRMCD,A,WRK)
     *
     *
     *	 Note: I have figured out what msg does.
     *	 It is actually a sum of bit flags as follows
     *	   1 = don't check/warn for 1-d problems
     *	   2 = don't check analytic gradients
     *	   4 = don't check analytic hessians
     *	   8 = don't print start and end info
     *	  16 = print at every iteration
     *	 Using msg=9 is absolutely minimal
     *	 I think we always check gradients and hessians
     */

    optif9(n, n, x, (fcn_p) fcn, (fcn_p) Cd1fcn, (d2fcn_p) Cd2fcn,
	   state, typsiz, fscale, method, iexp, &msg, ndigit, itnlim,
	   iagflg, iahflg, dlt, gradtl, stepmx, steptol, xpls, &fpls,
	   gpls, &code, a, wrk, &itncnt);

    if (msg < 0)
	opterror(msg);
    if (code != 0 && (omsg&8) == 0)
	optcode(code);

    if (want_hessian) {
	PROTECT(value = allocVector(VECSXP, 6));
	PROTECT(names = allocVector(STRSXP, 6));
	fdhess(n, xpls, fpls, (fcn_p) fcn, state, a, n, &wrk[0], &wrk[n],
	       ndigit, typsiz);
	for (i = 0; i < n; i++)
	    for (j = 0; j < i; j++)
		a[i + j * n] = a[j + i * n];
    }
    else {
	PROTECT(value = allocVector(VECSXP, 5));
	PROTECT(names = allocVector(STRSXP, 5));
    }
    k = 0;

    SET_STRING_ELT(names, k, mkChar("minimum"));
    SET_VECTOR_ELT(value, k, ScalarReal(fpls));
    k++;

    SET_STRING_ELT(names, k, mkChar("estimate"));
    SET_VECTOR_ELT(value, k, allocVector(REALSXP, n));
    for (i = 0; i < n; i++)
	REAL(VECTOR_ELT(value, k))[i] = xpls[i];
    k++;

    SET_STRING_ELT(names, k, mkChar("gradient"));
    SET_VECTOR_ELT(value, k, allocVector(REALSXP, n));
    for (i = 0; i < n; i++)
	REAL(VECTOR_ELT(value, k))[i] = gpls[i];
    k++;

    if (want_hessian) {
	SET_STRING_ELT(names, k, mkChar("hessian"));
	SET_VECTOR_ELT(value, k, allocMatrix(REALSXP, n, n));
	for (i = 0; i < n * n; i++)
	    REAL(VECTOR_ELT(value, k))[i] = a[i];
	k++;
    }

    SET_STRING_ELT(names, k, mkChar("code"));
    SET_VECTOR_ELT(value, k, allocVector(INTSXP, 1));
    INTEGER(VECTOR_ELT(value, k))[0] = code;
    k++;

    /* added by Jim K Lindsey */
    SET_STRING_ELT(names, k, mkChar("iterations"));
    SET_VECTOR_ELT(value, k, allocVector(INTSXP, 1));
    INTEGER(VECTOR_ELT(value, k))[0] = itncnt;
    k++;

    setAttrib(value, R_NamesSymbol, names);
    UNPROTECT(3);
    return value;
}
