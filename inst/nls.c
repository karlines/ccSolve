/*
 *  Routines used in calculating least squares solutions in a
 *  nonlinear model using compiled code
 *  Based on the file nls.c in R base, R-version 3.1.1 
 *  implementation: karline soetaert
 */

#include <R.h>
#include <Rdefines.h>
#include "ccSolve.h"

#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif

/*
 * get the list element named str. names is the name attribute of list
 */

static SEXP
getListElement(SEXP list, SEXP names, const char *str)
{
    SEXP elmt = (SEXP) NULL;
    const char *tempChar;
    int i;

    for (i = 0; i < length(list); i++) {
	tempChar = CHAR(STRING_ELT(names, i)); /* ASCII only */
	if( strcmp(tempChar,str) == 0) {
	    elmt = VECTOR_ELT(list, i);
	    break;
	}
    }
    return elmt;
}

/*
 * put some convergence-related information into list
 */
static SEXP
ConvInfoMsg(char* msg, int iter, int whystop, double fac,
	    double minFac, int maxIter, double convNew)
{
    const char *nms[] = {"isConv", "finIter", "finTol",
			 "stopCode", "stopMessage",  ""};
    SEXP ans;
    PROTECT(ans = mkNamed(VECSXP, nms));

    SET_VECTOR_ELT(ans, 0, ScalarLogical(whystop == 0)); /* isConv */
    SET_VECTOR_ELT(ans, 1, ScalarInteger(iter));	 /* finIter */
    SET_VECTOR_ELT(ans, 2, ScalarReal   (convNew));	 /* finTol */
    SET_VECTOR_ELT(ans, 3, ScalarInteger(whystop));      /* stopCode */
    SET_VECTOR_ELT(ans, 4, mkString(msg));               /* stopMessage */

    UNPROTECT(1);
    return ans;
}


/*
 *  call to nls_iter from R --- .Call("nls_iter", m, control, doTrace)
 *  where m and control are nlsModel and nlsControl objects
 *             doTrace is a logical value.
 *  m is modified; the return value is a "convergence-information" list.
 */
SEXP
cc_nls_iter(SEXP fn, SEXP initfunc, SEXP Pars, SEXP Nvar, SEXP Ndata, SEXP data, 
  SEXP control, SEXP doTraceArg)
{
    double dev, fac, minFac, tolerance, newDev, convNew = -1./*-Wall*/;
    int i, j, maxIter, hasConverged, nPars, doTrace, evaltotCnt = -1, warnOnly, printEval;
    SEXP tmp, conv, incr, deviance, setPars, getPars, pars, newPars, trace;

    C_init_dat_type *initializer;
    C_nls_type *fnls = NULL;
    
    if (!inherits(fn, "NativeSymbol")) 
       error("'fn' is not a compiled function");
    fnls = (C_nls_type *) R_ExternalPtrAddr(fn);  

    if (!inherits(initfunc, "NativeSymbol")) 
       error("'initfunc' is not a compiled function");

    doTrace = asLogical(doTraceArg);

    if(!isNewList(control))
	    error("'control' must be a list");

    PROTECT(tmp = getAttrib(control, R_NamesSymbol));

    conv = getListElement(control, tmp, "maxiter");
    if(conv == NULL || !isNumeric(conv))
	error("'%s' absent", "control$maxiter");
    maxIter = asInteger(conv);

    conv = getListElement(control, tmp, "tol");
    if(conv == NULL || !isNumeric(conv))
	error("'%s' absent", "control$tol");
    tolerance = asReal(conv);

    conv = getListElement(control, tmp, "minFactor");
    if(conv == NULL || !isNumeric(conv))
	error("'%s' absent", "control$minFactor");
    minFac = asReal(conv);

    conv = getListElement(control, tmp, "warnOnly");
    if(conv == NULL || !isLogical(conv))
	error("'%s' absent", "control$warnOnly");
    warnOnly = asLogical(conv);

    conv = getListElement(control, tmp, "printEval");
    if(conv == NULL || !isLogical(conv))
	error("'%s' absent", "control$printEval");
    printEval = asLogical(conv);

#define CONV_INFO_MSG(_STR_, _I_)					\
	ConvInfoMsg(_STR_, i, _I_, fac, minFac, maxIter, convNew)

#define NON_CONV_FINIS(_ID_, _MSG_)		\
    if(warnOnly) {				\
	warning(_MSG_);				\
	return CONV_INFO_MSG(_MSG_, _ID_);      \
    }						\
    else					\
	error(_MSG_);

#define NON_CONV_FINIS_1(_ID_, _MSG_, _A1_)	\
    if(warnOnly) {				\
	char msgbuf[1000];			\
	warning(_MSG_, _A1_);			\
	snprintf(msgbuf, 1000, _MSG_, _A1_);	\
	return CONV_INFO_MSG(msgbuf, _ID_);	\
    }						\
    else					\
	error(_MSG_, _A1_);

#define NON_CONV_FINIS_2(_ID_, _MSG_, _A1_, _A2_)	\
    if(warnOnly) {					\
	char msgbuf[1000];				\
	warning(_MSG_, _A1_, _A2_);			\
	snprintf(msgbuf, 1000, _MSG_, _A1_, _A2_);	\
	return CONV_INFO_MSG(msgbuf, _ID_);		\
    }							\
    else						\
	error(_MSG_, _A1_, _A2_);

                    
    fac = 1.0;
    hasConverged = FALSE;

    PROTECT(newPars = allocVector(REALSXP, nPars));
    if(printEval)
	evaltotCnt = 1;
    for (i = 0; i < maxIter; i++) {
	SEXP newIncr;
	int evalCnt = -1;
	if((convNew = asReal(eval(conv, R_GlobalEnv))) < tolerance) {
	    hasConverged = TRUE;
	    break;
	}
	PROTECT(newIncr = eval(incr, R_GlobalEnv));

	if(printEval)
	    evalCnt = 1;

	while(fac >= minFac) {
	    if(printEval) {
		Rprintf("  It. %3d, fac= %11.6g, eval (no.,total): (%2d,%3d):",
			i+1, fac, evalCnt, evaltotCnt);
		evalCnt++;
		evaltotCnt++;
	    }
	    for(j = 0; j < nPars; j++)
		REAL(newPars)[j] = REAL(pars)[j] + fac * REAL(newIncr)[j];

	    PROTECT(tmp = lang2(setPars, newPars));
	    if (asLogical(eval(tmp, R_GlobalEnv))) { /* singular gradient */
		UNPROTECT(11);

		NON_CONV_FINIS(1, _("singular gradient"));
	    }
	    UNPROTECT(1);

	    newDev = asReal(eval(deviance, R_GlobalEnv));
	    if(printEval)
		Rprintf(" new dev = %g\n", newDev);
	    if(newDev <= dev) {
		dev = newDev;
		fac = MIN(2*fac, 1);
		tmp = newPars;
		newPars = pars;
		pars = tmp;
		break;
	    }
	    fac /= 2.;
	}
	UNPROTECT(1);
	if( fac < minFac ) {
	    UNPROTECT(9);
	    NON_CONV_FINIS_2(2,
			     _("step factor %g reduced below 'minFactor' of %g"),
			     fac, minFac);
	}
	if(doTrace) eval(trace, R_GlobalEnv);
    }

    UNPROTECT(9);
    if(!hasConverged) {
	NON_CONV_FINIS_1(3,
			 _("number of iterations exceeded maximum of %d"),
			 maxIter);
    }
    /* else */
initializer = (C_init_dat_type *) R_ExternalPtrAddr(Initfunc);
initializer(Initnlsdata);

    return CONV_INFO_MSG(_("converged"), 0);
}
#undef CONV_INFO_MSG
#undef NON_CONV_FINIS
#undef NON_CONV_FINIS_1
#undef NON_CONV_FINIS_2




void Initnlsdata(int *N, double *data)
{
  int i, Nforcs;

  nvars = INTEGER(Nvars)[0];   
  ndata = INTEGER(Ndata)[0];   
  if ((*N) > nvars - 1)
    {
      warning("Number of variables passed to solver, %i; number in DLL, %i\n",
      Nforcs, *N);
      PROBLEM "Confusion over the number of variables"
      ERROR;
    } 
  else {
      for (i = 0; i < ndata; i++) data[i] = REAL(Data)[i + (*N - 1)*nvars];
  }
}

