
/* type definitions for ccSolve, some come from R_ext/Applic.h               */
#include <R_ext/Applic.h>

typedef double optimfn(int, double *, void *);
typedef void optimgr(int, double *, double *, void *);

typedef void C_func_type(int*, double*, double*, double  *, int*);
typedef void C_func_type2(double *, double*, double *, int*);
typedef void C_jac_type(int*, double*, double*, double *, int*);

/*  from Applic.h
void vmmin(int n, double *b, double *Fmin,
	   optimfn fn, optimgr gr, int maxit, int trace,
	   int *mask, double abstol, double reltol, int nREPORT,
	   void *ex, int *fncount, int *grcount, int *fail);
void nmmin(int n, double *Bvec, double *X, double *Fmin, optimfn fn,
	   int *fail, double abstol, double intol, void *ex,
	   double alpha, double bet, double gamm, int trace,
	   int *fncount, int maxit);
void cgmin(int n, double *Bvec, double *X, double *Fmin,
	   optimfn fn, optimgr gr,
	   int *fail, double abstol, double intol, void *ex,
	   int type, int trace, int *fncount, int *grcount, int maxit);
void lbfgsb(int n, int m, double *x, double *l, double *u, int *nbd,
	    double *Fmin, optimfn fn, optimgr gr, int *fail, void *ex,
	    double factr, double pgtol, int *fncount, int *grcount,
	    int maxit, char *msg, int trace, int nREPORT);
void samin(int n, double *pb, double *yb, optimfn fn, int maxit,
	   int tmax, double ti, int trace, void *ex);
*/

/* use with cc_nls */
typedef void C_init_dat_type(void  (int *, double *), int*);
void Initstdat (int *N, double *pdata);

extern double R_zeroin2(double ax, double bx, double fa, double fb, 
	  double (*f)(double x, void *info), void *info, 
	  double *Tol, int *Maxit);

void cc_fcn_lmder(int *m, int *n, double *par, double *fvec, 
  double *fjac, int *ldfjac, int *iflag);
void F77_NAME(lmder)(void (*fcn_lmder)(int *m, int *n, double *par, double *fvec, double *fjac, int *ldfjac, int *iflag),
                     int *m, int *n, double *par, double *fvec,
                     double *fjac, int *ldfjac,
                     double *ftol, double *ptol, double *gtol,
                     int *maxfev, double *diag, int *mode,
                     double *factor, int *nprint, int *info,
                     int *nfev, int *njev, int *ipvt,
                     double *qtf, double *wa1, double *wa2,
                     double *wa3, double *wa4);
char *fcn_message(char*, int, int, int);
void cc_fcn_lmdif(int *m, int *n, double *par, double *fvec, int *iflag);
void F77_NAME(lmdif)(void (*fcn_lmdif)(int *m, int *n, double *par, double *fvec, int *iflag),
                     int *m, int *n, double *par, double *fvec,
                     double *ftol, double *ptol, double *gtol, int *maxfev,
                     double *epsfcn, double *diag, int *mode, double *factor,
                     int *nprint, int *info, int *nfev, double *fjac,
                     int *ldfjac, int *ipvt, double *qtf,
                     double *wa1, double *wa2, double *wa3, double *wa4);

/* from minpack.lm */
void transpose(double *x, int nrx, int ncx, double *y);
void matprod  (double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z);
void crossprod(double *x, int nrx, int ncx, double *y, int nry, int ncy, double *z);

SEXP getListElement(SEXP list, char *str);

double *real_vector(int n);
int  *int_vector(int n);


