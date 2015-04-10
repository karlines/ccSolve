#include <R.h>

char *fcn_message(char *msg, int info, int n, int nit)
{
    if      (info == 1)
        sprintf(msg, "Relative error in the sum of squares is at most `ftol'.");
    else if (info == 2)
        sprintf(msg, "Relative error between `par' and the solution is at most `ptol'.");
    else if (info == 3)
        sprintf(msg, "Conditions for `info = 1' and `info = 2' both hold.");
    else if (info == 4)
        sprintf(msg, "The cosine of the angle between `fvec' and any column of the Jacobian is at most `gtol' in absolute value.");
    else if (info == 5)
        sprintf(msg, "Number of calls to `fcn' has reached or exceeded `maxfev' == %d.", n);
    else if (info == 6)
        sprintf(msg, "`ftol' is too small. No further reduction in the sum of squares is possible.");
    else if (info == 7)
        sprintf(msg, "`ptol' is too small. No further improvement in the approximate solution `par' is possible.");
    else if (info == 8)
        sprintf(msg, "`gtol' is too small. `fvec' is orthogonal to the columns of the Jacobian to machine precision.");
    else if (info < 0)
      sprintf(msg, "Number of iterations has reached `maxiter' == %d.", nit);
    else if (info == 0)
        sprintf(msg, "Improper input parameters.");
    return msg;
}
