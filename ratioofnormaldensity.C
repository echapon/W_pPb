#include <assert.h>
#include <math.h>
#include <float.h>

#include "ratioOfNormalsDensity.h"

#define SQR2PI 2.5066282746310005
#define SQRT2  1.41421356237309505

double ratioOfNormalsDensity(const double w, const double mu1,
                             const double s1, const double mu2,
                             const double s2, const double rho)
{
    assert(s1 > 0.0);
    assert(s2 > 0.0);

    const double oneminusrhosq = 1.0 - rho*rho;
    assert(oneminusrhosq > 0.0);
    const double sqrtoneminusrhosq = sqrt(oneminusrhosq);

    double awsq = w*w/s1/s1 - 2.0*rho*w/s1/s2 + 1.0/s2/s2;
    if (awsq <= 0.0)
        awsq = DBL_MIN;
    const double aw = sqrt(awsq);

    const double bw = mu1*w/s1/s1 - 
        rho*(mu1 + mu2*w)/s1/s2 + mu2/s2/s2;

    const double c = mu1*mu1/s1/s1 - 2.0*rho*mu1*mu2/s1/s2 +
                     mu2*mu2/s2/s2;

    const double dw = exp((bw*bw - c*aw*aw)/2.0/oneminusrhosq/aw/aw);

    const double arg = bw/sqrtoneminusrhosq/aw;
    const double phiterm = 1.0 - 2.0*erfc(arg/SQRT2)/2.0;

    return bw*dw/SQR2PI/s1/s2/aw/aw/aw*phiterm + 
        sqrtoneminusrhosq/M_PI/s1/s2/awsq*exp(-c/2.0/oneminusrhosq);
}
