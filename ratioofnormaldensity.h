#ifndef RATIOOFNORMALSDENSITY_H_
#define RATIOOFNORMALSDENSITY_H_

/**
// This function calculates the density of a variable which is a ratio
// of two variables each of which is distributed according to a Gaussian
// (normal) distribution.
//
// Formulae used in the implementation of this function come from
// D.V. Hinkley, "On the ratio of two correlated normal random variables",
// Biometrika, v. 56, p. 635 (1969).
//
// Function arguments are as follows:
//
//   x   -- Coordinate at which the density value is desired.
//
//   mu1 -- Mean of the Gaussian in the numerator.
//
//   s1  -- Standard deviation of the Gaussian in the numerator.
//          Must be positive.
//
//   mu2 -- Mean of the Gaussian in the denominator.
//
//   s2  -- Standard deviation of the Gaussian in the denominator.
//          Must be positive.
//
//   rho -- The correlation coefficient between the two Gaussians.
//          Its magnitude must be less than 1.
//
// Igor Volobouev, November 2013
*/
double ratioOfNormalsDensity(double x, double mu1, double s1,
                             double mu2, double s2, double rho);

#endif // RATIOOFNORMALSDENSITY_H_
