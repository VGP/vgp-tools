/*****************************************************************************************\
*                                                                                         *
*  Distribution generator data abstraction                                                *
*     One can create a distribution generator for a number of parameterized distribution  *
*     types and then generate events with that distribution                               *
*                                                                                         *
*  Author:  Gene Myers                                                                    *
*  Date  :  January 2007                                                                  *
*                                                                                         *
\*****************************************************************************************/

#ifndef _SR_CDF
#define _SR_CDF

#include "gene_core.h"
#include "cdf.h"

//  My own implementation of erand/eseed:
//     The current state of a random sequence generator is held in a uint64.
//     Myseed generates an initial state from a given seed value (such as the PID of the process)
//     Myrand advances the state of the generator and returns a pseudo-random value in [0,1)

double myrand(uint64 *state);
uint64 myseed(uint32 seedval);

//  A CDF is a random number generator for a particular distribution, unlike myrand which
//    gnenerates numbers from the uniform distribution over [0,1].  One can create a generator
//    from the currently implemented list below:
//
//    Distribution     PDF                             Mean      Variance     Domain
//
//    Normal(m,s)      ...                             m         s^2          R (all reals)
//    Exponential(a)   ae^(-ax)                        1/a       1/a          R >= 0
//    Poisson(a)       a^k e^-a / k!                   a         a            N (non-neg integers) 
//    Geometric(p)     (1-p)^k-1 p                     1/p       (1-p)/p^2    N
//    Binomial(n,p)    (n choose k) p^k (1-p)^(n-k)    np        np(1-p)      N in [0,n]
//    Uniform(l,h)     1/(h-l)                         (l+h)/2   (h-l)^2/12   R in [l,h]
//    FairCoin(n)      1/n                             n/2       n^2/12       N in [0,n)
//    Weighted(n,*w)   w[k]/(sum w)                    ...       ...          N in [0,n)

typedef void CDF;

CDF *Normal_CDF(double mean, double stdev);
CDF *Exponential_CDF(double a);
CDF *Poisson_CDF(double a);
CDF *Geometric_CDF(double p);
CDF *Binomial_CDF(int n, double p);
CDF *Uniform_CDF(double low, double hgh);
CDF *Fair_Coin_CDF(int n);
CDF *Weighted_Coin_CDF(int n, double *weight);

// Free all the space associated with generator cdf.

void Free_CDF(CDF *cdf);

// Generate the next random value for cdf.  The return value is double even if the
//   domain is a subset of the integers.

double Sample_CDF(CDF *cdf);

// Seed cdf with the specified value seedval.

void    Seed_CDF(CDF *cdf, uint32 seedval);

// Return a pointer to the generator of a cdf

uint64 *CDF_Generator(CDF *cdf);

//  Normally the current state of each CDF is independent of any other CDF objects.
//    One can link a CDF to another with Link_CDF so that CDF slave shares the same undelying
//    random number generator as CDF master.  One can later restore a CDF slave to using
//    its own random number generator by calling Unlink_CDF.

void   Link_CDF(CDF *master, CDF *slave);
void   Unlink_CDF(CDF *slave);

#endif
