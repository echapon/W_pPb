#ifndef stattest_h
#define stattest_h

#include <vector>
#include "ratioofnormaldensity.h"

using namespace std;

class stattest
{
   private:
      int m_nbins;
      int m_nsyst;
      vector<double> m_data; // size nbins
      vector<double> m_pred; // size nbins
      vector< vector<double> > m_systs; // size [nsyst][nbins]
      vector< vector<double> > m_cov; // size [nbins][nbins]

   public:
      stattest(int nbins, int nsyst)
      {
         m_nbins = nbins; m_nsyst = nsyst;
         m_data.resize(nbins);
         m_pred.resize(nbins);
         m_cov.resize(nbins);
         m_systs.resize(nsyst);
         for (int i=0; i<nsyst; i++) m_systs[i].resize(nbins);
         for (int i=0; i<nbins; i++) m_cov[i].resize(nbins);
      };
      ~stattest() {};

      // setters
      void set_data(vector<double> data);
      void set_pred(vector<double> pred);
      void set_pulls(vector<double> pulls);
      void set_systs(vector< vector<double> > systs);
      void set_cov(vector< vector<double> > cov);

      // stat tests
      double chi2(double *pulls);
      double log_poisson_likelihood(double *pulls);
      double log_rationormal_likelihood(double *pulls);
}

#endif // #ifndef stattest_h
