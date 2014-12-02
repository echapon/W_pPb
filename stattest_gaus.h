#ifndef stattest_gaus_h
#define stattest_gaus_h

#include <vector>

using namespace std;

class stattest_gaus
{
   private:
      int m_nbins;
      int m_nsyst;
      vector<double> m_data; // size nbins
      vector<double> m_data_stat; // size nbins
      vector<double> m_pred; // size nbins
      vector< vector<double> > m_systs; // size [nsyst][nbins]
      vector< vector<double> > m_cov; // size [nbins][nbins]
      vector< vector<double> > m_invcov; // size [nbins][nbins]

   public:
      stattest_gaus(int nbins, int nsyst)
      {
         m_nbins = nbins; m_nsyst = nsyst;
         m_data.resize(nbins);
         m_data_stat.resize(nbins);
         m_pred.resize(nbins);
         m_cov.resize(nbins);
         m_invcov.resize(nbins);
         m_systs.resize(nsyst);
         for (int i=0; i<nsyst; i++) m_systs[i].resize(nbins);
         for (int i=0; i<nbins; i++) m_invcov[i].resize(nbins);
         for (int i=0; i<nbins; i++) m_cov[i].resize(nbins);
      };
      ~stattest_gaus() {};

      // setters
      void set_data(vector<double> data);
      void set_data_stat(vector<double> data_stat);
      void set_pred(vector<double> pred);
      void set_systs(vector< vector<double> > systs);

      // stat tests
      double chi2(double *pulls=NULL);
};

#endif // #ifndef stattest_gaus_h
