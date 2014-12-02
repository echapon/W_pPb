#include "stattest_gaus.h"
#include <math.h>

using namespace std;

void stattest_gaus::set_data(vector<double> data)
{
   if (data.size() != m_nbins)
   {
      cout << "Error, wrong size for data vector" << endl;
      return;
   }
   m_data = data;
}

void stattest_gaus::set_data_stat(vector<double> data_stat)
{
   if (data_stat.size() != m_nbins)
   {
      cout << "Error, wrong size for data_stat vector" << endl;
      return;
   }
   m_data_stat = data_stat;
}

void stattest_gaus::set_pred(vector<double> pred)
{
   if (pred.size() != m_nbins)
   {
      cout << "Error, wrong size for pred vector" << endl;
      return;
   }
   m_pred = pred;
}

void stattest_gaus::set_systs(vector< vector<double> > systs)
{
   // set the systematics vector
   if (systs.size() != m_nsyst)
   {
      cout << "Error, wrong size for systs vector" << endl;
      return;
   }
   for (int i=0; i<systs.size(); i++)
      if (systs[i].size() != m_nbins)
      {
         cout << "Error, wrong size for systs vector" << endl;
         return;
      }
   m_systs = systs;

   // set the covariance matrix from there
   // cov = diag(m_data_stat^2) + sum_k ^t(fk) fk
   for (int i=0; i<m_nbins; i++)
      for (int j=0; j<m_nbins; j++)
      {
         if(i==j) m_cov[i][j] = pow(m_datastat[i],2);
         for (int k=0; k<m_nsyst; k++)
         {
            m_cov[i][j] += m_systs[k][i]*m_systs[k][j];
         }
      }

   // invert the covariance matrix and store it
   TMatrixDSym tm_cov(m_nbins), tm_invcov(m_nbins);
   for (int i=0; i<m_nbins; i++)
      for (int j=0; j<m_nbins; j++)
         tm_cov[i][j] = m_cov[i][j];

   tm_invconv = tm_cov.Invert();

   for (int i=0; i<m_nbins; i++)
      for (int j=0; j<m_nbins; j++)
         m_invcov[i][j] = tm_invcov[i][j];
}

double stattest_gaus::chi2(double *pulls)
{
   double val=0;

   // effect of systs
   vector<double> data_mod; data_mod = m_data;
   if (pulls != NULL)
   {
      for (int j=0; j<m_nbins; j++)
      {
         double fact=0.;
         for (int i=0; i<m_nsyst; i++)
            fact+=m_systs[i][j]*pulls[i];
         data_mod[j]*=(1+fact);
      }
   }

   for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
         val += (data_mod[i]-m_pred[i])*m_invcov[i][j]*(data_mod[j]-m_pred[j]);

   if (pulls != NULL)
      for (int i=0; i<m_nsysts; i++)
         val += pow(pulls[i],2);

   return val;
}
