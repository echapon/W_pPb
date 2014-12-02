#include "stattest.h"

using namespace std;

void stattest::set_data(vector<double> data)
{
   if (data.size() != m_nbins)
   {
      cout << "Error, wrong size for data vector" << endl;
      return;
   }
   m_data = data;
}

void stattest::set_pred(vector<double> pred)
{
   if (pred.size() != m_nbins)
   {
      cout << "Error, wrong size for pred vector" << endl;
      return;
   }
   m_pred = pred;
}

void stattest::set_systs(vector< vector<double> > systs)
{
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
}

void stattest::set_cov(vector< vector<double> > cov)
{
   if (cov.size() != m_nbins)
   {
      cout << "Error, wrong size for cov vector" << endl;
      return;
   }
   for (int i=0; i<cov.size(); i++)
      if (cov[i].size() != m_nbins)
      {
         cout << "Error, wrong size for cov vector" << endl;
         return;
      }
   m_cov = cov;
}

double stattest::chi2(vector<double> data, vector<double> prediction, vector< vector<double> > covariance)
{
   double val=0;
   if (data.size() != prediction.size() || data.size() != covariance.size()) return -1;
   int n = data.size();

   for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
         chi2 += (data[i]-prediction[i])*covariance[i][j]*(data[j]-prediction[j]);

   return chi2;
}

double stattest::log_poisson_likelihood(vector<double> data, vector<double> prediction, vector<double> pulls;)
{
   double val=0;
   if (data.size() != prediction.size()) return -1;
   int n = data.size();
   int nsyst = pulls.size();

   for (int i=0; i<n; i++)
      val += 2*(prediction[i]-data[i] - data[i]*log(prediction[i]/data[i]));

   for (int j=0; j<nsyst; j++)
      val += pulls[j]*pulls[j];

   return val;
}

double stattest::log_rationormal_likelihood(vector<double> data, vector<double> mu1, vector<double> s2, vector<double> mu2, vector<double> s2, vector<double> pulls;)
{
   double val=0;
   if (data.size() != prediction.size()) return -1;
   int n = data.size();
   int nsyst = pulls.size();

   for (int i=0; i<n; i++)
      val += 2*log(ratioOfNormalsDensity(data[i],mu1[i],s1[i],mu2[i],s2[i],rho[i]);

   for (int j=0; j<nsyst; j++)
      val += pulls[j]*pulls[j];

   return val;
}
