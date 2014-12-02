#include "statchannel.h"

#include "math.h"

using namespace std;


void statchannel::read(int n, const char* filename)
{
   m_nbins = n;
   ifstream file;
   file.open(filename);
   double x;

   getline(file,m_label);
   m_bins.resize(m_nbins+1);
   for (int i=0; i<m_nbins+1; i++) file >> m_bins[i];
   m_yields.resize(m_nbins);
   for (int i=0; i<m_nbins; i++) file >> m_yields[i];
   m_staterr.resize(m_nbins);
   for (int i=0; i<m_nbins; i++) file >> m_staterr[i];
   m_th.resize(m_nbins);
   for (int i=0; i<m_nbins; i++) file >> m_th[i];
   m_therr.resize(m_nbins);
   for (int i=0; i<m_nbins; i++) file >> m_therr[i];
   file >> x;
   while (!file.eof())
   {
      vector<double> tmpvec;
      tmpvec.resize(m_nbins);
      for (int i=0; i<m_nbins; i++) {tmpvec[i] = x; file >> x;};
      m_eff.push_back(tmpvec);
      for (int i=0; i<m_nbins; i++) {tmpvec[i] = x; file >> x;};
      m_efferr.push_back(tmpvec);
   }
}

void statchannel::print()
{
   cout << m_label << endl;
   cout << "bin edges: " << endl;
   for (int i=0; i<m_nbins+1; i++) cout << m_bins[i] << " ";
   cout << endl;
   cout << "yields: " << endl;
   for (int i=0; i<m_nbins; i++) cout << m_yields[i] << " ";
   cout << endl;
   cout << "statistical error on yields: " << endl;
   for (int i=0; i<m_nbins; i++) cout << m_staterr[i] << " ";
   cout << endl;
   cout << "theory: " << endl;
   for (int i=0; i<m_nbins; i++) cout << m_th[i] << " ";
   cout << endl;
   cout << "theory error: " << endl;
   for (int i=0; i<m_nbins; i++) cout << m_therr[i] << " ";
   cout << endl;
   for (int j=0; j<m_eff.size(); j++)
   {
      cout << "eff. " << j << ": " << endl;
      for (int i=0; i<m_nbins; i++) cout << m_eff[j][i] << " ";
      cout << endl;
      cout << "error on eff. " << j << ": " << endl;
      for (int i=0; i<m_nbins; i++) cout << m_efferr[j][i] << " ";
      cout << endl;
   }
   cout << endl;
}

TGraphErrors* statchannel::graph(const statchannel chanp, const statchannel chanm, asym mode, bool isth, bool dosyst)
{
   int nbins = chanp.get_nbins();
   vector<double> bins = chanp.get_bins();
   vector<double> th_p = chanp.get_th();
   vector<double> th_m = chanm.get_th();
   vector<double> therr_p = chanp.get_therr();
   vector<double> therr_m = chanm.get_therr();

   if (!isth)
   {
      th_p = chanp.get_yields();
      th_m = chanm.get_yields();
      therr_p = chanp.get_staterr();
      therr_m = chanm.get_staterr();
      vector< vector<double> > eff_p = chanp.get_eff();
      vector< vector<double> > eff_m = chanm.get_eff();
      vector< vector<double> > efferr_p = chanp.get_efferr();
      vector< vector<double> > efferr_m = chanm.get_efferr();

      // now compute the errors and the yields taking into account the efficiencies
      bool correl_pm = mode!=CH;
      int neff = eff_p.size();
      for (int i=0; i<nbins; i++)
      {
         // first, compute the yields from all the efficiencies
         for (int ieff=0; ieff<neff; ieff++)
         {
            th_p[i] /= eff_p[ieff][i];
            th_m[i] /= eff_m[ieff][i];
            // scale also the stat uncertainties
            therr_p[i] /= eff_p[ieff][i];
            therr_m[i] /= eff_m[ieff][i];
         }

         // now the error on the yield. we add everything in quadrature.
         // first the statistical error
         if (!dosyst || mode==YIELDS_PM)
         {
            therr_p[i] = pow(therr_p[i],2);
            therr_m[i] = pow(therr_m[i],2);
         }

         // now the systematics (be careful with correlations)
         else
         {
            therr_p[i] = 0; therr_m[i] = 0;
            for (int ieff=0; ieff<neff; ieff++)
            {
               if (correl_pm || efferr_p[ieff][i]>0)
                  therr_p[i] += pow(th_p[i]*efferr_p[ieff][i]/eff_p[ieff][i],2);
               if (correl_pm || efferr_m[ieff][i]>0)
                  therr_m[i] += pow(th_m[i]*efferr_m[ieff][i]/eff_m[ieff][i],2);
            }

            // // do not forget the luminosity
            // if (mode==YIELDS || mode==YIELDS_PM)
            // {
            //    therr_p[i] += pow(th_p[i]*LUMIERR,2);
            //    therr_m[i] += pow(th_m[i]*LUMIERR,2);
            // }
         }

         // at last, take the square root of the sum of the squares.
         therr_p[i] = sqrt(therr_p[i]);
         therr_m[i] = sqrt(therr_m[i]);
      }
   }

   double *x, *y, *ex, *ey;
   x = new double[nbins];
   y = new double[nbins];
   ex = new double[nbins];
   ey = new double[nbins];

   for (int i=0; i<nbins; i++)
   {
      // if we are not plotting the yields or the charge asymmetry, we only need nbins/2
      if (i>=(nbins/2) && !(mode==YIELDS || mode==YIELDS_PM || mode==CH)) break;

      int jp = (nbins/2)+i;
      int jm = (nbins/2)-i-1;

      // x axis: eta value
      int ix = (mode==YIELDS || mode==YIELDS_PM || mode==CH) ? i : jp;
      x[i] = (bins[ix+1]+bins[ix])/2.;
      ex[i] = (bins[ix+1]-bins[ix])/2.;
      // ex[i] = 0.05;

      // define some shorter names for later
      double Aa = th_p[i];
      double Ba = th_m[i];
      double dAa = therr_p[i];
      double dBa = therr_m[i];
      double A = th_p[jm]; double B = th_p[jp];
      double C = th_m[jm]; double D = th_m[jp];
      double dA = therr_p[jm]; double dB = therr_p[jp];
      double dC = therr_m[jm]; double dD = therr_m[jp];

      switch (mode)
      {
         case YIELDS : // W+ yields only
            {
               y[i] = th_p[i]/(bins[ix+1]-bins[ix]);
               ey[i] = therr_p[i]/(bins[ix+1]-bins[ix]);
               break;
            }

         case YIELDS_PM : // full yields (W+ + W-)
            {
               y[i] = (th_p[i]+th_m[i])/(bins[ix+1]-bins[ix]);
               // for the systematic error we need to be a bit careful with correlations
               if (!isth && dosyst)
               {
                  ey[i] = 0;
                  vector< vector<double> > eff_p = chanp.get_eff();
                  vector< vector<double> > eff_m = chanm.get_eff();
                  vector< vector<double> > efferr_p = chanp.get_efferr();
                  vector< vector<double> > efferr_m = chanm.get_efferr();
                  int neff = eff_p.size();
                  for (int ieff=0; ieff<neff; ieff++)
                  {
                     if (efferr_p[ieff][i]>0)
                        ey[i] += pow(th_p[i]*efferr_p[ieff][i]/eff_p[ieff][i],2) + pow(th_m[i]*efferr_m[ieff][i]/eff_m[ieff][i],2);
                     else
                        ey[i] += pow(th_p[i]*efferr_p[ieff][i]/eff_p[ieff][i] + th_m[i]*efferr_m[ieff][i]/eff_m[ieff][i],2);

                  }

                  // // do not forget the luminosity
                  // ey[i] += pow(LUMIERR*(th_p[i]+th_m[i]),2);

                  ey[i] = sqrt(ey[i]);
               }
               else
                  ey[i] = sqrt(pow(therr_p[i],2)+pow(therr_m[i],2));
               ey[i] = ey[i]/(bins[ix+1]-bins[ix]);
               break;
            }

         case CH : // charge asymmetry
            {
               y[i] = (Aa-Ba)/(Aa+Ba);
               ey[i] = sqrt((4*Ba*Ba/pow(Aa+Ba,4))*dAa*dAa + (4*Aa*Aa/pow(Aa+Ba,4)*dBa*dBa));
               break;
            }

         case A1p : // A1+
            {
               y[i] = A/B;
               ey[i] = fabs(y[i])*sqrt(pow(dA/A,2) + pow(dB/B,2));
               break;
            }

         case A1m : // A1-
            {
               y[i] = C/D;
               ey[i] = fabs(y[i])*sqrt(pow(dC/C,2) + pow(dD/D,2));
               break;
            }

         case A2 : // A2
            {
               y[i] = (A-B)/(C-D);
               ey[i] = fabs(y[i])*sqrt((dA*dA+dB*dB)/pow(A-B,2) + (dC*dC+dD*dD)/pow(C-D,2));
               break;
            }

         case A3 : // A3
            {
               y[i] = (A+C)/(B+D);
               ey[i] = fabs(y[i])*sqrt((dA*dA+dC*dC)/pow(A+C,2) + (dB*dB+dD*dD)/pow(B+D,2));
               break;
            }

         case A4 : // A4
            {
               y[i] = (A+C-B-D)/(total(nbins,th_p)+total(nbins,th_m));
               ey[i] = fabs(y[i])*sqrt((dA*dA+dB*dB+dC*dC+dD*dD)/pow(A+C-B-D,2) + (total2(nbins,therr_p)+total2(nbins,therr_m))/pow(total(nbins,th_p)+total(nbins,th_m),2));
               break;
            }

         default :
            {
               cout << "statchannel::graph_th: Error, unknown asymmetry type!"  << endl;
               return NULL;
            }
      }
   }

   if (!(mode==YIELDS || mode==YIELDS_PM || mode==CH)) nbins = nbins/2;
   TGraphErrors *gr = new TGraphErrors(nbins,x,y,ex,ey);
   return gr;
}

double statchannel::total(int n, vector<double> v)
{
   double tot=0;
   for (int i=0; i<n; i++) tot += v[i];
   return tot;
}

double statchannel::total2(int n, vector<double> v)
{
   double tot=0;
   for (int i=0; i<n; i++) tot += v[i]*v[i];
   return tot;
}

// TH2F* statchannel::covmat(const statchannel chanp, const statchannel chanm, asym mode, bool dosyst, bool dolumi)
// {
//    int nbins = chanp.get_nbins();
//    vector<double> bins = chanp.get_bins();
//    vector<double> th_p = chanp.get_th();
//    vector<double> th_m = chanm.get_th();
//    vector<double> therr_p = chanp.get_therr();
//    vector<double> therr_m = chanm.get_therr();

//    double *x, *y, *ex, *ey;
//    x = new double[nbins];
//    y = new double[nbins];
//    ex = new double[nbins];
//    ey = new double[nbins];

//    for (int i=0; i<nbins; i++)
//    {
//       // if we are not plotting the yields or the charge asymmetry, we only need nbins/2
//       if (i>=(nbins/2) && !(mode==YIELDS || mode==YIELDS_PM || mode==CH)) break;

//       int jp = (nbins/2)+i;
//       int jm = (nbins/2)-i-1;

//       // x axis: eta value
//       int ix = (mode==YIELDS || mode==YIELDS_PM || mode==CH) ? i : jp;
//       x[i] = (bins[ix+1]+bins[ix])/2.;
//       ex[i] = (bins[ix+1]-bins[ix])/2.;
//       // ex[i] = 0.05;

//       // define some shorter names for later
//       double Aa = th_p[i];
//       double Ba = th_m[i];
//       double dAa = therr_p[i];
//       double dBa = therr_m[i];
//       double A = th_p[jm]; double B = th_p[jp];
//       double C = th_m[jm]; double D = th_m[jp];
//       double dA = therr_p[jm]; double dB = therr_p[jp];
//       double dC = therr_m[jm]; double dD = therr_m[jp];

//       switch (mode)
//       {
//          case YIELDS : // W+ yields only
//             {
//                y[i] = th_p[i]/(bins[ix+1]-bins[ix]);
//                ey[i] = therr_p[i]/(bins[ix+1]-bins[ix]);
//                break;
//             }

//          case YIELDS_PM : // full yields (W+ + W-)
//             {
//                y[i] = (th_p[i]+th_m[i])/(bins[ix+1]-bins[ix]);
//                // for the systematic error we need to be a bit careful with correlations
//                if (!isth && dosyst)
//                {
//                   ey[i] = 0;
//                   vector< vector<double> > eff_p = chanp.get_eff();
//                   vector< vector<double> > eff_m = chanm.get_eff();
//                   vector< vector<double> > efferr_p = chanp.get_efferr();
//                   vector< vector<double> > efferr_m = chanm.get_efferr();
//                   int neff = eff_p.size();
//                   for (int ieff=0; ieff<neff; ieff++)
//                   {
//                      if (efferr_p[ieff][i]>0)
//                         ey[i] += pow(th_p[i]*efferr_p[ieff][i]/eff_p[ieff][i],2) + pow(th_m[i]*efferr_m[ieff][i]/eff_m[ieff][i],2);
//                      else
//                         ey[i] += pow(th_p[i]*efferr_p[ieff][i]/eff_p[ieff][i] + th_m[i]*efferr_m[ieff][i]/eff_m[ieff][i],2);

//                   }

//                   // // do not forget the luminosity
//                   // ey[i] += pow(LUMIERR*(th_p[i]+th_m[i]),2);

//                   ey[i] = sqrt(ey[i]);
//                }
//                else
//                   ey[i] = sqrt(pow(therr_p[i],2)+pow(therr_m[i],2));
//                ey[i] = ey[i]/(bins[ix+1]-bins[ix]);
//                break;
//             }

//          case CH : // charge asymmetry
//             {
//                y[i] = (Aa-Ba)/(Aa+Ba);
//                ey[i] = sqrt((4*Ba*Ba/pow(Aa+Ba,4))*dAa*dAa + (4*Aa*Aa/pow(Aa+Ba,4)*dBa*dBa));
//                break;
//             }

//          case A1p : // A1+
//             {
//                y[i] = A/B;
//                ey[i] = fabs(y[i])*sqrt(pow(dA/A,2) + pow(dB/B,2));
//                break;
//             }

//          case A1m : // A1-
//             {
//                y[i] = C/D;
//                ey[i] = fabs(y[i])*sqrt(pow(dC/C,2) + pow(dD/D,2));
//                break;
//             }

//          case A2 : // A2
//             {
//                y[i] = (A-B)/(C-D);
//                ey[i] = fabs(y[i])*sqrt((dA*dA+dB*dB)/pow(A-B,2) + (dC*dC+dD*dD)/pow(C-D,2));
//                break;
//             }

//          case A3 : // A3
//             {
//                y[i] = (A+C)/(B+D);
//                ey[i] = fabs(y[i])*sqrt((dA*dA+dC*dC)/pow(A+C,2) + (dB*dB+dD*dD)/pow(B+D,2));
//                break;
//             }

//          case A4 : // A4
//             {
//                y[i] = (A+C-B-D)/(total(nbins,th_p)+total(nbins,th_m));
//                ey[i] = fabs(y[i])*sqrt((dA*dA+dB*dB+dC*dC+dD*dD)/pow(A+C-B-D,2) + (total2(nbins,therr_p)+total2(nbins,therr_m))/pow(total(nbins,th_p)+total(nbins,th_m),2));
//                break;
//             }

//          default :
//             {
//                cout << "statchannel::graph_th: Error, unknown asymmetry type!"  << endl;
//                return NULL;
//             }
//       }
//    }

//    if (!(mode==YIELDS || mode==YIELDS_PM || mode==CH)) nbins = nbins/2;
//    TGraphErrors *gr = new TGraphErrors(nbins,x,y,ex,ey);
//    return gr;
// }
