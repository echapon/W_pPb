#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <iomanip>

#include "TRandom3.h"
#include "TGraph2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TMatrixDSym.h"
#include "TMath.h"
#include "TH2.h"
#include "combine.h"

#define NSET_CT10 26
#define NSET_EPS09 15

#define DOEFFS 1.
#define DOLUMI 1.
#define DOCT10 1.
#define DOEPS09 1.
#define DOSCALE 1.
#define DOSTAT 1.

#define APb 208

#define NPE 100000

#include "statchannel.h"

using namespace std;

double quad(double a, double b, double c, double x);
double chi2(int nbins, double *ni, double *dni, double *ti, vector<double> Sk);
double scalemod(double xsecs[7], double muf, double mur);
TGraph2D* scalegraph(double xsecs[7]);
double pdfmod(int nsets, double xsec0, double *xsecsp, double *xsecsm, double *pulls);
double mysqrt(double a);
double reverse(int nbins, double *array);
double reverse(int nbins1, int nbins2, double **array);
double total(int n, double *array);
double total(int n, vector<double> array);
double init(int n, double* array);
double init(int n1, int n2, double** array);

int main(int argc, const char** argv)
{
   if (argc != 19)
   {
      cout << "Usage: " << argv[0] << " nbins1 data_plus1.txt data_minus1.txt th_plus_nuc1.txt th_minus_nuc1.txt th_plus_nonuc1.txt th_minus_nonuc1.txt th_scale_plus1.txt th_minus_nonuc1.txt nbins2 data_plus2.txt data_minus2.txt th_plus_nuc2.txt th_minus_nuc2.txt th_plus_nonuc2.txt th_minus_nonuc2.txt th_scale_plus2.txt th_minus_nonuc2.txt" << endl;
      return 0;
   }

   // initialize random generator
   gRandom->SetSeed();

   TFile *fout = new TFile("cls.root","RECREATE");

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // 1. read in 2 (+, -) * 4 (exp, pdf uncert ct10, pdf uncert eps09, scale uncert) files
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   statchannel chan1p, chan1m, chan2m, chan2p;
   int nbins1;
   int nsyst1;
   int nbins2;
   int nsyst2;

   // channel 1
   nbins1 = atoi(argv[1]);
   chan1p.read(nbins1,argv[2]);
   chan1m.read(nbins1,argv[3]);
   // chan1p.print();
   // chan1m.print();
   nsyst1 = chan1p.get_eff().size();
   if (chan1m.get_eff().size() != nsyst1)
   {
      cout << "Error, " << argv[2] << " and " << argv[3] << " have a different number of efficiencies" << endl;
      return -1;
   }

   vector<double> bins1 = chan1p.get_bins();
   vector<double> yields1p = chan1p.get_yields();
   vector<double> yields1m = chan1m.get_yields();
   vector<double> staterr1p = chan1p.get_staterr();
   vector<double> staterr1m = chan1m.get_staterr();
   vector< vector<double> > eff1p = chan1p.get_eff();
   vector< vector<double> > eff1m = chan1m.get_eff();
   vector< vector<double> > efferr1p = chan1p.get_efferr();
   vector< vector<double> > efferr1m = chan1m.get_efferr();

   // channel 2
   nbins2 = atoi(argv[10]);
   chan2p.read(nbins2,argv[11]);
   chan2m.read(nbins2,argv[12]);
   // chan2p.print();
   // chan2m.print();
   nsyst2 = chan2p.get_eff().size();
   if (chan2m.get_eff().size() != nsyst2)
   {
      cout << "Error, " << argv[11] << " and " << argv[12] << " have a different number of efficiencies" << endl;
      return -1;
   }

   vector<double> bins2 = chan2p.get_bins();
   vector<double> yields2p = chan2p.get_yields();
   vector<double> yields2m = chan2m.get_yields();
   vector<double> staterr2p = chan2p.get_staterr();
   vector<double> staterr2m = chan2m.get_staterr();
   vector< vector<double> > eff2p = chan2p.get_eff();
   vector< vector<double> > eff2m = chan2m.get_eff();
   vector< vector<double> > efferr2p = chan2p.get_efferr();
   vector< vector<double> > efferr2m = chan2m.get_efferr();

   // some declarations
   double *etamin1 = new double[nbins1];
   double *etamax1 = new double[nbins1];
   double *etaavg1 = new double[nbins1];
   double *xsec_ch1p_nuc01 = new double[nbins1];
   double **xsec_ch1p_nucp1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_ch1p_nucp1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch1p_nucm1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_ch1p_nucm1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch1p_nucp1_eps09 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_ch1p_nucp1_eps09[ibin] = new double[NSET_EPS09];
   double **xsec_ch1p_nucm1_eps09 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_ch1p_nucm1_eps09[ibin] = new double[NSET_EPS09];
   double *xsec_ch1p_nonuc01 = new double[nbins1];
   double **xsec_ch1p_nonucp1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_ch1p_nonucp1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch1p_nonucm1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_ch1p_nonucm1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch1p_scale = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_ch1p_scale[ibin] = new double[7];
   double *xsec_ch1m_nuc01 = new double[nbins1];
   double **xsec_ch1m_nucp1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_ch1m_nucp1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch1m_nucm1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_ch1m_nucm1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch1m_nucp1_eps09 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_ch1m_nucp1_eps09[ibin] = new double[NSET_EPS09];
   double **xsec_ch1m_nucm1_eps09 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_ch1m_nucm1_eps09[ibin] = new double[NSET_EPS09];
   double *xsec_ch1m_nonuc01 = new double[nbins1];
   double **xsec_ch1m_nonucp1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_ch1m_nonucp1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch1m_nonucm1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_ch1m_nonucm1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch1m_scale = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_ch1m_scale[ibin] = new double[7];
   double *etamin2 = new double[nbins2];
   double *etamax2 = new double[nbins2];
   double *etaavg2 = new double[nbins2];
   double *xsec_ch2p_nuc01 = new double[nbins2];
   double **xsec_ch2p_nucp1_ct10 = new double*[nbins2];
   for (int ibin=0; ibin<nbins2; ibin++)
      xsec_ch2p_nucp1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch2p_nucm1_ct10 = new double*[nbins2];
   for (int ibin=0; ibin<nbins2; ibin++)
      xsec_ch2p_nucm1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch2p_nucp1_eps09 = new double*[nbins2];
   for (int ibin=0; ibin<nbins2; ibin++)
      xsec_ch2p_nucp1_eps09[ibin] = new double[NSET_EPS09];
   double **xsec_ch2p_nucm1_eps09 = new double*[nbins2];
   for (int ibin=0; ibin<nbins2; ibin++)
      xsec_ch2p_nucm1_eps09[ibin] = new double[NSET_EPS09];
   double *xsec_ch2p_nonuc01 = new double[nbins2];
   double **xsec_ch2p_nonucp1_ct10 = new double*[nbins2];
   for (int ibin=0; ibin<nbins2; ibin++)
      xsec_ch2p_nonucp1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch2p_nonucm1_ct10 = new double*[nbins2];
   for (int ibin=0; ibin<nbins2; ibin++)
      xsec_ch2p_nonucm1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch2p_scale = new double*[nbins2];
   for (int ibin=0; ibin<nbins2; ibin++)
      xsec_ch2p_scale[ibin] = new double[7];
   double *xsec_ch2m_nuc01 = new double[nbins2];
   double **xsec_ch2m_nucp1_ct10 = new double*[nbins2];
   for (int ibin=0; ibin<nbins2; ibin++)
      xsec_ch2m_nucp1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch2m_nucm1_ct10 = new double*[nbins2];
   for (int ibin=0; ibin<nbins2; ibin++)
      xsec_ch2m_nucm1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch2m_nucp1_eps09 = new double*[nbins2];
   for (int ibin=0; ibin<nbins2; ibin++)
      xsec_ch2m_nucp1_eps09[ibin] = new double[NSET_EPS09];
   double **xsec_ch2m_nucm1_eps09 = new double*[nbins2];
   for (int ibin=0; ibin<nbins2; ibin++)
      xsec_ch2m_nucm1_eps09[ibin] = new double[NSET_EPS09];
   double *xsec_ch2m_nonuc01 = new double[nbins2];
   double **xsec_ch2m_nonucp1_ct10 = new double*[nbins2];
   for (int ibin=0; ibin<nbins2; ibin++)
      xsec_ch2m_nonucp1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch2m_nonucm1_ct10 = new double*[nbins2];
   for (int ibin=0; ibin<nbins2; ibin++)
      xsec_ch2m_nonucm1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_ch2m_scale = new double*[nbins2];
   for (int ibin=0; ibin<nbins2; ibin++)
      xsec_ch2m_scale[ibin] = new double[7];

   // read the files
   ifstream file_ch1p_nuc(argv[4]);
   ifstream file_ch1m_nuc(argv[5]);
   ifstream file_ch1p_nonuc(argv[6]);
   ifstream file_ch1m_nonuc(argv[7]);
   ifstream file_ch1p_scale(argv[8]);
   ifstream file_ch1m_scale(argv[9]);
   ifstream file_ch2p_nuc(argv[13]);
   ifstream file_ch2m_nuc(argv[14]);
   ifstream file_ch2p_nonuc(argv[15]);
   ifstream file_ch2m_nonuc(argv[16]);
   ifstream file_ch2p_scale(argv[17]);
   ifstream file_ch2m_scale(argv[18]);
   for (int ibin=0; ibin<nbins1; ibin++)
   {
      file_ch1p_nuc >> etamin1[ibin] >> etamax1[ibin] >> etaavg1[ibin] >> xsec_ch1p_nuc01[ibin];
      for (int iset=0; iset<NSET_CT10; iset++)
         file_ch1p_nuc >> xsec_ch1p_nucp1_ct10[ibin][iset] >> xsec_ch1p_nucm1_ct10[ibin][iset];
      for (int iset=0; iset<NSET_EPS09; iset++)
         file_ch1p_nuc >> xsec_ch1p_nucp1_eps09[ibin][iset] >> xsec_ch1p_nucm1_eps09[ibin][iset];

      double dummy;
      file_ch1p_nonuc >> dummy >> dummy >> dummy >> xsec_ch1p_nonuc01[ibin];
      for (int iset=0; iset<NSET_CT10; iset++)
         file_ch1p_nonuc >> xsec_ch1p_nonucp1_ct10[ibin][iset] >> xsec_ch1p_nonucm1_ct10[ibin][iset];

      file_ch1p_scale >> dummy >> dummy >> dummy;
      for (int iset=0; iset<7; iset++)
         file_ch1p_scale >> xsec_ch1p_scale[ibin][iset];

      file_ch1m_nuc >> dummy >> dummy >> dummy >> xsec_ch1m_nuc01[ibin];
      for (int iset=0; iset<NSET_CT10; iset++)
         file_ch1m_nuc >> xsec_ch1m_nucp1_ct10[ibin][iset] >> xsec_ch1m_nucm1_ct10[ibin][iset];
      for (int iset=0; iset<NSET_EPS09; iset++)
         file_ch1m_nuc >> xsec_ch1m_nucp1_eps09[ibin][iset] >> xsec_ch1m_nucm1_eps09[ibin][iset];

      file_ch1m_nonuc >> dummy >> dummy >> dummy >> xsec_ch1m_nonuc01[ibin];
      for (int iset=0; iset<NSET_CT10; iset++)
         file_ch1m_nonuc >> xsec_ch1m_nonucp1_ct10[ibin][iset] >> xsec_ch1m_nonucm1_ct10[ibin][iset];

      file_ch1m_scale >> dummy >> dummy >> dummy;
      for (int iset=0; iset<7; iset++)
         file_ch1m_scale >> xsec_ch1m_scale[ibin][iset];
   }
   for (int ibin=0; ibin<nbins2; ibin++)
   {
      file_ch2p_nuc >> etamin2[ibin] >> etamax2[ibin] >> etaavg2[ibin] >> xsec_ch2p_nuc01[ibin];
      for (int iset=0; iset<NSET_CT10; iset++)
         file_ch2p_nuc >> xsec_ch2p_nucp1_ct10[ibin][iset] >> xsec_ch2p_nucm1_ct10[ibin][iset];
      for (int iset=0; iset<NSET_EPS09; iset++)
         file_ch2p_nuc >> xsec_ch2p_nucp1_eps09[ibin][iset] >> xsec_ch2p_nucm1_eps09[ibin][iset];

      double dummy;
      file_ch2p_nonuc >> dummy >> dummy >> dummy >> xsec_ch2p_nonuc01[ibin];
      for (int iset=0; iset<NSET_CT10; iset++)
         file_ch2p_nonuc >> xsec_ch2p_nonucp1_ct10[ibin][iset] >> xsec_ch2p_nonucm1_ct10[ibin][iset];

      file_ch2p_scale >> dummy >> dummy >> dummy;
      for (int iset=0; iset<7; iset++)
         file_ch2p_scale >> xsec_ch2p_scale[ibin][iset];

      file_ch2m_nuc >> dummy >> dummy >> dummy >> xsec_ch2m_nuc01[ibin];
      for (int iset=0; iset<NSET_CT10; iset++)
         file_ch2m_nuc >> xsec_ch2m_nucp1_ct10[ibin][iset] >> xsec_ch2m_nucm1_ct10[ibin][iset];
      for (int iset=0; iset<NSET_EPS09; iset++)
         file_ch2m_nuc >> xsec_ch2m_nucp1_eps09[ibin][iset] >> xsec_ch2m_nucm1_eps09[ibin][iset];

      file_ch2m_nonuc >> dummy >> dummy >> dummy >> xsec_ch2m_nonuc01[ibin];
      for (int iset=0; iset<NSET_CT10; iset++)
         file_ch2m_nonuc >> xsec_ch2m_nonucp1_ct10[ibin][iset] >> xsec_ch2m_nonucm1_ct10[ibin][iset];

      file_ch2m_scale >> dummy >> dummy >> dummy;
      for (int iset=0; iset<7; iset++)
         file_ch2m_scale >> xsec_ch2m_scale[ibin][iset];
   }

   // theory is not in the same eta convention as data!!! reverse it
   reverse(nbins1,xsec_ch1p_nuc01);
   reverse(nbins1,NSET_CT10,xsec_ch1p_nucp1_ct10);
   reverse(nbins1,NSET_CT10,xsec_ch1p_nucm1_ct10);
   reverse(nbins1,NSET_EPS09,xsec_ch1p_nucp1_eps09);
   reverse(nbins1,NSET_EPS09,xsec_ch1p_nucm1_eps09);
   reverse(nbins1,xsec_ch1p_nonuc01);
   reverse(nbins1,NSET_CT10,xsec_ch1p_nonucp1_ct10);
   reverse(nbins1,NSET_CT10,xsec_ch1p_nonucm1_ct10);
   reverse(nbins1,7,xsec_ch1p_scale);
   reverse(nbins1,xsec_ch1m_nuc01);
   reverse(nbins1,NSET_CT10,xsec_ch1m_nucp1_ct10);
   reverse(nbins1,NSET_CT10,xsec_ch1m_nucm1_ct10);
   reverse(nbins1,NSET_EPS09,xsec_ch1m_nucp1_eps09);
   reverse(nbins1,NSET_EPS09,xsec_ch1m_nucm1_eps09);
   reverse(nbins1,xsec_ch1m_nonuc01);
   reverse(nbins1,NSET_CT10,xsec_ch1m_nonucp1_ct10);
   reverse(nbins1,NSET_CT10,xsec_ch1m_nonucm1_ct10);
   reverse(nbins1,7,xsec_ch1m_scale);
   reverse(nbins2,xsec_ch2p_nuc01);
   reverse(nbins2,NSET_CT10,xsec_ch2p_nucp1_ct10);
   reverse(nbins2,NSET_CT10,xsec_ch2p_nucm1_ct10);
   reverse(nbins2,NSET_EPS09,xsec_ch2p_nucp1_eps09);
   reverse(nbins2,NSET_EPS09,xsec_ch2p_nucm1_eps09);
   reverse(nbins2,xsec_ch2p_nonuc01);
   reverse(nbins2,NSET_CT10,xsec_ch2p_nonucp1_ct10);
   reverse(nbins2,NSET_CT10,xsec_ch2p_nonucm1_ct10);
   reverse(nbins2,7,xsec_ch2p_scale);
   reverse(nbins2,xsec_ch2m_nuc01);
   reverse(nbins2,NSET_CT10,xsec_ch2m_nucp1_ct10);
   reverse(nbins2,NSET_CT10,xsec_ch2m_nucm1_ct10);
   reverse(nbins2,NSET_EPS09,xsec_ch2m_nucp1_eps09);
   reverse(nbins2,NSET_EPS09,xsec_ch2m_nucm1_eps09);
   reverse(nbins2,xsec_ch2m_nonuc01);
   reverse(nbins2,NSET_CT10,xsec_ch2m_nonucp1_ct10);
   reverse(nbins2,NSET_CT10,xsec_ch2m_nonucm1_ct10);
   reverse(nbins2,7,xsec_ch2m_scale);

   // // check xsecs scales
   // for (int ibin=0; ibin<nbins1; ibin++)
   // {
   //    cout << "bin " << ibin << ": ";
   //    for (int iset=0; iset<7; iset++)
   //       cout << xsec_mup_scale[ibin][iset] << " ";
   //    cout << endl;
   // }
   // for (int ibin=0; ibin<nbins1; ibin++)
   // {
   //    cout << "bin " << ibin << ": ";
   //    for (int iset=0; iset<7; iset++)
   //       cout << xsec_mum_scale[ibin][iset] << " ";
   //    cout << endl;
   // }

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // compute quantities for data
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   double yp_data_ch1[nbins1], ym_data_ch1[nbins1], yt_data_ch1[2*nbins1];
   double A1p_data_ch1[nbins1/2], A1m_data_ch1[nbins1/2], A3_data_ch1[nbins1/2], A4_data_ch1[nbins1/2];;
   double Ch_data_ch1[nbins1];
   double yp_data_ch2[nbins2], ym_data_ch2[nbins2], yt_data_ch2[2*nbins2];
   double A1p_data_ch2[nbins2/2], A1m_data_ch2[nbins2/2], A3_data_ch2[nbins2/2], A4_data_ch2[nbins1/2];
   double Ch_data_ch2[nbins2];
   double yp_data_ch12[nbins1+nbins2], ym_data_ch12[nbins1+nbins2], yt_data_ch12[2*(nbins1+nbins2)];
   double A1p_data_ch12[(nbins1+nbins2)/2], A1m_data_ch12[(nbins1+nbins2)/2], A3_data_ch12[(nbins1+nbins2)/2], A4_data_ch12[(nbins1+nbins2)/2];
   double Ch_data_ch12[nbins1+nbins2];
   for (int ibin=0; ibin<nbins1/2; ibin++)
   {
      int jp = (nbins1/2)+ibin;
      int jm = (nbins1/2)-ibin-1;

      double A_data_ch1 = yields1p[jm]; double B_data_ch1 = yields1p[jp];
      double C_data_ch1 = yields1m[jm]; double D_data_ch1 = yields1m[jp];

      for (int ieff=0; ieff<eff1p.size(); ieff++)
      {
         A_data_ch1 /= eff1p[ieff][jm];
         B_data_ch1 /= eff1p[ieff][jp];
         C_data_ch1 /= eff1m[ieff][jm];
         D_data_ch1 /= eff1m[ieff][jp];
      }

      yp_data_ch1[jm] = A_data_ch1;
      yp_data_ch1[jp] = B_data_ch1;
      ym_data_ch1[jm] = C_data_ch1;
      ym_data_ch1[jp] = D_data_ch1;
      yt_data_ch1[jm] = A_data_ch1;
      yt_data_ch1[jp] = B_data_ch1;
      yt_data_ch1[jm+nbins1] = C_data_ch1;
      yt_data_ch1[jp+nbins1] = D_data_ch1;
      Ch_data_ch1[jp] = (B_data_ch1-D_data_ch1)/(B_data_ch1+D_data_ch1);
      Ch_data_ch1[jm] = (A_data_ch1-C_data_ch1)/(A_data_ch1+C_data_ch1);
      A1p_data_ch1[ibin] = A_data_ch1/B_data_ch1;
      A1m_data_ch1[ibin] = C_data_ch1/D_data_ch1;
      A3_data_ch1[ibin] = (A_data_ch1+C_data_ch1)/(B_data_ch1+D_data_ch1);
   }
   for (int ibin=0; ibin<nbins1/2; ibin++)
   {
      int jp = (nbins1/2)+ibin;
      int jm = (nbins1/2)-ibin-1;
      double A_data_ch1 = yields1p[jm]; double B_data_ch1 = yields1p[jp];
      double C_data_ch1 = yields1m[jm]; double D_data_ch1 = yields1m[jp];
      for (int ieff=0; ieff<eff1p.size(); ieff++)
      {
         A_data_ch1 /= eff1p[ieff][jm];
         B_data_ch1 /= eff1p[ieff][jp];
         C_data_ch1 /= eff1m[ieff][jm];
         D_data_ch1 /= eff1m[ieff][jp];
      }
      A4_data_ch1[ibin] = (A_data_ch1+C_data_ch1-B_data_ch1-D_data_ch1)/(total(nbins1,yp_data_ch1)+total(nbins1,ym_data_ch1));
   }
   for (int ibin=0; ibin<nbins2/2; ibin++)
   {
      int jp = (nbins2/2)+ibin;
      int jm = (nbins2/2)-ibin-1;

      double A_data_ch2 = yields2p[jm]; double B_data_ch2 = yields2p[jp];
      double C_data_ch2 = yields2m[jm]; double D_data_ch2 = yields2m[jp];

      for (int ieff=0; ieff<eff2p.size(); ieff++)
      {
         A_data_ch2 /= eff2p[ieff][jm];
         B_data_ch2 /= eff2p[ieff][jp];
         C_data_ch2 /= eff2m[ieff][jm];
         D_data_ch2 /= eff2m[ieff][jp];
      }

      yp_data_ch2[jm] = A_data_ch2;
      yp_data_ch2[jp] = B_data_ch2;
      ym_data_ch2[jm] = C_data_ch2;
      ym_data_ch2[jp] = D_data_ch2;
      yt_data_ch2[jm] = A_data_ch2;
      yt_data_ch2[jp] = B_data_ch2;
      yt_data_ch2[jm+nbins2] = C_data_ch2;
      yt_data_ch2[jp+nbins2] = D_data_ch2;
      Ch_data_ch2[jp] = (B_data_ch2-D_data_ch2)/(B_data_ch2+D_data_ch2);
      Ch_data_ch2[jm] = (A_data_ch2-C_data_ch2)/(A_data_ch2+C_data_ch2);
      A1p_data_ch2[ibin] = A_data_ch2/B_data_ch2;
      A1m_data_ch2[ibin] = C_data_ch2/D_data_ch2;
      A3_data_ch2[ibin] = (A_data_ch2+C_data_ch2)/(B_data_ch2+D_data_ch2);
   }
   for (int ibin=0; ibin<nbins2/2; ibin++)
   {
      int jp = (nbins2/2)+ibin;
      int jm = (nbins2/2)-ibin-1;
      double A_data_ch2 = yields2p[jm]; double B_data_ch2 = yields2p[jp];
      double C_data_ch2 = yields2m[jm]; double D_data_ch2 = yields2m[jp];
      for (int ieff=0; ieff<eff2p.size(); ieff++)
      {
         A_data_ch2 /= eff2p[ieff][jm];
         B_data_ch2 /= eff2p[ieff][jp];
         C_data_ch2 /= eff2m[ieff][jm];
         D_data_ch2 /= eff2m[ieff][jp];
      }
      A4_data_ch2[ibin] = (A_data_ch2+C_data_ch2-B_data_ch2-D_data_ch2)/(total(nbins2,yp_data_ch2)+total(nbins2,ym_data_ch2));
   }
   cat(nbins1,yp_data_ch1,nbins2,yp_data_ch2,yp_data_ch12);
   cat(nbins1,ym_data_ch1,nbins2,ym_data_ch2,ym_data_ch12);
   cat(2*nbins1,yt_data_ch1,2*nbins2,yt_data_ch2,yt_data_ch12);
   cat(nbins1,Ch_data_ch1,nbins2,Ch_data_ch2,Ch_data_ch12);
   cat(nbins1/2,A1p_data_ch1,nbins2/2,A1p_data_ch2,A1p_data_ch12);
   cat(nbins1/2,A1m_data_ch1,nbins2/2,A1m_data_ch2,A1m_data_ch12);
   cat(nbins1/2,A3_data_ch1,nbins2/2,A3_data_ch2,A3_data_ch12);
   cat(nbins1/2,A4_data_ch1,nbins2/2,A4_data_ch2,A4_data_ch12);

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // compute mean expectation on raw yields from theory
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   double *avg_yp_nuc_ch1 = new double[nbins1];
   double *avg_ym_nuc_ch1 = new double[nbins1];
   double *avg_yt_nuc_ch1 = new double[2*nbins1];
   double *avg_ch_nuc_ch1 = new double[nbins1];
   double *avg_A1p_nuc_ch1 = new double[nbins1/2];
   double *avg_A1m_nuc_ch1 = new double[nbins1/2];
   double *avg_A3_nuc_ch1 = new double[nbins1/2];
   double *avg_A4_nuc_ch1 = new double[nbins1/2];
   double *avg_yp_nonuc_ch1 = new double[nbins1];
   double *avg_ym_nonuc_ch1 = new double[nbins1];
   double *avg_yt_nonuc_ch1 = new double[2*nbins1];
   double *avg_ch_nonuc_ch1 = new double[nbins1];
   double *avg_A1p_nonuc_ch1 = new double[nbins1/2];
   double *avg_A1m_nonuc_ch1 = new double[nbins1/2];
   double *avg_A3_nonuc_ch1 = new double[nbins1/2];
   double *avg_A4_nonuc_ch1 = new double[nbins1/2];
   double *avg_yp_dat_ch1 = new double[nbins1];
   double *avg_ym_dat_ch1 = new double[nbins1];
   double *avg_yt_dat_ch1 = new double[2*nbins1];
   double *avg_ch_dat_ch1 = new double[nbins1];
   double *avg_A1p_dat_ch1 = new double[nbins1/2];
   double *avg_A1m_dat_ch1 = new double[nbins1/2];
   double *avg_A3_dat_ch1 = new double[nbins1/2];
   double *avg_A4_dat_ch1 = new double[nbins1/2];
   double *avg_yp_nuc_ch2 = new double[nbins2];
   double *avg_ym_nuc_ch2 = new double[nbins2];
   double *avg_yt_nuc_ch2 = new double[2*nbins2];
   double *avg_ch_nuc_ch2 = new double[nbins2];
   double *avg_A1p_nuc_ch2 = new double[nbins2/2];
   double *avg_A1m_nuc_ch2 = new double[nbins2/2];
   double *avg_A3_nuc_ch2 = new double[nbins2/2];
   double *avg_A4_nuc_ch2 = new double[nbins2/2];
   double *avg_yp_nonuc_ch2 = new double[nbins2];
   double *avg_ym_nonuc_ch2 = new double[nbins2];
   double *avg_yt_nonuc_ch2 = new double[2*nbins2];
   double *avg_ch_nonuc_ch2 = new double[nbins2];
   double *avg_A1p_nonuc_ch2 = new double[nbins2/2];
   double *avg_A1m_nonuc_ch2 = new double[nbins2/2];
   double *avg_A3_nonuc_ch2 = new double[nbins2/2];
   double *avg_A4_nonuc_ch2 = new double[nbins2/2];
   double *avg_yp_dat_ch2 = new double[nbins2];
   double *avg_ym_dat_ch2 = new double[nbins2];
   double *avg_yt_dat_ch2 = new double[2*nbins2];
   double *avg_ch_dat_ch2 = new double[nbins2];
   double *avg_A1p_dat_ch2 = new double[nbins2/2];
   double *avg_A1m_dat_ch2 = new double[nbins2/2];
   double *avg_A3_dat_ch2 = new double[nbins2/2];
   double *avg_A4_dat_ch2 = new double[nbins2/2];
   double *avg_yp_nuc_ch12 = new double[(nbins1+nbins2)];
   double *avg_ym_nuc_ch12 = new double[(nbins1+nbins2)];
   double *avg_yt_nuc_ch12 = new double[2*(nbins1+nbins2)];
   double *avg_ch_nuc_ch12 = new double[(nbins1+nbins2)];
   double *avg_A1p_nuc_ch12 = new double[(nbins1+nbins2)/2];
   double *avg_A1m_nuc_ch12 = new double[(nbins1+nbins2)/2];
   double *avg_A3_nuc_ch12 = new double[(nbins1+nbins2)/2];
   double *avg_A4_nuc_ch12 = new double[(nbins1+nbins2)/2];
   double *avg_yp_nonuc_ch12 = new double[(nbins1+nbins2)];
   double *avg_ym_nonuc_ch12 = new double[(nbins1+nbins2)];
   double *avg_yt_nonuc_ch12 = new double[2*(nbins1+nbins2)];
   double *avg_ch_nonuc_ch12 = new double[(nbins1+nbins2)];
   double *avg_A1p_nonuc_ch12 = new double[(nbins1+nbins2)/2];
   double *avg_A1m_nonuc_ch12 = new double[(nbins1+nbins2)/2];
   double *avg_A3_nonuc_ch12 = new double[(nbins1+nbins2)/2];
   double *avg_A4_nonuc_ch12 = new double[(nbins1+nbins2)/2];
   double *avg_yp_dat_ch12 = new double[(nbins1+nbins2)];
   double *avg_ym_dat_ch12 = new double[(nbins1+nbins2)];
   double *avg_yt_dat_ch12 = new double[2*(nbins1+nbins2)];
   double *avg_ch_dat_ch12 = new double[(nbins1+nbins2)];
   double *avg_A1p_dat_ch12 = new double[(nbins1+nbins2)/2];
   double *avg_A1m_dat_ch12 = new double[(nbins1+nbins2)/2];
   double *avg_A3_dat_ch12 = new double[(nbins1+nbins2)/2];
   double *avg_A4_dat_ch12 = new double[(nbins1+nbins2)/2];

   for (int ibin=0; ibin<nbins1/2; ibin++)
   {
      int jp = (nbins1/2)+ibin;
      int jm = (nbins1/2)-ibin-1;

      double A_nuc = xsec_ch1p_nuc01[jm]*APb*LUMI*1e-6; double B_nuc = xsec_ch1p_nuc01[jp]*APb*LUMI*1e-6;
      double C_nuc = xsec_ch1m_nuc01[jm]*APb*LUMI*1e-6; double D_nuc = xsec_ch1m_nuc01[jp]*APb*LUMI*1e-6;
      double A_nonuc = xsec_ch1p_nonuc01[jm]*APb*LUMI*1e-6; double B_nonuc = xsec_ch1p_nonuc01[jp]*APb*LUMI*1e-6;
      double C_nonuc = xsec_ch1m_nonuc01[jm]*APb*LUMI*1e-6; double D_nonuc = xsec_ch1m_nonuc01[jp]*APb*LUMI*1e-6;
      double A_dat = yields1p[jm]/(etamax1[jm]-etamin1[jm]); double B_dat = yields1p[jp]/(etamax1[jp]-etamin1[jp]);
      double C_dat = yields1m[jm]/(etamax1[jm]-etamin1[jm]); double D_dat = yields1m[jp]/(etamax1[jp]-etamin1[jp]);
      for (int ieff=0; ieff<eff1p.size(); ieff++)
      {
         A_dat /= eff1p[ieff][jm];
         B_dat /= eff1p[ieff][jp];
         C_dat /= eff1m[ieff][jm];
         D_dat /= eff1m[ieff][jp];
      }

      avg_yp_nuc_ch1[jm] = A_nuc;
      avg_yp_nuc_ch1[jp] = B_nuc;
      avg_ym_nuc_ch1[jm] = C_nuc;
      avg_ym_nuc_ch1[jp] = D_nuc;
      avg_yt_nuc_ch1[jm] = A_nuc;
      avg_yt_nuc_ch1[jp] = B_nuc;
      avg_yt_nuc_ch1[jm+nbins1] = C_nuc;
      avg_yt_nuc_ch1[jp+nbins1] = D_nuc;
      avg_ch_nuc_ch1[jp] = (B_nuc-D_nuc)/(B_nuc+D_nuc);
      avg_ch_nuc_ch1[jm] = (A_nuc-C_nuc)/(A_nuc+C_nuc);
      avg_A1p_nuc_ch1[ibin] = A_nuc/B_nuc;
      avg_A1m_nuc_ch1[ibin] = C_nuc/D_nuc;
      avg_A3_nuc_ch1[ibin] = (A_nuc+C_nuc)/(B_nuc+D_nuc);
      avg_yp_nonuc_ch1[jm] = A_nonuc;
      avg_yp_nonuc_ch1[jp] = B_nonuc;
      avg_ym_nonuc_ch1[jm] = C_nonuc;
      avg_ym_nonuc_ch1[jp] = D_nonuc;
      avg_yt_nonuc_ch1[jm] = A_nonuc;
      avg_yt_nonuc_ch1[jp] = B_nonuc;
      avg_yt_nonuc_ch1[jm+nbins1] = C_nonuc;
      avg_yt_nonuc_ch1[jp+nbins1] = D_nonuc;
      avg_ch_nonuc_ch1[jp] = (B_nonuc-D_nonuc)/(B_nonuc+D_nonuc);
      avg_ch_nonuc_ch1[jm] = (A_nonuc-C_nonuc)/(A_nonuc+C_nonuc);
      avg_A1p_nonuc_ch1[ibin] = A_nonuc/B_nonuc;
      avg_A1m_nonuc_ch1[ibin] = C_nonuc/D_nonuc;
      avg_A3_nonuc_ch1[ibin] = (A_nonuc+C_nonuc)/(B_nonuc+D_nonuc);
      avg_yp_dat_ch1[jm] = A_dat;
      avg_yp_dat_ch1[jp] = B_dat;
      avg_ym_dat_ch1[jm] = C_dat;
      avg_ym_dat_ch1[jp] = D_dat;
      avg_yt_dat_ch1[jm] = A_dat;
      avg_yt_dat_ch1[jp] = B_dat;
      avg_yt_dat_ch1[jm+nbins1] = C_dat;
      avg_yt_dat_ch1[jp+nbins1] = D_dat;
      avg_ch_dat_ch1[jp] = (B_dat-D_dat)/(B_dat+D_dat);
      avg_ch_dat_ch1[jm] = (A_dat-C_dat)/(A_dat+C_dat);
      avg_A1p_dat_ch1[ibin] = A_dat/B_dat;
      avg_A1m_dat_ch1[ibin] = C_dat/D_dat;
      avg_A3_dat_ch1[ibin] = (A_dat+C_dat)/(B_dat+D_dat);
      cout << jm << " " << avg_yp_dat_ch1[jm] << " " << avg_ym_dat_ch1[jm] << endl;
      cout << jp << " " << avg_yp_dat_ch1[jp] << " " << avg_ym_dat_ch1[jp] << endl;
      cout << "A3 " << ibin << " " << avg_A3_dat_ch1[ibin] << endl;
   }
   for (int ibin=0; ibin<nbins1/2; ibin++)
   {
      int jp = (nbins1/2)+ibin;
      int jm = (nbins1/2)-ibin-1;
      double A_nuc = xsec_ch1p_nuc01[jm]*APb*LUMI*1e-6; double B_nuc = xsec_ch1p_nuc01[jp]*APb*LUMI*1e-6;
      double C_nuc = xsec_ch1m_nuc01[jm]*APb*LUMI*1e-6; double D_nuc = xsec_ch1m_nuc01[jp]*APb*LUMI*1e-6;
      double A_nonuc = xsec_ch1p_nonuc01[jm]*APb*LUMI*1e-6; double B_nonuc = xsec_ch1p_nonuc01[jp]*APb*LUMI*1e-6;
      double C_nonuc = xsec_ch1m_nonuc01[jm]*APb*LUMI*1e-6; double D_nonuc = xsec_ch1m_nonuc01[jp]*APb*LUMI*1e-6;
      double A_dat = yields1p[jm]/(etamax1[jm]-etamin1[jm]); double B_dat = yields1p[jp]/(etamax1[jp]-etamin1[jp]);
      double C_dat = yields1m[jm]/(etamax1[jm]-etamin1[jm]); double D_dat = yields1m[jp]/(etamax1[jp]-etamin1[jp]);
      for (int ieff=0; ieff<eff1p.size(); ieff++)
      {
         A_dat /= eff1p[ieff][jm];
         B_dat /= eff1p[ieff][jp];
         C_dat /= eff1m[ieff][jm];
         D_dat /= eff1m[ieff][jp];
      }
      avg_A4_nuc_ch1[ibin] = (A_nuc+C_nuc-B_nuc-D_nuc)/(total(nbins1,avg_yp_nuc_ch1)+total(nbins1,avg_ym_nuc_ch1));
      avg_A4_nonuc_ch1[ibin] = (A_nonuc+C_nonuc-B_nonuc-D_nonuc)/(total(nbins1,avg_yp_nonuc_ch1)+total(nbins1,avg_ym_nonuc_ch1));
      avg_A4_dat_ch1[ibin] = (A_dat+C_dat-B_dat-D_dat)/(total(nbins1,avg_yp_dat_ch1)+total(nbins1,avg_ym_dat_ch1));
   }

   for (int ibin=0; ibin<nbins2/2; ibin++)
   {
      int jp = (nbins2/2)+ibin;
      int jm = (nbins2/2)-ibin-1;

      double A_nuc = xsec_ch2p_nuc01[jm]*APb*LUMI*1e-6; double B_nuc = xsec_ch2p_nuc01[jp]*APb*LUMI*1e-6;
      double C_nuc = xsec_ch2m_nuc01[jm]*APb*LUMI*1e-6; double D_nuc = xsec_ch2m_nuc01[jp]*APb*LUMI*1e-6;
      double A_nonuc = xsec_ch2p_nonuc01[jm]*APb*LUMI*1e-6; double B_nonuc = xsec_ch2p_nonuc01[jp]*APb*LUMI*1e-6;
      double C_nonuc = xsec_ch2m_nonuc01[jm]*APb*LUMI*1e-6; double D_nonuc = xsec_ch2m_nonuc01[jp]*APb*LUMI*1e-6;
      double A_dat = yields2p[jm]/(etamax2[jm]-etamin2[jm]); double B_dat = yields2p[jp]/(etamax2[jp]-etamin2[jp]);
      double C_dat = yields2m[jm]/(etamax2[jm]-etamin2[jm]); double D_dat = yields2m[jp]/(etamax2[jp]-etamin2[jp]);
      for (int ieff=0; ieff<eff2p.size(); ieff++)
      {
         A_dat /= eff2p[ieff][jm];
         B_dat /= eff2p[ieff][jp];
         C_dat /= eff2m[ieff][jm];
         D_dat /= eff2m[ieff][jp];
      }

      avg_yp_nuc_ch2[jm] = A_nuc;
      avg_yp_nuc_ch2[jp] = B_nuc;
      avg_ym_nuc_ch2[jm] = C_nuc;
      avg_ym_nuc_ch2[jp] = D_nuc;
      avg_yt_nuc_ch2[jm] = A_nuc;
      avg_yt_nuc_ch2[jp] = B_nuc;
      avg_yt_nuc_ch2[jm+nbins2] = C_nuc;
      avg_yt_nuc_ch2[jp+nbins2] = D_nuc;
      avg_ch_nuc_ch2[jp] = (B_nuc-D_nuc)/(B_nuc+D_nuc);
      avg_ch_nuc_ch2[jm] = (A_nuc-C_nuc)/(A_nuc+C_nuc);
      avg_A1p_nuc_ch2[ibin] = A_nuc/B_nuc;
      avg_A1m_nuc_ch2[ibin] = C_nuc/D_nuc;
      avg_A3_nuc_ch2[ibin] = (A_nuc+C_nuc)/(B_nuc+D_nuc);
      avg_yp_nonuc_ch2[jm] = A_nonuc;
      avg_yp_nonuc_ch2[jp] = B_nonuc;
      avg_ym_nonuc_ch2[jm] = C_nonuc;
      avg_ym_nonuc_ch2[jp] = D_nonuc;
      avg_yt_nonuc_ch2[jm] = A_nonuc;
      avg_yt_nonuc_ch2[jp] = B_nonuc;
      avg_yt_nonuc_ch2[jm+nbins2] = C_nonuc;
      avg_yt_nonuc_ch2[jp+nbins2] = D_nonuc;
      avg_ch_nonuc_ch2[jp] = (B_nonuc-D_nonuc)/(B_nonuc+D_nonuc);
      avg_ch_nonuc_ch2[jm] = (A_nonuc-C_nonuc)/(A_nonuc+C_nonuc);
      avg_A1p_nonuc_ch2[ibin] = A_nonuc/B_nonuc;
      avg_A1m_nonuc_ch2[ibin] = C_nonuc/D_nonuc;
      avg_A3_nonuc_ch2[ibin] = (A_nonuc+C_nonuc)/(B_nonuc+D_nonuc);
      avg_yp_dat_ch2[jm] = A_dat;
      avg_yp_dat_ch2[jp] = B_dat;
      avg_ym_dat_ch2[jm] = C_dat;
      avg_ym_dat_ch2[jp] = D_dat;
      avg_yt_dat_ch2[jm] = A_dat;
      avg_yt_dat_ch2[jp] = B_dat;
      avg_yt_dat_ch2[jm+nbins2] = C_dat;
      avg_yt_dat_ch2[jp+nbins2] = D_dat;
      avg_ch_dat_ch2[jp] = (B_dat-D_dat)/(B_dat+D_dat);
      avg_ch_dat_ch2[jm] = (A_dat-C_dat)/(A_dat+C_dat);
      avg_A1p_dat_ch2[ibin] = A_dat/B_dat;
      avg_A1m_dat_ch2[ibin] = C_dat/D_dat;
      avg_A3_dat_ch2[ibin] = (A_dat+C_dat)/(B_dat+D_dat);
   }
   for (int ibin=0; ibin<nbins2/2; ibin++)
   {
      int jp = (nbins2/2)+ibin;
      int jm = (nbins2/2)-ibin-1;
      double A_nuc = xsec_ch2p_nuc01[jm]*APb*LUMI*1e-6; double B_nuc = xsec_ch2p_nuc01[jp]*APb*LUMI*1e-6;
      double C_nuc = xsec_ch2m_nuc01[jm]*APb*LUMI*1e-6; double D_nuc = xsec_ch2m_nuc01[jp]*APb*LUMI*1e-6;
      double A_nonuc = xsec_ch2p_nonuc01[jm]*APb*LUMI*1e-6; double B_nonuc = xsec_ch2p_nonuc01[jp]*APb*LUMI*1e-6;
      double C_nonuc = xsec_ch2m_nonuc01[jm]*APb*LUMI*1e-6; double D_nonuc = xsec_ch2m_nonuc01[jp]*APb*LUMI*1e-6;
      double A_dat = yields2p[jm]/(etamax2[jm]-etamin2[jm]); double B_dat = yields2p[jp]/(etamax2[jp]-etamin2[jp]);
      double C_dat = yields2m[jm]/(etamax2[jm]-etamin2[jm]); double D_dat = yields2m[jp]/(etamax2[jp]-etamin2[jp]);
      for (int ieff=0; ieff<eff2p.size(); ieff++)
      {
         A_dat /= eff2p[ieff][jm];
         B_dat /= eff2p[ieff][jp];
         C_dat /= eff2m[ieff][jm];
         D_dat /= eff2m[ieff][jp];
      }
      avg_A4_nuc_ch2[ibin] = (A_nuc+C_nuc-B_nuc-D_nuc)/(total(nbins2,avg_yp_nuc_ch2)+total(nbins2,avg_ym_nuc_ch2));
      avg_A4_nonuc_ch2[ibin] = (A_nonuc+C_nonuc-B_nonuc-D_nonuc)/(total(nbins2,avg_yp_nonuc_ch2)+total(nbins2,avg_ym_nonuc_ch2));
      avg_A4_dat_ch2[ibin] = (A_dat+C_dat-B_dat-D_dat)/(total(nbins2,avg_yp_dat_ch2)+total(nbins2,avg_ym_dat_ch2));
   }

   cat(nbins1,avg_yp_nuc_ch1,nbins2,avg_yp_nuc_ch2,avg_yp_nuc_ch12);
   cat(nbins1,avg_ym_nuc_ch1,nbins2,avg_ym_nuc_ch2,avg_ym_nuc_ch12);
   cat(nbins1,avg_ch_nuc_ch1,nbins2,avg_ch_nuc_ch2,avg_ch_nuc_ch12);
   cat(2*nbins1,avg_yt_nuc_ch1,2*nbins2,avg_yt_nuc_ch2,avg_yt_nuc_ch12);
   cat(nbins1/2,avg_A1p_nuc_ch1,nbins2/2,avg_A1p_nuc_ch2,avg_A1p_nuc_ch12);
   cat(nbins1/2,avg_A1m_nuc_ch1,nbins2/2,avg_A1m_nuc_ch2,avg_A1m_nuc_ch12);
   cat(nbins1/2,avg_A3_nuc_ch1,nbins2/2,avg_A3_nuc_ch2,avg_A3_nuc_ch12);
   cat(nbins1/2,avg_A4_nuc_ch1,nbins2/2,avg_A4_nuc_ch2,avg_A4_nuc_ch12);
   cat(nbins1,avg_yp_nonuc_ch1,nbins2,avg_yp_nonuc_ch2,avg_yp_nonuc_ch12);
   cat(nbins1,avg_ym_nonuc_ch1,nbins2,avg_ym_nonuc_ch2,avg_ym_nonuc_ch12);
   cat(nbins1,avg_ch_nonuc_ch1,nbins2,avg_ch_nonuc_ch2,avg_ch_nonuc_ch12);
   cat(2*nbins1,avg_yt_nonuc_ch1,2*nbins2,avg_yt_nonuc_ch2,avg_yt_nonuc_ch12);
   cat(nbins1/2,avg_A1p_nonuc_ch1,nbins2/2,avg_A1p_nonuc_ch2,avg_A1p_nonuc_ch12);
   cat(nbins1/2,avg_A1m_nonuc_ch1,nbins2/2,avg_A1m_nonuc_ch2,avg_A1m_nonuc_ch12);
   cat(nbins1/2,avg_A3_nonuc_ch1,nbins2/2,avg_A3_nonuc_ch2,avg_A3_nonuc_ch12);
   cat(nbins1/2,avg_A4_nonuc_ch1,nbins2/2,avg_A4_nonuc_ch2,avg_A4_nonuc_ch12);
   cat(nbins1,avg_yp_dat_ch1,nbins2,avg_yp_dat_ch2,avg_yp_dat_ch12);
   cat(nbins1,avg_ym_dat_ch1,nbins2,avg_ym_dat_ch2,avg_ym_dat_ch12);
   cat(nbins1,avg_ch_dat_ch1,nbins2,avg_ch_dat_ch2,avg_ch_dat_ch12);
   cat(2*nbins1,avg_yt_dat_ch1,2*nbins2,avg_yt_dat_ch2,avg_yt_dat_ch12);
   cat(nbins1/2,avg_A1p_dat_ch1,nbins2/2,avg_A1p_dat_ch2,avg_A1p_dat_ch12);
   cat(nbins1/2,avg_A1m_dat_ch1,nbins2/2,avg_A1m_dat_ch2,avg_A1m_dat_ch12);
   cat(nbins1/2,avg_A3_dat_ch1,nbins2/2,avg_A3_dat_ch2,avg_A3_dat_ch12);
   cat(nbins1/2,avg_A4_dat_ch1,nbins2/2,avg_A4_dat_ch2,avg_A4_dat_ch12);

   // print out central predictions from theory
   cout << "central prediction W+ (w/ nuclear effects, channel 1)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
      cout << xsec_ch1p_nuc01[ibin]*APb*LUMI*1e-6*(etamax1[ibin]-etamin1[ibin]) << " ";
   cout << endl;
   cout << "central prediction W- (w/ nuclear effects, channel 1)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
      cout << xsec_ch1m_nuc01[ibin]*APb*LUMI*1e-6*(etamax1[ibin]-etamin1[ibin]) << " ";
   cout << endl;
   cout << "central prediction charge asymmetry (w/ nuclear effects, channel 1)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
      cout << avg_ch_nuc_ch1[ibin] << " ";
   cout << endl;
   cout << "central prediction A1- asymmetry (w/ nuclear effects, channel 1)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
      cout << avg_A1m_nuc_ch1[ibin] << " ";
   cout << endl;
   cout << "central prediction A1+ asymmetry (w/ nuclear effects, channel 1)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
      cout << avg_A1p_nuc_ch1[ibin] << " ";
   cout << endl;
   cout << "central prediction A3 asymmetry (w/ nuclear effects, channel 1)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
      cout << avg_A3_nuc_ch1[ibin] << " ";
   cout << endl;
   cout << "central prediction A4 asymmetry (w/ nuclear effects, channel 1)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
      cout << avg_A4_nuc_ch1[ibin] << " ";
   cout << endl;
   cout << "central prediction W+ (w/o nuclear effects, channel 1)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
      cout << xsec_ch1p_nonuc01[ibin]*APb*LUMI*1e-6*(etamax1[ibin]-etamin1[ibin]) << " ";
   cout << endl;
   cout << "central prediction W- (w/o nuclear effects, channel 1)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
      cout << xsec_ch1m_nonuc01[ibin]*APb*LUMI*1e-6*(etamax1[ibin]-etamin1[ibin]) << " ";
   cout << endl;
   cout << "central prediction charge asymmetry (w/o nuclear effects, channel 1)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
      cout << avg_ch_nonuc_ch1[ibin] << " ";
   cout << endl;
   cout << "central prediction A1- asymmetry (w/o nuclear effects, channel 1)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
      cout << avg_A1m_nonuc_ch1[ibin] << " ";
   cout << endl;
   cout << "central prediction A1+ asymmetry (w/o nuclear effects, channel 1)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
      cout << avg_A1p_nonuc_ch1[ibin] << " ";
   cout << endl;
   cout << "central prediction A3 asymmetry (w/o nuclear effects, channel 1)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
      cout << avg_A3_nonuc_ch1[ibin] << " ";
   cout << endl;
   cout << "central prediction A4 asymmetry (w/o nuclear effects, channel 1)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
      cout << avg_A4_nonuc_ch1[ibin] << " ";
   cout << endl;
   cout << "central prediction W+ (w/ nuclear effects, channel 2)" << endl;
   for (int ibin=0; ibin<nbins2; ibin++)
      cout << xsec_ch2p_nuc01[ibin]*APb*LUMI*1e-6*(etamax2[ibin]-etamin2[ibin]) << " ";
   cout << endl;
   cout << "central prediction W- (w/ nuclear effects, channel 2)" << endl;
   for (int ibin=0; ibin<nbins2; ibin++)
      cout << xsec_ch2m_nuc01[ibin]*APb*LUMI*1e-6*(etamax2[ibin]-etamin2[ibin]) << " ";
   cout << endl;
   cout << "central prediction charge asymmetry (w/ nuclear effects, channel 2)" << endl;
   for (int ibin=0; ibin<nbins2; ibin++)
      cout << avg_ch_nuc_ch2[ibin] << " ";
   cout << endl;
   cout << "central prediction A1- asymmetry (w/ nuclear effects, channel 2)" << endl;
   for (int ibin=0; ibin<nbins2/2; ibin++)
      cout << avg_A1m_nuc_ch2[ibin] << " ";
   cout << endl;
   cout << "central prediction A1+ asymmetry (w/ nuclear effects, channel 2)" << endl;
   for (int ibin=0; ibin<nbins2/2; ibin++)
      cout << avg_A1p_nuc_ch2[ibin] << " ";
   cout << endl;
   cout << "central prediction A3 asymmetry (w/ nuclear effects, channel 2)" << endl;
   for (int ibin=0; ibin<nbins2/2; ibin++)
      cout << avg_A3_nuc_ch2[ibin] << " ";
   cout << endl;
   cout << "central prediction A4 asymmetry (w/ nuclear effects, channel 2)" << endl;
   for (int ibin=0; ibin<nbins2/2; ibin++)
      cout << avg_A4_nuc_ch2[ibin] << " ";
   cout << endl;
   cout << "central prediction W+ (w/o nuclear effects, channel 2)" << endl;
   for (int ibin=0; ibin<nbins2; ibin++)
      cout << xsec_ch2p_nonuc01[ibin]*APb*LUMI*1e-6*(etamax2[ibin]-etamin2[ibin]) << " ";
   cout << endl;
   cout << "central prediction W- (w/o nuclear effects, channel 2)" << endl;
   for (int ibin=0; ibin<nbins2; ibin++)
      cout << xsec_ch2m_nonuc01[ibin]*APb*LUMI*1e-6*(etamax2[ibin]-etamin2[ibin]) << " ";
   cout << endl;
   cout << "central prediction charge asymmetry (w/o nuclear effects, channel 2)" << endl;
   for (int ibin=0; ibin<nbins2; ibin++)
      cout << avg_ch_nonuc_ch2[ibin] << " ";
   cout << endl;
   cout << "central prediction A1- asymmetry (w/o nuclear effects, channel 2)" << endl;
   for (int ibin=0; ibin<nbins2/2; ibin++)
      cout << avg_A1m_nonuc_ch2[ibin] << " ";
   cout << endl;
   cout << "central prediction A1+ asymmetry (w/o nuclear effects, channel 2)" << endl;
   for (int ibin=0; ibin<nbins2/2; ibin++)
      cout << avg_A1p_nonuc_ch2[ibin] << " ";
   cout << endl;
   cout << "central prediction A3 asymmetry (w/o nuclear effects, channel 2)" << endl;
   for (int ibin=0; ibin<nbins2/2; ibin++)
      cout << avg_A3_nonuc_ch2[ibin] << " ";
   cout << endl;
   cout << "central prediction A4 asymmetry (w/o nuclear effects, channel 2)" << endl;
   for (int ibin=0; ibin<nbins2/2; ibin++)
      cout << avg_A4_nonuc_ch2[ibin] << " ";
   cout << endl;


   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // declaration of everything going to the ttree, and creation of branches
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////

   TTree *tr = new TTree("PEtree","pseudo-experiments");
   double mur, muf;
   double yp_penuc_ch1[nbins1], ym_penuc_ch1[nbins1], yt_penuc_ch1[2*nbins1];
   double A1p_penuc_ch1[nbins1/2], A1m_penuc_ch1[nbins1/2], A3_penuc_ch1[nbins1/2], A4_penuc_ch1[nbins1/2];
   double Ch_penuc_ch1[nbins1];
   double yp_penonuc_ch1[nbins1], ym_penonuc_ch1[nbins1], yt_penonuc_ch1[2*nbins1];
   double A1p_penonuc_ch1[nbins1/2], A1m_penonuc_ch1[nbins1/2], A3_penonuc_ch1[nbins1/2], A4_penonuc_ch1[nbins1/2];
   double Ch_penonuc_ch1[nbins1];
   double yp_pedat_ch1[nbins1], ym_pedat_ch1[nbins1], yt_pedat_ch1[2*nbins1];
   double A1p_pedat_ch1[nbins1/2], A1m_pedat_ch1[nbins1/2], A3_pedat_ch1[nbins1/2], A4_pedat_ch1[nbins1/2];
   double Ch_pedat_ch1[nbins1];
   double yp_penuc_ch2[nbins2], ym_penuc_ch2[nbins2], yt_penuc_ch2[2*nbins2];
   double A1p_penuc_ch2[nbins2/2], A1m_penuc_ch2[nbins2/2], A3_penuc_ch2[nbins2/2], A4_penuc_ch2[nbins2/2];
   double Ch_penuc_ch2[nbins2];
   double yp_penonuc_ch2[nbins2], ym_penonuc_ch2[nbins2], yt_penonuc_ch2[2*nbins2];
   double A1p_penonuc_ch2[nbins2/2], A1m_penonuc_ch2[nbins2/2], A3_penonuc_ch2[nbins2/2], A4_penonuc_ch2[nbins2/2];
   double Ch_penonuc_ch2[nbins2];
   double yp_pedat_ch2[nbins2], ym_pedat_ch2[nbins2], yt_pedat_ch2[2*nbins2];
   double A1p_pedat_ch2[nbins2/2], A1m_pedat_ch2[nbins2/2], A3_pedat_ch2[nbins2/2], A4_pedat_ch2[nbins2/2];
   double Ch_pedat_ch2[nbins2];
   double yp_penuc_ch12[(nbins1+nbins2)], ym_penuc_ch12[(nbins1+nbins2)], yt_penuc_ch12[2*(nbins1+nbins2)];
   double A1p_penuc_ch12[(nbins1+nbins2)/2], A1m_penuc_ch12[(nbins1+nbins2)/2], A3_penuc_ch12[(nbins1+nbins2)/2], A4_penuc_ch12[(nbins1+nbins2)/2];
   double Ch_penuc_ch12[(nbins1+nbins2)];
   double yp_penonuc_ch12[(nbins1+nbins2)], ym_penonuc_ch12[(nbins1+nbins2)], yt_penonuc_ch12[2*(nbins1+nbins2)];
   double A1p_penonuc_ch12[(nbins1+nbins2)/2], A1m_penonuc_ch12[(nbins1+nbins2)/2], A3_penonuc_ch12[(nbins1+nbins2)/2], A4_penonuc_ch12[(nbins1+nbins2)/2];
   double Ch_penonuc_ch12[(nbins1+nbins2)];
   double yp_pedat_ch12[(nbins1+nbins2)], ym_pedat_ch12[(nbins1+nbins2)], yt_pedat_ch12[2*(nbins1+nbins2)];
   double A1p_pedat_ch12[(nbins1+nbins2)/2], A1m_pedat_ch12[(nbins1+nbins2)/2], A3_pedat_ch12[(nbins1+nbins2)/2], A4_pedat_ch12[(nbins1+nbins2)/2];
   double Ch_pedat_ch12[(nbins1+nbins2)];

   tr->Branch("mur",&mur,"mur/D");
   tr->Branch("muf",&muf,"muf/D");
   tr->Branch("yp_penuc_ch12",yp_penuc_ch12,Form("yp_penuc_ch12[%i]/D",(nbins1+nbins2)));
   tr->Branch("ym_penuc_ch12",ym_penuc_ch12,Form("ym_penuc_ch12[%i]/D",(nbins1+nbins2)));
   tr->Branch("yt_penuc_ch12",yt_penuc_ch12,Form("yt_penuc_ch12[%i]/D",2*(nbins1+nbins2)));
   tr->Branch("A1p_penuc_ch12",A1p_penuc_ch12,Form("A1p_penuc_ch12[%i]/D",(nbins1+nbins2)/2));
   tr->Branch("A1m_penuc_ch12",A1m_penuc_ch12,Form("A1m_penuc_ch12[%i]/D",(nbins1+nbins2)/2));
   tr->Branch("A3_penuc_ch12",A3_penuc_ch12,Form("A3_penuc_ch12[%i]/D",(nbins1+nbins2)/2));
   tr->Branch("A4_penuc_ch12",A4_penuc_ch12,Form("A4_penuc_ch12[%i]/D",(nbins1+nbins2)/2));
   tr->Branch("Ch_penuc_ch12",Ch_penuc_ch12,Form("Ch_penuc_ch12[%i]/D",(nbins1+nbins2)));
   tr->Branch("yp_penonuc_ch12",yp_penonuc_ch12,Form("yp_penonuc_ch12[%i]/D",(nbins1+nbins2)));
   tr->Branch("ym_penonuc_ch12",ym_penonuc_ch12,Form("ym_penonuc_ch12[%i]/D",(nbins1+nbins2)));
   tr->Branch("yt_penonuc_ch12",yt_penonuc_ch12,Form("yt_penonuc_ch12[%i]/D",2*(nbins1+nbins2)));
   tr->Branch("A1p_penonuc_ch12",A1p_penonuc_ch12,Form("A1p_penonuc_ch12[%i]/D",(nbins1+nbins2)/2));
   tr->Branch("A1m_penonuc_ch12",A1m_penonuc_ch12,Form("A1m_penonuc_ch12[%i]/D",(nbins1+nbins2)/2));
   tr->Branch("A3_penonuc_ch12",A3_penonuc_ch12,Form("A3_penonuc_ch12[%i]/D",(nbins1+nbins2)/2));
   tr->Branch("A4_penonuc_ch12",A4_penonuc_ch12,Form("A4_penonuc_ch12[%i]/D",(nbins1+nbins2)/2));
   tr->Branch("Ch_penonuc_ch12",Ch_penonuc_ch12,Form("Ch_penonuc_ch12[%i]/D",(nbins1+nbins2)));
   tr->Branch("yp_pedat_ch12",yp_pedat_ch12,Form("yp_pedat_ch12[%i]/D",(nbins1+nbins2)));
   tr->Branch("ym_pedat_ch12",ym_pedat_ch12,Form("ym_pedat_ch12[%i]/D",(nbins1+nbins2)));
   tr->Branch("yt_pedat_ch12",yt_pedat_ch12,Form("yt_pedat_ch12[%i]/D",2*(nbins1+nbins2)));
   tr->Branch("A1p_pedat_ch12",A1p_pedat_ch12,Form("A1p_pedat_ch12[%i]/D",(nbins1+nbins2)/2));
   tr->Branch("A1m_pedat_ch12",A1m_pedat_ch12,Form("A1m_pedat_ch12[%i]/D",(nbins1+nbins2)/2));
   tr->Branch("A3_pedat_ch12",A3_pedat_ch12,Form("A3_pedat_ch12[%i]/D",(nbins1+nbins2)/2));
   tr->Branch("A4_pedat_ch12",A4_pedat_ch12,Form("A4_pedat_ch12[%i]/D",(nbins1+nbins2)/2));
   tr->Branch("Ch_pedat_ch12",Ch_pedat_ch12,Form("Ch_pedat_ch12[%i]/D",(nbins1+nbins2)));

   // chi2's
   double chi2_yp_penuc_thnuc;
   double chi2_ym_penuc_thnuc;
   double chi2_Ch_penuc_thnuc;
   double chi2_A1p_penuc_thnuc;
   double chi2_A1m_penuc_thnuc;
   double chi2_A3_penuc_thnuc;
   double chi2_A4_penuc_thnuc;
   double chi2_yp_penuc_thnonuc;
   double chi2_ym_penuc_thnonuc;
   double chi2_Ch_penuc_thnonuc;
   double chi2_A1p_penuc_thnonuc;
   double chi2_A1m_penuc_thnonuc;
   double chi2_A3_penuc_thnonuc;
   double chi2_A4_penuc_thnonuc;
   double chi2_yp_penonuc_thnuc;
   double chi2_ym_penonuc_thnuc;
   double chi2_Ch_penonuc_thnuc;
   double chi2_A1p_penonuc_thnuc;
   double chi2_A1m_penonuc_thnuc;
   double chi2_A3_penonuc_thnuc;
   double chi2_A4_penonuc_thnuc;
   double chi2_yp_penonuc_thnonuc;
   double chi2_ym_penonuc_thnonuc;
   double chi2_Ch_penonuc_thnonuc;
   double chi2_A1p_penonuc_thnonuc;
   double chi2_A1m_penonuc_thnonuc;
   double chi2_A3_penonuc_thnonuc;
   double chi2_A4_penonuc_thnonuc;
   tr->Branch("chi2_yp_penuc_thnuc",&chi2_yp_penuc_thnuc,"chi2_yp_penuc_thnuc/D");
   tr->Branch("chi2_ym_penuc_thnuc",&chi2_ym_penuc_thnuc,"chi2_ym_penuc_thnuc/D");
   tr->Branch("chi2_Ch_penuc_thnuc",&chi2_Ch_penuc_thnuc,"chi2_Ch_penuc_thnuc/D");
   tr->Branch("chi2_A1p_penuc_thnuc",&chi2_A1p_penuc_thnuc,"chi2_A1p_penuc_thnuc/D");
   tr->Branch("chi2_A1m_penuc_thnuc",&chi2_A1m_penuc_thnuc,"chi2_A1m_penuc_thnuc/D");
   tr->Branch("chi2_A3_penuc_thnuc",&chi2_A3_penuc_thnuc,"chi2_A3_penuc_thnuc/D");
   tr->Branch("chi2_A4_penuc_thnuc",&chi2_A4_penuc_thnuc,"chi2_A4_penuc_thnuc/D");
   tr->Branch("chi2_yp_penuc_thnonuc",&chi2_yp_penuc_thnonuc,"chi2_yp_penuc_thnonuc/D");
   tr->Branch("chi2_ym_penuc_thnonuc",&chi2_ym_penuc_thnonuc,"chi2_ym_penuc_thnonuc/D");
   tr->Branch("chi2_Ch_penuc_thnonuc",&chi2_Ch_penuc_thnonuc,"chi2_Ch_penuc_thnonuc/D");
   tr->Branch("chi2_A1p_penuc_thnonuc",&chi2_A1p_penuc_thnonuc,"chi2_A1p_penuc_thnonuc/D");
   tr->Branch("chi2_A1m_penuc_thnonuc",&chi2_A1m_penuc_thnonuc,"chi2_A1m_penuc_thnonuc/D");
   tr->Branch("chi2_A3_penuc_thnonuc",&chi2_A3_penuc_thnonuc,"chi2_A3_penuc_thnonuc/D");
   tr->Branch("chi2_A4_penuc_thnonuc",&chi2_A4_penuc_thnonuc,"chi2_A4_penuc_thnonuc/D");
   tr->Branch("chi2_yp_penonuc_thnuc",&chi2_yp_penonuc_thnuc,"chi2_yp_penonuc_thnuc/D");
   tr->Branch("chi2_ym_penonuc_thnuc",&chi2_ym_penonuc_thnuc,"chi2_ym_penonuc_thnuc/D");
   tr->Branch("chi2_Ch_penonuc_thnuc",&chi2_Ch_penonuc_thnuc,"chi2_Ch_penonuc_thnuc/D");
   tr->Branch("chi2_A1p_penonuc_thnuc",&chi2_A1p_penonuc_thnuc,"chi2_A1p_penonuc_thnuc/D");
   tr->Branch("chi2_A1m_penonuc_thnuc",&chi2_A1m_penonuc_thnuc,"chi2_A1m_penonuc_thnuc/D");
   tr->Branch("chi2_A3_penonuc_thnuc",&chi2_A3_penonuc_thnuc,"chi2_A3_penonuc_thnuc/D");
   tr->Branch("chi2_A4_penonuc_thnuc",&chi2_A4_penonuc_thnuc,"chi2_A4_penonuc_thnuc/D");
   tr->Branch("chi2_yp_penonuc_thnonuc",&chi2_yp_penonuc_thnonuc,"chi2_yp_penonuc_thnonuc/D");
   tr->Branch("chi2_ym_penonuc_thnonuc",&chi2_ym_penonuc_thnonuc,"chi2_ym_penonuc_thnonuc/D");
   tr->Branch("chi2_Ch_penonuc_thnonuc",&chi2_Ch_penonuc_thnonuc,"chi2_Ch_penonuc_thnonuc/D");
   tr->Branch("chi2_A1p_penonuc_thnonuc",&chi2_A1p_penonuc_thnonuc,"chi2_A1p_penonuc_thnonuc/D");
   tr->Branch("chi2_A1m_penonuc_thnonuc",&chi2_A1m_penonuc_thnonuc,"chi2_A1m_penonuc_thnonuc/D");
   tr->Branch("chi2_A3_penonuc_thnonuc",&chi2_A3_penonuc_thnonuc,"chi2_A3_penonuc_thnonuc/D");
   tr->Branch("chi2_A4_penonuc_thnonuc",&chi2_A4_penonuc_thnonuc,"chi2_A4_penonuc_thnonuc/D");

   // get covariance matrices
   // w/ and w/o nuc effects, all observables, ch1 only since we only want the combo
   TFile *fth = new TFile("pes_th.root");
   TH2F *hcov2_yp_nuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_yp_nuc_ch1");
   TH2F *hcov2_yp_nonuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_yp_nonuc_ch1");
   TH2F *hcov2_ym_nuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_ym_nuc_ch1");
   TH2F *hcov2_ym_nonuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_ym_nonuc_ch1");
   TH2F *hcov2_yt_nuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_yt_nuc_ch1");
   TH2F *hcov2_yt_nonuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_yt_nonuc_ch1");
   TH2F *hcov2_ch_nuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_ch_nuc_ch1");
   TH2F *hcov2_ch_nonuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_ch_nonuc_ch1");
   TH2F *hcov2_A1p_nuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_A1p_nuc_ch1");
   TH2F *hcov2_A1p_nonuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_A1p_nonuc_ch1");
   TH2F *hcov2_A1m_nuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_A1m_nuc_ch1");
   TH2F *hcov2_A1m_nonuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_A1m_nonuc_ch1");
   TH2F *hcov2_A3_nuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_A3_nuc_ch1");
   TH2F *hcov2_A3_nonuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_A3_nonuc_ch1");
   TH2F *hcov2_A4_nuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_A4_nuc_ch1");
   TH2F *hcov2_A4_nonuc_ch1 = (TH2F*) fth->Get("matrices/hcov2_A4_nonuc_ch1");

   // we also need the covariance matrices for experimental uncertainties (stat+syst+lumi)
   TFile *fexp = new TFile("pes_exp_muacc.root");
   fexp->cd("matrices");
   TH2F *hcov2_yp_dat_ch12 = (TH2F*) fexp->Get("matrices/hcov2_yp_dat_ch12");
   TH2F *hcov2_ym_dat_ch12 = (TH2F*) fexp->Get("matrices/hcov2_ym_dat_ch12");
   TH2F *hcov2_yt_dat_ch12 = (TH2F*) fexp->Get("matrices/hcov2_yt_dat_ch12");
   TH2F *hcov2_ch_dat_ch12 = (TH2F*) fexp->Get("matrices/hcov2_ch_dat_ch12");
   TH2F *hcov2_A1p_dat_ch12 = (TH2F*) fexp->Get("matrices/hcov2_A1p_dat_ch12");
   TH2F *hcov2_A1m_dat_ch12 = (TH2F*) fexp->Get("matrices/hcov2_A1m_dat_ch12");
   TH2F *hcov2_A3_dat_ch12 = (TH2F*) fexp->Get("matrices/hcov2_A3_dat_ch12");
   TH2F *hcov2_A4_dat_ch12 = (TH2F*) fexp->Get("matrices/hcov2_A4_dat_ch12");
   

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // pseudo-experiment loop
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   for (int ipe=0; ipe<NPE; ipe++)
   {
      if (ipe%(NPE/10)==0) cout << ipe << "/" << NPE << endl;
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      // 2. throw pulls:
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      // - nbins pulls per eff / SF (or 2*nbins if +/- are independent)
      //  - trig
      //  - id+iso
      //  - reco
      //  - single mu eff
      //  - qcd
      // - 1 pull for lumi
      // - 26 pulls for CT10 pdf
      // - 15 pulls for eps pdf
      // - 2 pulls for scale (if pull<0 then divide by 1+1/-pull instead of multiplying by 1+pull. retry if scalemod(pull1, pull2)=0)
      vector< vector<double> > pull1_eff_plus;
      vector< vector<double> > pull1_eff_minus;
      vector< vector<double> > pull1_eff_all;
      vector< vector<double> > pull2_eff_plus;
      vector< vector<double> > pull2_eff_minus;
      vector< vector<double> > pull2_eff_all;
      double pull_lumi;
      double *pull_ct10 = new double[NSET_CT10];
      double *pull_eps09 = new double[NSET_EPS09];
      double pull_scale[2];
      bool dosysts = true;//(ipe>0);

      // pulls for efficiencies: channel 1
      for (int ieff=0; ieff<eff1p.size(); ieff++)
      {
         vector<double> tmpvec, tmpvec2;
         if (efferr1p[ieff][0]>0) 
            // in this case no covariance between muon charges: one pull for each charge
         {
            for (int ibin=0; ibin<nbins1; ibin++)
            {
               double thepull = !(ibin>0 && eff1p[ieff][ibin]==eff1p[ieff][ibin-1] && efferr1p[ieff][ibin]==efferr1p[ieff][ibin-1]) ? gRandom->Gaus(0,dosysts*DOEFFS) : tmpvec.back();
               double thepull2 = !(ibin>0 && eff1m[ieff][ibin]==eff1m[ieff][ibin-1] && efferr1m[ieff][ibin]==efferr1m[ieff][ibin-1]) ? gRandom->Gaus(0,dosysts*DOEFFS) : tmpvec2.back();
               tmpvec.push_back(thepull);
               tmpvec2.push_back(thepull2);
            }
            pull1_eff_all.push_back(tmpvec);
            pull1_eff_plus.push_back(tmpvec);
            pull1_eff_minus.push_back(tmpvec2);
         }
         else // in this case only one pull for both plus and minus charges
         {
            for (int ibin=0; ibin<nbins1; ibin++)
            {
               double thepull = !(ibin>0 && eff1p[ieff][ibin]==eff1p[ieff][ibin-1] && efferr1p[ieff][ibin]==efferr1p[ieff][ibin-1]) ? gRandom->Gaus(0,dosysts*DOEFFS) : tmpvec.back();
               tmpvec.push_back(thepull);
            }
            pull1_eff_all.push_back(tmpvec);
            pull1_eff_plus.push_back(tmpvec);
            pull1_eff_minus.push_back(tmpvec);
         }
      }

      // pulls for efficiencies: channel 2
      for (int ieff=0; ieff<eff2p.size(); ieff++)
      {
         vector<double> tmpvec, tmpvec2;
         if (efferr2p[ieff][0]>0) 
            // in this case no covariance between muon charges: one pull for each charge
         {
            for (int ibin=0; ibin<nbins2; ibin++)
            {
               double thepull = !(ibin>0 && eff2p[ieff][ibin]==eff2p[ieff][ibin-1] && efferr2p[ieff][ibin]==efferr2p[ieff][ibin-1]) ? gRandom->Gaus(0,dosysts*DOEFFS) : tmpvec.back();
               double thepull2 = !(ibin>0 && eff2m[ieff][ibin]==eff2m[ieff][ibin-1] && efferr2m[ieff][ibin]==efferr2m[ieff][ibin-1]) ? gRandom->Gaus(0,dosysts*DOEFFS) : tmpvec2.back();
               // check also if the efficiency is the same in the other channel
               for (int ieff1=0; ieff1<eff1p.size(); ieff1++)
                  if (eff2p[ieff][ibin]==eff1p[ieff1][ibin] && efferr2p[ieff][ibin]==efferr1p[ieff1][ibin]) thepull = pull1_eff_plus[ieff1][ibin];
               for (int ieff1=0; ieff1<eff1m.size(); ieff1++)
                  if (eff2m[ieff][ibin]==eff1m[ieff1][ibin] && efferr2m[ieff][ibin]==efferr1m[ieff1][ibin]) thepull2 = pull1_eff_minus[ieff1][ibin];
               tmpvec.push_back(thepull);
               tmpvec2.push_back(thepull2);
            }
            pull2_eff_all.push_back(tmpvec);
            pull2_eff_plus.push_back(tmpvec);
            pull2_eff_minus.push_back(tmpvec2);
         }
         else // in this case only one pull for both plus and minus charges
         {
            for (int ibin=0; ibin<nbins2; ibin++)
            {
               double thepull = !(ibin>0 && eff2p[ieff][ibin]==eff2p[ieff][ibin-1] && efferr2p[ieff][ibin]==efferr2p[ieff][ibin-1]) ? gRandom->Gaus(0,dosysts*DOEFFS) : tmpvec.back();
               // check also if the efficiency is the same in the other channel
               for (int ieff1=0; ieff1<eff1p.size(); ieff1++)
                  if (eff2p[ieff][ibin]==eff1p[ieff1][ibin] && efferr2p[ieff][ibin]==efferr1p[ieff1][ibin]) thepull = pull1_eff_all[ieff1][ibin];
               tmpvec.push_back(thepull);
            }
            pull2_eff_all.push_back(tmpvec);
            pull2_eff_plus.push_back(tmpvec);
            pull2_eff_minus.push_back(tmpvec);
         }
      }

      // pull for lumi
      pull_lumi = gRandom->Gaus(0,dosysts*DOLUMI);

      // pulls for scale. Be careful to check that we are in the allowed range for scales.
      double tmp[2] = {0};
      bool ok=false;
      double dummyxsecs[7] = {1,1,1,1,1,1,1};
      while (!ok)
      {
         tmp[0] = pow(2,gRandom->Gaus(0,dosysts*DOSCALE));
         tmp[1] = pow(2,gRandom->Gaus(0,dosysts*DOSCALE));
         ok = (scalemod(dummyxsecs,tmp[0],tmp[1])!=0);
      }
      for (int iscale=0; iscale<2; iscale++)
         pull_scale[iscale] = tmp[iscale];
      mur = tmp[0];
      muf = tmp[1];

      // pulls for ct10
      for (int iset=0; iset<NSET_CT10; iset++)
         pull_ct10[iset] = gRandom->Gaus(0,dosysts*DOCT10);

      // pulls for EPS09
      for (int iset=0; iset<NSET_EPS09; iset++)
         pull_eps09[iset] = gRandom->Gaus(0,dosysts*DOEPS09);




      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      // 3. compute new model
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      // - 1 function for the scale variation
      // - 1 function for effect of pdf uncert on CT10
      // - 1 function for effect of pdf uncert on CT10+EPS09
      // - multiply modified theory by product (1+sk) where sk are the pulls for the different effs
      double *smeared_nuc_ch1plus = new double[nbins1];
      double *smeared_nuc_ch1minus = new double[nbins1];
      double *smeared_nonuc_ch1plus = new double[nbins1];
      double *smeared_nonuc_ch1minus = new double[nbins1];
      double *smeared_nuc_ch2plus = new double[nbins2];
      double *smeared_nuc_ch2minus = new double[nbins2];
      double *smeared_nonuc_ch2plus = new double[nbins2];
      double *smeared_nonuc_ch2minus = new double[nbins2];

      for (int ibin=0; ibin<nbins1; ibin++)
      {
         smeared_nuc_ch1plus[ibin] = xsec_ch1p_nuc01[ibin]*APb*LUMI*(1+pull_lumi*LUMIERR)*1e-6*(etamax1[ibin]-etamin1[ibin]) * scalemod(xsec_ch1p_scale[ibin],pull_scale[0],pull_scale[1]);
         smeared_nuc_ch1plus[ibin] *= pdfmod(NSET_CT10,xsec_ch1p_nuc01[ibin],xsec_ch1p_nucp1_ct10[ibin],xsec_ch1p_nucm1_ct10[ibin],pull_ct10);
         smeared_nuc_ch1plus[ibin] *= pdfmod(NSET_EPS09,xsec_ch1p_nuc01[ibin],xsec_ch1p_nucp1_eps09[ibin],xsec_ch1p_nucm1_eps09[ibin],pull_eps09);
         smeared_nonuc_ch1plus[ibin] = xsec_ch1p_nonuc01[ibin]*APb*LUMI*(1+pull_lumi*LUMIERR)*1e-6*(etamax1[ibin]-etamin1[ibin]) * scalemod(xsec_ch1p_scale[ibin],pull_scale[0],pull_scale[1]);
         smeared_nonuc_ch1plus[ibin] *= pdfmod(NSET_CT10,xsec_ch1p_nonuc01[ibin],xsec_ch1p_nonucp1_ct10[ibin],xsec_ch1p_nonucm1_ct10[ibin],pull_ct10);
         smeared_nuc_ch1minus[ibin] = xsec_ch1m_nuc01[ibin]*APb*LUMI*(1+pull_lumi*LUMIERR)*1e-6*(etamax1[ibin]-etamin1[ibin]) * scalemod(xsec_ch1m_scale[ibin],pull_scale[0],pull_scale[1]);
         smeared_nuc_ch1minus[ibin] *= pdfmod(NSET_CT10,xsec_ch1m_nuc01[ibin],xsec_ch1m_nucp1_ct10[ibin],xsec_ch1m_nucm1_ct10[ibin],pull_ct10);
         smeared_nuc_ch1minus[ibin] *= pdfmod(NSET_EPS09,xsec_ch1m_nuc01[ibin],xsec_ch1m_nucp1_eps09[ibin],xsec_ch1m_nucm1_eps09[ibin],pull_eps09);
         smeared_nonuc_ch1minus[ibin] = xsec_ch1m_nonuc01[ibin]*APb*LUMI*(1+pull_lumi*LUMIERR)*1e-6*(etamax1[ibin]-etamin1[ibin]) * scalemod(xsec_ch1m_scale[ibin],pull_scale[0],pull_scale[1]);
         smeared_nonuc_ch1minus[ibin] *= pdfmod(NSET_CT10,xsec_ch1m_nonuc01[ibin],xsec_ch1m_nonucp1_ct10[ibin],xsec_ch1m_nonucm1_ct10[ibin],pull_ct10);
      }
      for (int ibin=0; ibin<nbins2; ibin++)
      {
         smeared_nuc_ch2plus[ibin] = xsec_ch2p_nuc01[ibin]*APb*LUMI*(1+pull_lumi*LUMIERR)*1e-6*(etamax2[ibin]-etamin2[ibin]) * scalemod(xsec_ch2p_scale[ibin],pull_scale[0],pull_scale[1]);
         smeared_nuc_ch2plus[ibin] *= pdfmod(NSET_CT10,xsec_ch2p_nuc01[ibin],xsec_ch2p_nucp1_ct10[ibin],xsec_ch2p_nucm1_ct10[ibin],pull_ct10);
         smeared_nuc_ch2plus[ibin] *= pdfmod(NSET_EPS09,xsec_ch2p_nuc01[ibin],xsec_ch2p_nucp1_eps09[ibin],xsec_ch2p_nucm1_eps09[ibin],pull_eps09);
         smeared_nonuc_ch2plus[ibin] = xsec_ch2p_nonuc01[ibin]*APb*LUMI*(1+pull_lumi*LUMIERR)*1e-6*(etamax2[ibin]-etamin2[ibin]) * scalemod(xsec_ch2p_scale[ibin],pull_scale[0],pull_scale[1]);
         smeared_nonuc_ch2plus[ibin] *= pdfmod(NSET_CT10,xsec_ch2p_nonuc01[ibin],xsec_ch2p_nonucp1_ct10[ibin],xsec_ch2p_nonucm1_ct10[ibin],pull_ct10);
         smeared_nuc_ch2minus[ibin] = xsec_ch2m_nuc01[ibin]*APb*LUMI*(1+pull_lumi*LUMIERR)*1e-6*(etamax2[ibin]-etamin2[ibin]) * scalemod(xsec_ch2m_scale[ibin],pull_scale[0],pull_scale[1]);
         smeared_nuc_ch2minus[ibin] *= pdfmod(NSET_CT10,xsec_ch2m_nuc01[ibin],xsec_ch2m_nucp1_ct10[ibin],xsec_ch2m_nucm1_ct10[ibin],pull_ct10);
         smeared_nuc_ch2minus[ibin] *= pdfmod(NSET_EPS09,xsec_ch2m_nuc01[ibin],xsec_ch2m_nucp1_eps09[ibin],xsec_ch2m_nucm1_eps09[ibin],pull_eps09);
         smeared_nonuc_ch2minus[ibin] = xsec_ch2m_nonuc01[ibin]*APb*LUMI*(1+pull_lumi*LUMIERR)*1e-6*(etamax2[ibin]-etamin2[ibin]) * scalemod(xsec_ch2m_scale[ibin],pull_scale[0],pull_scale[1]);
         smeared_nonuc_ch2minus[ibin] *= pdfmod(NSET_CT10,xsec_ch2m_nonuc01[ibin],xsec_ch2m_nonucp1_ct10[ibin],xsec_ch2m_nonucm1_ct10[ibin],pull_ct10);
      }



      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      // 4. throw gaussian distributed new yields from new model. Width = nominal stat uncert * sqrt(modified model / nominal model)
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      double *pe_nuc_ch1plus = new double[nbins1];
      double *pe_nuc_ch1minus = new double[nbins1];
      double *pe_nonuc_ch1plus = new double[nbins1];
      double *pe_nonuc_ch1minus = new double[nbins1];
      double *pe_dat_ch1plus = new double[nbins1];
      double *pe_dat_ch1minus = new double[nbins1];
      double *pe_nuc_ch2plus = new double[nbins2];
      double *pe_nuc_ch2minus = new double[nbins2];
      double *pe_nonuc_ch2plus = new double[nbins2];
      double *pe_nonuc_ch2minus = new double[nbins2];
      double *pe_dat_ch2plus = new double[nbins2];
      double *pe_dat_ch2minus = new double[nbins2];

      for (int ibin=0; ibin<nbins1; ibin++)
      {
         double raw_nuc_ch1plus = smeared_nuc_ch1plus[ibin];
         double raw_nonuc_ch1plus = smeared_nonuc_ch1plus[ibin];
         double raw_dat_ch1plus = yields1p[ibin]/(1+pull_lumi*LUMIERR);
         double raw_nuc_ch1minus = smeared_nuc_ch1minus[ibin];
         double raw_nonuc_ch1minus = smeared_nonuc_ch1minus[ibin];
         double raw_dat_ch1minus = yields1m[ibin]/(1+pull_lumi*LUMIERR);

         for (int ieff=0; ieff<eff1p.size(); ieff++)
         {
            raw_nuc_ch1plus *= eff1p[ieff][ibin];
            raw_nonuc_ch1plus *= eff1p[ieff][ibin];
            raw_nuc_ch1minus *= eff1m[ieff][ibin];
            raw_nonuc_ch1minus *= eff1m[ieff][ibin];
         }

         pe_nuc_ch1plus[ibin] = gRandom->Gaus(raw_nuc_ch1plus,DOSTAT*staterr1p[ibin]*sqrt(raw_nuc_ch1plus/yields1p[ibin]));
         pe_nonuc_ch1plus[ibin] = gRandom->Gaus(raw_nonuc_ch1plus,DOSTAT*staterr1p[ibin]*sqrt(raw_nonuc_ch1plus/yields1p[ibin]));
         pe_dat_ch1plus[ibin] = gRandom->Gaus(raw_dat_ch1plus,DOSTAT*staterr1p[ibin]*sqrt(raw_dat_ch1plus/yields1p[ibin]));
         pe_nuc_ch1minus[ibin] = gRandom->Gaus(raw_nuc_ch1minus,DOSTAT*staterr1m[ibin]*sqrt(raw_nuc_ch1minus/yields1m[ibin]));
         pe_nonuc_ch1minus[ibin] = gRandom->Gaus(raw_nonuc_ch1minus,DOSTAT*staterr1m[ibin]*sqrt(raw_nonuc_ch1minus/yields1m[ibin]));
         pe_dat_ch1minus[ibin] = gRandom->Gaus(raw_dat_ch1minus,DOSTAT*staterr1m[ibin]*sqrt(raw_dat_ch1minus/yields1m[ibin]));

         for (int ieff=0; ieff<eff1p.size(); ieff++)
            if (efferr1p[ieff][ibin]>0)
            {
               pe_nuc_ch1plus[ibin] /= (eff1p[ieff][ibin] + efferr1p[ieff][ibin]*pull1_eff_plus[ieff][ibin]);
               pe_nonuc_ch1plus[ibin] /= (eff1p[ieff][ibin] + efferr1p[ieff][ibin]*pull1_eff_plus[ieff][ibin]);
               pe_dat_ch1plus[ibin] /= (eff1p[ieff][ibin] + efferr1p[ieff][ibin]*pull1_eff_plus[ieff][ibin]);
               pe_nuc_ch1minus[ibin] /= (eff1m[ieff][ibin] + efferr1m[ieff][ibin]*pull1_eff_minus[ieff][ibin]);
               pe_nonuc_ch1minus[ibin] /= (eff1m[ieff][ibin] + efferr1m[ieff][ibin]*pull1_eff_minus[ieff][ibin]);
               pe_dat_ch1minus[ibin] /= (eff1m[ieff][ibin] + efferr1m[ieff][ibin]*pull1_eff_minus[ieff][ibin]);
            }
            else
            {
               pe_nuc_ch1plus[ibin] /= (eff1p[ieff][ibin] + efferr1p[ieff][ibin]*pull1_eff_all[ieff][ibin]);
               pe_nonuc_ch1plus[ibin] /= (eff1p[ieff][ibin] + efferr1p[ieff][ibin]*pull1_eff_all[ieff][ibin]);
               pe_dat_ch1plus[ibin] /= (eff1p[ieff][ibin] + efferr1p[ieff][ibin]*pull1_eff_all[ieff][ibin]);
               pe_nuc_ch1minus[ibin] /= (eff1m[ieff][ibin] + efferr1m[ieff][ibin]*pull1_eff_all[ieff][ibin]);
               pe_nonuc_ch1minus[ibin] /= (eff1m[ieff][ibin] + efferr1m[ieff][ibin]*pull1_eff_all[ieff][ibin]);
               pe_dat_ch1minus[ibin] /= (eff1m[ieff][ibin] + efferr1m[ieff][ibin]*pull1_eff_all[ieff][ibin]);
            }
         pe_nuc_ch1plus[ibin] /= (etamax1[ibin]-etamin1[ibin]);
         pe_nonuc_ch1plus[ibin] /= (etamax1[ibin]-etamin1[ibin]);
         pe_dat_ch1plus[ibin] /= (etamax1[ibin]-etamin1[ibin]);
         pe_nuc_ch1minus[ibin] /= (etamax1[ibin]-etamin1[ibin]);
         pe_nonuc_ch1minus[ibin] /= (etamax1[ibin]-etamin1[ibin]);
         pe_dat_ch1minus[ibin] /= (etamax1[ibin]-etamin1[ibin]);
      }

      for (int ibin=0; ibin<nbins2; ibin++)
      {
         double raw_nuc_ch2plus = smeared_nuc_ch2plus[ibin];
         double raw_nonuc_ch2plus = smeared_nonuc_ch2plus[ibin];
         double raw_dat_ch2plus = yields2p[ibin]/(1+pull_lumi*LUMIERR);
         double raw_nuc_ch2minus = smeared_nuc_ch2minus[ibin];
         double raw_nonuc_ch2minus = smeared_nonuc_ch2minus[ibin];
         double raw_dat_ch2minus = yields2m[ibin]/(1+pull_lumi*LUMIERR);
         for (int ieff=0; ieff<eff2p.size(); ieff++)
         {
            raw_dat_ch2plus /= eff2p[ieff][ibin];
            raw_dat_ch2minus /= eff2m[ieff][ibin];
         }

         for (int ieff=0; ieff<eff2p.size(); ieff++)
            if (efferr2p[ieff][ibin]>0)
            {
               raw_nuc_ch2plus *= (eff2p[ieff][ibin] + efferr2p[ieff][ibin]*pull2_eff_plus[ieff][ibin]);
               raw_nonuc_ch2plus *= (eff2p[ieff][ibin] + efferr2p[ieff][ibin]*pull2_eff_plus[ieff][ibin]);
               raw_dat_ch2plus *= (eff2p[ieff][ibin] + efferr2p[ieff][ibin]*pull2_eff_plus[ieff][ibin]);
               raw_nuc_ch2minus *= (eff2m[ieff][ibin] + efferr2m[ieff][ibin]*pull2_eff_minus[ieff][ibin]);
               raw_nonuc_ch2minus *= (eff2m[ieff][ibin] + efferr2m[ieff][ibin]*pull2_eff_minus[ieff][ibin]);
               raw_dat_ch2minus *= (eff2m[ieff][ibin] + efferr2m[ieff][ibin]*pull2_eff_minus[ieff][ibin]);
            }
            else
            {
               raw_nuc_ch2plus *= (eff2p[ieff][ibin] + efferr2p[ieff][ibin]*pull2_eff_all[ieff][ibin]);
               raw_nonuc_ch2plus *= (eff2p[ieff][ibin] + efferr2p[ieff][ibin]*pull2_eff_all[ieff][ibin]);
               raw_dat_ch2plus *= (eff2p[ieff][ibin] + efferr2p[ieff][ibin]*pull2_eff_all[ieff][ibin]);
               raw_nuc_ch2minus *= (eff2m[ieff][ibin] + efferr2m[ieff][ibin]*pull2_eff_all[ieff][ibin]);
               raw_nonuc_ch2minus *= (eff2m[ieff][ibin] + efferr2m[ieff][ibin]*pull2_eff_all[ieff][ibin]);
               raw_dat_ch2minus *= (eff2m[ieff][ibin] + efferr2m[ieff][ibin]*pull2_eff_all[ieff][ibin]);
            }

         pe_nuc_ch2plus[ibin] = gRandom->Gaus(raw_nuc_ch2plus,DOSTAT*staterr2p[ibin]*sqrt(raw_nuc_ch2plus/yields2p[ibin]));
         pe_nonuc_ch2plus[ibin] = gRandom->Gaus(raw_nonuc_ch2plus,DOSTAT*staterr2p[ibin]*sqrt(raw_nonuc_ch2plus/yields2p[ibin]));
         pe_dat_ch2plus[ibin] = gRandom->Gaus(raw_dat_ch2plus,DOSTAT*staterr2p[ibin]*sqrt(raw_dat_ch2plus/yields2p[ibin]));
         pe_nuc_ch2minus[ibin] = gRandom->Gaus(raw_nuc_ch2minus,DOSTAT*staterr2m[ibin]*sqrt(raw_nuc_ch2minus/yields2m[ibin]));
         pe_nonuc_ch2minus[ibin] = gRandom->Gaus(raw_nonuc_ch2minus,DOSTAT*staterr2m[ibin]*sqrt(raw_nonuc_ch2minus/yields2m[ibin]));
         pe_dat_ch2minus[ibin] = gRandom->Gaus(raw_dat_ch2minus,DOSTAT*staterr2m[ibin]*sqrt(raw_dat_ch2minus/yields2m[ibin]));

         for (int ieff=0; ieff<eff2p.size(); ieff++)
         {
            pe_nuc_ch2plus[ibin] /= eff2p[ieff][ibin];
            pe_nonuc_ch2plus[ibin] /= eff2p[ieff][ibin];
            pe_dat_ch2plus[ibin] /= eff2p[ieff][ibin];
            pe_nuc_ch2minus[ibin] /= eff2m[ieff][ibin];
            pe_nonuc_ch2minus[ibin] /= eff2m[ieff][ibin];
            pe_dat_ch2minus[ibin] /= eff2m[ieff][ibin];
         }
         pe_nuc_ch2plus[ibin] /= (etamax2[ibin]-etamin2[ibin]);
         pe_nonuc_ch2plus[ibin] /= (etamax2[ibin]-etamin2[ibin]);
         pe_dat_ch2plus[ibin] /= (etamax2[ibin]-etamin2[ibin]);
         pe_nuc_ch2minus[ibin] /= (etamax2[ibin]-etamin2[ibin]);
         pe_nonuc_ch2minus[ibin] /= (etamax2[ibin]-etamin2[ibin]);
         pe_dat_ch2minus[ibin] /= (etamax2[ibin]-etamin2[ibin]);
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      // 5. compute derived quantities (asymetries...) for both pseudo-data and models
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      for (int ibin=0; ibin<nbins1/2; ibin++)
      {
         int jp = (nbins1/2)+ibin;
         int jm = (nbins1/2)-ibin-1;

         double A_penuc_ch1 = pe_nuc_ch1plus[jm]; double B_penuc_ch1 = pe_nuc_ch1plus[jp];
         double C_penuc_ch1 = pe_nuc_ch1minus[jm]; double D_penuc_ch1 = pe_nuc_ch1minus[jp];
         double A_penonuc_ch1 = pe_nonuc_ch1plus[jm]; double B_penonuc_ch1 = pe_nonuc_ch1plus[jp];
         double C_penonuc_ch1 = pe_nonuc_ch1minus[jm]; double D_penonuc_ch1 = pe_nonuc_ch1minus[jp];
         double A_pedat_ch1 = pe_dat_ch1plus[jm]; double B_pedat_ch1 = pe_dat_ch1plus[jp];
         double C_pedat_ch1 = pe_dat_ch1minus[jm]; double D_pedat_ch1 = pe_dat_ch1minus[jp];

         yp_penuc_ch1[jp] = pe_nuc_ch1plus[jp];
         yp_penuc_ch1[jm] = pe_nuc_ch1plus[jm];
         ym_penuc_ch1[jp] = pe_nuc_ch1minus[jp];
         ym_penuc_ch1[jm] = pe_nuc_ch1minus[jm];
         yt_penuc_ch1[jp] = pe_nuc_ch1plus[jp];
         yt_penuc_ch1[jm] = pe_nuc_ch1plus[jm];
         yt_penuc_ch1[nbins1+jp] = pe_nuc_ch1minus[jp];
         yt_penuc_ch1[nbins1+jm] = pe_nuc_ch1minus[jm];
         Ch_penuc_ch1[jp] = (B_penuc_ch1-D_penuc_ch1)/(B_penuc_ch1+D_penuc_ch1);
         Ch_penuc_ch1[jm] = (A_penuc_ch1-C_penuc_ch1)/(A_penuc_ch1+C_penuc_ch1);
         A1p_penuc_ch1[ibin] = A_penuc_ch1/B_penuc_ch1;
         A1m_penuc_ch1[ibin] = C_penuc_ch1/D_penuc_ch1;
         A3_penuc_ch1[ibin] = (A_penuc_ch1+C_penuc_ch1)/(B_penuc_ch1+D_penuc_ch1);

         yp_penonuc_ch1[jp] = pe_nonuc_ch1plus[jp];
         yp_penonuc_ch1[jm] = pe_nonuc_ch1plus[jm];
         ym_penonuc_ch1[jp] = pe_nonuc_ch1minus[jp];
         ym_penonuc_ch1[jm] = pe_nonuc_ch1minus[jm];
         yt_penonuc_ch1[jp] = pe_nonuc_ch1plus[jp];
         yt_penonuc_ch1[jm] = pe_nonuc_ch1plus[jm];
         yt_penonuc_ch1[nbins1+jp] = pe_nonuc_ch1minus[jp];
         yt_penonuc_ch1[nbins1+jm] = pe_nonuc_ch1minus[jm];
         Ch_penonuc_ch1[jp] = (B_penonuc_ch1-D_penonuc_ch1)/(B_penonuc_ch1+D_penonuc_ch1);
         Ch_penonuc_ch1[jm] = (A_penonuc_ch1-C_penonuc_ch1)/(A_penonuc_ch1+C_penonuc_ch1);
         A1p_penonuc_ch1[ibin] = A_penonuc_ch1/B_penonuc_ch1;
         A1m_penonuc_ch1[ibin] = C_penonuc_ch1/D_penonuc_ch1;
         A3_penonuc_ch1[ibin] = (A_penonuc_ch1+C_penonuc_ch1)/(B_penonuc_ch1+D_penonuc_ch1);

         yp_pedat_ch1[jp] = pe_dat_ch1plus[jp];
         yp_pedat_ch1[jm] = pe_dat_ch1plus[jm];
         ym_pedat_ch1[jp] = pe_dat_ch1minus[jp];
         ym_pedat_ch1[jm] = pe_dat_ch1minus[jm];
         yt_pedat_ch1[jp] = pe_dat_ch1plus[jp];
         yt_pedat_ch1[jm] = pe_dat_ch1plus[jm];
         yt_pedat_ch1[nbins1+jp] = pe_dat_ch1minus[jp];
         yt_pedat_ch1[nbins1+jm] = pe_dat_ch1minus[jm];
         Ch_pedat_ch1[jp] = (B_pedat_ch1-D_pedat_ch1)/(B_pedat_ch1+D_pedat_ch1);
         Ch_pedat_ch1[jm] = (A_pedat_ch1-C_pedat_ch1)/(A_pedat_ch1+C_pedat_ch1);
         A1p_pedat_ch1[ibin] = A_pedat_ch1/B_pedat_ch1;
         A1m_pedat_ch1[ibin] = C_pedat_ch1/D_pedat_ch1;
         A3_pedat_ch1[ibin] = (A_pedat_ch1+C_pedat_ch1)/(B_pedat_ch1+D_pedat_ch1);
      }
      for (int ibin=0; ibin<nbins1/2; ibin++)
      {
         int jp = (nbins1/2)+ibin;
         int jm = (nbins1/2)-ibin-1;
         double A_penuc_ch1 = pe_nuc_ch1plus[jm]; double B_penuc_ch1 = pe_nuc_ch1plus[jp];
         double C_penuc_ch1 = pe_nuc_ch1minus[jm]; double D_penuc_ch1 = pe_nuc_ch1minus[jp];
         double A_penonuc_ch1 = pe_nonuc_ch1plus[jm]; double B_penonuc_ch1 = pe_nonuc_ch1plus[jp];
         double C_penonuc_ch1 = pe_nonuc_ch1minus[jm]; double D_penonuc_ch1 = pe_nonuc_ch1minus[jp];
         double A_pedat_ch1 = pe_dat_ch1plus[jm]; double B_pedat_ch1 = pe_dat_ch1plus[jp];
         double C_pedat_ch1 = pe_dat_ch1minus[jm]; double D_pedat_ch1 = pe_dat_ch1minus[jp];
         A4_penuc_ch1[ibin] = (A_penuc_ch1+C_penuc_ch1-B_penuc_ch1-D_penuc_ch1)/(total(nbins1,pe_nuc_ch1plus)+total(nbins1,pe_nuc_ch1minus));
         A4_penonuc_ch1[ibin] = (A_penonuc_ch1+C_penonuc_ch1-B_penonuc_ch1-D_penonuc_ch1)/(total(nbins1,pe_nonuc_ch1plus)+total(nbins1,pe_nonuc_ch1minus));
         A4_pedat_ch1[ibin] = (A_pedat_ch1+C_pedat_ch1-B_pedat_ch1-D_pedat_ch1)/(total(nbins1,pe_dat_ch1plus)+total(nbins1,pe_dat_ch1minus));
      }
      for (int ibin=0; ibin<nbins2/2; ibin++)
      {
         int jp = (nbins2/2)+ibin;
         int jm = (nbins2/2)-ibin-1;

         double A_penuc_ch2 = pe_nuc_ch2plus[jm]; double B_penuc_ch2 = pe_nuc_ch2plus[jp];
         double C_penuc_ch2 = pe_nuc_ch2minus[jm]; double D_penuc_ch2 = pe_nuc_ch2minus[jp];
         double A_penonuc_ch2 = pe_nonuc_ch2plus[jm]; double B_penonuc_ch2 = pe_nonuc_ch2plus[jp];
         double C_penonuc_ch2 = pe_nonuc_ch2minus[jm]; double D_penonuc_ch2 = pe_nonuc_ch2minus[jp];
         double A_pedat_ch2 = pe_dat_ch2plus[jm]; double B_pedat_ch2 = pe_dat_ch2plus[jp];
         double C_pedat_ch2 = pe_dat_ch2minus[jm]; double D_pedat_ch2 = pe_dat_ch2minus[jp];

         yp_penuc_ch2[jp] = pe_nuc_ch2plus[jp];
         yp_penuc_ch2[jm] = pe_nuc_ch2plus[jm];
         ym_penuc_ch2[jp] = pe_nuc_ch2minus[jp];
         ym_penuc_ch2[jm] = pe_nuc_ch2minus[jm];
         yt_penuc_ch2[jp] = pe_nuc_ch2plus[jp];
         yt_penuc_ch2[jm] = pe_nuc_ch2plus[jm];
         yt_penuc_ch2[nbins2+jp] = pe_nuc_ch2minus[jp];
         yt_penuc_ch2[nbins2+jm] = pe_nuc_ch2minus[jm];
         Ch_penuc_ch2[jp] = (B_penuc_ch2-D_penuc_ch2)/(B_penuc_ch2+D_penuc_ch2);
         Ch_penuc_ch2[jm] = (A_penuc_ch2-C_penuc_ch2)/(A_penuc_ch2+C_penuc_ch2);
         A1p_penuc_ch2[ibin] = A_penuc_ch2/B_penuc_ch2;
         A1m_penuc_ch2[ibin] = C_penuc_ch2/D_penuc_ch2;
         A3_penuc_ch2[ibin] = (A_penuc_ch2+C_penuc_ch2)/(B_penuc_ch2+D_penuc_ch2);

         yp_penonuc_ch2[jp] = pe_nonuc_ch2plus[jp];
         yp_penonuc_ch2[jm] = pe_nonuc_ch2plus[jm];
         ym_penonuc_ch2[jp] = pe_nonuc_ch2minus[jp];
         ym_penonuc_ch2[jm] = pe_nonuc_ch2minus[jm];
         yt_penonuc_ch2[jp] = pe_nonuc_ch2plus[jp];
         yt_penonuc_ch2[jm] = pe_nonuc_ch2plus[jm];
         yt_penonuc_ch2[nbins2+jp] = pe_nonuc_ch2minus[jp];
         yt_penonuc_ch2[nbins2+jm] = pe_nonuc_ch2minus[jm];
         Ch_penonuc_ch2[jp] = (B_penonuc_ch2-D_penonuc_ch2)/(B_penonuc_ch2+D_penonuc_ch2);
         Ch_penonuc_ch2[jm] = (A_penonuc_ch2-C_penonuc_ch2)/(A_penonuc_ch2+C_penonuc_ch2);
         A1p_penonuc_ch2[ibin] = A_penonuc_ch2/B_penonuc_ch2;
         A1m_penonuc_ch2[ibin] = C_penonuc_ch2/D_penonuc_ch2;
         A3_penonuc_ch2[ibin] = (A_penonuc_ch2+C_penonuc_ch2)/(B_penonuc_ch2+D_penonuc_ch2);

         yp_pedat_ch2[jp] = pe_dat_ch2plus[jp];
         yp_pedat_ch2[jm] = pe_dat_ch2plus[jm];
         ym_pedat_ch2[jp] = pe_dat_ch2minus[jp];
         ym_pedat_ch2[jm] = pe_dat_ch2minus[jm];
         yt_pedat_ch2[jp] = pe_dat_ch2plus[jp];
         yt_pedat_ch2[jm] = pe_dat_ch2plus[jm];
         yt_pedat_ch2[nbins2+jp] = pe_dat_ch2minus[jp];
         yt_pedat_ch2[nbins2+jm] = pe_dat_ch2minus[jm];
         Ch_pedat_ch2[jp] = (B_pedat_ch2-D_pedat_ch2)/(B_pedat_ch2+D_pedat_ch2);
         Ch_pedat_ch2[jm] = (A_pedat_ch2-C_pedat_ch2)/(A_pedat_ch2+C_pedat_ch2);
         A1p_pedat_ch2[ibin] = A_pedat_ch2/B_pedat_ch2;
         A1m_pedat_ch2[ibin] = C_pedat_ch2/D_pedat_ch2;
         A3_pedat_ch2[ibin] = (A_pedat_ch2+C_pedat_ch2)/(B_pedat_ch2+D_pedat_ch2);
      }
      for (int ibin=0; ibin<nbins2/2; ibin++)
      {
         int jp = (nbins2/2)+ibin;
         int jm = (nbins2/2)-ibin-1;
         double A_penuc_ch2 = pe_nuc_ch2plus[jm]; double B_penuc_ch2 = pe_nuc_ch2plus[jp];
         double C_penuc_ch2 = pe_nuc_ch2minus[jm]; double D_penuc_ch2 = pe_nuc_ch2minus[jp];
         double A_penonuc_ch2 = pe_nonuc_ch2plus[jm]; double B_penonuc_ch2 = pe_nonuc_ch2plus[jp];
         double C_penonuc_ch2 = pe_nonuc_ch2minus[jm]; double D_penonuc_ch2 = pe_nonuc_ch2minus[jp];
         double A_pedat_ch2 = pe_dat_ch2plus[jm]; double B_pedat_ch2 = pe_dat_ch2plus[jp];
         double C_pedat_ch2 = pe_dat_ch2minus[jm]; double D_pedat_ch2 = pe_dat_ch2minus[jp];
         A4_penuc_ch2[ibin] = (A_penuc_ch2+C_penuc_ch2-B_penuc_ch2-D_penuc_ch2)/(total(nbins2,pe_nuc_ch2plus)+total(nbins2,pe_nuc_ch2minus));
         A4_penonuc_ch2[ibin] = (A_penonuc_ch2+C_penonuc_ch2-B_penonuc_ch2-D_penonuc_ch2)/(total(nbins2,pe_nonuc_ch2plus)+total(nbins2,pe_nonuc_ch2minus));
         A4_pedat_ch2[ibin] = (A_pedat_ch2+C_pedat_ch2-B_pedat_ch2-D_pedat_ch2)/(total(nbins2,pe_dat_ch2plus)+total(nbins2,pe_dat_ch2minus));
      }

      cat(nbins1,yp_penuc_ch1,nbins2,yp_penuc_ch2,yp_penuc_ch12);
      cat(nbins1,ym_penuc_ch1,nbins2,ym_penuc_ch2,ym_penuc_ch12);
      cat(nbins1,Ch_penuc_ch1,nbins2,Ch_penuc_ch2,Ch_penuc_ch12);
      cat(2*nbins1,yt_penuc_ch1,2*nbins2,yt_penuc_ch2,yt_penuc_ch12);
      cat(nbins1/2,A1p_penuc_ch1,nbins2/2,A1p_penuc_ch2,A1p_penuc_ch12);
      cat(nbins1/2,A1m_penuc_ch1,nbins2/2,A1m_penuc_ch2,A1m_penuc_ch12);
      cat(nbins1/2,A3_penuc_ch1,nbins2/2,A3_penuc_ch2,A3_penuc_ch12);
      cat(nbins1/2,A4_penuc_ch1,nbins2/2,A4_penuc_ch2,A4_penuc_ch12);
      cat(nbins1,yp_penonuc_ch1,nbins2,yp_penonuc_ch2,yp_penonuc_ch12);
      cat(nbins1,ym_penonuc_ch1,nbins2,ym_penonuc_ch2,ym_penonuc_ch12);
      cat(nbins1,Ch_penonuc_ch1,nbins2,Ch_penonuc_ch2,Ch_penonuc_ch12);
      cat(2*nbins1,yt_penonuc_ch1,2*nbins2,yt_penonuc_ch2,yt_penonuc_ch12);
      cat(nbins1/2,A1p_penonuc_ch1,nbins2/2,A1p_penonuc_ch2,A1p_penonuc_ch12);
      cat(nbins1/2,A1m_penonuc_ch1,nbins2/2,A1m_penonuc_ch2,A1m_penonuc_ch12);
      cat(nbins1/2,A3_penonuc_ch1,nbins2/2,A3_penonuc_ch2,A3_penonuc_ch12);
      cat(nbins1/2,A4_penonuc_ch1,nbins2/2,A4_penonuc_ch2,A4_penonuc_ch12);
      cat(nbins1,yp_pedat_ch1,nbins2,yp_pedat_ch2,yp_pedat_ch12);
      cat(nbins1,ym_pedat_ch1,nbins2,ym_pedat_ch2,ym_pedat_ch12);
      cat(nbins1,Ch_pedat_ch1,nbins2,Ch_pedat_ch2,Ch_pedat_ch12);
      cat(2*nbins1,yt_pedat_ch1,2*nbins2,yt_pedat_ch2,yt_pedat_ch12);
      cat(nbins1/2,A1p_pedat_ch1,nbins2/2,A1p_pedat_ch2,A1p_pedat_ch12);
      cat(nbins1/2,A1m_pedat_ch1,nbins2/2,A1m_pedat_ch2,A1m_pedat_ch12);
      cat(nbins1/2,A3_pedat_ch1,nbins2/2,A3_pedat_ch2,A3_pedat_ch12);
      cat(nbins1/2,A4_pedat_ch1,nbins2/2,A4_pedat_ch2,A4_pedat_ch12);

      // compute chi2's
      chi2_yp_penuc_thnuc = combine_and_chi2(nbins1+nbins2,yp_penuc_ch12,"yieldsWp","yieldsWp_muon.dta",hcov2_yp_dat_ch12,hcov2_yp_nuc_ch1,true);
      chi2_ym_penuc_thnuc = combine_and_chi2(nbins1+nbins2,ym_penuc_ch12,"yieldsWm","yieldsWm_muon.dta",hcov2_ym_dat_ch12,hcov2_ym_nuc_ch1,true);
      chi2_Ch_penuc_thnuc = combine_and_chi2(nbins1+nbins2,Ch_penuc_ch12,"CA","CA_muon.dta",hcov2_ch_dat_ch12,hcov2_ch_nuc_ch1,true);
      chi2_A1p_penuc_thnuc = combine_and_chi2((nbins1+nbins2)/2,A1p_penuc_ch12,"A1p","A1p_muon.dta",hcov2_A1p_dat_ch12,hcov2_A1p_nuc_ch1,true);
      chi2_A1m_penuc_thnuc = combine_and_chi2((nbins1+nbins2)/2,A1m_penuc_ch12,"A1m","A1m_muon.dta",hcov2_A1m_dat_ch12,hcov2_A1m_nuc_ch1,true);
      chi2_A3_penuc_thnuc = combine_and_chi2((nbins1+nbins2)/2,A3_penuc_ch12,"A3","A3_muon.dta",hcov2_A3_dat_ch12,hcov2_A3_nuc_ch1,true);
      chi2_A4_penuc_thnuc = combine_and_chi2((nbins1+nbins2)/2,A4_penuc_ch12,"A4","A4_muon.dta",hcov2_A4_dat_ch12,hcov2_A4_nuc_ch1,true);
      chi2_yp_penuc_thnonuc = combine_and_chi2(nbins1+nbins2,yp_penuc_ch12,"yieldsWp","yieldsWp_muon.dta",hcov2_yp_dat_ch12,hcov2_yp_nonuc_ch1,false);
      chi2_ym_penuc_thnonuc = combine_and_chi2(nbins1+nbins2,ym_penuc_ch12,"yieldsWm","yieldsWm_muon.dta",hcov2_ym_dat_ch12,hcov2_ym_nonuc_ch1,false);
      chi2_Ch_penuc_thnonuc = combine_and_chi2(nbins1+nbins2,Ch_penuc_ch12,"CA","CA_muon.dta",hcov2_ch_dat_ch12,hcov2_ch_nonuc_ch1,false);
      chi2_A1p_penuc_thnonuc = combine_and_chi2((nbins1+nbins2)/2,A1p_penuc_ch12,"A1p","A1p_muon.dta",hcov2_A1p_dat_ch12,hcov2_A1p_nonuc_ch1,false);
      chi2_A1m_penuc_thnonuc = combine_and_chi2((nbins1+nbins2)/2,A1m_penuc_ch12,"A1m","A1m_muon.dta",hcov2_A1m_dat_ch12,hcov2_A1m_nonuc_ch1,false);
      chi2_A3_penuc_thnonuc = combine_and_chi2((nbins1+nbins2)/2,A3_penuc_ch12,"A3","A3_muon.dta",hcov2_A3_dat_ch12,hcov2_A3_nonuc_ch1,false);
      chi2_A4_penuc_thnonuc = combine_and_chi2((nbins1+nbins2)/2,A4_penuc_ch12,"A4","A4_muon.dta",hcov2_A4_dat_ch12,hcov2_A4_nonuc_ch1,false);
      chi2_yp_penonuc_thnuc = combine_and_chi2(nbins1+nbins2,yp_penonuc_ch12,"yieldsWp","yieldsWp_muon.dta",hcov2_yp_dat_ch12,hcov2_yp_nuc_ch1,true);
      chi2_ym_penonuc_thnuc = combine_and_chi2(nbins1+nbins2,ym_penonuc_ch12,"yieldsWm","yieldsWm_muon.dta",hcov2_ym_dat_ch12,hcov2_ym_nuc_ch1,true);
      chi2_Ch_penonuc_thnuc = combine_and_chi2(nbins1+nbins2,Ch_penonuc_ch12,"CA","CA_muon.dta",hcov2_ch_dat_ch12,hcov2_ch_nuc_ch1,true);
      chi2_A1p_penonuc_thnuc = combine_and_chi2((nbins1+nbins2)/2,A1p_penonuc_ch12,"A1p","A1p_muon.dta",hcov2_A1p_dat_ch12,hcov2_A1p_nuc_ch1,true);
      chi2_A1m_penonuc_thnuc = combine_and_chi2((nbins1+nbins2)/2,A1m_penonuc_ch12,"A1m","A1m_muon.dta",hcov2_A1m_dat_ch12,hcov2_A1m_nuc_ch1,true);
      chi2_A3_penonuc_thnuc = combine_and_chi2((nbins1+nbins2)/2,A3_penonuc_ch12,"A3","A3_muon.dta",hcov2_A3_dat_ch12,hcov2_A3_nuc_ch1,true);
      chi2_A4_penonuc_thnuc = combine_and_chi2((nbins1+nbins2)/2,A4_penonuc_ch12,"A4","A4_muon.dta",hcov2_A4_dat_ch12,hcov2_A4_nuc_ch1,true);
      chi2_yp_penonuc_thnonuc = combine_and_chi2(nbins1+nbins2,yp_penonuc_ch12,"yieldsWp","yieldsWp_muon.dta",hcov2_yp_dat_ch12,hcov2_yp_nonuc_ch1,false);
      chi2_ym_penonuc_thnonuc = combine_and_chi2(nbins1+nbins2,ym_penonuc_ch12,"yieldsWm","yieldsWm_muon.dta",hcov2_ym_dat_ch12,hcov2_ym_nonuc_ch1,false);
      chi2_Ch_penonuc_thnonuc = combine_and_chi2(nbins1+nbins2,Ch_penonuc_ch12,"CA","CA_muon.dta",hcov2_ch_dat_ch12,hcov2_ch_nonuc_ch1,false);
      chi2_A1p_penonuc_thnonuc = combine_and_chi2((nbins1+nbins2)/2,A1p_penonuc_ch12,"A1p","A1p_muon.dta",hcov2_A1p_dat_ch12,hcov2_A1p_nonuc_ch1,false);
      chi2_A1m_penonuc_thnonuc = combine_and_chi2((nbins1+nbins2)/2,A1m_penonuc_ch12,"A1m","A1m_muon.dta",hcov2_A1m_dat_ch12,hcov2_A1m_nonuc_ch1,false);
      chi2_A3_penonuc_thnonuc = combine_and_chi2((nbins1+nbins2)/2,A3_penonuc_ch12,"A3","A3_muon.dta",hcov2_A3_dat_ch12,hcov2_A3_nonuc_ch1,false);
      chi2_A4_penonuc_thnonuc = combine_and_chi2((nbins1+nbins2)/2,A4_penonuc_ch12,"A4","A4_muon.dta",hcov2_A4_dat_ch12,hcov2_A4_nonuc_ch1,false);


      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      // 8. fill tree
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      tr->Fill();

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      // 9. delete stuff
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      delete pull_ct10;
      delete pull_eps09;
      delete smeared_nuc_ch1plus;
      delete smeared_nuc_ch1minus;
      delete smeared_nonuc_ch1plus;
      delete smeared_nonuc_ch1minus;
      delete pe_nuc_ch1plus;
      delete pe_nuc_ch1minus;
      delete pe_nonuc_ch1plus;
      delete pe_nonuc_ch1minus;
      delete pe_dat_ch1plus;
      delete pe_dat_ch1minus;
      delete smeared_nuc_ch2plus;
      delete smeared_nuc_ch2minus;
      delete smeared_nonuc_ch2plus;
      delete smeared_nonuc_ch2minus;
      delete pe_nuc_ch2plus;
      delete pe_nuc_ch2minus;
      delete pe_nonuc_ch2plus;
      delete pe_nonuc_ch2minus;
      delete pe_dat_ch2plus;
      delete pe_dat_ch2minus;
   }

   fout->Write();
   fout->Close();
   fexp->Close();
   fth->Close();
}

double quad(double a, double b, double c, double x)
// a = f(-1), b = f(0), c = f(1), return f(x)
{
   double p0 = b;
   double p1 = (c-a)/2.;
   double p2 = c-p0-p1;
   return p0 + p1*x + p2*x*x;
}

double chi2(int nbins, double *ni, double *dni, double *ti, vector<double> Sk)
{
   double chi2=0;
   // cout << nbins << " " << Sk.size() << endl;
   for (int i=0; i<nbins; i++)
   {
      //     double prod = 1;
      //     for (int k=0; k<Sk.size(); k++)
      //       prod *= (1+fik[i][k]*Sk[k]);
      //     double tiprime = ti[i]/prod;
      //     chi2 += pow(ni[i]-tiprime,2)/pow(dni[i],2);
      chi2 += pow(ni[i]-ti[i],2)/pow(dni[i],2);
   }

   for (int k=0; k<Sk.size(); k++)
      chi2 += pow(Sk[k],2);

   return chi2;
}

double scalemod(double xsecs[7], double muf, double mur)
{
   double x[7] = {1,0.5,1,0.5,1,2,2};
   double y[7] = {1,1,0.5,0.5,2,1,2};
   TGraph2D *gr = new TGraph2D(7,x,y,xsecs);
   double val = gr->Interpolate(muf,mur)/xsecs[0];
   delete gr;
   return val;
}

TGraph2D* scalegraph(double xsecs[7])
{
   double x[7] = {1,0.5,1,0.5,1,2,2};
   double y[7] = {1,1,0.5,0.5,2,1,2};
   double mod[7];
   for (int i=0; i<7; i++) mod[i] = xsecs[i]/xsecs[0]-1;
   TGraph2D *gr = new TGraph2D(7,x,y,mod);
   return gr;
}

double pdfmod(int nsets, double xsec0, double *xsecsp, double *xsecsm, double *pulls)
{
   double fact = 1.;
   for (int iset=0; iset<nsets; iset++)
   {
      fact *= quad(xsecsm[iset]/xsec0,1.,xsecsp[iset]/xsec0,pulls[iset]);
   }
   return fact;
}

double mysqrt(double a)
{
   return a>0 ? sqrt(a) : -sqrt(-a);
   // return a;
}

double reverse(int nbins, double *array)
{
   for (int ibin=0; ibin<nbins/2; ibin++)
   {
      double buf = array[ibin];
      array[ibin] = array[nbins-1-ibin];
      array[nbins-1-ibin] = buf;
   }
}

double reverse(int nbins1, int nbins2, double **array)
{
   for (int iset=0; iset<nbins2; iset++)
      for (int ibin=0; ibin<nbins1/2; ibin++)
      {
         double buf = array[ibin][iset];
         array[ibin][iset] = array[nbins1-1-ibin][iset];
         array[nbins1-1-ibin][iset] = buf;
      }
}

double total(int n, double *array)
{
   double tot=0;
   for (int i=0; i<n; i++) tot += array[i];
   return tot;
}

double total(int n, vector<double> array)
{
   double tot=0;
   for (int i=0; i<n; i++) tot += array[i];
   return tot;
}

double init(int n, double* array)
{
   for (int ibin=0; ibin<n; ibin++)
      array[ibin] = 0.;
}

double init(int n1, int n2, double** array)
{
   for (int ibin=0; ibin<n1; ibin++)
      for (int jbin=0; jbin<n2; jbin++)
         array[ibin][jbin] = 0.;
}
