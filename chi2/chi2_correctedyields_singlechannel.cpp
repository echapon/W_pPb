#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <math.h>

#include "TRandom3.h"
#include "TGraph2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TMatrixDSym.h"
#include "TMath.h"
#include "TH2.h"

#define NSET_CT10 26
#define NSET_EPS09 15

#define DOEFFS 0.
#define DOLUMI 0.
#define DOCT10 0.
#define DOEPS09 0.
#define DOSCALE 0.
#define DOSTAT 1.

#define APb 208

#define NPE 10000

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

int main(int argc, const char** argv)
{
   if (argc != 10)
   {
      cout << "Usage: " << argv[0] << " nbins data_plus.txt data_minus.txt th_plus_nuc.txt th_minus_nuc.txt th_plus_nonuc.txt th_minus_nonuc.txt th_scale_plus.txt th_minus_nonuc.txt" << endl;
      return 0;
   }

   // initialize random generator
   gRandom->SetSeed();

   TFile *fout = new TFile("pes.root","RECREATE");

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // 1. read in 2 (+, -) * 4 (exp, pdf uncert ct10, pdf uncert eps09, scale uncert) files
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   statchannel chan1p, chan1m;
   int nbins1;
   int nsyst1;

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

   vector<double> bins = chan1p.get_bins();
   vector<double> yields1p = chan1p.get_yields();
   vector<double> yields1m = chan1m.get_yields();
   vector<double> staterr1p = chan1p.get_staterr();
   vector<double> staterr1m = chan1m.get_staterr();
   vector< vector<double> > eff1p = chan1p.get_eff();
   vector< vector<double> > eff1m = chan1m.get_eff();
   vector< vector<double> > efferr1p = chan1p.get_efferr();
   vector< vector<double> > efferr1m = chan1m.get_efferr();

   // some declarations
   double *etamin1 = new double[nbins1];
   double *etamax1 = new double[nbins1];
   double *etaavg1 = new double[nbins1];
   double *xsec_mup_nuc01 = new double[nbins1];
   double **xsec_mup_nucp1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_mup_nucp1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_mup_nucm1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_mup_nucm1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_mup_nucp1_eps09 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_mup_nucp1_eps09[ibin] = new double[NSET_EPS09];
   double **xsec_mup_nucm1_eps09 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_mup_nucm1_eps09[ibin] = new double[NSET_EPS09];
   double *xsec_mup_nonuc01 = new double[nbins1];
   double **xsec_mup_nonucp1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_mup_nonucp1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_mup_nonucm1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_mup_nonucm1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_mup_scale = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_mup_scale[ibin] = new double[7];
   double *xsec_mum_nuc01 = new double[nbins1];
   double **xsec_mum_nucp1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_mum_nucp1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_mum_nucm1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_mum_nucm1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_mum_nucp1_eps09 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_mum_nucp1_eps09[ibin] = new double[NSET_EPS09];
   double **xsec_mum_nucm1_eps09 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_mum_nucm1_eps09[ibin] = new double[NSET_EPS09];
   double *xsec_mum_nonuc01 = new double[nbins1];
   double **xsec_mum_nonucp1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_mum_nonucp1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_mum_nonucm1_ct10 = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_mum_nonucm1_ct10[ibin] = new double[NSET_CT10];
   double **xsec_mum_scale = new double*[nbins1];
   for (int ibin=0; ibin<nbins1; ibin++)
      xsec_mum_scale[ibin] = new double[7];

   // read the files
   ifstream file_mup_nuc(argv[4]);
   ifstream file_mum_nuc(argv[5]);
   ifstream file_mup_nonuc(argv[6]);
   ifstream file_mum_nonuc(argv[7]);
   ifstream file_mup_scale(argv[8]);
   ifstream file_mum_scale(argv[9]);
   for (int ibin=0; ibin<nbins1; ibin++)
   {
      file_mup_nuc >> etamin1[ibin] >> etamax1[ibin] >> etaavg1[ibin] >> xsec_mup_nuc01[ibin];
      for (int iset=0; iset<NSET_CT10; iset++)
         file_mup_nuc >> xsec_mup_nucp1_ct10[ibin][iset] >> xsec_mup_nucm1_ct10[ibin][iset];
      for (int iset=0; iset<NSET_EPS09; iset++)
         file_mup_nuc >> xsec_mup_nucp1_eps09[ibin][iset] >> xsec_mup_nucm1_eps09[ibin][iset];

      double dummy;
      file_mup_nonuc >> dummy >> dummy >> dummy >> xsec_mup_nonuc01[ibin];
      for (int iset=0; iset<NSET_CT10; iset++)
         file_mup_nonuc >> xsec_mup_nonucp1_ct10[ibin][iset] >> xsec_mup_nonucm1_ct10[ibin][iset];

      file_mup_scale >> dummy >> dummy >> dummy;
      for (int iset=0; iset<7; iset++)
         file_mup_scale >> xsec_mup_scale[ibin][iset];

      file_mum_nuc >> dummy >> dummy >> dummy >> xsec_mum_nuc01[ibin];
      for (int iset=0; iset<NSET_CT10; iset++)
         file_mum_nuc >> xsec_mum_nucp1_ct10[ibin][iset] >> xsec_mum_nucm1_ct10[ibin][iset];
      for (int iset=0; iset<NSET_EPS09; iset++)
         file_mum_nuc >> xsec_mum_nucp1_eps09[ibin][iset] >> xsec_mum_nucm1_eps09[ibin][iset];

      file_mum_nonuc >> dummy >> dummy >> dummy >> xsec_mum_nonuc01[ibin];
      for (int iset=0; iset<NSET_CT10; iset++)
         file_mum_nonuc >> xsec_mum_nonucp1_ct10[ibin][iset] >> xsec_mum_nonucm1_ct10[ibin][iset];

      file_mum_scale >> dummy >> dummy >> dummy;
      for (int iset=0; iset<7; iset++)
         file_mum_scale >> xsec_mum_scale[ibin][iset];
   }

   // theory is not in the same eta convention as data!!! reverse it
   reverse(nbins1,xsec_mup_nuc01);
   reverse(nbins1,NSET_CT10,xsec_mup_nucp1_ct10);
   reverse(nbins1,NSET_CT10,xsec_mup_nucm1_ct10);
   reverse(nbins1,NSET_EPS09,xsec_mup_nucp1_eps09);
   reverse(nbins1,NSET_EPS09,xsec_mup_nucm1_eps09);
   reverse(nbins1,xsec_mup_nonuc01);
   reverse(nbins1,NSET_CT10,xsec_mup_nonucp1_ct10);
   reverse(nbins1,NSET_CT10,xsec_mup_nonucm1_ct10);
   reverse(nbins1,7,xsec_mup_scale);
   reverse(nbins1,xsec_mum_nuc01);
   reverse(nbins1,NSET_CT10,xsec_mum_nucp1_ct10);
   reverse(nbins1,NSET_CT10,xsec_mum_nucm1_ct10);
   reverse(nbins1,NSET_EPS09,xsec_mum_nucp1_eps09);
   reverse(nbins1,NSET_EPS09,xsec_mum_nucm1_eps09);
   reverse(nbins1,xsec_mum_nonuc01);
   reverse(nbins1,NSET_CT10,xsec_mum_nonucp1_ct10);
   reverse(nbins1,NSET_CT10,xsec_mum_nonucm1_ct10);
   reverse(nbins1,7,xsec_mum_scale);

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
   double yp_data[nbins1], ym_data[nbins1], yt_data[2*nbins1];
   double A1p_data[nbins1/2], A1m_data[nbins1/2], A3_data[nbins1/2];
   double Ch_data[nbins1];
   for (int ibin=0; ibin<nbins1/2; ibin++)
   {
      int jp = (nbins1/2)+ibin;
      int jm = (nbins1/2)-ibin-1;

      double A_data = yields1p[jm]; double B_data = yields1p[jp];
      double C_data = yields1m[jm]; double D_data = yields1m[jp];

      for (int ieff=0; ieff<eff1p.size(); ieff++)
      {
         A_data /= eff1p[ieff][jm];
         B_data /= eff1p[ieff][jp];
         C_data /= eff1m[ieff][jm];
         D_data /= eff1m[ieff][jp];
      }

      yp_data[jm] = A_data;
      yp_data[jp] = B_data;
      ym_data[jm] = C_data;
      ym_data[jp] = D_data;
      yt_data[jm] = A_data;
      yt_data[jp] = B_data;
      yt_data[jm+nbins1] = C_data;
      yt_data[jp+nbins1] = D_data;
      Ch_data[jp] = (B_data-D_data)/(B_data+D_data);
      Ch_data[jm] = (A_data-C_data)/(A_data+C_data);
      A1p_data[ibin] = A_data/B_data;
      A1m_data[ibin] = C_data/D_data;
      A3_data[ibin] = (A_data+C_data)/(B_data+D_data);
   }

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // compute mean expectation on raw yields from theory
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   double *avg_yp_nuc = new double[nbins1];
   double *avg_ym_nuc = new double[nbins1];
   double *avg_yt_nuc = new double[2*nbins1];
   double *avg_ch_nuc = new double[nbins1];
   double *avg_A1p_nuc = new double[nbins1/2];
   double *avg_A1m_nuc = new double[nbins1/2];
   double *avg_A3_nuc = new double[nbins1/2];
   double *avg_yp_nonuc = new double[nbins1];
   double *avg_ym_nonuc = new double[nbins1];
   double *avg_yt_nonuc = new double[2*nbins1];
   double *avg_ch_nonuc = new double[nbins1];
   double *avg_A1p_nonuc = new double[nbins1/2];
   double *avg_A1m_nonuc = new double[nbins1/2];
   double *avg_A3_nonuc = new double[nbins1/2];

   for (int ibin=0; ibin<nbins1/2; ibin++)
   {
      int jp = (nbins1/2)+ibin;
      int jm = (nbins1/2)-ibin-1;

      double A_nuc = xsec_mup_nuc01[jm]*APb*LUMI*1e-6*(etamax1[jm]-etamin1[jm]); double B_nuc = xsec_mup_nuc01[jp]*APb*LUMI*1e-6*(etamax1[jp]-etamin1[jp]);
      double C_nuc = xsec_mum_nuc01[jm]*APb*LUMI*1e-6*(etamax1[jm]-etamin1[jm]); double D_nuc = xsec_mum_nuc01[jp]*APb*LUMI*1e-6*(etamax1[jp]-etamin1[jp]);
      double A_nonuc = xsec_mup_nonuc01[jm]*APb*LUMI*1e-6*(etamax1[jm]-etamin1[jm]); double B_nonuc = xsec_mup_nonuc01[jp]*APb*LUMI*1e-6*(etamax1[jp]-etamin1[jp]);
      double C_nonuc = xsec_mum_nonuc01[jm]*APb*LUMI*1e-6*(etamax1[jm]-etamin1[jm]); double D_nonuc = xsec_mum_nonuc01[jp]*APb*LUMI*1e-6*(etamax1[jp]-etamin1[jp]);

      avg_yp_nuc[jm] = A_nuc;
      avg_yp_nuc[jp] = B_nuc;
      avg_ym_nuc[jm] = C_nuc;
      avg_ym_nuc[jp] = D_nuc;
      avg_yt_nuc[jm] = A_nuc;
      avg_yt_nuc[jp] = B_nuc;
      avg_yt_nuc[jm+nbins1] = C_nuc;
      avg_yt_nuc[jp+nbins1] = D_nuc;
      avg_ch_nuc[jp] = (B_nuc-D_nuc)/(B_nuc+D_nuc);
      avg_ch_nuc[jm] = (A_nuc-C_nuc)/(A_nuc+C_nuc);
      avg_A1p_nuc[ibin] = A_nuc/B_nuc;
      avg_A1m_nuc[ibin] = C_nuc/D_nuc;
      avg_A3_nuc[ibin] = (A_nuc+C_nuc)/(B_nuc+D_nuc);
      avg_yp_nonuc[jm] = A_nonuc;
      avg_yp_nonuc[jp] = B_nonuc;
      avg_ym_nonuc[jm] = C_nonuc;
      avg_ym_nonuc[jp] = D_nonuc;
      avg_yt_nonuc[jm] = A_nonuc;
      avg_yt_nonuc[jp] = B_nonuc;
      avg_yt_nonuc[jm+nbins1] = C_nonuc;
      avg_yt_nonuc[jp+nbins1] = D_nonuc;
      avg_ch_nonuc[jp] = (B_nonuc-D_nonuc)/(B_nonuc+D_nonuc);
      avg_ch_nonuc[jm] = (A_nonuc-C_nonuc)/(A_nonuc+C_nonuc);
      avg_A1p_nonuc[ibin] = A_nonuc/B_nonuc;
      avg_A1m_nonuc[ibin] = C_nonuc/D_nonuc;
      avg_A3_nonuc[ibin] = (A_nonuc+C_nonuc)/(B_nonuc+D_nonuc);
   }

   // print out central predictions from theory
   cout << "central prediction W+ (w/ nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
      cout << xsec_mup_nuc01[ibin]*APb*LUMI*1e-6*(etamax1[ibin]-etamin1[ibin]) << " ";
   cout << endl;
   cout << "central prediction W- (w/ nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
      cout << xsec_mum_nuc01[ibin]*APb*LUMI*1e-6*(etamax1[ibin]-etamin1[ibin]) << " ";
   cout << endl;
   cout << "central prediction charge asymmetry (w/ nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
      cout << avg_ch_nuc[ibin] << " ";
   cout << endl;
   cout << "central prediction A1- asymmetry (w/ nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
      cout << avg_A1m_nuc[ibin] << " ";
   cout << endl;
   cout << "central prediction A1+ asymmetry (w/ nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
      cout << avg_A1p_nuc[ibin] << " ";
   cout << endl;
   cout << "central prediction A3 asymmetry (w/ nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
      cout << avg_A3_nuc[ibin] << " ";
   cout << endl;
   cout << "central prediction W+ (w/o nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
      cout << xsec_mup_nonuc01[ibin]*APb*LUMI*1e-6*(etamax1[ibin]-etamin1[ibin]) << " ";
   cout << endl;
   cout << "central prediction W- (w/o nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
      cout << xsec_mum_nonuc01[ibin]*APb*LUMI*1e-6*(etamax1[ibin]-etamin1[ibin]) << " ";
   cout << endl;
   cout << "central prediction charge asymmetry (w/o nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
      cout << avg_ch_nonuc[ibin] << " ";
   cout << endl;
   cout << "central prediction A1- asymmetry (w/o nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
      cout << avg_A1m_nonuc[ibin] << " ";
   cout << endl;
   cout << "central prediction A1+ asymmetry (w/o nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
      cout << avg_A1p_nonuc[ibin] << " ";
   cout << endl;
   cout << "central prediction A3 asymmetry (w/o nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
      cout << avg_A3_nonuc[ibin] << " ";
   cout << endl;


   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // declaration of everything going to the ttree, and creation of branches
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////

   TTree *tr = new TTree("PEtree","pseudo-experiments");
   double mur, muf;
   double yp_penuc[nbins1], ym_penuc[nbins1], yt_penuc[2*nbins1];
   double A1p_penuc[nbins1/2], A1m_penuc[nbins1/2], A3_penuc[nbins1/2];
   double Ch_penuc[nbins1];
   double yp_penonuc[nbins1], ym_penonuc[nbins1], yt_penonuc[2*nbins1];
   double A1p_penonuc[nbins1/2], A1m_penonuc[nbins1/2], A3_penonuc[nbins1/2];
   double Ch_penonuc[nbins1];

   tr->Branch("mur",&mur,"mur/D");
   tr->Branch("muf",&muf,"muf/D");
   tr->Branch("yp_penuc",yp_penuc,Form("yp_penuc[%i]/D",nbins1));
   tr->Branch("ym_penuc",ym_penuc,Form("ym_penuc[%i]/D",nbins1));
   tr->Branch("yt_penuc",yt_penuc,Form("yt_penuc[%i]/D",2*nbins1));
   tr->Branch("A1p_penuc",A1p_penuc,Form("A1p_penuc[%i]/D",nbins1/2));
   tr->Branch("A1m_penuc",A1m_penuc,Form("A1m_penuc[%i]/D",nbins1/2));
   tr->Branch("A3_penuc",A3_penuc,Form("A3_penuc[%i]/D",nbins1/2));
   tr->Branch("Ch_penuc",Ch_penuc,Form("Ch_penuc[%i]/D",nbins1));
   tr->Branch("yp_penonuc",yp_penonuc,Form("yp_penonuc[%i]/D",nbins1));
   tr->Branch("ym_penonuc",ym_penonuc,Form("ym_penonuc[%i]/D",nbins1));
   tr->Branch("yt_penonuc",yt_penonuc,Form("yt_penonuc[%i]/D",2*nbins1));
   tr->Branch("A1p_penonuc",A1p_penonuc,Form("A1p_penonuc[%i]/D",nbins1/2));
   tr->Branch("A1m_penonuc",A1m_penonuc,Form("A1m_penonuc[%i]/D",nbins1/2));
   tr->Branch("A3_penonuc",A3_penonuc,Form("A3_penonuc[%i]/D",nbins1/2));
   tr->Branch("Ch_penonuc",Ch_penonuc,Form("Ch_penonuc[%i]/D",nbins1));
   tr->Branch("etamin1",etamin1,Form("etamin1[%i]/D",nbins1));
   tr->Branch("etamax1",etamax1,Form("etamax1[%i]/D",nbins1));
   tr->Branch("etaavg1",etaavg1,Form("etaavg1[%i]/D",nbins1));

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // covariance matrices
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   double **cm_yp_nuc = new double*[nbins1];
   double **cm_ym_nuc = new double*[nbins1];
   double **cm_yt_nuc = new double*[2*nbins1];
   double **cm_ch_nuc = new double*[nbins1];
   double **cm_A1p_nuc = new double*[nbins1/2];
   double **cm_A1m_nuc = new double*[nbins1/2];
   double **cm_A3_nuc = new double*[nbins1/2];
   double **cm_yp_nonuc = new double*[nbins1];
   double **cm_ym_nonuc = new double*[nbins1];
   double **cm_yt_nonuc = new double*[2*nbins1];
   double **cm_ch_nonuc = new double*[nbins1];
   double **cm_A1p_nonuc = new double*[nbins1/2];
   double **cm_A1m_nonuc = new double*[nbins1/2];
   double **cm_A3_nonuc = new double*[nbins1/2];
   for (int ibin=0; ibin<2*nbins1; ibin++)
   {
      cm_yt_nuc[ibin] = new double[2*nbins1];
      cm_yt_nonuc[ibin] = new double[2*nbins1];
   }
   for (int ibin=0; ibin<nbins1; ibin++)
   {
      cm_yp_nuc[ibin] = new double[nbins1];
      cm_ym_nuc[ibin] = new double[nbins1];
      cm_ch_nuc[ibin] = new double[nbins1];
      cm_yp_nonuc[ibin] = new double[nbins1];
      cm_ym_nonuc[ibin] = new double[nbins1];
      cm_ch_nonuc[ibin] = new double[nbins1];
   }
   for (int ibin=0; ibin<nbins1/2; ibin++)
   {
      cm_A1p_nuc[ibin] = new double[nbins1/2];
      cm_A1m_nuc[ibin] = new double[nbins1/2];
      cm_A3_nuc[ibin] = new double[nbins1/2];
      cm_A1p_nonuc[ibin] = new double[nbins1/2];
      cm_A1m_nonuc[ibin] = new double[nbins1/2];
      cm_A3_nonuc[ibin] = new double[nbins1/2];
   }

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
      vector< vector<double> > pull_eff_plus;
      vector< vector<double> > pull_eff_minus;
      vector< vector<double> > pull_eff_all;
      double pull_lumi;
      double *pull_ct10 = new double[NSET_CT10];
      double *pull_eps09 = new double[NSET_EPS09];
      double pull_scale[2];
      bool dosysts = (ipe>0);

      // pulls for efficiencies
      for (int ieff=0; ieff<eff1p.size(); ieff++)
      {
         vector<double> tmpvec, tmpvec2;
         if (efferr1p[ieff][0]>0) 
            // in this case no covariance between muon charges: one pull for each charge
         {
            for (int ibin=0; ibin<nbins1; ibin++)
            {
               tmpvec.push_back(gRandom->Gaus(0,dosysts*DOEFFS));
               tmpvec2.push_back(gRandom->Gaus(0,dosysts*DOEFFS));
            }
            pull_eff_all.push_back(tmpvec);
            pull_eff_plus.push_back(tmpvec);
            pull_eff_minus.push_back(tmpvec2);
         }
         else // in this case only one pull for both plus and minus charges
         {
            for (int ibin=0; ibin<nbins1; ibin++)
            {
               tmpvec.push_back(gRandom->Gaus(0,dosysts*DOEFFS));
            }
            pull_eff_all.push_back(tmpvec);
            pull_eff_plus.push_back(tmpvec);
            pull_eff_minus.push_back(tmpvec);
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
      double *smeared_nuc_muplus = new double[nbins1];
      double *smeared_nuc_muminus = new double[nbins1];
      double *smeared_nonuc_muplus = new double[nbins1];
      double *smeared_nonuc_muminus = new double[nbins1];

      for (int ibin=0; ibin<nbins1; ibin++)
      {
         smeared_nuc_muplus[ibin] = xsec_mup_nuc01[ibin]*APb*LUMI*(1+pull_lumi*LUMIERR)*1e-6*(etamax1[ibin]-etamin1[ibin]) * scalemod(xsec_mup_scale[ibin],pull_scale[0],pull_scale[1]);
         smeared_nuc_muplus[ibin] *= pdfmod(NSET_CT10,xsec_mup_nuc01[ibin],xsec_mup_nucp1_ct10[ibin],xsec_mup_nucm1_ct10[ibin],pull_ct10);
         smeared_nuc_muplus[ibin] *= pdfmod(NSET_EPS09,xsec_mup_nuc01[ibin],xsec_mup_nucp1_eps09[ibin],xsec_mup_nucm1_eps09[ibin],pull_eps09);
         smeared_nonuc_muplus[ibin] = xsec_mup_nonuc01[ibin]*APb*LUMI*(1+pull_lumi*LUMIERR)*1e-6*(etamax1[ibin]-etamin1[ibin]) * scalemod(xsec_mup_scale[ibin],pull_scale[0],pull_scale[1]);
         smeared_nonuc_muplus[ibin] *= pdfmod(NSET_CT10,xsec_mup_nonuc01[ibin],xsec_mup_nonucp1_ct10[ibin],xsec_mup_nonucm1_ct10[ibin],pull_ct10);
         smeared_nuc_muminus[ibin] = xsec_mum_nuc01[ibin]*APb*LUMI*(1+pull_lumi*LUMIERR)*1e-6*(etamax1[ibin]-etamin1[ibin]) * scalemod(xsec_mum_scale[ibin],pull_scale[0],pull_scale[1]);
         smeared_nuc_muminus[ibin] *= pdfmod(NSET_CT10,xsec_mum_nuc01[ibin],xsec_mum_nucp1_ct10[ibin],xsec_mum_nucm1_ct10[ibin],pull_ct10);
         smeared_nuc_muminus[ibin] *= pdfmod(NSET_EPS09,xsec_mum_nuc01[ibin],xsec_mum_nucp1_eps09[ibin],xsec_mum_nucm1_eps09[ibin],pull_eps09);
         smeared_nonuc_muminus[ibin] = xsec_mum_nonuc01[ibin]*APb*LUMI*(1+pull_lumi*LUMIERR)*1e-6*(etamax1[ibin]-etamin1[ibin]) * scalemod(xsec_mum_scale[ibin],pull_scale[0],pull_scale[1]);
         smeared_nonuc_muminus[ibin] *= pdfmod(NSET_CT10,xsec_mum_nonuc01[ibin],xsec_mum_nonucp1_ct10[ibin],xsec_mum_nonucm1_ct10[ibin],pull_ct10);
      }



      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      // 4. throw gaussian distributed new yields from new model. Width = nominal stat uncert * sqrt(modified model / nominal model)
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      double *pe_nuc_muplus = new double[nbins1];
      double *pe_nuc_muminus = new double[nbins1];
      double *pe_nonuc_muplus = new double[nbins1];
      double *pe_nonuc_muminus = new double[nbins1];

      for (int ibin=0; ibin<nbins1; ibin++)
      {
         double raw_nuc_muplus = smeared_nuc_muplus[ibin];
         double raw_nonuc_muplus = smeared_nonuc_muplus[ibin];
         double raw_nuc_muminus = smeared_nuc_muminus[ibin];
         double raw_nonuc_muminus = smeared_nonuc_muminus[ibin];

         for (int ieff=0; ieff<eff1p.size(); ieff++)
            if (efferr1p[ieff][ibin]>0)
            {
               raw_nuc_muplus *= (eff1p[ieff][ibin] + efferr1p[ieff][ibin]*pull_eff_plus[ieff][ibin]);
               raw_nonuc_muplus *= (eff1p[ieff][ibin] + efferr1p[ieff][ibin]*pull_eff_plus[ieff][ibin]);
               raw_nuc_muminus *= (eff1m[ieff][ibin] + efferr1m[ieff][ibin]*pull_eff_minus[ieff][ibin]);
               raw_nonuc_muminus *= (eff1m[ieff][ibin] + efferr1m[ieff][ibin]*pull_eff_minus[ieff][ibin]);
            }
            else
            {
               raw_nuc_muplus *= (eff1p[ieff][ibin] + efferr1p[ieff][ibin]*pull_eff_all[ieff][ibin]);
               raw_nonuc_muplus *= (eff1p[ieff][ibin] + efferr1p[ieff][ibin]*pull_eff_all[ieff][ibin]);
               raw_nuc_muminus *= (eff1m[ieff][ibin] + efferr1m[ieff][ibin]*pull_eff_all[ieff][ibin]);
               raw_nonuc_muminus *= (eff1m[ieff][ibin] + efferr1m[ieff][ibin]*pull_eff_all[ieff][ibin]);
            }

         pe_nuc_muplus[ibin] = gRandom->Gaus(raw_nuc_muplus,DOSTAT*staterr1p[ibin]*sqrt(raw_nuc_muplus/yields1p[ibin]));
         pe_nonuc_muplus[ibin] = gRandom->Gaus(raw_nonuc_muplus,DOSTAT*staterr1p[ibin]*sqrt(raw_nonuc_muplus/yields1p[ibin]));
         pe_nuc_muminus[ibin] = gRandom->Gaus(raw_nuc_muminus,DOSTAT*staterr1m[ibin]*sqrt(raw_nuc_muminus/yields1m[ibin]));
         pe_nonuc_muminus[ibin] = gRandom->Gaus(raw_nonuc_muminus,DOSTAT*staterr1m[ibin]*sqrt(raw_nonuc_muminus/yields1m[ibin]));

         for (int ieff=0; ieff<eff1p.size(); ieff++)
         {
            pe_nuc_muplus[ibin] /= eff1p[ieff][ibin];
            pe_nonuc_muplus[ibin] /= eff1p[ieff][ibin];
            pe_nuc_muminus[ibin] /= eff1m[ieff][ibin];
            pe_nonuc_muminus[ibin] /= eff1m[ieff][ibin];
         }
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      // 5. compute derived quantities (asymetries...) for both pseudo-data and models
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      for (int ibin=0; ibin<nbins1/2; ibin++)
      {
         int jp = (nbins1/2)+ibin;
         int jm = (nbins1/2)-ibin-1;

         double A_penuc = pe_nuc_muplus[jm]; double B_penuc = pe_nuc_muplus[jp];
         double C_penuc = pe_nuc_muminus[jm]; double D_penuc = pe_nuc_muminus[jp];
         double A_penonuc = pe_nonuc_muplus[jm]; double B_penonuc = pe_nonuc_muplus[jp];
         double C_penonuc = pe_nonuc_muminus[jm]; double D_penonuc = pe_nonuc_muminus[jp];

         double A_nuc = smeared_nuc_muplus[jm]; double B_nuc = smeared_nuc_muplus[jp];
         double C_nuc = smeared_nuc_muminus[jm]; double D_nuc = smeared_nuc_muminus[jp];
         double A_nonuc = smeared_nonuc_muplus[jm]; double B_nonuc = smeared_nonuc_muplus[jp];
         double C_nonuc = smeared_nonuc_muminus[jm]; double D_nonuc = smeared_nonuc_muminus[jp];

         yp_penuc[jp] = pe_nuc_muplus[jp];
         yp_penuc[jm] = pe_nuc_muplus[jm];
         ym_penuc[jp] = pe_nuc_muminus[jp];
         ym_penuc[jm] = pe_nuc_muminus[jm];
         yt_penuc[jp] = pe_nuc_muplus[jp];
         yt_penuc[jm] = pe_nuc_muplus[jm];
         yt_penuc[nbins1+jp] = pe_nuc_muminus[jp];
         yt_penuc[nbins1+jm] = pe_nuc_muminus[jm];
         Ch_penuc[jp] = (B_penuc-D_penuc)/(B_penuc+D_penuc);
         Ch_penuc[jm] = (A_penuc-C_penuc)/(A_penuc+C_penuc);
         A1p_penuc[ibin] = A_penuc/B_penuc;
         A1m_penuc[ibin] = C_penuc/D_penuc;
         A3_penuc[ibin] = (A_penuc+C_penuc)/(B_penuc+D_penuc);

         yp_penonuc[jp] = pe_nonuc_muplus[jp];
         yp_penonuc[jm] = pe_nonuc_muplus[jm];
         ym_penonuc[jp] = pe_nonuc_muminus[jp];
         ym_penonuc[jm] = pe_nonuc_muminus[jm];
         yt_penonuc[jp] = pe_nonuc_muplus[jp];
         yt_penonuc[jm] = pe_nonuc_muplus[jm];
         yt_penonuc[nbins1+jp] = pe_nonuc_muminus[jp];
         yt_penonuc[nbins1+jm] = pe_nonuc_muminus[jm];
         Ch_penonuc[jp] = (B_penonuc-D_penonuc)/(B_penonuc+D_penonuc);
         Ch_penonuc[jm] = (A_penonuc-C_penonuc)/(A_penonuc+C_penonuc);
         A1p_penonuc[ibin] = A_penonuc/B_penonuc;
         A1m_penonuc[ibin] = C_penonuc/D_penonuc;
         A3_penonuc[ibin] = (A_penonuc+C_penonuc)/(B_penonuc+D_penonuc);
      }

      // fill covariance matrices
      for (int ibin=0; ibin<2*nbins1; ibin++)
         for (int jbin=0; jbin<2*nbins1; jbin++)
         {
            cm_yt_nuc[ibin][jbin] += (yt_penuc[ibin]-avg_yt_nuc[ibin])*(yt_penuc[jbin]-avg_yt_nuc[jbin]);
            cm_yt_nonuc[ibin][jbin] += (yt_penonuc[ibin]-avg_yt_nonuc[ibin])*(yt_penonuc[jbin]-avg_yt_nonuc[jbin]);
         }
      for (int ibin=0; ibin<nbins1; ibin++)
         for (int jbin=0; jbin<nbins1; jbin++)
         {
            cm_yp_nuc[ibin][jbin] += (yp_penuc[ibin]-avg_yp_nuc[ibin])*(yp_penuc[jbin]-avg_yp_nuc[jbin]);
            cm_ym_nuc[ibin][jbin] += (ym_penuc[ibin]-avg_ym_nuc[ibin])*(ym_penuc[jbin]-avg_ym_nuc[jbin]);
            cm_ch_nuc[ibin][jbin] += (Ch_penuc[ibin]-avg_ch_nuc[ibin])*(Ch_penuc[jbin]-avg_ch_nuc[jbin]);
            cm_yp_nonuc[ibin][jbin] += (yp_penonuc[ibin]-avg_yp_nonuc[ibin])*(yp_penonuc[jbin]-avg_yp_nonuc[jbin]);
            cm_ym_nonuc[ibin][jbin] += (ym_penonuc[ibin]-avg_ym_nonuc[ibin])*(ym_penonuc[jbin]-avg_ym_nonuc[jbin]);
            cm_ch_nonuc[ibin][jbin] += (Ch_penonuc[ibin]-avg_ch_nonuc[ibin])*(Ch_penonuc[jbin]-avg_ch_nonuc[jbin]);
         }
      for (int ibin=0; ibin<nbins1/2; ibin++)
         for (int jbin=0; jbin<nbins1/2; jbin++)
         {
            cm_A1p_nuc[ibin][jbin] += (A1p_penuc[ibin]-avg_A1p_nuc[ibin])*(A1p_penuc[jbin]-avg_A1p_nuc[jbin]);
            cm_A1m_nuc[ibin][jbin] += (A1m_penuc[ibin]-avg_A1m_nuc[ibin])*(A1m_penuc[jbin]-avg_A1m_nuc[jbin]);
            cm_A3_nuc[ibin][jbin] += (A3_penuc[ibin]-avg_A3_nuc[ibin])*(A3_penuc[jbin]-avg_A3_nuc[jbin]);
            cm_A1p_nonuc[ibin][jbin] += (A1p_penonuc[ibin]-avg_A1p_nonuc[ibin])*(A1p_penonuc[jbin]-avg_A1p_nonuc[jbin]);
            cm_A1m_nonuc[ibin][jbin] += (A1m_penonuc[ibin]-avg_A1m_nonuc[ibin])*(A1m_penonuc[jbin]-avg_A1m_nonuc[jbin]);
            cm_A3_nonuc[ibin][jbin] += (A3_penonuc[ibin]-avg_A3_nonuc[ibin])*(A3_penonuc[jbin]-avg_A3_nonuc[jbin]);
         }




      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      // 8. fill tree
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      tr->Fill();

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      // 9. delete stuff
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      delete pull_ct10;
      delete pull_eps09;
      delete smeared_nuc_muplus;
      delete smeared_nuc_muminus;
      delete smeared_nonuc_muplus;
      delete smeared_nonuc_muminus;
      delete pe_nuc_muplus;
      delete pe_nuc_muminus;
      delete pe_nonuc_muplus;
      delete pe_nonuc_muminus;
   }


   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // compute covariance matrices
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////

   for (int ibin=0; ibin<2*nbins1; ibin++)
      for (int jbin=0; jbin<2*nbins1; jbin++)
      {
         cm_yt_nuc[ibin][jbin] = cm_yt_nuc[ibin][jbin]/((double) NPE-1);
         cm_yt_nonuc[ibin][jbin] = cm_yt_nonuc[ibin][jbin]/((double) NPE-1);
      }
   for (int ibin=0; ibin<nbins1; ibin++)
      for (int jbin=0; jbin<nbins1; jbin++)
      {
         cm_yp_nuc[ibin][jbin] = cm_yp_nuc[ibin][jbin]/((double) NPE-1);
         cm_ym_nuc[ibin][jbin] = cm_ym_nuc[ibin][jbin]/((double) NPE-1);
         cm_ch_nuc[ibin][jbin] = cm_ch_nuc[ibin][jbin]/((double) NPE-1);
         cm_yp_nonuc[ibin][jbin] = cm_yp_nonuc[ibin][jbin]/((double) NPE-1);
         cm_ym_nonuc[ibin][jbin] = cm_ym_nonuc[ibin][jbin]/((double) NPE-1);
         cm_ch_nonuc[ibin][jbin] = cm_ch_nonuc[ibin][jbin]/((double) NPE-1);
      }
   for (int ibin=0; ibin<nbins1/2; ibin++)
      for (int jbin=0; jbin<nbins1/2; jbin++)
      {
         cm_A1p_nuc[ibin][jbin] = cm_A1p_nuc[ibin][jbin]/((double) NPE-1);
         cm_A1m_nuc[ibin][jbin] = cm_A1m_nuc[ibin][jbin]/((double) NPE-1);
         cm_A3_nuc[ibin][jbin] = cm_A3_nuc[ibin][jbin]/((double) NPE-1);
         cm_A1p_nonuc[ibin][jbin] = cm_A1p_nonuc[ibin][jbin]/((double) NPE-1);
         cm_A1m_nonuc[ibin][jbin] = cm_A1m_nonuc[ibin][jbin]/((double) NPE-1);
         cm_A3_nonuc[ibin][jbin] = cm_A3_nonuc[ibin][jbin]/((double) NPE-1);
      }

   // ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // // print covariance matrices
   // ///////////////////////////////////////////////////////////////////////////////////////////////////////////

   // cout << "covariance matrix for total yield (with nuclear effects)" << endl;
   // for (int ibin=0; ibin<2*nbins1; ibin++)
   // {
   //    for (int jbin=0; jbin<2*nbins1; jbin++)
   //       cout << cm_yt_nuc[ibin][jbin] << "\t";
   //    cout << endl;
   // }
   // cout << endl;
   // cout << "covariance matrix for total yield (without nuclear effects)" << endl;
   // for (int ibin=0; ibin<2*nbins1; ibin++)
   // {
   //    for (int jbin=0; jbin<2*nbins1; jbin++)
   //       cout << cm_yt_nonuc[ibin][jbin] << "\t";
   //    cout << endl;
   // }
   // cout << endl;
   // cout << "covariance matrix for mu+ yield (with nuclear effects)" << endl;
   // for (int ibin=0; ibin<nbins1; ibin++)
   // {
   //    for (int jbin=0; jbin<nbins1; jbin++)
   //       cout << cm_yp_nuc[ibin][jbin] << "\t";
   //    cout << endl;
   // }
   // cout << endl;
   // cout << "covariance matrix for mu+ yield (without nuclear effects)" << endl;
   // for (int ibin=0; ibin<nbins1; ibin++)
   // {
   //    for (int jbin=0; jbin<nbins1; jbin++)
   //       cout << cm_yp_nonuc[ibin][jbin] << "\t";
   //    cout << endl;
   // }
   // cout << endl;
   // cout << "covariance matrix for mu- yield (with nuclear effects)" << endl;
   // for (int ibin=0; ibin<nbins1; ibin++)
   // {
   //    for (int jbin=0; jbin<nbins1; jbin++)
   //       cout << cm_ym_nuc[ibin][jbin] << "\t";
   //    cout << endl;
   // }
   // cout << endl;
   // cout << "covariance matrix for mu- yield (without nuclear effects)" << endl;
   // for (int ibin=0; ibin<nbins1; ibin++)
   // {
   //    for (int jbin=0; jbin<nbins1; jbin++)
   //       cout << cm_ym_nonuc[ibin][jbin] << "\t";
   //    cout << endl;
   // }
   // cout << endl;
   // cout << "covariance matrix for charge asymmetry (with nuclear effects)" << endl;
   // for (int ibin=0; ibin<nbins1; ibin++)
   // {
   //    for (int jbin=0; jbin<nbins1; jbin++)
   //       cout << cm_ch_nuc[ibin][jbin] << "\t";
   //    cout << endl;
   // }
   // cout << endl;
   // cout << "covariance matrix for charge asymmetry (without nuclear effects)" << endl;
   // for (int ibin=0; ibin<nbins1; ibin++)
   // {
   //    for (int jbin=0; jbin<nbins1; jbin++)
   //       cout << cm_ch_nonuc[ibin][jbin] << "\t";
   //    cout << endl;
   // }
   // cout << endl;
   // cout << "covariance matrix for A1+ asymmetry (with nuclear effects)" << endl;
   // for (int ibin=0; ibin<nbins1/2; ibin++)
   // {
   //    for (int jbin=0; jbin<nbins1/2; jbin++)
   //       cout << cm_A1p_nuc[ibin][jbin] << "\t";
   //    cout << endl;
   // }
   // cout << endl;
   // cout << "covariance matrix for A1+ asymmetry (without nuclear effects)" << endl;
   // for (int ibin=0; ibin<nbins1/2; ibin++)
   // {
   //    for (int jbin=0; jbin<nbins1/2; jbin++)
   //       cout << cm_A1p_nonuc[ibin][jbin] << "\t";
   //    cout << endl;
   // }
   // cout << endl;
   // cout << "covariance matrix for A1- asymmetry (with nuclear effects)" << endl;
   // for (int ibin=0; ibin<nbins1/2; ibin++)
   // {
   //    for (int jbin=0; jbin<nbins1/2; jbin++)
   //       cout << cm_A1m_nuc[ibin][jbin] << "\t";
   //    cout << endl;
   // }
   // cout << endl;
   // cout << "covariance matrix for A1- asymmetry (without nuclear effects)" << endl;
   // for (int ibin=0; ibin<nbins1/2; ibin++)
   // {
   //    for (int jbin=0; jbin<nbins1/2; jbin++)
   //       cout << cm_A1m_nonuc[ibin][jbin] << "\t";
   //    cout << endl;
   // }
   // cout << endl;
   // cout << "covariance matrix for A3 asymmetry (with nuclear effects)" << endl;
   // for (int ibin=0; ibin<nbins1/2; ibin++)
   // {
   //    for (int jbin=0; jbin<nbins1/2; jbin++)
   //       cout << cm_A3_nuc[ibin][jbin] << "\t";
   //    cout << endl;
   // }
   // cout << endl;
   // cout << "covariance matrix for A3 asymmetry (without nuclear effects)" << endl;
   // for (int ibin=0; ibin<nbins1/2; ibin++)
   // {
   //    for (int jbin=0; jbin<nbins1/2; jbin++)
   //       cout << cm_A3_nonuc[ibin][jbin] << "\t";
   //    cout << endl;
   // }
   // cout << endl;

   fout->mkdir("matrices","Correlation and covariance matrices");
   fout->cd("matrices");
   TH2F *hcov_yp_nuc = new TH2F("hcov_yp_nuc","hcov_yp_nuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcov_yp_nonuc = new TH2F("hcov_yp_nonuc","hcov_yp_nonuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcov_ym_nuc = new TH2F("hcov_ym_nuc","hcov_ym_nuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcov_ym_nonuc = new TH2F("hcov_ym_nonuc","hcov_ym_nonuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcov_yt_nuc = new TH2F("hcov_yt_nuc","hcov_yt_nuc",2*nbins1,0,2*nbins1,2*nbins1,0,2*nbins1);
   TH2F *hcov_yt_nonuc = new TH2F("hcov_yt_nonuc","hcov_yt_nonuc",2*nbins1,0,2*nbins1,2*nbins1,0,2*nbins1);
   TH2F *hcov_ch_nuc = new TH2F("hcov_ch_nuc","hcov_ch_nuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcov_ch_nonuc = new TH2F("hcov_ch_nonuc","hcov_ch_nonuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcov_A1p_nuc = new TH2F("hcov_A1p_nuc","hcov_A1p_nuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcov_A1p_nonuc = new TH2F("hcov_A1p_nonuc","hcov_A1p_nonuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcov_A1m_nuc = new TH2F("hcov_A1m_nuc","hcov_A1m_nuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcov_A1m_nonuc = new TH2F("hcov_A1m_nonuc","hcov_A1m_nonuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcov_A3_nuc = new TH2F("hcov_A3_nuc","hcov_A3_nuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcov_A3_nonuc = new TH2F("hcov_A3_nonuc","hcov_A3_nonuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcov2_yp_nuc = new TH2F("hcov2_yp_nuc","hcov2_yp_nuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcov2_yp_nonuc = new TH2F("hcov2_yp_nonuc","hcov2_yp_nonuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcov2_ym_nuc = new TH2F("hcov2_ym_nuc","hcov2_ym_nuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcov2_ym_nonuc = new TH2F("hcov2_ym_nonuc","hcov2_ym_nonuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcov2_yt_nuc = new TH2F("hcov2_yt_nuc","hcov2_yt_nuc",2*nbins1,0,2*nbins1,2*nbins1,0,2*nbins1);
   TH2F *hcov2_yt_nonuc = new TH2F("hcov2_yt_nonuc","hcov2_yt_nonuc",2*nbins1,0,2*nbins1,2*nbins1,0,2*nbins1);
   TH2F *hcov2_ch_nuc = new TH2F("hcov2_ch_nuc","hcov2_ch_nuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcov2_ch_nonuc = new TH2F("hcov2_ch_nonuc","hcov2_ch_nonuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcov2_A1p_nuc = new TH2F("hcov2_A1p_nuc","hcov2_A1p_nuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcov2_A1p_nonuc = new TH2F("hcov2_A1p_nonuc","hcov2_A1p_nonuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcov2_A1m_nuc = new TH2F("hcov2_A1m_nuc","hcov2_A1m_nuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcov2_A1m_nonuc = new TH2F("hcov2_A1m_nonuc","hcov2_A1m_nonuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcov2_A3_nuc = new TH2F("hcov2_A3_nuc","hcov2_A3_nuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcov2_A3_nonuc = new TH2F("hcov2_A3_nonuc","hcov2_A3_nonuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcor_yp_nuc = new TH2F("hcor_yp_nuc","hcor_yp_nuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcor_yp_nonuc = new TH2F("hcor_yp_nonuc","hcor_yp_nonuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcor_ym_nuc = new TH2F("hcor_ym_nuc","hcor_ym_nuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcor_ym_nonuc = new TH2F("hcor_ym_nonuc","hcor_ym_nonuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcor_yt_nuc = new TH2F("hcor_yt_nuc","hcor_yt_nuc",2*nbins1,0,2*nbins1,2*nbins1,0,2*nbins1);
   TH2F *hcor_yt_nonuc = new TH2F("hcor_yt_nonuc","hcor_yt_nonuc",2*nbins1,0,2*nbins1,2*nbins1,0,2*nbins1);
   TH2F *hcor_ch_nuc = new TH2F("hcor_ch_nuc","hcor_ch_nuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcor_ch_nonuc = new TH2F("hcor_ch_nonuc","hcor_ch_nonuc",nbins1,0,nbins1,nbins1,0,nbins1);
   TH2F *hcor_A1p_nuc = new TH2F("hcor_A1p_nuc","hcor_A1p_nuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcor_A1p_nonuc = new TH2F("hcor_A1p_nonuc","hcor_A1p_nonuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcor_A1m_nuc = new TH2F("hcor_A1m_nuc","hcor_A1m_nuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcor_A1m_nonuc = new TH2F("hcor_A1m_nonuc","hcor_A1m_nonuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcor_A3_nuc = new TH2F("hcor_A3_nuc","hcor_A3_nuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);
   TH2F *hcor_A3_nonuc = new TH2F("hcor_A3_nonuc","hcor_A3_nonuc",nbins1/2,0,nbins1/2,nbins1/2,0,nbins1/2);

   for (int ibin=0; ibin<nbins1; ibin++)
      for (int jbin=0; jbin<nbins1; jbin++)
      {
         hcov_yp_nuc->Fill(ibin,jbin,mysqrt(cm_yp_nuc[ibin][jbin]));
         hcov_yp_nonuc->Fill(ibin,jbin,mysqrt(cm_yp_nonuc[ibin][jbin]));
         hcov_ym_nuc->Fill(ibin,jbin,mysqrt(cm_ym_nuc[ibin][jbin]));
         hcov_ym_nonuc->Fill(ibin,jbin,mysqrt(cm_ym_nonuc[ibin][jbin]));
         hcov_ch_nuc->Fill(ibin,jbin,mysqrt(cm_ch_nuc[ibin][jbin]));
         hcov_ch_nonuc->Fill(ibin,jbin,mysqrt(cm_ch_nonuc[ibin][jbin]));
         hcov2_yp_nuc->Fill(ibin,jbin,cm_yp_nuc[ibin][jbin]);
         hcov2_yp_nonuc->Fill(ibin,jbin,cm_yp_nonuc[ibin][jbin]);
         hcov2_ym_nuc->Fill(ibin,jbin,cm_ym_nuc[ibin][jbin]);
         hcov2_ym_nonuc->Fill(ibin,jbin,cm_ym_nonuc[ibin][jbin]);
         hcov2_ch_nuc->Fill(ibin,jbin,cm_ch_nuc[ibin][jbin]);
         hcov2_ch_nonuc->Fill(ibin,jbin,cm_ch_nonuc[ibin][jbin]);
         hcor_yp_nuc->Fill(ibin,jbin,cm_yp_nuc[ibin][jbin]/sqrt(cm_yp_nuc[ibin][ibin]*cm_yp_nuc[jbin][jbin]));
         hcor_yp_nonuc->Fill(ibin,jbin,cm_yp_nonuc[ibin][jbin]/sqrt(cm_yp_nonuc[ibin][ibin]*cm_yp_nonuc[jbin][jbin]));
         hcor_ym_nuc->Fill(ibin,jbin,cm_ym_nuc[ibin][jbin]/sqrt(cm_ym_nuc[ibin][ibin]*cm_ym_nuc[jbin][jbin]));
         hcor_ym_nonuc->Fill(ibin,jbin,cm_ym_nonuc[ibin][jbin]/sqrt(cm_ym_nonuc[ibin][ibin]*cm_ym_nonuc[jbin][jbin]));
         hcor_ch_nuc->Fill(ibin,jbin,cm_ch_nuc[ibin][jbin]/sqrt(cm_ch_nuc[ibin][ibin]*cm_ch_nuc[jbin][jbin]));
         hcor_ch_nonuc->Fill(ibin,jbin,cm_ch_nonuc[ibin][jbin]/sqrt(cm_ch_nonuc[ibin][ibin]*cm_ch_nonuc[jbin][jbin]));
      }
   for (int ibin=0; ibin<2*nbins1; ibin++)
      for (int jbin=0; jbin<2*nbins1; jbin++)
      {
         hcov_yt_nuc->Fill(ibin,jbin,mysqrt(cm_yt_nuc[ibin][jbin]));
         hcov_yt_nonuc->Fill(ibin,jbin,mysqrt(cm_yt_nonuc[ibin][jbin]));
         hcov2_yt_nuc->Fill(ibin,jbin,cm_yt_nuc[ibin][jbin]);
         hcov2_yt_nonuc->Fill(ibin,jbin,cm_yt_nonuc[ibin][jbin]);
         hcor_yt_nuc->Fill(ibin,jbin,cm_yt_nuc[ibin][jbin]/sqrt(cm_yt_nuc[ibin][ibin]*cm_yt_nuc[jbin][jbin]));
         hcor_yt_nonuc->Fill(ibin,jbin,cm_yt_nonuc[ibin][jbin]/sqrt(cm_yt_nonuc[ibin][ibin]*cm_yt_nonuc[jbin][jbin]));
      }
   for (int ibin=0; ibin<nbins1/2; ibin++)
      for (int jbin=0; jbin<nbins1/2; jbin++)
      {
         hcov_A1p_nuc->Fill(ibin,jbin,mysqrt(cm_A1p_nuc[ibin][jbin]));
         hcov_A1p_nonuc->Fill(ibin,jbin,mysqrt(cm_A1p_nonuc[ibin][jbin]));
         hcov_A1m_nuc->Fill(ibin,jbin,mysqrt(cm_A1m_nuc[ibin][jbin]));
         hcov_A1m_nonuc->Fill(ibin,jbin,mysqrt(cm_A1m_nonuc[ibin][jbin]));
         hcov_A3_nuc->Fill(ibin,jbin,mysqrt(cm_A3_nuc[ibin][jbin]));
         hcov_A3_nonuc->Fill(ibin,jbin,mysqrt(cm_A3_nonuc[ibin][jbin]));
         hcov2_A1p_nuc->Fill(ibin,jbin,cm_A1p_nuc[ibin][jbin]);
         hcov2_A1p_nonuc->Fill(ibin,jbin,cm_A1p_nonuc[ibin][jbin]);
         hcov2_A1m_nuc->Fill(ibin,jbin,cm_A1m_nuc[ibin][jbin]);
         hcov2_A1m_nonuc->Fill(ibin,jbin,cm_A1m_nonuc[ibin][jbin]);
         hcov2_A3_nuc->Fill(ibin,jbin,cm_A3_nuc[ibin][jbin]);
         hcov2_A3_nonuc->Fill(ibin,jbin,cm_A3_nonuc[ibin][jbin]);
         hcor_A1p_nuc->Fill(ibin,jbin,cm_A1p_nuc[ibin][jbin]/sqrt(cm_A1p_nuc[ibin][ibin]*cm_A1p_nuc[jbin][jbin]));
         hcor_A1p_nonuc->Fill(ibin,jbin,cm_A1p_nonuc[ibin][jbin]/sqrt(cm_A1p_nonuc[ibin][ibin]*cm_A1p_nonuc[jbin][jbin]));
         hcor_A1m_nuc->Fill(ibin,jbin,cm_A1m_nuc[ibin][jbin]/sqrt(cm_A1m_nuc[ibin][ibin]*cm_A1m_nuc[jbin][jbin]));
         hcor_A1m_nonuc->Fill(ibin,jbin,cm_A1m_nonuc[ibin][jbin]/sqrt(cm_A1m_nonuc[ibin][ibin]*cm_A1m_nonuc[jbin][jbin]));
         hcor_A3_nuc->Fill(ibin,jbin,cm_A3_nuc[ibin][jbin]/sqrt(cm_A3_nuc[ibin][ibin]*cm_A3_nuc[jbin][jbin]));
         hcor_A3_nonuc->Fill(ibin,jbin,cm_A3_nonuc[ibin][jbin]/sqrt(cm_A3_nonuc[ibin][ibin]*cm_A3_nonuc[jbin][jbin]));
      }

   fout->cd();
   fout->mkdir("scales","scale variations");
   fout->cd("scales");

   for (int ibin=0; ibin<nbins1; ibin++)
   {
      TGraph2D *grplus = scalegraph(xsec_mup_scale[ibin]); grplus->SetName(Form("grplus%i",ibin)); grplus->Write();
      TGraph2D *grminus = scalegraph(xsec_mum_scale[ibin]); grminus->SetName(Form("grminus%i",ibin)); grminus->Write();
   }

   fout->Write();
   fout->Close();


   // now compute chi2 between data and theory, using the covariance matrices
   // first create the ROOT objects
   TMatrixDSym tcov_yp_nuc(nbins1), tcov_yp_nonuc(nbins1);
   TMatrixDSym tcov_ym_nuc(nbins1), tcov_ym_nonuc(nbins1);
   TMatrixDSym tcov_yt_nuc(2*nbins1), tcov_yt_nonuc(2*nbins1);
   TMatrixDSym tcov_ch_nuc(nbins1), tcov_ch_nonuc(nbins1);
   TMatrixDSym tcov_A1p_nuc(nbins1/2), tcov_A1p_nonuc(nbins1/2);
   TMatrixDSym tcov_A1m_nuc(nbins1/2), tcov_A1m_nonuc(nbins1/2);
   TMatrixDSym tcov_A3_nuc(nbins1/2), tcov_A3_nonuc(nbins1/2);

   for (int ibin=0; ibin<nbins1; ibin++)
      for (int jbin=0; jbin<nbins1; jbin++)
      {
         tcov_yp_nuc[ibin][jbin] = cm_yp_nuc[ibin][jbin];
         tcov_yp_nonuc[ibin][jbin] = cm_yp_nonuc[ibin][jbin];
         tcov_ym_nuc[ibin][jbin] = cm_ym_nuc[ibin][jbin];
         tcov_ym_nonuc[ibin][jbin] = cm_ym_nonuc[ibin][jbin];
         tcov_ch_nuc[ibin][jbin] = cm_ch_nuc[ibin][jbin];
         tcov_ch_nonuc[ibin][jbin] = cm_ch_nonuc[ibin][jbin];
      }
   for (int ibin=0; ibin<2*nbins1; ibin++)
      for (int jbin=0; jbin<2*nbins1; jbin++)
      {
         tcov_yt_nuc[ibin][jbin] = cm_yt_nuc[ibin][jbin];
         tcov_yt_nonuc[ibin][jbin] = cm_yt_nonuc[ibin][jbin];
      }
   for (int ibin=0; ibin<nbins1/2; ibin++)
      for (int jbin=0; jbin<nbins1/2; jbin++)
      {
         tcov_A1p_nuc[ibin][jbin] = cm_A1p_nuc[ibin][jbin];
         tcov_A1p_nonuc[ibin][jbin] = cm_A1p_nonuc[ibin][jbin];
         tcov_A1m_nuc[ibin][jbin] = cm_A1m_nuc[ibin][jbin];
         tcov_A1m_nonuc[ibin][jbin] = cm_A1m_nonuc[ibin][jbin];
         tcov_A3_nuc[ibin][jbin] = cm_A3_nuc[ibin][jbin];
         tcov_A3_nonuc[ibin][jbin] = cm_A3_nonuc[ibin][jbin];
      }

   // now invert the matrices
   TMatrixDSym tinvcov_yp_nuc = tcov_yp_nuc.Invert();
   TMatrixDSym tinvcov_yp_nonuc = tcov_yp_nonuc.Invert();
   TMatrixDSym tinvcov_ym_nuc = tcov_ym_nuc.Invert();
   TMatrixDSym tinvcov_ym_nonuc = tcov_ym_nonuc.Invert();
   TMatrixDSym tinvcov_yt_nuc = tcov_yt_nuc.Invert();
   TMatrixDSym tinvcov_yt_nonuc = tcov_yt_nonuc.Invert();
   TMatrixDSym tinvcov_ch_nuc = tcov_ch_nuc.Invert();
   TMatrixDSym tinvcov_ch_nonuc = tcov_ch_nonuc.Invert();
   TMatrixDSym tinvcov_A1p_nuc = tcov_A1p_nuc.Invert();
   TMatrixDSym tinvcov_A1p_nonuc = tcov_A1p_nonuc.Invert();
   TMatrixDSym tinvcov_A1m_nuc = tcov_A1m_nuc.Invert();
   TMatrixDSym tinvcov_A1m_nonuc = tcov_A1m_nonuc.Invert();
   TMatrixDSym tinvcov_A3_nuc = tcov_A3_nuc.Invert();
   TMatrixDSym tinvcov_A3_nonuc = tcov_A3_nonuc.Invert();

   // actually compute the chi2's
   double chi2_yp_nuc=0, chi2_yp_nonuc=0;
   double chi2_ym_nuc=0, chi2_ym_nonuc=0;
   double chi2_yt_nuc=0, chi2_yt_nonuc=0;
   double chi2_ch_nuc=0, chi2_ch_nonuc=0;
   double chi2_A1p_nuc=0, chi2_A1p_nonuc=0;
   double chi2_A1m_nuc=0, chi2_A1m_nonuc=0;
   double chi2_A3_nuc=0, chi2_A3_nonuc=0;
   double chi2_yp_thnuc=0, chi2_yp_thnonuc=0;
   double chi2_ym_thnuc=0, chi2_ym_thnonuc=0;
   double chi2_yt_thnuc=0, chi2_yt_thnonuc=0;
   double chi2_ch_thnuc=0, chi2_ch_thnonuc=0;
   double chi2_A1p_thnuc=0, chi2_A1p_thnonuc=0;
   double chi2_A1m_thnuc=0, chi2_A1m_thnonuc=0;
   double chi2_A3_thnuc=0, chi2_A3_thnonuc=0;

   for (int ibin=0; ibin<nbins1; ibin++)
      for (int jbin=0; jbin<nbins1; jbin++)
      {
         chi2_yp_nuc += (yp_data[ibin]-avg_yp_nuc[ibin])*tinvcov_yp_nuc[ibin][jbin]*(yp_data[jbin]-avg_yp_nuc[jbin]);
         chi2_yp_nonuc += (yp_data[ibin]-avg_yp_nonuc[ibin])*tinvcov_yp_nonuc[ibin][jbin]*(yp_data[jbin]-avg_yp_nonuc[jbin]);
         chi2_ym_nuc += (ym_data[ibin]-avg_ym_nuc[ibin])*tinvcov_ym_nuc[ibin][jbin]*(ym_data[jbin]-avg_ym_nuc[jbin]);
         chi2_ym_nonuc += (ym_data[ibin]-avg_ym_nonuc[ibin])*tinvcov_ym_nonuc[ibin][jbin]*(ym_data[jbin]-avg_ym_nonuc[jbin]);
         chi2_ch_nuc += (Ch_data[ibin]-avg_ch_nuc[ibin])*tinvcov_ch_nuc[ibin][jbin]*(Ch_data[jbin]-avg_ch_nuc[jbin]);
         chi2_ch_nonuc += (Ch_data[ibin]-avg_ch_nonuc[ibin])*tinvcov_ch_nonuc[ibin][jbin]*(Ch_data[jbin]-avg_ch_nonuc[jbin]);
         chi2_yp_thnuc += (avg_yp_nonuc[ibin]-avg_yp_nuc[ibin])*tinvcov_yp_nuc[ibin][jbin]*(avg_yp_nonuc[jbin]-avg_yp_nuc[jbin]);
         chi2_yp_thnonuc += (avg_yp_nonuc[ibin]-avg_yp_nuc[ibin])*tinvcov_yp_nonuc[ibin][jbin]*(avg_yp_nonuc[jbin]-avg_yp_nuc[jbin]);
         chi2_ym_thnuc += (avg_ym_nonuc[ibin]-avg_ym_nuc[ibin])*tinvcov_ym_nuc[ibin][jbin]*(avg_ym_nonuc[jbin]-avg_ym_nuc[jbin]);
         chi2_ym_thnonuc += (avg_ym_nonuc[ibin]-avg_ym_nuc[ibin])*tinvcov_ym_nonuc[ibin][jbin]*(avg_ym_nonuc[jbin]-avg_ym_nuc[jbin]);
         chi2_ch_thnuc += (avg_ch_nonuc[ibin]-avg_ch_nuc[ibin])*tinvcov_ch_nuc[ibin][jbin]*(avg_ch_nonuc[jbin]-avg_ch_nuc[jbin]);
         chi2_ch_thnonuc += (avg_ch_nonuc[ibin]-avg_ch_nuc[ibin])*tinvcov_ch_nonuc[ibin][jbin]*(avg_ch_nonuc[jbin]-avg_ch_nuc[jbin]);
      }

   for (int ibin=0; ibin<2*nbins1; ibin++)
      for (int jbin=0; jbin<2*nbins1; jbin++)
      {
         chi2_yt_nuc += (yt_data[ibin]-avg_yt_nuc[ibin])*tinvcov_yt_nuc[ibin][jbin]*(yt_data[jbin]-avg_yt_nuc[jbin]);
         chi2_yt_nonuc += (yt_data[ibin]-avg_yt_nonuc[ibin])*tinvcov_yt_nonuc[ibin][jbin]*(yt_data[jbin]-avg_yt_nonuc[jbin]);
         chi2_yt_thnuc += (avg_yt_nonuc[ibin]-avg_yt_nuc[ibin])*tinvcov_yt_nuc[ibin][jbin]*(avg_yt_nonuc[jbin]-avg_yt_nuc[jbin]);
         chi2_yt_thnonuc += (avg_yt_nonuc[ibin]-avg_yt_nuc[ibin])*tinvcov_yt_nonuc[ibin][jbin]*(avg_yt_nonuc[jbin]-avg_yt_nuc[jbin]);
      }

   for (int ibin=0; ibin<nbins1/2; ibin++)
      for (int jbin=0; jbin<nbins1/2; jbin++)
      {
         chi2_A1p_nuc += (A1p_data[ibin]-avg_A1p_nuc[ibin])*tinvcov_A1p_nuc[ibin][jbin]*(A1p_data[jbin]-avg_A1p_nuc[jbin]);
         chi2_A1p_nonuc += (A1p_data[ibin]-avg_A1p_nonuc[ibin])*tinvcov_A1p_nonuc[ibin][jbin]*(A1p_data[jbin]-avg_A1p_nonuc[jbin]);
         chi2_A1m_nuc += (A1m_data[ibin]-avg_A1m_nuc[ibin])*tinvcov_A1m_nuc[ibin][jbin]*(A1m_data[jbin]-avg_A1m_nuc[jbin]);
         chi2_A1m_nonuc += (A1m_data[ibin]-avg_A1m_nonuc[ibin])*tinvcov_A1m_nonuc[ibin][jbin]*(A1m_data[jbin]-avg_A1m_nonuc[jbin]);
         chi2_A3_nuc += (A3_data[ibin]-avg_A3_nuc[ibin])*tinvcov_A3_nuc[ibin][jbin]*(A3_data[jbin]-avg_A3_nuc[jbin]);
         chi2_A3_nonuc += (A3_data[ibin]-avg_A3_nonuc[ibin])*tinvcov_A3_nonuc[ibin][jbin]*(A3_data[jbin]-avg_A3_nonuc[jbin]);
         chi2_A1p_thnuc += (avg_A1p_nonuc[ibin]-avg_A1p_nuc[ibin])*tinvcov_A1p_nuc[ibin][jbin]*(avg_A1p_nonuc[jbin]-avg_A1p_nuc[jbin]);
         chi2_A1p_thnonuc += (avg_A1p_nonuc[ibin]-avg_A1p_nuc[ibin])*tinvcov_A1p_nonuc[ibin][jbin]*(avg_A1p_nonuc[jbin]-avg_A1p_nuc[jbin]);
         chi2_A1m_thnuc += (avg_A1m_nonuc[ibin]-avg_A1m_nuc[ibin])*tinvcov_A1m_nuc[ibin][jbin]*(avg_A1m_nonuc[jbin]-avg_A1m_nuc[jbin]);
         chi2_A1m_thnonuc += (avg_A1m_nonuc[ibin]-avg_A1m_nuc[ibin])*tinvcov_A1m_nonuc[ibin][jbin]*(avg_A1m_nonuc[jbin]-avg_A1m_nuc[jbin]);
         chi2_A3_thnuc += (avg_A3_nonuc[ibin]-avg_A3_nuc[ibin])*tinvcov_A3_nuc[ibin][jbin]*(avg_A3_nonuc[jbin]-avg_A3_nuc[jbin]);
         chi2_A3_thnonuc += (avg_A3_nonuc[ibin]-avg_A3_nuc[ibin])*tinvcov_A3_nonuc[ibin][jbin]*(avg_A3_nonuc[jbin]-avg_A3_nuc[jbin]);
      }

   for (int ibin=0; ibin<nbins1; ibin++) cout << yp_data[ibin] << " ";
   cout << endl;
   for (int ibin=0; ibin<nbins1; ibin++) cout << avg_yp_nuc[ibin] << " ";
   cout << endl;

   // print the results
   cout << "chi2 for W+ (assuming nuclear effects): " << chi2_yp_nuc << " (prob. " << TMath::Prob(chi2_yp_nuc,nbins1) << ") (exp. " << chi2_yp_thnuc << ", prob. " << TMath::Prob(chi2_yp_thnuc,nbins1) << ")" << endl;
   cout << "chi2 for W+ (assuming no nuclear effects): " << chi2_yp_nonuc << " (prob. " << TMath::Prob(chi2_yp_nonuc,nbins1) << ") (exp. " << chi2_yp_thnonuc << ", prob. " << TMath::Prob(chi2_yp_thnonuc,nbins1) << ")" << endl;
   cout << "chi2 for W- (assuming nuclear effects): " << chi2_ym_nuc << " (prob. " << TMath::Prob(chi2_ym_nuc,nbins1) << ") (exp. " << chi2_ym_thnuc << ", prob. " << TMath::Prob(chi2_ym_thnuc,nbins1) << ")" << endl;
   cout << "chi2 for W- (assuming no nuclear effects): " << chi2_ym_nonuc << " (prob. " << TMath::Prob(chi2_ym_nonuc,nbins1) << ") (exp. " << chi2_ym_thnonuc << ", prob. " << TMath::Prob(chi2_ym_thnonuc,nbins1) << ")" << endl;
   cout << "chi2 for W+ and W- (assuming nuclear effects): " << chi2_yt_nuc << " (prob. " << TMath::Prob(chi2_yt_nuc,2*nbins1) << ") (exp. " << chi2_yt_thnuc << ", prob. " << TMath::Prob(chi2_yt_thnuc,2*nbins1) << ")" << endl;
   cout << "chi2 for W+ and W- (assuming no nuclear effects): " << chi2_yt_nonuc << " (prob. " << TMath::Prob(chi2_yt_nonuc,2*nbins1) << ") (exp. " << chi2_yt_thnonuc << ", prob. " << TMath::Prob(chi2_yt_thnonuc,2*nbins1) << ")" << endl;
   cout << "chi2 for charge asymmetry (assuming nuclear effects): " << chi2_ch_nuc << " (prob. " << TMath::Prob(chi2_ch_nuc,nbins1) << ") (exp. " << chi2_ch_thnuc << ", prob. " << TMath::Prob(chi2_ch_thnuc,nbins1) << ")" << endl;
   cout << "chi2 for charge asymmetry (assuming no nuclear effects): " << chi2_ch_nonuc << " (prob. " << TMath::Prob(chi2_ch_nonuc,nbins1) << ") (exp. " << chi2_ch_thnonuc << ", prob. " << TMath::Prob(chi2_ch_thnonuc,nbins1) << ")" << endl;
   cout << "chi2 for A1+ asymmetry (assuming nuclear effects): " << chi2_A1p_nuc << " (prob. " << TMath::Prob(chi2_A1p_nuc,nbins1/2) << ") (exp. " << chi2_A1p_thnuc << ", prob. " << TMath::Prob(chi2_A1p_thnuc,nbins1/2) << ")" << endl;
   cout << "chi2 for A1+ asymmetry (assuming no nuclear effects): " << chi2_A1p_nonuc << " (prob. " << TMath::Prob(chi2_A1p_nonuc,nbins1/2) << ") (exp. " << chi2_A1p_thnonuc << ", prob. " << TMath::Prob(chi2_A1p_thnonuc,nbins1/2) << ")" << endl;
   cout << "chi2 for A1- asymmetry (assuming nuclear effects): " << chi2_A1m_nuc << " (prob. " << TMath::Prob(chi2_A1m_nuc,nbins1/2) << ") (exp. " << chi2_A1m_thnuc << ", prob. " << TMath::Prob(chi2_A1m_thnuc,nbins1/2) << ")" << endl;
   cout << "chi2 for A1- asymmetry (assuming no nuclear effects): " << chi2_A1m_nonuc << " (prob. " << TMath::Prob(chi2_A1m_nonuc,nbins1/2) << ") (exp. " << chi2_A1m_thnonuc << ", prob. " << TMath::Prob(chi2_A1m_thnonuc,nbins1/2) << ")" << endl;
   cout << "chi2 for A3 asymmetry (assuming nuclear effects): " << chi2_A3_nuc << " (prob. " << TMath::Prob(chi2_A3_nuc,nbins1/2) << ") (exp. " << chi2_A3_thnuc << ", prob. " << TMath::Prob(chi2_A3_thnuc,nbins1/2) << ")" << endl;
   cout << "chi2 for A3 asymmetry (assuming no nuclear effects): " << chi2_A3_nonuc << " (prob. " << TMath::Prob(chi2_A3_nonuc,nbins1/2) << ") (exp. " << chi2_A3_thnonuc << ", prob. " << TMath::Prob(chi2_A3_thnonuc,nbins1/2) << ")" << endl;
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
