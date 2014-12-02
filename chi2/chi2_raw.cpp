#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <math.h>

#include "TRandom3.h"
#include "TGraph2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphErrors.h"
#include "TMatrixDSym.h"
#include "TMath.h"

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
   double dyp_data[nbins1], dym_data[nbins1], dyt_data[2*nbins1];
   double dA1p_data[nbins1/2], dA1m_data[nbins1/2], dA3_data[nbins1/2];
   double dCh_data[nbins1];
   for (int ibin=0; ibin<nbins1/2; ibin++)
   {
      int jp = (nbins1/2)+ibin;
      int jm = (nbins1/2)-ibin-1;

      double A_data = yields1p[jm]; double B_data = yields1p[jp];
      double C_data = yields1m[jm]; double D_data = yields1m[jp];
      double dA_data = staterr1p[jm]; double dB_data = staterr1p[jp];
      double dC_data = staterr1m[jm]; double dD_data = staterr1m[jp];

      yp_data[jp] = yields1p[jp];
      yp_data[jm] = yields1p[jm];
      ym_data[jp] = yields1m[jp];
      ym_data[jm] = yields1m[jm];
      yt_data[jp] = yields1p[jp];
      yt_data[jm] = yields1p[jm];
      yt_data[jp+nbins1] = yields1m[jp];
      yt_data[jm+nbins1] = yields1m[jm];
      Ch_data[jp] = (B_data-D_data)/(B_data+D_data);
      Ch_data[jm] = (A_data-C_data)/(A_data+C_data);
      A1p_data[ibin] = A_data/B_data;
      A1m_data[ibin] = C_data/D_data;
      A3_data[ibin] = (A_data+C_data)/(B_data+D_data);

      dyp_data[jm] = dA_data;
      dyp_data[jp] = dB_data;
      dym_data[jm] = dC_data;
      dym_data[jp] = dD_data;
      dyt_data[jm] = dA_data;
      dyt_data[jp] = dB_data;
      dyt_data[jm+nbins1] = dC_data;
      dyt_data[jp+nbins1] = dD_data;
      dCh_data[jp] = sqrt((4*D_data*D_data/pow(B_data+D_data,4))*dB_data*dB_data + (4*B_data*B_data/pow(B_data+D_data,4)*dD_data*dD_data));
      dCh_data[jm] = sqrt((4*C_data*C_data/pow(A_data+C_data,4))*dA_data*dA_data + (4*A_data*A_data/pow(A_data+C_data,4)*dC_data*dC_data));
      dA1p_data[ibin] = fabs(A1p_data[ibin])*sqrt(pow(dA_data/A_data,2) + pow(dB_data/B_data,2));
      dA1m_data[ibin] = fabs(A1m_data[ibin])*sqrt(pow(dC_data/C_data,2) + pow(dD_data/D_data,2));
      dA3_data[ibin] = fabs(A3_data[ibin])*sqrt((dA_data*dA_data+dC_data*dC_data)/pow(A_data+C_data,2) + (dB_data*dB_data+dD_data*dD_data)/pow(B_data+D_data,2));
   }

   TGraphErrors *gyields_exp_statonly_1 = statchannel::graph(chan1p, chan1m, YIELDS_PM, false, false);  gyields_exp_statonly_1->SetName("gyields_exp_statonly_1"); gyields_exp_statonly_1->Write();
   TGraphErrors *gyieldsp_exp_statonly_1 = statchannel::graph(chan1p, chan1m, YIELDS, false, false);  gyieldsp_exp_statonly_1->SetName("gyieldsp_exp_statonly_1"); gyieldsp_exp_statonly_1->Write();
   TGraphErrors *gyieldsm_exp_statonly_1 = statchannel::graph(chan1m, chan1p, YIELDS, false, false);  gyieldsm_exp_statonly_1->SetName("gyieldsm_exp_statonly_1"); gyieldsm_exp_statonly_1->Write();
   TGraphErrors *gch_exp_statonly_1 = statchannel::graph(chan1p, chan1m, CH, false, false); gch_exp_statonly_1->SetName("gch_exp_statonly_1"); gch_exp_statonly_1->Write();
   TGraphErrors *gA1p_exp_statonly_1 = statchannel::graph(chan1p, chan1m, A1p, false, false); gA1p_exp_statonly_1->SetName("gA1p_exp_statonly_1"); gA1p_exp_statonly_1->Write();
   TGraphErrors *gA1m_exp_statonly_1 = statchannel::graph(chan1p, chan1m, A1m, false, false); gA1m_exp_statonly_1->SetName("gA1m_exp_statonly_1"); gA1m_exp_statonly_1->Write();
   TGraphErrors *gA3_exp_statonly_1 = statchannel::graph(chan1p, chan1m, A3, false, false); gA3_exp_statonly_1->SetName("gA3_exp_statonly_1"); gA3_exp_statonly_1->Write();


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

      for (int ieff=0; ieff<eff1p.size(); ieff++)
      {
         A_nuc *= eff1p[ieff][jm];
         B_nuc *= eff1p[ieff][jp];
         C_nuc *= eff1m[ieff][jm];
         D_nuc *= eff1m[ieff][jp];
         A_nonuc *= eff1p[ieff][jm];
         B_nonuc *= eff1p[ieff][jp];
         C_nonuc *= eff1m[ieff][jm];
         D_nonuc *= eff1m[ieff][jp];
      }

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
   cout << "central prediction W+ (w/o nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
      cout << xsec_mup_nonuc01[ibin]*APb*LUMI*1e-6*(etamax1[ibin]-etamin1[ibin]) << " ";
   cout << endl;
   cout << "central prediction W- (w/o nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
      cout << xsec_mum_nonuc01[ibin]*APb*LUMI*1e-6*(etamax1[ibin]-etamin1[ibin]) << " ";
   cout << endl;


   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // declaration of everything going to the ttree, and creation of branches
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////

   TTree *tr = new TTree("PEtree","pseudo-experiments");

   // chi2's
   double chi2_yp_data_thnuc;
   double chi2_ym_data_thnuc;
   double chi2_yt_data_thnuc;
   double chi2_Ch_data_thnuc;
   double chi2_A1p_data_thnuc;
   double chi2_A1m_data_thnuc;
   double chi2_A3_data_thnuc;
   double chi2_yp_data_thnonuc;
   double chi2_ym_data_thnonuc;
   double chi2_yt_data_thnonuc;
   double chi2_Ch_data_thnonuc;
   double chi2_A1p_data_thnonuc;
   double chi2_A1m_data_thnonuc;
   double chi2_A3_data_thnonuc;
   double chi2_yp_penuc_thnuc;
   double chi2_ym_penuc_thnuc;
   double chi2_yt_penuc_thnuc;
   double chi2_Ch_penuc_thnuc;
   double chi2_A1p_penuc_thnuc;
   double chi2_A1m_penuc_thnuc;
   double chi2_A3_penuc_thnuc;
   double chi2_yp_penuc_thnonuc;
   double chi2_ym_penuc_thnonuc;
   double chi2_yt_penuc_thnonuc;
   double chi2_Ch_penuc_thnonuc;
   double chi2_A1p_penuc_thnonuc;
   double chi2_A1m_penuc_thnonuc;
   double chi2_A3_penuc_thnonuc;
   double chi2_yp_penonuc_thnuc;
   double chi2_ym_penonuc_thnuc;
   double chi2_yt_penonuc_thnuc;
   double chi2_Ch_penonuc_thnuc;
   double chi2_A1p_penonuc_thnuc;
   double chi2_A1m_penonuc_thnuc;
   double chi2_A3_penonuc_thnuc;
   double chi2_yp_penonuc_thnonuc;
   double chi2_ym_penonuc_thnonuc;
   double chi2_yt_penonuc_thnonuc;
   double chi2_Ch_penonuc_thnonuc;
   double chi2_A1p_penonuc_thnonuc;
   double chi2_A1m_penonuc_thnonuc;
   double chi2_A3_penonuc_thnonuc;
   tr->Branch("chi2_yp_data_thnuc",&chi2_yp_data_thnuc,"chi2_yp_data_thnuc/D");
   tr->Branch("chi2_ym_data_thnuc",&chi2_ym_data_thnuc,"chi2_ym_data_thnuc/D");
   tr->Branch("chi2_yt_data_thnuc",&chi2_yt_data_thnuc,"chi2_yt_data_thnuc/D");
   tr->Branch("chi2_Ch_data_thnuc",&chi2_Ch_data_thnuc,"chi2_Ch_data_thnuc/D");
   tr->Branch("chi2_A1p_data_thnuc",&chi2_A1p_data_thnuc,"chi2_A1p_data_thnuc/D");
   tr->Branch("chi2_A1m_data_thnuc",&chi2_A1m_data_thnuc,"chi2_A1m_data_thnuc/D");
   tr->Branch("chi2_A3_data_thnuc",&chi2_A3_data_thnuc,"chi2_A3_data_thnuc/D");
   tr->Branch("chi2_yp_data_thnonuc",&chi2_yp_data_thnonuc,"chi2_yp_data_thnonuc/D");
   tr->Branch("chi2_ym_data_thnonuc",&chi2_ym_data_thnonuc,"chi2_ym_data_thnonuc/D");
   tr->Branch("chi2_yt_data_thnonuc",&chi2_yt_data_thnonuc,"chi2_yt_data_thnonuc/D");
   tr->Branch("chi2_Ch_data_thnonuc",&chi2_Ch_data_thnonuc,"chi2_Ch_data_thnonuc/D");
   tr->Branch("chi2_A1p_data_thnonuc",&chi2_A1p_data_thnonuc,"chi2_A1p_data_thnonuc/D");
   tr->Branch("chi2_A1m_data_thnonuc",&chi2_A1m_data_thnonuc,"chi2_A1m_data_thnonuc/D");
   tr->Branch("chi2_A3_data_thnonuc",&chi2_A3_data_thnonuc,"chi2_A3_data_thnonuc/D");
   tr->Branch("chi2_yp_penuc_thnuc",&chi2_yp_penuc_thnuc,"chi2_yp_penuc_thnuc/D");
   tr->Branch("chi2_ym_penuc_thnuc",&chi2_ym_penuc_thnuc,"chi2_ym_penuc_thnuc/D");
   tr->Branch("chi2_yt_penuc_thnuc",&chi2_yt_penuc_thnuc,"chi2_yt_penuc_thnuc/D");
   tr->Branch("chi2_Ch_penuc_thnuc",&chi2_Ch_penuc_thnuc,"chi2_Ch_penuc_thnuc/D");
   tr->Branch("chi2_A1p_penuc_thnuc",&chi2_A1p_penuc_thnuc,"chi2_A1p_penuc_thnuc/D");
   tr->Branch("chi2_A1m_penuc_thnuc",&chi2_A1m_penuc_thnuc,"chi2_A1m_penuc_thnuc/D");
   tr->Branch("chi2_A3_penuc_thnuc",&chi2_A3_penuc_thnuc,"chi2_A3_penuc_thnuc/D");
   tr->Branch("chi2_yp_penuc_thnonuc",&chi2_yp_penuc_thnonuc,"chi2_yp_penuc_thnonuc/D");
   tr->Branch("chi2_ym_penuc_thnonuc",&chi2_ym_penuc_thnonuc,"chi2_ym_penuc_thnonuc/D");
   tr->Branch("chi2_yt_penuc_thnonuc",&chi2_yt_penuc_thnonuc,"chi2_yt_penuc_thnonuc/D");
   tr->Branch("chi2_Ch_penuc_thnonuc",&chi2_Ch_penuc_thnonuc,"chi2_Ch_penuc_thnonuc/D");
   tr->Branch("chi2_A1p_penuc_thnonuc",&chi2_A1p_penuc_thnonuc,"chi2_A1p_penuc_thnonuc/D");
   tr->Branch("chi2_A1m_penuc_thnonuc",&chi2_A1m_penuc_thnonuc,"chi2_A1m_penuc_thnonuc/D");
   tr->Branch("chi2_A3_penuc_thnonuc",&chi2_A3_penuc_thnonuc,"chi2_A3_penuc_thnonuc/D");
   tr->Branch("chi2_yp_penonuc_thnuc",&chi2_yp_penonuc_thnuc,"chi2_yp_penonuc_thnuc/D");
   tr->Branch("chi2_ym_penonuc_thnuc",&chi2_ym_penonuc_thnuc,"chi2_ym_penonuc_thnuc/D");
   tr->Branch("chi2_yt_penonuc_thnuc",&chi2_yt_penonuc_thnuc,"chi2_yt_penonuc_thnuc/D");
   tr->Branch("chi2_Ch_penonuc_thnuc",&chi2_Ch_penonuc_thnuc,"chi2_Ch_penonuc_thnuc/D");
   tr->Branch("chi2_A1p_penonuc_thnuc",&chi2_A1p_penonuc_thnuc,"chi2_A1p_penonuc_thnuc/D");
   tr->Branch("chi2_A1m_penonuc_thnuc",&chi2_A1m_penonuc_thnuc,"chi2_A1m_penonuc_thnuc/D");
   tr->Branch("chi2_A3_penonuc_thnuc",&chi2_A3_penonuc_thnuc,"chi2_A3_penonuc_thnuc/D");
   tr->Branch("chi2_yp_penonuc_thnonuc",&chi2_yp_penonuc_thnonuc,"chi2_yp_penonuc_thnonuc/D");
   tr->Branch("chi2_ym_penonuc_thnonuc",&chi2_ym_penonuc_thnonuc,"chi2_ym_penonuc_thnonuc/D");
   tr->Branch("chi2_yt_penonuc_thnonuc",&chi2_yt_penonuc_thnonuc,"chi2_yt_penonuc_thnonuc/D");
   tr->Branch("chi2_Ch_penonuc_thnonuc",&chi2_Ch_penonuc_thnonuc,"chi2_Ch_penonuc_thnonuc/D");
   tr->Branch("chi2_A1p_penonuc_thnonuc",&chi2_A1p_penonuc_thnonuc,"chi2_A1p_penonuc_thnonuc/D");
   tr->Branch("chi2_A1m_penonuc_thnonuc",&chi2_A1m_penonuc_thnonuc,"chi2_A1m_penonuc_thnonuc/D");
   tr->Branch("chi2_A3_penonuc_thnonuc",&chi2_A3_penonuc_thnonuc,"chi2_A3_penonuc_thnonuc/D");


   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // correlation matrices
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
      vector<double> Sk_nuc, Sk_nonuc;
      bool dosysts = (ipe>0);

      // pulls for efficiencies
      for (int ieff=0; ieff<eff1p.size(); ieff++)
      {
         vector<double> tmpvec, tmpvec2;
         if (eff1p[ieff][0]>0) 
            // in this case no correlation between muon charges: one pull for each charge
         {
            for (int ibin=0; ibin<nbins1; ibin++)
            {
               Sk_nonuc.push_back(gRandom->Gaus(0,dosysts*DOEFFS));
               tmpvec.push_back(Sk_nonuc.back());
               Sk_nonuc.push_back(gRandom->Gaus(0,dosysts*DOEFFS));
               tmpvec2.push_back(Sk_nonuc.back());
            }
            pull_eff_plus.push_back(tmpvec);
            pull_eff_minus.push_back(tmpvec2);
         }
         else // in this case only one pull for both plus and minus charges
         {
            for (int ibin=0; ibin<nbins1; ibin++)
            {
               Sk_nonuc.push_back(gRandom->Gaus(0,dosysts*DOEFFS));
               tmpvec.push_back(Sk_nonuc.back());
            }
            pull_eff_all.push_back(tmpvec);
         }
      }

      // pull for lumi
      pull_lumi = gRandom->Gaus(0,dosysts*DOLUMI);
      Sk_nonuc.push_back(pull_lumi);

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
      {
         Sk_nonuc.push_back(tmp[iscale]);
         pull_scale[iscale] = tmp[iscale];
      }

      // pulls for ct10
      for (int iset=0; iset<NSET_CT10; iset++)
      {
         Sk_nonuc.push_back(gRandom->Gaus(0,dosysts*DOCT10));
         pull_ct10[iset] = Sk_nonuc.back();
      }

      // pulls for EPS09
      Sk_nuc = Sk_nonuc;
      for (int iset=0; iset<NSET_EPS09; iset++)
      {
         Sk_nuc.push_back(gRandom->Gaus(0,dosysts*DOEPS09));
         pull_eps09[iset] = Sk_nuc.back();
      }




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
         // smeared_nonuc_muplus[ibin] *= pdfmod(NSET_EPS09,xsec_mup_nonuc01[ibin],xsec_mup_nonucp1_eps09[ibin],xsec_mup_nonucm1_eps09[ibin],pull_eps09);
         smeared_nuc_muminus[ibin] = xsec_mum_nuc01[ibin]*APb*LUMI*(1+pull_lumi*LUMIERR)*1e-6*(etamax1[ibin]-etamin1[ibin]) * scalemod(xsec_mum_scale[ibin],pull_scale[0],pull_scale[1]);
         smeared_nuc_muminus[ibin] *= pdfmod(NSET_CT10,xsec_mum_nuc01[ibin],xsec_mum_nucp1_ct10[ibin],xsec_mum_nucm1_ct10[ibin],pull_ct10);
         smeared_nuc_muminus[ibin] *= pdfmod(NSET_EPS09,xsec_mum_nuc01[ibin],xsec_mum_nucp1_eps09[ibin],xsec_mum_nucm1_eps09[ibin],pull_eps09);
         smeared_nonuc_muminus[ibin] = xsec_mum_nonuc01[ibin]*APb*LUMI*(1+pull_lumi*LUMIERR)*1e-6*(etamax1[ibin]-etamin1[ibin]) * scalemod(xsec_mum_scale[ibin],pull_scale[0],pull_scale[1]);
         smeared_nonuc_muminus[ibin] *= pdfmod(NSET_CT10,xsec_mum_nonuc01[ibin],xsec_mum_nonucp1_ct10[ibin],xsec_mum_nonucm1_ct10[ibin],pull_ct10);
         // smeared_nonuc_muminus[ibin] *= pdfmod(NSET_EPS09,xsec_mum_nonuc01[ibin],xsec_mum_nonucp1_eps09[ibin],xsec_mum_nonucm1_eps09[ibin],pull_eps09);
         for (int ieff=0; ieff<eff1p.size(); ieff++)
            if (eff1p[ieff][ibin]>0)
            {
               smeared_nuc_muplus[ibin] *= (eff1p[ieff][ibin] + efferr1p[ieff][ibin]*pull_eff_plus[ieff][ibin]);
               smeared_nonuc_muplus[ibin] *= (eff1p[ieff][ibin] + efferr1p[ieff][ibin]*pull_eff_plus[ieff][ibin]);
               smeared_nuc_muminus[ibin] *= (eff1m[ieff][ibin] + efferr1m[ieff][ibin]*pull_eff_minus[ieff][ibin]);
               smeared_nonuc_muminus[ibin] *= (eff1m[ieff][ibin] + efferr1m[ieff][ibin]*pull_eff_minus[ieff][ibin]);
            }
            else
            {
               smeared_nuc_muplus[ibin] *= (eff1p[ieff][ibin] + efferr1p[ieff][ibin]*pull_eff_all[ieff][ibin]);
               smeared_nonuc_muplus[ibin] *= (eff1p[ieff][ibin] + efferr1p[ieff][ibin]*pull_eff_all[ieff][ibin]);
               smeared_nuc_muminus[ibin] *= (eff1m[ieff][ibin] + efferr1m[ieff][ibin]*pull_eff_all[ieff][ibin]);
               smeared_nonuc_muminus[ibin] *= (eff1m[ieff][ibin] + efferr1m[ieff][ibin]*pull_eff_all[ieff][ibin]);
            }
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
         pe_nuc_muplus[ibin] = gRandom->Gaus(smeared_nuc_muplus[ibin],DOSTAT*staterr1p[ibin]*sqrt(smeared_nuc_muplus[ibin]/yields1p[ibin]));
         pe_nonuc_muplus[ibin] = gRandom->Gaus(smeared_nonuc_muplus[ibin],DOSTAT*staterr1p[ibin]*sqrt(smeared_nonuc_muplus[ibin]/yields1p[ibin]));
         pe_nuc_muminus[ibin] = gRandom->Gaus(smeared_nuc_muminus[ibin],DOSTAT*staterr1m[ibin]*sqrt(smeared_nuc_muminus[ibin]/yields1m[ibin]));
         pe_nonuc_muminus[ibin] = gRandom->Gaus(smeared_nonuc_muminus[ibin],DOSTAT*staterr1m[ibin]*sqrt(smeared_nonuc_muminus[ibin]/yields1m[ibin]));
      }

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      // 5. compute derived quantities (asymetries...) for both pseudo-data and models
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      double yp_penuc[nbins1], ym_penuc[nbins1], yt_penuc[2*nbins1];
      double A1p_penuc[nbins1/2], A1m_penuc[nbins1/2], A3_penuc[nbins1/2];
      double Ch_penuc[nbins1];
      double yp_penonuc[nbins1], ym_penonuc[nbins1], yt_penonuc[2*nbins1];
      double A1p_penonuc[nbins1/2], A1m_penonuc[nbins1/2], A3_penonuc[nbins1/2];
      double Ch_penonuc[nbins1];
      double yp_nuc[nbins1], ym_nuc[nbins1], yt_nuc[2*nbins1];
      double A1p_nuc[nbins1/2], A1m_nuc[nbins1/2], A3_nuc[nbins1/2];
      double Ch_nuc[nbins1];
      double yp_nonuc[nbins1], ym_nonuc[nbins1], yt_nonuc[2*nbins1];
      double A1p_nonuc[nbins1/2], A1m_nonuc[nbins1/2], A3_nonuc[nbins1/2];
      double Ch_nonuc[nbins1];
      double dyp_nuc[nbins1], dym_nuc[nbins1], dyt_nuc[2*nbins1];
      double dA1p_nuc[nbins1/2], dA1m_nuc[nbins1/2], dA3_nuc[nbins1/2];
      double dCh_nuc[nbins1];
      double dyp_nonuc[nbins1], dym_nonuc[nbins1], dyt_nonuc[2*nbins1];
      double dA1p_nonuc[nbins1/2], dA1m_nonuc[nbins1/2], dA3_nonuc[nbins1/2];
      double dCh_nonuc[nbins1];

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

         double dA_nuc = staterr1p[jm]*sqrt(smeared_nuc_muplus[jm]/yields1p[jm]); 
         double dB_nuc = staterr1p[jp]*sqrt(smeared_nuc_muplus[jp]/yields1p[jp]);
         double dC_nuc = staterr1m[jm]*sqrt(smeared_nuc_muminus[jm]/yields1m[jm]); 
         double dD_nuc = staterr1m[jp]*sqrt(smeared_nuc_muminus[jp]/yields1m[jp]);
         double dA_nonuc = staterr1p[jm]*sqrt(pe_nonuc_muplus[jm]/yields1p[jm]); 
         double dB_nonuc = staterr1p[jp]*sqrt(pe_nonuc_muplus[jp]/yields1p[jp]);
         double dC_nonuc = staterr1m[jm]*sqrt(pe_nonuc_muminus[jm]/yields1m[jm]); 
         double dD_nonuc = staterr1m[jp]*sqrt(pe_nonuc_muminus[jp]/yields1m[jp]);

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

         yp_nuc[jp] = smeared_nuc_muplus[jp];
         yp_nuc[jm] = smeared_nuc_muplus[jm];
         ym_nuc[jp] = smeared_nuc_muminus[jp];
         ym_nuc[jm] = smeared_nuc_muminus[jm];
         yt_nuc[jp] = smeared_nuc_muplus[jp];
         yt_nuc[jm] = smeared_nuc_muplus[jm];
         yt_nuc[nbins1+jp] = smeared_nuc_muminus[jp];
         yt_nuc[nbins1+jm] = smeared_nuc_muminus[jm];
         Ch_nuc[jp] = (B_nuc-D_nuc)/(B_nuc+D_nuc);
         Ch_nuc[jm] = (A_nuc-C_nuc)/(A_nuc+C_nuc);
         A1p_nuc[ibin] = A_nuc/B_nuc;
         A1m_nuc[ibin] = C_nuc/D_nuc;
         A3_nuc[ibin] = (A_nuc+C_nuc)/(B_nuc+D_nuc);

         yp_nonuc[jp] = smeared_nonuc_muplus[jp];
         yp_nonuc[jm] = smeared_nonuc_muplus[jm];
         ym_nonuc[jp] = smeared_nonuc_muminus[jp];
         ym_nonuc[jm] = smeared_nonuc_muminus[jm];
         yt_nonuc[jp] = smeared_nonuc_muplus[jp];
         yt_nonuc[jm] = smeared_nonuc_muplus[jm];
         yt_nonuc[nbins1+jp] = smeared_nonuc_muminus[jp];
         yt_nonuc[nbins1+jm] = smeared_nonuc_muminus[jm];
         Ch_nonuc[jp] = (B_nonuc-D_nonuc)/(B_nonuc+D_nonuc);
         Ch_nonuc[jm] = (A_nonuc-C_nonuc)/(A_nonuc+C_nonuc);
         A1p_nonuc[ibin] = A_nonuc/B_nonuc;
         A1m_nonuc[ibin] = C_nonuc/D_nonuc;
         A3_nonuc[ibin] = (A_nonuc+C_nonuc)/(B_nonuc+D_nonuc);

         dyp_nonuc[jm] = dA_nonuc;
         dyp_nonuc[jp] = dB_nonuc;
         dym_nonuc[jm] = dC_nonuc;
         dym_nonuc[jp] = dD_nonuc;
         dyt_nonuc[jm] = dA_nonuc;
         dyt_nonuc[jp] = dB_nonuc;
         dyt_nonuc[nbins1+jm] = dC_nonuc;
         dyt_nonuc[nbins1+jp] = dD_nonuc;
         dCh_nonuc[jp] = sqrt((4*D_nonuc*D_nonuc/pow(B_nonuc+D_nonuc,4))*dB_nonuc*dB_nonuc + (4*B_nonuc*B_nonuc/pow(B_nonuc+D_nonuc,4)*dD_nonuc*dD_nonuc));
         dCh_nonuc[jm] = sqrt((4*C_nonuc*C_nonuc/pow(A_nonuc+C_nonuc,4))*dA_nonuc*dA_nonuc + (4*A_nonuc*A_nonuc/pow(A_nonuc+C_nonuc,4)*dC_nonuc*dC_nonuc));
         dA1p_nonuc[ibin] = fabs(A1p_penonuc[ibin])*sqrt(pow(dA_nonuc/A_penonuc,2) + pow(dB_nonuc/B_penonuc,2));
         dA1m_nonuc[ibin] = fabs(A1m_penonuc[ibin])*sqrt(pow(dC_nonuc/C_penonuc,2) + pow(dD_nonuc/D_penonuc,2));
         dA3_nonuc[ibin] = fabs(A3_penonuc[ibin])*sqrt((dA_nonuc*dA_nonuc+dC_nonuc*dC_nonuc)/pow(A_penonuc+C_penonuc,2) + (dB_nonuc*dB_nonuc+dD_nonuc*dD_nonuc)/pow(B_penonuc+D_penonuc,2));

         dyp_nuc[jm] = dA_nuc;
         dyp_nuc[jp] = dB_nuc;
         dym_nuc[jm] = dC_nuc;
         dym_nuc[jp] = dD_nuc;
         dyt_nuc[jm] = dA_nuc;
         dyt_nuc[jp] = dB_nuc;
         dyt_nuc[nbins1+jm] = dC_nuc;
         dyt_nuc[nbins1+jp] = dD_nuc;
         dCh_nuc[jp] = sqrt((4*D_nuc*D_nuc/pow(B_nuc+D_nuc,4))*dB_nuc*dB_nuc + (4*B_nuc*B_nuc/pow(B_nuc+D_nuc,4)*dD_nuc*dD_nuc));
         dCh_nuc[jm] = sqrt((4*C_nuc*C_nuc/pow(A_nuc+C_nuc,4))*dA_nuc*dA_nuc + (4*A_nuc*A_nuc/pow(A_nuc+C_nuc,4)*dC_nuc*dC_nuc));
         dA1p_nuc[ibin] = fabs(A1p_penuc[ibin])*sqrt(pow(dA_nuc/A_penuc,2) + pow(dB_nuc/B_penuc,2));
         dA1m_nuc[ibin] = fabs(A1m_penuc[ibin])*sqrt(pow(dC_nuc/C_penuc,2) + pow(dD_nuc/D_penuc,2));
         dA3_nuc[ibin] = fabs(A3_penuc[ibin])*sqrt((dA_nuc*dA_nuc+dC_nuc*dC_nuc)/pow(A_penonuc+C_penonuc,2) + (dB_nuc*dB_nuc+dD_nuc*dD_nuc)/pow(B_penonuc+D_penonuc,2));
      }

      // fill correlation matrices
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
      // 6. compute new chi2 from eq 4.5 of CDF doc, for +/- yields and derived quantities
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      chi2_ym_data_thnuc = chi2(nbins1,ym_data,dym_data,ym_nuc,Sk_nuc);
      chi2_yp_data_thnuc = chi2(nbins1,yp_data,dyp_data,yp_nuc,Sk_nuc);
      chi2_yt_data_thnuc = chi2(nbins1,yt_data,dyt_data,yt_nuc,Sk_nuc);
      chi2_Ch_data_thnuc = chi2(nbins1,Ch_data,dCh_data,Ch_nuc,Sk_nuc);
      chi2_A1p_data_thnuc = chi2(nbins1/2,A1p_data,dA1p_data,A1p_nuc,Sk_nuc);
      chi2_A1m_data_thnuc = chi2(nbins1/2,A1m_data,dA1m_data,A1m_nuc,Sk_nuc);
      chi2_A3_data_thnuc = chi2(nbins1/2,A3_data,dA3_data,A3_nuc,Sk_nuc);
      chi2_ym_data_thnonuc = chi2(nbins1,ym_data,dym_data,ym_nonuc,Sk_nonuc);
      chi2_yp_data_thnonuc = chi2(nbins1,yp_data,dyp_data,yp_nonuc,Sk_nonuc);
      chi2_yt_data_thnonuc = chi2(nbins1,yt_data,dyt_data,yt_nonuc,Sk_nonuc);
      chi2_Ch_data_thnonuc = chi2(nbins1,Ch_data,dCh_data,Ch_nonuc,Sk_nonuc);
      chi2_A1p_data_thnonuc = chi2(nbins1/2,A1p_data,dA1p_data,A1p_nonuc,Sk_nonuc);
      chi2_A1m_data_thnonuc = chi2(nbins1/2,A1m_data,dA1m_data,A1m_nonuc,Sk_nonuc);
      chi2_A3_data_thnonuc = chi2(nbins1/2,A3_data,dA3_data,A3_nonuc,Sk_nonuc);
      chi2_ym_penuc_thnuc = chi2(nbins1,ym_penuc,dym_nuc,ym_nuc,Sk_nuc);
      chi2_yp_penuc_thnuc = chi2(nbins1,yp_penuc,dyp_nuc,yp_nuc,Sk_nuc);
      chi2_yt_penuc_thnuc = chi2(nbins1,yt_penuc,dyt_nuc,yt_nuc,Sk_nuc);
      chi2_Ch_penuc_thnuc = chi2(nbins1,Ch_penuc,dCh_nuc,Ch_nuc,Sk_nuc);
      chi2_A1p_penuc_thnuc = chi2(nbins1/2,A1p_penuc,dA1p_nuc,A1p_nuc,Sk_nuc);
      chi2_A1m_penuc_thnuc = chi2(nbins1/2,A1m_penuc,dA1m_nuc,A1m_nuc,Sk_nuc);
      chi2_A3_penuc_thnuc = chi2(nbins1/2,A3_penuc,dA3_nuc,A3_nuc,Sk_nuc);
      chi2_ym_penuc_thnonuc = chi2(nbins1,ym_penuc,dym_nuc,ym_nonuc,Sk_nonuc);
      chi2_yp_penuc_thnonuc = chi2(nbins1,yp_penuc,dyp_nuc,yp_nonuc,Sk_nonuc);
      chi2_yt_penuc_thnonuc = chi2(nbins1,yt_penuc,dyt_nuc,yt_nonuc,Sk_nonuc);
      chi2_Ch_penuc_thnonuc = chi2(nbins1,Ch_penuc,dCh_nuc,Ch_nonuc,Sk_nonuc);
      chi2_A1p_penuc_thnonuc = chi2(nbins1/2,A1p_penuc,dA1p_nuc,A1p_nonuc,Sk_nonuc);
      chi2_A1m_penuc_thnonuc = chi2(nbins1/2,A1m_penuc,dA1m_nuc,A1m_nonuc,Sk_nonuc);
      chi2_A3_penuc_thnonuc = chi2(nbins1/2,A3_penuc,dA3_nuc,A3_nonuc,Sk_nonuc);
      chi2_ym_penonuc_thnuc = chi2(nbins1,ym_penonuc,dym_nonuc,ym_nuc,Sk_nuc);
      chi2_yp_penonuc_thnuc = chi2(nbins1,yp_penonuc,dyp_nonuc,yp_nuc,Sk_nuc);
      chi2_yt_penonuc_thnuc = chi2(nbins1,yt_penonuc,dyt_nonuc,yt_nuc,Sk_nuc);
      chi2_Ch_penonuc_thnuc = chi2(nbins1,Ch_penonuc,dCh_nonuc,Ch_nuc,Sk_nuc);
      chi2_A1p_penonuc_thnuc = chi2(nbins1/2,A1p_penonuc,dA1p_nonuc,A1p_nuc,Sk_nuc);
      chi2_A1m_penonuc_thnuc = chi2(nbins1/2,A1m_penonuc,dA1m_nonuc,A1m_nuc,Sk_nuc);
      chi2_A3_penonuc_thnuc = chi2(nbins1/2,A3_penonuc,dA3_nonuc,A3_nuc,Sk_nuc);
      chi2_ym_penonuc_thnonuc = chi2(nbins1,ym_penonuc,dym_nonuc,ym_nonuc,Sk_nonuc);
      chi2_yp_penonuc_thnonuc = chi2(nbins1,yp_penonuc,dyp_nonuc,yp_nonuc,Sk_nonuc);
      chi2_yt_penonuc_thnonuc = chi2(nbins1,yt_penonuc,dyt_nonuc,yt_nonuc,Sk_nonuc);
      chi2_Ch_penonuc_thnonuc = chi2(nbins1,Ch_penonuc,dCh_nonuc,Ch_nonuc,Sk_nonuc);
      chi2_A1p_penonuc_thnonuc = chi2(nbins1/2,A1p_penonuc,dA1p_nonuc,A1p_nonuc,Sk_nonuc);
      chi2_A1m_penonuc_thnonuc = chi2(nbins1/2,A1m_penonuc,dA1m_nonuc,A1m_nonuc,Sk_nonuc);
      chi2_A3_penonuc_thnonuc = chi2(nbins1/2,A3_penonuc,dA3_nonuc,A3_nonuc,Sk_nonuc);


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
   // compute correlation matrices
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////

   for (int ibin=0; ibin<2*nbins1; ibin++)
      for (int jbin=0; jbin<2*nbins1; jbin++)
      {
         cm_yt_nuc[ibin][jbin] = mysqrt(cm_yt_nuc[ibin][jbin]/((double) NPE-1));
         cm_yt_nonuc[ibin][jbin] = mysqrt(cm_yt_nonuc[ibin][jbin]/((double) NPE-1));
      }
   for (int ibin=0; ibin<nbins1; ibin++)
      for (int jbin=0; jbin<nbins1; jbin++)
      {
         cm_yp_nuc[ibin][jbin] = mysqrt(cm_yp_nuc[ibin][jbin]/((double) NPE-1));
         cm_ym_nuc[ibin][jbin] = mysqrt(cm_ym_nuc[ibin][jbin]/((double) NPE-1));
         cm_ch_nuc[ibin][jbin] = mysqrt(cm_ch_nuc[ibin][jbin]/((double) NPE-1));
         cm_yp_nonuc[ibin][jbin] = mysqrt(cm_yp_nonuc[ibin][jbin]/((double) NPE-1));
         cm_ym_nonuc[ibin][jbin] = mysqrt(cm_ym_nonuc[ibin][jbin]/((double) NPE-1));
         cm_ch_nonuc[ibin][jbin] = mysqrt(cm_ch_nonuc[ibin][jbin]/((double) NPE-1));
      }
   for (int ibin=0; ibin<nbins1/2; ibin++)
      for (int jbin=0; jbin<nbins1/2; jbin++)
      {
         cm_A1p_nuc[ibin][jbin] = mysqrt(cm_A1p_nuc[ibin][jbin]/((double) NPE-1));
         cm_A1m_nuc[ibin][jbin] = mysqrt(cm_A1m_nuc[ibin][jbin]/((double) NPE-1));
         cm_A3_nuc[ibin][jbin] = mysqrt(cm_A3_nuc[ibin][jbin]/((double) NPE-1));
         cm_A1p_nonuc[ibin][jbin] = mysqrt(cm_A1p_nonuc[ibin][jbin]/((double) NPE-1));
         cm_A1m_nonuc[ibin][jbin] = mysqrt(cm_A1m_nonuc[ibin][jbin]/((double) NPE-1));
         cm_A3_nonuc[ibin][jbin] = mysqrt(cm_A3_nonuc[ibin][jbin]/((double) NPE-1));
      }

   ///////////////////////////////////////////////////////////////////////////////////////////////////////////
   // print correlation matrices
   ///////////////////////////////////////////////////////////////////////////////////////////////////////////

   cout << "correlation matrix for total yield (with nuclear effects)" << endl;
   for (int ibin=0; ibin<2*nbins1; ibin++)
   {
      for (int jbin=0; jbin<2*nbins1; jbin++)
         cout << cm_yt_nuc[ibin][jbin] << "\t";
      cout << endl;
   }
   cout << endl;
   cout << "correlation matrix for total yield (without nuclear effects)" << endl;
   for (int ibin=0; ibin<2*nbins1; ibin++)
   {
      for (int jbin=0; jbin<2*nbins1; jbin++)
         cout << cm_yt_nonuc[ibin][jbin] << "\t";
      cout << endl;
   }
   cout << endl;
   cout << "correlation matrix for mu+ yield (with nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
   {
      for (int jbin=0; jbin<nbins1; jbin++)
         cout << cm_yp_nuc[ibin][jbin] << "\t";
      cout << endl;
   }
   cout << endl;
   cout << "correlation matrix for mu+ yield (without nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
   {
      for (int jbin=0; jbin<nbins1; jbin++)
         cout << cm_yp_nonuc[ibin][jbin] << "\t";
      cout << endl;
   }
   cout << endl;
   cout << "correlation matrix for mu- yield (with nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
   {
      for (int jbin=0; jbin<nbins1; jbin++)
         cout << cm_ym_nuc[ibin][jbin] << "\t";
      cout << endl;
   }
   cout << endl;
   cout << "correlation matrix for mu- yield (without nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
   {
      for (int jbin=0; jbin<nbins1; jbin++)
         cout << cm_ym_nonuc[ibin][jbin] << "\t";
      cout << endl;
   }
   cout << endl;
   cout << "correlation matrix for charge asymmetry (with nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
   {
      for (int jbin=0; jbin<nbins1; jbin++)
         cout << cm_ch_nuc[ibin][jbin] << "\t";
      cout << endl;
   }
   cout << endl;
   cout << "correlation matrix for charge asymmetry (without nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1; ibin++)
   {
      for (int jbin=0; jbin<nbins1; jbin++)
         cout << cm_ch_nonuc[ibin][jbin] << "\t";
      cout << endl;
   }
   cout << endl;
   cout << "correlation matrix for A1+ asymmetry (with nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
   {
      for (int jbin=0; jbin<nbins1/2; jbin++)
         cout << cm_A1p_nuc[ibin][jbin] << "\t";
      cout << endl;
   }
   cout << endl;
   cout << "correlation matrix for A1+ asymmetry (without nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
   {
      for (int jbin=0; jbin<nbins1/2; jbin++)
         cout << cm_A1p_nonuc[ibin][jbin] << "\t";
      cout << endl;
   }
   cout << endl;
   cout << "correlation matrix for A1- asymmetry (with nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
   {
      for (int jbin=0; jbin<nbins1/2; jbin++)
         cout << cm_A1m_nuc[ibin][jbin] << "\t";
      cout << endl;
   }
   cout << endl;
   cout << "correlation matrix for A1- asymmetry (without nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
   {
      for (int jbin=0; jbin<nbins1/2; jbin++)
         cout << cm_A1m_nonuc[ibin][jbin] << "\t";
      cout << endl;
   }
   cout << endl;
   cout << "correlation matrix for A3 asymmetry (with nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
   {
      for (int jbin=0; jbin<nbins1/2; jbin++)
         cout << cm_A3_nuc[ibin][jbin] << "\t";
      cout << endl;
   }
   cout << endl;
   cout << "correlation matrix for A3 asymmetry (without nuclear effects)" << endl;
   for (int ibin=0; ibin<nbins1/2; ibin++)
   {
      for (int jbin=0; jbin<nbins1/2; jbin++)
         cout << cm_A3_nonuc[ibin][jbin] << "\t";
      cout << endl;
   }
   cout << endl;

   fout->Write();
   fout->Close();


   // now compute chi2 between data and theory, using the correlation matrices
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

   for (int ibin=0; ibin<nbins1; ibin++)
      for (int jbin=0; jbin<nbins1; jbin++)
      {
         chi2_yp_nuc += (yp_data[ibin]-avg_yp_nuc[ibin])*tinvcov_yp_nuc[ibin][jbin]*(yp_data[jbin]-avg_yp_nuc[jbin]);
         chi2_yp_nonuc += (yp_data[ibin]-avg_yp_nonuc[ibin])*tinvcov_yp_nonuc[ibin][jbin]*(yp_data[jbin]-avg_yp_nonuc[jbin]);
         chi2_ym_nuc += (ym_data[ibin]-avg_ym_nuc[ibin])*tinvcov_ym_nuc[ibin][jbin]*(ym_data[jbin]-avg_ym_nuc[jbin]);
         chi2_ym_nonuc += (ym_data[ibin]-avg_ym_nonuc[ibin])*tinvcov_ym_nonuc[ibin][jbin]*(ym_data[jbin]-avg_ym_nonuc[jbin]);
         chi2_ch_nuc += (Ch_data[ibin]-avg_ch_nuc[ibin])*tinvcov_ch_nuc[ibin][jbin]*(Ch_data[jbin]-avg_ch_nuc[jbin]);
         chi2_ch_nonuc += (Ch_data[ibin]-avg_ch_nonuc[ibin])*tinvcov_ch_nonuc[ibin][jbin]*(Ch_data[jbin]-avg_ch_nonuc[jbin]);
      }

   for (int ibin=0; ibin<2*nbins1; ibin++)
      for (int jbin=0; jbin<2*nbins1; jbin++)
      {
         chi2_yt_nuc += (yt_data[ibin]-avg_yt_nuc[ibin])*tinvcov_yt_nuc[ibin][jbin]*(yt_data[jbin]-avg_yt_nuc[jbin]);
         chi2_yt_nonuc += (yt_data[ibin]-avg_yt_nonuc[ibin])*tinvcov_yt_nonuc[ibin][jbin]*(yt_data[jbin]-avg_yt_nonuc[jbin]);
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
      }

   // print the results
   cout << "chi2 for W+ (assuming nuclear effects): " << chi2_yp_nuc << " (prob. " << TMath::Prob(chi2_yp_nuc,nbins1) << ")" << endl;
   cout << "chi2 for W+ (assuming no nuclear effects): " << chi2_yp_nonuc << " (prob. " << TMath::Prob(chi2_yp_nonuc,nbins1) << ")" << endl;
   cout << "chi2 for W- (assuming nuclear effects): " << chi2_ym_nuc << " (prob. " << TMath::Prob(chi2_ym_nuc,nbins1) << ")" << endl;
   cout << "chi2 for W- (assuming no nuclear effects): " << chi2_ym_nonuc << " (prob. " << TMath::Prob(chi2_ym_nonuc,nbins1) << ")" << endl;
   cout << "chi2 for W+ and W- (assuming nuclear effects): " << chi2_yt_nuc << " (prob. " << TMath::Prob(chi2_yt_nuc,2*nbins1) << ")" << endl;
   cout << "chi2 for W+ and W- (assuming no nuclear effects): " << chi2_yt_nonuc << " (prob. " << TMath::Prob(chi2_yt_nonuc,2*nbins1) << ")" << endl;
   cout << "chi2 for charge asymmetry (assuming nuclear effects): " << chi2_ch_nuc << " (prob. " << TMath::Prob(chi2_ch_nuc,nbins1) << ")" << endl;
   cout << "chi2 for charge asymmetry (assuming no nuclear effects): " << chi2_ch_nonuc << " (prob. " << TMath::Prob(chi2_ch_nonuc,nbins1) << ")" << endl;
   cout << "chi2 for A1+ asymmetry (assuming nuclear effects): " << chi2_A1p_nuc << " (prob. " << TMath::Prob(chi2_A1p_nuc,nbins1/2) << ")" << endl;
   cout << "chi2 for A1+ asymmetry (assuming no nuclear effects): " << chi2_A1p_nonuc << " (prob. " << TMath::Prob(chi2_A1p_nonuc,nbins1/2) << ")" << endl;
   cout << "chi2 for A1- asymmetry (assuming nuclear effects): " << chi2_A1m_nuc << " (prob. " << TMath::Prob(chi2_A1m_nuc,nbins1/2) << ")" << endl;
   cout << "chi2 for A1- asymmetry (assuming no nuclear effects): " << chi2_A1m_nonuc << " (prob. " << TMath::Prob(chi2_A1m_nonuc,nbins1/2) << ")" << endl;
   cout << "chi2 for A3 asymmetry (assuming nuclear effects): " << chi2_A3_nuc << " (prob. " << TMath::Prob(chi2_A3_nuc,nbins1/2) << ")" << endl;
   cout << "chi2 for A3 asymmetry (assuming no nuclear effects): " << chi2_A3_nonuc << " (prob. " << TMath::Prob(chi2_A3_nonuc,nbins1/2) << ")" << endl;
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
   // return a>0 ? sqrt(a) : -sqrt(-a);
   return a;
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
