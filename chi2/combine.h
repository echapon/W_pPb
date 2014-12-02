#ifndef combine_h
#define combine_h

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "TGraphErrors.h"
#include "TFile.h"
#include "TH2.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TRandom.h"

#include "lumierror.h"

using namespace std;

void combine_blue(double val1, double err1, double val2, double err2);
void cat(int n1, double *arr1, int n2, double* arr2, double *out);
TGraphErrors *combine(TGraphErrors *g1, TGraphErrors *g2);
void combine_blue(TGraphErrors *gr12, TH2F *cov12, TGraphErrors *grcombo, TH2F *covcombo);
void combine_blue(TGraphErrors *gr1_stat, TGraphErrors *gr1_syst, TGraphErrors *gr2_stat, TGraphErrors *gr2_syst, TGraphErrors *gr12_stat, TGraphErrors *gr12_syst);
void combine_blue_nocor(const char* infilename, const char* outfilename);
void combine_blue(const char* infilename, const char* outfilename, bool build_covmat_from_graphs);
double chi2(unsigned int nbins, double *y, const char* label, const char* input_dta, TH2F *hcov, bool donuc);
double chi2(TGraphErrors *grstat, const char* label, const char* input_dta, TH2F *hcov, bool donuc);
void chi2(const char* file_graphs, const char* file_matrices_exp, const char* file_matrices_th, int nchan=1);
void chi2(const char* file_graphs, const char* file_matrices_th, int nchan=1);
double combine_and_chi2(int nbins, double *y12, const char* label, const char* input_dta, TH2F *hcov_exp, TH2F *hcov_th, bool donuc);
TH1F* chi2_expected(TGraphErrors *gr_stat, TGraphErrors *gr_syst, TGraph *grth1, TGraph *grth2, TH2F *hcov_th1, TH2F *hcov_th2, TH2F *hcov_exp, int npe=10000);
TH1F* chi2_expected(TGraphErrors *gr_stat, TGraphErrors *gr_syst, const char* label, TH2F *hcov_th1, TH2F *hcov_th2, TH2F *hcov_exp, bool donuc, int npe=10000);
TH1F* chi2_expected(const char* file_exp, const char* file_th, const char* label, bool pred_nuc, int npe=10000);

#endif // ifndef combine_h
