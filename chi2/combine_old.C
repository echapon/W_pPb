#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "TGraphErrors.h"
#include "TFile.h"
#include "TH2.h"
#include "TMatrixDSym.h"

using namespace std;

void combine_blue(TGraphErrors *gr1_stat, TGraphErrors *gr1_syst, TGraphErrors *gr2_stat, TGraphErrors *gr2_syst, TGraphErrors *gr12_stat, TGraphErrors *gr12_syst)
{
   int nbins = gr1_stat->GetN();
   for (int i=0; i<nbins; i++)
   {
      double x1, y1, x2, y2;
      gr1_stat->GetPoint(i,x1,y1);
      gr2_stat->GetPoint(i,x2,y2);
      if (x1 != x2) cout << "Warning, in bin " << i << ": x1 = " << x1 << ", x2 = " << x2 << endl;
      double ex1 = gr1_stat->GetErrorX(i);
      double ex2 = gr2_stat->GetErrorX(i);
      double ey1_stat = gr1_stat->GetErrorY(i);
      double ey1_syst = gr1_syst->GetErrorY(i);
      double ey2_stat = gr2_stat->GetErrorY(i);
      double ey2_syst = gr2_syst->GetErrorY(i);

      double ey2_1 = ey1_stat*ey1_stat+ey1_syst*ey1_syst;
      double ey2_2 = ey2_stat*ey2_stat+ey2_syst*ey2_syst;

      double wt1 = ey2_2/(ey2_1+ey2_2);
      double wt2 = ey2_1/(ey2_1+ey2_2);

      double y12 = wt1*y1 + wt2*y2;
      double ey12_stat = sqrt(wt1*wt1*ey1_stat*ey1_stat + wt2*wt2*ey2_stat*ey2_stat);
      double ey12_syst = sqrt(wt1*wt1*ey1_syst*ey1_syst + wt2*wt2*ey2_syst*ey2_syst);
      cout << "bin " << i << ": combine " << 
         y1 << " +- " << ey1_stat << " +- " << ey1_syst << " (weight " << wt1 << "), " << 
         y2 << " +- " << ey2_stat << " +- " << ey2_syst << " (weight " << wt2 << ") => " << 
         y12 << " +- " << ey12_stat << " +- " << ey12_syst << endl;

      gr12_stat->SetPoint(i,x1,y12);
      gr12_syst->SetPoint(i,x1,y12);
      gr12_stat->SetPointError(i,ex1,ey12_stat);
      gr12_syst->SetPointError(i,ex1,ey12_syst);
   }
}

void combine_blue(const char* infilename, const char* outfilename)
{
   TFile *fin = new TFile(infilename);
   TFile *fout = new TFile(outfilename,"RECREATE");

   vector<string> names;
   names.push_back("yields");
   names.push_back("yieldsp");
   names.push_back("yieldsm");
   names.push_back("ch");
   names.push_back("A1p");
   names.push_back("A1m");
   names.push_back("A3");
   names.push_back("A4");

   for (int i=0; i<names.size(); i++)
   {
      cout << "Combining " << names[i] << endl;
      TGraphErrors *gr1_stat = (TGraphErrors*) fin->Get(Form("g%s_exp_statonly_1",names[i].c_str()));
      TGraphErrors *gr1_syst = (TGraphErrors*) fin->Get(Form("g%s_exp_1",names[i].c_str()));
      TGraphErrors *gr2_stat = (TGraphErrors*) fin->Get(Form("g%s_exp_statonly_2",names[i].c_str()));
      TGraphErrors *gr2_syst = (TGraphErrors*) fin->Get(Form("g%s_exp_2",names[i].c_str()));
      int nbins = gr1_stat->GetN();
      TGraphErrors *gr12_stat = new TGraphErrors(nbins);
      TGraphErrors *gr12_syst = new TGraphErrors(nbins);

      fout->cd();
      combine_blue(gr1_stat, gr1_syst, gr2_stat, gr2_syst, gr12_stat, gr12_syst);
      fout->cd();
      gr12_stat->SetName(Form("g%s_exp_statonly_1",names[i].c_str())); gr12_stat->Write();
      gr12_syst->SetName(Form("g%s_exp_1",names[i].c_str())); gr12_syst->Write();
   }

   fout->Write();
   fout->Close();
   fin->Close();
}

double chi2(TGraphErrors *grstat, const char* label, const char* input_dta, TH2F *hcov, bool donuc)
{
   int nbins = grstat->GetN();
   double *y = grstat->GetY();

   ifstream file(input_dta);

   double dummy;
   double *A_CT10 = new double[nbins];
   double *A_EPS09 = new double[nbins];

   for (int ibin=0; ibin<nbins; ibin++)
   {
      int ibin2=ibin;
      if (label=="CA" || label=="yieldsWp" || label=="yieldsWm") ibin2 = nbins-1-ibin;
      file >> dummy >> dummy >> A_CT10[ibin2] >> dummy >> dummy >> A_EPS09[ibin2] >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;

      if (label == "yieldsWp" || label == "yieldsWm")
      {
         double deta = 1.;
         // if (channel_number==1 && (ibin2==0 || ibin2==nbins-1)) deta = 0.4;
         A_CT10[ibin2] *= 208.*deta*LUMI*1e-6;
         A_EPS09[ibin2] *= 208.*deta*LUMI*1e-6;
      }
   }

   double *th = donuc ? A_EPS09 : A_CT10;

   TMatrixDSym tcov(nbins);
   for (int ibin=0; ibin<nbins; ibin++)
      for (int jbins=0; jbin<nbins; jbin++)
         tcov[ibin][jbin] = hcov->GetBinContent(ibin+1,jbin+1);

   TMatrixDSym tinvcov = tcov.Invert();

   double chi2;
   for (int ibin=0; ibin<nbins; ibin++)
      for (int jbin=0; jbin<nbins; jbin++)
         chi2 += (y[ibin]-th[ibin])*tinvcov[ibin][jbin]*(y[jbin]-th[jbin]);

   return chi2;
}

void chi2(const char* file_graphs, const char* file_matrices)
{
   TFile *fgr = new TFile(file_graphs);
   TFile *fmat = new TFile(file_matrices);

   vector<string> names;
   names.push_back("yieldsp");
   names.push_back("yieldsm");
   names.push_back("ch");
   names.push_back("A1p");
   names.push_back("A1m");
   names.push_back("A3");
   names.push_back("A4");

   for (int i=0; i<names.size(); i++)
   {
      TGraphErrors *gr = (TGraphErrors*) fgr->Get(Form("g%s_exp_1",names[i].c_str()));
      string label,label2;
      if (string(names[i])==string("A1p")) {label="A1p"; label2="A1p";}
      else if (string(names[i])==string("A1m")) {label="A1m"; label2="A1m";}
      else if (string(names[i])==string("ch")) {label="CA"; label2="ch";}
      else if (string(names[i])==string("A3")) {label="A3"; label2="A3";}
      else if (string(names[i])==string("A4")) {label="A4"; label2="A4";}
      else if (string(names[i])==string("yieldsp")) {label="yieldsWp"; label2="yp";}
      else if (string(names[i])==string("yieldsm")) {label="yieldsWm"; label2="ym";}
      string input_dta = string(label) + string("_elect.dta");
      TH2F *hnuc = (TH2F*) fmat->Get(Form("matrices/hcov_%s_"));
   }
}
