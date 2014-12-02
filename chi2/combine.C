#ifndef combine_C
#define combine_C

#include "combine.h"
#include "build_covmat.C"

void combine_blue(double val1, double err1, double val2, double err2)
{
   double w1 = err2*err2/(err1*err1+err2*err2);
   double w2 = err1*err1/(err1*err1+err2*err2);

   cout << w1*val1+w2*val2 << "\\pm" << err1*err2/sqrt(err1*err1+err2*err2) << endl;
}

void cat(int n1, double *arr1, int n2, double* arr2, double *out)
{
   for (int i=0; i<n1; i++)
      out[i] = arr1[i];
   for (int i=0; i<n2; i++)
      out[n1+i] = arr2[i];
}

TGraphErrors *combine(TGraphErrors *g1, TGraphErrors *g2)
{
   int nbins1 = g1->GetN();
   int nbins2 = g2->GetN();

   double *x1 = g1->GetX();
   double *y1 = g1->GetY();
   double *ex1 = g1->GetEX();
   double *ey1 = g1->GetEY();
   double *x2 = g2->GetX();
   double *y2 = g2->GetY();
   double *ex2 = g2->GetEX();
   double *ey2 = g2->GetEY();

   double *x = new double[nbins1+nbins2];
   cat(nbins1,x1,nbins2,x2,x);
   double *y = new double[nbins1+nbins2];
   cat(nbins1,y1,nbins2,y2,y);
   double *ex = new double[nbins1+nbins2];
   cat(nbins1,ex1,nbins2,ex2,ex);
   double *ey = new double[nbins1+nbins2];
   cat(nbins1,ey1,nbins2,ey2,ey);

   return new TGraphErrors(nbins1+nbins2,x,y,ex,ey);
}

void combine_blue(TGraphErrors *gr12, TH2F *cov12, TGraphErrors *grcombo, TH2F *covcombo)
{
   unsigned int nbins = gr12->GetN()/2;
   double* x = gr12->GetX();
   double* y = gr12->GetY();
   double* ex = gr12->GetEX();
   double* ey = gr12->GetEY();

   TMatrixD U(2*nbins,nbins);
   TMatrixD Ut(nbins,2*nbins);
   TMatrixD M(2*nbins,2*nbins);
   TMatrixD Minv(2*nbins,2*nbins);
   TMatrixD UMU(nbins,nbins);
   TMatrixD UMUinv(nbins,nbins);
   TMatrixD UM(nbins,2*nbins);
   TMatrixD lambda(nbins,2*nbins);

   // U matrix
   for (unsigned int i=0; i<2*nbins; i++)
      for (unsigned int j=0; j<nbins; j++)
      {
         U[i][j]=0;
         if (i==j || i==j+nbins) U[i][j]=1;
      }
   Ut = Ut.Transpose(U);

   // M matrix
   for (unsigned int i=0; i<2*nbins; i++)
      for (unsigned int j=0; j<2*nbins; j++)
         M[i][j] = cov12->GetBinContent(i+1,j+1);

   // Minv
   Minv = M;
   Minv.Invert();

   // UMU
   UMU = Ut*Minv*U;
   UMUinv = UMU;
   UMUinv.Invert();
   UM = Ut*Minv;

   // weights lambda
   for (unsigned int alpha=0; alpha<nbins; alpha++)
      for (unsigned int i=0; i<2*nbins; i++)
      {
         lambda[alpha][i]=0;
         for (unsigned int beta=0; beta<nbins; beta++)
            lambda[alpha][i] += UMUinv[alpha][beta]*UM[beta][i];
      }

   for (unsigned int alpha=0; alpha<nbins; alpha++)
      for (unsigned int beta=0; beta<nbins; beta++)
         covcombo->SetBinContent(alpha+1,beta+1,UMUinv[alpha][beta]);

   for (unsigned int alpha=0; alpha<nbins; alpha++)
   {
      double newx = x[alpha];
      double newy = 0;
      for (unsigned int i=0; i<2*nbins; i++) newy += lambda[alpha][i] * y[i];
      double newex = ex[alpha];
      double newey = 0;
      for (unsigned int i=0; i<2*nbins; i++) newey += pow(lambda[alpha][i]*ey[i],2);
      newey = sqrt(newey);
      grcombo->SetPoint(alpha,newx,newy);
      grcombo->SetPointError(alpha,newex,newey);
      // print the combination and check if the combined value is in between the measurements in the two channels or not
      cout << "bin " << alpha << ": " << y[alpha] << " +/- " << ey[alpha] << ", " << y[alpha+nbins] << " +/- " << ey[alpha+nbins] << " -> " << newy << " +/- " << newey << " (" << (newy-y[alpha])/(y[alpha+nbins]-y[alpha]) << ")" << endl;
      // print the weights
      cout << "weight: ";
      for (unsigned int i=0; i<2*nbins; i++) cout << lambda[alpha][i] << " ";
      cout << endl;
   }
}

void combine_blue(TGraphErrors *gr1_stat, TGraphErrors *gr1_syst, TGraphErrors *gr2_stat, TGraphErrors *gr2_syst, TGraphErrors *gr12_stat, TGraphErrors *gr12_syst)
{
   unsigned int nbins = gr1_stat->GetN();
   for (unsigned int i=0; i<nbins; i++)
   {
      double x1, y1, x2, y2;
      gr1_stat->GetPoint(i,x1,y1);
      gr2_stat->GetPoint(i,x2,y2);
      if (x1 != x2) cout << "Warning, in bin " << i << ": x1 = " << x1 << ", x2 = " << x2 << endl;
      double ex1 = gr1_stat->GetErrorX(i);
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

void combine_blue_nocor(const char* infilename, const char* outfilename)
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

   for (unsigned int i=0; i<names.size(); i++)
   {
      cout << "Combining " << names[i] << endl;
      TGraphErrors *gr1_stat = (TGraphErrors*) fin->Get(Form("g%s_exp_statonly_1",names[i].c_str()));
      TGraphErrors *gr1_syst = (TGraphErrors*) fin->Get(Form("g%s_exp_1",names[i].c_str()));
      TGraphErrors *gr2_stat = (TGraphErrors*) fin->Get(Form("g%s_exp_statonly_2",names[i].c_str()));
      TGraphErrors *gr2_syst = (TGraphErrors*) fin->Get(Form("g%s_exp_2",names[i].c_str()));
      unsigned int nbins = gr1_stat->GetN();
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

void combine_blue(const char* infilename, const char* outfilename, bool build_covmat_from_graphs)
{
   TFile *fin = new TFile(infilename);
   TFile *fout = new TFile(outfilename,"RECREATE");
   fout->mkdir("matrices");

   vector<string> names;
   names.push_back("yieldsp");
   names.push_back("yieldsm");
   names.push_back("yields");
   names.push_back("ch");
   names.push_back("A1p");
   names.push_back("A1m");
   names.push_back("A3");
   names.push_back("A4");

   for (unsigned int i=0; i<names.size(); i++)
   {
      cout << "Combining " << names[i] << endl;
      TGraphErrors *gr_stat1 = (TGraphErrors*) fin->Get(Form("g%s_exp_statonly_1",names[i].c_str()));
      TGraphErrors *gr_syst1 = (TGraphErrors*) fin->Get(Form("g%s_exp_1",names[i].c_str()));
      TGraphErrors *gr_stat2 = (TGraphErrors*) fin->Get(Form("g%s_exp_statonly_2",names[i].c_str()));
      TGraphErrors *gr_syst2 = (TGraphErrors*) fin->Get(Form("g%s_exp_2",names[i].c_str()));
      TGraphErrors *gr_stat = combine(gr_stat1,gr_stat2);
      TGraphErrors *gr_syst = combine(gr_syst1,gr_syst2);
      unsigned int nbins = gr_stat1->GetN();
      TGraphErrors *gr12_stat = new TGraphErrors(nbins);
      TGraphErrors *gr12_syst = new TGraphErrors(nbins);

      string label=names[i];
      if (label==string("yieldsp")) label="yp";
      else if (label==string("yieldsm")) label="ym";
      else if (label==string("yields")) label="yt";
      TH2F *hcov = (TH2F*) fin->Get(Form("matrices/hcov2_%s_dat_ch12",label.c_str()));
      if (!hcov || build_covmat_from_graphs) hcov = build_covmat(gr_stat,gr_syst,false);
      fout->cd("matrices");
      TH2F *hcovdat12 = new TH2F(Form("hcov2_%s_dat_ch1",label.c_str()),Form("hcov2_%s_dat_ch1",label.c_str()),nbins,0,nbins,nbins,0,nbins);

      fout->cd();
      combine_blue(gr_stat, hcov, gr12_stat, hcovdat12);
      combine_blue(gr_syst, hcov, gr12_syst, hcovdat12);
      fout->cd();
      gr12_stat->SetName(Form("g%s_exp_statonly_1",names[i].c_str())); gr12_stat->Write();
      gr12_syst->SetName(Form("g%s_exp_1",names[i].c_str())); gr12_syst->Write();
   }

   fout->Write();
   fout->Close();
   fin->Close();
}

double chi2(unsigned int nbins, double *y, const char* label, const char* input_dta, TH2F *hcov, bool donuc)
{
   ifstream file(input_dta);

   double dummy;
   double *A_CT10 = new double[nbins];
   double *A_EPS09 = new double[nbins];

   for (unsigned int ibin=0; ibin<nbins; ibin++)
   {
      int ibin2=ibin;
      if (string(label)==string("CA") || string(label)==string("yieldsWp") || string(label)==string("yieldsWm")) ibin2 = nbins-1-ibin;
      file >> dummy >> dummy >> A_CT10[ibin2] >> dummy >> dummy >> A_EPS09[ibin2] >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;

      if (string(label) == "yieldsWp" || string(label) == "yieldsWm")
      {
         double deta = 1.;
         // if (channel_number==1 && (ibin2==0 || ibin2==nbins-1)) deta = 0.4;
         A_CT10[ibin2] *= 208.*deta*LUMI*1e-6;
         A_EPS09[ibin2] *= 208.*deta*LUMI*1e-6;
      }
   }

   double *th = donuc ? A_EPS09 : A_CT10;

   TMatrixDSym tcov(nbins);
   for (unsigned int ibin=0; ibin<nbins; ibin++)
      for (unsigned int jbin=0; jbin<nbins; jbin++)
         tcov[ibin][jbin] = hcov->GetBinContent(ibin+1,jbin+1);
   TMatrixDSym tinvcov = tcov.Invert();

   // cout << label << ", th: ";
   // for (unsigned int ibin=0; ibin<nbins; ibin++)
   //    cout << th[ibin] << " ";
   // cout << endl;
   // cout << label << ", exp: ";
   // for (unsigned int ibin=0; ibin<nbins; ibin++)
   //    cout << y[ibin] << " ";
   // cout << endl;

   double chi2=0;
   for (unsigned int ibin=0; ibin<nbins; ibin++)
      for (unsigned int jbin=0; jbin<nbins; jbin++)
         chi2 += (y[ibin]-th[ibin])*tinvcov[ibin][jbin]*(y[jbin]-th[jbin]);

   return chi2;
}

double chi2(TGraphErrors *grstat, const char* label, const char* input_dta, TH2F *hcov, bool donuc)
{
   unsigned int nbins = grstat->GetN();
   double *y = grstat->GetY();
   return chi2(nbins, y, label, input_dta, hcov, donuc);
}

void chi2(const char* file_graphs, const char* file_matrices_exp, const char* file_matrices_th, int nchan)
{
   TFile *fgr = new TFile(file_graphs);
   TFile *fmat_exp = new TFile(file_matrices_exp);
   TFile *fmat_th = new TFile(file_matrices_th);

   TCanvas *c1 = new TCanvas();
   TCanvas *c2 = new TCanvas();
   TCanvas *c3 = new TCanvas();

   vector<string> names;
   names.push_back("yieldsp");
   names.push_back("yieldsm");
   names.push_back("ch");
   names.push_back("A1p");
   names.push_back("A1m");
   names.push_back("A3");
   names.push_back("A4");

   // print the chi2's in TeX format into a file
   ofstream texfile("table_chi2.tex");
   texfile << "\% table for channel 1" << endl;
   texfile << "\\begin{tabular}{|l|c|c||c|c|}" << endl;
   texfile << "\\hline" << endl;
   texfile << " & \\multicolumn{2}{c||}{No nuclear effects} & \\multicolumn{2}{c|}{Nuclear effects} \\\\" << endl;
   texfile << "Asymmetry & $\\chi^2$ & Prob.($\\chi^2$) & $\\chi^2$ & Prob.($\\chi^2$) \\\\" << endl;
   texfile << "\\hline" << endl;

   for (unsigned int i=0; i<names.size(); i++)
   {
      TGraphErrors *gr = (TGraphErrors*) fgr->Get(Form("g%s_exp_%i",names[i].c_str(),nchan));
      string label,label2;
      if (string(names[i])==string("A1p")) {label="A1p"; label2="A1p";}
      else if (string(names[i])==string("A1m")) {label="A1m"; label2="A1m";}
      else if (string(names[i])==string("ch")) {label="CA"; label2="ch";}
      else if (string(names[i])==string("A3")) {label="A3"; label2="A3";}
      else if (string(names[i])==string("A4")) {label="A4"; label2="A4";}
      else if (string(names[i])==string("yieldsp")) {label="yieldsWp"; label2="yp";}
      else if (string(names[i])==string("yieldsm")) {label="yieldsWm"; label2="ym";}
      TH2F *hcov_nuc_th = (TH2F*) fmat_th->Get(Form("matrices/hcov2_%s_nuc_ch%i",label2.c_str(),nchan));
      TH2F *hcov_nonuc_th = (TH2F*) fmat_th->Get(Form("matrices/hcov2_%s_nonuc_ch%i",label2.c_str(),nchan));
      TH2F *hcov_exp = (TH2F*) fmat_exp->Get(Form("matrices/hcov2_%s_dat_ch%i",label2.c_str(),nchan));
      TH2F *hcov_nuc_all = new TH2F(*hcov_nuc_th);
      hcov_nuc_all->Add(hcov_exp);
      TH2F *hcov_nonuc_all = new TH2F(*hcov_nonuc_th);
      hcov_nonuc_all->Add(hcov_exp);
      string input_dta = label + string("_muon.dta");
      double chi2_nuc = chi2(gr, label.c_str(), input_dta.c_str(), hcov_nuc_all, true);
      double chi2_nonuc = chi2(gr, label.c_str(), input_dta.c_str(), hcov_nonuc_all, false);
      int nbins = gr->GetN();
      cout << names[i] << ": nuclear effects -> " << chi2_nuc << " (prob. " << TMath::Prob(chi2_nuc,nbins) << "), no nuclear effects -> " << chi2_nonuc << " (prob. " << TMath::Prob(chi2_nonuc,nbins) << ")" << endl;
      texfile << label << " & " << chi2_nonuc << " & " << TMath::Prob(chi2_nonuc,nbins) << " & " << chi2_nuc << " & " << TMath::Prob(chi2_nuc,nbins) << "\\\\" << endl;

      if (i==5)
      {
         c1->cd();
         hcov_exp->Draw("COLZ TEXT");
         c2->cd();
         hcov_nuc_th->Draw("COLZ TEXT");
         c3->cd();
         hcov_nuc_all->Draw("COLZ TEXT");
      }
   }
   texfile << "\\hline" << endl;
   texfile << "\\end{tabular}" << endl << endl << endl;
   texfile.close();
}

void chi2(const char* file_graphs, const char* file_matrices_th, int nchan)
{
   TFile *fgr = new TFile(file_graphs);
   TFile *fmat_th = new TFile(file_matrices_th);

   TCanvas *c1 = new TCanvas();
   TCanvas *c2 = new TCanvas();
   TCanvas *c3 = new TCanvas();

   vector<string> names;
   names.push_back("yieldsp");
   names.push_back("yieldsm");
   names.push_back("ch");
   names.push_back("A1p");
   names.push_back("A1m");
   names.push_back("A3");
   names.push_back("A4");

   // print the chi2's in TeX format into a file
   ofstream texfile("table_chi2.tex");
   texfile << "\% table for channel 1" << endl;
   texfile << "\\begin{tabular}{|l|c|c||c|c|}" << endl;
   texfile << "\\hline" << endl;
   texfile << " & \\multicolumn{2}{c||}{No nuclear effects} & \\multicolumn{2}{c|}{Nuclear effects} \\\\" << endl;
   texfile << "Asymmetry & $\\chi^2$ & Prob.(\\%) & $\\chi^2$ & Prob.(\\%) \\\\" << endl;
   texfile << "\\hline" << endl;

   for (unsigned int i=0; i<names.size(); i++)
   {
      TGraphErrors *gr = (TGraphErrors*) fgr->Get(Form("g%s_exp_%i",names[i].c_str(),nchan));
      TGraphErrors *grstat = (TGraphErrors*) fgr->Get(Form("g%s_exp_statonly_%i",names[i].c_str(),nchan));
      string label,label2;
      if (string(names[i])==string("A1p")) {label="A1p"; label2="A1p";}
      else if (string(names[i])==string("A1m")) {label="A1m"; label2="A1m";}
      else if (string(names[i])==string("ch")) {label="CA"; label2="ch";}
      else if (string(names[i])==string("A3")) {label="A3"; label2="A3";}
      else if (string(names[i])==string("A4")) {label="A4"; label2="A4";}
      else if (string(names[i])==string("yieldsp")) {label="yieldsWp"; label2="yp";}
      else if (string(names[i])==string("yieldsm")) {label="yieldsWm"; label2="ym";}
      TH2F *hcov_nuc_th = (TH2F*) fmat_th->Get(Form("matrices/hcov2_%s_nuc_ch%i",label2.c_str(),nchan));
      TH2F *hcov_nonuc_th = (TH2F*) fmat_th->Get(Form("matrices/hcov2_%s_nonuc_ch%i",label2.c_str(),nchan));
      TH2F *hcov_exp = (TH2F*) build_covmat(grstat,gr,(label2=="yp"||label2=="ym")); // We want the luminosity uncertainty in the covariance matrix only for yields
      TH2F *hcov_nuc_all = new TH2F(*hcov_nuc_th);
      hcov_nuc_all->Add(hcov_exp);
      TH2F *hcov_nonuc_all = new TH2F(*hcov_nonuc_th);
      hcov_nonuc_all->Add(hcov_exp);
      string input_dta = label + string("_muon.dta");
      double chi2_nuc = chi2(gr, label.c_str(), input_dta.c_str(), hcov_nuc_all, true);
      double chi2_nonuc = chi2(gr, label.c_str(), input_dta.c_str(), hcov_nonuc_all, false);
      int nbins = gr->GetN();
      cout << names[i] << ": nuclear effects -> " << chi2_nuc << " (prob. " << TMath::Prob(chi2_nuc,nbins) << "), no nuclear effects -> " << chi2_nonuc << " (prob. " << TMath::Prob(chi2_nonuc,nbins) << ")" << endl;
      texfile << label << " & " << chi2_nonuc << " & " << 100.*TMath::Prob(chi2_nonuc,nbins) << " & " << chi2_nuc << " & " << 100.*TMath::Prob(chi2_nuc,nbins) << "\\\\" << endl;

      if (i==5)
      {
         c1->cd();
         hcov_exp->Draw("COLZ TEXT");
         c2->cd();
         hcov_nuc_th->Draw("COLZ TEXT");
         c3->cd();
         hcov_nuc_all->Draw("COLZ TEXT");
      }
   }
   texfile << "\\hline" << endl;
   texfile << "\\end{tabular}" << endl << endl << endl;
   texfile.close();
}

double combine_and_chi2(int nbins, double *y12, const char* label, const char* input_dta, TH2F *hcov_exp, TH2F *hcov_th, bool donuc)
{
   TGraphErrors *gr = new TGraphErrors(nbins,y12,y12);
   TGraphErrors *grcombo = new TGraphErrors(nbins/2);
   TH2F *hcov_exp_combo = new TH2F(Form("h%i",gRandom->Integer(pow(2,31)-1)),Form("h%i",gRandom->Integer(pow(2,31)-1)),nbins/2,0,nbins/2,nbins/2,0,nbins/2);
   combine_blue(gr,hcov_exp,grcombo,hcov_exp_combo);

   TH2F *hcov_expth = new TH2F(Form("h%i",gRandom->Integer(pow(2,31)-1)),Form("h%i",gRandom->Integer(pow(2,31)-1)),nbins/2,0,nbins/2,nbins/2,0,nbins/2);
   for (int ibin=1; ibin<nbins/2+1; ibin++)
   {
      for (int jbin=1; jbin<nbins/2+1; jbin++)
      {
         hcov_expth->SetBinContent(ibin,jbin,hcov_exp_combo->GetBinContent(ibin,jbin)+hcov_th->GetBinContent(ibin,jbin));
      }
   }

   double thechi2 = chi2(grcombo,label,input_dta,hcov_expth,donuc);
   delete gr;
   delete grcombo;
   delete hcov_exp_combo;
   delete hcov_expth;

   return thechi2;
}

TH1F* chi2_expected(TGraphErrors *gr_stat, TGraphErrors *gr_syst, TGraph *grth1, TGraph *grth2, TH2F *hcov_th1, TH2F *hcov_th2, TH2F *hcov_exp, int npe)
{
   int nbins = gr_stat->GetN();

   TH1F *hpe = new TH1F("hpe","hpe",npe/100,-5*nbins,5*nbins);
   TH2F *hpe2 = new TH2F("hpe2","hpe2",nbins,0,nbins,100,0,5000);

   double *y = gr_stat->GetY();
   double *dy_stat = gr_stat->GetEY();
   double *dy_syst = gr_syst->GetEY();
   double *y_th1 = grth1->GetY();
   double *y_th2 = grth2->GetY();

   TMatrixDSym tcov1(nbins), tcov2(nbins);
   for (unsigned int ibin=0; ibin<nbins; ibin++)
      for (unsigned int jbin=0; jbin<nbins; jbin++)
      {
         tcov1[ibin][jbin] = hcov_th1->GetBinContent(ibin+1,jbin+1)+hcov_exp->GetBinContent(ibin+1,jbin+1);
         tcov2[ibin][jbin] = hcov_th2->GetBinContent(ibin+1,jbin+1)+hcov_exp->GetBinContent(ibin+1,jbin+1);
         // tcov[ibin][jbin] = hcov_exp->GetBinContent(ibin+1,jbin+1);
      }
   TMatrixDSym tinvcov1 = tcov1.Invert();
   TMatrixDSym tinvcov2 = tcov2.Invert();

   for (int ipe=0; ipe<npe; ipe++)
   {
      double chi2_1=0;
      double chi2_2=0;

      double *ype = new double[nbins];
      for (unsigned int ibin=0; ibin<nbins; ibin++)
      {
         // ype[ibin] = y_th[ibin] * gRandom->Gaus(1.,dy_stat[ibin]/y[ibin]) * gRandom->Gaus(1.,dy_syst[ibin]/y[ibin]);
         ype[ibin] = gRandom->Gaus(y[ibin],sqrt(pow(dy_stat[ibin],2)+pow(dy_syst[ibin],2)));
         hpe2->Fill(ibin,ype[ibin]);
      }

      for (unsigned int ibin=0; ibin<nbins; ibin++)
         for (unsigned int jbin=0; jbin<nbins; jbin++)
         {
            chi2_1 += (y_th1[ibin]-ype[ibin])*tinvcov1[ibin][jbin]*(y_th1[jbin]-ype[jbin]);
            chi2_2 += (y_th2[ibin]-ype[ibin])*tinvcov2[ibin][jbin]*(y_th2[jbin]-ype[jbin]);
         }

      hpe->Fill(chi2_1-chi2_2);
   }

   double q[1], probSum[1];
   probSum[0] = 0.5;
   hpe->GetQuantiles(1,q,probSum);
   cout << "median: " << q[0] << endl;
   cout << "mean: " << hpe->GetMean() << endl;
   hpe2->Draw("COLZ");
   return hpe;
}

TH1F* chi2_expected(TGraphErrors *gr_stat, TGraphErrors *gr_syst, const char* label, TH2F *hcov_th1, TH2F *hcov_th2, TH2F *hcov_exp, bool donuc, int npe)
{
   string input_dta = string(label) + string("_muon.dta");
   ifstream file(input_dta.c_str());
   int nbins = gr_stat->GetN();

   double dummy;
   double *A_CT10 = new double[nbins];
   double *A_EPS09 = new double[nbins];

   for (unsigned int ibin=0; ibin<nbins; ibin++)
   {
      int ibin2=ibin;
      if (string(label)==string("CA") || string(label)==string("yieldsWm") || string(label)==string("yieldsWp")) 
      {
         ibin2 = nbins-1-ibin;
      }
      file >> dummy >> dummy >> A_CT10[ibin2] >> dummy >> dummy >> A_EPS09[ibin2] >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy;

      if (string(label) == "yieldsWp" || string(label) == "yieldsWm")
      {
         double deta = 1.;
         // if (channel_number==1 && (ibin2==0 || ibin2==nbins-1)) deta = 0.4;
         A_CT10[ibin2] *= 208.*deta*LUMI*1e-6;
         A_EPS09[ibin2] *= 208.*deta*LUMI*1e-6;
      }
   }

   // double *th = donuc ? A_EPS09 : A_CT10;
   TGraph *gth1 = new TGraph(nbins,A_EPS09,A_EPS09);
   TGraph *gth2 = new TGraph(nbins,A_CT10,A_CT10);

   for (unsigned int ibin=0; ibin<nbins; ibin++)
   {
      gr_stat->SetPoint(ibin,gr_stat->GetX()[ibin],donuc ? gth1->GetY()[ibin] : gth2->GetY()[ibin]);
   }

   return chi2_expected(gr_stat, gr_syst, gth1, gth2, hcov_th1, hcov_th2, hcov_exp, npe);
}

TH1F* chi2_expected(const char* file_exp, const char* file_th, const char* label, bool pred_nuc, int npe)
{
   TFile *fexp = new TFile(file_exp);
   TFile *fth = new TFile(file_th);

   TGraphErrors *gr_stat = (TGraphErrors*) fexp->Get(Form("g%s_exp_statonly_1",label));
   TGraphErrors *gr_syst = (TGraphErrors*) fexp->Get(Form("g%s_exp_1",label));

   string label_cov = string(label);
   if (label_cov=="yieldsp") label_cov = string("yp");
   else if (label_cov=="yieldsm") label_cov = string("ym");

   string label_dta = string(label);
   if (label_dta=="yieldsp") label_dta = string("yieldsWp");
   else if (label_dta=="yieldsm") label_dta = string("yieldsWm");
   else if (label_dta=="ch") label_dta = string("CA");

   // string label_covnuc = cov_nuc ? "nuc" : "nonuc";

   TH2F *hcov_th1 = (TH2F*) fth->Get(Form("matrices/hcov2_%s_nuc_ch1",label_cov.c_str()));
   TH2F *hcov_th2 = (TH2F*) fth->Get(Form("matrices/hcov2_%s_nonuc_ch1",label_cov.c_str()));
   TH2F *hcov_exp = (TH2F*) fexp->Get(Form("matrices/hcov2_%s_dat_ch1",label_cov.c_str()));

   return chi2_expected(gr_stat, gr_syst, label_dta.c_str(), hcov_th1, hcov_th2, hcov_exp, pred_nuc, npe);
}

#endif // ifndef combine_C
