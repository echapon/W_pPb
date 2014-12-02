#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include <fstream>
#include <string>
#include "lumierror.h"

EColor gMyColor = kBlack;
// EColor gMyColor = kBlue;
EColor gMyColor1 = kRed;
EColor gMyColor2 = kBlue;
EColor gColorCT10 = kRed;
EColor gColorCT10_fill = kYellow;
EColor gColorEPS09 = kGreen+2;
int gMyMarker = 22;
int gMyMarker1 = 20;
int gMyMarker2 = 21;

float gTextSize = 0.05;

const char *channeltext = "W #rightarrow #it{l} + #nu";
const char *channeltext_p = "W^{+} #rightarrow #it{l}^{+} + #nu";
const char *channeltext_m = "W^{-} #rightarrow #it{l}^{-} + #nu";
// const char *channeltext = "W #rightarrow e + #nu";
// const char *channeltext_p = "w^{+} #rightarrow e^{+} + #nu";
// const char *channeltext_m = "w^{-} #rightarrow e^{-} + #nu";
const char *channeltext1 = "W #rightarrow #mu + #nu";
const char *channeltext1_p = "W^{+} #rightarrow #mu^{+} + #nu";
const char *channeltext1_m = "W^{-} #rightarrow #mu^{-} + #nu";
// const char *channeltext1_p = "W^{+} #rightarrow e^{+} + #nu (iso)";
// const char *channeltext1_m = "W^{-} #rightarrow e^{-} + #nu (iso)";
const char *channeltext2 = "W #rightarrow e + #nu";
const char *channeltext2_p = "W^{+} #rightarrow e^{+} + #nu";
const char *channeltext2_m = "W^{-} #rightarrow e^{-} + #nu";
// const char *channeltext2_p = "W^{+} #rightarrow e^{+} + #nu (no iso)";
// const char *channeltext2_m = "W^{-} #rightarrow e^{-} + #nu (no iso)";

void revertX(TGraphErrors *gr)
{
   int nbins = gr->GetN();
   for (int i=0; i<nbins; i++)
   {
      double x = gr->GetX()[i];
      double y = gr->GetY()[i];
      gr->SetPoint(i,-x,y);
   }
}

void shiftX(TGraphErrors *gr, double dx)
{
   int nbins = gr->GetN();
   for (int i=0; i<nbins; i++)
   {
      double x = gr->GetX()[i];
      double y = gr->GetY()[i];
      gr->SetPoint(i,x+dx,y);
   }
}

void setErrorX(TGraphErrors *gr, double dx)
{
   int nbins = gr->GetN();
   for (int i=0; i<nbins; i++)
   {
      double dy = gr->GetEY()[i];
      gr->SetPointError(i,dx,dy);
   }
}

void setErrorY(TGraphAsymmErrors *gr, double dy)
{
   int nbins = gr->GetN();
   for (int i=0; i<nbins; i++)
   {
      double dx = gr->GetEXhigh()[i];
      gr->SetPointError(i,dx,dx,dy,dy);
   }
}

void revertX(TGraphAsymmErrors *gr)
{
   int nbins = gr->GetN();
   for (int i=0; i<nbins; i++)
   {
      double x = gr->GetX()[i];
      double y = gr->GetY()[i];
      gr->SetPoint(i,-x,y);
   }
}

void scaleY(TGraphErrors *gr, double scale)
{
   int nbins = gr->GetN();
   for (int i=0; i<nbins; i++)
   {
      double x = gr->GetX()[i];
      double y = gr->GetY()[i];
      double ex = gr->GetEX()[i];
      double ey = gr->GetEY()[i];
      gr->SetPoint(i,x,y*scale);
      gr->SetPointError(i,ex,ey*scale);
   }
}

void scaleY(TGraphAsymmErrors *gr, double scale)
{
   int nbins = gr->GetN();
   for (int i=0; i<nbins; i++)
   {
      double x = gr->GetX()[i];
      double y = gr->GetY()[i];
      double exl = gr->GetEXlow()[i];
      double eyl = gr->GetEYlow()[i];
      double exh = gr->GetEXhigh()[i];
      double eyh = gr->GetEYhigh()[i];
      gr->SetPoint(i,x,y*scale);
      gr->SetPointError(i,exl,exh,eyl*scale,eyh*scale);
   }
}

TGraphAsymmErrors* theory_errors(const char* gbasename, int channel_number, int nbins, const char* pdfname)
{
   string label;
   if (string(gbasename)==string("gA1p")) label="A1p";
   else if (string(gbasename)==string("gA1m")) label="A1m";
   else if (string(gbasename)==string("gch")) label="CA";
   else if (string(gbasename)==string("gA3")) label="A3";
   else if (string(gbasename)==string("gA4")) label="A4";
   else if (string(gbasename)==string("gyieldsp")) label="yieldsWp";
   else if (string(gbasename)==string("gyieldsm")) label="yieldsWm";

   string channel = (channel_number==1) ? "muon" : "elect";
   // string channel = (channel_number==1) ? "wide" : "elect";

   string filename = string(label) + string("_") + string(channel) + string(".dta");
   ifstream file(filename.c_str());

   double dummy;
   double *y_median = new double[nbins];
   double *A_CT10 = new double[nbins];
   double *errA_CT10_CT10_up = new double[nbins];
   double *errA_CT10_CT10_down = new double[nbins];
   double *A_EPS09 = new double[nbins];
   double *errA_EPS09_CT10_up = new double[nbins];
   double *errA_EPS09_CT10_down = new double[nbins];
   double *errA_EPS09_EPS09_up = new double[nbins];
   double *errA_EPS09_EPS09_down = new double[nbins];
   double *errA_EPS09_scales_up = new double[nbins];
   double *errA_EPS09_scales_down = new double[nbins];

   for (int ibin=0; ibin<nbins; ibin++)
   {
      int ibin2=ibin;
      if (label=="CA" || label=="yieldsWp" || label=="yieldsWm") ibin2 = nbins-1-ibin;
      file >> dummy >> y_median[ibin] >> A_CT10[ibin2] >> errA_CT10_CT10_up[ibin2] >> errA_CT10_CT10_down[ibin2] >> A_EPS09[ibin2] >> errA_EPS09_CT10_up[ibin2] >> errA_EPS09_CT10_down[ibin2] >> errA_EPS09_EPS09_up[ibin2] >> errA_EPS09_EPS09_down[ibin2] >> errA_EPS09_scales_up[ibin2] >> errA_EPS09_scales_down[ibin2];

      if (label == "yieldsWp" || label == "yieldsWm")
      {
         double deta = 1.;
         // if (channel_number==1 && (ibin2==0 || ibin2==nbins-1)) deta = 0.4;
         A_CT10[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_CT10_CT10_up[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_CT10_CT10_down[ibin2] *= 208.*deta*LUMI*1e-6;
         A_EPS09[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_EPS09_CT10_up[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_EPS09_CT10_down[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_EPS09_EPS09_up[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_EPS09_EPS09_down[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_EPS09_scales_up[ibin2] *= 208.*deta*LUMI*1e-6;
         errA_EPS09_scales_down[ibin2] *= 208.*deta*LUMI*1e-6;
      }
   }

   double *erry = new double[nbins];
   for (int ibin=0; ibin<nbins; ibin++) erry[ibin] = 0.25;
   if (channel_number==1) 
   {
      if (label=="yieldsWp" || label=="yieldsWm" || label=="CA") erry[0]=0.2; 
      erry[nbins-1]=0.2;
   }

   double *A = new double[nbins];
   double *dA_up = new double[nbins];
   double *dA_down = new double[nbins];

   if (string(pdfname) == "CT10")
   {
      for (int ibin=0; ibin<nbins; ibin++)
      {
         A[ibin] = A_CT10[ibin];
         dA_up[ibin] = sqrt(pow(errA_CT10_CT10_up[ibin],2)+pow(errA_EPS09_scales_up[ibin],2));
         dA_down[ibin] = sqrt(pow(errA_CT10_CT10_down[ibin],2)+pow(errA_EPS09_scales_down[ibin],2));
      }
   }
   else
   {
      for (int ibin=0; ibin<nbins; ibin++)
      {
         A[ibin] = A_EPS09[ibin];
         dA_up[ibin] = sqrt(pow(errA_EPS09_EPS09_up[ibin],2)+pow(errA_EPS09_scales_up[ibin],2)+pow(errA_EPS09_CT10_up[ibin],2));
         dA_down[ibin] = sqrt(pow(errA_EPS09_EPS09_down[ibin],2)+pow(errA_EPS09_scales_down[ibin],2)+pow(errA_EPS09_CT10_down[ibin],2));
      }
   }

   for (int ibin=0; ibin<nbins; ibin++)
      cout << ibin << " " << A[ibin] << " +" << dA_up[ibin] << " -" << dA_down[ibin] << endl;

   TGraphAsymmErrors *gr = new TGraphAsymmErrors(nbins,y_median,A,erry,erry,dA_down,dA_up);

   return gr;
}

void plot_graph(const char* fname_cteq="graph.root", const char* fname_eps="graph.root", const char *gbasename="gyields", int channel_number=1, const char* xlabel="#eta_{lab}", const char* ylabel="asymetry", double xmin=-2.5, double xmax=2.5, double ymin=0, double ymax=-1)
{
   gStyle->SetOptTitle(0);
   TFile *fcteq = new TFile(fname_cteq);
   TFile *feps = new TFile(fname_eps);

   TString gname_stat = Form("%s_exp_statonly_%i",gbasename,channel_number);
   TString gname_syst = Form("%s_exp_%i",gbasename,channel_number);

   TGraphErrors *gexp_stat = (TGraphErrors*) fcteq->Get(gname_stat);
   TGraphErrors *gexp_syst = (TGraphErrors*) fcteq->Get(gname_syst);
   int nbins = gexp_stat->GetN();
   TGraphAsymmErrors *gth_cteq = theory_errors(gbasename, channel_number, nbins, "CT10");
   TGraphAsymmErrors *gth_eps = theory_errors(gbasename, channel_number, nbins, "EPS09");

   gth_cteq->SetLineColor(gColorCT10);
   gth_cteq->SetLineWidth(2);
   gth_cteq->SetFillColor(gColorCT10);
   gth_cteq->SetFillStyle(3002);
   gth_cteq->Draw("AL3");
   if (xmin<xmax) gth_cteq->GetXaxis()->SetRangeUser(xmin,xmax);
   if (ymin<ymax) gth_cteq->GetYaxis()->SetRangeUser(xmin,xmax);
   gth_cteq->GetXaxis()->SetTitle(xlabel);
   gth_cteq->GetYaxis()->SetTitle(ylabel);
   gth_cteq->Draw("AL3");

   gth_eps->SetLineColor(gColorEPS09);
   gth_eps->SetLineStyle(9);
   gth_eps->SetLineWidth(2);
   gth_eps->SetFillColor(gColorEPS09);
   gth_eps->SetFillStyle(3003);
   gth_eps->Draw("L3");

   gexp_syst->SetLineColor(gMyColor);
   gexp_syst->SetMarkerColor(gMyColor);
   gexp_syst->SetMarkerStyle(8);
   gexp_syst->SetFillStyle(0);
   gexp_syst->SetLineWidth(2);
   gexp_syst->Draw("P5");

   gexp_stat->SetLineColor(gMyColor);
   gexp_stat->SetMarkerColor(gMyColor);
   gexp_stat->SetMarkerStyle(8);
   gexp_stat->SetLineWidth(2);
   gexp_stat->Draw("PZ");

   TH1F *dummy = new TH1F();
   TLegend *tleg = new TLegend(0.34,0.2,0.66,0.4);
   tleg->AddEntry(dummy,"CMS Preliminary","");
   tleg->AddEntry(dummy,"pPb data L=34.6 nb^{-1}","");
   tleg->AddEntry(gexp_syst,channeltext,"p");
   tleg->AddEntry(gth_cteq,"CT10","lf");
   tleg->AddEntry(gth_eps,"EPS09","lf");
   tleg->SetFillColor(kWhite);
   tleg->SetBorderSize(0);
   tleg->Draw();
}

void plot_graph_2chan(const char* fname_cteq1="graph_cteq_1.root", const char* fname_eps1="graph_eps_1.root", const char* fname_2="graph_2.root", const char *gbasename="gyields", const char* xlabel="#eta_{lab}", const char* ylabel="asymetry", double xmin=-2.5, double xmax=2.5, double ymin=0, double ymax=-1)
{
   gStyle->SetOptTitle(0);
   gStyle->SetEndErrorSize(10);
   TFile *fcteq1 = new TFile(fname_cteq1);
   TFile *feps1 = new TFile(fname_eps1);
   TFile *f2 = new TFile(fname_2);

   int channel_number = 1;

   TString gname_th = Form("%s_th_%i",gbasename,channel_number);
   TString gname_stat1 = Form("%s_exp_statonly_%i",gbasename,channel_number);
   TString gname_syst1 = Form("%s_exp_%i",gbasename,channel_number);
   channel_number=1;
   TString gname_stat2 = Form("%s_exp_statonly_%i",gbasename,channel_number);
   TString gname_syst2 = Form("%s_exp_%i",gbasename,channel_number);

   TGraphErrors *gexp_stat1 = (TGraphErrors*) fcteq1->Get(gname_stat1);
   TGraphErrors *gexp_syst1 = (TGraphErrors*) fcteq1->Get(gname_syst1);
   TGraphErrors *gexp_stat2 = (TGraphErrors*) f2->Get(gname_stat2);
   TGraphErrors *gexp_syst2 = (TGraphErrors*) f2->Get(gname_syst2);
   int nbins = gexp_stat1->GetN();
   TGraphErrors *gth_cteq = (TGraphErrors*) fcteq1->Get(gname_th);
   TGraphErrors *gth_eps = (TGraphErrors*) feps1->Get(gname_th);

   gth_cteq->SetLineColor(gColorCT10);
   gth_cteq->SetLineWidth(2);
   gth_cteq->SetFillColor(gColorCT10);
   gth_cteq->SetFillStyle(3002);
   gth_cteq->Draw("AL3");
   if (xmin<xmax) gth_cteq->GetXaxis()->SetRangeUser(xmin,xmax);
   if (ymin<ymax) gth_cteq->GetYaxis()->SetRangeUser(ymin,ymax);
   gth_cteq->GetXaxis()->SetTitle(xlabel);
   gth_cteq->GetYaxis()->SetTitle(ylabel);
   gth_cteq->Draw("AL3");

   gth_eps->SetLineColor(gColorEPS09);
   gth_eps->SetLineStyle(9);
   gth_eps->SetLineWidth(2);
   gth_eps->SetFillColor(gColorPES09);
   gth_eps->SetFillStyle(3003);
   gth_eps->Draw("L3");

   gexp_syst1->SetLineColor(gMyColor1);
   gexp_syst1->SetMarkerColor(gMyColor1);
   gexp_syst1->SetMarkerStyle(8);
   gexp_syst1->SetFillStyle(0);
   gexp_syst1->SetLineWidth(2);
   gexp_syst1->Draw("P5");

   gexp_stat1->SetLineColor(gMyColor1);
   gexp_stat1->SetMarkerColor(gMyColor1);
   gexp_stat1->SetMarkerStyle(8);
   gexp_stat1->SetLineWidth(2);
   gexp_stat1->Draw("PZ");

   if (twochan)
   {
      gexp_syst2->SetLineColor(gMyColor2);
      gexp_syst2->SetMarkerColor(gMyColor2);
      gexp_syst2->SetMarkerStyle(8);
      gexp_syst2->SetFillStyle(0);
      gexp_syst2->SetLineWidth(2);
      gexp_syst2->Draw("P5");

      gexp_stat2->SetLineColor(gMyColor2);
      gexp_stat2->SetMarkerColor(gMyColor2);
      gexp_stat2->SetMarkerStyle(8);
      gexp_stat2->SetLineWidth(2);
      gexp_stat2->Draw("PZ");
   }

   TH1F *dummy = new TH1F();
   TLegend *tleg = new TLegend(0.34,0.2,0.66,0.4);
   tleg->AddEntry(dummy,"CMS Preliminary","");
   tleg->AddEntry(dummy,"pPb data L=34.6 nb^{-1}","");
   tleg->AddEntry(gexp_syst1,channeltext1,"p");
   tleg->AddEntry(gexp_syst2,channeltext2,"p");
   tleg->AddEntry(gth_cteq,"CT10","lf");
   tleg->AddEntry(gth_eps,"EPS09","lf");
   tleg->SetFillColor(kWhite);
   tleg->SetBorderSize(0);
   tleg->Draw();
}

void plot_graph_1file(const char* fname="graph.root", const char *gbasename="gyields", int channel_number=0, bool reverteta=false, double scale=1.)
{
   gStyle->SetOptTitle(0);
   gStyle->SetEndErrorSize(6);
   gStyle->SetHatchesLineWidth(2);
   gStyle->SetHatchesSpacing(1.7);
   TFile *f = new TFile(fname);

   TCanvas *c1 = new TCanvas("c1","c1",600,600);

   string xlabel, ylabel;
   double xmin, xmax, ymin, ymax;

   // x-axis
   xlabel = string("#eta_{lab}");
   if (string(gbasename)=="gyieldsp" || string(gbasename)=="gyieldsm" || string(gbasename)=="gch")
   {
      xmin = -4.1;
      xmax = 4.1;
   }
   else
   {
      xmin = 0.;
      xmax = 2.5;
   }

   bool twochan = (channel_number==0 || channel_number==3);
   if (twochan) channel_number = 1;

   // y-axis range
   if (string(gbasename)=="gyieldsp" || string(gbasename)=="gyieldsm")
   {
      ymin = 0.;
      ymax = 140.;
      if (string(gbasename)=="gyieldsp") ylabel = string("d#sigma (W^{+}#rightarrow#it{l}^{+}#nu) / d#eta_{lab} [nb]");
      else ylabel = string("d#sigma (W^{-}#rightarrow#it{l}^{-}#nu) / d#eta_{lab} [nb]");
   }
   else if (string(gbasename)=="gch")
   {
      ymin = -0.4;
      ymax = 0.4;
      ylabel = string("(N^{+}-N^{-})/(N^{+}+N^{-})");
   }
   else if (string(gbasename)=="gA1m")
   {
      ymin = 0.6;
      ymax = 1.4;
      ylabel = string("N^{-}(+#eta_{lab})/N^{-}(-#eta_{lab})");
   }
   else if (string(gbasename)=="gA1p")
   {
      ymin = 0.8;
      ymax = 3.2;
      ylabel = string("N^{+}(+#eta_{lab})/N^{+}(-#eta_{lab})");
   }
   else if (string(gbasename)=="gA3")
   {
      ymin = 0.8;
      ymax = (!twochan && channel_number==1) ? 1.8 : 1.9;
      ylabel = string("N(+#eta_{lab})/N(-#eta_{lab})");
   }
   else if (string(gbasename)=="gA4")
   {
      ymin = -0.01;
      ymax = 0.05;
      ylabel = string("A_{4}");
   }

   TString gname_stat1 = Form("%s_exp_statonly_%i",gbasename,channel_number);
   TString gname_syst1 = Form("%s_exp_%i",gbasename,channel_number);
   TString gname_stat2 = Form("%s_exp_statonly_%i",gbasename,2);
   TString gname_syst2 = Form("%s_exp_%i",gbasename,2);

   TGraphErrors *gexp_stat1 = (TGraphErrors*) f->Get(gname_stat1);
   TGraphErrors *gexp_syst1 = (TGraphErrors*) f->Get(gname_syst1);
   TGraphErrors *gexp_stat2 = (TGraphErrors*) f->Get(gname_stat2);
   TGraphErrors *gexp_syst2 = (TGraphErrors*) f->Get(gname_syst2);
   // TFile *fc = new TFile("combo_muacc_20140219.root");
   int nbins = gexp_stat1->GetN();
   // TGraphAsymmErrors *gth_cteq = theory_errors(gbasename, twochan ? 1 : channel_number, nbins, "CT10");
   // TGraphAsymmErrors *gth_eps = theory_errors(gbasename, twochan ? 1 : channel_number, nbins, "EPS09");
   TGraphAsymmErrors *gth_cteq = theory_errors(gbasename, twochan ? 1 : channel_number, 12, "CT10");
   TGraphAsymmErrors *gth_eps = theory_errors(gbasename, twochan ? 1 : channel_number, 12, "EPS09");

   if (reverteta)
   {
      revertX(gexp_stat1);
      revertX(gexp_syst1);
      if (twochan)
      {
         revertX(gexp_stat2);
         revertX(gexp_syst2);
      }
      revertX(gth_cteq);
      revertX(gth_eps);
   }

   scaleY(gexp_stat1,scale);
   scaleY(gexp_syst1,scale);
   setErrorX(gexp_stat1,0);
   setErrorX(gexp_syst1,0);
   gexp_stat1->SetMarkerSize(1.2);
   gexp_syst1->SetMarkerSize(1.2);
   if (twochan)
   {
      scaleY(gexp_stat2,scale);
      scaleY(gexp_syst2,scale);
      setErrorX(gexp_stat2,0);
      setErrorX(gexp_syst2,0);
      gexp_stat2->SetMarkerSize(1.2);
      gexp_syst2->SetMarkerSize(1.2);
      shiftX(gexp_stat1,-0.05);
      shiftX(gexp_syst1,-0.05);
      shiftX(gexp_stat2,0.05);
      shiftX(gexp_syst2,0.05);
   }
   else
   {
      scaleY(gth_cteq,scale);
      scaleY(gth_eps,scale);
   }

   gth_cteq->SetLineColor(gColorCT10);
   gth_cteq->SetMarkerSize(0);
   gth_cteq->SetLineWidth(2);
   gth_cteq->SetFillColor(gColorCT10_fill);
   // gth_cteq->SetFillStyle(3002);
   gth_cteq->SetFillStyle(1001);
   // gth_cteq->Draw("AL3");
   gth_cteq->Draw("A2");
   gth_cteq->GetXaxis()->SetLimits(xmin,xmax);
   gth_cteq->GetYaxis()->SetRangeUser(ymin,ymax);
   gth_cteq->GetXaxis()->SetTitle(xlabel.c_str());
   gth_cteq->GetXaxis()->SetTitleSize(gTextSize);
   gth_cteq->GetXaxis()->SetLabelSize(gTextSize);
   gth_cteq->GetYaxis()->SetTitle(ylabel.c_str());
   gth_cteq->GetYaxis()->SetTitleSize(gTextSize);
   gth_cteq->GetYaxis()->SetLabelSize(gTextSize);
   // gth_cteq->Draw("AL3");
   gth_cteq->Draw("A2");
   TGraphAsymmErrors *gth_cteq2 = new TGraphAsymmErrors(*gth_cteq);
   setErrorY(gth_cteq2,0);
   gth_cteq2->Draw("Z");

   gth_eps->SetLineColor(gColorEPS09);
   gth_eps->SetLineStyle(9);
   gth_eps->SetLineWidth(2);
   gth_eps->SetFillColor(gColorEPS09);
   gth_eps->SetFillStyle(3345);
   // gth_eps->Draw("L3");
   gth_eps->Draw("2");
   TGraphAsymmErrors *gth_eps2 = new TGraphAsymmErrors(*gth_eps);
   setErrorY(gth_eps2,0);
   gth_eps2->Draw("Z");

   gexp_syst1->SetLineColor((twochan) ? gMyColor1 : gMyColor);
   gexp_syst1->SetMarkerColor((twochan) ? gMyColor1 : gMyColor);
   gexp_syst1->SetMarkerStyle(gMyMarker1);
   gexp_syst1->SetFillStyle(0);
   gexp_syst1->SetLineWidth(2);
   gexp_syst1->Draw("||");

   gexp_stat1->SetLineColor((twochan) ? gMyColor1 : gMyColor);
   gexp_stat1->SetMarkerColor((twochan) ? gMyColor1 : gMyColor);
   gexp_stat1->SetMarkerStyle(gMyMarker1);
   gexp_stat1->SetLineWidth(2);
   gexp_stat1->Draw("PZ");

   if (twochan)
   {
      gexp_syst2->SetLineColor((twochan) ? gMyColor2 : gMyColor);
      gexp_syst2->SetMarkerColor((twochan) ? gMyColor2 : gMyColor);
      gexp_syst2->SetMarkerStyle(gMyMarker2);
      gexp_syst2->SetFillStyle(0);
      gexp_syst2->SetLineWidth(2);
      gexp_syst2->Draw("||");

      gexp_stat2->SetLineColor((twochan) ? gMyColor2 : gMyColor);
      gexp_stat2->SetMarkerColor((twochan) ? gMyColor2 : gMyColor);
      gexp_stat2->SetMarkerStyle(gMyMarker2);
      gexp_stat2->SetLineWidth(2);
      gexp_stat2->Draw("PZ");
   }

   TH1F *dummy = new TH1F();
   double xl=0.5,yl=0.18,dx=0.3,dy=0.28;
   if (string(gbasename)=="gA1m" || string(gbasename)=="gA1p" || string(gbasename)=="gA3" || string(gbasename)=="gA4")
   {
      xl=0.21;
      yl=0.62;
   }
   if (string(gbasename)=="gyieldsm") xl=0.4;
   if (string(gbasename)=="gA1m") xl=0.5;
   if ((string(gbasename)=="gch"||string(gbasename)=="gyieldsp")&&!reverteta) xl = 0.2;
   TLegend *tleg = new TLegend(xl,yl,xl+dx,yl+dy);
   tleg->SetTextFont(42);
   tleg->SetTextSize(gTextSize*0.8);
   tleg->AddEntry(dummy,"CMS Preliminary","");
   tleg->AddEntry(dummy,"pPb #sqrt{s_{NN}} = 5.02 TeV","");
   tleg->AddEntry(dummy,"L = 34.6 nb^{-1}","");
   string legend1 = twochan ? string(channeltext1) : string("Data"); //string(channeltext);
   if (string(gbasename)=="gyieldsp" || string(gbasename)=="gA1p") legend1 = twochan ? string(channeltext1_p) : string("Data"); //string(channeltext_p);
   else if (string(gbasename)=="gyieldsm" || string(gbasename)=="gA1m") legend1 = twochan ? string(channeltext1_m) : string("Data"); //string(channeltext_m);
   tleg->AddEntry(gexp_syst1,legend1.c_str(),"p");
   if (twochan) 
   {
      string legend2 = string(channeltext2);
      if (string(gbasename)=="gyieldsp" || string(gbasename)=="gA1p") legend2 = string(channeltext2_p);
      else if (string(gbasename)=="gyieldsm" || string(gbasename)=="gA1m") legend2 = string(channeltext2_m);
      tleg->AddEntry(gexp_syst2,legend2.c_str(),"p");
   }
   else
   {
      tleg->AddEntry(gth_cteq,"CT10","lf");
      tleg->AddEntry(gth_eps,"EPS09","lf");
   }
   tleg->SetFillColor(kWhite);
   tleg->SetBorderSize(0);
   tleg->Draw();

   if (string(gbasename)=="gyieldsp" || string(gbasename)=="gyieldsm")
   {
      double txl=-1.8, tyl=115., tdx=2., tdy=20.;
      if (!reverteta) txl = 0.3;
      // if (string(gbasename)=="gyieldsp") txl = 0.4;
      TPaveText *tp = new TPaveText(txl,tyl,txl+tdx,tyl+tdy);
      tp->AddText("Luminosity uncertainty: 3.5%");
      tp->AddText("p_{T}^{l} > 25 GeV/c");
      tp->SetFillColor(kWhite);
      tp->SetBorderSize(0);
      tp->SetTextFont(42);
      tp->SetTextSize(gTextSize*0.8);
      tp->Draw();
   }

   if (string(gbasename)=="gyieldsp" || string(gbasename)=="gyieldsm" || string(gbasename)=="gA1p" || string(gbasename)=="gA1m" || string(gbasename)=="gA3")
   {
      double txl=0.38, tyl=60, tdx=1.,tdy=20.;
      if (string(gbasename)=="gyieldsm") tyl=95;
      if (!reverteta) txl = -1.8;
      if (string(gbasename)=="gA1p") {txl = 0.38; tdx = 0.5; tyl = 1.75; tdy = 0.2;}
      if (string(gbasename)=="gA1m") {txl = 0.25; tdx = 0.5; tyl = 1.2; tdy = 0.1;}
      if (string(gbasename)=="gA3") {txl = 0.25; tdx = 0.5; tyl = 1.2; tdy = 0.1;}
      TPaveText *tp = new TPaveText(txl,tyl,txl+tdx,tyl+tdy);
      string legend;
      if (string(gbasename)=="gyieldsp" || string(gbasename)=="gA1p") legend = string(channeltext_p);
      else if (string(gbasename)=="gyieldsm" || string(gbasename)=="gA1m") legend = string(channeltext_m);
      else legend = string(channeltext);
      tp->AddText(legend.c_str());
      tp->SetFillColor(kWhite);
      tp->SetBorderSize(0);
      tp->SetTextFont(42);
      tp->SetTextSize(gTextSize);
      tp->SetTextAlign(13);
      tp->Draw();
   }

   gPad->SaveAs(Form("%s.pdf",gbasename));
   gPad->SaveAs(Form("%s.png",gbasename));
   gPad->SaveAs(Form("%s.root",gbasename));
   gPad->SaveAs(Form("%s.C",gbasename));
}

void plot_graph_1file_all(const char* file, int chan, bool reverse_eta=true, double scale=1./34.62)
{
   plot_graph_1file(file,"gA3",chan);
   plot_graph_1file(file,"gch",chan,reverse_eta);
   plot_graph_1file(file,"gA1p",chan);
   plot_graph_1file(file,"gA1m",chan);
   plot_graph_1file(file,"gA4",chan);
   plot_graph_1file(file,"gyieldsp",chan,reverse_eta,scale);
   plot_graph_1file(file,"gyieldsm",chan,reverse_eta,scale);
}
