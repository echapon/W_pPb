#include "TFile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TColor.h"
#include <fstream>
#include <iostream>
#include <string>
#include "lumierror.h"

#include "CMS_lumi.C"
writeExtraText = false;       // if extra text
int iPeriod = 99;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 99 = pPb

using namespace std;

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

bool doAlice = false;

float gTextSize = 0.05;

const char *channeltext = "W #rightarrow #font[12]{l} + #nu";
const char *channeltext_p = "W^{+} #rightarrow #font[12]{l}^{+} + #nu";
const char *channeltext_m = "W^{-} #rightarrow #font[12]{l}^{-} + #nu";
// const char *channeltext = "W #rightarrow e + #nu";
// const char *channeltext_p = "w^{+} #rightarrow e^{+} + #nu";
// const char *channeltext_m = "w^{-} #rightarrow e^{-} + #nu";
const char *channeltext1 = "W #rightarrow #mu + #nu";
const char *channeltext1_p = "W^{+} #rightarrow #mu^{+} + #nu";
const char *channeltext1_m = "W^{-} #rightarrow #mu^{-} + #nu";
// const char *channeltext1 = "W #rightarrow  + #font[12]{l}#nu (new)";
// const char *channeltext1_p = "W^{+} #rightarrow #font[12]{l}^{+} + #nu (new)";
// const char *channeltext1_m = "W^{-} #rightarrow #font[12]{l}^{-} + #nu (new)";
const char *channeltext2 = "W #rightarrow e + #nu";
const char *channeltext2_p = "W^{+} #rightarrow e^{+} + #nu";
const char *channeltext2_m = "W^{-} #rightarrow e^{-} + #nu";
// const char *channeltext2 = "W #rightarrow  + #font[12]{l}#nu (old)";
// const char *channeltext2_p = "W^{+} #rightarrow #font[12]{l}^{+} + #nu (old)";
// const char *channeltext2_m = "W^{-} #rightarrow #font[12]{l}^{-} + #nu (old)";
const char *channeltext3 = "W #rightarrow #mu + #nu (40-100%)";
const char *channeltext3_p = "W^{+} #rightarrow #mu^{+} + #nu (40-100%)";
const char *channeltext3_m = "W^{-} #rightarrow #mu^{-} + #nu (40-100%)";

double acc_alice_plus[2] = {1.2285,1.5920};
double acc_alice_plus_err[2] = {0.0035,0.0174};
double acc_alice_minus[2] = {1.4495,1.7864};
double acc_alice_minus_err[2] = {0.0063,0.0109};

void toto(const char* tut)
{
   cout << tut << endl;
}

TGraphErrors* divideGraphs(TGraphErrors *gr1, TGraphErrors *gr2)
{
   int nbins = gr1->GetN();

   double *x1 = gr1->GetX();
   double *ex1 = gr2->GetEX();
   double *y1 = gr1->GetY();
   double *ey1 = gr1->GetEY();
   double *y2 = gr2->GetY();
   double *ey2 = gr2->GetEY();

   double *y = new double[nbins];
   double *ey = new double[nbins];

   for (int i=0; i<nbins; i++)
   {
      y[i] = y1[i]/y2[i];
      // set error of denominator to 0
      ey2[i] = 0;

      ey[i] = sqrt(pow(ey1[i],2)+pow(ey2[i]*y1[i]/y2[i],2))/y2[i];
   }

   TGraphErrors *gr = new TGraphErrors(nbins,x1,y,ex1,ey);
   return gr;
}

TGraphErrors* convert(TGraphAsymmErrors* gr)
{
   int n = gr->GetN();
   double *x = gr->GetX();
   double *y = gr->GetY();
   double *ex = new double[n];
   double *ey = new double[n];
   for (int i=0; i<n; i++)
   {
      ex[i] = gr->GetErrorX(i);
      ey[i] = gr->GetErrorY(i);
   }

   TGraphErrors *grout = new TGraphErrors(n,x,y,ex,ey);
   return grout;
}

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

void setErrorY(TGraphErrors *gr, double dy)
{
   int nbins = gr->GetN();
   for (int i=0; i<nbins; i++)
   {
      double dx = gr->GetEX()[i];
      gr->SetPointError(i,dx,dy);
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

TGraphAsymmErrors* theory_errors(const char* gbasename, int channel_number, int nbins, const char* pdfname, bool doct10err=true)
{
   string label;
   if (string(gbasename)==string("gA1p")) label="A1p";
   else if (string(gbasename)==string("gA1m")) label="A1m";
   else if (string(gbasename)==string("gch")) label="CA";
   else if (string(gbasename)==string("gA3")) label="A3";
   else if (string(gbasename)==string("gA4")) label="A4";
   else if (string(gbasename)==string("gyieldsp")) label="yieldsWp";
   else if (string(gbasename)==string("gyieldsm")) label="yieldsWm";

   string channel = (channel_number==1 || !(nbins=10||nbins==5)) ? "muon" : "elect";

   string ext = (nbins==10||nbins==5) ? string(".dta") : string("_4bins.dta");

   string filename = string(label) + string("_") + string(channel) + ext;
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
   for (int ibin=0; ibin<nbins; ibin++) erry[ibin] = (nbins==10||nbins==5) ? 0.25 : 0.5;
   if (channel_number==1) 
   {
      if (label=="yieldsWp" || label=="yieldsWm" || label=="CA") erry[0] = (nbins==10||nbins==5) ? 0.2 : 0.7; 
      erry[nbins-1] = (nbins==10||nbins==5) ? 0.2 : 0.7;
   }

   double *A = new double[nbins];
   double *dA_up = new double[nbins];
   double *dA_down = new double[nbins];

   if (string(pdfname) == "CT10")
   {
      for (int ibin=0; ibin<nbins; ibin++)
      {
         A[ibin] = A_CT10[ibin];
         dA_up[ibin] = doct10err ? sqrt(pow(errA_CT10_CT10_up[ibin],2)+pow(errA_EPS09_scales_up[ibin],2)) : 0;
         dA_down[ibin] = doct10err ? sqrt(pow(errA_CT10_CT10_down[ibin],2)+pow(errA_EPS09_scales_down[ibin],2)) : 0;
      }
   }
   else
   {
      for (int ibin=0; ibin<nbins; ibin++)
      {
         A[ibin] = A_EPS09[ibin];
         dA_up[ibin] = doct10err ? sqrt(pow(errA_EPS09_EPS09_up[ibin],2)+pow(errA_EPS09_scales_up[ibin],2)+pow(errA_EPS09_CT10_up[ibin],2)) : fabs(errA_EPS09_EPS09_up[ibin]);
         dA_down[ibin] = doct10err ? sqrt(pow(errA_EPS09_EPS09_down[ibin],2)+pow(errA_EPS09_scales_down[ibin],2)+pow(errA_EPS09_CT10_down[ibin],2)) : fabs(errA_EPS09_EPS09_down[ibin]);
      }
   }

   for (int ibin=0; ibin<nbins; ibin++)
      cout << ibin << " " << y_median[ibin] << " " << A[ibin] << " +" << dA_up[ibin] << " -" << dA_down[ibin] << endl;

   TGraphAsymmErrors *gr = new TGraphAsymmErrors(nbins,y_median,A,erry,erry,dA_down,dA_up);

   return gr;
}

TGraphErrors* theory_alice_plus()
{
   double x[2] = {-3.245,3.245};
   double ex[2] = {0.75,0.75};
   double y[2] = {102.2,18.3}; // values "stolen" from the ALICE plot
   for (int i=0; i<2; i++) y[i]=y[i]*LUMI/acc_alice_plus[i]/(2*ex[i]);
   double ey[2] = {0,0};

   TGraphErrors *gr = new TGraphErrors(2,x,y,ex,ey);
   return gr;
}

TGraphErrors* theory_alice_minus()
{
   double x[2] = {-3.245,3.245};
   double ex[2] = {0.75,0.75};
   double y[2] = {84.5,73.6}; // values "stolen" from the ALICE plot
   for (int i=0; i<2; i++) y[i]=y[i]*LUMI/acc_alice_minus[i]/(2*ex[i]);
   double ey[2] = {0,0};

   TGraphErrors *gr = new TGraphErrors(2,x,y,ex,ey);
   return gr;
}

TGraphErrors* data_alice_plus_stat()
{
   double x[2] = {-3.245,3.245};
   double ex[2] = {0.75,0.75};
   double y[2] = {86.2,12.8}; 
   double ey[2] = {5.4,2.0};
   for (int i=0; i<2; i++) 
   {
      y[i]=y[i]*LUMI/acc_alice_plus[i]/(2*ex[i]);
      ey[i]=ey[i]*(ex[i])*LUMI/acc_alice_plus[i]/(2*ex[i]);
   }

   TGraphErrors *gr = new TGraphErrors(2,x,y,ex,ey);
   return gr;
}

TGraphErrors* data_alice_plus_syst()
{
   double x[2] = {-3.245,3.245};
   double ex[2] = {0.75,0.75};
   double y[2] = {86.2,12.8}; 
   double ey[2] = {8.6,2.8};
   for (int i=0; i<2; i++) 
   {
      y[i]=y[i]*LUMI/acc_alice_plus[i]/(2*ex[i]);
      ey[i]=sqrt(pow(ey[i]*(ex[i])*LUMI/acc_alice_plus[i],2)+pow(y[i]*acc_alice_minus_err[i]/acc_alice_minus[i],2))/(2*ex[i]);
   }

   TGraphErrors *gr = new TGraphErrors(2,x,y,ex,ey);
   return gr;
}

TGraphErrors* data_alice_minus_stat()
{
   double x[2] = {-3.245,3.245};
   double ex[2] = {0.75,0.75};
   double y[2] = {82.9,75.9};
   double ey[2] = {5.6,5.4};
   for (int i=0; i<2; i++) 
   {
      y[i]=y[i]*LUMI/acc_alice_minus[i]/(2*ex[i]);
      ey[i]=ey[i]*(ex[i])*LUMI/acc_alice_minus[i]/(2*ex[i]);
   }

   TGraphErrors *gr = new TGraphErrors(2,x,y,ex,ey);
   return gr;
}

TGraphErrors* data_alice_minus_syst()
{
   double x[2] = {-3.245,3.245};
   double ex[2] = {0.75,0.75};
   double y[2] = {82.9,75.9};
   double ey[2] = {14.8,10.4};
   for (int i=0; i<2; i++) 
   {
      y[i]=y[i]*LUMI/acc_alice_minus[i]/(2*ex[i]);
      ey[i]=sqrt(pow(ey[i]*(ex[i])*LUMI/acc_alice_minus[i],2)+pow(y[i]*acc_alice_minus_err[i]/acc_alice_minus[i],2))/(2*ex[i]);
   }

   TGraphErrors *gr = new TGraphErrors(2,x,y,ex,ey);
   return gr;
}

void plot_graph_1file(const char* fname="graph.root", const char *gbasename="gyields", int channel_number=0, bool reverteta=false, double scale=1.)
{
   gStyle->SetOptTitle(0);
   gStyle->SetEndErrorSize(6);
   gStyle->SetHatchesLineWidth(2);
   gStyle->SetHatchesSpacing(1.7);
   TFile *f = new TFile(fname);

   float lTextSize = gTextSize;

   bool twopads = (channel_number != 0 && (string(gbasename)=="gyieldsm" || string(gbasename)=="gyieldsp"));

   TCanvas *c1 = new TCanvas("c1","c1",600,600);
   c1->cd();
	TPad *pad1 = new TPad("pad1","pad1",0,twopads ? 0.26 : 0,1,1);
	c1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3*0.96);
   if (twopads)
   {
      pad1->SetBottomMargin(0.04);
      pad2->SetTopMargin(0);
      pad2->SetFillColor(0);
      pad2->SetFillStyle(0);
      pad2->SetBottomMargin(gStyle->GetPadBottomMargin()/0.3);
      pad1->SetTopMargin(gStyle->GetPadTopMargin()/0.7);
      pad2->SetGridy();
      pad1->Draw();
      pad1->cd();

      // adapt text size
      lTextSize *= 1./0.7;
   }
   else
   {
      c1->cd();
   }

   string xlabel, ylabel;
   double xmin, xmax, ymin, ymax;

   // x-axis
   xlabel = string("#eta_{lab}");
   if (string(gbasename)=="gyieldsp" || string(gbasename)=="gyieldsm" || string(gbasename)=="gch")
   {
      xmin = !doAlice ? -2.5 : -4.1;
      xmax = !doAlice ? 2.5 : 4.1;
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
      if (string(gbasename)=="gyieldsp") ylabel = string("d#sigma (W^{+}#rightarrow#font[12]{l}^{+}#nu) / d#eta_{lab} [nb]");
      else ylabel = string("d#sigma (W^{-}#rightarrow#font[12]{l}^{-}#nu) / d#eta_{lab} [nb]");
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
   TGraphErrors *galice_stat = (string(gbasename)=="gyieldsp") ? data_alice_plus_stat() : data_alice_minus_stat();
   TGraphErrors *galice_syst = (string(gbasename)=="gyieldsp") ? data_alice_plus_syst() : data_alice_minus_syst();
   // TFile *fc = new TFile("combo_muacc_20140219.root");
   int nbins = gexp_stat1->GetN();
   TGraphAsymmErrors *gth_cteq = theory_errors(gbasename, twochan ? 1 : channel_number, nbins, "CT10");
   TGraphAsymmErrors *gth_eps = theory_errors(gbasename, twochan ? 1 : channel_number, nbins, "EPS09");
   // TGraphAsymmErrors *gth_cteq = theory_errors(gbasename, twochan ? 1 : channel_number, 12, "CT10");
   // TGraphAsymmErrors *gth_eps = theory_errors(gbasename, twochan ? 1 : channel_number, 12, "EPS09");
   TGraphErrors *galice_th = (string(gbasename)=="gyieldsp") ? theory_alice_plus() : theory_alice_minus();

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
      revertX(galice_stat);
      revertX(galice_syst);
      revertX(galice_th);
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
   scaleY(galice_stat,scale);
   scaleY(galice_syst,scale);
   scaleY(galice_th,scale);
   setErrorX(galice_stat,0);
   setErrorX(galice_syst,0);
   galice_stat->SetMarkerSize(1.2);
   galice_syst->SetMarkerSize(1.2);

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
   if (twopads) gth_cteq->GetXaxis()->SetTitle("");
   gth_cteq->GetXaxis()->SetNdivisions(505);
   gth_cteq->GetXaxis()->SetTitleSize(lTextSize);
   gth_cteq->GetXaxis()->SetLabelSize(lTextSize);
   gth_cteq->GetYaxis()->SetTitle(ylabel.c_str());
   gth_cteq->GetYaxis()->SetTitleSize(lTextSize);
   gth_cteq->GetYaxis()->SetLabelSize(lTextSize);
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

   if (doAlice)
   {
      galice_th->SetLineColor(gColorCT10);
      galice_th->SetMarkerSize(0);
      galice_th->SetLineWidth(2);
      galice_th->SetFillColor(gColorCT10);
      galice_th->SetFillStyle(1001);
      // galice_th->Draw("L3");
      galice_th->Draw("2");
      cout << galice_th->GetX()[0] << " " << galice_th->GetY()[0] << endl;
   }

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

   if (doAlice)
   {
      galice_syst->SetLineColor(kGreen);
      galice_syst->SetMarkerColor(kGreen);
      galice_syst->SetMarkerStyle(kFullStar);
      galice_syst->SetFillStyle(0);
      galice_syst->SetLineWidth(2);
      galice_syst->Draw("||");

      galice_stat->SetLineColor(kGreen);
      galice_stat->SetMarkerColor(kGreen);
      galice_stat->SetMarkerStyle(kFullStar);
      galice_stat->SetLineWidth(2);
      galice_stat->Draw("PZ");
   }

   double xl=0.4,yl=0.18,dx=0.3,dy=0.15;
   if (string(gbasename)=="gA1m" || string(gbasename)=="gA1p" || string(gbasename)=="gA3" || string(gbasename)=="gA4")
   {
      xl=0.21;
      yl=0.62;
   }
   if (string(gbasename)=="gyieldsm") xl=0.3;
   if (twopads) {yl=0.1; dy*=1./0.7;}
   if (string(gbasename)=="gA1m") xl=0.5;
   if (string(gbasename)=="gA3") yl=0.45;
   if ((string(gbasename)=="gch"||string(gbasename)=="gyieldsp")&&!reverteta) xl = 0.2;
   // double shiftx = dx/2.;

   TLegend *tleg2 = new TLegend(xl,yl,xl+dx,yl+dy);
   tleg2->SetTextFont(42);
   tleg2->SetTextSize(lTextSize*0.8);
   string legend1 = twochan ? string(channeltext1) : string("Data"); //string(channeltext);
   if (string(gbasename)=="gyieldsp" || string(gbasename)=="gA1p") legend1 = twochan ? string(channeltext1_p) : string("Data"); //string(channeltext_p);
   else if (string(gbasename)=="gyieldsm" || string(gbasename)=="gA1m") legend1 = twochan ? string(channeltext1_m) : string("Data"); //string(channeltext_m);
   tleg2->AddEntry(gexp_syst1,legend1.c_str(),"p");
   if (twochan) 
   {
      string legend2 = string(channeltext2);
      if (string(gbasename)=="gyieldsp" || string(gbasename)=="gA1p") legend2 = string(channeltext2_p);
      else if (string(gbasename)=="gyieldsm" || string(gbasename)=="gA1m") legend2 = string(channeltext2_m);
      tleg2->AddEntry(gexp_syst2,legend2.c_str(),"p");
   }
   else
   {
      // if (string(gbasename)=="gyieldsp"||string(gbasename)=="gyieldsm")
      // {
      //    tleg2->AddEntry(gth_cteq,"208#times(CT10)","lf");
      //    tleg2->AddEntry(gth_eps,"208#times(CT10+EPS09)","lf");
      // }
      // else
      // {
      tleg2->AddEntry(gth_cteq,"CT10","lf");
      tleg2->AddEntry(gth_eps,"CT10+EPS09","lf");
      // }
   }
   if (doAlice) tleg2->AddEntry(galice_syst,"ALICE data", "p");
   tleg2->SetFillColor(kWhite);
   tleg2->SetBorderSize(0);
   tleg2->Draw();

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
      tp->SetTextSize(lTextSize*0.8);
      tp->SetFillStyle(0);
      tp->Draw();
   }

   if (string(gbasename)=="gyieldsp" || string(gbasename)=="gyieldsm" || string(gbasename)=="gA1p" || string(gbasename)=="gA1m" || string(gbasename)=="gA3")
   {
      double txl=0.38, tyl=60, tdx=1.,tdy=20.;
      if (string(gbasename)=="gyieldsm") tyl=95;
      if (!reverteta) txl = -1.8;
      if (string(gbasename)=="gA1p") {txl = 0.38; tdx = 0.5; tyl = 1.75; tdy = 0.2;}
      if (string(gbasename)=="gA1m") {txl = 0.25; tdx = 0.5; tyl = 1.2; tdy = 0.1;}
      if (string(gbasename)=="gA3") {txl = 0.25; tdx = 0.5; tyl = 1.5; tdy = 0.1;}
      TPaveText *tp = new TPaveText(txl,tyl,txl+tdx,tyl+tdy);
      string legend;
      if (string(gbasename)=="gyieldsp" || string(gbasename)=="gA1p") legend = string(channeltext_p);
      else if (string(gbasename)=="gyieldsm" || string(gbasename)=="gA1m") legend = string(channeltext_m);
      else legend = string(channeltext);
      tp->AddText(legend.c_str());
      tp->SetFillColor(kWhite);
      tp->SetBorderSize(0);
      tp->SetTextFont(42);
      tp->SetTextSize(lTextSize);
      tp->SetTextAlign(13);
      tp->Draw();
   }


   if (twopads)
   {
      // adapt text size
      lTextSize *= 0.7/0.3;

      // fix pad 1
      gth_cteq->GetYaxis()->SetTitleOffset(gStyle->GetTitleYOffset()*0.7);
      gth_cteq->GetXaxis()->SetLabelSize(0);
      pad1->Update();

      c1->cd();
      pad2->SetFrameFillStyle(4000);
      pad2->Draw();
      pad2->cd();
      TH1F *haxis = new TH1F(*(gth_cteq->GetHistogram()));
      haxis->GetYaxis()->SetRangeUser(0.75,1.35);
      haxis->GetYaxis()->SetTitle("Ratio to CT10");
      haxis->GetXaxis()->SetTitle("#eta_{lab}");
      haxis->GetXaxis()->SetTitleSize(lTextSize);
      haxis->GetXaxis()->SetLabelSize(lTextSize);
      haxis->GetYaxis()->SetTitleSize(lTextSize);
      haxis->GetYaxis()->SetLabelSize(lTextSize);
      haxis->GetYaxis()->SetTitleOffset(gStyle->GetTitleYOffset()*0.3);
      haxis->GetYaxis()->SetNdivisions(503);
      haxis->Draw();

      TGraphErrors *gratio_stat = divideGraphs(gexp_stat1, convert(gth_cteq));
      TGraphErrors *gratio_syst = divideGraphs(gexp_syst1, convert(gth_cteq));
      gratio_stat->SetLineColor((twochan) ? gMyColor1 : gMyColor);
      gratio_stat->SetMarkerColor((twochan) ? gMyColor1 : gMyColor);
      gratio_stat->SetMarkerStyle(gMyMarker1);
      gratio_stat->SetLineWidth(2);
      gratio_syst->SetLineColor((twochan) ? gMyColor1 : gMyColor);
      gratio_syst->SetMarkerColor((twochan) ? gMyColor1 : gMyColor);
      gratio_syst->SetMarkerStyle(gMyMarker1);
      gratio_syst->SetFillStyle(0);
      gratio_syst->SetLineWidth(2);
      TGraphErrors *gratio_ct10err = divideGraphs(convert(gth_cteq), convert(gth_cteq));
      setErrorX(gratio_syst,0);
      setErrorX(gratio_stat,0);
      gratio_ct10err->SetLineColor(gColorCT10);
      gratio_ct10err->SetMarkerSize(0);
      gratio_ct10err->SetLineWidth(2);
      gratio_ct10err->SetFillColor(gColorCT10_fill);
      gratio_ct10err->SetFillStyle(1001);
      gratio_ct10err->Draw("2");

      TGraphErrors *gratio_ct10err2 = new TGraphErrors(*gratio_ct10err);
      setErrorY(gratio_ct10err2,0);
      gratio_ct10err2->Draw("Z");

      TGraphAsymmErrors *gth_eps_noct10err = theory_errors(gbasename, twochan ? 1 : channel_number, nbins, "EPS09",false);
      if (!twochan) scaleY(gth_eps_noct10err,scale);
      revertX(gth_eps_noct10err);
      TGraphErrors *gratio_eps = divideGraphs(convert(gth_eps_noct10err), convert(gth_cteq));
      for (int i=0; i<gratio_eps->GetN(); i++) cout << gratio_eps->GetX()[i] << " " << gratio_eps->GetY()[i] << "+-" << gratio_eps->GetEY()[i] << endl;
      gratio_eps->SetLineColor(gColorEPS09);
      gratio_eps->SetLineStyle(9);
      gratio_eps->SetLineWidth(2);
      gratio_eps->SetFillColor(gColorEPS09);
      gratio_eps->SetFillStyle(3345);
      gratio_eps->Draw("2");

      TGraphErrors *gratio_eps2 = new TGraphErrors(*gratio_eps);
      setErrorY(gratio_eps2,0);
      gratio_eps2->Draw("Z");

      gratio_syst->Draw("||");
      gratio_stat->Draw("PZ");

      pad1->cd();
      pad1->Draw();
   }

   int iPos = (string(gbasename)=="gch") ? 11 : 33;
   CMS_lumi( gPad, iPeriod, iPos );

   c1->cd();
   c1->Update();
   c1->RedrawAxis();
   pad1->RedrawAxis();
   pad2->RedrawAxis();

   gPad->SaveAs(Form("%s.pdf",gbasename));
   gPad->SaveAs(Form("%s.png",gbasename));
   gPad->SaveAs(Form("%s.root",gbasename));
   gPad->SaveAs(Form("%s.C",gbasename));
}

void plot_graph_1file_all(const char* filename, int chan, bool reverse_eta=true, double scale=1./34.62)
{
   plot_graph_1file(filename,"gA3",chan);
   plot_graph_1file(filename,"gch",chan,reverse_eta);
   plot_graph_1file(filename,"gA1p",chan);
   plot_graph_1file(filename,"gA1m",chan);
   plot_graph_1file(filename,"gA4",chan);
   plot_graph_1file(filename,"gyieldsp",chan,reverse_eta,scale);
   plot_graph_1file(filename,"gyieldsm",chan,reverse_eta,scale);
}

void plot_graph_2file(const char* fname1="graph.root", const char* fname2="graph.root", const char *gbasename="gyields", bool reverteta=false, double scale=1., bool threegraphs=false)
{
   gStyle->SetOptTitle(0);
   gStyle->SetEndErrorSize(6);
   gStyle->SetHatchesLineWidth(2);
   gStyle->SetHatchesSpacing(1.7);
   TFile *f1 = new TFile(fname1);
   TFile *f2 = new TFile(fname2);

   float lTextSize = gTextSize;
   bool twopads = false;

   TCanvas *c1 = new TCanvas("c1","c1",600,600);
   c1->cd();
	TPad *pad1 = new TPad("pad1","pad1",0,0.26,1,1);
	c1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3*0.96);
   if (twopads)
   {
      pad1->SetBottomMargin(0.04);
      pad2->SetTopMargin(0);
      pad2->SetFillColor(0);
      pad2->SetFillStyle(0);
      pad2->SetBottomMargin(gStyle->GetPadBottomMargin()/0.3);
      pad1->SetTopMargin(gStyle->GetPadTopMargin()/0.7);
      pad2->SetGridy();
      pad1->Draw();
      pad1->cd();

      // adapt text size
      lTextSize *= 1./0.7;
   }
   else
   {
      c1->cd();
   }

   string xlabel, ylabel;
   double xmin, xmax, ymin, ymax;

   // x-axis
   xlabel = string("#eta_{lab}");
   if (string(gbasename)=="gyieldsp" || string(gbasename)=="gyieldsm" || string(gbasename)=="gch")
   {
      xmin = !doAlice ? -2.5 : -4.1;
      xmax = !doAlice ? 2.5 : 4.1;
   }
   else
   {
      xmin = 0.;
      xmax = 2.5;
   }

   int channel_number=0;
   bool twochan = (channel_number==0 || channel_number==3);
   if (twochan) channel_number = 1;

   // y-axis range
   if (string(gbasename)=="gyieldsp" || string(gbasename)=="gyieldsm")
   {
      ymin = 0.;
      ymax = 140.;
      if (string(gbasename)=="gyieldsp") ylabel = string("d#sigma (W^{+}#rightarrow#font[12]{l}^{+}#nu) / d#eta_{lab} [nb]");
      else ylabel = string("d#sigma (W^{-}#rightarrow#font[12]{l}^{-}#nu) / d#eta_{lab} [nb]");
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
   TString gname_stat2 = Form("%s_exp_statonly_%i",gbasename,threegraphs ? 2 : channel_number);
   TString gname_syst2 = Form("%s_exp_%i",gbasename,threegraphs ? 2 : channel_number);

   TGraphErrors *gexp_stat1 = (TGraphErrors*) f1->Get(gname_stat1);
   TGraphErrors *gexp_syst1 = (TGraphErrors*) f1->Get(gname_syst1);
   TGraphErrors *gexp_stat2 = (TGraphErrors*) f2->Get(gname_stat2);
   TGraphErrors *gexp_syst2 = (TGraphErrors*) f2->Get(gname_syst2);
   TGraphErrors *gexp_stat3 = (TGraphErrors*) f1->Get(gname_stat2);
   TGraphErrors *gexp_syst3 = (TGraphErrors*) f1->Get(gname_syst2);
   TGraphErrors *galice_stat = (string(gbasename)=="gyieldsp") ? data_alice_plus_stat() : data_alice_minus_stat();
   TGraphErrors *galice_syst = (string(gbasename)=="gyieldsp") ? data_alice_plus_syst() : data_alice_minus_syst();
   // TFile *fc = new TFile("combo_muacc_20140219.root");
   int nbins = gexp_stat1->GetN();
   TGraphAsymmErrors *gth_cteq = theory_errors(gbasename, twochan ? 1 : channel_number, nbins, "CT10");
   TGraphAsymmErrors *gth_eps = theory_errors(gbasename, twochan ? 1 : channel_number, nbins, "EPS09");
   // TGraphAsymmErrors *gth_cteq = theory_errors(gbasename, twochan ? 1 : channel_number, 12, "CT10");
   // TGraphAsymmErrors *gth_eps = theory_errors(gbasename, twochan ? 1 : channel_number, 12, "EPS09");
   TGraphErrors *galice_th = (string(gbasename)=="gyieldsp") ? theory_alice_plus() : theory_alice_minus();

   if (reverteta)
   {
      revertX(gexp_stat1);
      revertX(gexp_syst1);
      if (twochan)
      {
         revertX(gexp_stat2);
         revertX(gexp_syst2);
         revertX(gexp_stat3);
         revertX(gexp_syst3);
      }
      revertX(gth_cteq);
      revertX(gth_eps);
      revertX(galice_stat);
      revertX(galice_syst);
      revertX(galice_th);
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
      scaleY(gexp_stat3,scale);
      scaleY(gexp_syst3,scale);
      setErrorX(gexp_stat2,0);
      setErrorX(gexp_syst2,0);
      setErrorX(gexp_stat3,0);
      setErrorX(gexp_syst3,0);
      gexp_stat2->SetMarkerSize(1.2);
      gexp_syst2->SetMarkerSize(1.2);
      gexp_stat3->SetMarkerSize(1.2);
      gexp_syst3->SetMarkerSize(1.2);
      shiftX(gexp_stat1,-0.05);
      shiftX(gexp_syst1,-0.05);
      shiftX(gexp_stat2,0.05);
      shiftX(gexp_syst2,0.05);
      shiftX(gexp_stat3,0.05);
      shiftX(gexp_syst3,0.05);
   }
   else
   {
      scaleY(gth_cteq,scale);
      scaleY(gth_eps,scale);
   }
   scaleY(galice_stat,scale);
   scaleY(galice_syst,scale);
   scaleY(galice_th,scale);
   setErrorX(galice_stat,0);
   setErrorX(galice_syst,0);
   galice_stat->SetMarkerSize(1.2);
   galice_syst->SetMarkerSize(1.2);

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
   if (twopads) gth_cteq->GetXaxis()->SetTitle("");
   gth_cteq->GetXaxis()->SetTitleSize(lTextSize);
   gth_cteq->GetXaxis()->SetLabelSize(lTextSize);
   gth_cteq->GetYaxis()->SetTitle(ylabel.c_str());
   gth_cteq->GetYaxis()->SetTitleSize(lTextSize);
   gth_cteq->GetYaxis()->SetLabelSize(lTextSize);
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

   if (doAlice)
   {
      galice_th->SetLineColor(gColorCT10);
      galice_th->SetMarkerSize(0);
      galice_th->SetLineWidth(2);
      galice_th->SetFillColor(gColorCT10);
      galice_th->SetFillStyle(1001);
      // galice_th->Draw("L3");
      galice_th->Draw("2");
      cout << galice_th->GetX()[0] << " " << galice_th->GetY()[0] << endl;
   }

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

      gexp_syst3->SetLineColor(kGreen);
      gexp_syst3->SetMarkerColor(kGreen);
      gexp_syst3->SetMarkerStyle(34);
      gexp_syst3->SetFillStyle(0);
      gexp_syst3->SetLineWidth(2);
      if (threegraphs) gexp_syst3->Draw("||");

      gexp_stat3->SetLineColor(kGreen);
      gexp_stat3->SetMarkerColor(kGreen);
      gexp_stat3->SetMarkerStyle(34);
      gexp_stat3->SetLineWidth(2);
      if (threegraphs) gexp_stat3->Draw("PZ");
   }

   if (doAlice)
   {
      galice_syst->SetLineColor(kGreen);
      galice_syst->SetMarkerColor(kGreen);
      galice_syst->SetMarkerStyle(kFullStar);
      galice_syst->SetFillStyle(0);
      galice_syst->SetLineWidth(2);
      galice_syst->Draw("||");

      galice_stat->SetLineColor(kGreen);
      galice_stat->SetMarkerColor(kGreen);
      galice_stat->SetMarkerStyle(kFullStar);
      galice_stat->SetLineWidth(2);
      galice_stat->Draw("PZ");
   }

   double xl=0.5,yl=0.18,dx=0.3,dy=0.15;
   if (string(gbasename)=="gA1m" || string(gbasename)=="gA1p" || string(gbasename)=="gA3" || string(gbasename)=="gA4")
   {
      xl=0.21;
      yl=0.62;
   }
   if (string(gbasename)=="gyieldsm") xl=0.3;
   if (twopads) {yl=0.1; dy*=1./0.7;}
   if (string(gbasename)=="gA1m") xl=0.6;
   if (string(gbasename)=="gA3") yl=0.45;
   if ((string(gbasename)=="gch"||string(gbasename)=="gyieldsp")&&!reverteta) xl = 0.2;
   // double shiftx = dx/2.;

   TLegend *tleg2 = new TLegend(xl,yl,xl+dx,yl+dy);
   tleg2->SetTextFont(42);
   tleg2->SetTextSize(lTextSize*0.8);
   string legend1 = twochan ? string(channeltext1) : string("Data"); //string(channeltext);
   if (string(gbasename)=="gyieldsp" || string(gbasename)=="gA1p") legend1 = twochan ? string(channeltext1_p) : string("Data"); //string(channeltext_p);
   else if (string(gbasename)=="gyieldsm" || string(gbasename)=="gA1m") legend1 = twochan ? string(channeltext1_m) : string("Data"); //string(channeltext_m);
   tleg2->AddEntry(gexp_syst1,legend1.c_str(),"p");
   if (twochan) 
   {
      string legend2 = string(channeltext2);
      if (string(gbasename)=="gyieldsp" || string(gbasename)=="gA1p") legend2 = string(channeltext2_p);
      else if (string(gbasename)=="gyieldsm" || string(gbasename)=="gA1m") legend2 = string(channeltext2_m);
      string legend3 = string(channeltext3);
      if (string(gbasename)=="gyieldsp" || string(gbasename)=="gA1p") legend3 = string(channeltext3_p);
      else if (string(gbasename)=="gyieldsm" || string(gbasename)=="gA1m") legend3 = string(channeltext3_m);
      tleg2->AddEntry(gexp_syst2,legend2.c_str(),"p");
      if (threegraphs) tleg2->AddEntry(gexp_syst3,legend3.c_str(),"p");
   }
   else
   {
      tleg2->AddEntry(gth_cteq,"CT10","lf");
      tleg2->AddEntry(gth_eps,"CT10+EPS09","lf");
   }
   if (doAlice) tleg2->AddEntry(galice_syst,"ALICE data", "p");
   tleg2->SetFillColor(kWhite);
   tleg2->SetBorderSize(0);
   tleg2->Draw();

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
      tp->SetTextSize(lTextSize*0.8);
      tp->SetFillStyle(0);
      tp->Draw();
   }

   if (string(gbasename)=="gyieldsp" || string(gbasename)=="gyieldsm" || string(gbasename)=="gA1p" || string(gbasename)=="gA1m" || string(gbasename)=="gA3")
   {
      double txl=0.38, tyl=60, tdx=1.,tdy=20.;
      if (string(gbasename)=="gyieldsm") tyl=95;
      if (!reverteta) txl = -1.8;
      if (string(gbasename)=="gA1p") {txl = 0.38; tdx = 0.5; tyl = 1.75; tdy = 0.2;}
      if (string(gbasename)=="gA1m") {txl = 0.25; tdx = 0.5; tyl = 1.2; tdy = 0.1;}
      if (string(gbasename)=="gA3") {txl = 0.25; tdx = 0.5; tyl = 1.5; tdy = 0.1;}
      TPaveText *tp = new TPaveText(txl,tyl,txl+tdx,tyl+tdy);
      string legend;
      if (string(gbasename)=="gyieldsp" || string(gbasename)=="gA1p") legend = string(channeltext_p);
      else if (string(gbasename)=="gyieldsm" || string(gbasename)=="gA1m") legend = string(channeltext_m);
      else legend = string(channeltext);
      tp->AddText(legend.c_str());
      tp->SetFillColor(kWhite);
      tp->SetBorderSize(0);
      tp->SetTextFont(42);
      tp->SetTextSize(lTextSize);
      tp->SetTextAlign(13);
      tp->Draw();
   }


   if (twopads)
   {
      // adapt text size
      lTextSize *= 0.7/0.3;

      // fix pad 1
      gth_cteq->GetYaxis()->SetTitleOffset(gStyle->GetTitleYOffset()*0.7);
      gth_cteq->GetXaxis()->SetLabelSize(0);
      pad1->Update();

      c1->cd();
      pad2->SetFrameFillStyle(4000);
      pad2->Draw();
      pad2->cd();
      TH1F *haxis = new TH1F(*(gth_cteq->GetHistogram()));
      haxis->GetYaxis()->SetRangeUser(0.75,1.35);
      haxis->GetYaxis()->SetTitle("Ratio to CT10");
      haxis->GetXaxis()->SetTitle("#eta_{lab}");
      haxis->GetXaxis()->SetTitleSize(lTextSize);
      haxis->GetXaxis()->SetLabelSize(lTextSize);
      haxis->GetYaxis()->SetTitleSize(lTextSize);
      haxis->GetYaxis()->SetLabelSize(lTextSize);
      haxis->GetYaxis()->SetTitleOffset(gStyle->GetTitleYOffset()*0.3);
      haxis->GetYaxis()->SetNdivisions(503);
      haxis->Draw();

      TGraphErrors *gratio_stat = divideGraphs(gexp_stat1, convert(gth_cteq));
      TGraphErrors *gratio_syst = divideGraphs(gexp_syst1, convert(gth_cteq));
      gratio_stat->SetLineColor((twochan) ? gMyColor1 : gMyColor);
      gratio_stat->SetMarkerColor((twochan) ? gMyColor1 : gMyColor);
      gratio_stat->SetMarkerStyle(gMyMarker1);
      gratio_stat->SetLineWidth(2);
      gratio_syst->SetLineColor((twochan) ? gMyColor1 : gMyColor);
      gratio_syst->SetMarkerColor((twochan) ? gMyColor1 : gMyColor);
      gratio_syst->SetMarkerStyle(gMyMarker1);
      gratio_syst->SetFillStyle(0);
      gratio_syst->SetLineWidth(2);
      TGraphErrors *gratio_ct10err = divideGraphs(convert(gth_cteq), convert(gth_cteq));
      setErrorX(gratio_syst,0);
      setErrorX(gratio_stat,0);
      gratio_ct10err->SetLineColor(gColorCT10);
      gratio_ct10err->SetMarkerSize(0);
      gratio_ct10err->SetLineWidth(2);
      gratio_ct10err->SetFillColor(gColorCT10_fill);
      gratio_ct10err->SetFillStyle(1001);
      gratio_ct10err->Draw("2");

      TGraphErrors *gratio_ct10err2 = new TGraphErrors(*gratio_ct10err);
      setErrorY(gratio_ct10err2,0);
      gratio_ct10err2->Draw("Z");

      TGraphAsymmErrors *gth_eps_noct10err = theory_errors(gbasename, twochan ? 1 : channel_number, nbins, "EPS09",false);
      if (!twochan) scaleY(gth_eps_noct10err,scale);
      revertX(gth_eps_noct10err);
      TGraphErrors *gratio_eps = divideGraphs(convert(gth_eps_noct10err), convert(gth_cteq));
      for (int i=0; i<gratio_eps->GetN(); i++) cout << gratio_eps->GetX()[i] << " " << gratio_eps->GetY()[i] << "+-" << gratio_eps->GetEY()[i] << endl;
      gratio_eps->SetLineColor(gColorEPS09);
      gratio_eps->SetLineStyle(9);
      gratio_eps->SetLineWidth(2);
      gratio_eps->SetFillColor(gColorEPS09);
      gratio_eps->SetFillStyle(3345);
      gratio_eps->Draw("2");

      TGraphErrors *gratio_eps2 = new TGraphErrors(*gratio_eps);
      setErrorY(gratio_eps2,0);
      gratio_eps2->Draw("Z");

      gratio_syst->Draw("||");
      gratio_stat->Draw("PZ");

      pad1->cd();
      pad1->Draw();
   }

   int iPos = (string(gbasename)=="gch") ? 11 : 33;
   CMS_lumi( c1, iPeriod, iPos );

   c1->cd();
   c1->Update();

   gPad->SaveAs(Form("%s.pdf",gbasename));
   gPad->SaveAs(Form("%s.png",gbasename));
   gPad->SaveAs(Form("%s.root",gbasename));
   gPad->SaveAs(Form("%s.C",gbasename));
}

