#include <TGraphErrors.h>
#include <math.h>

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

void ratio_graphs(TGraphErrors *g1_stat, TGraphErrors *g1_syst, TGraphErrors *g2_stat, TGraphErrors *g2_syst, bool systIsTotErr=true)
{
   int nbins = g1_stat->GetN();

   double *x = new double[nbins]; x = g1_stat->GetX();
   double *ex = new double[nbins]; ex = g1_stat->GetEX();

   double *y = new double[nbins];
   double *ey_stat = new double[nbins];
   double *ey_syst = new double[nbins];

   for (int i=0; i<nbins; i++)
   {
      double y1 = g1_stat->GetY()[i];
      double y2 = g2_stat->GetY()[i];
      double ey1_stat = g1_stat->GetEY()[i];
      double ey2_stat = g2_stat->GetEY()[i];
      double ey1_syst = g1_syst->GetEY()[i];
      double ey2_syst = g2_syst->GetEY()[i];

      y[i] = (y1-y2)/(y1+y2);
      ey_stat[i] = sqrt(pow((2*y2)/pow(y1+y2,2),2)*pow(ey1_stat,2) + pow((-2*y1)/pow(y1+y2,2),2)*pow(ey2_stat,2));
      ey_syst[i] = sqrt(pow((2*y2)/pow(y1+y2,2),2)*pow(ey1_syst,2) + pow((-2*y1)/pow(y1+y2,2),2)*pow(ey2_syst,2));
      if (systIsTotErr) ey_syst[i] = sqrt(pow(ey_stat[i],2) + pow(ey_syst[i],2));
   }

   TGraphErrors *gr_stat = new TGraphErrors(nbins,x,y,ex,ey_stat);
   revertX(gr_stat);
   TGraphErrors *gr_syst = new TGraphErrors(nbins,x,y,ex,ey_syst);
   revertX(gr_syst);

   // now draw these graphs

   TCanvas *c1 = new TCanvas();
   c1->SetGridy();
   gr_syst->SetLineColor(kBlack);
   gr_syst->SetMarkerColor(kBlack);
   gr_syst->SetMarkerStyle(8);
   gr_syst->SetFillStyle(0);
   gr_syst->SetLineWidth(2);
   gr_syst->Draw("AP5");
   double xmin=-2.5,xmax=2.5;
   double ymin=-.2,ymax=.2;
   gr_syst->GetXaxis()->SetRangeUser(xmin,xmax);
   gr_syst->GetYaxis()->SetRangeUser(ymin,ymax);
   gr_syst->GetXaxis()->SetTitle("#eta");
   gr_syst->GetYaxis()->SetTitle("(N(#mu^{+})-N(e^{+})) / (N(#mu^{+})+N(e^{+}))");
   gr_syst->Draw("AP5");

   gr_stat->SetLineColor(kBlack);
   gr_stat->SetMarkerColor(kBlack);
   gr_stat->SetMarkerStyle(8);
   gr_stat->SetLineWidth(2);
   gr_stat->Draw("PZ");

   // now the pull histo
   TCanvas *c2 = new TCanvas();
   c2->SetGridy();
   c2->cd();

   double *xedges = new double[nbins+3];
   xedges[0] = -2.5;
   xedges[1] = x[0]-ex[0];
   for (int i=0; i<nbins; i++)
      xedges[i+2] = x[i]+ex[i];
   xedges[nbins+2]=2.5;
   TH1F *h1 = new TH1F("pulls","pulls",nbins+2,xedges);
   for (int i=0; i<nbins; i++)
      h1->SetBinContent(nbins-i+1,y[i]/ey_syst[i]);
   double xmin=-2.5,xmax=2.5;
   double ymin=-3,ymax=3;
   h1->GetXaxis()->SetRangeUser(xmin,xmax);
   h1->GetYaxis()->SetRangeUser(ymin,ymax);
   h1->GetXaxis()->SetTitle("#eta");
   h1->GetYaxis()->SetTitle("pull");
   h1->SetFillColor(kBlack);
   h1->Draw("B");
}
