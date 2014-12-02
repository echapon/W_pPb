#include "TH1F.h"
#include "TF1.h"
#include "RooStats/HypoTestCalculator.h"

TF1* tfchi2()
{
   TF1 *tfchi2 = new TF1("tfchi2","[0]*(pow(x,[1]/2.-1.)*exp(-x/2.))/(pow(2.,[1]/2.)*TMath::Gamma([1]/2.))",0,50);
   tfchi2->SetParName(0,"Norm");
   tfchi2->SetParName(1,"Ndf");
   tfchi2->SetParameter(0,500);
   tfchi2->SetParameter(1,10);
   return tfchi2;
}

void plotCL(TH1F *hist)
{
   hist->Scale(1./hist->Integral("width"));

   int med=0, m68=0, p68=0, m95=0, p95=0;

   double nbins = hist->GetNbinsX();
   double xmin = hist->GetBinLowEdge(1);
   double xmax = hist->GetBinLowEdge(nbins)+hist->GetBinWidth(nbins);
   TH1F *h68 = new TH1F("h68","h68",nbins,xmin,xmax);
   TH1F *h95 = new TH1F("h95","h95",nbins,xmin,xmax);

   double intg=0;

   for (int ibin=1; ibin<nbins+1; ibin++)
   {
      intg += hist->GetBinContent(ibin);
      
      if (!med && intg>=0.5) med=ibin;
      if (!m68 && intg>=0.16) m68=ibin;
      if (!p68 && intg>=0.84) p68=ibin;
      if (!m95 && intg>=0.025) m95=ibin;
      if (!p95 && intg>=0.975) p95=ibin;

      if (m68 && !p68) h68->SetBinContent(ibin,hist->GetBinContent(ibin));
      if (m95 && !p95) h95->SetBinContent(ibin,hist->GetBinContent(ibin));
   }

   double dmed = hist->GetBinCenter(med);
   double dm68 = hist->GetBinCenter(m68);
   double dp68 = hist->GetBinCenter(p68);
   double dm95 = hist->GetBinCenter(m95);
   double dp95 = hist->GetBinCenter(p95);

   cout << "median: " << dmed << endl;
   cout << "68%: [" << dm68 << "," << dp68 << "]" << endl;
   cout << "95%: [" << dm95 << "," << dp95 << "]" << endl;

   TLine *tl = new TLine();
   hist->Draw("hist");
   h95->SetFillColor(kYellow);
   h95->Draw("hist same");
   h68->SetFillColor(kGreen);
   h68->Draw("hist same");
   gPad->Update();
   double ymin = gPad->GetUymin();
   double ymax = gPad->GetUymax();
   tl->SetLineWidth(2);
   tl->DrawLine(dmed,ymin,dmed,ymax);
}

double plotCL(TH1F *histnuc, TH1F *histnonuc, double obs)
{
   histnuc->Scale(1./histnuc->Integral("width"));
   histnuc->SetLineColor(kRed);
   histnuc->SetLineWidth(2);
   histnonuc->Scale(1./histnonuc->Integral("width"));
   histnonuc->SetLineColor(kBlue);
   histnonuc->SetLineWidth(2);

   double nbins = histnonuc->GetNbinsX();
   double xmin = histnonuc->GetBinLowEdge(1);
   double xmax = histnonuc->GetBinLowEdge(nbins)+histnonuc->GetBinWidth(nbins);
   TH1F *h68 = new TH1F("h68","h68",nbins,xmin,xmax);
   TH1F *h95 = new TH1F("h95","h95",nbins,xmin,xmax);

   double amednonuc[5];
   double probnonuc[5] = {0.025,0.16,0.5,0.84,0.975};
   histnonuc->GetQuantiles(5,amednonuc,probnonuc);

   double dmed = amednonuc[2];
   double dm68 = amednonuc[1];
   double dp68 = amednonuc[3];
   double dm95 = amednonuc[0];
   double dp95 = amednonuc[4];

   for (int ibin=1; ibin<nbins+1; ibin++)
   {
      double bc = histnonuc->GetBinCenter(ibin);
      if (bc>=dm68&&bc<dp68) h68->SetBinContent(ibin,histnonuc->GetBinContent(ibin));
      if (bc>=dm95&&bc<dp95) h95->SetBinContent(ibin,histnonuc->GetBinContent(ibin));
   }

   double amednuc[1];
   double prob[1] = {0.5};
   histnuc->GetQuantiles(1,amednuc,prob);
   double mednuc = amednuc[0];

   cout << "No nuclear effects:" << endl;
   cout << "median: " << dmed << endl;
   cout << "68%: [" << dm68 << "," << dp68 << "]" << endl;
   cout << "95%: [" << dm95 << "," << dp95 << "]" << endl;
   cout << "Nuclear effects:" << endl;
   cout << "median: " << mednuc << endl;

   TLine *tl_mednuc = new TLine();
   TLine *tl_mednonuc = new TLine();
   TLine *tl_obs = new TLine();
   histnonuc->Draw("hist");
   histnonuc->GetXaxis()->SetTitle("#Delta#chi^{2}");
   histnonuc->GetYaxis()->SetTitle("Prob.");
   h95->SetFillColor(kYellow);
   h95->Draw("hist same");
   h68->SetFillColor(kGreen);
   h68->Draw("hist same");
   histnonuc->Draw("hist same");
   histnuc->Draw("hist same");
   gPad->Update();
   double ymin = gPad->GetUymin();
   double ymax = gPad->GetUymax();
   tl_mednonuc->SetLineWidth(2);
   tl_mednonuc->SetLineColor(kBlack);
   tl_mednonuc->SetLineStyle(2);
   tl_mednonuc->DrawLine(dmed,ymin,dmed,ymax);
   tl_mednuc->SetLineWidth(2);
   tl_mednuc->SetLineColor(kRed);
   tl_mednuc->SetLineStyle(3);
   tl_mednuc->DrawLine(mednuc,ymin,mednuc,ymax);
   tl_obs->SetLineColor(kBlack);
   tl_obs->SetLineWidth(2);
   tl_obs->DrawLine(obs,ymin,obs,ymax);

   TLegend *tleg = new TLegend(0.2,0.7,0.5,0.9);
   tleg->AddEntry(histnonuc,"CT10 #Delta#chi^{2}","l");
   tleg->AddEntry(h68,"CT10 #Delta#chi^{2} (68\% C.L.)","f");
   tleg->AddEntry(h95,"CT10 #Delta#chi^{2} (95\% C.L.)","f");
   tleg->AddEntry(histnuc,"EPS09 #Delta#chi^{2}","l");
   tleg->AddEntry(tl_obs,"Observed #Delta#chi^{2}","l");
   tleg->SetFillColor(0);
   tleg->SetBorderSize(0);
   tleg->Draw();

   double clnuc_obs = histnuc->Integral(histnuc->FindBin(obs),histnuc->GetNbinsX())/histnuc->Integral();
   double clnuc_exp = histnuc->Integral(histnuc->FindBin(dmed),histnuc->GetNbinsX())/histnuc->Integral();
   double clnonuc_obs = 1.-histnonuc->Integral(histnonuc->FindBin(obs),histnonuc->GetNbinsX())/histnonuc->Integral();
   double clnonuc_exp = 1.-histnonuc->Integral(histnonuc->FindBin(mednuc),histnonuc->GetNbinsX())/histnonuc->Integral();
   double cls_obs = clnuc_obs/clnonuc_obs;

   cout << "CL nuc:" << endl;
   cout << "Obs.: " << clnuc_obs << " or " << RooStats::PValueToSignificance(clnuc_obs) << " sigma";
   cout << " (exp. from CT10: " << clnuc_exp << " or " << RooStats::PValueToSignificance(clnuc_exp) << " sigma)" << endl;
   cout << "CL nonuc:" << endl;
   cout << "Obs.: " << clnonuc_obs << " or " << RooStats::PValueToSignificance(clnonuc_obs) << " sigma";
   cout << " (exp. from EPS09: " << clnonuc_exp << " or " << RooStats::PValueToSignificance(clnonuc_exp) << " sigma)" << endl;
   cout << "cls: " << cls_obs << endl;

   return clnonuc_obs;
}
