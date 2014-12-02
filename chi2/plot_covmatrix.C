#include <fstream>
#include <iostream>
#include "TH2.h"

using namespace std;

void plot_covmatrix(int nbins, const char* filename)
{
   TH2F *h2 = new TH2F("h2","h2",nbins,0,nbins,nbins,0,nbins);
   ifstream file(filename);
   double buf;
   for (int ibin=0; ibin<nbins; ibin++)
      for (int jbin=0; jbin<nbins; jbin++)
      {
         file >> buf;
         h2->Fill(ibin,jbin,buf);
      }

   h2->Draw("COLZ");
}

void plot_cormatrix(int nbins, const char* filename)
{
   TH2F *h2 = new TH2F("h2","h2",nbins,0,nbins,nbins,0,nbins);
   ifstream file(filename);
   double buf;
   for (int ibin=0; ibin<nbins; ibin++)
      for (int jbin=0; jbin<nbins; jbin++)
      {
         file >> buf;
         h2->Fill(ibin,jbin,buf);
      }

   TH2F *h2c = new TH2F("h2c","h2c",nbins,0,nbins,nbins,0,nbins);
   for (int ibin=1; ibin<nbins+1; ibin++)
      for (int jbin=1; jbin<nbins+1; jbin++)
      {
         h2c->SetBinContent(ibin,jbin,h2->GetBinContent(ibin,jbin)/sqrt(h2->GetBinContent(ibin,ibin))/sqrt(h2->GetBinContent(jbin,jbin)));
      }

   h2c->Draw("COLZ");
}
