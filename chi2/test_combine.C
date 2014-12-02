double mysqrt(double a)
{
   return a>0 ? sqrt(a) : -sqrt(-a);
   // return a;
}

void test_combine()
{
   gROOT->ProcessLine(".L combine.C");

   int nbins = 4;
   double x[4] = {1,2,3,4};
   // double y[4] = {10.5,9.5,13.5,14.};
   double y[4] = {2.36e3/34.62,2.22e3/34.62,1.96e3/34.62,2.20e3/34.62};
   double ex[4] = {0.1,0.1,0.1,0.1};
   double ey[4] = {sqrt(5.04e4)/34.62,sqrt(2.31e4)/34.62,sqrt(5.45e4)/34.62,sqrt(5.67e4)/34.62};
   // double ey[4] = {sqrt(0.0145),sqrt(0.016),sqrt(0.017),sqrt(0.026)};
   // double ey[4] = {1,3,3,3};
   TGraphErrors *gr12 = new TGraphErrors(nbins,x,y,ex,ey);

   TH2F *hcov12 = new TH2F("hcov12","hcov12",nbins,0,nbins,nbins,0,nbins);
   for (int i=0; i<nbins; i++)
      for (int j=0; j<nbins; j++)
         hcov12->SetBinContent(i+1,j+1,i==j ? ey[i]*ey[i] : 0);
   // hcov12->SetBinContent(1,3,0.45);
   // hcov12->SetBinContent(3,1,0.45);
   hcov12->SetBinContent(3,4,3.71e4/34.62/34.62);
   hcov12->SetBinContent(4,3,3.71e4/34.62/34.62);
   // hcov12->SetBinContent(3,4,8.96);
   // hcov12->SetBinContent(4,3,8.96);
   for (int i=0; i<nbins; i++)
   {
      for (int j=0; j<nbins; j++)
         cout << hcov12->GetBinContent(i+1,j+1) << " ";
      cout << endl;
   }

   TGraphErrors *grcombo = new TGraphErrors(nbins/2);
   TH2F *hcovcombo = new TH2F("hcovcombo","hcovcombo",nbins/2,0,nbins/2,nbins/2,0,nbins/2);
   combine_blue(gr12,hcov12,grcombo,hcovcombo);
   for (int i=0; i<nbins/2; i++)
   {
      for (int j=0; j<nbins/2; j++)
         cout << hcovcombo->GetBinContent(i+1,j+1) << " ";
      cout << endl;
   }
}
