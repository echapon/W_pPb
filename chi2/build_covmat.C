// double lumierr=0.035;
double lumierr=0.0362;
double puerr=0.;

TH2F *build_covmat(TGraphErrors *gstat, TGraphErrors *gsyst, bool docov=true)
{
   int n = gstat->GetN();

   TH2F *h2 = new TH2F("h2","h2",n,0,n,n,0,n);

   for (int i=0; i<n; i++)
   {
      double estat = gstat->GetEY()[i];
      double esyst = sqrt(pow(gsyst->GetEY()[i],2)-docov*pow(puerr*gsyst->GetY()[i],2));
      for (int j=0; j<n; j++)
      {
         double err=0;
         if (i==j) err = estat*estat+esyst*esyst;
         err += docov*(lumierr+puerr)*gstat->GetY()[i]*(lumierr+puerr)*gstat->GetY()[j];
         h2->SetBinContent(i+1,j+1,err);
      }
   }

   return h2;
}
