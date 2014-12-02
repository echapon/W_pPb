printgraph(const char* file, const char* graphname, const char* graphname_syst, bool reverseta=false, double scale=1./34.62)
{
   TFile * f = new TFile(file);
   TGraphErrors *gr = (TGraphErrors*) f->Get(graphname);
   TGraphErrors *grsyst = (TGraphErrors*) f->Get(graphname_syst);

   int nbins = gr->GetN();
   double *x = gr->GetX();
   if (reverseta)
      for (int i=0; i<nbins; i++)
         x[i] = gr->GetX()[nbins-i-1];
   double *y = gr->GetY();
   double *ey = gr->GetEY();
   double *ey_syst = grsyst->GetEY();

   for (int i=0; i<nbins; i++)
      cout << x[i] << ": " << y[i]*scale << "+-" << ey[i]*scale << "(stat.) +-" << ey_syst[i]*scale << "(syst.)" << endl;
}
