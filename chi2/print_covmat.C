void print_covmat(const char *file, const char *covmat)
{
   TFile *f = new TFile(file);
   TH2F *hcov = (TH2F*) f->Get(Form("matrices/%s",covmat));

   cout.setf(std::ios::fixed);
   cout << setprecision(1);
   cout << "\\begin{tabular}{c|cccccccccc}" << endl;
   cout << "|\\eta| bin & [-2.4,-2] & [-2,-1.5] & [-1.5,-1] & [-1,-0.5] & [-0.5,0] & [0,0.5] & [0.5,1] & [1,1.5] & [1.5,2] & [2,2.5]\\\\" << endl;
   for (int i=0; i<10; i++)
   {
      cout << i;
      for (int j=0; j<10; j++)
      {
         if (j<=i) cout << " & " << hcov->GetBinContent(i+1,j+1)/34.62/34.62;
         else cout << " & ";
      }
      cout << "\\\\" << endl;
   }
   cout << "\\end{tabular}" << endl;
   f->Close();
}
