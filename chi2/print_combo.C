#define NBINS 10
#include "../plot_graph.C"

void print_combo(const char *file_split, const char* file_combo, const char* gbasename)
{
   TFile *fs = new TFile(file_split);
   TFile *fc = new TFile(file_combo);

   TGraphErrors *g1_syst = (TGraphErrors*) fs->Get(Form("%s_exp_1",gbasename));
   TGraphErrors *g2_syst = (TGraphErrors*) fs->Get(Form("%s_exp_2",gbasename));
   TGraphErrors *gc_syst = (TGraphErrors*) fc->Get(Form("%s_exp_1",gbasename));
   TGraphErrors *g1_stat = (TGraphErrors*) fs->Get(Form("%s_exp_statonly_1",gbasename));
   TGraphErrors *g2_stat = (TGraphErrors*) fs->Get(Form("%s_exp_statonly_2",gbasename));
   TGraphErrors *gc_stat = (TGraphErrors*) fc->Get(Form("%s_exp_statonly_1",gbasename));

   TGraphAsymmErrors *gnuc = theory_errors(gbasename,1,NBINS,"EPS09");
   TGraphAsymmErrors *gnonuc = theory_errors(gbasename,1,NBINS,"CT10");

   cout.setf(std::ios::fixed);
   cout << setprecision(1);
   cout << "\\begin{center}" << endl;
   cout << "\\resizebox{\\textwidth}{!}{" << endl;
   cout << "\\begin{tabular}{|l|ccccc|}" << endl;
   cout << "\\hline" << endl;
   cout << "$\\frac{d\\sigma}{d\\eta} (\\eta)$ & $[-2.4,-2.0] $ & $[-2.0,-1.5]$& $[-1.5,-1.0]$ & $[-1.0,-0.5]$ & $[-0.5,0]$\\\\" << endl;
   cout << "\\hline" << endl;
   cout << "$\\Pgm^-$";
   for (int i=NBINS-1; i>NBINS/2-1; i--)
      cout << " & $" << g1_syst->GetY()[i]/34.62 << "\\pm" << g1_stat->GetEY()[i]/34.62 << "\\pm" << g1_syst->GetEY()[i]/34.62 << "$";
   cout << "\\\\" << endl;
   cout << "$\\Pe^-$";
   for (int i=NBINS-1; i>NBINS/2-1; i--)
      cout << " & $" << g2_syst->GetY()[i]/34.62 << "\\pm" << g2_stat->GetEY()[i]/34.62 << "\\pm" << g2_syst->GetEY()[i]/34.62 << "$";
   cout << "\\\\" << endl;
   cout << "$l^-$";
   for (int i=NBINS-1; i>NBINS/2-1; i--)
      cout << " & $" << gc_syst->GetY()[i]/34.62 << "\\pm" << gc_stat->GetEY()[i]/34.62 << "\\pm" << gc_syst->GetEY()[i]/34.62 << "$";
   cout << "\\\\" << endl;
   cout << "\\hline" << endl;
   cout << "$CT10+EPS09";
   for (int i=NBINS-1; i>NBINS/2-1; i--)
      cout << " & $" << gnuc->GetY()[i]/34.62 << "^{+" << gnuc->GetEYhigh()[i]/34.62 << "}_{-" << gnuc->GetEYlow()[i]/34.62 << "}$";
   cout << "\\\\" << endl;
   cout << "$CT10";
   for (int i=NBINS-1; i>NBINS/2-1; i--)
      cout << " & $" << gnonuc->GetY()[i]/34.62 << "^{+" << gnonuc->GetEYhigh()[i]/34.62 << "}_{-" << gnonuc->GetEYlow()[i]/34.62 << "}$";
   cout << "\\\\" << endl;
   cout << "\\hline" << endl;
   cout << "\\hline" << endl;
   cout << "\\hline" << endl;
   cout << "$\\frac{d\\sigma}{d\\eta} (\\eta)$ &$[0,0.5]$ & $[0.5,1.0]$ & $[1.0,1.5]$ & $[1.5,2.0]$& $[2.0,2.4]$ \\\\" << endl;
   cout << "$\\Pgm^-$";
   for (int i=NBINS/2-1; i>=0; i--)
      cout << " & $" << g1_syst->GetY()[i]/34.62 << "\\pm" << g1_stat->GetEY()[i]/34.62 << "\\pm" << g1_syst->GetEY()[i]/34.62 << "$";
   cout << "\\\\" << endl;
   cout << "$\\Pe^-$";
   for (int i=NBINS/2-1; i>=0; i--)
      cout << " & $" << g2_syst->GetY()[i]/34.62 << "\\pm" << g2_stat->GetEY()[i]/34.62 << "\\pm" << g2_syst->GetEY()[i]/34.62 << "$";
   cout << "\\\\" << endl;
   cout << "$l^-$";
   for (int i=NBINS/2-1; i>=0; i--)
      cout << " & $" << gc_syst->GetY()[i]/34.62 << "\\pm" << gc_stat->GetEY()[i]/34.62 << "\\pm" << gc_syst->GetEY()[i]/34.62 << "$";
   cout << "\\\\" << endl;
   cout << "\\hline" << endl;
   cout << "$CT10+EPS09";
   for (int i=NBINS/2-1; i>=0; i--)
      cout << " & $" << gnuc->GetY()[i]/34.62 << "^{+" << gnuc->GetEYhigh()[i]/34.62 << "}_{-" << gnuc->GetEYlow()[i]/34.62 << "}$";
   cout << "\\\\" << endl;
   cout << "$CT10";
   for (int i=NBINS/2-1; i>=0; i--)
      cout << " & $" << gnonuc->GetY()[i]/34.62 << "^{+" << gnonuc->GetEYhigh()[i]/34.62 << "}_{-" << gnonuc->GetEYlow()[i]/34.62 << "}$";
   cout << "\\\\" << endl;
   cout << "\\hline" << endl;
   cout << "\\end{tabular}}" << endl;
   cout << "\\end{center}" << endl;
   cout << "\\end{table}" << endl;
}
