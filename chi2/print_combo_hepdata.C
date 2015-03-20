#define NBINS 10
#include "../plot_graph.C"

void print_combo_hepdata(const char *file, const char* gbasename, int channel=1, double norm=1./34.62, int prec=1, int nbins=10)
{
   TFile *fs = new TFile(file);

   TGraphErrors *g1_syst = (TGraphErrors*) fs->Get(Form("%s_exp_%i",gbasename,channel));
   TGraphErrors *g1_stat = (TGraphErrors*) fs->Get(Form("%s_exp_statonly_%i",gbasename,channel));

   // TGraphAsymmErrors *gnuc = theory_errors(gbasename,1,NBINS,"EPS09");
   // TGraphAsymmErrors *gnonuc = theory_errors(gbasename,1,NBINS,"CT10");


   cout.setf(std::ios::fixed);
   cout << setprecision(prec);


   if (nbins==NBINS)
      for (int i=NBINS-1; i>=0; i--)
      {
         cout << -g1_syst->GetX()[i]-g1_syst->GetEX()[i] << " TO " << -g1_syst->GetX()[i]+g1_syst->GetEX()[i] << "; ";
         cout << g1_syst->GetY()[i]*norm << " +- " << g1_stat->GetEY()[i]*norm;
         cout << "(DSYS=" << g1_syst->GetEY()[i]*norm << ",DSYS=" << 0.035*g1_syst->GetY()[i]*norm << ":lumi);" << endl;
      }

   else
      for (int i=0; i<nbins; i++)
      {
         cout << g1_syst->GetX()[i]-g1_syst->GetEX()[i] << " TO " << g1_syst->GetX()[i]+g1_syst->GetEX()[i] << "; ";
         cout << g1_syst->GetY()[i]*norm << " +- " << g1_stat->GetEY()[i]*norm;
         cout << "(DSYS=" << g1_syst->GetEY()[i]*norm << ");" << endl;
      }

}
