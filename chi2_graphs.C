void chi2_graphs(const char* fname_cteq="graph.root", const char* fname_eps="graph.root", const char *gbasename="gyields", int channel_number=1, double lumierr=0)
{
   TFile *fcteq = new TFile(fname_cteq);
   TFile *feps = new TFile(fname_eps);

   TString gname_th = Form("%s_th_%i",gbasename,channel_number);
   TString gname_stat = Form("%s_exp_statonly_%i",gbasename,channel_number);
   TString gname_syst = Form("%s_exp_%i",gbasename,channel_number);

   TGraphErrors *gth_cteq = (TGraphErrors*) fcteq->Get(gname_th);
   TGraphErrors *gth_eps = (TGraphErrors*) feps->Get(gname_th);
   TGraphErrors *gexp_stat = (TGraphErrors*) fcteq->Get(gname_stat);
   TGraphErrors *gexp_syst = (TGraphErrors*) fcteq->Get(gname_syst);

   double chi2_cteq=0, chi2_eps=0;
   int nbins = gth_cteq->GetN();

   // build the covariance matrix
   TMatrixDSym tcov_cteq(nbins), tcov_eps(nbins);
   for (int i=0; i<nbins; i++)
   {
      double therr_cteq = gth_cteq->GetEY()[i];
      double therr_eps = gth_eps->GetEY()[i];
      double stat = gexp_stat->GetEY()[i];
      double syst = gexp_syst->GetEY()[i];

      for (int j=0; j<=i; j++)
      {
         if (i==j) 
         {
            tcov_cteq[i][j] = therr_cteq*therr_cteq+stat*stat+syst*syst;
            tcov_eps[i][j] = therr_eps*therr_eps+stat*stat+syst*syst;
         }
         else
         {
            tcov_cteq[i][j] = lumierr*gth_cteq->GetEY()[i]*lumierr*gth_cteq->GetEY()[j];
            tcov_eps[i][j] = lumierr*gth_eps->GetEY()[i]*lumierr*gth_eps->GetEY()[j];
         }
      }
   }

   TMatrixDSym tinvcov_cteq(nbins), tinvcov_eps(nbins);
   tinvcov_cteq = tcov_cteq.Invert();
   tinvcov_eps = tcov_eps.Invert();

   for (int i=0; i<nbins; i++)
   {
      double pred_cteq = gth_cteq->GetY()[i];
      double therr_cteq = gth_cteq->GetEY()[i];
      double pred_eps = gth_eps->GetY()[i];
      double therr_eps = gth_eps->GetEY()[i];
      double exp = gexp_stat->GetY()[i];
      double stat = gexp_stat->GetEY()[i];
      double syst = gexp_syst->GetEY()[i];

      chi2_cteq += pow(pred_cteq-exp,2)/(therr_cteq*therr_cteq+stat*stat+syst*syst);
      chi2_eps += pow(pred_eps-exp,2)/(therr_eps*therr_eps+stat*stat+syst*syst);
   }

   cout << "Old school chi2 (no correlations):" << endl;
   cout << "CTEQ: " << chi2_cteq << "/" << nbins << " = " << chi2_cteq/nbins << endl;
   cout << "EPS: " << chi2_eps << "/" << nbins << " = " << chi2_eps/nbins << endl;

   chi2_cteq=0; chi2_eps=0;
   for (int i=0; i<nbins; i++)
      for (int j=0; j<nbins; j++)
      {
         chi2_cteq += (gth_cteq->GetY()[i]-gexp_stat->GetY()[i])*tinvcov_cteq[i][j]*(gth_cteq->GetY()[j]-gexp_stat->GetY()[j]);
         chi2_eps += (gth_eps->GetY()[i]-gexp_stat->GetY()[i])*tinvcov_eps[i][j]*(gth_eps->GetY()[j]-gexp_stat->GetY()[j]);
      }

   cout << "Cool chi2:" << endl;
   cout << "CTEQ: " << chi2_cteq << "/" << nbins << " = " << chi2_cteq/nbins << endl;
   cout << "EPS: " << chi2_eps << "/" << nbins << " = " << chi2_eps/nbins << endl;
}

double chi2(TGraphErrors *g1_stat, TGraphErrors *g1_syst, TGraphErrors *g2_stat, TGraphErrors *g2_syst)
{
   double val=0.;

   for (int i=0; i<g1_stat->GetN(); i++)
   {
      double dval = pow(g1_stat->GetY()[i]-g2_stat->GetY()[i],2)/(pow(g1_stat->GetEY()[i],2)+pow(g1_syst->GetEY()[i],2)+pow(g2_stat->GetEY()[i],2)+pow(g2_syst->GetEY()[i],2));
      cout << "bin " << i << ": (" << g1_stat->GetY()[i]/34.62 << "-" << g2_stat->GetY()[i]/34.62 << ")^2/(" << g1_stat->GetEY()[i]/34.62 << "+" << g1_syst->GetEY()[i]/34.62 << "+" << g2_stat->GetEY()[i]/34.62 << "+" << g2_syst->GetEY()[i]/34.62 << ") = " <<  dval << endl;
      val += dval;
   }

   cout << "chi2 = " << val << " (prob = " << TMath::Prob(val,g1_stat->GetN()) << " -> " << TMath::ErfInverse(1-(0.5*(TMath::Prob(val,g1_stat->GetN())))) << " sigma)" << endl;
   return val;
}

double chi2_sum(TGraphErrors *g1_stat, TGraphErrors *g1_syst, TGraphErrors *g2_stat, TGraphErrors *g2_syst, double scale = 1./34.62)
{
   double s1=0,s2=0,st1=0,st2=0,sy1=0,sy2=0;

   for (int i=0; i<g1_stat->GetN(); i++)
   {
      // double deta = g2_stat->GetEX()[i]*2.;
      double deta = 0.5;
      s1 += g1_stat->GetY()[i]*scale*deta;
      s2 += g2_stat->GetY()[i]*scale*deta;
      st1 += pow(g1_stat->GetEY()[i]*scale*deta,2);
      st2 += pow(g2_stat->GetEY()[i]*scale*deta,2);
      sy1 += pow(g1_syst->GetEY()[i]*scale*deta,2);
      sy2 += pow(g2_syst->GetEY()[i]*scale*deta,2);
   }

   st1 = sqrt(st1);
   st2 = sqrt(st2);
   sy1 = sqrt(sy1);
   sy2 = sqrt(sy2);

   cout << "ch1: " << s1 << " +- " << st1 << " +- " << sy1 << endl;
   cout << "ch2: " << s2 << " +- " << st2 << " +- " << sy2 << endl;

   double val = pow(s1-s2,2)/(st1*st1+st2*st2+sy1*sy1+sy2*sy2);

   cout << "chi2 = " << val << " (prob = " << TMath::Prob(val,1) << " -> " << sqrt(2)*TMath::ErfcInverse(TMath::Prob(val,1)) << " sigma)" << endl;

   return val;
}

double chi2_sum(const char* fname="graph.root", const char* gbasename="gyieldsp", double scale = 1./34.62)
{
   TFile *f = new TFile(fname);
   TString gname_stat = Form("%s_exp_statonly_1",gbasename);
   TString gname_syst = Form("%s_exp_1",gbasename);
   TGraphErrors *gexp_stat1 = (TGraphErrors*) f->Get(gname_stat);
   TGraphErrors *gexp_syst1 = (TGraphErrors*) f->Get(gname_syst);
   gname_stat = Form("%s_exp_statonly_2",gbasename);
   gname_syst = Form("%s_exp_2",gbasename);
   TGraphErrors *gexp_stat2 = (TGraphErrors*) f->Get(gname_stat);
   TGraphErrors *gexp_syst2 = (TGraphErrors*) f->Get(gname_syst);

   return chi2_sum(gexp_stat1, gexp_syst1, gexp_stat2, gexp_syst2, scale);
}

double chi2(const char* fname="graph.root", const char* gbasename="gyieldsp")
{
   TFile *f = new TFile(fname);
   TString gname_stat = Form("%s_exp_statonly_1",gbasename);
   TString gname_syst = Form("%s_exp_1",gbasename);
   TGraphErrors *gexp_stat1 = (TGraphErrors*) f->Get(gname_stat);
   TGraphErrors *gexp_syst1 = (TGraphErrors*) f->Get(gname_syst);
   gname_stat = Form("%s_exp_statonly_2",gbasename);
   gname_syst = Form("%s_exp_2",gbasename);
   TGraphErrors *gexp_stat2 = (TGraphErrors*) f->Get(gname_stat);
   TGraphErrors *gexp_syst2 = (TGraphErrors*) f->Get(gname_syst);

   return chi2(gexp_stat1, gexp_syst1, gexp_stat2, gexp_syst2);
}
