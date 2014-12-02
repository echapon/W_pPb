#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

#include "TFile.h"
#include "TGraphErrors.h"

// #include "stattest_gaus.h"
#include "statchannel.h"

using namespace std;

int main(int argc, const char** argv)
{
   if (argc!=4 && argc!=7)
   {
      cout << "Usage:" << endl;
      cout << argv[0] << " n_etabins1 channel1_plus.txt channel1_minus.txt n_etabins2 channel2_plus.txt chennael2_minus.txt" << endl;
      return -1;
   }

   bool singlechan = argc==4;

   statchannel chan1p, chan1m, chan2p, chan2m;
   int nbins1, nbins2;
   int nsyst1, nsyst2;

   TFile *f = new TFile("graph.root","RECREATE");
   f->cd();

   // channel 1
   nbins1 = atoi(argv[1]);
   chan1p.read(nbins1,argv[2]);
   chan1m.read(nbins1,argv[3]);
   // chan1p.print();
   // chan1m.print();
   nsyst1 = chan1p.get_eff().size();
   if (chan1m.get_eff().size() != nsyst1)
   {
      cout << "Error, " << argv[2] << " and " << argv[3] << " have a different number of efficiencies" << endl;
      return -1;
   }

   // create the graphs
   TGraphErrors *gyields_th_1 = statchannel::graph(chan1p, chan1m, YIELDS_PM, true); gyields_th_1->SetName("gyields_th_1"); gyields_th_1->Write();
   TGraphErrors *gyields_exp_1 = statchannel::graph(chan1p, chan1m, YIELDS_PM, false, true);  gyields_exp_1->SetName("gyields_exp_1"); gyields_exp_1->Write();
   TGraphErrors *gyields_exp_statonly_1 = statchannel::graph(chan1p, chan1m, YIELDS_PM, false, false);  gyields_exp_statonly_1->SetName("gyields_exp_statonly_1"); gyields_exp_statonly_1->Write();
   TGraphErrors *gyieldsp_th_1 = statchannel::graph(chan1p, chan1m, YIELDS, true);  gyieldsp_th_1->SetName("gyieldsp_th_1"); gyieldsp_th_1->Write();
   TGraphErrors *gyieldsp_exp_1 = statchannel::graph(chan1p, chan1m, YIELDS, false, true);  gyieldsp_exp_1->SetName("gyieldsp_exp_1"); gyieldsp_exp_1->Write();
   TGraphErrors *gyieldsp_exp_statonly_1 = statchannel::graph(chan1p, chan1m, YIELDS, false, false);  gyieldsp_exp_statonly_1->SetName("gyieldsp_exp_statonly_1"); gyieldsp_exp_statonly_1->Write();
   TGraphErrors *gyieldsm_th_1 = statchannel::graph(chan1m, chan1p, YIELDS, true);  gyieldsm_th_1->SetName("gyieldsm_th_1"); gyieldsm_th_1->Write();
   TGraphErrors *gyieldsm_exp_1 = statchannel::graph(chan1m, chan1p, YIELDS, false, true);  gyieldsm_exp_1->SetName("gyieldsm_exp_1"); gyieldsm_exp_1->Write();
   TGraphErrors *gyieldsm_exp_statonly_1 = statchannel::graph(chan1m, chan1p, YIELDS, false, false);  gyieldsm_exp_statonly_1->SetName("gyieldsm_exp_statonly_1"); gyieldsm_exp_statonly_1->Write();
   TGraphErrors *gch_th_1 = statchannel::graph(chan1p, chan1m, CH, true); gch_th_1->SetName("gch_th_1"); gch_th_1->Write();
   TGraphErrors *gch_exp_1 = statchannel::graph(chan1p, chan1m, CH, false, true); gch_exp_1->SetName("gch_exp_1"); gch_exp_1->Write();
   TGraphErrors *gch_exp_statonly_1 = statchannel::graph(chan1p, chan1m, CH, false, false); gch_exp_statonly_1->SetName("gch_exp_statonly_1"); gch_exp_statonly_1->Write();
   TGraphErrors *gA1p_th_1 = statchannel::graph(chan1p, chan1m, A1p, true); gA1p_th_1->SetName("gA1p_th_1"); gA1p_th_1->Write();
   TGraphErrors *gA1p_exp_1 = statchannel::graph(chan1p, chan1m, A1p, false, true); gA1p_exp_1->SetName("gA1p_exp_1"); gA1p_exp_1->Write();
   TGraphErrors *gA1p_exp_statonly_1 = statchannel::graph(chan1p, chan1m, A1p, false, false); gA1p_exp_statonly_1->SetName("gA1p_exp_statonly_1"); gA1p_exp_statonly_1->Write();
   TGraphErrors *gA1m_th_1 = statchannel::graph(chan1p, chan1m, A1m, true); gA1m_th_1->SetName("gA1m_th_1"); gA1m_th_1->Write();
   TGraphErrors *gA1m_exp_1 = statchannel::graph(chan1p, chan1m, A1m, false, true); gA1m_exp_1->SetName("gA1m_exp_1"); gA1m_exp_1->Write();
   TGraphErrors *gA1m_exp_statonly_1 = statchannel::graph(chan1p, chan1m, A1m, false, false); gA1m_exp_statonly_1->SetName("gA1m_exp_statonly_1"); gA1m_exp_statonly_1->Write();
   TGraphErrors *gA2_th_1 = statchannel::graph(chan1p, chan1m, A2, true); gA2_th_1->SetName("gA2_th_1"); gA2_th_1->Write();
   TGraphErrors *gA2_exp_1 = statchannel::graph(chan1p, chan1m, A2, false, true); gA2_exp_1->SetName("gA2_exp_1"); gA2_exp_1->Write();
   TGraphErrors *gA2_exp_statonly_1 = statchannel::graph(chan1p, chan1m, A2, false, false); gA2_exp_statonly_1->SetName("gA2_exp_statonly_1"); gA2_exp_statonly_1->Write();
   TGraphErrors *gA3_th_1 = statchannel::graph(chan1p, chan1m, A3, true); gA3_th_1->SetName("gA3_th_1"); gA3_th_1->Write();
   TGraphErrors *gA3_exp_1 = statchannel::graph(chan1p, chan1m, A3, false, true); gA3_exp_1->SetName("gA3_exp_1"); gA3_exp_1->Write();
   TGraphErrors *gA3_exp_statonly_1 = statchannel::graph(chan1p, chan1m, A3, false, false); gA3_exp_statonly_1->SetName("gA3_exp_statonly_1"); gA3_exp_statonly_1->Write();
   TGraphErrors *gA4_th_1 = statchannel::graph(chan1p, chan1m, A4, true); gA4_th_1->SetName("gA4_th_1"); gA4_th_1->Write();
   TGraphErrors *gA4_exp_1 = statchannel::graph(chan1p, chan1m, A4, false, true); gA4_exp_1->SetName("gA4_exp_1"); gA4_exp_1->Write();
   TGraphErrors *gA4_exp_statonly_1 = statchannel::graph(chan1p, chan1m, A4, false, false); gA4_exp_statonly_1->SetName("gA4_exp_statonly_1"); gA4_exp_statonly_1->Write();


   // channel 2
   if (!singlechan)
   {
      nbins2 = atoi(argv[4]);
      chan2p.read(nbins2,argv[5]);
      chan2m.read(nbins2,argv[6]);
      // chan2p.print();
      // chan2m.print();
      nsyst2 = chan2p.get_eff().size();
      if (chan2m.get_eff().size() != nsyst2)
      {
         cout << "Error, " << argv[5] << " and " << argv[6] << " have a different number of efficiencies" << endl;
         return -1;
      }

      // create the graphs
      TGraphErrors *gyields_th_2 = statchannel::graph(chan2p, chan2m, YIELDS_PM, true); gyields_th_2->SetName("gyields_th_2"); gyields_th_2->Write();
      TGraphErrors *gyields_exp_2 = statchannel::graph(chan2p, chan2m, YIELDS_PM, false, true);  gyields_exp_2->SetName("gyields_exp_2"); gyields_exp_2->Write();
      TGraphErrors *gyields_exp_statonly_2 = statchannel::graph(chan2p, chan2m, YIELDS_PM, false, false);  gyields_exp_statonly_2->SetName("gyields_exp_statonly_2"); gyields_exp_statonly_2->Write();
      TGraphErrors *gyieldsp_th_2 = statchannel::graph(chan2p, chan2m, YIELDS, true);  gyieldsp_th_2->SetName("gyieldsp_th_2"); gyieldsp_th_2->Write();
      TGraphErrors *gyieldsp_exp_2 = statchannel::graph(chan2p, chan2m, YIELDS, false, true);  gyieldsp_exp_2->SetName("gyieldsp_exp_2"); gyieldsp_exp_2->Write();
      TGraphErrors *gyieldsp_exp_statonly_2 = statchannel::graph(chan2p, chan2m, YIELDS, false, false);  gyieldsp_exp_statonly_2->SetName("gyieldsp_exp_statonly_2"); gyieldsp_exp_statonly_2->Write();
      TGraphErrors *gyieldsm_th_2 = statchannel::graph(chan2m, chan2p, YIELDS, true);  gyieldsm_th_2->SetName("gyieldsm_th_2"); gyieldsm_th_2->Write();
      TGraphErrors *gyieldsm_exp_2 = statchannel::graph(chan2m, chan2p, YIELDS, false, true);  gyieldsm_exp_2->SetName("gyieldsm_exp_2"); gyieldsm_exp_2->Write();
      TGraphErrors *gyieldsm_exp_statonly_2 = statchannel::graph(chan2m, chan2p, YIELDS, false, false);  gyieldsm_exp_statonly_2->SetName("gyieldsm_exp_statonly_2"); gyieldsm_exp_statonly_2->Write();
      TGraphErrors *gch_th_2 = statchannel::graph(chan2p, chan2m, CH, true); gch_th_2->SetName("gch_th_2"); gch_th_2->Write();
      TGraphErrors *gch_exp_2 = statchannel::graph(chan2p, chan2m, CH, false, true); gch_exp_2->SetName("gch_exp_2"); gch_exp_2->Write();
      TGraphErrors *gch_exp_statonly_2 = statchannel::graph(chan2p, chan2m, CH, false, false); gch_exp_statonly_2->SetName("gch_exp_statonly_2"); gch_exp_statonly_2->Write();
      TGraphErrors *gA1p_th_2 = statchannel::graph(chan2p, chan2m, A1p, true); gA1p_th_2->SetName("gA1p_th_2"); gA1p_th_2->Write();
      TGraphErrors *gA1p_exp_2 = statchannel::graph(chan2p, chan2m, A1p, false, true); gA1p_exp_2->SetName("gA1p_exp_2"); gA1p_exp_2->Write();
      TGraphErrors *gA1p_exp_statonly_2 = statchannel::graph(chan2p, chan2m, A1p, false, false); gA1p_exp_statonly_2->SetName("gA1p_exp_statonly_2"); gA1p_exp_statonly_2->Write();
      TGraphErrors *gA1m_th_2 = statchannel::graph(chan2p, chan2m, A1m, true); gA1m_th_2->SetName("gA1m_th_2"); gA1m_th_2->Write();
      TGraphErrors *gA1m_exp_2 = statchannel::graph(chan2p, chan2m, A1m, false, true); gA1m_exp_2->SetName("gA1m_exp_2"); gA1m_exp_2->Write();
      TGraphErrors *gA1m_exp_statonly_2 = statchannel::graph(chan2p, chan2m, A1m, false, false); gA1m_exp_statonly_2->SetName("gA1m_exp_statonly_2"); gA1m_exp_statonly_2->Write();
      TGraphErrors *gA2_th_2 = statchannel::graph(chan2p, chan2m, A2, true); gA2_th_2->SetName("gA2_th_2"); gA2_th_2->Write();
      TGraphErrors *gA2_exp_2 = statchannel::graph(chan2p, chan2m, A2, false, true); gA2_exp_2->SetName("gA2_exp_2"); gA2_exp_2->Write();
      TGraphErrors *gA2_exp_statonly_2 = statchannel::graph(chan2p, chan2m, A2, false, false); gA2_exp_statonly_2->SetName("gA2_exp_statonly_2"); gA2_exp_statonly_2->Write();
      TGraphErrors *gA3_th_2 = statchannel::graph(chan2p, chan2m, A3, true); gA3_th_2->SetName("gA3_th_2"); gA3_th_2->Write();
      TGraphErrors *gA3_exp_2 = statchannel::graph(chan2p, chan2m, A3, false, true); gA3_exp_2->SetName("gA3_exp_2"); gA3_exp_2->Write();
      TGraphErrors *gA3_exp_statonly_2 = statchannel::graph(chan2p, chan2m, A3, false, false); gA3_exp_statonly_2->SetName("gA3_exp_statonly_2"); gA3_exp_statonly_2->Write();
      TGraphErrors *gA4_th_2 = statchannel::graph(chan2p, chan2m, A4, true); gA4_th_2->SetName("gA4_th_2"); gA4_th_2->Write();
      TGraphErrors *gA4_exp_2 = statchannel::graph(chan2p, chan2m, A4, false, true); gA4_exp_2->SetName("gA4_exp_2"); gA4_exp_2->Write();
      TGraphErrors *gA4_exp_statonly_2 = statchannel::graph(chan2p, chan2m, A4, false, false); gA4_exp_statonly_2->SetName("gA4_exp_statonly_2"); gA4_exp_statonly_2->Write();
   }

   // // yields per channel
   // // make the stattest objects
   // stattest_gaus statg1p(nbins1,nsyst1), statg1m(nbins1,nsyst1);
   // stattest_gaus statg2p(nbins2,nsyst2), statg2m(nbins2,nsyst2);
   // statg1p.set_data(chan1p.get_yields());
   // statg1p.set_data_stat(chan1p.get_staterr());
   // statg1p.set_pred(chan1p.get_th());
   // statg1m.set_data(chan1m.get_yields());
   // statg1m.set_data_stat(chan1m.get_staterr());
   // statg1m.set_pred(chan1m.get_th());

   f->Write();
   f->Close();
}

// NB
// objets stat à faire 
// (corréler les syst selon leur signe), ne pas oublier la lumi
// -> 1 par asymétrie
// -> 1 par canal (1p, 1m, 2p, 2m)
// -> 1p+1m, 2p+2m
// -> 1p+1m+2p+2m
//
// automatiser un peu tout ca
// -> une classe par canal
// -> une fonction qui prend l'objet stat, calcule le chi2 et fait la distrib de PE
//
// Creation des graphes
// - fct de statchannel
// - une fct par type de graphe (yields, C, A1p, A1m, A2, A3, A4)
// - retourne le graphe
