{
//=========Macro generated from canvas: c1/c1
//=========  (Tue Dec  2 14:50:33 2014) by ROOT version5.34/10
   TCanvas *c1 = new TCanvas("c1", "c1",0,23,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-0.4746836,-0.01962963,2.689873,0.05444444);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(0);
   c1->SetTickx(1);
   c1->SetTicky(1);
   c1->SetLeftMargin(0.15);
   c1->SetRightMargin(0.06);
   c1->SetTopMargin(0.06);
   c1->SetBottomMargin(0.13);
   c1->SetFrameLineColor(0);
   c1->SetFrameBorderMode(0);
   c1->SetFrameLineColor(0);
   c1->SetFrameBorderMode(0);
   
   TGraphAsymmErrors *grae = new TGraphAsymmErrors(5);
   grae->SetName("Graph0");
   grae->SetTitle("Graph");

   Int_t ci;   // for color index setting
   ci = TColor::GetColor("#ffff00");
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0);
   grae->SetPoint(0,0.25,0.00122301);
   grae->SetPointError(0,0.25,0.25,0.0006756227,0.0007680071);
   grae->SetPoint(1,0.75,0.00447037);
   grae->SetPointError(1,0.25,0.25,0.001788515,0.002344638);
   grae->SetPoint(2,1.25,0.0105736);
   grae->SetPointError(2,0.25,0.25,0.002406357,0.002716315);
   grae->SetPoint(3,1.75,0.0216333);
   grae->SetPointError(3,0.25,0.25,0.001954734,0.002083388);
   grae->SetPoint(4,2.2,0.0349209);
   grae->SetPointError(4,0.2,0.2,0.00166497,0.001340596);
   
   TH1F *Graph_Graph9 = new TH1F("Graph_Graph9","Graph",100,0,2.5);
   Graph_Graph9->SetMinimum(-0.01);
   Graph_Graph9->SetMaximum(0.05);
   Graph_Graph9->SetDirectory(0);
   Graph_Graph9->SetStats(0);
   Graph_Graph9->SetLineStyle(0);
   Graph_Graph9->SetMarkerStyle(20);
   Graph_Graph9->GetXaxis()->SetTitle("#eta_{lab}");
   Graph_Graph9->GetXaxis()->SetNdivisions(505);
   Graph_Graph9->GetXaxis()->SetLabelFont(42);
   Graph_Graph9->GetXaxis()->SetLabelOffset(0.01);
   Graph_Graph9->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph9->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph9->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph9->GetXaxis()->SetTitleFont(42);
   Graph_Graph9->GetYaxis()->SetTitle("A_{4}");
   Graph_Graph9->GetYaxis()->SetLabelFont(42);
   Graph_Graph9->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph9->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph9->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph9->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph9->GetYaxis()->SetTitleFont(42);
   Graph_Graph9->GetZaxis()->SetLabelFont(42);
   Graph_Graph9->GetZaxis()->SetLabelSize(0.045);
   Graph_Graph9->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph9);
   
   grae->Draw("a2");
   
   grae = new TGraphAsymmErrors(5);
   grae->SetName("Graph1");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#ffff00");
   grae->SetFillColor(ci);

   ci = TColor::GetColor("#ff0000");
   grae->SetLineColor(ci);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetMarkerSize(0);
   grae->SetPoint(0,0.25,0.00122301);
   grae->SetPointError(0,0.25,0.25,0,0);
   grae->SetPoint(1,0.75,0.00447037);
   grae->SetPointError(1,0.25,0.25,0,0);
   grae->SetPoint(2,1.25,0.0105736);
   grae->SetPointError(2,0.25,0.25,0,0);
   grae->SetPoint(3,1.75,0.0216333);
   grae->SetPointError(3,0.25,0.25,0,0);
   grae->SetPoint(4,2.2,0.0349209);
   grae->SetPointError(4,0.2,0.2,0,0);
   grae->Draw("z");
   
   grae = new TGraphAsymmErrors(5);
   grae->SetName("Graph2");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#009900");
   grae->SetFillColor(ci);
   grae->SetFillStyle(3345);

   ci = TColor::GetColor("#009900");
   grae->SetLineColor(ci);
   grae->SetLineStyle(9);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,0.25,-0.00169038);
   grae->SetPointError(0,0.25,0.25,0.001116978,0.001604194);
   grae->SetPoint(1,0.75,-0.00332414);
   grae->SetPointError(1,0.25,0.25,0.003372264,0.003863969);
   grae->SetPoint(2,1.25,-2.02198e-05);
   grae->SetPointError(2,0.25,0.25,0.005124405,0.004955424);
   grae->SetPoint(3,1.75,0.0103701);
   grae->SetPointError(3,0.25,0.25,0.005712242,0.005574079);
   grae->SetPoint(4,2.2,0.0243959);
   grae->SetPointError(4,0.2,0.2,0.005904851,0.005571725);
   
   TH1F *Graph_Graph10 = new TH1F("Graph_Graph10","Graph",100,0,2.64);
   Graph_Graph10->SetMinimum(-0.01036281);
   Graph_Graph10->SetMaximum(0.03363403);
   Graph_Graph10->SetDirectory(0);
   Graph_Graph10->SetStats(0);
   Graph_Graph10->SetLineStyle(0);
   Graph_Graph10->SetMarkerStyle(20);
   Graph_Graph10->GetXaxis()->SetLabelFont(42);
   Graph_Graph10->GetXaxis()->SetLabelOffset(0.01);
   Graph_Graph10->GetXaxis()->SetLabelSize(0.045);
   Graph_Graph10->GetXaxis()->SetTitleSize(0.055);
   Graph_Graph10->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph10->GetXaxis()->SetTitleFont(42);
   Graph_Graph10->GetYaxis()->SetLabelFont(42);
   Graph_Graph10->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph10->GetYaxis()->SetLabelSize(0.045);
   Graph_Graph10->GetYaxis()->SetTitleSize(0.055);
   Graph_Graph10->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph10->GetYaxis()->SetTitleFont(42);
   Graph_Graph10->GetZaxis()->SetLabelFont(42);
   Graph_Graph10->GetZaxis()->SetLabelSize(0.045);
   Graph_Graph10->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph10);
   
   grae->Draw("2");
   
   grae = new TGraphAsymmErrors(5);
   grae->SetName("Graph3");
   grae->SetTitle("Graph");

   ci = TColor::GetColor("#009900");
   grae->SetFillColor(ci);
   grae->SetFillStyle(3345);

   ci = TColor::GetColor("#009900");
   grae->SetLineColor(ci);
   grae->SetLineStyle(9);
   grae->SetLineWidth(2);
   grae->SetMarkerStyle(20);
   grae->SetPoint(0,0.25,-0.00169038);
   grae->SetPointError(0,0.25,0.25,0,0);
   grae->SetPoint(1,0.75,-0.00332414);
   grae->SetPointError(1,0.25,0.25,0,0);
   grae->SetPoint(2,1.25,-2.02198e-05);
   grae->SetPointError(2,0.25,0.25,0,0);
   grae->SetPoint(3,1.75,0.0103701);
   grae->SetPointError(3,0.25,0.25,0,0);
   grae->SetPoint(4,2.2,0.0243959);
   grae->SetPointError(4,0.2,0.2,0,0);
   grae->Draw("z");
   
   TGraphErrors *gre = new TGraphErrors(5);
   gre->SetName("gA4_exp_1");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetFillStyle(0);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,0.25,0.002897462);
   gre->SetPointError(0,0,0.003226835);
   gre->SetPoint(1,0.75,-0.0006229979);
   gre->SetPointError(1,0,0.002936148);
   gre->SetPoint(2,1.25,-0.00121096);
   gre->SetPointError(2,0,0.003464562);
   gre->SetPoint(3,1.75,0.01196637);
   gre->SetPointError(3,0,0.00396891);
   gre->SetPoint(4,2.2,0.0283991);
   gre->SetPointError(4,0,0.004234413);
   gre->Draw("||");
   
   gre = new TGraphErrors(5);
   gre->SetName("gA4_exp_statonly_1");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,0.25,0.002897462);
   gre->SetPointError(0,0,0.002438023);
   gre->SetPoint(1,0.75,-0.0006229979);
   gre->SetPointError(1,0,0.002404777);
   gre->SetPoint(2,1.25,-0.00121096);
   gre->SetPointError(2,0,0.00250501);
   gre->SetPoint(3,1.75,0.01196637);
   gre->SetPointError(3,0,0.002477774);
   gre->SetPoint(4,2.2,0.0283991);
   gre->SetPointError(4,0,0.002189718);
   
   TH1F *Graph_gA4_exp_statonly_15 = new TH1F("Graph_gA4_exp_statonly_15","Graph",100,0.055,2.395);
   Graph_gA4_exp_statonly_15->SetMinimum(-0.007146448);
   Graph_gA4_exp_statonly_15->SetMaximum(0.03401929);
   Graph_gA4_exp_statonly_15->SetDirectory(0);
   Graph_gA4_exp_statonly_15->SetStats(0);
   Graph_gA4_exp_statonly_15->SetLineStyle(0);
   Graph_gA4_exp_statonly_15->SetMarkerStyle(20);
   Graph_gA4_exp_statonly_15->GetXaxis()->SetLabelFont(42);
   Graph_gA4_exp_statonly_15->GetXaxis()->SetLabelOffset(0.01);
   Graph_gA4_exp_statonly_15->GetXaxis()->SetLabelSize(0.045);
   Graph_gA4_exp_statonly_15->GetXaxis()->SetTitleSize(0.055);
   Graph_gA4_exp_statonly_15->GetXaxis()->SetTitleOffset(1.1);
   Graph_gA4_exp_statonly_15->GetXaxis()->SetTitleFont(42);
   Graph_gA4_exp_statonly_15->GetYaxis()->SetLabelFont(42);
   Graph_gA4_exp_statonly_15->GetYaxis()->SetLabelOffset(0.01);
   Graph_gA4_exp_statonly_15->GetYaxis()->SetLabelSize(0.045);
   Graph_gA4_exp_statonly_15->GetYaxis()->SetTitleSize(0.055);
   Graph_gA4_exp_statonly_15->GetYaxis()->SetTitleOffset(1.4);
   Graph_gA4_exp_statonly_15->GetYaxis()->SetTitleFont(42);
   Graph_gA4_exp_statonly_15->GetZaxis()->SetLabelFont(42);
   Graph_gA4_exp_statonly_15->GetZaxis()->SetLabelSize(0.045);
   Graph_gA4_exp_statonly_15->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_gA4_exp_statonly_15);
   
   gre->Draw("pz");
   
   TPaveText *pt = new TPaveText(0.15,0.925,0.3,0.99,"brNDC");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextAlign(11);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   TText *text = pt->AddText("CMS pPb 34.6 nb^{-1}");
   pt->Draw();
   
   pt = new TPaveText(0.5,0.925,0.955,0.99,"brNDC");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextAlign(31);
   pt->SetTextFont(42);
   pt->SetTextSize(0.04);
   text = pt->AddText("#sqrt{s_{NN}} = 5.02 TeV");
   pt->Draw();
   
   TLegend *leg = new TLegend(0.21,0.62,0.51,0.77,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("gA4_exp_1","Data","p");
   entry->SetLineColor(1);
   entry->SetLineStyle(1);
   entry->SetLineWidth(1);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(20);
   entry->SetMarkerSize(1.2);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph0","CT10","lf");

   ci = TColor::GetColor("#ffff00");
   entry->SetFillColor(ci);
   entry->SetFillStyle(1001);

   ci = TColor::GetColor("#ff0000");
   entry->SetLineColor(ci);
   entry->SetLineStyle(1);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   entry=leg->AddEntry("Graph2","CT10+EPS09","lf");

   ci = TColor::GetColor("#009900");
   entry->SetFillColor(ci);
   entry->SetFillStyle(3345);

   ci = TColor::GetColor("#009900");
   entry->SetLineColor(ci);
   entry->SetLineStyle(9);
   entry->SetLineWidth(2);
   entry->SetMarkerColor(1);
   entry->SetMarkerStyle(21);
   entry->SetMarkerSize(1);
   entry->SetTextFont(42);
   leg->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
