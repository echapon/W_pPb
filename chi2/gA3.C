{
//=========Macro generated from canvas: c1/c1
//=========  (Tue Dec  2 14:50:33 2014) by ROOT version5.34/10
   TCanvas *c1 = new TCanvas("c1", "c1",0,23,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-0.4746836,0.6395062,2.689873,1.874074);
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
   grae->SetPoint(0,0.25,1.01143);
   grae->SetPointError(0,0.25,0.25,0.006224436,0.006999013);
   grae->SetPoint(1,0.75,1.04282);
   grae->SetPointError(1,0.25,0.25,0.0169939,0.02247715);
   grae->SetPoint(2,1.25,1.1075);
   grae->SetPointError(2,0.25,0.25,0.02518055,0.0292036);
   grae->SetPoint(3,1.75,1.25194);
   grae->SetPointError(3,0.25,0.25,0.0265282,0.02983749);
   grae->SetPoint(4,2.2,1.51393);
   grae->SetPointError(4,0.2,0.2,0.03276547,0.03201791);
   
   TH1F *Graph_Graph1 = new TH1F("Graph_Graph1","Graph",100,0,2.5);
   Graph_Graph1->SetMinimum(0.8);
   Graph_Graph1->SetMaximum(1.8);
   Graph_Graph1->SetDirectory(0);
   Graph_Graph1->SetStats(0);
   Graph_Graph1->SetLineStyle(0);
   Graph_Graph1->SetMarkerStyle(20);
   Graph_Graph1->GetXaxis()->SetTitle("#eta_{lab}");
   Graph_Graph1->GetXaxis()->SetNdivisions(505);
   Graph_Graph1->GetXaxis()->SetLabelFont(42);
   Graph_Graph1->GetXaxis()->SetLabelOffset(0.01);
   Graph_Graph1->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph1->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph1->GetXaxis()->SetTitleFont(42);
   Graph_Graph1->GetYaxis()->SetTitle("N(+#eta_{lab})/N(-#eta_{lab})");
   Graph_Graph1->GetYaxis()->SetLabelFont(42);
   Graph_Graph1->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph1->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph1->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph1->GetYaxis()->SetTitleFont(42);
   Graph_Graph1->GetZaxis()->SetLabelFont(42);
   Graph_Graph1->GetZaxis()->SetLabelSize(0.045);
   Graph_Graph1->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph1);
   
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
   grae->SetPoint(0,0.25,1.01143);
   grae->SetPointError(0,0.25,0.25,0,0);
   grae->SetPoint(1,0.75,1.04282);
   grae->SetPointError(1,0.25,0.25,0,0);
   grae->SetPoint(2,1.25,1.1075);
   grae->SetPointError(2,0.25,0.25,0,0);
   grae->SetPoint(3,1.75,1.25194);
   grae->SetPointError(3,0.25,0.25,0,0);
   grae->SetPoint(4,2.2,1.51393);
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
   grae->SetPoint(0,0.25,0.984845);
   grae->SetPointError(0,0.25,0.25,0.009826457,0.01433731);
   grae->SetPoint(1,0.75,0.969881);
   grae->SetPointError(1,0.25,0.25,0.02972037,0.03470488);
   grae->SetPoint(2,1.25,0.999806);
   grae->SetPointError(2,0.25,0.25,0.04877509,0.04803818);
   grae->SetPoint(3,1.75,1.11566);
   grae->SetPointError(3,0.25,0.25,0.06597888,0.06589705);
   grae->SetPoint(4,2.2,1.35065);
   grae->SetPointError(4,0.2,0.2,0.09201274,0.08930102);
   
   TH1F *Graph_Graph2 = new TH1F("Graph_Graph2","Graph",100,0,2.64);
   Graph_Graph2->SetMinimum(0.8901816);
   Graph_Graph2->SetMaximum(1.48993);
   Graph_Graph2->SetDirectory(0);
   Graph_Graph2->SetStats(0);
   Graph_Graph2->SetLineStyle(0);
   Graph_Graph2->SetMarkerStyle(20);
   Graph_Graph2->GetXaxis()->SetLabelFont(42);
   Graph_Graph2->GetXaxis()->SetLabelOffset(0.01);
   Graph_Graph2->GetXaxis()->SetLabelSize(0.045);
   Graph_Graph2->GetXaxis()->SetTitleSize(0.055);
   Graph_Graph2->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph2->GetXaxis()->SetTitleFont(42);
   Graph_Graph2->GetYaxis()->SetLabelFont(42);
   Graph_Graph2->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph2->GetYaxis()->SetLabelSize(0.045);
   Graph_Graph2->GetYaxis()->SetTitleSize(0.055);
   Graph_Graph2->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph2->GetYaxis()->SetTitleFont(42);
   Graph_Graph2->GetZaxis()->SetLabelFont(42);
   Graph_Graph2->GetZaxis()->SetLabelSize(0.045);
   Graph_Graph2->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph2);
   
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
   grae->SetPoint(0,0.25,0.984845);
   grae->SetPointError(0,0.25,0.25,0,0);
   grae->SetPoint(1,0.75,0.969881);
   grae->SetPointError(1,0.25,0.25,0,0);
   grae->SetPoint(2,1.25,0.999806);
   grae->SetPointError(2,0.25,0.25,0,0);
   grae->SetPoint(3,1.75,1.11566);
   grae->SetPointError(3,0.25,0.25,0,0);
   grae->SetPoint(4,2.2,1.35065);
   grae->SetPointError(4,0.2,0.2,0,0);
   grae->Draw("z");
   
   TGraphErrors *gre = new TGraphErrors(5);
   gre->SetName("gA3_exp_1");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetFillStyle(0);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,0.25,1.01792);
   gre->SetPointError(0,0,0.02869542);
   gre->SetPoint(1,0.75,0.9940261);
   gre->SetPointError(1,0,0.02621272);
   gre->SetPoint(2,1.25,0.9881206);
   gre->SetPointError(2,0,0.03210871);
   gre->SetPoint(3,1.75,1.13314);
   gre->SetPointError(3,0,0.04673086);
   gre->SetPoint(4,2.2,1.491612);
   gre->SetPointError(4,0,0.08419942);
   gre->Draw("||");
   
   gre = new TGraphErrors(5);
   gre->SetName("gA3_exp_statonly_1");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,0.25,1.01792);
   gre->SetPointError(0,0,0.02157377);
   gre->SetPoint(1,0.75,0.9940261);
   gre->SetPointError(1,0,0.02129124);
   gre->SetPoint(2,1.25,0.9881206);
   gre->SetPointError(2,0,0.02331959);
   gre->SetPoint(3,1.75,1.13314);
   gre->SetPointError(3,0,0.02951123);
   gre->SetPoint(4,2.2,1.491612);
   gre->SetPointError(4,0,0.04584313);
   
   TH1F *Graph_gA3_exp_statonly_11 = new TH1F("Graph_gA3_exp_statonly_11","Graph",100,0.055,2.395);
   Graph_gA3_exp_statonly_11->SetMinimum(0.9075355);
   Graph_gA3_exp_statonly_11->SetMaximum(1.594721);
   Graph_gA3_exp_statonly_11->SetDirectory(0);
   Graph_gA3_exp_statonly_11->SetStats(0);
   Graph_gA3_exp_statonly_11->SetLineStyle(0);
   Graph_gA3_exp_statonly_11->SetMarkerStyle(20);
   Graph_gA3_exp_statonly_11->GetXaxis()->SetLabelFont(42);
   Graph_gA3_exp_statonly_11->GetXaxis()->SetLabelOffset(0.01);
   Graph_gA3_exp_statonly_11->GetXaxis()->SetLabelSize(0.045);
   Graph_gA3_exp_statonly_11->GetXaxis()->SetTitleSize(0.055);
   Graph_gA3_exp_statonly_11->GetXaxis()->SetTitleOffset(1.1);
   Graph_gA3_exp_statonly_11->GetXaxis()->SetTitleFont(42);
   Graph_gA3_exp_statonly_11->GetYaxis()->SetLabelFont(42);
   Graph_gA3_exp_statonly_11->GetYaxis()->SetLabelOffset(0.01);
   Graph_gA3_exp_statonly_11->GetYaxis()->SetLabelSize(0.045);
   Graph_gA3_exp_statonly_11->GetYaxis()->SetTitleSize(0.055);
   Graph_gA3_exp_statonly_11->GetYaxis()->SetTitleOffset(1.4);
   Graph_gA3_exp_statonly_11->GetYaxis()->SetTitleFont(42);
   Graph_gA3_exp_statonly_11->GetZaxis()->SetLabelFont(42);
   Graph_gA3_exp_statonly_11->GetZaxis()->SetLabelSize(0.045);
   Graph_gA3_exp_statonly_11->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_gA3_exp_statonly_11);
   
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
   
   TLegend *leg = new TLegend(0.21,0.45,0.51,0.6,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("gA3_exp_1","Data","p");
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
   
   pt = new TPaveText(0.25,1.5,0.75,1.6,"br");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetTextAlign(13);
   pt->SetTextFont(42);
   pt->SetTextSize(0.05);
   text = pt->AddText("W #rightarrow #font[12]{l} + #nu");
   pt->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
