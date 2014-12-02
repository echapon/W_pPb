{
//=========Macro generated from canvas: c1/c1
//=========  (Tue Dec  2 14:50:33 2014) by ROOT version5.34/10
   TCanvas *c1 = new TCanvas("c1", "c1",0,23,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-0.4746836,0.4148148,2.689873,3.377778);
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
   grae->SetPoint(0,0.25,1.0283);
   grae->SetPointError(0,0.25,0.25,0.008340636,0.006600096);
   grae->SetPoint(1,0.75,1.1045);
   grae->SetPointError(1,0.25,0.25,0.01759299,0.02525716);
   grae->SetPoint(2,1.25,1.25445);
   grae->SetPointError(2,0.25,0.25,0.02682874,0.03849533);
   grae->SetPoint(3,1.75,1.60653);
   grae->SetPointError(3,0.25,0.25,0.04307044,0.05198622);
   grae->SetPoint(4,2.2,2.38443);
   grae->SetPointError(4,0.2,0.2,0.07661457,0.07231854);
   
   TH1F *Graph_Graph5 = new TH1F("Graph_Graph5","Graph",100,0,2.5);
   Graph_Graph5->SetMinimum(0.8);
   Graph_Graph5->SetMaximum(3.2);
   Graph_Graph5->SetDirectory(0);
   Graph_Graph5->SetStats(0);
   Graph_Graph5->SetLineStyle(0);
   Graph_Graph5->SetMarkerStyle(20);
   Graph_Graph5->GetXaxis()->SetTitle("#eta_{lab}");
   Graph_Graph5->GetXaxis()->SetNdivisions(505);
   Graph_Graph5->GetXaxis()->SetLabelFont(42);
   Graph_Graph5->GetXaxis()->SetLabelOffset(0.01);
   Graph_Graph5->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph5->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph5->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph5->GetXaxis()->SetTitleFont(42);
   Graph_Graph5->GetYaxis()->SetTitle("N^{+}(+#eta_{lab})/N^{+}(-#eta_{lab})");
   Graph_Graph5->GetYaxis()->SetLabelFont(42);
   Graph_Graph5->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph5->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph5->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph5->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph5->GetYaxis()->SetTitleFont(42);
   Graph_Graph5->GetZaxis()->SetLabelFont(42);
   Graph_Graph5->GetZaxis()->SetLabelSize(0.045);
   Graph_Graph5->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph5);
   
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
   grae->SetPoint(0,0.25,1.0283);
   grae->SetPointError(0,0.25,0.25,0,0);
   grae->SetPoint(1,0.75,1.1045);
   grae->SetPointError(1,0.25,0.25,0,0);
   grae->SetPoint(2,1.25,1.25445);
   grae->SetPointError(2,0.25,0.25,0,0);
   grae->SetPoint(3,1.75,1.60653);
   grae->SetPointError(3,0.25,0.25,0,0);
   grae->SetPoint(4,2.2,2.38443);
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
   grae->SetPoint(0,0.25,0.999177);
   grae->SetPointError(0,0.25,0.25,0.01123775,0.01533736);
   grae->SetPoint(1,0.75,1.02266);
   grae->SetPointError(1,0.25,0.25,0.03336155,0.0398917);
   grae->SetPoint(2,1.25,1.1311);
   grae->SetPointError(2,0.25,0.25,0.05871014,0.05983429);
   grae->SetPoint(3,1.75,1.44437);
   grae->SetPointError(3,0.25,0.25,0.09367862,0.09323451);
   grae->SetPoint(4,2.2,2.17478);
   grae->SetPointError(4,0.2,0.2,0.1643659,0.158067);
   
   TH1F *Graph_Graph6 = new TH1F("Graph_Graph6","Graph",100,0,2.64);
   Graph_Graph6->SetMinimum(0.8534485);
   Graph_Graph6->SetMaximum(2.467338);
   Graph_Graph6->SetDirectory(0);
   Graph_Graph6->SetStats(0);
   Graph_Graph6->SetLineStyle(0);
   Graph_Graph6->SetMarkerStyle(20);
   Graph_Graph6->GetXaxis()->SetLabelFont(42);
   Graph_Graph6->GetXaxis()->SetLabelOffset(0.01);
   Graph_Graph6->GetXaxis()->SetLabelSize(0.045);
   Graph_Graph6->GetXaxis()->SetTitleSize(0.055);
   Graph_Graph6->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph6->GetXaxis()->SetTitleFont(42);
   Graph_Graph6->GetYaxis()->SetLabelFont(42);
   Graph_Graph6->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph6->GetYaxis()->SetLabelSize(0.045);
   Graph_Graph6->GetYaxis()->SetTitleSize(0.055);
   Graph_Graph6->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph6->GetYaxis()->SetTitleFont(42);
   Graph_Graph6->GetZaxis()->SetLabelFont(42);
   Graph_Graph6->GetZaxis()->SetLabelSize(0.045);
   Graph_Graph6->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph6);
   
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
   grae->SetPoint(0,0.25,0.999177);
   grae->SetPointError(0,0.25,0.25,0,0);
   grae->SetPoint(1,0.75,1.02266);
   grae->SetPointError(1,0.25,0.25,0,0);
   grae->SetPoint(2,1.25,1.1311);
   grae->SetPointError(2,0.25,0.25,0,0);
   grae->SetPoint(3,1.75,1.44437);
   grae->SetPointError(3,0.25,0.25,0,0);
   grae->SetPoint(4,2.2,2.17478);
   grae->SetPointError(4,0.2,0.2,0,0);
   grae->Draw("z");
   
   TGraphErrors *gre = new TGraphErrors(5);
   gre->SetName("gA1p_exp_1");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetFillStyle(0);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,0.25,1.058523);
   gre->SetPointError(0,0,0.04102057);
   gre->SetPoint(1,0.75,1.03161);
   gre->SetPointError(1,0,0.03482058);
   gre->SetPoint(2,1.25,1.163431);
   gre->SetPointError(2,0,0.05161941);
   gre->SetPoint(3,1.75,1.629723);
   gre->SetPointError(3,0,0.0885695);
   gre->SetPoint(4,2.2,2.471921);
   gre->SetPointError(4,0,0.1901017);
   gre->Draw("||");
   
   gre = new TGraphErrors(5);
   gre->SetName("gA1p_exp_statonly_1");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,0.25,1.058523);
   gre->SetPointError(0,0,0.02967911);
   gre->SetPoint(1,0.75,1.03161);
   gre->SetPointError(1,0,0.0298508);
   gre->SetPoint(2,1.25,1.163431);
   gre->SetPointError(2,0,0.0372495);
   gre->SetPoint(3,1.75,1.629723);
   gre->SetPointError(3,0,0.05830586);
   gre->SetPoint(4,2.2,2.471921);
   gre->SetPointError(4,0,0.1112945);
   
   TH1F *Graph_gA1p_exp_statonly_13 = new TH1F("Graph_gA1p_exp_statonly_13","Graph",100,0.055,2.395);
   Graph_gA1p_exp_statonly_13->SetMinimum(0.8436134);
   Graph_gA1p_exp_statonly_13->SetMaximum(2.741361);
   Graph_gA1p_exp_statonly_13->SetDirectory(0);
   Graph_gA1p_exp_statonly_13->SetStats(0);
   Graph_gA1p_exp_statonly_13->SetLineStyle(0);
   Graph_gA1p_exp_statonly_13->SetMarkerStyle(20);
   Graph_gA1p_exp_statonly_13->GetXaxis()->SetLabelFont(42);
   Graph_gA1p_exp_statonly_13->GetXaxis()->SetLabelOffset(0.01);
   Graph_gA1p_exp_statonly_13->GetXaxis()->SetLabelSize(0.045);
   Graph_gA1p_exp_statonly_13->GetXaxis()->SetTitleSize(0.055);
   Graph_gA1p_exp_statonly_13->GetXaxis()->SetTitleOffset(1.1);
   Graph_gA1p_exp_statonly_13->GetXaxis()->SetTitleFont(42);
   Graph_gA1p_exp_statonly_13->GetYaxis()->SetLabelFont(42);
   Graph_gA1p_exp_statonly_13->GetYaxis()->SetLabelOffset(0.01);
   Graph_gA1p_exp_statonly_13->GetYaxis()->SetLabelSize(0.045);
   Graph_gA1p_exp_statonly_13->GetYaxis()->SetTitleSize(0.055);
   Graph_gA1p_exp_statonly_13->GetYaxis()->SetTitleOffset(1.4);
   Graph_gA1p_exp_statonly_13->GetYaxis()->SetTitleFont(42);
   Graph_gA1p_exp_statonly_13->GetZaxis()->SetLabelFont(42);
   Graph_gA1p_exp_statonly_13->GetZaxis()->SetLabelSize(0.045);
   Graph_gA1p_exp_statonly_13->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_gA1p_exp_statonly_13);
   
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
   TLegendEntry *entry=leg->AddEntry("gA1p_exp_1","Data","p");
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
   
   pt = new TPaveText(0.38,1.75,0.88,1.95,"br");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetTextAlign(13);
   pt->SetTextFont(42);
   pt->SetTextSize(0.05);
   text = pt->AddText("W^{+} #rightarrow #font[12]{l}^{+} + #nu");
   pt->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
