{
//=========Macro generated from canvas: c1/c1
//=========  (Tue Dec  2 14:50:33 2014) by ROOT version5.34/10
   TCanvas *c1 = new TCanvas("c1", "c1",0,23,600,600);
   gStyle->SetOptFit(1);
   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   c1->Range(-0.4746836,0.4716049,2.689873,1.459259);
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
   grae->SetPoint(0,0.25,0.99022);
   grae->SetPointError(0,0.25,0.25,0.00855494,0.01078848);
   grae->SetPoint(1,0.75,0.968122);
   grae->SetPointError(1,0.25,0.25,0.02487038,0.02671408);
   grae->SetPoint(2,1.25,0.942255);
   grae->SetPointError(2,0.25,0.25,0.03707417,0.03485948);
   grae->SetPoint(3,1.75,0.919187);
   grae->SetPointError(3,0.25,0.25,0.02794373,0.03604308);
   grae->SetPoint(4,2.2,0.924656);
   grae->SetPointError(4,0.2,0.2,0.02711127,0.03571177);
   
   TH1F *Graph_Graph7 = new TH1F("Graph_Graph7","Graph",100,0,2.5);
   Graph_Graph7->SetMinimum(0.6);
   Graph_Graph7->SetMaximum(1.4);
   Graph_Graph7->SetDirectory(0);
   Graph_Graph7->SetStats(0);
   Graph_Graph7->SetLineStyle(0);
   Graph_Graph7->SetMarkerStyle(20);
   Graph_Graph7->GetXaxis()->SetTitle("#eta_{lab}");
   Graph_Graph7->GetXaxis()->SetNdivisions(505);
   Graph_Graph7->GetXaxis()->SetLabelFont(42);
   Graph_Graph7->GetXaxis()->SetLabelOffset(0.01);
   Graph_Graph7->GetXaxis()->SetLabelSize(0.05);
   Graph_Graph7->GetXaxis()->SetTitleSize(0.05);
   Graph_Graph7->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph7->GetXaxis()->SetTitleFont(42);
   Graph_Graph7->GetYaxis()->SetTitle("N^{-}(+#eta_{lab})/N^{-}(-#eta_{lab})");
   Graph_Graph7->GetYaxis()->SetLabelFont(42);
   Graph_Graph7->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph7->GetYaxis()->SetLabelSize(0.05);
   Graph_Graph7->GetYaxis()->SetTitleSize(0.05);
   Graph_Graph7->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph7->GetYaxis()->SetTitleFont(42);
   Graph_Graph7->GetZaxis()->SetLabelFont(42);
   Graph_Graph7->GetZaxis()->SetLabelSize(0.045);
   Graph_Graph7->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph7);
   
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
   grae->SetPoint(0,0.25,0.99022);
   grae->SetPointError(0,0.25,0.25,0,0);
   grae->SetPoint(1,0.75,0.968122);
   grae->SetPointError(1,0.25,0.25,0,0);
   grae->SetPoint(2,1.25,0.942255);
   grae->SetPointError(2,0.25,0.25,0,0);
   grae->SetPoint(3,1.75,0.919187);
   grae->SetPointError(3,0.25,0.25,0,0);
   grae->SetPoint(4,2.2,0.924656);
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
   grae->SetPoint(0,0.25,0.966757);
   grae->SetPointError(0,0.25,0.25,0.01079338,0.01549479);
   grae->SetPoint(1,0.75,0.905839);
   grae->SetPointError(1,0.25,0.25,0.03091994,0.0344406);
   grae->SetPoint(2,1.25,0.853546);
   grae->SetPointError(2,0.25,0.25,0.04381241,0.04491883);
   grae->SetPoint(3,1.75,0.816165);
   grae->SetPointError(3,0.25,0.25,0.04722097,0.05274991);
   grae->SetPoint(4,2.2,0.818614);
   grae->SetPointError(4,0.2,0.2,0.05398694,0.05718776);
   
   TH1F *Graph_Graph8 = new TH1F("Graph_Graph8","Graph",100,0,2.64);
   Graph_Graph8->SetMinimum(0.7428646);
   Graph_Graph8->SetMaximum(1.004014);
   Graph_Graph8->SetDirectory(0);
   Graph_Graph8->SetStats(0);
   Graph_Graph8->SetLineStyle(0);
   Graph_Graph8->SetMarkerStyle(20);
   Graph_Graph8->GetXaxis()->SetLabelFont(42);
   Graph_Graph8->GetXaxis()->SetLabelOffset(0.01);
   Graph_Graph8->GetXaxis()->SetLabelSize(0.045);
   Graph_Graph8->GetXaxis()->SetTitleSize(0.055);
   Graph_Graph8->GetXaxis()->SetTitleOffset(1.1);
   Graph_Graph8->GetXaxis()->SetTitleFont(42);
   Graph_Graph8->GetYaxis()->SetLabelFont(42);
   Graph_Graph8->GetYaxis()->SetLabelOffset(0.01);
   Graph_Graph8->GetYaxis()->SetLabelSize(0.045);
   Graph_Graph8->GetYaxis()->SetTitleSize(0.055);
   Graph_Graph8->GetYaxis()->SetTitleOffset(1.4);
   Graph_Graph8->GetYaxis()->SetTitleFont(42);
   Graph_Graph8->GetZaxis()->SetLabelFont(42);
   Graph_Graph8->GetZaxis()->SetLabelSize(0.045);
   Graph_Graph8->GetZaxis()->SetTitleFont(42);
   grae->SetHistogram(Graph_Graph8);
   
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
   grae->SetPoint(0,0.25,0.966757);
   grae->SetPointError(0,0.25,0.25,0,0);
   grae->SetPoint(1,0.75,0.905839);
   grae->SetPointError(1,0.25,0.25,0,0);
   grae->SetPoint(2,1.25,0.853546);
   grae->SetPointError(2,0.25,0.25,0,0);
   grae->SetPoint(3,1.75,0.816165);
   grae->SetPointError(3,0.25,0.25,0,0);
   grae->SetPoint(4,2.2,0.818614);
   grae->SetPointError(4,0.2,0.2,0,0);
   grae->Draw("z");
   
   TGraphErrors *gre = new TGraphErrors(5);
   gre->SetName("gA1m_exp_1");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetFillStyle(0);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,0.25,0.9653543);
   gre->SetPointError(0,0,0.03883504);
   gre->SetPoint(1,0.75,0.9455861);
   gre->SetPointError(1,0,0.03890125);
   gre->SetPoint(2,1.25,0.8116163);
   gre->SetPointError(2,0,0.03897397);
   gre->SetPoint(3,1.75,0.7434323);
   gre->SetPointError(3,0,0.04405454);
   gre->SetPoint(4,2.2,0.8833141);
   gre->SetPointError(4,0,0.07076644);
   gre->Draw("||");
   
   gre = new TGraphErrors(5);
   gre->SetName("gA1m_exp_statonly_1");
   gre->SetTitle("Graph");
   gre->SetFillColor(1);
   gre->SetLineWidth(2);
   gre->SetMarkerStyle(20);
   gre->SetMarkerSize(1.2);
   gre->SetPoint(0,0.25,0.9653543);
   gre->SetPointError(0,0,0.03137784);
   gre->SetPoint(1,0.75,0.9455861);
   gre->SetPointError(1,0,0.03044051);
   gre->SetPoint(2,1.25,0.8116163);
   gre->SetPointError(2,0,0.02861161);
   gre->SetPoint(3,1.75,0.7434323);
   gre->SetPointError(3,0,0.02860939);
   gre->SetPoint(4,2.2,0.8833141);
   gre->SetPointError(4,0,0.03901511);
   
   TH1F *Graph_gA1m_exp_statonly_14 = new TH1F("Graph_gA1m_exp_statonly_14","Graph",100,0.055,2.395);
   Graph_gA1m_exp_statonly_14->SetMinimum(0.686632);
   Graph_gA1m_exp_statonly_14->SetMaximum(1.024923);
   Graph_gA1m_exp_statonly_14->SetDirectory(0);
   Graph_gA1m_exp_statonly_14->SetStats(0);
   Graph_gA1m_exp_statonly_14->SetLineStyle(0);
   Graph_gA1m_exp_statonly_14->SetMarkerStyle(20);
   Graph_gA1m_exp_statonly_14->GetXaxis()->SetLabelFont(42);
   Graph_gA1m_exp_statonly_14->GetXaxis()->SetLabelOffset(0.01);
   Graph_gA1m_exp_statonly_14->GetXaxis()->SetLabelSize(0.045);
   Graph_gA1m_exp_statonly_14->GetXaxis()->SetTitleSize(0.055);
   Graph_gA1m_exp_statonly_14->GetXaxis()->SetTitleOffset(1.1);
   Graph_gA1m_exp_statonly_14->GetXaxis()->SetTitleFont(42);
   Graph_gA1m_exp_statonly_14->GetYaxis()->SetLabelFont(42);
   Graph_gA1m_exp_statonly_14->GetYaxis()->SetLabelOffset(0.01);
   Graph_gA1m_exp_statonly_14->GetYaxis()->SetLabelSize(0.045);
   Graph_gA1m_exp_statonly_14->GetYaxis()->SetTitleSize(0.055);
   Graph_gA1m_exp_statonly_14->GetYaxis()->SetTitleOffset(1.4);
   Graph_gA1m_exp_statonly_14->GetYaxis()->SetTitleFont(42);
   Graph_gA1m_exp_statonly_14->GetZaxis()->SetLabelFont(42);
   Graph_gA1m_exp_statonly_14->GetZaxis()->SetLabelSize(0.045);
   Graph_gA1m_exp_statonly_14->GetZaxis()->SetTitleFont(42);
   gre->SetHistogram(Graph_gA1m_exp_statonly_14);
   
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
   
   TLegend *leg = new TLegend(0.5,0.62,0.8,0.77,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.04);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);
   TLegendEntry *entry=leg->AddEntry("gA1m_exp_1","Data","p");
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
   
   pt = new TPaveText(0.25,1.2,0.75,1.3,"br");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetTextAlign(13);
   pt->SetTextFont(42);
   pt->SetTextSize(0.05);
   text = pt->AddText("W^{-} #rightarrow #font[12]{l}^{-} + #nu");
   pt->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
