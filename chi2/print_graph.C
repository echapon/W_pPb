void print_graph(TGraphErrors *g1, TGraphErrors *g2, double scale=1., bool do_deta=true)
{
   double deta=1.;
   for (int i=0; i<g1->GetN(); i++)
   {
      if (do_deta) deta = 2*g1->GetEX()[i];
      cout << "[" << g1->GetX()[i]-g1->GetEX()[i] << "," << g1->GetX()[i]+g1->GetEX()[i] <<"]: ";
      cout << g1->GetY()[i]*scale*deta << "+-" << g1->GetEY()[i]*scale*deta << "+-" << g2->GetEY()[i]*scale*deta << endl;
   }
}
