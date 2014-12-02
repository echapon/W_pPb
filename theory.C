#include "Math/Interpolator.h"

double th_plus_cteq(double val)
{
   vector<double> x,y_plus_cteq;
   x.push_back(-3);
   x.push_back(-2.5);
   x.push_back(-2);
   x.push_back(-1.5);
   x.push_back(-1);
   x.push_back(-.5);
   x.push_back(0);
   x.push_back(.5);
   x.push_back(1);
   x.push_back(1.5);
   x.push_back(2);
   x.push_back(2.5);
   x.push_back(3);

   y_plus_cteq.push_back(126868.4);
   y_plus_cteq.push_back(254823.6);
   y_plus_cteq.push_back(361235.3);
   y_plus_cteq.push_back(418815.4);
   y_plus_cteq.push_back(443500.2);
   y_plus_cteq.push_back(456005.7);
   y_plus_cteq.push_back(468217.6);
   y_plus_cteq.push_back(484117.0);
   y_plus_cteq.push_back(502132.4);
   y_plus_cteq.push_back(509063.6);
   y_plus_cteq.push_back(473572.4);
   y_plus_cteq.push_back(360251.7);
   y_plus_cteq.push_back(188964.1);

   ROOT::Math::Interpolator inter_plus_cteq(x,y_plus_cteq, ROOT::Math::Interpolation::kCSPLINE);

   return inter_plus_cteq.Eval(val);
}

double th_minus_cteq(double val)
{
   vector<double> x,y_minus_cteq;
   x.push_back(-3);
   x.push_back(-2.5);
   x.push_back(-2);
   x.push_back(-1.5);
   x.push_back(-1);
   x.push_back(-.5);
   x.push_back(0);
   x.push_back(.5);
   x.push_back(1);
   x.push_back(1.5);
   x.push_back(2);
   x.push_back(2.5);
   x.push_back(3);

   y_minus_cteq.push_back(271489.0);
   y_minus_cteq.push_back(322681.6);
   y_minus_cteq.push_back(343223.2);
   y_minus_cteq.push_back(352120.6);
   y_minus_cteq.push_back(356207.5);
   y_minus_cteq.push_back(356162.1);
   y_minus_cteq.push_back(351002.9);
   y_minus_cteq.push_back(339307.4);
   y_minus_cteq.push_back(320784.2);
   y_minus_cteq.push_back(296742.8);
   y_minus_cteq.push_back(267804.6);
   y_minus_cteq.push_back(228431.7);
   y_minus_cteq.push_back(169198.6);

   ROOT::Math::Interpolator inter_minus_cteq(x,y_minus_cteq, ROOT::Math::Interpolation::kCSPLINE);

   return inter_minus_cteq.Eval(val);
}

double th_plus_eps(double val)
{
   vector<double> x,y_plus_eps;
   x.push_back(-3);
   x.push_back(-2.5);
   x.push_back(-2);
   x.push_back(-1.5);
   x.push_back(-1);
   x.push_back(-.5);
   x.push_back(0);
   x.push_back(.5);
   x.push_back(1);
   x.push_back(1.5);
   x.push_back(2);
   x.push_back(2.5);
   x.push_back(3);

   y_plus_eps.push_back(119925.9);
   y_plus_eps.push_back(249967.7);
   y_plus_eps.push_back(367369.4);
   y_plus_eps.push_back(433169.2);
   y_plus_eps.push_back(457427.7);
   y_plus_eps.push_back(460981.2);
   y_plus_eps.push_back(457752.1);
   y_plus_eps.push_back(457031.8);
   y_plus_eps.push_back(459302.0);
   y_plus_eps.push_back(454393.3);
   y_plus_eps.push_back(416368.0);
   y_plus_eps.push_back(312309.4);
   y_plus_eps.push_back(162876.1);

   ROOT::Math::Interpolator inter_plus_eps(x,y_plus_eps, ROOT::Math::Interpolation::kCSPLINE);

   return inter_plus_eps.Eval(val);
}

double th_minus_eps(double val)
{
   vector<double> x,y_minus_eps;
   x.push_back(-3);
   x.push_back(-2.5);
   x.push_back(-2);
   x.push_back(-1.5);
   x.push_back(-1);
   x.push_back(-.5);
   x.push_back(0);
   x.push_back(.5);
   x.push_back(1);
   x.push_back(1.5);
   x.push_back(2);
   x.push_back(2.5);
   x.push_back(3);

   y_minus_eps.push_back(271284.6);
   y_minus_eps.push_back(330370.5);
   y_minus_eps.push_back(356200.1);
   y_minus_eps.push_back(365296.1);
   y_minus_eps.push_back(365943.3);
   y_minus_eps.push_back(359035.5);
   y_minus_eps.push_back(344691.9);
   y_minus_eps.push_back(324350.0);
   y_minus_eps.push_back(298937.7);
   y_minus_eps.push_back(270459.5);
   y_minus_eps.push_back(239464.2);
   y_minus_eps.push_back(200540.7);
   y_minus_eps.push_back(146484.2);

   ROOT::Math::Interpolator inter_minus_eps(x,y_minus_eps, ROOT::Math::Interpolation::kCSPLINE);

   return inter_minus_eps.Eval(val);
}

void theory(double luminb = 34.621)
{
   double a = 208;
   double lumifb = luminb*1e-6*a;

   vector<double> myx;
   myx.push_back(-2.4);
   myx.push_back(-2.0);
   myx.push_back(-1.5);
   myx.push_back(-1.);
   myx.push_back(-.5);
   myx.push_back(-0);
   myx.push_back(.5);
   myx.push_back(1);
   myx.push_back(1.5);
   myx.push_back(2.);
   myx.push_back(2.4);

   TF1 *f_plus_eps = new TF1("f_plus_eps","th_plus_eps(x)",-3,3);
   TF1 *f_minus_eps = new TF1("f_minus_eps","th_minus_eps(x)",-3,3);
   TF1 *f_plus_cteq = new TF1("f_plus_cteq","th_plus_cteq(x)",-3,3);
   TF1 *f_minus_cteq = new TF1("f_minus_cteq","th_minus_cteq(x)",-3,3);

   cout << "EPS09, W+" << endl;
   for (int i=0; i<myx.size(); i++) cout << myx[i] << " ";
   cout << endl;
   for (int i=0; i<myx.size()-1; i++) cout << f_plus_eps->Integral(-myx[i+1]-0.465,-myx[i]-0.465)*lumifb << " ";
   cout << endl;
   cout << "EPS09, W-" << endl;
   for (int i=0; i<myx.size(); i++) cout << myx[i] << " ";
   cout << endl;
   for (int i=0; i<myx.size()-1; i++) cout << f_minus_eps->Integral(-myx[i+1]-0.465,-myx[i]-0.465)*lumifb << " ";
   cout << endl;
   cout << "CTEQ6.6, W+" << endl;
   for (int i=0; i<myx.size(); i++) cout << myx[i] << " ";
   cout << endl;
   for (int i=0; i<myx.size()-1; i++) cout << f_plus_cteq->Integral(-myx[i+1]-0.465,-myx[i]-0.465)*lumifb << " ";
   cout << endl;
   cout << "CTEQ6.6, W-" << endl;
   for (int i=0; i<myx.size(); i++) cout << myx[i] << " ";
   cout << endl;
   for (int i=0; i<myx.size()-1; i++) cout << f_minus_cteq->Integral(-myx[i+1]-0.465,-myx[i]-0.465)*lumifb << " ";
   cout << endl;
}
