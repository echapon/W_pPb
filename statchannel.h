#ifndef statchannel_h
#define statchannel_h

#include <vector>
#include <fstream>
#include <iostream>
#include <string>

#include "TGraphErrors.h"

#include "lumierror.h"

using namespace std;

enum asym {
   YIELDS=0,
   YIELDS_PM,
   CH,
   A1p,
   A1m,
   A2,
   A3,
   A4
};

class statchannel
{
   private:
      string m_label;
      int m_nbins;
      vector<double> m_bins;
      vector<double> m_yields;
      vector<double> m_staterr;
      vector< vector<double> > m_eff;
      vector< vector<double> > m_efferr;
      vector<double> m_th;
      vector<double> m_therr;

   public:
      statchannel() {};
      ~statchannel() {};

      void read(int n, const char* filename);
      void print();
      
      // getters
      string get_label() const {return m_label;};
      int get_nbins() const {return m_nbins;};
      vector<double> get_bins() const {return m_bins;};
      vector<double> get_yields() const {return m_yields;};
      vector<double> get_staterr() const {return m_staterr;};
      vector< vector<double> > get_eff() const {return m_eff;};
      vector< vector<double> > get_efferr() const {return m_efferr;};
      vector<double> get_th() const {return m_th;};
      vector<double> get_therr() const {return m_therr;};

      // static
      static TGraphErrors* graph(const statchannel chanp, const statchannel chanm, asym mode=YIELDS, bool isth=true, bool dosyst=false);
      static double total(int n, vector<double> v);
      static double total2(int n, vector<double> v);
};

#endif // #ifdef statchannel_h
