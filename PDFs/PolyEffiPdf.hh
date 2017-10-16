#ifndef POLYEFFI_PDF_HH
#define POLYEFFI_PDF_HH

#include "GooPdf.hh" 

class PolyEffiPdf : public GooPdf {
public:
  PolyEffiPdf (std::string n, Variable* _x, std::vector<Variable*> weights, Variable* x0 = 0, unsigned int lowestDegree = 0); 
  PolyEffiPdf (string n, vector<Variable*> obses, vector<Variable*> coeffs, vector<Variable*> offsets, unsigned int maxDegree); 
   __host__ fptype integrate (fptype lo, fptype hi) const; 
  //__host__ virtual bool hasAnalyticIntegral () const {return (1 == observables.size());} 
  __host__ fptype getCoefficient (int coef) const;
  __host__ virtual bool hasAnalyticIntegral () const {return true;}

private:
  Variable* center; 
};

#endif
