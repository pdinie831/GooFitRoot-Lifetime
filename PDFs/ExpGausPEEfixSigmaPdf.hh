#ifndef ExpGausPEEfixSigma_PDF_HH
#define ExpGausPEEfixSigma_PDF_HH

#include "GooPdf.hh" 

class ExpGausPEEfixSigmaPdf : public GooPdf {
public:
  ExpGausPEEfixSigmaPdf (std::string n, Variable* _x, Variable* _s, Variable* m,  Variable* t); 
   __host__ fptype integrate (fptype lo, fptype hi) const; 
   __host__ virtual bool hasAnalyticIntegral () const {return true;} 

private:


};

#endif
