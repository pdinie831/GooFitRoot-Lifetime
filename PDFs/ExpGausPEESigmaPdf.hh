#ifndef EXPGAUSPEESigma_PDF_HH
#define EXPGAUSPEESigma_PDF_HH

#include "GooPdf.hh" 

class ExpGausPEESigmaPdf : public GooPdf {
public:
  ExpGausPEESigmaPdf (std::string n, Variable* _x, Variable* _s, Variable* m,  Variable* t); 
   __host__ fptype integrate (fptype lo, fptype hi) const; 
   __host__ virtual bool hasAnalyticIntegral () const {return true;} 

private:


};

#endif
