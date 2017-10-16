#ifndef EXPGAUSPEESigmaB_PDF_HH
#define EXPGAUSPEESigmaB_PDF_HH

#include "GooPdf.hh" 

class ExpGausPEESigmaBPdf : public GooPdf {
public:
  ExpGausPEESigmaBPdf (std::string n, Variable* _x, Variable* _s, Variable* m,  Variable* t, Variable* lo,  Variable* hi); 
   __host__ fptype integrate (fptype lo, fptype hi) const; 
   __host__ virtual bool hasAnalyticIntegral () const {return true;} 

private:


};

#endif
