#ifndef EXPGAUSProd_PDF_HH
#define EXPGAUSProd_PDF_HH

#include "GooPdf.hh" 

class ExpGausProdPdf : public GooPdf {
public:
  ExpGausProdPdf (std::string n, Variable* _x, Variable* _s, Variable* m,  Variable* t, Variable* ss, Variable* ms,  Variable* ts); 
   __host__ fptype integrate (fptype lo, fptype hi) const; 
   __host__ virtual bool hasAnalyticIntegral () const {return true;} 

private:


};

#endif
