#ifndef ExpGausProdBResoCorr_PDF_HH
#define ExpGausProdBResoCorr_PDF_HH

#include "GooPdf.hh" 

class ExpGausProdBResoCorrPdf : public GooPdf {
public:
  ExpGausProdBResoCorrPdf (std::string n, Variable* _x, Variable* _s, Variable* m,  Variable* t, Variable* ss, Variable* ms,  Variable* ts, 
                   Variable* lo, Variable* hi, Variable* los, Variable* his, Variable* corrP1, Variable* corrP0); 
   __host__ fptype integrate (fptype lo, fptype hi) const; 
   __host__ virtual bool hasAnalyticIntegral () const {return true;} 

private:


};

#endif
