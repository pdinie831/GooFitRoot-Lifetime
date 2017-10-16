#ifndef EXPGausProdEffiB_PDF_HH
#define EXPGausProdEffiB_PDF_HH

#include "GooPdf.hh" 

class ExpGausProdEffiBPdf : public GooPdf {
public:
  ExpGausProdEffiBPdf (std::string n, Variable* _x, Variable* _s, Variable* m,  Variable* t, Variable* ss, Variable* ms,  Variable* ts, 
                   Variable* lo, Variable* hi, Variable* los, Variable* his); 
   __host__ fptype integrate (fptype lo, fptype hi) const; 
   __host__ virtual bool hasAnalyticIntegral () const {return true;} 

private:


};

#endif
