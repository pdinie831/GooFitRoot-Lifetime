#ifndef EXPGAUSWITHINT_PDF_HH
#define EXPGAUSWITHINT_PDF_HH

#include "GooPdf.hh" 

class ExpGausWithIntPdf : public GooPdf {
public:
  ExpGausWithIntPdf (std::string n, Variable* _x, Variable* m, Variable* s, Variable* t); 
   __host__ fptype integrate (fptype lo, fptype hi) const; 
   __host__ virtual bool hasAnalyticIntegral () const {return true;} 

private:

};

#endif
