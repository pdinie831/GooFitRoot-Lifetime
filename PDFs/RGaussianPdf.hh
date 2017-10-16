#ifndef RGAUSSIAN_PDF_HH
#define RGAUSSIAN_PDF_HH

#include "GooPdf.hh" 

class RGaussianPdf : public GooPdf {
public:
  RGaussianPdf (std::string n, Variable* _x, Variable* m, Variable* s); 
  __host__ fptype integrate (fptype lo, fptype hi) const; 
  __host__ virtual bool hasAnalyticIntegral () const {return true;} 



private:

};

#endif
