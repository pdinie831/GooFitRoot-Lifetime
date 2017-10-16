#ifndef GAMMA_PDF_HH
#define GAMMA_PDF_HH

#include "GooPdf.hh" 

class GammaPdf : public GooPdf {
public:
  GammaPdf (std::string n, Variable* _x, Variable* gamma, Variable* beta, Variable* mu); 
//  __host__ fptype integrate (fptype lo, fptype hi) const; 
//  __host__ virtual bool hasAnalyticIntegral () const {return true;} 



private:

};

#endif
