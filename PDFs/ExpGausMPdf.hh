#ifndef EXPGAUSM_PDF_HH
#define EXPGAUSM_PDF_HH

#include "GooPdf.hh" 

class ExpGausMPdf : public GooPdf {
public:
  ExpGausMPdf (std::string n, Variable* _x, Variable* s, Variable* m, Variable* t); 

private:

};

#endif
