#ifndef SIGMOIDGAUS_PDF_HH
#define SIGMOIDGAUS_PDF_HH

#include "GooPdf.hh" 

class SigmoidGausPdf : public GooPdf {
public:
  SigmoidGausPdf (std::string n, Variable* _x,Variable* p0, Variable* p1, Variable* p2, Variable* p3, Variable* p4,Variable* mean,Variable* sigma) ; 

private:

};

#endif
