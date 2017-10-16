#ifndef SIGMOIDBP_PDF_HH
#define SIGMOIDBP_PDF_HH

#include "GooPdf.hh" 

class SigmoidBpPdf : public GooPdf {
public:
  SigmoidBpPdf (std::string n, Variable* _x,Variable* p0, Variable* p1, Variable* p2, Variable* p3, Variable* p4, Variable* p5, Variable* p6, Variable* p7) ; 

private:

};

#endif
