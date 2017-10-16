#ifndef SIGMOID_PDF_HH
#define SIGMOID_PDF_HH

#include "GooPdf.hh" 

class SigmoidPdf : public GooPdf {
public:
  SigmoidPdf (std::string n, Variable* _x,Variable* p0, Variable* p1, Variable* p2, Variable* p3) ; 

private:

};

#endif
