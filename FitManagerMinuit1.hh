#ifndef FITMANAGER_MINUIT1_HH
#define FITMANAGER_MINUIT1_HH

#include "TMinuit.h" 
extern PdfBase* pdfPointer; 
extern int numPars; 

void FitFun(int &npar, double *gin, double &fun, double *fp, int iflag); 

class FitManager { 
public:
  FitManager (PdfBase* dat);
  ~FitManager ();
  void setMaxCalls (double mxc) {overrideCallLimit = mxc;}
  void setupMinuit ();
  void runMigrad (); 
  void runHesse (); 
  void runMinos (); 
  void runImprove (); 
  void runCommand (const char * command); 
  void fit (); 
  TMinuit* getMinuitObject () {return minuit;} 
  void getMinuitValues () const;
  TMinuit* minuit; 
private:
  double overrideCallLimit; 
};

#endif 
