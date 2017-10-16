//
//
//
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sys/stat.h>
#include "Riostream.h"
#include <map>
#include <string>
//#include <boost/algorithm/string.hpp>
//#include <boost/algorithm/string/trim.hpp>
#include <vector>
#include <math.h>
#include <TCint.h>
#include <TGenericClassInfo.h> 
#include <TH1.h>
#include <TH2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TSystem.h>
#include <TTree.h>
#include "TBranch.h"
#include <TApplication.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TStyle.h> 
#include <TF1.h>  
#include <TF2.h> 
#include <TLorentzVector.h>
#include <TVector3.h>
#include "TDSet.h"
#include "TChain.h"
#include <time.h> 
#include <TSystemDirectory.h>
#include <TLegend.h>
#include "FitManager.hh"
#include "UnbinnedDataSet.hh" 
#include "LandauPdf.hh" 
#include "NovosibirskPdf.hh"
#include "BifurGaussPdf.hh" 
#include "SimpleCheby2Pdf.hh" 
#include "GammaPdf.hh" 

#include "TRandom.h" 
#if HAVE_ROOT
#  include "Variable.hh"
#  include "TH1F.h"
//#  include "TH2F.h"			// unused?
#  include "TStyle.h"
#  include "TCanvas.h"
#else
#  include "fakeTH1F.h"
#endif

#include <sys/time.h>
#include <sys/times.h>
#include <iostream>
// GooFit stuff
#include "Variable.hh" 
#include "KinLimitBWPdf.hh" 
#include "ConvolutionPdf.hh"
#include "GaussianPdf.hh"
#include "ScaledGaussianPdf.hh"
#include "ArgusPdf.hh"
#include "AddPdf.hh"
#include "PolynomialPdf.hh" 
#include "FitManager.hh" 
#include "ExpGausPdf.hh" 
#include "ExpGausPEEPdf.hh" 
#include "ExpGausPEESigmaPdf.hh" 
#include "ExpPdf.hh" 
#include "ProdPdf.hh" 
#include "RGaussianPdf.hh"
#include "BifurGaussPdf.hh"
#include "ExpGausMPdf.hh" 
#include "ExpGausWithIntPdf.hh"
#include "ExpGausPEEfixSigmaPdf.hh" 
#include "ExpGausProdPdf.hh" 
#include "ExpGausProdBPdf.hh" 
// System stuff
#include <fstream> 
#include <sys/time.h>
#include <sys/times.h>

timeval startTime, stopTime, totalTime;
clock_t startCPU, stopCPU; 
tms startProc, stopProc; 

void fitTauSBModel();
using namespace std; 


int main (int argc, char** argv) {

   TApplication app("App",&argc, argv);

   fitTauSBModel(); 
//  app.Run() ;
   cout<<"esco..." <<endl;
   return 0 ;
}


void fitTauSBModel(){
//
//
   gettimeofday(&startTime, NULL);
   startCPU = times(&startProc);
//
   gROOT ->Reset();
   gROOT->SetStyle("Plain");
   gStyle->SetOptStat(000000);
   gStyle->SetOptFit(000000);
//
   float PlotLineWidth = 1.2;
//

   TCanvas* c1 = new TCanvas("c1","Mass PLOTS",200,10,900,780);
   TCanvas* c2 = new TCanvas("c2","cTau PLOTS",200,10,900,780);
   TCanvas* c3 = new TCanvas("c3","STau PLOTS",200,10,900,780);
//   TCanvas*c2 = new TCanvas("c2","PLOTS",200,10,900,780);
   //TPad* pad1 = (TPad*)c1->GetPad(0);
//   TPad* pad2 = (TPad*)c2->GetPad(0);
   //pad1->SetLeftMargin(0.15); 
//   pad2->SetLeftMargin(0.15); 





  Char_t    InputFileName[300] = "testproof2DCut-Run2016-SaraCuts-NewCutPDGMass.root";
  Char_t    InputTauBcTreeName[10]   = "TauBcTree";
  TFile*InputFile = TFile::Open(InputFileName,"READ","ROOT file");
  
   Char_t    OutFileName[300] =
   "testGoofit2DBc-2016.root";
   gSystem->Exec(Form("rm %s",OutFileName));
   TFile*OutFile = TFile::Open(OutFileName,"RECREATE");
  double xBcMass;
//  double xBcTau;
  double xBccTau;
  double xSBccTau;
//  double c_const       = 0.0299792458;
//     double XMinSign = 5.12;
//     double XMaxSign = 5.44;
//double XMinSign = 5.15;
//double XMaxSign = 5.4;
  double BcMass   = 6.273;
  double BcSigma  = 0.028;

  double NSigmaSB = 7.7;

  double XMinSign = 6.049000;
  double XMaxSign = 6.497000;
//   double XMinSBL = 4.879000;
//   double XMaxSBL = 5.159000;
//   double XMinSBR = 5.399000;
//   double XMaxSBR = 5.679000;
  double XMinSBL = BcMass -(6+NSigmaSB)*BcSigma;
  double XMaxSBL = BcMass - 6          *BcSigma;
  double XMinSBR = BcMass + 6          *BcSigma;
  double XMaxSBR = BcMass +(6+NSigmaSB)*BcSigma;
  double XMin = 0.01;
  double XMax = 0.1;
  double SXMin = 0.0003;
  double SXMax = 0.007;
  double XStepSign = 0.01;
  double XStepcTau = 0.001;
  double XStepScTau = 0.0001;
//  double XStepMinuit = 0.00001;
  double XHScale = 10;
//  
//  int    NFiness     = 900;
//  int    NFinessMass = 900;
//  
  double c_const       = 0.0299792458;
//

  Variable* xMass  = new Variable("xMass",XMinSign, XMaxSign); 
  xMass->numbins = (XMaxSign -XMinSign)/XStepSign;
  TH1F HxMass( "HxMass" , "Bc+ Mass"    ,  xMass->numbins, xMass->lowerlimit, xMass->upperlimit);
  TH1F pdfHist("pdfHist", "Bc+ Mass Fit",  XHScale*xMass->numbins, xMass->lowerlimit, xMass->upperlimit);
  TH1F sigHist("sigHist", "Bc+ Mass Fit",  XHScale*xMass->numbins, xMass->lowerlimit, xMass->upperlimit);
  TH1F bkgHist("bgkHist", "Bc+ Mass Fit",  XHScale*xMass->numbins, xMass->lowerlimit, xMass->upperlimit);
//
  Variable* xcTau  = new Variable("xcTau",XMin, XMax); 
  xcTau->numbins = (XMax -XMin)/XStepcTau;
  Variable* xScTau  = new Variable("xScTau",SXMin, SXMax); 
  xScTau->numbins = (SXMax -SXMin)/(XStepScTau);
// 5 xScTau->numbins = (XMax -XMin)/XStepcTau;

  std::cout<<"xMass ->numbins = "<<xMass ->numbins<<std::endl;
  std::cout<<"xcTau ->numbins = "<<xcTau ->numbins<<std::endl;
  std::cout<<"xScTau->numbins = "<<xScTau->numbins<<std::endl;
  TH1F HxcTau(    "HxcTau"   , "Bc+ cTau",          xcTau  ->numbins, xcTau ->lowerlimit,  xcTau ->upperlimit);
  TH1F HxcTauSB(  "HxcTauSB" , "Bc+ cTau SB",       xcTau  ->numbins, xcTau ->lowerlimit,  xcTau ->upperlimit);
  TH1F HxScTau(   "HxScTau"  , "Bc+ cTau Sigma",    xScTau ->numbins, xScTau->lowerlimit, xScTau ->upperlimit);
  TH1F HxScTauSB( "HxScTauSB", "Bc+ cTau Sigma SB", xScTau ->numbins, xScTau->lowerlimit, xScTau ->upperlimit);
//   TH1F pdf_cTau_Hist( "pdf_cTau_Hist" , "Bc+ cTau model   pdf",    xcTau ->numbins, xcTau ->lowerlimit, xcTau ->upperlimit);
//   TH1F sig_cTau_Hist( "sig_cTau_Hist" , "Bc+ cTau signal  pdf",    xcTau ->numbins, xcTau ->lowerlimit, xcTau ->upperlimit);
//   TH1F bkg_cTau_Hist( "bkg_cTau_Hist" , "Bc+ cTau bckg    pdf",    xcTau ->numbins, xcTau ->lowerlimit, xcTau ->upperlimit);
//   TH1F pdf_STau_Hist( "pdf_STau_Hist" , "Bc+ ScTau model  pdf",    xcTau ->numbins, xScTau ->lowerlimit, xScTau ->upperlimit);
//   TH1F sig_STau_Hist( "sig_STau_Hist" , "Bc+ ScTau signal pdf",    xcTau ->numbins, xScTau ->lowerlimit, xScTau ->upperlimit);
//   TH1F bkg_STau_Hist( "bkg_STau_Hist" , "Bc+ ScTau bckg   pdf",    xcTau ->numbins, xScTau ->lowerlimit, xScTau ->upperlimit);
// 
// //                                                                        xcTau ->numbins, xScTau ->lowerlimit, xScTau ->upperlimit);              
//   TH1F pdf_cTau_Hist2D( "pdf_cTau_Hist2D" , "Bc+ cTau model   pdf",    xcTau ->numbins, xcTau ->lowerlimit, xcTau ->upperlimit);
//   TH1F sig_cTau_Hist2D( "sig_cTau_Hist2D" , "Bc+ cTau signal  pdf",    xcTau ->numbins, xcTau ->lowerlimit, xcTau ->upperlimit);
//   TH1F bkg_cTau_Hist2D( "bkg_cTau_Hist2D" , "Bc+ cTau bckg    pdf",    xcTau ->numbins, xcTau ->lowerlimit, xcTau ->upperlimit);              

  TH2F pdf_cTauSTau_Hist2D( "pdf_cTauSTau_Hist2D" , "Bc+ cTau model   pdf",  XHScale*xcTau ->numbins, xcTau ->lowerlimit, xcTau ->upperlimit, XHScale*xScTau->numbins, xScTau ->lowerlimit, xScTau ->upperlimit);
  TH2F sig_cTauSTau_Hist2D( "sig_cTauSTau_Hist2D" , "Bc+ cTau model   sig",  XHScale*xcTau ->numbins, xcTau ->lowerlimit, xcTau ->upperlimit, XHScale*xScTau->numbins, xScTau ->lowerlimit, xScTau ->upperlimit);
  TH2F bkg_cTauSTau_Hist2D( "bkg_cTauSTau_Hist2D" , "Bc+ cTau model   bkg",  XHScale*xcTau ->numbins, xcTau ->lowerlimit, xcTau ->upperlimit, XHScale*xScTau->numbins, xScTau ->lowerlimit, xScTau ->upperlimit);
//  xMass->numbins = 400;
//  xcTau->numbins = 400;
//  xScTau->numbins = 400;
//
// Mass Spectrum
//
  Variable* mean   = new Variable("mean" ,6.2738, 6.2, 6.4);
  Variable* sigma  = new Variable("sigma",0.029 , 0.01, 1.);
//  Variable* mean   = new Variable("mean"  ,6.273,XStepMinuit, 6., 6.5);
//  Variable* sigma1 = new Variable("sigma1",0.028,XStepMinuit, 0., 1.);
//  Variable* sigma2 = new Variable("sigma2",0.030,XStepMinuit, 0., 1.);
//  Variable* sigma3 = new Variable("sigma3",0.060,XStepMinuit, 0., 1.);

  GaussianPdf* signalMass = new GaussianPdf("signalMass", xMass, mean, sigma);
//  GaussianPdf* gauss1 = new GaussianPdf("gauss1", xMass, mean, sigma1);
//  GaussianPdf* gauss2 = new GaussianPdf("gauss2", xMass, mean, sigma2);
//  GaussianPdf* gauss3 = new GaussianPdf("gauss3", xMass, mean, sigma3);

/*   Variable* meanBckgBc   = new Variable("meanBckgBc" ,5.360,0.00001, 5., 5.5);
  Variable* sigmaBckgBc  = new Variable("sigmaBckgBc",0.030,0.00001, 0. , 1. );
  Variable* meanBckgB0   = new Variable("meanBckgB0" ,5.090,0.00001, 5. , 5.2);
  Variable* sigmaBckgB0  = new Variable("sigmaBckgB0",0.025,0.00001, 0. , 1.);
 */  
//  Variable* meanBckgBc   = new Variable("meanBckgBc" ,5.37 , 5.3, 5.5);
//  Variable* sigmaBckgBc  = new Variable("sigmaBckgBc",0.033, 0.01 , 5. );
//  Variable* meanBckgB0   = new Variable("meanBckgB0" ,5.090, 5.0 , 5.25);
//  Variable* sigmaBckgB0  = new Variable("sigmaBckgB0",0.025, 0. , 3.);
//  Variable* wb1 = new Variable("wb1",1.75628e-01, 0., 1.);
//  Variable* wb2 = new Variable("wb2",1.17246e-01, 0., 1.);
//  Variable* wb3 = new Variable("wb3",0.1, 0., 1.);

//  RGaussianPdf* gaussBckgBc = new RGaussianPdf("gaussBckgBc", xMass, meanBckgBc, sigmaBckgBc);
//  RGaussianPdf* gaussBckgB0 = new RGaussianPdf("gaussBckgB0", xMass, meanBckgB0, sigmaBckgB0);
  

//  Variable* wg1 = new Variable("wg1",0.44, 0., 1.);
//  Variable* wg2 = new Variable("wg2",0.5 , 0., 1.);
//    Variable* signalYield = new Variable("signalYield"  , 1.15052e+03, 0.,0., 10000.);
//      Variable* bckgYield   = new Variable("bckgYield"  , 4.92744e+03, 0.,0., 20000.);
    Variable* signalYield = new Variable("signalYield", 1100, 0.1, 0., 10000.);
    Variable* bckgYield   = new Variable("bckgYield"  , 5000, 0.1, 0., 20000.);
//Variable* signalYield = new Variable("signalYield", 1.15052e+03);
//Variable* bckgYield   = new Variable("bckgYield"  , 4.92744e+03);

//   Variable* constaCoef = new Variable("constaCoef", 70, 0.001, 10, 100); 
//   Variable* linearCoef = new Variable("linearCoef", 0.1, 0.001, -0.35, 10.); 
//   Variable* secondCoef = new Variable("secondCoef", 0.1, 0.001, 0, 10);
//   Variable* thirdCoef  = new Variable("thirdCoef" , 0.1, 0.001, 0, 10);
 
//   Variable* constaCoef = new Variable("constaCoef", 70., XStepMinuit, 20., 1000); 
//   Variable* linearCoef = new Variable("linearCoef", 0.1, XStepMinuit, -3.5, 10.); 
//   Variable* secondCoef = new Variable("secondCoef", 0.1, XStepMinuit, 0, 10);
//   Variable* thirdCoef  = new Variable("thirdCoef" , 0.1, XStepMinuit, 0, 10);
  Variable* constaCoef = new Variable("constaCoef", 1. ,0.,1000. ); 
  Variable* linearCoef = new Variable("linearCoef", 0.001,0,10 ); 

  Variable* p0   = new Variable("p0",-2.32996e-01,0.001,-10.,10. ); 
  Variable* p1   = new Variable("p1", 0,0,-10,10 ); 
  Variable* xxmin = new Variable("xmin",XMinSign ); 
  Variable* xxmax = new Variable("xmax",XMaxSign); 
  SimpleCheby2Pdf* bckgMass  = new SimpleCheby2Pdf("bckgMass", xMass, p0, p1,xxmin,xxmax);
//  SimpleCheby2Pdf* Test  = new SimpleCheby2Pdf("Test", xMass, p0, p1,xxmin,xxmax);

// double  fullRange = XMaxSign - XMinSign;
// double  minScaled = -1. + 2. * (XMinSign - xminfull) / fullRange;
// 
// double  maxScaled = +1. - 2. * (xmaxfull - XMaxSign)) / fullRange; 

//  Variable* aslope     = new Variable("slope", -1.);
//  Variable* aslope     = new Variable("slope", 0.39, -10, 10);
//  Variable* apower     = new Variable("apower", 6, 0, 10);
//  Variable* apower     = new Variable("apower", 1.18, XStepMinuit, 0.9, 15.);
//  Variable* apower     = new Variable("apower", 1.18, XStepMinuit, 0.9, 6.);
//  Variable* apower     = new Variable("apower", 1.18, 0.001, 0.9, 5.);
//  Variable* treshold   = new Variable("treshold" ,5.168,XStepMinuit, 5.02, 6.);
//  Variable* treshold   = new Variable("treshold" ,5.33,0, 5.04, 6.);

 
 
//  std::vector<Variable*> weightsSignalMass;
//  weightsSignalMass.push_back(wg1);
//  weightsSignalMass.push_back(wg2);

//  std::vector<PdfBase*> compsSignalMass;
//  compsSignalMass.push_back(gauss1);
//  compsSignalMass.push_back(gauss2);
//  compsSignalMass.push_back(gauss3);

//  AddPdf signalMass("signalMass", weightsSignalMass, compsSignalMass); 
//  signalMass.addSpecialMask(PdfBase::ForceCommonNorm) ;
  

//  vector<Variable*> weightsPoly;
//  weightsPoly.push_back(constaCoef);
//  weightsPoly.push_back(linearCoef);
//  weightsPoly.push_back(secondCoef);
//  weightsPoly.push_back(thirdCoef);

  
//  PolynomialPdf* polyTmp = new PolynomialPdf("polyTmp", xMass, weightsPoly); 
//   std::vector<PdfBase*> compsPoly2;
//   compsSignalMass.push_back(polyTmp);
//   compsSignalMass.push_back(polyTmp);
//   
//   ProdPdf* poly      = new ProdPdf("poly"  ,compsPoly2 );
 
//  PolynomialPdf* poly = new PolynomialPdf("poly", xMass, weightsPoly); 
  
//  std::vector<Variable*> weightsBckgMass;
//  weightsBckgMass.push_back(wb1);
//  weightsBckgMass.push_back(wb2);
//weightsBckgMass.push_back(wb3);

//  ArgusPdf* argus = new  ArgusPdf("argus", xMass, treshold, aslope, true, apower);  

//  std::vector<PdfBase*> compsBckgMass;
//  compsBckgMass.push_back(gaussBckgBc);
//  compsBckgMass.push_back(gaussBckgB0);
//    compsBckgMass.push_back(argus);
//  compsBckgMass.push_back(poly);
//  compsBckgMass.push_back(Test);
 
//  AddPdf bckgMass("bckgMass", weightsBckgMass, compsBckgMass);
//  bckgMass.addSpecialMask(PdfBase::ForceCommonNorm) ;

  
//==============================================================================
//==============================================================================
//==============================================================================
// Lifetime
//==============================================================================
//==============================================================================
//==============================================================================

//   Variable* cTau     = new Variable("cTau"  ,1./(-1.638 *c_const), -1000., 1000.);
//   Variable* tauSB1   = new Variable("tauSB1",1./(-1.440 *c_const),-1000., 1000.);
//   Variable* tauSB2   = new Variable("tauSB2",1./(-1.600 *c_const),-1000., 1000.);

  Variable* cTau     = new Variable("cTau"  ,1./( 0.54  *c_const),0.01,0.00, 100.);
  Variable* tauSB1   = new Variable("tauSB1",1./( 0.113 *c_const),0.01,0.00, 1000.);
  Variable* tauSB2   = new Variable("tauSB2",1./( 0.487 *c_const),0.01,0.00, 1000.);

  Variable* meanRes      = new Variable("meanRes"        ,0,0,0,XMax);

  Variable* meanResBckg  = new Variable("meanResBckg"    ,0,0,0,XMax);

  Variable* sigmaRes     = new Variable("sigmaRes"      ,0.0004,0,2*SXMax);

  Variable* sigmaResBckg = new Variable("sigmaResBckg"  ,0.0004,0,2*SXMax);


  Variable* sigmaFix     = new Variable("sigmaFix"      ,0.001,0,0.01);
 
    Variable* meanGaussianErrSign    = new Variable( "meanGaussianErrSign"      ,1.72100e-03 ,0.0001, SXMax);
    Variable* sigmaGaussianErrorSign = new Variable( "sigmaGaussianErrorSign"   ,4.55213e-04 ,0.0001, SXMax);

    Variable* tauErrSign   = new Variable("tauErrSign",1.49315e+03,1000., 5100.);

    Variable* meanGaussianErrBckg1    = new Variable( "meanGaussianErrBckg1"	,1.78238e-03 ,0.0001, SXMax);
    Variable* sigmaGaussianErrorBckg1 = new Variable( "sigmaGaussianErrorBckg1"	,4.75277e-04 ,0.0001, SXMax);

    Variable* tauErrBckg1   = new Variable("tauErrBck1g",1.58241e+03,900., 8000.);

    Variable* meanGaussianErrBckg2    = new Variable( "meanGaussianErrBckg2"	,1.78238e-03 ,0.0001, SXMax);
    Variable* sigmaGaussianErrorBckg2 = new Variable( "sigmaGaussianErrorBckg2"	,4.75277e-04 ,0.0001, SXMax);

    Variable* tauErrBckg2   = new Variable("tauErrBckg2",1.58241e+03,1000., 5100.);


//Variable* meanGaussianErrSign    = new Variable( "meanGaussianErrSign"      ,1.72100e-03 ,0.00000, SXMax);
//Variable* sigmaGaussianErrorSign = new Variable( "sigmaGaussianErrorSign"   ,4.55213e-04 ,0.00001, SXMax);
//
//Variable* meanGaussianErrBckg    = new Variable( "meanGaussianErrBckg"      ,1.78238e-03 ,0.00000, SXMax);
//Variable* sigmaGaussianErrorBckg = new Variable( "sigmaGaussianErrorBckg"   ,4.75277e-04 ,0.00001, SXMax);
//
//Variable* tauErrSign   = new Variable("tauErrSign",1.49315e+03,100, 2100.);
//Variable* tauErrBckg   = new Variable("tauErrBckg",1.58241e+03,100, 2100.);

  Variable* gammaErrSign = new Variable( "gammaErrSign" ,5.3 ,0.00000, 10);
  Variable* gammaErrBckg = new Variable( "gammaErrBckg" ,5.3 ,0.00000, 10);

  Variable* betaErrSign  = new Variable( "betaErrSign"  ,0.0005 ,0.000001, 10);
  Variable* betaErrBckg  = new Variable( "betaErrBckg"  ,0.0005 ,0.000001, 10);

  Variable* muErrSign    = new Variable( "muErrSign"    ,0.0 );
  Variable* muErrBckg    = new Variable( "muErrBckg"    ,0.0 );

  

//  double unoParam =  1.         ;
//   double ef0Param =  5.06489e-02;
//   double ef1Param =  4.05836e-02; 
//   double ef2Param =  2.63305e-02; 
//   double ef3Param = -1.67310e+00; 
//   
	
	
  double ef0Param =  9.77208e-02;
  double ef1Param = -1.57397e-01;
//   double ef2Param =  4.45647e-01; 
//   double ef3Param = -7.82135e-01; 
  
  
  
  
//  Variable* ef0 = new Variable("ef0",  5.06489e-02); 
//  Variable* uno = new Variable("uno",  1); 
  Variable* ef0 = new Variable("ef0",  ef0Param); 
  Variable* ef1 = new Variable("ef1",  ef1Param); 
//  Variable* ef2 = new Variable("ef2",  ef2Param); 
//  Variable* ef3 = new Variable("ef3",  ef3Param); 

  vector<Variable*> coeffEffi;
  coeffEffi.push_back(ef0);
//  coeffEffi.push_back(uno);
  coeffEffi.push_back(ef1);
//  coeffEffi.push_back(ef2);
//  coeffEffi.push_back(ef3);

 
  PolynomialPdf* Effi = new PolynomialPdf("Effi", xcTau, coeffEffi); 
  
  Variable* XMinV     = new Variable("XMinV"  ,XMin ,0,0,1);
  Variable* XMaxV     = new Variable("XMaxV"  ,XMax ,0,0,1);
  Variable* SXMinV    = new Variable("SXMinV" ,SXMin,0,0,1);
  Variable* SXMaxV    = new Variable("SXMaxV" ,SXMax,0,0,1);

// ExpPdf* DecayBc     = new ExpPdf("DecayBc"	,  xcTau, cTau  );
// ExpPdf* pdfFitBckg1 = new ExpPdf("pdfFitBckg1", xcTau, tauSB1);
// ExpPdf* pdfFitBckg2 = new ExpPdf("pdfFitBckg2", xcTau, tauSB2);

//
//     ExpGausPEEPdf* DecayBc     = new ExpGausPEEPdf("DecayBc"    , xcTau, xScTau, meanRes    , cTau  );
//     ExpGausPEEPdf* pdfFitBckg1 = new ExpGausPEEPdf("pdfFitBckg1", xcTau, xScTau, meanResBckg, tauSB1);
//     ExpGausPEEPdf* pdfFitBckg2 = new ExpGausPEEPdf("pdfFitBckg2", xcTau, xScTau, meanResBckg, tauSB2);
//   ExpGausPdf* DecayBc     = new ExpGausPdf("DecayBc"    , xcTau, xScTau, meanRes    , cTau  );
//   ExpGausPdf* pdfFitBckg1 = new ExpGausPdf("pdfFitBckg1", xcTau, xScTau, meanResBckg, tauSB1);
//   ExpGausPdf* pdfFitBckg2 = new ExpGausPdf("pdfFitBckg2", xcTau, xScTau, meanResBckg, tauSB2);
 
     ExpGausPEEPdf* DecayBcProj	    = new ExpGausPEEPdf("DecayBcProj"	 , xcTau, xScTau, meanRes    , cTau  );
     ExpGausPEEPdf* pdfFitBckg1Proj = new ExpGausPEEPdf("pdfFitBckg1Proj", xcTau, xScTau, meanResBckg, tauSB1);
     ExpGausPEEPdf* pdfFitBckg2Proj = new ExpGausPEEPdf("pdfFitBckg2Proj", xcTau, xScTau, meanResBckg, tauSB2);
 
    ExpGausProdBPdf* DecayBc  = new ExpGausProdBPdf("DecayBc"    , xcTau, xScTau, meanRes    , cTau  , sigmaGaussianErrorSign, meanGaussianErrSign,tauErrSign,
    XMinV,XMaxV,SXMinV,SXMaxV);
    ExpGausProdBPdf* pdfFitBckg1 = new ExpGausProdBPdf("pdfFitBckg1", xcTau, xScTau, meanResBckg, tauSB1, sigmaGaussianErrorBckg1, meanGaussianErrBckg1,tauErrBckg1,
    XMinV,XMaxV,SXMinV,SXMaxV);
    ExpGausProdBPdf* pdfFitBckg2 = new ExpGausProdBPdf("pdfFitBckg2", xcTau, xScTau, meanResBckg, tauSB2, sigmaGaussianErrorBckg2, meanGaussianErrBckg2,tauErrBckg2,
    XMinV,XMaxV,SXMinV,SXMaxV);


//     ExpGausProdPdf* DecayBc  = new ExpGausProdPdf("DecayBc"    , xcTau, xScTau, meanRes    , cTau  , sigmaGaussianErrorSign, meanGaussianErrSign,tauErrSign);
//     ExpGausProdPdf* pdfFitBckg1 = new ExpGausProdPdf("pdfFitBckg1", xcTau, xScTau, meanResBckg, tauSB1, sigmaGaussianErrorBckg, meanGaussianErrBckg,tauErrBckg);
//     ExpGausProdPdf* pdfFitBckg2 = new ExpGausProdPdf("pdfFitBckg2", xcTau, xScTau, meanResBckg, tauSB2, sigmaGaussianErrorBckg, meanGaussianErrBckg,tauErrBckg);



//  ExpGausPEEfixSigmaPdf* DecayBc     = new ExpGausPEEfixSigmaPdf("DecayBc"    , xcTau, sigmaFix, meanRes    , cTau  );
//  ExpGausPEEfixSigmaPdf* pdfFitBckg1 = new ExpGausPEEfixSigmaPdf("pdfFitBckg1", xcTau, sigmaFix, meanResBckg, tauSB1);
//  ExpGausPEEfixSigmaPdf* pdfFitBckg2 = new ExpGausPEEfixSigmaPdf("pdfFitBckg2", xcTau, sigmaFix, meanResBckg, tauSB2);

//ExpGausMPdf* DecayBc	   = new ExpGausMPdf("DecayBc"	  , xcTau, xScTau,  meanRes    ,  cTau  );
//ExpGausMPdf* pdfFitBckg1 = new ExpGausMPdf("pdfFitBckg1", xcTau, xScTau,  meanResBckg,  tauSB1);
//ExpGausMPdf* pdfFitBckg2 = new ExpGausMPdf("pdfFitBckg2", xcTau, xScTau,  meanResBckg,  tauSB2);
  
//  ExpGausPdf* DecayBc     = new ExpGausPdf("DecayBc"    , xcTau, meanRes    , sigmaRes, cTau  );
//  ExpGausPdf* pdfFitBckg1 = new ExpGausPdf("pdfFitBckg1", xcTau, meanResBckg, sigmaResBckg, tauSB1);
//  ExpGausPdf* pdfFitBckg2 = new ExpGausPdf("pdfFitBckg2", xcTau, meanResBckg, sigmaResBckg, tauSB2);
//  ExpGausPdf* DecayBc	    = new ExpGausPdf("DecayBc"	  , xcTau, meanGaussianErrSign    , sigmaGaussianErrorSign, cTau  );
//  ExpGausPdf* pdfFitBckg1 = new ExpGausPdf("pdfFitBckg1", xcTau, meanGaussianErrBckg, sigmaGaussianErrorBckg, tauSB1);
//  ExpGausPdf* pdfFitBckg2 = new ExpGausPdf("pdfFitBckg2", xcTau, meanGaussianErrBckg, sigmaGaussianErrorBckg, tauSB2);

 
// Res Models   

//  ExpGausPdf* ExpGauSign = new ExpGausPdf("ExpGauSign" , xScTau, meanGaussianErrSign, sigmaGaussianErrorSign, tauErrSign);
//  ExpGausPdf* ExpGauBckg = new ExpGausPdf("ExpGauBckg" , xScTau, meanGaussianErrBckg, sigmaGaussianErrorBckg, tauErrBckg);

//    ExpGausPEESigmaPdf* ExpGauSign = new ExpGausPEESigmaPdf("ExpGauSig" , xScTau,  sigmaGaussianErrorSign, meanGaussianErrSign,tauErrSign);
//    ExpGausPEESigmaPdf* ExpGauBckg = new ExpGausPEESigmaPdf("ExpGauBckg", xScTau,  sigmaGaussianErrorBckg, meanGaussianErrBckg1,tauErrBckg1);

//    ExpGausWithIntPdf* ExpGauSign = new ExpGausWithIntPdf("ExpGauSig" , xScTau, meanGaussianErrSign, sigmaGaussianErrorSign, tauErrSign);
//    ExpGausWithIntPdf* ExpGauBckg = new ExpGausWithIntPdf("ExpGauBckg", xScTau, meanGaussianErrBckg, sigmaGaussianErrorBckg, tauErrBckg);

//    GammaPdf* GammaSign = new GammaPdf("GammaSig" , xScTau, gammaErrSign, betaErrSign, muErrSign);
//    GammaPdf* GammaBckg = new GammaPdf("GammaBckg", xScTau, gammaErrBckg, betaErrBckg, muErrBckg);

//    GaussianPdf* ExpGauSign = new GaussianPdf("ExpGauSign", xScTau, meanGaussianErrSign, sigmaGaussianErrorSign);
//    GaussianPdf* ExpGauBckg = new GaussianPdf("ExpGauBckg", xScTau, meanGaussianErrBckg, sigmaGaussianErrorBckg);
    
//  DecayBc     ->addSpecialMask(PdfBase::ForceSeparateNorm); 
//  pdfFitBckg1 ->addSpecialMask(PdfBase::ForceSeparateNorm); 
//  Effi        ->addSpecialMask(PdfBase::ForceSeparateNorm); 
  
  std::vector<PdfBase*> compspdfFitBc;
  compspdfFitBc.push_back(DecayBc);
  compspdfFitBc.push_back(Effi);

  ProdPdf* pdfFitBc     = new ProdPdf("pdfFitBc"  , compspdfFitBc);

  std::vector<PdfBase*> compspdfFitBcProj;
  compspdfFitBcProj.push_back(DecayBcProj);
  compspdfFitBcProj.push_back(Effi);

  ProdPdf* pdfFitBcProj     = new ProdPdf("pdfFitBcProj"  , compspdfFitBcProj);

  Variable* b1 = new Variable("b1",0.5,0.001,0., 1.);
  std::vector<Variable*> weightspdfFitBckg;
  weightspdfFitBckg.push_back(b1);

  std::vector<PdfBase*> compspdfFitBckgAdd;
  compspdfFitBckgAdd.push_back(pdfFitBckg1);
  compspdfFitBckgAdd.push_back(pdfFitBckg2);
  
  AddPdf pdfFitBckg("pdfFitBckg", weightspdfFitBckg, compspdfFitBckgAdd);

  std::vector<PdfBase*> compspdfFitBckgAddProj;
  compspdfFitBckgAddProj.push_back(pdfFitBckg1Proj);
  compspdfFitBckgAddProj.push_back(pdfFitBckg2Proj);
  
  AddPdf pdfFitBckgProj("pdfFitBckgProj", weightspdfFitBckg, compspdfFitBckgAddProj);
//  pdfFitBckg.addSpecialMask(PdfBase::ForceSeparateNorm) ;
//  AddPdf pdfFitBckgAdd("pdfFitBckgAdd", weightspdfFitBckg, compspdfFitBckgAdd);
  
 
//  std::vector<PdfBase*> compspdfFitBckg;
//compspdfFitBckg.push_back(pdfFitBckg1);
//    compspdfFitBckg.push_back(&pdfFitBckgAdd);
//  compspdfFitBckg.push_back(Effi);
 
//  ProdPdf* pdfFitBc     = new ProdPdf("pdfFitBc"  , compspdfFitBc);
//  ProdPdf* pdfFitBckg   = new ProdPdf("pdfFitBckg", compspdfFitBckg);


///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//
// 2DFit
//  
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

/* 
   signalMass->setIntegrationFineness(NFinessMass);
   bckgMass->setIntegrationFineness(NFinessMass);
   DecayBc    ->setIntegrationFineness(NFiness);
   pdfFitBckg1->setIntegrationFineness(NFiness);
   pdfFitBckg2->setIntegrationFineness(NFiness);
   ExpGauSign->setIntegrationFineness(NFiness);
   ExpGauBckg->setIntegrationFineness(NFiness);
   pdfFitBc->setIntegrationFineness(NFiness);
   pdfFitBckg.setIntegrationFineness(NFiness);
 */

/*   std::vector<PdfBase*> compsSignalTau;
  compsSignalTau.push_back(pdfFitBc);
  compsSignalTau.push_back(ExpGauSign);
//  compsSignalLife.push_back(DecayBc);
 
  std::vector<PdfBase*> compsBckgTau;
  compsBckgTau.push_back(&pdfFitBckg);
  compsSignalTau.push_back(ExpGauBckg);

  ProdPdf* signalTau = new ProdPdf("signalTau", compsSignalTau);
  ProdPdf* bckgTau   = new ProdPdf("bckgTau  ", compsBckgTau);
 */
  std::vector<PdfBase*> compsSignalLife;
 
  compsSignalLife.push_back(signalMass);
  compsSignalLife.push_back(pdfFitBc);
//  compsSignalLife.push_back(ExpGauSign);
//  compsSignalLife.push_back(ExpGauBckg);
 // compsSignalLife.push_back(GammaSign);/
//  compsSignalLife.push_back(signalTau);
//  compsSignalLife.push_back(DecayBc);
//   signalMass     ->addSpecialMask(PdfBase::ForceSeparateNorm); 
//   pdfFitBc       ->addSpecialMask(PdfBase::ForceSeparateNorm); 
//   ExpGauSign     ->addSpecialMask(PdfBase::ForceSeparateNorm); 
// 


  std::vector<PdfBase*> compsBckgLife;
  compsBckgLife.push_back(bckgMass);
  compsBckgLife.push_back(&pdfFitBckg);
//  compsSignalLife.push_back(ExpGauBckg);
  //compsSignalLife.push_back(GammaBckg);
  //compsSignalLife.push_back(ExpGauSign);
//  compsBckgLife.push_back(pdfFitBckg1);
//  compsBckgLife.push_back(bckgTau);
//compsBckgLife.push_back(&pdfFitBckgAdd);

//  bckgMass     ->addSpecialMask(PdfBase::ForceSeparateNorm); 
//  pdfFitBckg.addSpecialMask(PdfBase::ForceSeparateNorm); 
//  ExpGauBckg   ->addSpecialMask(PdfBase::ForceSeparateNorm); 

  ProdPdf* signalLife = new ProdPdf("signalLife", compsSignalLife);
  ProdPdf* bckgLife   = new ProdPdf("bckgLife  ", compsBckgLife);
//  signalLife->addSpecialMask(PdfBase::ForceCommonNorm) ;
//  bckgLife->addSpecialMask(PdfBase::ForceCommonNorm) ;

  std::vector<Variable*> weightsYield;
  weightsYield.push_back(signalYield);
  weightsYield.push_back(bckgYield);
  

  std::vector<PdfBase*> compsModel;
  
  
  compsModel.push_back(signalLife);
//  compsModel.push_back(gaussBckgBc);
  compsModel.push_back(bckgLife);
  
  
//  compsModel.push_back(poly);
//  compsModel.push_back(argus);
  AddPdf model("model", weightsYield, compsModel); 
// 
//   AddPdf model1("model1", weightsYield, compsModel); 
//   
//   std::vector<PdfBase*> compsModel1;
//   compsModel1.push_back(&model1);
//   compsModel1.push_back(ExpGauBckg);
//   
//   ProdPdf model(" model  ", compsModel1);
  
  
//  model.addSpecialMask(PdfBase::ForceSeparateNorm) ;
  

//  model.addSpecialMask(PdfBase::ForceCommonNorm) ;
 
//
// These are used for Plots....
//
  std::vector<PdfBase*> compsMass;
  compsMass.push_back(signalMass);
  compsMass.push_back(bckgMass);
  AddPdf modelMass("modelMass"  , weightsYield, compsMass); 
  
  std::vector<PdfBase*> compscTau;
//  compscTau.push_back(pdfFitBckg1);
//  compscTau.push_back(pdfFitBcProj);
//  compscTau.push_back(&pdfFitBckgProj);
  compscTau.push_back(pdfFitBc);
  compscTau.push_back(&pdfFitBckg);
  AddPdf model_cTau("model_cTau", weightsYield, compscTau); 
  
  //std::vector<PdfBase*> compsSTau;
  //   compsSTau.push_back(ExpGauSign);
  //   compsSTau.push_back(ExpGauBckg);
//compsSTau.push_back(GammaSign);
//compsSTau.push_back(GammaBckg);
  //AddPdf model_STau("model_cTau", weightsYield, compsSTau); 
  
  

//
// Data
//
  vector<Variable*> dataVec;
  
  dataVec.push_back(xMass);
  dataVec.push_back(xcTau);
  dataVec.push_back(xScTau);
  UnbinnedDataSet* dataLife = new UnbinnedDataSet(dataVec);

//  UnbinnedDataSet* dataSLife = new UnbinnedDataSet(xScTau);
//
  if (!InputFile)
   {
     cout<<"File:"<<InputFileName<<" not found!!!"<<endl;
    exit(1);
   }
   InputFile->ls();
   
   TTree *TauBcTree    = (TTree*)InputFile->Get(InputTauBcTreeName);
   if(!TauBcTree ){
     cout<<"TTree cTau Data: "<< InputTauBcTreeName <<" not found!!!"<<endl;
     exit(1);
   }else{
     cout<<"TTree cTau Data: "<< InputTauBcTreeName <<" OK FOUND!!!"<<endl;
   }  
    
   TauBcTree->SetBranchAddress("xBcMass",&xBcMass);
//   TauBcTree->SetBranchAddress("xBcTau" ,&xBcTau);
   TauBcTree->SetBranchAddress("xBccTau",&xBccTau);
   TauBcTree->SetBranchAddress("xSBccTau",&xSBccTau);
   int nentries = (int)TauBcTree->GetEntries();
   for (Int_t i=0;i<nentries;i++) {
    TauBcTree->GetEntry(i);
    if(xBccTau>XMin&&xBccTau<XMax&&xSBccTau>SXMin&&xSBccTau<SXMax){
     if(xBcMass>XMinSign&&xBcMass<XMaxSign){
      xMass->value  = xBcMass;
      xcTau->value  = xBccTau;
      xScTau->value = xSBccTau; 
      dataLife->addEvent();
      HxMass.Fill(xBcMass);
      HxcTau.Fill(xBccTau);
      HxScTau.Fill(xSBccTau);
//      printf("Mass = %f xcTau = %f sigma = %f \n", xBcMass, xBccTau, xSBccTau);     
     } 
     if(xBcMass>XMinSBL&&xBcMass<XMaxSBL){
      HxcTauSB.Fill(xBccTau);
      HxScTauSB.Fill(xSBccTau);
//      dataSLife->addEvent();
     } 
     if(xBcMass>XMinSBR&&xBcMass<XMaxSBR){
      HxcTauSB.Fill(xBccTau);
      HxScTauSB.Fill(xSBccTau);
//      dataSLife->addEvent();
     } 
    } 
   }
//================================================================================
//================================================================================
///FIT
//================================================================================
//================================================================================
    model.setData(dataLife);
    FitManager fitter(&model);
//ExpGauBckg->setData(dataSLife);
//FitManager fitter(ExpGauBckg);
  fitter.setMaxCalls(20000);
  std::cout<<"========================================================" <<std::endl;
  std::cout<<"======       	  Start Fit 	  		======"<<std::endl;
  std::cout<<"======       	  Start Fit 	  		======"<<std::endl;
  std::cout<<"======       	  Start Fit 	  		======"<<std::endl;
  std::cout<<"========================================================" <<std::endl;
  fitter.fit(); 
//  fitter.runMigrad(); 
  //TMinuit * Minuit = fitter.getMinuitObject();
  //Minuit->SetPrintLevel(1);
//  Minuit->mnmigr();
//  Minuit->mnhess();
//  Minuit->mnmigr();
  fitter.getMinuitValues(); 
  std::cout<<"========================================================"<<std::endl;
  std::cout<<"======	   	   End  Fit 	  		======"<<std::endl;
  std::cout<<"======	   	   End  Fit 	  		======"<<std::endl;
  std::cout<<"======	   	   End  Fit 	  		======"<<std::endl;
  std::cout<<"========================================================"<<std::endl;
  
//================================================================================
//================================================================================
///FIT
//================================================================================
//================================================================================
//exit(1);
//std::cout<<" warning!!!" <<std::endl;
 //================================================================================
///PLOT

// Mass
  UnbinnedDataSet gridMass(xMass);
  double totalDataMass = 0; 
  double NStep = XHScale*xMass->numbins;
  for (int i = 0; i < NStep; ++i) {
    double step = (xMass->upperlimit - xMass->lowerlimit)/NStep;
    xMass->value = xMass->lowerlimit + (i + 0.5) * step;
    gridMass.addEvent(); 
   totalDataMass++; 
  }

  modelMass.setData(&gridMass);
  vector<vector<double> > pdfValsMass;
  modelMass.getCompProbsAtDataPoints(pdfValsMass); 
  double totalPdfMass = 0; 
  for (int i = 0; i < gridMass.getNumEvents(); ++i) {
    gridMass.loadEvent(i); 
    pdfHist.Fill(xMass->value, pdfValsMass[0][i]);
    sigHist.Fill(xMass->value, pdfValsMass[1][i]);
    bkgHist.Fill(xMass->value, pdfValsMass[2][i]);
    totalPdfMass += pdfValsMass[0][i]; 
  }
  
  
  pdfHist.Scale((signalYield->value+bckgYield->value)/pdfHist.Integral()*XHScale);
  sigHist.Scale(signalYield->value/sigHist.Integral()*XHScale);
  bkgHist.Scale(bckgYield->value/bkgHist.Integral()*XHScale);
  std::cout<<"Signal Yield = "<< signalYield->value<<std::endl;
  std::cout<<"Bckg   Yield = "<< bckgYield->value<<std::endl;
  std::cout<<"(SB    Yield  = "<<HxcTauSB.GetEntries() <<")"<<std::endl;
  std::cout<<"Tot   Yield  = "<< signalYield->value+bckgYield->value<<std::endl;
//--------------------------------------------------  
// Tau


  int NIntegral = 1;
  

  vector<Variable*> dataPlot2D;
  dataPlot2D.push_back(xcTau);
  dataPlot2D.push_back(xScTau);
//  UnbinnedDataSet grid_cTau(dataPlot);
  UnbinnedDataSet grid_cTau2D(dataPlot2D);
//  UnbinnedDataSet grid_STau(dataPlotS);
  
//  bool first = true;
//  UnbinnedDataSet grid_cTau(xcTau);
//  double totalData_cTau = 0; 
  NStep  = XHScale*xcTau->numbins;
//  double NSStep   = XHScale*xcTau->numbins;
  double NSStep2D = XHScale*NIntegral*xScTau->numbins;
  double step  = (xcTau->upperlimit - xcTau->lowerlimit)/NStep;
//  double sstep = (xScTau->upperlimit - xScTau->lowerlimit)/NSStep;
  double sstep2D = (xScTau->upperlimit - xScTau->lowerlimit)/NSStep2D;
  for (int i = 0; i < NStep; ++i) {
    xcTau->value  = xcTau->lowerlimit  + (i + 0.5) * step;
//    grid_cTau.addEvent(); 
//    totalData_cTau++; 
//    xScTau->value = xScTau->lowerlimit + (i + 0.5) * sstep;
//    cout<<"X = "<<xcTau->value<<" sx = "<<xScTau->value<<endl;
//    grid_cTau.addEvent(); 
//    xcTau2D->value  = xcTau2D ->lowerlimit + (i + 0.5) * step;
//   cout<<"======================================     \n"<<endl;
//   cout<<"======================================     \n"<<endl;
//   cout<<"======================================     \n"<<endl;
    for (int ii = 0; ii < NSStep2D; ++ii) {
     xScTau->value = xScTau->lowerlimit + (ii + 0.5) * sstep2D;
//     xScTau2D->value = xScTau2D->lowerlimit + (ii + 0.5) * sstep2D;
//    cout<<"X = "<<xcTau->value<<" sx = "<<xScTau->value<<endl;
     grid_cTau2D.addEvent(); 
//     if (first) grid_STau.addEvent();
    }
//    first = false;
  }

//  model_cTau.setData(&grid_cTau);
//   vector<vector<double> > pdfVals_cTau;
//   model_cTau.getCompProbsAtDataPoints(pdfVals_cTau); 
//   double totalPdf_cTau = 0; 
//   for (int i = 0; i < grid_cTau.getNumEvents(); ++i) {
//     grid_cTau.loadEvent(i); 
//     pdf_cTau_Hist.Fill(xcTau->value , pdfVals_cTau[0][i]);
//     sig_cTau_Hist.Fill(xcTau->value , pdfVals_cTau[1][i]);
//     bkg_cTau_Hist.Fill(xcTau->value , pdfVals_cTau[2][i]);
//     totalPdf_cTau += pdfVals_cTau[0][i]; 
//   }
  
//   double pdf_cTau_Integral2D = 0;
//   double sig_cTau_Integral2D = 0;
//   double bkg_cTau_Integral2D = 0;
//  int NStep2D = NStep*NIntegral;
  model_cTau.setData(&grid_cTau2D);
  vector<vector<double> > pdfVals_cTau2D;
  model_cTau.getCompProbsAtDataPoints(pdfVals_cTau2D); 
  for (int i = 0; i < grid_cTau2D.getNumEvents(); ++i) {
    grid_cTau2D.loadEvent(i); 
    pdf_cTauSTau_Hist2D.Fill(xcTau->value ,xScTau->value , pdfVals_cTau2D[0][i]);
    sig_cTauSTau_Hist2D.Fill(xcTau->value ,xScTau->value , pdfVals_cTau2D[1][i]);
    bkg_cTauSTau_Hist2D.Fill(xcTau->value ,xScTau->value , pdfVals_cTau2D[2][i]);
//     if (i%NStep2D == 1 && i>0){
//      pdf_cTau_Hist2D.Fill(xcTau->value , pdf_cTau_Integral2D/step);
//      sig_cTau_Hist2D.Fill(xcTau->value , sig_cTau_Integral2D/step);
//      bkg_cTau_Hist2D.Fill(xcTau->value , bkg_cTau_Integral2D/step);
// //     cout<<"Int = "<<bkg_cTau_Integral2D<<endl;
//      pdf_cTau_Integral2D=0;
//      sig_cTau_Integral2D=0;
//      bkg_cTau_Integral2D=0;
// //     exit(0);
//     }else{
// //     cout<<"X = "<<xcTau->value<<" sx = "<<xScTau->value<<" NStep2D = "<<NStep2D<<endl;
//      pdf_cTau_Integral2D =+ pdfVals_cTau2D[0][i]*sstep2D; 
//      sig_cTau_Integral2D =+ pdfVals_cTau2D[1][i]*sstep2D; 
//      bkg_cTau_Integral2D =+ pdfVals_cTau2D[2][i]*sstep2D;
// //     cout<<"Int = "<<pdfVals_cTau2D[2][i]<<endl;
//     } 
  }

  TH1D * pdf_cTauSTau_X = pdf_cTauSTau_Hist2D.ProjectionX("pdf_cTauSTau_X");
  TH1D * pdf_cTauSTau_Y = pdf_cTauSTau_Hist2D.ProjectionY("pdf_cTauSTau_Y");

  TH1D * sig_cTauSTau_X = sig_cTauSTau_Hist2D.ProjectionX("sig_cTauSTau_X");
  TH1D * sig_cTauSTau_Y = sig_cTauSTau_Hist2D.ProjectionY("sig_cTauSTau_Y");

  TH1D * bkg_cTauSTau_X = bkg_cTauSTau_Hist2D.ProjectionX("bkg_cTauSTau_X");
  TH1D * bkg_cTauSTau_Y = bkg_cTauSTau_Hist2D.ProjectionY("bkg_cTauSTau_Y");
/*   vector<Variable*> dataSPlot;
  dataSPlot.push_back(xScTau);
  dataSPlot.push_back(xcTau);
  UnbinnedDataSet grid_STau(dataSPlot);
//  UnbinnedDataSet grid_cTau(xcTau);
  double totalData_STau = 0; 
  NStep  = XHScale*xScTau->numbins;
  double NSStep = XHScale*xcTau->numbins;
  for (int i = 0; i < NSStep; ++i) {
    totalData_STau++; 
    double sstep = (xScTau->upperlimit - xScTau->lowerlimit)/NSStep;
    xScTau->value = xScTau->lowerlimit + (i + 0.5) * sstep;
    grid_STau.addEvent(); 
  }
 */

//
// STau
//
//   model_STau.setData(&grid_cTau);
//   vector<vector<double> > pdfVals_STau;
//   model_STau.getCompProbsAtDataPoints(pdfVals_STau); 
//   double totalPdf_STau = 0; 
//   for (int i = 0; i < grid_cTau.getNumEvents(); ++i) {
//     grid_cTau.loadEvent(i); 
//     pdf_STau_Hist.Fill(xScTau->value, pdfVals_STau[0][i]);
//     sig_STau_Hist.Fill(xScTau->value, pdfVals_STau[1][i]);
//     bkg_STau_Hist.Fill(xScTau->value, pdfVals_STau[2][i]);
//     totalPdf_STau += pdfVals_STau[0][i]; 
//   }
  
//
// Models plot  
  pdf_cTauSTau_X->Scale((signalYield->value+bckgYield->value)/pdf_cTauSTau_X->Integral()*XHScale);
  pdf_cTauSTau_Y->Scale((signalYield->value+bckgYield->value)/pdf_cTauSTau_Y->Integral()*XHScale);
  sig_cTauSTau_X->Scale((signalYield->value)/sig_cTauSTau_X->Integral()*XHScale);
  sig_cTauSTau_Y->Scale((signalYield->value)/sig_cTauSTau_Y->Integral()*XHScale);
  bkg_cTauSTau_X->Scale((bckgYield->value)/bkg_cTauSTau_X->Integral()*XHScale);
  bkg_cTauSTau_Y->Scale((bckgYield->value)/bkg_cTauSTau_Y->Integral()*XHScale);

//   pdf_cTau_Hist.Scale((signalYield->value+bckgYield->value)/pdf_cTau_Hist.Integral()*XHScale);
//   sig_cTau_Hist.Scale(signalYield->value/sig_cTau_Hist.Integral()*XHScale);
//   bkg_cTau_Hist.Scale(HxcTauSB.GetEntries()/bkg_cTau_Hist.Integral()*XHScale);
//   
//   pdf_cTau_Hist2D.Scale((signalYield->value+bckgYield->value)/pdf_cTau_Hist2D.Integral()*XHScale);
//   sig_cTau_Hist2D.Scale(signalYield->value/sig_cTau_Hist2D.Integral()*XHScale);
//   bkg_cTau_Hist2D.Scale(HxcTauSB.GetEntries()/bkg_cTau_Hist2D.Integral()*XHScale);
//   
//   pdf_STau_Hist.Scale((signalYield->value+bckgYield->value)/pdf_STau_Hist.Integral()*XHScale);
//   sig_STau_Hist.Scale(signalYield->value/sig_STau_Hist.Integral()*XHScale);
//   bkg_STau_Hist.Scale(HxScTauSB.GetEntries()/bkg_STau_Hist.Integral()*XHScale);
//  bkg_cTau_Hist.Scale(bckgYield->value/bkg_cTau_Hist.Integral()*XHScale);
 
//   for (int i = 0; i < xMass->numbins; ++i) {
//     double val = pdfHist.GetBinContent(i+1); 
//     val /= totalPdf; 
//     val *= totalData;
//     pdfHist.SetBinContent(i+1, val); 
//     val = sigHist.GetBinContent(i+1); 
//     val /= totalPdf; 
//     val *= sigFrac->value; 
//     val *= totalData;
//     sigHist.SetBinContent(i+1, val); 
//     val = bkgHist.GetBinContent(i+1); 
//     val /= totalPdf; 
//     val *= (1.0 - sigFrac->value);
//     val *= totalData;
//     bkgHist.SetBinContent(i+1, val); 
//   }
  c1->cd();
    TLegend* leg_sign = new TLegend(0.30,0.75,0.90,0.90);
    leg_sign->SetTextSize(0.025) ;
    leg_sign->SetTextAlign(31);
    leg_sign->SetBorderSize(0.);
    leg_sign->SetFillStyle(0);
    leg_sign->SetHeader("B^{+} mass spectrum     ");
    if(signalYield->error!=0){
      leg_sign->AddEntry(&HxMass ,Form( "Yield_{Sign} =    %5.0f  #pm %5.0f",signalYield->value,signalYield->error),"");
    }else{
      leg_sign->AddEntry(&HxMass ,Form( "Yield_{Sign} =    %5.0f Fixed",signalYield->value),"");
    }
    if(bckgYield->error!=0){
      leg_sign->AddEntry(&HxMass ,Form( "Yield_{Bckg} =   %5.0f  #pm %5.0f",bckgYield->value,bckgYield->error),"");
    }else{
      leg_sign->AddEntry(&HxMass ,Form( "Yield_{Bckg} =   %5.0f Fixed",bckgYield->value),"");
    }
    
    if(mean->error!=0){
     leg_sign->AddEntry(&HxMass ,Form( "M_{B_{c}^{+}} =   %5.4f  #pm %5.4f",mean->value,mean->error),"");
    }else{
     leg_sign->AddEntry(&HxMass ,Form( "M_{B_{c}^{+}} =   %5.4f Fixed",mean->value),"");
     }
    if(sigma->error!=0){
     leg_sign->AddEntry(&HxMass ,Form( "#sigma_{B_{c}^{+}} =   %5.4f  #pm %5.4f",sigma->value,sigma->error),"");
    }else{
     leg_sign->AddEntry(&HxMass ,Form( "#sigma_{B_{c}^{+}} =   %5.4f Fixed",sigma->value),"");
    }
  HxMass.SetMarkerStyle(8);
  HxMass.SetMarkerSize(0.5);
  HxMass.SetTitle("");
  HxMass.Draw("E1"); 
  leg_sign->Draw("same");
//  HxMass.Draw("p"); 
  pdfHist.SetLineColor(kBlue);
  pdfHist.SetLineWidth(3); 
  pdfHist.Draw("same"); 
  sigHist.SetLineColor(kMagenta);
  sigHist.SetLineStyle(kDashed); 
  sigHist.SetLineWidth(2); 
  sigHist.Draw("same"); 
  bkgHist.SetLineColor(kRed);
  bkgHist.SetLineStyle(kDashed); 
  bkgHist.SetLineWidth(2); 
  bkgHist.Draw("same"); 
  HxMass.Write();
  pdfHist.Write();
  sigHist.Write();
  bkgHist.Write();
//  
  c2->cd();
  c2->SetLogy();
  TLegend* leg_pdfSB = new TLegend(0.60,0.65,0.90,0.90);
  leg_pdfSB->SetTextAlign(12);
  leg_pdfSB->SetHeader("B^{+} proper time Fit Projections");
  leg_pdfSB->SetTextSize(0.025) ;
  leg_pdfSB->SetBorderSize(0.);
  leg_pdfSB->SetFillStyle(0);
  leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[4]{#scale[1.5]{#tau}_{B^{+}}  =  %5.3f #pm %5.3f     }",1/(c_const*cTau->value),cTau->error/((c_const*cTau->value)*(cTau->value)))   ,"");
  if( b1->error!=0){
      leg_pdfSB->AddEntry(&HxcTau ,Form( "b1   =  %5.3f #pm %5.3f     ",b1->value,b1->error)   ,"");
  }else{      
      leg_pdfSB->AddEntry(&HxcTau ,Form( "b1   =  %5.3f     Fixed     ",b1->value)   ,"");
  }   
  if( tauSB1->error!=0){
      leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB1} =  %5.3f #pm %5.3f     }",1/(c_const*tauSB1->value),tauSB1->error/((c_const*tauSB1->value)*(tauSB1->value)))   ,"");
  }else{      
      leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB1} =  %5.3f	Fixed	  }",1/(c_const*tauSB1->value))   ,"");
  }   
  if( tauSB2->error!=0){
      leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB2} =  %5.3f #pm %5.3f     }",1/(c_const*tauSB2->value),tauSB2->error/((c_const*tauSB2->value)*(tauSB2->value)))   ,"");
  }else{      
      leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB2} =  %5.3f	Fixed	  }",1/(c_const*tauSB2->value))   ,"");
  }   
  leg_pdfSB->AddEntry(&HxcTau ,"#color[4]{#scale[0.8]{- Fit model	     }}"   ,"");
  leg_pdfSB->AddEntry(&HxcTau ,"#color[6]{#scale[0.8]{- Signal model	     }}"   ,"");
  leg_pdfSB->AddEntry(&HxcTau ,"#color[2]{#scale[0.8]{- Background model on SB}}"   ,"");
  HxcTau.SetMarkerStyle(8);
  HxcTau.SetMarkerSize(0.5);
  HxcTau.SetTitle("");
  HxcTau.Draw("E1");
  HxcTauSB.SetMarkerStyle(8);
  HxcTauSB.SetMarkerSize(0.5);
  HxcTauSB.SetMarkerColor(2);
  HxcTauSB.Draw("same,E1");
  leg_pdfSB->Draw("same");
  pdf_cTauSTau_X->SetLineColor(kBlue);
  pdf_cTauSTau_X->SetLineWidth(PlotLineWidth);
  pdf_cTauSTau_X->Draw("same");
  sig_cTauSTau_X->SetLineColor(kMagenta);
  sig_cTauSTau_X->SetLineWidth(PlotLineWidth);
  sig_cTauSTau_X->SetLineStyle(kDashed);
  sig_cTauSTau_X->Draw("same");
  bkg_cTauSTau_X->SetLineColor(kRed);
  bkg_cTauSTau_X->SetLineWidth(PlotLineWidth);
  bkg_cTauSTau_X->SetLineStyle(kDashed);
  bkg_cTauSTau_X->Draw("same");
 
//   pdf_cTau_Hist2D.SetLineColor(kBlue);
//   pdf_cTau_Hist2D.SetLineWidth(3); 
//   pdf_cTau_Hist2D.Draw("same"); 
//   sig_cTau_Hist2D.SetLineColor(kMagenta);
//   sig_cTau_Hist2D.SetLineStyle(kDashed); 
//   sig_cTau_Hist2D.SetLineWidth(2); 
//   sig_cTau_Hist2D.Draw("same"); 
//   bkg_cTau_Hist2D.SetLineColor(kRed);
//   bkg_cTau_Hist2D.SetLineStyle(kDashed); 
//   bkg_cTau_Hist2D.SetLineWidth(2); 
//   bkg_cTau_Hist2D.Draw("same"); 
  
//   pdf_cTau_Hist.SetLineColor(kBlue);
//   pdf_cTau_Hist.SetLineWidth(3); 
//   pdf_cTau_Hist.Draw("same"); 
//   sig_cTau_Hist.SetLineColor(kMagenta);
//   sig_cTau_Hist.SetLineStyle(kDashed); 
//   sig_cTau_Hist.SetLineWidth(2); 
//   sig_cTau_Hist.Draw("same"); 
//   bkg_cTau_Hist.SetLineColor(kRed);
//   bkg_cTau_Hist.SetLineStyle(kDashed); 
//   bkg_cTau_Hist.SetLineWidth(2); 
//   bkg_cTau_Hist.Draw("same"); 
  
  c3->cd();
  TLegend* leg_pdfResolution = new TLegend(0.40,0.65,0.90,0.90);
  leg_pdfResolution->SetTextAlign(12);
  leg_pdfResolution->SetHeader("B_{c}^{+} resolution Fit Projections");
  leg_pdfResolution->SetTextSize(0.025) ;
  leg_pdfResolution->SetBorderSize(0.);
  leg_pdfResolution->SetFillStyle(0);
//  leg_pdfResolution->AddEntry(&HxcTau ,Form( "#color[4]{#scale[1.5]{#tau}_{B^{+}}  =  %5.3f #pm %5.3f     }",1/(c_const*cTau->value),cTau->error/((c_const*cTau->value)*(cTau->value)))   ,"");
//   if( b1->error!=0){
//       leg_pdfResolution->AddEntry(&HxcTau ,Form( "b1   =  %5.3f #pm %5.3f     ",b1->value,b1->error)   ,"");
//   }else{      
//       leg_pdfResolution->AddEntry(&HxcTau ,Form( "b1   =  %5.3f     Fixed     ",b1->value)   ,"");
//   }   
//   if( tauResolution1->error!=0){
//       leg_pdfResolution->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{Resolution1} =  %5.3f #pm %5.3f     }",1/(c_const*tauResolution1->value),tauResolution1->error/((c_const*tauResolution1->value)*(tauResolution1->value)))   ,"");
//   }else{      
//       leg_pdfResolution->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{Resolution1} =  %5.3f	Fixed	  }",1/(c_const*tauResolution1->value))   ,"");
//   }   
//   if( tauResolution2->error!=0){
//       leg_pdfResolution->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{Resolution2} =  %5.3f #pm %5.3f     }",1/(c_const*tauResolution2->value),tauResolution2->error/((c_const*tauResolution2->value)*(tauResolution2->value)))   ,"");
//   }else{      
//       leg_pdfResolution->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{Resolution2} =  %5.3f	Fixed	  }",1/(c_const*tauResolution2->value))   ,"");
//   }   
  leg_pdfResolution->AddEntry(&HxScTau ,"#color[4]{#scale[0.8]{- Pdf model Resolution	     }}"   ,"");
  leg_pdfResolution->AddEntry(&HxScTau ,"#color[6]{#scale[0.8]{- Signal model Resolution     }}"   ,"");
  leg_pdfResolution->AddEntry(&HxScTau ,"#color[2]{#scale[0.8]{- Background model  Resolution (on SB events)}}"   ,"");
  HxScTau.Draw("E1");
  leg_pdfResolution->Draw("same");
  HxScTauSB.SetMarkerStyle(8);
  HxScTauSB.SetMarkerSize(0.5);
  HxScTauSB.SetMarkerColor(kRed);
  HxScTauSB.Draw("same,E1");
  pdf_cTauSTau_Y->SetLineColor(kBlue);
  pdf_cTauSTau_Y->SetLineWidth(PlotLineWidth);
  pdf_cTauSTau_Y->Draw("same");
  sig_cTauSTau_Y->SetLineColor(kMagenta);
  sig_cTauSTau_Y->SetLineWidth(PlotLineWidth);
  sig_cTauSTau_Y->SetLineStyle(kDashed);
  sig_cTauSTau_Y->Draw("same");
  bkg_cTauSTau_Y->SetLineColor(kRed);
  bkg_cTauSTau_Y->SetLineWidth(PlotLineWidth);
  bkg_cTauSTau_Y->SetLineStyle(kDashed);
  bkg_cTauSTau_Y->Draw("same");
  
  
  
  
  
  HxcTau.Write();
  HxcTauSB.Write();
  HxScTau.Write();
  HxScTauSB.Write();
  pdf_cTauSTau_X->Write();
  pdf_cTauSTau_Y->Write();
  sig_cTauSTau_X->Write();
  sig_cTauSTau_Y->Write();
  bkg_cTauSTau_X->Write();
  bkg_cTauSTau_Y->Write();
//   pdf_cTau_Hist2D.Write();
//   sig_cTau_Hist2D.Write();
//   bkg_cTau_Hist2D.Write();
//   bkg_cTau_Hist2D.Write();
//   pdf_STau_Hist.Write();
//   sig_STau_Hist.Write();
//   bkg_STau_Hist.Write();
  c1->Write();
  c2->Write();
  c3->Write();
  char PDFNameMass[50] = "Bc-Mass-2016.pdf";
  char PDFNamecTau[50] = "Bc-cTau-2016.pdf";
  char PDFNameReso[50] = "Bc-Reso-2016.pdf";
  char testo[130] ;
  sprintf(testo,"mv %s %s.tmp",PDFNameMass,PDFNameMass);
  gSystem->Exec(testo);
  sprintf(testo,"mv %s %s.tmp",PDFNamecTau,PDFNamecTau);
  gSystem->Exec(testo);
  sprintf(testo,"mv %s %s.tmp",PDFNameReso,PDFNameReso);
  gSystem->Exec(testo);

  std::cout<<"Tau    = "<<1/(c_const*cTau->value)<<"+/-"<<cTau->error/((c_const*cTau->value)*(cTau->value))<<std::endl;
  std::cout<<"TauSB1 = "<<1/(c_const*tauSB1->value)<<"+/-"<<tauSB1->error/((c_const*tauSB1->value)*(tauSB1->value))<<std::endl;
  std::cout<<"TauSB2 = "<<1/(c_const*tauSB2->value)<<"+/-"<<tauSB2->error/((c_const*tauSB2->value)*(tauSB2->value))<<std::endl;
  
  c1->Print(PDFNameMass);
  c2->Print(PDFNamecTau);
  c3->Print(PDFNameReso);
  OutFile->Close();
  std::cout<<"==========================================" <<std::endl;
  std::cout<<"==========================================" <<std::endl;

  stopCPU = times(&stopProc);
  gettimeofday(&stopTime, NULL);
  // Print total minimization time
  double myCPU = stopCPU - startCPU;
  double totalCPU = myCPU; 

  timersub(&stopTime, &startTime, &totalTime);
  std::cout << "Wallclock time  : " << totalTime.tv_sec + totalTime.tv_usec/1000000.0 << " seconds." << std::endl;
  std::cout << "CPU time: " << (myCPU / CLOCKS_PER_SEC) << std::endl; 
  std::cout << "Total CPU time: " << (totalCPU / CLOCKS_PER_SEC) << std::endl; 
  myCPU = stopProc.tms_utime - startProc.tms_utime;
  std::cout << "Processor time: " << (myCPU / CLOCKS_PER_SEC) << std::endl;
  std::cout<<"==========================================" <<std::endl;
  std::cout<<"==========================================" <<std::endl;
}
