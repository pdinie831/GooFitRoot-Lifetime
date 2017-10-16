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
#include "ExpPdf.hh" 
#include "ProdPdf.hh" 
#include "RGaussianPdf.hh"
#include "BifurGaussPdf.hh"
#include "ExpGausMPdf.hh" 
#include "ExpGausWithIntPdf.hh"
#include "ExpGausPEEfixSigmaPdf.hh" 
#include "ExpGausProdBPdf.hh"
#include "ExpGausProdEffiBPdf.hh"
#include "ExpGausPEESigmaBPdf.hh" 
#include "PolyEffiPdf.hh" 
#include "ErfcPolyPdf.hh"
#include "SigmoidBpPdf.hh"
#include "SigmoidGausPdf.hh"
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


  gettimeofday(&startTime, NULL);
  startCPU = times(&startProc);

   gROOT ->Reset();
   gROOT->SetStyle("Plain");
   gStyle->SetOptStat(000000);
   gStyle->SetOptFit(000000);


   TCanvas* c1 = new TCanvas("c1","Mass PLOTS",200,10,900,780);
   TCanvas* c2 = new TCanvas("c2","cTau PLOTS",200,10,900,780);
   TCanvas* c3 = new TCanvas("c3","STau PLOTS",200,10,900,780);
    Char_t    InputFileName[300] = "test-ntupla.root";
  Char_t    InputTauBpTreeName[10]   = "TauBpTree";
  TFile*InputFile = TFile::Open(InputFileName,"READ","ROOT file");
  
   Char_t    OutFileName[300] =  "testGoofit3DBp-2016.root";
   gSystem->Exec(Form("rm %s",OutFileName));
   TFile*OutFile = TFile::Open(OutFileName,"RECREATE");
   
   float PlotLineWidth = 1.2;
   float MarkerSize    = 0.35;
 
  double xBpMass;
//  double xBpTau;
  double xBpcTau;
  double xSBpcTau;
//  double c_const       = 0.0299792458;
//     double XMinSign = 5.12;
//     double XMaxSign = 5.44;
  double XMinSign = 5.15;
  double XMaxSign = 5.40;
  double BpMass   = 5.279;
//double BpSigma  = 0.020;
  double BpSigma  = 0.022;

  double NSigmaSB = 3.6;
  double BiasSB   = 6;
 
//      double XMinSign = 5.12;
//      double XMaxSign = 5.44;
//double XMinSign = 5.15;
//double XMaxSign = 5.40;
//   double XMinSBL = 4.879000;
//   double XMaxSBL = 5.159000;
//   double XMinSBR = 5.399000;
//   double XMaxSBR = 5.679000;
  double XMinSBL = BpMass -(BiasSB+NSigmaSB)*BpSigma;
  double XMaxSBL = BpMass - BiasSB          *BpSigma;
  double XMinSBR = BpMass + BiasSB          *BpSigma;
  double XMaxSBR = BpMass +(BiasSB+NSigmaSB)*BpSigma;
  double XMin = 0.006;
  double XMax = 0.34;
  double SXMin = 0.0003;
  double SXMax = 0.005;
//  double SXMax = 0.0079;
//  double SXMax = 0.0079;
  double XStepSign = 0.001;
  double XStepcTau = 0.001;
  double XStepScTau = 0.00002;
  double XStepMinuit = 0.00001;
  double XHScale = 4;
  
  double c_const       = 0.0299792458;

  printf("(xBpMass>%8f && xBpMass<%8f || xBpMass>%8f && xBpMass<%8f)",XMinSBL,XMaxSBL,XMinSBR,XMaxSBR);

  Variable* xMass  = new Variable("xMass",XMinSign, XMaxSign); 
  xMass->numbins = (XMaxSign -XMinSign)/XStepSign;
  TH1F HxMass( "HxMass" , "B^{+} Mass"    ,	     xMass->numbins, xMass->lowerlimit, xMass->upperlimit);
  TH1F pdfHist("pdfHist", "B^{+} Mass Fit",  XHScale*xMass->numbins, xMass->lowerlimit, xMass->upperlimit);
  TH1F sigHist("sigHist", "B^{+} Mass Fit",  XHScale*xMass->numbins, xMass->lowerlimit, xMass->upperlimit);
  TH1F bkgHist("bgkHist", "B^{+} Mass Fit",  XHScale*xMass->numbins, xMass->lowerlimit, xMass->upperlimit);

  Variable* xcTau  = new Variable("xcTau",XMin, XMax); 
  xcTau->numbins   = (XMax -XMin)/XStepcTau;
  Variable* xScTau = new Variable("xScTau",SXMin, SXMax); 
  xScTau->numbins = (SXMax -SXMin)/(XStepScTau);
  std::cout<<"xMass ->numbins = "<<xMass ->numbins<<std::endl;
  std::cout<<"xcTau ->numbins = "<<xcTau ->numbins<<std::endl;
  std::cout<<"xScTau->numbins = "<<xScTau->numbins<<std::endl;
  std::cout<<"xScTau->numbins = "<<xScTau->numbins<<std::endl;
 
  TH1F HxcTau(    "HxcTau"   , "B^{+} cTau",          xcTau  ->numbins, xcTau ->lowerlimit,  xcTau ->upperlimit);
  TH1F HxcTauSB(  "HxcTauSB" , "B^{+} cTau SB",       xcTau  ->numbins, xcTau ->lowerlimit,  xcTau ->upperlimit);
  TH1F HxScTau(   "HxScTau"  , "B^{+} cTau Sigma",    xScTau  ->numbins, xScTau->lowerlimit, xScTau ->upperlimit);
  TH1F HxScTauSB( "HxScTauSB", "B^{+} cTau Sigma SB", xScTau  ->numbins, xScTau->lowerlimit, xScTau ->upperlimit);
  TH2F pdf_cTauSTau_Hist2D( "pdf_cTauSTau_Hist2D" , "Bc+ cTau model   pdf",  XHScale*xcTau ->numbins, xcTau ->lowerlimit, xcTau ->upperlimit, XHScale*xScTau->numbins, xScTau ->lowerlimit, xScTau ->upperlimit);
  TH2F sig_cTauSTau_Hist2D( "sig_cTauSTau_Hist2D" , "Bc+ cTau model   sig",  XHScale*xcTau ->numbins, xcTau ->lowerlimit, xcTau ->upperlimit, XHScale*xScTau->numbins, xScTau ->lowerlimit, xScTau ->upperlimit);
  TH2F bkg_cTauSTau_Hist2D( "bkg_cTauSTau_Hist2D" , "Bc+ cTau model   bkg",  XHScale*xcTau ->numbins, xcTau ->lowerlimit, xcTau ->upperlimit, XHScale*xScTau->numbins, xScTau ->lowerlimit, xScTau ->upperlimit);
//
// Mass Spectrum
//
  Variable* mean   = new Variable("mean"  ,5.2785,XStepMinuit, 5., 5.5);
  Variable* sigma1 = new Variable("sigma1",0.0139,XStepMinuit, 0., 1.);
  Variable* sigma2 = new Variable("sigma2",0.0228,XStepMinuit, 0., 1.);
  Variable* sigma3 = new Variable("sigma3",0.0601,XStepMinuit, 0., 1.);

//   GaussianPdf* gauss1 = new GaussianPdf("gauss1", xMass, mean, sigma1);
//   GaussianPdf* gauss2 = new GaussianPdf("gauss2", xMass, mean, sigma2);
//   GaussianPdf* gauss3 = new GaussianPdf("gauss3", xMass, mean, sigma3);
  RGaussianPdf* gauss1 = new RGaussianPdf("gauss1", xMass, mean, sigma1);
  RGaussianPdf* gauss2 = new RGaussianPdf("gauss2", xMass, mean, sigma2);
  RGaussianPdf* gauss3 = new RGaussianPdf("gauss3", xMass, mean, sigma3);

/*   Variable* meanBckgBp   = new Variable("meanBckgBp" ,5.360,0.00001, 5., 5.5);
  Variable* sigmaBckgBp  = new Variable("sigmaBckgBp",0.030,0.00001, 0. , 1. );
  Variable* meanBckgB0   = new Variable("meanBckgB0" ,5.090,0.00001, 5. , 5.2);
  Variable* sigmaBckgB0  = new Variable("sigmaBckgB0",0.025,0.00001, 0. , 1.);
 */  
//   Variable* meanBckgBp   = new Variable("meanBckgBp" ,5.37 ,0,  5.2, 5.7);
//   Variable* sigmaBckgBp  = new Variable("sigmaBckgBp",0.033,0,  0.013 , 0.05 );
//   Variable* meanBckgB0   = new Variable("meanBckgB0" ,5.090,0,  5.0 , 5.15);
//   Variable* sigmaBckgB0  = new Variable("sigmaBckgB0",0.025,0,  0.01 , 0.05);
  Variable* meanBckgBp   = new Variable("meanBckgBp" ,5.355,0, 5.2, 5.7);
  Variable* sigmaBckgBp  = new Variable("sigmaBckgBp",0.020,0, 0.013 , 0.05 );
  Variable* meanBckgB0   = new Variable("meanBckgB0" ,5.010,0, 5.0 , 5.15);
  Variable* sigmaBckgB0  = new Variable("sigmaBckgB0",0.029,0, 0.02 , 0.04);
  Variable* wb1 = new Variable("wb1",1.75628e-01, 0., 1.);
  Variable* wb2 = new Variable("wb2",1.17246e-01, 0., 1.);
  Variable* wb3 = new Variable("wb3",0.1, 0., 1.);

  RGaussianPdf* gaussBckgBp = new RGaussianPdf("gaussBckgBp", xMass, meanBckgBp, sigmaBckgBp);
  RGaussianPdf* gaussBckgB0 = new RGaussianPdf("gaussBckgB0", xMass, meanBckgB0, sigmaBckgB0);
  

  Variable* wg1 = new Variable("wg1",0.44, 0., 1.);
  Variable* wg2 = new Variable("wg2",0.5 , 0., 1.);
  Variable* signalYield = new Variable("signalYield",8.88646e+05,      100000., 2000000.);
  Variable* bckgYield   = new Variable("bckgYield"  ,1.25481e+05 ,      10000., 400000.);

//   Variable* constaCoef = new Variable("constaCoef", 70, 0.001, 10, 100); 
//   Variable* linearCoef = new Variable("linearCoef", 0.1, 0.001, -0.35, 10.); 
//   Variable* secondCoef = new Variable("secondCoef", 0.1, 0.001, 0, 10);
//   Variable* thirdCoef  = new Variable("thirdCoef" , 0.1, 0.001, 0, 10);
 
//   Variable* constaCoef = new Variable("constaCoef", 70., XStepMinuit, 20., 1000); 
//   Variable* linearCoef = new Variable("linearCoef", 0.1, XStepMinuit, -3.5, 10.); 
//   Variable* secondCoef = new Variable("secondCoef", 0.1, XStepMinuit, 0, 10);
//   Variable* thirdCoef  = new Variable("thirdCoef" , 0.1, XStepMinuit, 0, 10);
  Variable* constaCoef = new Variable("constaCoef", 1. ,XStepMinuit,0.,1000. ); 
  Variable* linearCoef = new Variable("linearCoef", 0.001,XStepMinuit,0,10 ); 

  Variable* p0   = new Variable("p0", -6.86109e-01,-10.,10. ); 
  Variable* p1   = new Variable("p1", 0,0,-10,10 ); 
  Variable* VMinSign = new Variable("VMinSign",XMinSign ); 
  Variable* VMaxSign = new Variable("VMaxSign",XMaxSign ); 
  SimpleCheby2Pdf* SimpleCheby2  = new SimpleCheby2Pdf("SimpleCheby2", xMass, p0, p1,VMinSign,VMaxSign);

// double  fullRange = XMaxSign - XMinSign;
// double  minScaled = -1. + 2. * (XMinSign - xminfull) / fullRange;
// 
// double  maxScaled = +1. - 2. * (xmaxfull - XMaxSign)) / fullRange; 

//  Variable* aslope     = new Variable("slope", -1.);
  //Variable* aslope     = new Variable("slope", 0.39, -10, 10);
  //Variable* apower     = new Variable("apower", 6, 0, 10);
//  Variable* apower     = new Variable("apower", 1.18, XStepMinuit, 0.9, 15.);
//  Variable* apower     = new Variable("apower", 1.18, XStepMinuit, 0.9, 6.);
//  Variable* apower     = new Variable("apower", 1.18, 0.001, 0.9, 5.);
  //Variable* treshold   = new Variable("treshold" ,5.168,XStepMinuit, 5.02, 6.);
//  Variable* treshold   = new Variable("treshold" ,5.33,0, 5.04, 6.);

 
 
  std::vector<Variable*> weightsSignalMass;
  weightsSignalMass.push_back(wg1);
//  weightsSignalMass.push_back(wg2); 

  std::vector<PdfBase*> compsSignalMass;
  compsSignalMass.push_back(gauss1);
  compsSignalMass.push_back(gauss2);
//  compsSignalMass.push_back(gauss3);

  AddPdf signalMass("signalMass", weightsSignalMass, compsSignalMass); 
//  signalMass.addSpecialMask(PdfBase::ForceCommonNorm) ;
  

  vector<Variable*> weightsPoly;
  weightsPoly.push_back(constaCoef);
//  weightsPoly.push_back(linearCoef);
//  weightsPoly.push_back(secondCoef);
//  weightsPoly.push_back(thirdCoef);

  
//  PolynomialPdf* polyTmp = new PolynomialPdf("polyTmp", xMass, weightsPoly); 
//   std::vector<PdfBase*> compsPoly2;
//   compsSignalMass.push_back(polyTmp);
//   compsSignalMass.push_back(polyTmp);
//   
//   ProdPdf* poly      = new ProdPdf("poly"  ,compsPoly2 );
 
  PolynomialPdf* poly = new PolynomialPdf("poly", xMass, weightsPoly); 
  
  std::vector<Variable*> weightsBckgMass;
  weightsBckgMass.push_back(wb1);
  weightsBckgMass.push_back(wb2);
//weightsBckgMass.push_back(wb3);

//  ArgusPdf* argus = new  ArgusPdf("argus", xMass, treshold, aslope, true, apower);  

  std::vector<PdfBase*> compsBckgMass;
  compsBckgMass.push_back(gaussBckgBp);
  compsBckgMass.push_back(gaussBckgB0);
//    compsBckgMass.push_back(argus);
//  compsBckgMass.push_back(poly);
  compsBckgMass.push_back(SimpleCheby2);
 
  AddPdf bckgMass("bckgMass", weightsBckgMass, compsBckgMass);
//  bckgMass.addSpecialMask(PdfBase::ForceCommonNorm) ;

  
//==============================================================================
//==============================================================================
//==============================================================================
// Lifetime
//==============================================================================
//==============================================================================
//==============================================================================

     Variable* cTau     = new Variable("cTau"  ,1./(1.638 *c_const), 10., 40.);
//     Variable* tauSB1   = new Variable("tauSB1",1./(1.620 *c_const), 0., 1000.);
//     Variable* tauSB2   = new Variable("tauSB2",1./(0.400 *c_const), 0., 1000.);

     Variable* tauSB1   = new Variable("tauSB1",1/3.11894e-03, 0., 2000.);
     Variable* tauSB2   = new Variable("tauSB2",1/4.69903e-02, 0., 1000.);
     Variable* tauSB3   = new Variable("tauSB3",1/5.42673e-03, 0., 1000.);

//Variable* cTau     = new Variable("cTau"  ,1./( 1.638 *c_const),0.00, 1000.);
//Variable* tauSB1   = new Variable("tauSB1",1./( 1.440 *c_const),0.00, 1000.);
//Variable* tauSB2   = new Variable("tauSB2",1./( 1.600 *c_const),0.00, 1000.);

  Variable* meanResSign    = new Variable("meanResSign",0);

  Variable* meanResBckg1 = new Variable("meanResBckg1" ,0.);
  Variable* meanResBckg2 = new Variable("meanResBckg2" ,0.);
  Variable* meanResBckg3 = new Variable("meanResBckg3" ,0.);

  Variable* sigmaRes     = new Variable("sigmaRes"      ,0.0003,0,001);

  Variable* sigmaResBckg = new Variable("sigmaResBckg"  ,0.0003,0,001);



//  Variable* meanLandauErrSign      = new Variable( "meanLandauErrSign"        ,0.0015 ,SXMin, SXMax);
//  Variable* sigmaLandauErrorSign   = new Variable( "sigmaLandauErrorSign"     ,0.0002 ,SXMin, SXMax);

//  Variable* meanLandauErrBckg      = new Variable( "meanLandauErrBckg"        ,0.0015 ,SXMin, SXMax);
//  Variable* sigmaLandauErrorBckg   = new Variable( "sigmaLandauErrorBckg"     ,0.0002 ,SXMin, SXMax);
 
//  Variable* meanGaussianErrSign    = new Variable( "meanGaussianErrSign"      ,0.0013  ,SXMin, SXMax);
//  Variable* sigmaGaussianErrorSign = new Variable( "sigmaGaussianErrorSign"   ,0.0003  ,0.00001, SXMax);

//  Variable* meanGaussianErrBckg1   = new Variable( "meanGaussianErrBckg1"      ,0.0013  ,0.,SXMin, SXMax);
//  Variable* sigmaGaussianErrorBckg1= new Variable( "sigmaGaussianErrorBckg1"   ,0.0003  ,0.,0.00001, SXMax);

//  Variable* meanBifurGErrSign      = new Variable( "meanBifurGErrSign"        ,0.0015 ,0., SXMax);
//  Variable* sigmaLBifurGErrSign    = new Variable( "sigmaLBifurGErrSign"      ,0.0003 ,0.00001, SXMax);
//  Variable* sigmaRBifurGErrSign    = new Variable( "sigmaRBifurGErrSign"      ,0.0009 ,0.00001, SXMax);

//  Variable* meanBifurGErrBckg      = new Variable( "meanBifurGErrBckg"	      ,0.0015 ,0., SXMax);
//  Variable* sigmaLBifurGErrBckg    = new Variable( "sigmaLBifurGErrBckg"      ,0.0003 ,0.00001, SXMax);
//  Variable* sigmaRBifurGErrBckg    = new Variable( "sigmaRBifurGErrBckg"      ,0.0009 ,0.00001, SXMax);
  
//   Variable* tauErrSign   = new Variable("tauErrSign",2100,0.,0., 10000.);
//   Variable* tauErrBckg1  = new Variable("tauErrBckg1",2100,0.,0., 10000.);

  Variable* meanGaussianErrSign1    = new Variable( "meanGaussianErrSign1"     ,1.72100e-03  ,0.0001, SXMax);
  Variable* meanGaussianErrSign2    = new Variable( "meanGaussianErrSign2"     ,1.72100e-03  ,0.0001, SXMax);
  Variable* meanGaussianErrSign3    = new Variable( "meanGaussianErrSign3"     ,1.72100e-03  ,0.0001, SXMax);

  Variable* sigmaGaussianErrorSign1 = new Variable( "sigmaGaussianErrorSign1"  ,4.55213e-04  ,0.0001, 0.001);
  Variable* sigmaGaussianErrorSign2 = new Variable( "sigmaGaussianErrorSign2"  ,4.55213e-04  ,0.0001, 0.001);
  Variable* sigmaGaussianErrorSign3 = new Variable( "sigmaGaussianErrorSign3"  ,4.55213e-04  ,0.0001, 0.001);

  Variable* meanGaussianErrBckg1    = new Variable( "meanGaussianErrBckg1"      ,2.00197e-03 ,0.0001, SXMax);
  Variable* meanGaussianErrBckg2    = new Variable( "meanGaussianErrBckg2"      ,1.51338e-03 ,0.00005, 0.005);
  Variable* meanGaussianErrBckg3    = new Variable( "meanGaussianErrBckg3"      ,1.8-03      ,0.0001, SXMax);

  Variable* sigmaGaussianErrorBckg1 = new Variable( "sigmaGaussianErrorBckg1"   ,4.28955e-04,0.0001, 0.001);
  Variable* sigmaGaussianErrorBckg2 = new Variable( "sigmaGaussianErrorBckg2"   ,3.67023e-04,0.0001, 0.001);
  Variable* sigmaGaussianErrorBckg3 = new Variable( "sigmaGaussianErrorBckg3"   , 2.8e-04   ,0.0001, 0.001);


  Variable* tauErrSign1  = new Variable("tauErrSign1" ,1.49315e+03,1000., 20000.);
  Variable* tauErrSign2  = new Variable("tauErrSign2" ,1.49315e+03,1000., 11000.);

  Variable* tauErrBckg1  = new Variable("tauErrBckg1",1/3.73273e-04, 2000., 10000);
  Variable* tauErrBckg2  = new Variable("tauErrBckg2",1/3.72029e-04, 2000., 18000.);
  Variable* tauErrBckg3  = new Variable("tauErrBckg3",1/3.71851e-04, 2000., 20000.);

  Variable* wt1 = new Variable("wt1",0.03, 0., 1.);
  Variable* wt2 = new Variable("wt2",0.05, 0., 1.);
//GenTrue break 5 sigma 0.006-0.35
         double ef0Param = 1.66895e-01       ;
         double ef1Param = 3.02444e+02       ;
         double ef2Param =-1.34418e-01       ;
         double ef3Param =-1.88642e-01       ;
         double ef4Param = 2.93237e+00       ;
         double ef5Param =-2.28278e+01       ;
         double ef6Param = 6.77957e+01       ;
         double ef7Param =-6.83260e+01       ;


//  Variable* ef0 = new Variable("ef0",  5.06489e-02); 
//  Variable* uno = new Variable("uno",  1); 
  Variable* ef0 = new Variable("ef0",  ef0Param); 
  Variable* ef1 = new Variable("ef1",  ef1Param); 
  Variable* ef2 = new Variable("ef2",  ef2Param); 
  Variable* ef3 = new Variable("ef3",  ef3Param); 
  Variable* ef4 = new Variable("ef4",  ef4Param); 
  Variable* ef5 = new Variable("ef5",  ef5Param); 
  Variable* ef6 = new Variable("ef6",  ef6Param); 
  Variable* ef7 = new Variable("ef7",  ef7Param); 

//  coeffEffi.push_back(uno);
//   vector<Variable*> coeffEffi;
//   coeffEffi.push_back(ef0);
//   coeffEffi.push_back(ef1);
//   coeffEffi.push_back(ef2);
//  coeffEffi.push_back(ef3);

 
//   PolyEffiPdf* Effi = new PolyEffiPdf("Effi", xcTau, coeffEffi);
// PolynomialPdf* Effi = new PolynomialPdf("Effi", xcTau, coeffEffi);

//  ErfcPolyPdf  *Effi = new ErfcPolyPdf("Effi", xcTau, ef0,ef1,ef2,ef3,ef4);
      SigmoidBpPdf  *Effi = new SigmoidBpPdf("Effi", xcTau, ef0,ef1,ef2,ef3,ef4,ef5,ef6,ef7);
//      SigmoidGausPdf  *Effi = new SigmoidGausPdf("Effi", xcTau, ef0,ef1,ef2,ef3,ef4,ef5,ef6);


     Variable* XMinV     = new Variable("XMinV"  ,XMin ,0,0,1);
     Variable* XMaxV     = new Variable("XMaxV"  ,XMax ,0,0,1);
     Variable* SXMinV    = new Variable("SXMinV" ,SXMin,0,0,1);
     Variable* SXMaxV    = new Variable("SXMaxV" ,SXMax,0,0,1);
//  Effi->setParameterConstantness(true); 
//  pdfFitBckg1 ->addSpecialMask(PdfBase::ForceSeparateNorm); 
//  Effi        ->addSpecialMask(PdfBase::ForceSeparateNorm); 
   
//    ExpGausProdBPdf* DecayBp	= new ExpGausProdBPdf("DecayBp"    , xcTau, xScTau, meanResSign   , cTau  , sigmaGaussianErrorSign, meanGaussianErrSign,tauErrSign1,
//    XMinV,XMaxV,SXMinV,SXMaxV);
      ExpGausProdBPdf* DecayBp1	= new ExpGausProdBPdf("DecayBp1"    , xcTau, xScTau, meanResSign   , cTau  , sigmaGaussianErrorSign1, meanGaussianErrSign1,tauErrSign2,
      XMinV,XMaxV,SXMinV,SXMaxV);
//  ExpGausProdEffiBPdf* DecayBp1	= new ExpGausProdEffiBPdf("DecayBp1"    , xcTau, xScTau, meanResSign   , cTau  , sigmaGaussianErrorSign, meanGaussianErrSign,tauErrSign1,
//  XMinV,XMaxV,SXMinV,SXMaxV);
//     ExpGausProdBPdf* DecayBp2	= new ExpGausProdBPdf("DecayBp2"    , xcTau, xScTau, meanResSign   , cTau  , sigmaGaussianErrorSign2, meanGaussianErrSign2,tauErrSign2,
//     XMinV,XMaxV,SXMinV,SXMaxV);
      ExpGausProdBPdf* DecayBp2	= new ExpGausProdBPdf("DecayBp2"    , xcTau, xScTau, meanResSign   , cTau  , sigmaGaussianErrorSign2, meanGaussianErrSign2,tauErrSign2,
      XMinV,XMaxV,SXMinV,SXMaxV);
      ExpGausProdBPdf* DecayBp3 = new ExpGausProdBPdf("DecayBp3"    , xcTau, xScTau, meanResSign   , cTau  , sigmaGaussianErrorSign3, meanGaussianErrSign3,tauErrSign2,
      XMinV,XMaxV,SXMinV,SXMaxV);
//  ExpGausProdEffiBPdf* DecayBp2	= new ExpGausProdEffiBPdf("DecayBp2"    , xcTau, xScTau, meanResSign   , cTau  , sigmaGaussianErrorSign2, meanGaussianErrSign2,tauErrSign2,
//  XMinV,XMaxV,SXMinV,SXMaxV);
//  ExpGausProdEffiBPdf* pdfFitBckg1 = new ExpGausProdEffiBPdf("pdfFitBckg1", xcTau, xScTau, meanResBckg, tauSB1, sigmaGaussianErrorBckg1, meanGaussianErrBckg1,tauErrBckg1,
//  XMinV,XMaxV,SXMinV,SXMaxV);
//        ExpGausProdBPdf* pdfFitBckg1 = new ExpGausProdBPdf("pdfFitBckg1", xcTau, xScTau, meanResBckg, tauSB1, sigmaGaussianErrorBckg1, meanGaussianErrBckg1,tauErrBckg1,
//        XMinV,XMaxV,SXMinV,SXMaxV);
       ExpGausProdBPdf* pdfFitBckg1 = new ExpGausProdBPdf("pdfFitBckg1", xcTau, xScTau, meanResBckg2, tauSB1, sigmaGaussianErrorBckg2, meanGaussianErrBckg2,tauErrBckg2,
       XMinV,XMaxV,SXMinV,SXMaxV);
       ExpGausProdBPdf* pdfFitBckg2 = new ExpGausProdBPdf("pdfFitBckg2", xcTau, xScTau, meanResBckg2, tauSB2, sigmaGaussianErrorBckg2, meanGaussianErrBckg2,tauErrBckg2,
       XMinV,XMaxV,SXMinV,SXMaxV);
       ExpGausProdBPdf* pdfFitBckg3 = new ExpGausProdBPdf("pdfFitBckg3", xcTau, xScTau, meanResBckg2, tauSB3, sigmaGaussianErrorBckg2, meanGaussianErrBckg2,tauErrBckg3,
       XMinV,XMaxV,SXMinV,SXMaxV);
//        ExpGausProdBPdf* pdfFitBckg3 = new ExpGausProdBPdf("pdfFitBckg3", xcTau, xScTau, meanResBckg3, tauSB3, sigmaGaussianErrorBckg3, meanGaussianErrBckg3,tauErrBckg3,
//        XMinV,XMaxV,SXMinV,SXMaxV);
 // ExpGausProdEffiBPdf* pdfFitBckg2 = new ExpGausProdEffiBPdf("pdfFitBckg2", xcTau, xScTau, meanResBckg, tauSB2, sigmaGaussianErrorBckg1, meanGaussianErrBckg1,tauErrBckg1,
 // XMinV,XMaxV,SXMinV,SXMaxV);
//       ExpGausProdBPdf* pdfFitBckg2 = new ExpGausProdBPdf("pdfFitBckg2", xcTau, xScTau, meanResBckg, tauSB2, sigmaGaussianErrorBckg2, meanGaussianErrBckg2,tauErrBckg2,
//       XMinV,XMaxV,SXMinV,SXMaxV);
//  ExpGausProdBPdf* pdfFitBckg2 = new ExpGausProdBPdf("pdfFitBckg2", xcTau, xScTau, meanResBckg, tauSB2, sigmaGaussianErrorBckg1, meanGaussianErrBckg1,tauErrBckg2,
//  XMinV,XMaxV,SXMinV,SXMaxV);

// ExpPdf* DecayBp     = new ExpPdf("DecayBp"    ,  xcTau, cTau  );
// ExpPdf* pdfFitBckg1 = new ExpPdf("pdfFitBckg1", xcTau, tauSB1);
// ExpPdf* pdfFitBckg2 = new ExpPdf("pdfFitBckg2", xcTau, tauSB2);


//  ExpGausPEEPdf* DecayBp     = new ExpGausPEEPdf("DecayBp"    , xcTau, xScTau, meanResSign   , cTau  );
//  ExpGausPEEPdf* pdfFitBckg1 = new ExpGausPEEPdf("pdfFitBckg1", xcTau, xScTau, meanResBckg, tauSB1);
//  ExpGausPEEPdf* pdfFitBckg2 = new ExpGausPEEPdf("pdfFitBckg2", xcTau, xScTau, meanResBckg, tauSB2);

  
//  ExpGausPdf* DecayBp	    = new ExpGausPdf("DecayBp"	, xcTau, meanResSign   , sigmaRes, cTau  );
//  ExpGausPdf* pdfFitBckg1 = new ExpGausPdf("pdfFitBckg1", xcTau, meanResBckg, sigmaResBckg, tauSB1);
//  ExpGausPdf* pdfFitBckg2 = new ExpGausPdf("pdfFitBckg2", xcTau, meanResBckg, sigmaResBckg, tauSB2);

//    std::vector<PdfBase*> compspdfFitBp1;
//    compspdfFitBp1.push_back(DecayBp1);
//    compspdfFitBp1.push_back(Effi);
//    ProdPdf* pdfFitBp1	 = new ProdPdf("pdfFitBp1"  , compspdfFitBp1);
// 
//    std::vector<PdfBase*> compspdfFitBp2;
//    compspdfFitBp2.push_back(DecayBp2);
//    compspdfFitBp2.push_back(Effi);
//    ProdPdf* pdfFitBp2	 = new ProdPdf("pdfFitBp2"  , compspdfFitBp2);
// 
//   std::vector<PdfBase*> compspdfDecayBpAdd;
//   compspdfDecayBpAdd.push_back(pdfFitBp1);
//   compspdfDecayBpAdd.push_back(pdfFitBp2);
//   AddPdf *pdfFitBp = new AddPdf("pdfFitBp", weightsSignalMass, compspdfDecayBpAdd);
// 
   std::vector<Variable*> weightsSignalTau;
   weightsSignalTau.push_back(wt1);
//   weightsSignalTau.push_back(wt2); 
   std::vector<PdfBase*> compspdfDecayBpAdd;
   compspdfDecayBpAdd.push_back(DecayBp1);
   compspdfDecayBpAdd.push_back(DecayBp2);
//   compspdfDecayBpAdd.push_back(DecayBp3);
   AddPdf DecayBp("DecayBp", weightsSignalMass, compspdfDecayBpAdd);
//   AddPdf DecayBp("DecayBp", weightsSignalTau, compspdfDecayBpAdd);
//  AddPdf* pdfFitBp =  new AddPdf("pdfFitBp", weightsSignalMass, compspdfDecayBpAdd);
  
 
//  DecayBp.addSpecialMask(PdfBase::ForceSeparateNorm); 
 
 
 
    std::vector<PdfBase*> compspdfFitBp;
    compspdfFitBp.push_back(&DecayBp);
//    compspdfFitBp.push_back(DecayBp1);
    compspdfFitBp.push_back(Effi);
    ProdPdf* pdfFitBp	 = new ProdPdf("pdfFitBp"  , compspdfFitBp);

  Variable* b1 = new Variable("b1",5.38694e-02,0., 1.);
  Variable* b2 = new Variable("b2",8.97145e-01,0., 1.);
  std::vector<Variable*> weightspdfFitBckg;
//  weightspdfFitBckg.push_back(b1);
  weightspdfFitBckg.push_back(b2);

  
  std::vector<PdfBase*> compspdfFitBckgAdd;
//  compspdfFitBckgAdd.push_back(pdfFitBckg1);
  compspdfFitBckgAdd.push_back(pdfFitBckg2);
  compspdfFitBckgAdd.push_back(pdfFitBckg3);
//  AddPdf pdfFitBckgAdd("pdfFitBckgAdd", weightspdfFitBckg, compspdfFitBckgAdd);
//  AddPdf pdfFitBckgTmp("pdfFitBckgTmp",weightsBckgMass , compspdfFitBckgAdd);
  AddPdf pdfFitBckgTmp("pdfFitBckgTmp", weightspdfFitBckg, compspdfFitBckgAdd);
//      AddPdf pdfFitBckg("pdfFitBckg", weightspdfFitBckg, compspdfFitBckgAdd);
  
   std::vector<PdfBase*> compspdfFitBpBckg;
   compspdfFitBpBckg.push_back(&pdfFitBckgTmp);
//   compspdfFitBpBckg.push_back(pdfFitBckg1);
   compspdfFitBpBckg.push_back(Effi);
   ProdPdf pdfFitBckg("pdfFitBckg"  , compspdfFitBpBckg);
 

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
//
// 2DFit
//  
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
  

  std::vector<PdfBase*> compsSignalLife;
  compsSignalLife.push_back(&signalMass);
  compsSignalLife.push_back(pdfFitBp);
//  compsSignalLife.push_back( GaussianErrorSign);
  //compsSignalLife.push_back(ExpGauSign);
//  compsSignalLife.push_back( LandauErrorSign);
  ///compsSignalLife.push_back(BifurGErrorSign);
//  compsSignalLife.push_back(DecayBp);

  std::vector<PdfBase*> compsBckgLife;
  compsBckgLife.push_back(&bckgMass);
  compsBckgLife.push_back(&pdfFitBckg);
//  compsBckgLife.push_back(pdfFitBckg1);
//  compsSignalLife.push_back( GaussianErrorBckg);
  //compsSignalLife.push_back(ExpGauBckg);
//  compsSignalLife.push_back( LandauErrorBckg);
  //compsSignalLife.push_back(BifurGErrorBckg);
//compsBckgLife.push_back(&pdfFitBckgAdd);

  ProdPdf* signalLife = new ProdPdf("signalLife", compsSignalLife);
  ProdPdf* bckgLife   = new ProdPdf("bckgLife  ", compsBckgLife);

  std::vector<Variable*> weightsYield;
  weightsYield.push_back(signalYield);
  weightsYield.push_back(bckgYield);
  

  std::vector<PdfBase*> compsModel;
  
  
  compsModel.push_back(signalLife);
//  compsModel.push_back(gaussBckgBp);
  compsModel.push_back(bckgLife);
  
  
//  compsModel.push_back(poly);
//  compsModel.push_back(argus);
  AddPdf model("model", weightsYield, compsModel); 

//  model.addSpecialMask(PdfBase::ForceCommonNorm) ;

//
// These are used for Plots....
//
  std::vector<PdfBase*> compsMass;
  compsMass.push_back(&signalMass);
  compsMass.push_back(&bckgMass);
  AddPdf modelMass("modelMass", weightsYield, compsMass); 
  
  std::vector<PdfBase*> compscTau;
  compscTau.push_back(pdfFitBp);
  compscTau.push_back(&pdfFitBckg);
  AddPdf model_cTau("model_cTau", weightsYield, compscTau); 
  
//  std::vector<PdfBase*> compsSTau;
//  compsSTau.push_back(LandauErrorSign);
//  compsSTau.push_back(LandauErrorBckg);
//  compsSTau.push_back(BifurGErrorSign);
//  compsSTau.push_back(BifurGErrorBckg);
//  compsSTau.push_back(GaussianErrorSign);
//  compsSTau.push_back(GaussianErrorBckg);
//  compsSTau.push_back(ExpGauSign);
//  compsSTau.push_back(ExpGauBckg);
//  AddPdf model_STau("model_cTau", weightsYield, compsSTau); 
  
  

//
// Data
//
  vector<Variable*> dataVec;
  
  dataVec.push_back(xMass);
  dataVec.push_back(xcTau);
  dataVec.push_back(xScTau);
  UnbinnedDataSet* dataLife = new UnbinnedDataSet(dataVec);
//
  if (!InputFile)
   {
     cout<<"File:"<<InputFileName<<" not found!!!"<<endl;
    exit(1);
   }
   InputFile->ls();
   
   TTree *TauBpTree    = (TTree*)InputFile->Get(InputTauBpTreeName);
   if(!TauBpTree ){
     cout<<"TTree cTau Data: "<< InputTauBpTreeName <<" not found!!!"<<endl;
     exit(1);
   }else{
     cout<<"TTree cTau Data: "<< InputTauBpTreeName <<" OK FOUND!!!"<<endl;
   }  
    
   TauBpTree->SetBranchAddress("xBpMass",&xBpMass);
//   TauBpTree->SetBranchAddress("xBpTau" ,&xBpTau);
   TauBpTree->SetBranchAddress("xBpcTau",&xBpcTau);
   TauBpTree->SetBranchAddress("xSBpcTau",&xSBpcTau);
   int nentries = (int)TauBpTree->GetEntries();
   
   for (Int_t i=0;i<nentries;i++) { 
    TauBpTree->GetEntry(i);
    if(xBpcTau>XMin&&xBpcTau<XMax&&xSBpcTau>SXMin&&xSBpcTau<SXMax){
     if(xBpMass>XMinSign&&xBpMass<XMaxSign){
      xMass->value = xBpMass;
      xcTau->value = xBpcTau;
      xScTau->value = xSBpcTau; 
      dataLife->addEvent();
      HxMass.Fill(xBpMass);
      HxcTau.Fill(xBpcTau);
      HxScTau.Fill(xSBpcTau);
     } 
     if(xBpMass>XMinSBL&&xBpMass<XMaxSBL){
      HxcTauSB.Fill(xBpcTau);
      HxScTauSB.Fill(xSBpcTau);
     } 
     if(xBpMass>XMinSBR&&xBpMass<XMaxSBR){
      HxcTauSB.Fill(xBpcTau);
      HxScTauSB.Fill(xSBpcTau);
     }  
    } 
   }
   cout<<"***********************************"<<endl;
   cout<<"***********************************\n"<<endl;
   cout<<"TauBpTree   Entries = "<<nentries<<endl;
   cout<<"Mass        Entries = "<<HxMass.GetEntries()<<endl;
   cout<<"SideBand    Entries = "<<HxcTauSB.GetEntries()<<endl;
   cout<<"\n***********************************"<<endl;
   cout<<"***********************************"<<endl;
//================================================================================
//================================================================================
///FIT
//================================================================================
//================================================================================
//  double arglis[10];
//  int ierflg= 0;
  model.setData(dataLife);
  FitManager fitter(&model);
  fitter.setMaxCalls(8000);
  cout<<"                  ===*** Start Fit ***=== "<<endl;
  cout<<"                  ===*** Start Fit ***=== "<<endl;
  cout<<"                  ===*** Start Fit ***=== "<<endl;
  fitter.setupMinuit();
  TMinuit * Minuit = fitter.getMinuitObject();  
// 
//       Minuit->FixParameter(7);
//      fitter.runCommand("MIGRAD");
   
//      Minuit->Release(7);
   fitter.runCommand("MIGRAD");
 
//  fitter.runCommand("MINOS");
   fitter.runCommand("HESSE");
//  fitter.fit();   
//  Minuit->SetPrintLevel(1);
//  Minuit->mnmigr();
//  Minuit->mnhess();
//  Minuit->mnmigr();
  fitter.getMinuitValues(); 
  cout<<"		   ===***  End  Fit ***=== "<<endl;
  cout<<"		   ===***  End  Fit ***=== "<<endl;
  cout<<"		   ===***  End  Fit ***=== "<<endl;
  
//================================================================================
//================================================================================
///FIT
//================================================================================
//================================================================================

//================================================================================
///PLOT

//XHScale=10;
// Mass
  UnbinnedDataSet gridMass(xMass);
  double totalDataMass = 0; 
  double NStepMass = XHScale*xMass->numbins;
  for (int i = 0; i < NStepMass; ++i) {
    double step = (xMass->upperlimit - xMass->lowerlimit)/NStepMass;
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


//XHScale=1;
  int NIntegral = 1;
  

//  vector<Variable*> dataPlot;
//  dataPlot.push_back(xcTau);
//  dataPlot.push_back(xScTau);

  vector<Variable*> dataPlot2D;
  dataPlot2D.push_back(xcTau);
  dataPlot2D.push_back(xScTau);

//  vector<Variable*> dataPlotS;
//  dataPlotS.push_back(xScTau);
  
  
//  UnbinnedDataSet grid_cTau(dataPlot);
  UnbinnedDataSet grid_cTau2D(dataPlot2D);
//  UnbinnedDataSet grid_STau(dataPlotS);
  
//  bool first = true;
//  UnbinnedDataSet grid_cTau(xcTau);
//  double totalData_cTau = 0; 
  double NStep  = XHScale*xcTau->numbins;
//  double NSStep   = XHScale*xcTau->numbins;
  double NSStep2D = NIntegral*XHScale*xScTau->numbins;
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
//  vector<vector<double> > pdfVals_cTau;
//  model_cTau.getCompProbsAtDataPoints(pdfVals_cTau); 
//  double totalPdf_cTau = 0; 
//   for (int i = 0; i < grid_cTau.getNumEvents(); ++i) {
//     grid_cTau.loadEvent(i); 
//     pdf_cTau_Hist.Fill(xcTau->value , pdfVals_cTau[0][i]);
//     sig_cTau_Hist.Fill(xcTau->value , pdfVals_cTau[1][i]);
//     bkg_cTau_Hist.Fill(xcTau->value , pdfVals_cTau[2][i]);
//     totalPdf_cTau += pdfVals_cTau[0][i]; 
//   }

//  double pdf_cTau_Integral2D = 0;
//  double sig_cTau_Integral2D = 0;
//  double bkg_cTau_Integral2D = 0;
//model.setData(&grid_cTau2D);
//vector<vector<double> > pdfVals_cTau2D;
//model.getCompProbsAtDataPoints(pdfVals_cTau2D);
     model_cTau.setData(&grid_cTau2D);
     vector<vector<double> > pdfVals_cTau2D;
     model_cTau.getCompProbsAtDataPoints(pdfVals_cTau2D);
  for (int i = 0; i < grid_cTau2D.getNumEvents(); ++i) {
    grid_cTau2D.loadEvent(i); 
    pdf_cTauSTau_Hist2D.Fill(xcTau->value ,xScTau->value , pdfVals_cTau2D[0][i]);
    sig_cTauSTau_Hist2D.Fill(xcTau->value ,xScTau->value , pdfVals_cTau2D[1][i]);
    bkg_cTauSTau_Hist2D.Fill(xcTau->value ,xScTau->value , pdfVals_cTau2D[2][i]);
//     if (i%int(NSStep2D) == 1 && i>0){
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
//   model_STau.setData(&grid_STau);
//   vector<vector<double> > pdfVals_STau;
//   model_STau.getCompProbsAtDataPoints(pdfVals_STau); 
//   double totalPdf_STau = 0;  
//   for (int i = 0; i < grid_STau.getNumEvents(); ++i) {
//     grid_STau.loadEvent(i); 
//     pdf_STau_Hist.Fill(xScTau->value, pdfVals_STau[0][i]);
//     sig_STau_Hist.Fill(xScTau->value, pdfVals_STau[1][i]);
//     bkg_STau_Hist.Fill(xScTau->value, pdfVals_STau[2][i]);
//     totalPdf_STau += pdfVals_STau[0][i]; 
//   }
  
//
// Models plot  
//   pdf_cTau_Hist.Scale((signalYield->value+bckgYield->value)/pdf_cTau_Hist.Integral()*XHScale);
//   sig_cTau_Hist.Scale(signalYield->value/sig_cTau_Hist.Integral()*XHScale);
//   bkg_cTau_Hist.Scale(HxcTauSB.GetEntries()/bkg_cTau_Hist.Integral()*XHScale);
// 
//   pdf_cTau_Hist2D.Scale((signalYield->value+bckgYield->value)/pdf_cTau_Hist2D.Integral()*XHScale);
//   sig_cTau_Hist2D.Scale(signalYield->value/sig_cTau_Hist2D.Integral()*XHScale);
//   bkg_cTau_Hist2D.Scale(HxcTauSB.GetEntries()/bkg_cTau_Hist2D.Integral()*XHScale);
    
  pdf_cTauSTau_X->Scale((signalYield->value+bckgYield->value)/pdf_cTauSTau_X->Integral()*XHScale);
  sig_cTauSTau_X->Scale((signalYield->value)/sig_cTauSTau_X->Integral()*XHScale);
  bkg_cTauSTau_X->Scale((HxcTauSB.GetEntries())/bkg_cTauSTau_X->Integral()*XHScale);
//  bkg_cTauSTau_X->Scale((bckgYield->value)/bkg_cTauSTau_X->Integral()*XHScale);
  pdf_cTauSTau_Y->Scale((signalYield->value+bckgYield->value)/pdf_cTauSTau_Y->Integral()*XHScale);
  sig_cTauSTau_Y->Scale((signalYield->value)/sig_cTauSTau_Y->Integral()*XHScale);
  bkg_cTauSTau_Y->Scale((HxcTauSB.GetEntries())/bkg_cTauSTau_Y->Integral()*XHScale);
//  bkg_cTauSTau_Y->Scale((bckgYield->value)/bkg_cTauSTau_Y->Integral()*XHScale);
  
  

//   sig_STau_Hist.Scale(signalYield->value/sig_STau_Hist.Integral()*XHScale);
//   bkg_STau_Hist.Scale(HxScTauSB.GetEntries()/bkg_STau_Hist.Integral()*XHScale);
//   pdf_STau_Hist.Scale((signalYield->value+bckgYield->value)/pdf_STau_Hist.Integral()*XHScale);
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

  double sigmaw    = sigma1->value*wg1->value+ (1-wg1->value)*sigma2->value;
  double sigmawErr = sqrt(sigma1->error*wg1->value*sigma1->error*wg1->value+ (1-wg1->value)*sigma2->error*(1-wg1->value)*sigma2->error);
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
     leg_sign->AddEntry(&HxMass ,Form( "M_{B^{+}} =   %5.4f  #pm %5.4f",mean->value,mean->error),"");
    }else{
     leg_sign->AddEntry(&HxMass ,Form( "M_{B^{+}} =   %5.4f Fixed",mean->value),"");
     }
    if(sigma1->error!=0){
     leg_sign->AddEntry(&HxMass ,Form( "#sigma_{B^{+}} =   %5.4f  #pm %5.4f",sigmaw,sigmawErr),"");
    }else{
     leg_sign->AddEntry(&HxMass ,Form( "#sigma_{B^{+}} =   %5.4f Fixed",sigmaw),"");
    }
  HxMass.SetMarkerStyle(8);
  HxMass.SetMarkerSize(MarkerSize);
  HxMass.SetTitle("");
  HxMass.Draw("E1"); 
  leg_sign->Draw("same");
//  HxMass.Draw("p"); 
  pdfHist.SetLineColor(kBlue);
  pdfHist.SetLineWidth(PlotLineWidth); 
  pdfHist.Draw("same"); 
  sigHist.SetLineColor(kMagenta);
  sigHist.SetLineStyle(kDashed); 
  sigHist.SetLineWidth(PlotLineWidth); 
  sigHist.Draw("same"); 
  bkgHist.SetLineColor(kRed);
  bkgHist.SetLineStyle(kDashed); 
  bkgHist.SetLineWidth(PlotLineWidth); 
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
  if( b2->error!=0){
      leg_pdfSB->AddEntry(&HxcTau ,Form( "b2   =  %5.3f #pm %5.3f     ",b2->value,b2->error)   ,"");
  }else{      
      leg_pdfSB->AddEntry(&HxcTau ,Form( "b2   =  %5.3f     Fixed     ",b2->value)   ,"");
  }   
  if( tauSB1->error!=0){
      leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB1} =  %5.3f #pm %5.3f     }",1/(tauSB1->value),tauSB1->error/((tauSB1->value)*(tauSB1->value)))   ,"");
  }else{      
      leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB1} =  %5.3f	 Fixed	  }",1/(tauSB1->value))   ,"");
  }   
  if( tauSB2->error!=0){
      leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB2} =  %5.3f #pm %5.3f     }",1/(tauSB2->value),tauSB2->error/((tauSB2->value)*(tauSB2->value)))   ,"");
  }else{      
      leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB2} =  %5.3f	 Fixed	  }",1/(tauSB2->value))   ,"");
  }   
  if( tauSB3->error!=0){
      leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB2} =  %5.3f #pm %5.3f     }",1/(tauSB3->value),tauSB3->error/((tauSB3->value)*(tauSB3->value)))   ,"");
  }else{      
      leg_pdfSB->AddEntry(&HxcTau ,Form( "#color[2]{#scale[1.5]{#tau}_{SB2} =  %5.3f Fixed     }",1/(tauSB3->value))   ,"");
  }   
  leg_pdfSB->AddEntry(&HxcTau ,"#color[4]{#scale[0.8]{- Fit model	     }}"   ,"");
  leg_pdfSB->AddEntry(&HxcTau ,"#color[6]{#scale[0.8]{- Signal model	     }}"   ,"");
  leg_pdfSB->AddEntry(&HxcTau ,"#color[2]{#scale[0.8]{- Background model on SB}}"   ,"");
  HxcTau.SetMarkerStyle(8);
  HxcTau.SetMarkerSize(MarkerSize);
  HxcTau.SetTitle("");
  HxcTau.Draw("E1");
  HxcTauSB.SetMarkerStyle(8);
  HxcTauSB.SetMarkerSize(0.5);
  HxcTauSB.SetMarkerColor(kRed);
  HxcTauSB.Draw("same,E1");
  leg_pdfSB->Draw("same");
/*   pdf_cTau_Hist2D.SetLineColor(kBlue);
  pdf_cTau_Hist2D.SetLineWidth(3); 
  pdf_cTau_Hist2D.Draw("same"); 
  sig_cTau_Hist2D.SetLineColor(kMagenta);
  sig_cTau_Hist2D.SetLineStyle(kDashed); 
  sig_cTau_Hist2D.SetLineWidth(2); 
  sig_cTau_Hist2D.Draw("same"); 
  bkg_cTau_Hist2D.SetLineColor(kRed);
  bkg_cTau_Hist2D.SetLineStyle(kDashed); 
  bkg_cTau_Hist2D.SetLineWidth(2); 
  bkg_cTau_Hist2D.Draw("same"); 
 */  
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
  leg_pdfResolution->SetHeader("B^{+} resolution Fit Projections");
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
  leg_pdfResolution->AddEntry(&HxScTau ,"#color[4]{#scale[0.8]{- Pdf  Resolution model	     }}"   ,"");
  leg_pdfResolution->AddEntry(&HxScTau ,"#color[6]{#scale[0.8]{- Signal  Resolution  model   }}"   ,"");
  leg_pdfResolution->AddEntry(&HxScTau ,"#color[2]{#scale[0.8]{- Background Resolution model (on SB events)}}"   ,"");
  HxScTau.Draw("E1");
  leg_pdfResolution->Draw("same");
  HxScTauSB.SetMarkerStyle(8);
  HxScTauSB.SetMarkerSize(MarkerSize);
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


//   pdf_STau_Hist.SetLineColor(kBlue);
//   pdf_STau_Hist.SetLineWidth(2); 
//   pdf_STau_Hist.Draw("same"); 
//   sig_STau_Hist.SetLineColor(kMagenta);
//   sig_STau_Hist.SetLineStyle(kDashed); 
//   sig_STau_Hist.SetLineWidth(2); 
//   sig_STau_Hist.Draw("same"); 
//   bkg_STau_Hist.SetLineColor(kRed);
//   bkg_STau_Hist.SetLineStyle(kDashed); 
//   bkg_STau_Hist.SetLineWidth(2); 
//   bkg_STau_Hist.Draw("same"); 
  
  HxcTau.Write();
  HxcTauSB.Write();
  HxScTau.Write();
  HxScTauSB.Write();
  pdf_cTauSTau_Hist2D.Write();
  sig_cTauSTau_Hist2D.Write();
  bkg_cTauSTau_Hist2D.Write();
//   pdf_cTau_Hist2D.Write();
//   sig_cTau_Hist2D.Write();
//   bkg_cTau_Hist2D.Write();
//   pdf_cTau_Hist.Write();
//   sig_cTau_Hist.Write();
//   bkg_cTau_Hist.Write();
//   pdf_STau_Hist.Write();
//   sig_STau_Hist.Write();
//   bkg_STau_Hist.Write();
  pdf_cTauSTau_X->Write();
  pdf_cTauSTau_Y->Write();
  sig_cTauSTau_X->Write();
  sig_cTauSTau_Y->Write();
  bkg_cTauSTau_X->Write();
  bkg_cTauSTau_Y->Write();
  c1->Write();
  c2->Write();
  c3->Write();
  char PDFNameMass[50] = "Bp-Mass-2016.pdf";
  char PDFNamecTau[50] = "Bp-cTau-2016.pdf";
  char PDFNameReso[50] = "Bp-Reso-2016.pdf";
  char testo[130] ;
  sprintf(testo,"mv %s %s.tmp",PDFNameMass,PDFNameMass);
  gSystem->Exec(testo);
  sprintf(testo,"mv %s %s.tmp",PDFNamecTau,PDFNamecTau);
  gSystem->Exec(testo);
  sprintf(testo,"mv %s %s.tmp",PDFNameReso,PDFNameReso);
  gSystem->Exec(testo);
  
  std::cout<<"Tau  [ps]= "<<1/(c_const*cTau->value)<<"+/-"<<cTau->error/((c_const*cTau->value)*(cTau->value))<<std::endl;
  std::cout<<"1/TauSB1 = "<<1/(tauSB1->value)<<"+/-"<<tauSB1->error/((tauSB1->value)*(tauSB1->value))<<std::endl;
  std::cout<<"1/TauSB2 = "<<1/(tauSB2->value)<<"+/-"<<tauSB2->error/((tauSB2->value)*(tauSB2->value))<<std::endl;
  std::cout<<"1/TauSB3 = "<<1/(tauSB3->value)<<"+/-"<<tauSB3->error/((tauSB3->value)*(tauSB3->value))<<std::endl;

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
