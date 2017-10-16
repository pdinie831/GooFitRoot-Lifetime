//g++ -o testGoofit testGoofit.cc   `root-config --cflags --libs` -lProof -lRooFit
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

TCanvas* c1; 

   c1 = new TCanvas("c1","GPU PLOTS",200,10,900,780);
//   TCanvas*c2 = new TCanvas("c2","PLOTS",200,10,900,780);
   //TPad* pad1 = (TPad*)c1->GetPad(0);
//   TPad* pad2 = (TPad*)c2->GetPad(0);
   //pad1->SetLeftMargin(0.15); 
//   pad2->SetLeftMargin(0.15); 





  Char_t    InputFileName[300] = "testproof2DCut-Bp-SaraCuts-NewCut.root";
  Char_t    InputTauBpTreeName[10]   = "TauBpTree";
  TFile*InputFile = TFile::Open(InputFileName,"READ","ROOT file");
  
   Char_t    OutFileName[300] =
   "testGoofit.root";
   gSystem->Exec(Form("rm %s",OutFileName));
   TFile*OutFile = TFile::Open(OutFileName,"RECREATE");
  double xBpMass;
  double xBpTau;
  double xBpcTau;
//  double c_const       = 0.0299792458;
    double XMinSign = 5.119;
    double XMaxSign = 5.439;
//double XMinSign = 5.00;
//double XMaxSign = 5.5;
  double XMin = 0.0;
  double XMax = 0.44;
  double XStep = 0.001;
  double XHScale = 100;
  
  Variable* xMass  = new Variable("xMass",XMinSign, XMaxSign); 
  xMass->numbins = (XMaxSign -XMinSign)/XStep;
  Variable* mean   = new Variable("mean"  ,5.278,0.001, 5., 5.5);
  Variable* sigma1 = new Variable("sigma1",0.014,0.001, 0., 1.);
  Variable* sigma2 = new Variable("sigma2",0.030,0.001, 0., 1.);
  Variable* sigma3 = new Variable("sigma3",0.060,0.001, 0., 1.);

  GaussianPdf* gauss1 = new GaussianPdf("gauss1", xMass, mean, sigma1);
  GaussianPdf* gauss2 = new GaussianPdf("gauss2", xMass, mean, sigma2);
  GaussianPdf* gauss3 = new GaussianPdf("gauss3", xMass, mean, sigma3);

/*   Variable* meanBckgBp   = new Variable("meanBckgBp" ,5.360,0.00001, 5., 5.5);
  Variable* sigmaBckgBp  = new Variable("sigmaBckgBp",0.030,0.00001, 0. , 1. );
  Variable* meanBckgB0   = new Variable("meanBckgB0" ,5.090,0.00001, 5. , 5.2);
  Variable* sigmaBckgB0  = new Variable("sigmaBckgB0",0.025,0.00001, 0. , 1.);
 */  
  Variable* meanBckgBp   = new Variable("meanBckgBp" ,5.360,0.001, 5.3, 5.5);
  Variable* sigmaBckgBp  = new Variable("sigmaBckgBp",0.030,0.001, 0.01 , 5. );
  Variable* meanBckgB0   = new Variable("meanBckgB0" ,5.090,0.001, 5.0 , 5.25);
  Variable* sigmaBckgB0  = new Variable("sigmaBckgB0",0.025,0.001, 0. , 3.);
  Variable* wb1 = new Variable("wb1",0.1,0.001, 0., 1.);
  Variable* wb2 = new Variable("wb2",0.1,0.001, 0., 1.);
  Variable* wb3 = new Variable("wb3",0.1,0.001, 0., 1.);

  GaussianPdf* gaussBckgBp = new GaussianPdf("gaussBckgBp", xMass, meanBckgBp, sigmaBckgBp);
  GaussianPdf* gaussBckgB0 = new GaussianPdf("gaussBckgB0", xMass, meanBckgB0, sigmaBckgB0);
  

  Variable* wg1 = new Variable("wg1",0.2,0.001, 0., 1.);
  Variable* wg2 = new Variable("wg2",0.5,0.001, 0., 1.);
  Variable* signalYield = new Variable("signalYield",500000,0.001, 0., 1000000.);
  Variable* bckgYield   = new Variable("bckgYield"  , 30000,0.001, 0., 1000000.);

  Variable* constaCoef = new Variable("constaCoef", 70, 0.001, 10, 100); 
  Variable* linearCoef = new Variable("linearCoef", 0.1, 0.001, -0.35, 10.); 
  Variable* secondCoef = new Variable("secondCoef", 0.1, 0.001, 0, 10);
  Variable* thirdCoef  = new Variable("thirdCoef" , 0.1, 0.001, 0, 10);
 
//  Variable* aslope     = new Variable("slope", -1.);
  Variable* aslope     = new Variable("slope", 0.39, 0.001, -0.3, 1.);
//  Variable* apower     = new Variable("apower", 0.5);
  Variable* apower     = new Variable("apower", 1.18, 0.001, 0.9, 5.);
  Variable* treshold   = new Variable("treshold" ,5.168,0.001, 5.15, 5.5);

  TH1F HxMass( "HxMass" , "B+ Mass"    ,          xMass->numbins, xMass->lowerlimit, xMass->upperlimit);
  TH1F pdfHist("pdfHist", "B+ Mass Fit",  XHScale*xMass->numbins, xMass->lowerlimit, xMass->upperlimit);
  TH1F sigHist("sigHist", "B+ Mass Fit",  XHScale*xMass->numbins, xMass->lowerlimit, xMass->upperlimit);
  TH1F bkgHist("bgkHist", "B+ Mass Fit",  XHScale*xMass->numbins, xMass->lowerlimit, xMass->upperlimit);
 
 
  std::vector<Variable*> weightsGauss;
  weightsGauss.push_back(wg1);
//  weightsGauss.push_back(wg2);

  std::vector<PdfBase*> compsSignal;
  compsSignal.push_back(gauss1);
  compsSignal.push_back(gauss2);
//  compsSignal.push_back(gauss3);

  AddPdf signal("signal", weightsGauss, compsSignal); 
  

  vector<Variable*> weightsPoly;
  weightsPoly.push_back(constaCoef);
  weightsPoly.push_back(linearCoef);
//  weightsPoly.push_back(secondCoef);
//  weightsPoly.push_back(thirdCoef);

  PolynomialPdf* poly = new PolynomialPdf("poly", xMass, weightsPoly); 
  
  std::vector<Variable*> weightsBckg;
  weightsBckg.push_back(wb1);
//  weightsBckg.push_back(wb2);
  weightsBckg.push_back(wb3);

  ArgusPdf* argus = new  ArgusPdf("argus", xMass, treshold, aslope, true, apower);  
  
  std::vector<PdfBase*> compsBckg;
  compsBckg.push_back(gaussBckgBp);
//  compsBckg.push_back(gaussBckgB0);
  compsBckg.push_back(argus);
  compsBckg.push_back(poly);
  
  AddPdf bckgSignTot("bckgSignTot", weightsBckg, compsBckg); 
  
  std::vector<Variable*> weightsYield;
  weightsYield.push_back(signalYield);
  weightsYield.push_back(bckgYield);
  

  std::vector<PdfBase*> compsModel;
  
  
  compsModel.push_back(&signal);
//  compsModel.push_back(gaussBckgBp);
  compsModel.push_back(&bckgSignTot);
//  compsModel.push_back(poly);
//  compsModel.push_back(argus);
  AddPdf model("model", weightsYield, compsModel); 

  UnbinnedDataSet* dataMass = new UnbinnedDataSet(xMass);

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
   TauBpTree->SetBranchAddress("xBpTau" ,&xBpTau);
   TauBpTree->SetBranchAddress("xBpcTau",&xBpcTau);
   int nentries = (int)TauBpTree->GetEntries();
   for (Int_t i=0;i<nentries;i++) {
    TauBpTree->GetEntry(i);
    if(xBpMass>XMinSign&&xBpMass<XMaxSign){
     if(xBpcTau>XMin&&xBpcTau<XMax){
      xMass->value = xBpMass;
      dataMass->addEvent(); 
      HxMass.Fill(xMass->value);
     } 
    } 
   }

///FIT
  model.setData(dataMass);
  FitManager fitter(&model);
  fitter.setMaxCalls(20000);
  fitter.fit(); 
  TMinuit * Minuit = fitter.getMinuitObject();
//  Minuit->SetPrintLevel(1);
//  Minuit->mnmigr();
  Minuit->mnhess();
//  Minuit->mnmigr();
  fitter.getMinuitValues(); 
  
///FIT

///PLOT
  UnbinnedDataSet grid(xMass);
  double totalData = 0; 
  double NStep = XHScale*xMass->numbins;
  for (int i = 0; i < NStep; ++i) {
    double step = (xMass->upperlimit - xMass->lowerlimit)/NStep;
    xMass->value = xMass->lowerlimit + (i + 0.5) * step;
    grid.addEvent(); 
   totalData++; 
  }

  model.setData(&grid);
  vector<vector<double> > pdfVals;
  model.getCompProbsAtDataPoints(pdfVals); 
  double totalPdf = 0; 
  for (int i = 0; i < grid.getNumEvents(); ++i) {
    grid.loadEvent(i); 
    pdfHist.Fill(xMass->value, pdfVals[0][i]);
    sigHist.Fill(xMass->value, pdfVals[1][i]);
    bkgHist.Fill(xMass->value, pdfVals[2][i]);
    totalPdf += pdfVals[0][i]; 
  }
  
  
  pdfHist.Scale((signalYield->value+bckgYield->value)/pdfHist.Integral()*XHScale);
  sigHist.Scale(signalYield->value/sigHist.Integral()*XHScale);
  bkgHist.Scale(bckgYield->value/bkgHist.Integral()*XHScale);
  std::cout<<"Signal Yield = "<< signalYield->value<<std::endl;
  std::cout<<"Bckg   Yield = "<< bckgYield->value<<std::endl;
  std::cout<<"Tot   Yield  = "<< signalYield->value+bckgYield->value<<std::endl;
  
  
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
    if(sigma1->error!=0){
     leg_sign->AddEntry(&HxMass ,Form( "#sigma_{B_{c}^{+}} =   %5.4f  #pm %5.4f",sigma1->value,sigma1->error),"");
    }else{
     leg_sign->AddEntry(&HxMass ,Form( "#sigma_{B_{c}^{+}} =   %5.4f Fixed",sigma1->value),"");
    }
  HxMass.SetMarkerStyle(8);
  HxMass.SetMarkerSize(0.5);
  HxMass.Draw(); 
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
  c1->Write();
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
