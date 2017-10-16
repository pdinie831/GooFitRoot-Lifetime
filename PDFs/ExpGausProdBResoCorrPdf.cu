#include "ExpGausProdBResoCorrPdf.hh"

EXEC_TARGET fptype device_ExpGausProdBResoCorr (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x         = evt[indices[2+indices[0]]]; 
  fptype sigmaTrue = evt[indices[3+indices[0]]]; 
  fptype mean      = p[indices[1]];
  fptype lambda    = p[indices[2]];
  fptype sigmas    = p[indices[3]];
  fptype means     = p[indices[4]];
  fptype lambdas   = p[indices[5]];
  fptype lo        = p[indices[6]];
  fptype hi        = p[indices[7]];
  fptype los       = p[indices[8]];
  fptype his       = p[indices[9]];
  fptype corrP1    = p[indices[10]];
  fptype corrP0    = p[indices[11]];
  
//  fptype sigma = p[indices[2]];
//  fptype alpha = p[indices[3]];
//  sigmaM = sigma;
    fptype sigma = (corrP1*x+corrP0)*sigmaTrue;
    if(sigma<0){
     printf("Error: ExpGausProdBResoCorrPdf sigma<0!!! ==> corrP1=%f  corrP0=%f \n", corrP1,corrP0);
     return 0;
    }
    fptype ret = 0.5*lambda; 
    fptype exparg = ret * (2*mean + lambda*sigma*sigma - 2*x);
    fptype erfarg = (mean + lambda*sigma*sigma - x) / (sigma * 1.4142135623);
    ret *= EXP(exparg); 
    ret *= ERFC(erfarg); 

    fptype y = sigma;
    fptype rets = 0.5*lambdas; 
    fptype expargs = rets * (2*means + lambdas*sigmas*sigmas - 2*y);
    fptype erfargs = (means + lambdas*sigmas*sigmas - y) / (sigmas * 1.4142135623);
//fptype exparg = lambda * (lambda*sigma*sigma/2. + mean-x);
//fptype erfarg = ((mean-x)/(sigma*sigma) + lambda)*sigma /  1.4142135623;
//
    rets *= EXP(expargs); 
    rets *= ERFC(erfargs); 
    ret=rets*ret;
  if (ret<=0){
    printf("Error: ExpGausProdBResoCorrPdf<=0!!! ==> x=%f  sigma=%f mean=%f lambda=%f sigmas=%f means=%f lambdas=%fret=%f\n", x, sigma, mean,lambda , sigmas, means,lambdas ,ret);
    printf("Error: ExpGausProdBResoCorrPdf<=0!!! ==> lo=%f  hi=%f los=%f his=%f \n", lo,hi,los,his);
    printf("Error: ExpGausProdBResoCorrPdf<=0!!! ==> corrP1=%f  corrP0=%f \n", corrP1,corrP0);
     return 0;
  }   
//
//    fptype lo = 0.01;
//    fptype hi = 0.1 ;
//fptype lo = x;
//fptype hi = x+0.00000000001 ;
//
  fptype intg_hi = 0.5*(1-EXP(lambda*(lambda*sigma*sigma/2.+mean - hi))*
                   erfc((lambda+(mean-hi)/(sigma*sigma))*sigma/1.4142135623)+ erf((hi-mean)/(sigma*1.4142135623)));   
  fptype intg_lo = 0.5*(1-EXP(lambda*(lambda*sigma*sigma/2.+mean - lo))*
                   erfc((lambda+(mean-lo)/(sigma*sigma))*sigma/1.4142135623)+ erf((lo-mean)/(sigma*1.4142135623)));  

//  fptype los = 0.0003;
//  fptype his = 0.007 ;
  fptype intg_his = 0.5*(1-EXP(lambdas*(lambdas*sigmas*sigmas/2.+means - his))*
                   erfc((lambdas+(means-his)/(sigmas*sigmas))*sigmas/1.4142135623)+ erf((his-means)/(sigmas*1.4142135623)));   
  fptype intg_los = 0.5*(1-EXP(lambdas*(lambdas*sigmas*sigmas/2.+means - los))*
                   erfc((lambdas+(means-los)/(sigmas*sigmas))*sigmas/1.4142135623)+ erf((los-means)/(sigmas*1.4142135623)));  
		    
//  
//   fptype u_hi = lambda * (hi - mean);
//   fptype v_hi = lambda * sigma;
//   fptype expa_hi = -u_hi+ v_hi*v_hi*0.5+LOG(0.5*( 1+ERF((u_hi-v_hi*v_hi)/(v_hi*1.4142135623))) );
//   fptype intg_hi1 = (0.5*(1+ERF(u_hi/(v_hi*1.4142135623))) - EXP(expa_hi));
// //
//   fptype u_lo = lambda * (lo - mean);
//   fptype v_lo = lambda * sigma;
//   fptype expa_lo = -u_lo+ v_lo*v_lo*0.5+LOG(0.5*( 1+ERF((u_lo-v_lo*v_lo)/(v_lo*1.4142135623))) );
//   fptype intg_lo1 = (0.5*(1+ERF(u_lo/(v_lo*1.4142135623))) - EXP(expa_lo));
//
 
 
   fptype integral = fabs(intg_hi-intg_lo)*fabs(intg_his-intg_los);
//  fptype integral1 = fabs(intg_hi1-intg_lo1);
   ret = ret/integral;
 
//if ((0 == THREADIDX) && (0 == BLOCKIDX)){
//  printf("\n\nExpGausProdBResoCorrPdf x=%f  sigma=%f mean=%f lambda=%f ret=%f integral =%f\n", x, sigma, mean,lambda , ret, integral);
//  printf("ExpGausProdBResoCorrPdf x=%f  sigma=%f mean=%f lambda=%f ret=%f integral1=%f\n", x, sigma, mean,lambda , ret, integral1);
//}  
 
/*  if ((0 == THREADIDX) && (0 == BLOCKIDX)){
  printf("ExpGausProdBResoCorrPdf->integrate  sigma=%f mean=%f lambda=%f Integral=%f\n",    sigma, mean,lambda , integral);
  printf("ExpGausProdBResoCorrPdf->integrate  sigma=%f mean=%f lambda=%f IntegralLO=%f\n",  sigma, mean,lambda , intg_lo);
  printf("ExpGausProdBResoCorrPdf->integrate  sigma=%f mean=%f lambda=%f IntegralHI=%f\n",  sigma, mean,lambda , intg_hi);
  printf("ExpGausProdBResoCorrPdf->integrate  sigma=%f mean=%f lambda=%f EXP(expa_hi)=%f\n",  sigma, mean,lambda , EXP(expa_hi));
  printf("ExpGausProdBResoCorrPdf->integrate  sigma=%f mean=%f lambda=%f EXP(expa_lo)=%f\n",  sigma, mean,lambda , EXP(expa_lo));
 }
 */// printf("ExpGausProdBResoCorrPdf->host_indices  host_indices0=%d host_indices1=%d host_indices2=%d\n",  host_indices[0], host_indices[1],host_indices[2]);
// printf("ExpGausProdBResoCorrPdf->     indices  indices0=%d indices1=%d indices2=%d\n",  indices[0], indices[1],indices[2]);
//  if ((0 == THREADIDX) && (0 == BLOCKIDX)){
//   printf("ExpGausProdBResoCorrPdf x=%f  sigma=%f mean=%f lambda=%f ret=%f integral=%f\n", x, sigma, mean,lambda , ret, integral);
//    printf("ExpGausProdBResoCorrPdf ind0 =%d  evt0 =%f \n",indices[0] ,evt[indices[0]] );
//    printf("ExpGausProdBResoCorrPdf ind1 =%d  evt1 =%f \n",indices[1] ,evt[indices[1]] );
//    printf("ExpGausProdBResoCorrPdf ind2 =%d  evt2 =%f \n",indices[2] ,evt[indices[2]] );
//    printf("ExpGausProdBResoCorrPdf ind3 =%d  evt3 =%f \n",indices[3] ,evt[indices[3]] );
//    printf("ExpGausProdBResoCorrPdf ind4 =%d  evt4 =%f \n",indices[4] ,evt[indices[4]] );
//    printf("ExpGausProdBResoCorrPdf ind5 =%d  evt5 =%f \n",indices[5] ,evt[indices[5]] );
//    printf("ExpGausProdBResoCorrPdf ind6 =%d  evt6 =%f \n",indices[6] ,evt[indices[6]] );
//    printf("ExpGausProdBResoCorrPdf ind7 =%d  evt7 =%f \n",indices[7] ,evt[indices[7]] );
//    printf("ExpGausProdBResoCorrPdf ind8 =%d  evt8 =%f \n",indices[8] ,evt[indices[8]] );
//    printf("ExpGausProdBResoCorrPdf ind9 =%d  evt9 =%f \n",indices[9] ,evt[indices[9]] );
//    printf("ExpGausProdBResoCorrPdf ind10=%d  evt10=%f \n",indices[10],evt[indices[10]] );
//    printf("ExpGausProdBResoCorrPdf ind11=%d  evt11=%f \n",indices[11],evt[indices[11]] );
//    printf("ExpGausProdBResoCorrPdf ind12=%d  evt12=%f \n",indices[12],evt[indices[12]] );
//    printf("ExpGausProdBResoCorrPdf ind13=%d  evt13=%f \n",indices[13],evt[indices[13]] );
//    printf("ExpGausProdBResoCorrPdf ind14=%d  evt14=%f \n",indices[14],evt[indices[14]] );
// }
//  return 0; 
//
  return ret;
}

MEM_DEVICE device_function_ptr ptr_to_ExpGausProdBResoCorr = device_ExpGausProdBResoCorr; 

 __host__ ExpGausProdBResoCorrPdf::ExpGausProdBResoCorrPdf (std::string n, Variable* _x, Variable* _s,  Variable* mean, Variable* lambda, Variable* sigmas, Variable* means, Variable* lambdas, 
                                             Variable* lo, Variable* hi, Variable* los, Variable* his, Variable* corrP1, Variable* corrP0) 
  : GooPdf(_x, n)
{
//   registerObservable(_x); //already registered!!!
   registerObservable(_s);
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(mean));
  pindices.push_back(registerParameter(lambda));
  pindices.push_back(registerParameter(sigmas));
  pindices.push_back(registerParameter(means));
  pindices.push_back(registerParameter(lambdas));
  pindices.push_back(registerParameter(lo));
  pindices.push_back(registerParameter(hi));
  pindices.push_back(registerParameter(los));
  pindices.push_back(registerParameter(his));
  pindices.push_back(registerParameter(corrP1));
  pindices.push_back(registerParameter(corrP0));
  GET_FUNCTION_ADDR(ptr_to_ExpGausProdBResoCorr);
  initialise(pindices); 
}

 __host__ fptype ExpGausProdBResoCorrPdf::integrate (fptype lo, fptype hi) const {
// printf("integratexxx\n");
/*  unsigned int* indices = host_indices+parameters; 
 fptype sigmaM = 0.0017; 
 fptype mean = host_params[indices[1]]  ;
 fptype lambda  = host_params[indices[2]]  ;
 fptype u_hi = lambda * (hi - mean);
 fptype v_hi = lambda * sigmaM;
// fptype v_hi = lambda * u_hi;
 fptype expa_hi = -u_hi+ v_hi*v_hi*0.5+LOG(0.5*( 1+ERF((u_hi-v_hi*v_hi)/(v_hi*1.4142135623))) );
 fptype intg_hi = (0.5*(1+ERF(u_hi/(v_hi*1.4142135623))) - EXP(expa_hi));

 fptype u_lo = lambda * (lo - mean);
 fptype v_lo = lambda * sigmaM;
// fptype v_lo = lambda * u_lo;
 fptype expa_lo = -u_lo+ v_lo*v_lo*0.5+LOG(0.5*( 1+ERF((u_lo-v_lo*v_lo)/(v_lo*1.4142135623))) );
 fptype intg_lo = (0.5*(1+ERF(u_lo/(v_lo*1.4142135623))) - EXP(expa_lo));
 
// printf("ExpGausProdBResoCorrPdf->integrate  sigma=%f mean=%f lambda=%f Integral=%f\n",  sigmaM, mean,lambda , (intg_hi-intg_lo));
// printf("ExpGausProdBResoCorrPdf->host_indices  host_indices0=%d host_indices1=%d host_indices2=%d\n",  host_indices[0], host_indices[1],host_indices[2]);
// printf("ExpGausProdBResoCorrPdf->     indices  indices0=%d indices1=%d indices2=%d\n",  indices[0], indices[1],indices[2]);
 return (intg_hi-intg_lo);
 */
  return 1.;
 }
