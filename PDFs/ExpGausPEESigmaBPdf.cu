#include "ExpGausPEESigmaBPdf.hh"

EXEC_TARGET fptype device_ExpGausPEESigmaB (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x     = evt[indices[2+indices[0]]]; 
//  fptype sigma = evt[indices[3+indices[0]]]; 
  fptype sigma = p[indices[1]]; 
  fptype mean  = p[indices[2]];
  fptype tau   = p[indices[3]];
  fptype lo    = p[indices[4]];
  fptype hi    = p[indices[5]];
  
//  fptype sigma = p[indices[2]];
//  fptype alpha = p[indices[3]];
//  sigmaM = sigma;
  fptype ret = 0.5*tau; 
//  fptype exparg = ret * (2*mean + tau*sigma*sigma - 2*x);
//  fptype erfarg = (mean + tau*sigma*sigma - x) / (sigma * 1.4142135623);
  fptype exparg = tau * (tau*sigma*sigma/2. + mean-x);
  fptype erfarg = ((mean-x)/(sigma*sigma) + tau)*sigma /  1.4142135623;
//
  ret *= EXP(exparg); 
//  ret *= (1- ERF(erfarg)); 
  ret *= ERFC(erfarg); 
  if (ret<=0){
    printf("Error: ExpGausPEESigmaBPdf<=0!!! ==> x=%f  sigma=%f mean=%f tau=%f ret=%f\n", x, sigma, mean,tau , ret);
  }   
//
//    fptype lo = 0.0003;
//    fptype hi = 0.007 ;
//fptype lo = x;
//fptype hi = x+0.00000000001 ;
//
  fptype intg_hi = 0.5*(1-EXP(tau*(tau*sigma*sigma/2.+mean - hi))*
                   erfc((tau+(mean-hi)/(sigma*sigma))*sigma/1.4142135623)+ erf((hi-mean)/(sigma*1.4142135623)));   
  fptype intg_lo = 0.5*(1-EXP(tau*(tau*sigma*sigma/2.+mean - lo))*
                   erfc((tau+(mean-lo)/(sigma*sigma))*sigma/1.4142135623)+ erf((lo-mean)/(sigma*1.4142135623)));  
		    
//  
//   fptype u_hi = tau * (hi - mean);
//   fptype v_hi = tau * sigma;
//   fptype expa_hi = -u_hi+ v_hi*v_hi*0.5+LOG(0.5*( 1+ERF((u_hi-v_hi*v_hi)/(v_hi*1.4142135623))) );
//   fptype intg_hi1 = (0.5*(1+ERF(u_hi/(v_hi*1.4142135623))) - EXP(expa_hi));
// //
//   fptype u_lo = tau * (lo - mean);
//   fptype v_lo = tau * sigma;
//   fptype expa_lo = -u_lo+ v_lo*v_lo*0.5+LOG(0.5*( 1+ERF((u_lo-v_lo*v_lo)/(v_lo*1.4142135623))) );
//   fptype intg_lo1 = (0.5*(1+ERF(u_lo/(v_lo*1.4142135623))) - EXP(expa_lo));
//
 
 
  fptype integral = fabs(intg_hi-intg_lo);
//  fptype integral1 = fabs(intg_hi1-intg_lo1);
   ret = ret/(integral);
 
//if ((0 == THREADIDX) && (0 == BLOCKIDX)){
//  printf("\n\nExpGausPEESigmaBPdf x=%f  sigma=%f mean=%f tau=%f ret=%f integral =%f\n", x, sigma, mean,tau , ret, integral);
//  printf("ExpGausPEESigmaBPdf x=%f  sigma=%f mean=%f tau=%f ret=%f integral1=%f\n", x, sigma, mean,tau , ret, integral1);
//}  
 
/*  if ((0 == THREADIDX) && (0 == BLOCKIDX)){
  printf("ExpGausPEESigmaBPdf->integrate  sigma=%f mean=%f tau=%f Integral=%f\n",    sigma, mean,tau , integral);
  printf("ExpGausPEESigmaBPdf->integrate  sigma=%f mean=%f tau=%f IntegralLO=%f\n",  sigma, mean,tau , intg_lo);
  printf("ExpGausPEESigmaBPdf->integrate  sigma=%f mean=%f tau=%f IntegralHI=%f\n",  sigma, mean,tau , intg_hi);
  printf("ExpGausPEESigmaBPdf->integrate  sigma=%f mean=%f tau=%f EXP(expa_hi)=%f\n",  sigma, mean,tau , EXP(expa_hi));
  printf("ExpGausPEESigmaBPdf->integrate  sigma=%f mean=%f tau=%f EXP(expa_lo)=%f\n",  sigma, mean,tau , EXP(expa_lo));
 }
 */// printf("ExpGausPEESigmaBPdf->host_indices  host_indices0=%d host_indices1=%d host_indices2=%d\n",  host_indices[0], host_indices[1],host_indices[2]);
// printf("ExpGausPEESigmaBPdf->     indices  indices0=%d indices1=%d indices2=%d\n",  indices[0], indices[1],indices[2]);
//  if ((0 == THREADIDX) && (0 == BLOCKIDX)){
//   printf("ExpGausPEESigmaBPdf x=%f  sigma=%f mean=%f tau=%f ret=%f integral=%f\n", x, sigma, mean,tau , ret, integral);
//    printf("ExpGausPEESigmaBPdf ind0 =%d  evt0 =%f \n",indices[0] ,evt[indices[0]] );
//    printf("ExpGausPEESigmaBPdf ind1 =%d  evt1 =%f \n",indices[1] ,evt[indices[1]] );
//    printf("ExpGausPEESigmaBPdf ind2 =%d  evt2 =%f \n",indices[2] ,evt[indices[2]] );
//    printf("ExpGausPEESigmaBPdf ind3 =%d  evt3 =%f \n",indices[3] ,evt[indices[3]] );
//    printf("ExpGausPEESigmaBPdf ind4 =%d  evt4 =%f \n",indices[4] ,evt[indices[4]] );
//    printf("ExpGausPEESigmaBPdf ind5 =%d  evt5 =%f \n",indices[5] ,evt[indices[5]] );
//    printf("ExpGausPEESigmaBPdf ind6 =%d  evt6 =%f \n",indices[6] ,evt[indices[6]] );
//    printf("ExpGausPEESigmaBPdf ind7 =%d  evt7 =%f \n",indices[7] ,evt[indices[7]] );
//    printf("ExpGausPEESigmaBPdf ind8 =%d  evt8 =%f \n",indices[8] ,evt[indices[8]] );
//    printf("ExpGausPEESigmaBPdf ind9 =%d  evt9 =%f \n",indices[9] ,evt[indices[9]] );
//    printf("ExpGausPEESigmaBPdf ind10=%d  evt10=%f \n",indices[10],evt[indices[10]] );
//    printf("ExpGausPEESigmaBPdf ind11=%d  evt11=%f \n",indices[11],evt[indices[11]] );
//    printf("ExpGausPEESigmaBPdf ind12=%d  evt12=%f \n",indices[12],evt[indices[12]] );
//    printf("ExpGausPEESigmaBPdf ind13=%d  evt13=%f \n",indices[13],evt[indices[13]] );
//    printf("ExpGausPEESigmaBPdf ind14=%d  evt14=%f \n",indices[14],evt[indices[14]] );
// }
//  return 0; 
//
  return ret;
}

MEM_DEVICE device_function_ptr ptr_to_ExpGausPEESigmaB = device_ExpGausPEESigmaB; 

 __host__ ExpGausPEESigmaBPdf::ExpGausPEESigmaBPdf (std::string n, Variable* _x, Variable* _s,  Variable* mean, Variable* tau, Variable* lo, Variable* hi  ) 
  : GooPdf(_x, n)
{
//   registerObservable(_x); //already registered!!!
//   registerObservable(_s);
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(_s));
  pindices.push_back(registerParameter(mean));
//  pindices.push_back(registerParameter(sigma));
  pindices.push_back(registerParameter(tau));
  pindices.push_back(registerParameter(lo));
  pindices.push_back(registerParameter(hi));
  GET_FUNCTION_ADDR(ptr_to_ExpGausPEESigmaB);
  initialise(pindices); 
}

 __host__ fptype ExpGausPEESigmaBPdf::integrate (fptype lo, fptype hi) const {
// printf("integratexxx\n");
/*  unsigned int* indices = host_indices+parameters; 
 fptype sigmaM = 0.0017; 
 fptype mean = host_params[indices[1]]  ;
 fptype tau  = host_params[indices[2]]  ;
 fptype u_hi = tau * (hi - mean);
 fptype v_hi = tau * sigmaM;
// fptype v_hi = tau * u_hi;
 fptype expa_hi = -u_hi+ v_hi*v_hi*0.5+LOG(0.5*( 1+ERF((u_hi-v_hi*v_hi)/(v_hi*1.4142135623))) );
 fptype intg_hi = (0.5*(1+ERF(u_hi/(v_hi*1.4142135623))) - EXP(expa_hi));

 fptype u_lo = tau * (lo - mean);
 fptype v_lo = tau * sigmaM;
// fptype v_lo = tau * u_lo;
 fptype expa_lo = -u_lo+ v_lo*v_lo*0.5+LOG(0.5*( 1+ERF((u_lo-v_lo*v_lo)/(v_lo*1.4142135623))) );
 fptype intg_lo = (0.5*(1+ERF(u_lo/(v_lo*1.4142135623))) - EXP(expa_lo));
 
// printf("ExpGausPEESigmaBPdf->integrate  sigma=%f mean=%f tau=%f Integral=%f\n",  sigmaM, mean,tau , (intg_hi-intg_lo));
// printf("ExpGausPEESigmaBPdf->host_indices  host_indices0=%d host_indices1=%d host_indices2=%d\n",  host_indices[0], host_indices[1],host_indices[2]);
// printf("ExpGausPEESigmaBPdf->     indices  indices0=%d indices1=%d indices2=%d\n",  indices[0], indices[1],indices[2]);
 return (intg_hi-intg_lo);
 */
  return 1.;
 }
