#include "ExpGausWithIntPdf.hh"

EXEC_TARGET fptype device_ExpGausWithInt (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x     = evt[indices[2 + indices[0]]]; 
  fptype mean  = p[indices[1]];
  fptype sigma = p[indices[2]];
  fptype tau = p[indices[3]];

  fptype ret = 0.5*tau; 
  fptype exparg = ret * (2*mean + tau*sigma*sigma - 2*x);
  fptype erfarg = (mean + tau*sigma*sigma - x) / (sigma * 1.4142135623);

  ret *= EXP(exparg); 
  ret *= ERFC(erfarg); 
//  if ((0 == THREADIDX) && (0 == BLOCKIDX))
//  printf("ExpGausWithIntPdf x=%f  sigma=%f mean=%f tau=%f ret=%f\n", x, sigma, mean,alpha , ret);
//  return 0.;

  if (ret<=0) {
       printf("ExpGausWithIntPdf <=0!!!: x = %f , mean = %f sigma = %f alpha = %f\n",x,mean,sigma,tau);
//       return 0.;
  }
  fptype lo = 0.0003;
  fptype hi = 0.007;
  fptype intg_hi = 0.5*(1-EXP(tau*(tau*sigma*sigma/2.+mean - hi))*
                   erfc((tau+(mean-hi)/(sigma*sigma))*sigma/1.4142135623)+ erf((hi-mean)/(sigma*1.4142135623)));   
  fptype intg_lo = 0.5*(1-EXP(tau*(tau*sigma*sigma/2.+mean - lo))*
                   erfc((tau+(mean-lo)/(sigma*sigma))*sigma/1.4142135623)+ erf((lo-mean)/(sigma*1.4142135623)));   
//   fptype u_hi = tau * (hi - mean);
//   fptype v_hi = tau * sigma;
//   fptype expa_hi = -u_hi+ v_hi*v_hi*0.5+LOG(0.5*( 1+ERF((u_hi-v_hi*v_hi)/(v_hi*1.4142135623))) );
//   fptype intg_hi = (0.5*(1+ERF(u_hi/(v_hi*1.4142135623))) - EXP(expa_hi));
// 
//   fptype u_lo = tau * (lo - mean);
//   fptype v_lo = tau * sigma;
//   fptype expa_lo = -u_lo+ v_lo*v_lo*0.5+LOG(0.5*( 1+ERF((u_lo-v_lo*v_lo)/(v_lo*1.4142135623))) );
//   fptype intg_lo = (0.5*(1+ERF(u_lo/(v_lo*1.4142135623))) - EXP(expa_lo));

  fptype integral = (intg_hi-intg_lo);
  return ret/integral; 
}

MEM_DEVICE device_function_ptr ptr_to_ExpGausWithInt = device_ExpGausWithInt; 

ExpGausWithIntPdf::ExpGausWithIntPdf (std::string n, Variable* _x, Variable* mean, Variable* sigma, Variable* tau) 
  : GooPdf(_x, n)
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(mean));
  pindices.push_back(registerParameter(sigma));
  pindices.push_back(registerParameter(tau));
  GET_FUNCTION_ADDR(ptr_to_ExpGausWithInt);
  initialise(pindices); 
}

 __host__ fptype ExpGausWithIntPdf::integrate (fptype lo, fptype hi) const {
/*  unsigned int* indices = host_indices+parameters; 
 fptype mean = host_params[indices[1]]  ;
 fptype sigma= host_params[indices[2]]; 
 fptype tau  = host_params[indices[3]]  ;
 fptype u_hi = tau * (hi - mean);
 fptype v_hi = tau *sigma;
// fptype v_hi = tau * u_hi;
 fptype expa_hi = -u_hi+ v_hi*v_hi*0.5+LOG(0.5*( 1+ERF((u_hi-v_hi*v_hi)/(v_hi*1.4142135623))) );
 fptype intg_hi = (0.5*(1+ERF(u_hi/(v_hi*1.4142135623))) - EXP(expa_hi));

 fptype u_lo = tau * (lo - mean);
 fptype v_lo = tau * sigma;
// fptype v_lo = tau * u_lo;
 fptype expa_lo = -u_lo+ v_lo*v_lo*0.5+LOG(0.5*( 1+ERF((u_lo-v_lo*v_lo)/(v_lo*1.4142135623))) );
 fptype intg_lo = (0.5*(1+ERF(u_lo/(v_lo*1.4142135623))) - EXP(expa_lo));
 return (intg_hi-intg_lo);
 */
   return 1.;
 }

