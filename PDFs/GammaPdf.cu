#include "GammaPdf.hh"

EXEC_TARGET fptype device_Gamma (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x       = evt[indices[2 + indices[0]]];
  fptype gamma   = p[indices[1]];
  fptype beta    = p[indices[2]];
  fptype mu      = p[indices[3]];
  if ((x<mu) || (gamma<=0) || (beta <=0)) {
       printf("GammaPdf illegal parameter values x = %f , gamma = %f beta = %f\n",x,gamma,beta);
       return 0.;
  }
  fptype arg1    = x-mu;
  fptype arg2    = gamma-1;
  fptype ret     = 0;
  
  if ((arg1) < 0) {
     ret = 0.0;
  } else if (arg1 == 0) {
   if (gamma == 1) {
    ret = 1.0/beta;
   } else {
    ret = 0.0;
   }
  } else if (gamma == 1) {
    ret = EXP(-(arg1)/beta)/beta;
  } else {
    ret = EXP((arg2) * log((arg1)/beta) - (arg1)/beta - lgamma(gamma))/beta;
  }

//  fptype ret = POW(arg1,arg2)*EXP((-arg1)/beta)/(tgamma(gamma)*POW(beta,gamma));

//  if ((0 == THREADIDX) && (0 == BLOCKIDX))
//   printf("device_Gamma x=%f gamma=%f beta=%f mu=%f  ret=%f\n", x, gamma, beta,mu , ret);

  if (ret<=0) {
       printf("GammaPdf <=0!!!: x = %f , gamma = %f beta = %f mu = %f\n",x,gamma,beta,mu);
       return 0.;
  }
  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_Gamma = device_Gamma; 

__host__ GammaPdf::GammaPdf (std::string n, Variable* _x, Variable* gamma, Variable* beta, Variable* mu  ) 
  : GooPdf(_x, n) 
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(gamma));
  pindices.push_back(registerParameter(beta));
  pindices.push_back(registerParameter(mu));
  GET_FUNCTION_ADDR(ptr_to_Gamma);
  initialise(pindices); 
}

//__host__ fptype GammaPdf::integrate (fptype lo, fptype hi) const {
  
//   unsigned int* indices = host_indices+parameters; 
//   fptype p0 = host_params[indices[1]];
//   fptype p1 = host_params[indices[2]];
//   
//   fptype xmin =host_params[indices[3]]; 
//   fptype xmax =host_params[indices[4]];
//   
//   fptype xMinL  = -1+2*(xmin-lo)/(xmax-xmin); 
//   fptype xMaxL  =  1+2*(hi-xmax)/(xmax-xmin); 
// 
// //  fptype  xdelta1 = xMaxL - xMinL;
// //  fptype  xdelta2 = xMaxL*xMaxL - xMinL*xMinL;
// //  fptype  xdelta3 = xMaxL*xMaxL*xMaxL - xMinL*xMinL*xMinL;
//   
//   
//   fptype xint = 0.5*(hi-lo)*((xMaxL + p0*xMaxL*xMaxL/2+p1*(2*xMaxL*xMaxL*xMaxL/3-xMaxL))- (xMinL +  p0*xMinL*xMinL/2+p1*(2*xMinL*xMinL*xMinL/3-xMinL)));
// 
// //   printf("integ_Gamma p0=%f p1=%f xmin=%f xmax=%f lo=%f hi=%f int=%f\n",  p0, p1,xmin,xmax, lo,hi,xint);
// 
//   return xint; 
//    return 1;
//}

