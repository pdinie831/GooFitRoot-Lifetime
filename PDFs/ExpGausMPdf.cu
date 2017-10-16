#include "ExpGausMPdf.hh"

EXEC_TARGET fptype device_ExpGausM (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x     = evt[indices[2 + indices[0]]]; 
  fptype mean  = p[indices[1]];
//  fptype sigma = p[indices[2]];
  fptype sigma = evt[indices[3 + indices[0]]];;
//  fptype alpha = p[indices[3]];
  fptype alpha = p[indices[2]];

  fptype ret = 0.5*alpha; 
  fptype exparg = ret * (2*mean + alpha*sigma*sigma - 2*x);
  fptype erfarg = (mean + alpha*sigma*sigma - x) / (sigma * 1.4142135623);

  ret *= EXP(exparg); 
  ret *= ERFC(erfarg); 
//  if ((0 == THREADIDX) && (0 == BLOCKIDX))
//  printf("ExpGausMPdf x=%f  sigma=%f mean=%f tau=%f ret=%f\n", x, sigma, mean,alpha , ret);
//  return 0.;

  if (ret<=0) {
       printf("ExpGausMPdf <=0!!!: x = %f , mean = %f sigma = %f alpha = %f\n",x,mean,sigma,alpha);
//       return 0.;
  }
  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_ExpGausM = device_ExpGausM; 

ExpGausMPdf::ExpGausMPdf (std::string n, Variable* _x, Variable* sigma, Variable* mean, Variable* tau) 
  : GooPdf(0, n)
{
   registerObservable(_x); //already registered!!!
   registerObservable(sigma);
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(mean));
//  pindices.push_back(registerParameter(sigma));
  pindices.push_back(registerParameter(tau));
  GET_FUNCTION_ADDR(ptr_to_ExpGausM);
  initialise(pindices); 
}


