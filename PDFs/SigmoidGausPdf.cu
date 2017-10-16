#include "SigmoidGausPdf.hh"

EXEC_TARGET fptype device_SigmoidGaus (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x     = evt[indices[2 + indices[0]]]; 
  fptype p0 = p[indices[1]];
  fptype p1 = p[indices[2]];
  fptype p2 = p[indices[3]];
  fptype p3 = p[indices[4]];
  fptype p4    = p[indices[5]];
  fptype mean  = p[indices[6]];
  fptype sigma = p[indices[7]];

  fptype ret = p0/(1+p3+ p1*exp(-p2*x))+p4*EXP(-0.5*(x-mean)*(x-mean)/(sigma*sigma));

//  if ((0 == THREADIDX) && (0 == BLOCKIDX)){
//  printf("SigmoidGaus x=%f  sigma=%f mean=%f tau=%f ret=%f\n", x, sigma, mean,alpha , ret);
//  } 

  if (ret<=0) {
       printf("SigmoidGaus <=0!!!: x = %f , p0 = %f p1 = %f f\n",x,p0,p1,p2,p3,p4,mean,sigma);
       return 0.;
  }
  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_SigmoidGaus = device_SigmoidGaus; 

__host__ SigmoidGausPdf::SigmoidGausPdf (std::string n, Variable* _x, Variable* p0, Variable* p1, Variable* p2, Variable* p3, Variable* p4,
Variable* mean, Variable* sigma) 
  : GooPdf(_x, n)
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(p0));
  pindices.push_back(registerParameter(p1));
  pindices.push_back(registerParameter(p2));
  pindices.push_back(registerParameter(p3));
  pindices.push_back(registerParameter(p4));
  pindices.push_back(registerParameter(mean));
  pindices.push_back(registerParameter(sigma));
  GET_FUNCTION_ADDR(ptr_to_SigmoidGaus);
  initialise(pindices); 
}


