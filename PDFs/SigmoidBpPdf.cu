#include "SigmoidBpPdf.hh"

EXEC_TARGET fptype device_SigmoidBp (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x     = evt[indices[2 + indices[0]]]; 
  fptype p0 = p[indices[1]];
  fptype p1 = p[indices[2]];
  fptype p2 = p[indices[3]];
  fptype p3 = p[indices[4]];
  fptype p4 = p[indices[5]];
  fptype p5 = p[indices[6]];
  fptype p6 = p[indices[7]];
  fptype p7 = p[indices[8]];

  if (x<=0){
       printf("SigmoidBp not defined for x<=0!!!: x = %f",x);
       return 0.;
  }
  fptype ret = p2+p3*x+ p4*x*x + p5*x*x*x + p6*x*x*x*x + p7*x*x*x*x*x + p0/(1+exp(-p1*x));
//  fptype ret = p0/(1+p3+ p1*exp(-p2*x));

//  if ((0 == THREADIDX) && (0 == BLOCKIDX)){
//  printf("SigmoidBp x=%f  sigma=%f mean=%f tau=%f ret=%f\n", x, sigma, mean,alpha , ret);
//  } 

  if (ret<=0) {
       printf("SigmoidBp <=0!!!: ret = %f x = %f , p0 = %f p1 = %f p2 = %f p3 = %f p4 = %f p5 = %f p6 = %f p7 = %f\n",ret,x,p0,p1,p2,p3,p4,p5,p6,p7);
       return 0.;
  }
  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_SigmoidBp = device_SigmoidBp; 

__host__ SigmoidBpPdf::SigmoidBpPdf (std::string n, Variable* _x, Variable* p0, Variable* p1, Variable* p2, Variable* p3, Variable* p4, Variable* p5, Variable* p6, Variable* p7) 
  : GooPdf(_x, n)
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(p0));
  pindices.push_back(registerParameter(p1));
  pindices.push_back(registerParameter(p2));
  pindices.push_back(registerParameter(p3));
  pindices.push_back(registerParameter(p4));
  pindices.push_back(registerParameter(p5));
  pindices.push_back(registerParameter(p6));
  pindices.push_back(registerParameter(p7));
  GET_FUNCTION_ADDR(ptr_to_SigmoidBp);
  initialise(pindices); 
}


