#include "SigmoidPdf.hh"

EXEC_TARGET fptype device_Sigmoid (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x     = evt[indices[2 + indices[0]]]; 
  fptype p0 = p[indices[1]];
  fptype p1 = p[indices[2]];
  fptype p2 = p[indices[3]];
  fptype p3 = p[indices[4]];

  fptype ret = p0/(1+p3+ p1*exp(-p2*x));

//  if ((0 == THREADIDX) && (0 == BLOCKIDX)){
//  printf("Sigmoid x=%f  sigma=%f mean=%f tau=%f ret=%f\n", x, sigma, mean,alpha , ret);
//  } 

  if (ret<=0) {
       printf("Sigmoid <=0!!!: x = %f , p0 = %f p1 = %f f\n",x,p0,p1,p2,p3);
       return 0.;
  }
  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_Sigmoid = device_Sigmoid; 

__host__ SigmoidPdf::SigmoidPdf (std::string n, Variable* _x, Variable* p0, Variable* p1, Variable* p2, Variable* p3) 
  : GooPdf(_x, n)
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(p0));
  pindices.push_back(registerParameter(p1));
  pindices.push_back(registerParameter(p2));
  pindices.push_back(registerParameter(p3));
  GET_FUNCTION_ADDR(ptr_to_Sigmoid);
  initialise(pindices); 
}


