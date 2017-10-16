#include "ErfcPolyPdf.hh"

EXEC_TARGET fptype device_ErfcPoly (fptype* evt, fptype* p, unsigned int* indices) {
  fptype x     = evt[indices[2 + indices[0]]]; 
  fptype p0 = p[indices[1]];
  fptype p1 = p[indices[2]];
  fptype p2 = p[indices[3]];
  fptype p3 = p[indices[4]];
  fptype p4 = p[indices[5]];

  fptype erfarg = p0+p1*x+p2*x*x+p3*x*x*x+p4*x*x*x*x;

  fptype ret = ERFC(erfarg); 
//  if ((0 == THREADIDX) && (0 == BLOCKIDX)){
//  printf("ErfcPoly x=%f  sigma=%f mean=%f tau=%f ret=%f\n", x, sigma, mean,alpha , ret);
//  } 

  if (ret<=0) {
       printf("ErfcPoly <=0!!!: x = %f , p0 = %f p1 = %f p2 = %f p3 = %f p4 = %f\n",x,p0,p1,p2,p3,p4);
//       return 0.;
  }
  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_ErfcPoly = device_ErfcPoly; 

__host__ ErfcPolyPdf::ErfcPolyPdf (std::string n, Variable* _x, Variable* p0, Variable* p1, Variable* p2, Variable* p3, Variable* p4) 
  : GooPdf(_x, n)
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(p0));
  pindices.push_back(registerParameter(p1));
  pindices.push_back(registerParameter(p2));
  pindices.push_back(registerParameter(p3));
  pindices.push_back(registerParameter(p4));
  GET_FUNCTION_ADDR(ptr_to_ErfcPoly);
  initialise(pindices); 
}


