#include "SimpleCheby2Pdf.hh"

EXEC_TARGET fptype device_SimpleCheby2 (fptype* evt, fptype* p, unsigned int* indices) {
  fptype xx   = evt[indices[2 + indices[0]]];
  fptype p0   = p[indices[1]];
  fptype p1   = p[indices[2]];
  fptype xmin = p[indices[3]];
  fptype xmax = p[indices[4]];
  fptype x    = -1+2*(xx-xmin)/(xmax-xmin); 

  fptype ret = 1+p0*x+p1*(2*x*x-1);

  //if ((0 == THREADIDX) && (0 == BLOCKIDX)) cuPrintf("SimpleCheby2 Values %f %i %i %f %f %i\n", x, indices[1], indices[2], p0, p1, callnumber); 
  //cuPrintf("device_SimpleCheby2 %f %i %i %f %f %i %p %f\n", x, indices[1], indices[2], p0, p1, callnumber, indices, ret); 
  //if ((0 == THREADIDX) && (0 == BLOCKIDX))
  //printf("device_SimpleCheby2 %f %f %f %i %f\n", x, p0, p1, callnumber, ret);     
//  if ((0 == THREADIDX) && (0 == BLOCKIDX))
//   printf("device_SimpleCheby2 x=%f p0=%f p1=%f xmin=%f xmax=%f ret=%f\n", x, p0, p1,xmin,xmax, ret);

  return ret; 
}

MEM_DEVICE device_function_ptr ptr_to_SimpleCheby2 = device_SimpleCheby2; 

__host__ SimpleCheby2Pdf::SimpleCheby2Pdf (std::string n, Variable* _x, Variable* p0, Variable* p1,Variable* xmin, Variable* xmax ) 
  : GooPdf(_x, n) 
{
  std::vector<unsigned int> pindices;
  pindices.push_back(registerParameter(p0));
  pindices.push_back(registerParameter(p1));
  pindices.push_back(registerParameter(xmin));
  pindices.push_back(registerParameter(xmax));
  GET_FUNCTION_ADDR(ptr_to_SimpleCheby2);
  initialise(pindices); 
}

__host__ fptype SimpleCheby2Pdf::integrate (fptype lo, fptype hi) const {
  
  unsigned int* indices = host_indices+parameters; 
  fptype p0 = host_params[indices[1]];
  fptype p1 = host_params[indices[2]];
  
  fptype xmin =host_params[indices[3]]; 
  fptype xmax =host_params[indices[4]];
  
  fptype xMinL  = -1+2*(xmin-lo)/(xmax-xmin); 
  fptype xMaxL  =  1+2*(hi-xmax)/(xmax-xmin); 

//  fptype  xdelta1 = xMaxL - xMinL;
//  fptype  xdelta2 = xMaxL*xMaxL - xMinL*xMinL;
//  fptype  xdelta3 = xMaxL*xMaxL*xMaxL - xMinL*xMinL*xMinL;
  
  
  fptype xint = 0.5*(hi-lo)*((xMaxL + p0*xMaxL*xMaxL/2+p1*(2*xMaxL*xMaxL*xMaxL/3-xMaxL))- (xMinL +  p0*xMinL*xMinL/2+p1*(2*xMinL*xMinL*xMinL/3-xMinL)));

//   printf("integ_SimpleCheby2 p0=%f p1=%f xmin=%f xmax=%f lo=%f hi=%f int=%f\n",  p0, p1,xmin,xmax, lo,hi,xint);

  return xint; 
}

