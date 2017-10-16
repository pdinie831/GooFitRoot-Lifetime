setenv PATH /usr/local/cuda-8.0/bin:${PATH}
if (  ($?LD_LIBRARY_PATH)  ) then
setenv LD_LIBRARY_PATH /usr/local/cuda-8.0/lib64:${LD_LIBRARY_PATH}
else
setenv LD_LIBRARY_PATH /usr/local/cuda-8.0/lib64
endif 
setenv LD_LIBRARY_PATH $PWD/../../rootstuff:${LD_LIBRARY_PATH}
setenv ROOTSYS /usr/local/root
