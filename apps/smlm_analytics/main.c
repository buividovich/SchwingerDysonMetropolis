#include <stdio.h>
#include <stdlib.h>

#include <largeN_QFT_parameters.h>

#include "parameters.h"
#include "recursion.h"

//TODO Plan of implementation:
// 1. Fast implementation for any LS using uint-packing, without buffering
// 2. Implement buffering
// 3. Implement fast operations for LS==2 using #defines/ inline functions

int main(int argc, char *argv[])
{
 ansi_colors = 1;
 print_errors_to_stdout = 1;
 parse_command_line_options(argc, argv);
 init_parameters();
 print_parameters();
 init_recursion();
 
 DECLARE_AND_MALLOC(Gpq, double, LS);
 DECLARE_AND_MALLOC( sU, double, (mmax+1));
 DECLARE_AND_MALLOC( sL, double, (mmax+1));
 
 logs_Write(0, "\nExpansion coefficients for G2: ");
 
 for(int io=0; io<=mmax; io++)
 {
  double U = 0.0;
  double L = 0.0;
  for(uint p=0; p<LS; p++)
  {
   Gpq[p]  = G(p*LS + (LS-p)%LS, 1, io);
   U      += Gpq[p];
   L      += Gpq[p]*cos(2.0*M_PI*(double)p/(double)LS);
  };
  sU[io]   = pow(lambda, (double)(io+1))*U + (io>0? sU[io-1] : 0.0);
  sL[io]   = pow(lambda, (double)(io+1))*L + (io>0? sL[io-1] : 0.0);
  logs_Write(0, "\t Order %02i: U = %+2.4lf\tL = %+2.4lf", io, sU[io], sL[io]);
 };
 logs_Write(0, ""); 
 
 SAFE_FREE(Gpq);
 SAFE_FREE(sU);
 SAFE_FREE(sL);
 free_recursion();
 return EXIT_SUCCESS;
}
