#include <stdio.h>
#include <stdlib.h>

#include <largeN_QFT_parameters.h>

#include "parameters.h"
#include "recursion.h"

//TODO Plan of implementation:
// 4. Instead of recursion, we can do 
// a) looping 
// b) save the result only for conserved momenta
// c) inline funcs for m=1, separate calculation for m=2

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
 
 //Opening the output file
 FILE* f = NULL;
 if(out_file!=NULL)
  f = (append_mode? fopen(out_file, "a") : fopen(out_file, "w"));
 if(f==NULL && out_file!=NULL)
  logs_WriteError("Could not open the file %s for %s", out_file, (append_mode? "writing" : "appending"));
 
 //testG4(2);
 logs_Write(0, "\nExpansion coefficients for G2: ");
 
 for(int io=0; io<=mmax; io++)
 {
  double U = 0.0;
  double L = 0.0;
  getG2(Gpq, io);
  for(uint p=0; p<LS; p++)
  {
   U      += Gpq[p];
   L      += Gpq[p]*cos(2.0*M_PI*(double)p/(double)LS);
  };
  sU[io]   = pow(lambda, (double)(io+1))*U + (io>0? sU[io-1] : 0.0);
  sL[io]   = pow(lambda, (double)(io+1))*L + (io>0? sL[io-1] : 0.0);
  logs_Write(0, "\t Order %02i: U = %+2.6lf\tL = %+2.6lf", io, sU[io], sL[io]);
  if(f!=NULL)
   fprintf(f, "%2.4lf % 4i %2.4lf %2.4lf %2.4lf\n", lambda, io, 1.0/(double)(io + 1), sU[io], sL[io]);
 };
 logs_Write(0, ""); 
 
 print_buffer_usage();
 
 if(f!=NULL)
 {
  fprintf(f, "%2.4lf % 5i %2.4lf %2.4lf %2.4lf\n", lambda, 9999, 0.0, 1.0, mean_link);
  fclose(f);        
 };
 
 SAFE_FREE(out_file);
 
 SAFE_FREE(Gpq);
 SAFE_FREE(sU);
 SAFE_FREE(sL);
 free_recursion();
 return EXIT_SUCCESS;
}
