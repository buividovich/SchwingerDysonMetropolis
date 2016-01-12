#include <stdio.h>
#include <stdlib.h>

#include <largeN_QFT_parameters.h>

#include "common_parameters.h"
#include "recursion_commons.h"

#include "generalized_sc.h"

//TODO Plan of implementation:
// 4. Instead of recursion, we can do 
// a) looping 
// b) save the result only for conserved momenta
// c) inline funcs for m=1, separate calculation for m=2

int main(int argc, char *argv[])
{
 ansi_colors = 1;
 print_errors_to_stdout = 0;
 parse_common_command_line_options(argc, argv);
 init_common_parameters();
 print_common_parameters();
 
 init_generalized_sc();
 run_generalized_sc();
 test_generalized_sc();
 
 //Calculating the final answers and saving to the file
 logs_Write(0, "\nExpansion results (increasing order): ");
 DECLARE_AND_MALLOC( sU, double, (mmax+1));
 DECLARE_AND_MALLOC( sL, double, (mmax+1));
 for(int m=0; m<=mmax; m++)
 {
  double U = 0.0; double L = 0.0;
  for(uint p=0; p<LS; p++)
  {
   uint P = LS*((LS-p)%LS) + p;
   double Gpq = get_generalized_sc(m, 1, P);
   U += Gpq;
   L += Gpq*cos(2.0*M_PI*(double)p/(double)LS);
  };
  sU[m] = pow(lambda, (double)(m+1))*U + (m>0? sU[m-1] : 0.0);
  sL[m] = pow(lambda, (double)(m+1))*L + (m>0? sL[m-1] : 0.0);
  logs_Write(0, "\t Order %02i: U = %+2.6E\tL = %+2.6E", m, sU[m], sL[m]);
  //if(f!=NULL)
  // fprintf(f, "%2.4lf % 4i %2.4lf %2.4lf %2.4lf\n", lambda, io, 1.0/(double)(io + 1), sU[io], sL[io]);
 };
 
 free_generalized_sc();
 
 /*
 //Opening the output file
 FILE* f = NULL;
 if(out_file!=NULL)
  f = (append_mode? fopen(out_file, "a") : fopen(out_file, "w"));
 if(f==NULL && out_file!=NULL)
  logs_WriteError("Could not open the file %s for %s", out_file, (append_mode? "writing" : "appending"));
 
 logs_Write(0, "\nExpansion coefficients for G2: ");
 
 
 logs_Write(0, ""); 
 
 if(f!=NULL)
 {
  fprintf(f, "%2.4lf % 5i %2.4lf %2.4lf %2.4lf\n", lambda, 9999, 0.0, 1.0, 0.0); //TODO: here should go the exact answer
  fclose(f);        
 };
 
 SAFE_FREE(out_file); */
 
 SAFE_FREE(sU);
 SAFE_FREE(sL);
 return EXIT_SUCCESS;
}
