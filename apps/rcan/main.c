#include <stdio.h>
#include <stdlib.h>

#include <largeN_QFT_parameters.h>

#include "common_parameters.h"
#include "recursion.h"

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
 if(consistency_check)
  test_recursion();
 
 //Calculating the final answers and saving to the file
 char* out_file1 = NULL;
 sprintf_append(&out_file1, "G:\\LAT\\sd_metropolis\\data\\rcan\\mscan_l%2.4lf.dat", lambda);
 FILE* f = NULL;
 if(out_file1!=NULL)
  f = (append_mode? fopen(out_file1, "a") : fopen(out_file1, "w"));
 if(f==NULL)
  logs_WriteErrorAndTerminate("Could not open the file %s for %s", out_file1, (append_mode? "appending" : "writing"));
 
  
 logs_Write(0, "\nExpansion results (increasing order), saving to file %s", out_file1);
 DECLARE_AND_MALLOC(Gxy, double, (mmax+1)*LS);
 get_Gxy_generalized_sc(Gxy);
 fprintf(f, "%2.6E ", meff_sq);
 for(int m=0; m<=mmax; m++)
 {
  //fprintf(f, "%02i %2.4E ", m, 1.0/(double)(m+1));
  //for(int x=0; x<LS; x++)
  // fprintf(f, "%+2.4E ", Gxy[m*LS + x]);
  char* rstr = NULL;
  for(int x=0; x<LS; x++)
   sprintf_append(&rstr, "G%i%i = %+2.4E,\t", 0, x, Gxy[m*LS + x]); 
  if(LS==2)
  { 
   double G01_exact = (lambda<4.0? 1.0 - 0.125*lambda : 2.0/lambda);
   double err = (Gxy[m*LS + 1]/Gxy[m*LS + 0] - G01_exact)/G01_exact;
   sprintf_append(&rstr, " err = %+2.4E", err);
   fprintf(f, "%+2.4E ", err);
  }; 
  logs_Write(0, "\t Order %02i: %s ", m, rstr);
  //fprintf(f, "\n");
 };
 fprintf(f, "\n");
 
 fclose(f);
 free_recursion();
  
 SAFE_FREE(Gxy);
 return EXIT_SUCCESS;
}
