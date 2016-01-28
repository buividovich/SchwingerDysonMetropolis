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
 
 //init_lattice_constants();
 //alloc_recursion_packed();
 
 //return EXIT_SUCCESS;
 
 //TODO: may be the sequences in the address book are always containing the sequences for previous n, 
 //Is it unnecessary to initialize n times?
 //TODO: still count multiplicities? For Gamma coefficients... Is it useful at all?
  
 //init_generalized_sc();
 //run_generalized_sc();
 init_stereographic();
  run_stereographic();
 if(consistency_check)
  test_recursion();
 //print_Gmn(); 
 
 //Calculating the final answers and saving to the file
 char* out_file1 = NULL;
 sprintf_append(&out_file1, "G:\\LAT\\sd_metropolis\\data\\rcan\\stereo_conv_l%2.4lf.dat", lambda);
 FILE* f = NULL;
 if(out_file1!=NULL)
  f = (append_mode? fopen(out_file1, "a") : fopen(out_file1, "w"));
 if(f==NULL)
  logs_WriteErrorAndTerminate("Could not open the file %s for %s", out_file1, (append_mode? "appending" : "writing"));
 
  
 DECLARE_AND_MALLOC(Gxy, double, (mmax+1)*LS);
 DECLARE_AND_MALLOC( Gx, double, (mmax+1)   );
 get_Gxy_stereographic(Gxy, Gx);
 
 double G01_exact = (lambda<4.0? 1.0 - 0.125*lambda : 2.0/lambda);
 char G01exstr[64];
 if(LS==2)
  sprintf(G01exstr, "%2.6lf", G01_exact);
 else
  sprintf(G01exstr, "not known for LS=%i", LS);
 
 logs_Write(0, "\nExpansion results (increasing order), exact result for G01 is %s, saving to file %s", G01exstr, out_file1);

 for(int mmax1=0; mmax1<=mmax; mmax1++)
 {
  char* rstr = NULL;
  for(int x=0; x<LS; x++)
   sprintf_append(&rstr, "G%i%i = %+2.6lf,\t", 0, x, Gxy[mmax1*LS + x]); 
  if(LS==2)
  { 
   double err = (Gxy[mmax1*LS + 1]/Gxy[mmax1*LS + 0] - G01_exact)/G01_exact;
   sprintf_append(&rstr, " err = %+2.4E,\t", err);
   fprintf(f, "%02i %2.4E %2.4E %2.4E %2.4E\n", mmax1, 1.0/(double)(mmax1+1), Gxy[mmax1*LS + 1], err, Gx[mmax1]);
  }; 
  sprintf_append(&rstr, " Gx = %+2.4E", Gx[mmax1]);
  logs_Write(0, "\t Order %02i: %s ", mmax1, rstr);
 };

 fclose(f);
 
 free_recursion();
  
 SAFE_FREE(Gxy);
 SAFE_FREE(Gx);
 return EXIT_SUCCESS;
}

/*
  
 logs_Write(0, "\nExpansion results (increasing order), exact result for G01 is %s, saving to file %s", G01exstr, out_file1);
 fprintf(f, "%2.6E ", meff_sq);
 for(int mmax1=0; mmax1<=mmax; mmax1++)
 {
  //fprintf(f, "%02i %2.4E ", m, 1.0/(double)(m+1));
  //for(int x=0; x<LS; x++)
  // fprintf(f, "%+2.4E ", Gxy[m*LS + x]);
  char* rstr = NULL;
  for(int x=0; x<LS; x++)
   sprintf_append(&rstr, "G%i%i = %+2.6lf,\t", 0, x, Gxy[mmax1*LS + x]); 
  if(LS==2)
  { 
   double err = (Gxy[mmax1*LS + 1]/Gxy[mmax1*LS + 0] - G01_exact)/G01_exact;
   sprintf_append(&rstr, " err = %+2.4E,\t", err);
   fprintf(f, "%+2.4E ", err);
  }; 
  sprintf_append(&rstr, " Gx = %+2.4E", Gx[mmax1]);
  logs_Write(0, "\t Order %02i: %s ", mmax1, rstr);
  //fprintf(f, "\n");
 };
 fprintf(f, "\n");
 
fclose(f);*/
