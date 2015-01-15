#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "parameters.h"
#include "my_stack.h"
#include "metropolis.h"
#include "statistics.h"

int main(int argc, char *argv[])
{
 int imc;
 ansi_colors = 1;
 print_errors_to_stdout = 0;
 logs_noise_level = 1;
 
 parse_command_line_options(argc, argv);
 init_parameters();
 print_parameters();
 
 init_stack();
 init_mc_stat();
 init_observable_stat();
 
 int* ns_history = NULL;
 if(ns_history_file!=NULL)
  SAFE_MALLOC(ns_history, int, nmc);
 
 overflow_count = 0;
 
 for(imc=0; imc<nmc; imc++)
 {
  if(mode==0)
   aac   += mc_step_sc();
  if(mode==1)
   aac   += mc_step_wc();
   
  gather_observable_stat();
  gather_mc_stat();
  
  if(ns_history!=NULL)
   ns_history[imc] = ns;
  logs_Write((imc%100000==0? 1 : 2), "step %i of %i, \t ns = %i \t stop[ns] = %i \t asign[ns] = %i", imc, nmc, ns, stop[ns], asign[ns]);
 };
 
 process_mc_stat();
 if(mode==0)
  process_observable_stat_sc();
 if(mode==1)
  process_observable_stat_wc(); 
 
 free_mc_stat();
 free_observable_stat();
 free_stack();
 
 if(ns_history!=NULL)
 {
  FILE* fnsh = fopen(ns_history_file, "wb");
  int snmc = fwrite (ns_history, sizeof(int), nmc, fnsh);
  if(snmc!=nmc)
   logs_WriteError("Could not save the history of ns to the file %s, only %i elements out of %i saved", ns_history_file, snmc, nmc);
  fclose(fnsh); 
  SAFE_FREE(ns_history);
 };
  
 return EXIT_SUCCESS;
}
