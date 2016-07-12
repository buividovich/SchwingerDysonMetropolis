#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <clue_logs.h>
#include <rand_num_generators.h>

#include <sd_metropolis.h>
#include <lattice_propagator.h>

#include "parameters.h"
#include "actions.h"
#include "statistics.h"

#include <fcntl.h>
#include <unistd.h>
#include <limits.h>

#ifdef TESTS
 #include "test_stereo_vertex.h"
#endif

//A good set ob parameters: --lambda 1.31 --alpha 0.015 --average-seq-len 1.5 --average-num-seq 1.1

int main(int argc, char *argv[])
{
 ansi_colors = 1;
 print_errors_to_stderr = 0;
  
 logs_Write(0, "\n\nDIAGRAMMATIC MONTE-CARLO FOR SU(N) SIGMA MODEL IN THE PLANAR LIMIT - BASED ON THE STEREOGRAPHIC PROJECTION\n");
 
 init_rand(time(NULL));
 
 parse_command_line_options(argc, argv);
 init_actions();
 init_parameters();

 init_metropolis();
 t_observable_stat* obs_stat = init_observable_stat();
 t_stack_stat* Xstack_stat = init_stack_statistics(max_stack_nel);
 
 print_parameters();
 
 logs_Write(0, "Starting thermalization process for %i MC steps", therm_mc_steps);
 for(int imc=0; imc<therm_mc_steps; imc++)
  metropolis_step(imc);
 
 logs_Write(0, "Starting production run for %i MC steps", prod_mc_steps);
 for(int imc=0; imc<prod_mc_steps; imc++)
 {
  metropolis_step(imc);
  if(imc%mc_interval==0)
  {
   gather_observable_stat(obs_stat); 
   gather_stack_statistics(Xstack_stat, &X);
  }; 
 };
 
 print_stack_statistics(Xstack_stat);
  
 char prefix[500];
 largeN_QFT_prefix(prefix);
 process_mc_stat(prefix, 1);
 process_observable_stat(obs_stat);
 
 free_observable_stat(obs_stat);
 free_stack_statistics(Xstack_stat);
 SAFE_FREE(obs_stat);
 SAFE_FREE(Xstack_stat);
 
 free_metropolis();
 free_actions();
 free_lat_propagator(&P);
 free_parameters();
 
 return 0;
}
