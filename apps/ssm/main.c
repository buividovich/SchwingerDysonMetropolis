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

/*void lambda_descent(int tune_mc_steps, int tune_iterations)
{
 find_cc_NN_minimum(param_tuning_accuracy, NULL);
 int prev_logs_noise_level = logs_noise_level;
 logs_noise_level = -1;
 
 for(int iter=0; iter<tune_iterations; iter++)
 {
  double mnA = tune_cc_NN_minimum(1.0E-4, tune_mc_steps);
  while(mnA<0.9)
  {
   lambda *= 0.9;
   init_metropolis();
   for(int imc=0; imc<tune_mc_steps; imc++)
    metropolis_step(imc);
   process_mc_stat("", 0);
   mnA = mean_nA;
   logs_Write(-1, "At lambda = %2.4E, <nA> = %2.4E", lambda, mean_nA); 
  };
  lambda /= 0.9;
 }; 
 logs_noise_level = prev_logs_noise_level;    
}*/

int main(int argc, char *argv[])
{
 ansi_colors = 1;
 print_errors_to_stderr = 0;
 
 logs_Write(0, "\n\nDIAGRAMMATIC MONTE-CARLO FOR SU(N) SIGMA MODEL IN THE PLANAR LIMIT - BASED ON SC EXPANSION\n");
 
 init_rand(time(NULL));
 
 parse_command_line_options(argc, argv);
 init_actions();
 init_parameters();

 init_metropolis();
 init_observable_stat();
 
 print_parameters();
 
 //lambda_descent(1000000, 20);
 
 logs_Write(0, "Starting thermalization process for %i MC steps", therm_mc_steps);
 for(int imc=0; imc<therm_mc_steps; imc++)
  metropolis_step(imc);
 
 logs_Write(0, "Starting production run for %i MC steps", prod_mc_steps);
 for(int imc=0; imc<prod_mc_steps; imc++)
 {
  metropolis_step(imc);
  if(imc%mc_interval==0)
   gather_observable_stat(); 
 };
  
 char prefix[500];
 largeN_QFT_prefix(prefix);
 process_mc_stat(prefix, 1);
 process_observable_stat();
 
 double cc_new = 1.5*cc;
 double NN_new = 1.2*NN;
 logs_WriteWarning("\n Predicted <nA> for cc = %2.4E, NN = %2.4E:\t %2.4E\n",     cc,     NN, mean_nA_prediction(    NN,     cc, action_create_amplitude(NULL)));
 logs_WriteWarning("\n Predicted <nA> for cc = %2.4E, NN = %2.4E:\t %2.4E\n", cc_new, NN_new, mean_nA_prediction(NN_new, cc_new, action_create_amplitude(NULL)));
 
 free_observable_stat();
 free_metropolis();
 free_actions();
 free_lat_propagator(&P);
 
 return 0;
}
