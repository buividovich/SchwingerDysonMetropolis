#include "sd_metropolis_parameters.h"

//Parameters of the stochastic process itself
int      max_recursion_depth      = 10000;   //Max. possible depth of the sequence of states
double   p_plus                   = 0.5;     //Probability of the "Forward" move
int      p_plus_tuning            = 1;       //If 1 (default), p_plus will be automatically tuned to maximize acceptance
int      p_plus_tuning_interval   = 10000;   //p_plus will be updated in this interval
int      number_mc_steps          = 100000;  //Number of MC steps for production runs
int      mc_reporting_interval    = 40000;   //The interval at which the current state of the MC process is reported
int      exit_upon_overflow       = 0;       //Whether to stop the MC process once the history or the stack overflow is detected

void print_metropolis_parameters()
{
 logs_Write(0, "\tPARAMETERS OF THE METROPOLIS ALGORITHM");
 logs_WriteParameter(0,                  "Max. recursion depth",     "%i",   max_recursion_depth);
 logs_WriteParameter(0,           "Probability of forward move", "%2.4lf",   p_plus);
 if(p_plus_tuning)
  logs_WriteParameter(0,               "p_plus tuning interval",     "%i",   p_plus_tuning_interval); 
 logs_WriteParameter(0,    "Action upon state/history overflow",     "%s",   (exit_upon_overflow? "Stop MC process" : "Reset MC process") );
 logs_WriteParameter(0,                   "Number of MC steps ",     "%i (%i%% of INT_MAX)",    number_mc_steps, (int)(round(100.0*(double)number_mc_steps/(double)INT_MAX)) );
 logs_WriteParameter(0,           "MC state reporting interval",     "%i",   mc_reporting_interval);
 logs_WriteParameter(0,                      "Logs noise level",     "%i",   logs_noise_level);
}
