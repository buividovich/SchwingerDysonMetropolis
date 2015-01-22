#include "sd_metropolis_parameters.h"

//Parameters of the stochastic process itself
int      max_recursion_depth      = 10000;   //Max. possible depth of the sequence of states
double   p_plus                   = 0.5;     //Probability of the "Forward" move
int      prod_mc_steps            = 10;      //Number of MC steps for production runs
int      therm_mc_steps           = 0;       //Number of MC steps for thermalization
int      mc_interval              = 1;       //Interval between the successive presumably uncorrelated measurements
int      mc_reporting_interval    = 1;       //The interval at which the current state of the MC process is reported
char*    mc_stat_file             = NULL;    //File for the quantities characterizing the MC process itself
char*    ns_history_file          = NULL;    //File for saving the MC history of ns

void print_metropolis_parameters()
{
 logs_Write(0, "\tPARAMETERS OF THE METROPOLIS ALGORITHM");
 logs_WriteParameter(                  "Max. recursion depth",     "%i",   max_recursion_depth);
 logs_WriteParameter(           "Probability of forward move", "%2.4lf",   p_plus);
 logs_WriteParameter(     "Number of MC steps for production",     "%i",   prod_mc_steps);
 logs_WriteParameter( "Number of MC steps for thermalization",     "%i",   therm_mc_steps);
 logs_WriteParameter(         "Interval between measurements",     "%i",   mc_interval);
 logs_WriteParameter(           "MC state reporting interval",     "%i",   mc_reporting_interval);
 logs_WriteParameter(                      "Logs noise level",     "%i",   logs_noise_level);
 if(mc_stat_file!=NULL)
  logs_WriteParameter(                   "MC statistics file",     "%s",   mc_stat_file);
 if(ns_history_file!=NULL)
  logs_WriteParameter(         "Recursion depth history file",     "%s",   ns_history_file);
}

