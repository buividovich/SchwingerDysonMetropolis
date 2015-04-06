#include "sd_metropolis_parameters.h"

//Parameters of the stochastic process itself
int      max_recursion_depth      = 10000;   //Max. possible depth of the sequence of states
double   p_plus                   = 0.5;     //Probability of the "Forward" move
int      p_plus_tuning            = 1;       //If 1 (default), p_plus will be automatically tuned to maximize acceptance
int      p_plus_tuning_interval   = 10000;   //p_plus will be updated in this interval
int      prod_mc_steps            = 10;      //Number of MC steps for production runs
int      therm_mc_steps           = 0;       //Number of MC steps for thermalization
int      mc_interval              = 1;       //Interval between the successive presumably uncorrelated measurements
int      mc_reporting_interval    = 1000;    //The interval at which the current state of the MC process is reported
char*    mc_stat_file             = NULL;    //File for the quantities characterizing the MC process itself
char*    action_stat_file         = NULL;    //Statistics on different actions
char*    ns_history_file          = NULL;    //File for saving the MC history of ns
double   io_sleep_time            = 0.0;   //Time to sleep during append_str_to_file
int      io_write_attempts        = 1;      //Attempts to write to the file if it is locked

void print_metropolis_parameters()
{
 if(io_sleep_time > 0.0 || io_write_attempts>1)
 {    
  logs_Write(0, "\tIO PARAMETERS");
  logs_WriteParameter(0,                "I/O sleeping interval", "%2.4lf", io_sleep_time);
  logs_WriteParameter(0,                   "I/O write attempts",     "%i", io_write_attempts);
 }; 
 logs_Write(0, "\tPARAMETERS OF THE METROPOLIS ALGORITHM");
 logs_WriteParameter(0,                  "Max. recursion depth",     "%i",   max_recursion_depth);
 logs_WriteParameter(0,           "Probability of forward move", "%2.4lf",   p_plus);
 if(p_plus_tuning)
  logs_WriteParameter(0,               "p_plus tuning interval",     "%i",   p_plus_tuning_interval);
 logs_WriteParameter(0,     "Number of MC steps for production",     "%i (%i%% of INT_MAX)",    prod_mc_steps, (int)(round(100.0*(double)prod_mc_steps/(double)INT_MAX)) );
 logs_WriteParameter(0, "Number of MC steps for thermalization",     "%i (%i%% of INT_MAX)",   therm_mc_steps, (int)(round(100.0*(double)therm_mc_steps/(double)INT_MAX)) );
 logs_WriteParameter(0,         "Interval between measurements",     "%i",   mc_interval);
 logs_WriteParameter(0,           "MC state reporting interval",     "%i",   mc_reporting_interval);
 logs_WriteParameter(0,                      "Logs noise level",     "%i",   logs_noise_level);
 if(mc_stat_file!=NULL)
  logs_WriteParameter(0,                   "MC statistics file",     "%s",   mc_stat_file);
 if(action_stat_file!=NULL)
  logs_WriteParameter(0,               "Action statistics file",     "%s",   action_stat_file); 
 if(ns_history_file!=NULL)
  logs_WriteParameter(0,         "Recursion depth history file",     "%s",   ns_history_file);
}

