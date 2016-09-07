#include "sd_metropolis_statistics.h"

double                   maxnA    = 0.0; //Max. value of nA
double                     anA    = 0.0; //Expectation value of nA[x] - it is necessary in order to restore the correct weight
double                     dnA    = 0.0; //Expectation value of its square - necessary to estimate the error!
double                   msign    = 0.0; //Expectation value of the sign - also necessary for recovering the correct normalization
int                        aac    = 0;
double                     ans    = 0.0;
int                        nmc    = 0;
int*            action_counter    = NULL;
int*                ns_history    = NULL; //MC history of sequence lengths

//Variables for tuning of p_plus
double  ppl_trial_min             = 0.00;
double  ppl_trial_max             = 1.00;
double  dppl                      = 0.111111;
int     n_ppl_trials              = 10;
double* acc_estimates             = NULL;
int     pplus_tuning_ndata        = 0;

//Variables for estimating return times
int     prev_return_time          = 0;
int     n_returns                 = 0;
double  mean_return_time          = 0.0;

//These are the variables which are only set by process_mc_stat
double acceptance_rate            = 0.0;
double mean_recursion_depth       = 0.0;
double mean_nA                    = 0.0;
double mean_sign                  = 0.0;

void init_metropolis_statistics()
{
 maxnA    =     0.0;
 aac      =     0.0;
 ans      =     0.0;
 anA      =     0.0;
 dnA      =     0.0;
 msign    =     0.0;
 ans      =     0.0;
 nmc      =     0;
 
 SAFE_MALLOC_IF_NULL(action_counter, int, action_collection_size);
 for(int i=0; i<action_collection_size; i++)
  action_counter[i] = 0;
 
 //These are the variables which are only set by process_mc_stat
 acceptance_rate      = 0.0;
 mean_recursion_depth = 0.0;
 mean_nA              = 0.0;
 mean_sign            = 0.0;
 
 //Variables for estimating return times
 prev_return_time     = 0;
 n_returns            = 0;
 mean_return_time     = 0.0;
 
 if(p_plus_tuning)
  init_pplus_tuning();
}

void free_metropolis_statistics()
{
 free_pplus_tuning();
 SAFE_FREE(action_counter);
}

void gather_mc_stat()
{
 ans   += (double)ns;
 anA   += nA[ns];
 dnA   += SQR(nA[ns]);
 
 maxnA  = MAX(nA[ns], maxnA);
 msign += asign[ns];
 
 if(p_plus_tuning)
  collect_pplus_tuning_data();
 if(p_plus_tuning && nmc>0 && nmc%p_plus_tuning_interval==0)
  tune_pplus();  
  
 //Updating the statistics of return time
 if(ns==0)
  n_returns ++;
  
 nmc   ++;
}

void process_mc_stat()
{
 //Summarizing the post-run properties of the MC process
 acceptance_rate      = (double)aac/(double)nmc;
 mean_recursion_depth =         ans/(double)nmc;
 mean_nA              =         anA/(double)nmc;
 mean_sign            =       msign/(double)nmc;
 mean_return_time     = (double)nmc/(double)n_returns;
 
 logs_Write(0, "\nSTATISTICS ON THE MC PROCESS (over %i steps): ",       nmc);
 logs_WriteParameter(0,                                   "Acceptance rate",  "%2.4lf",            acceptance_rate);
 logs_WriteParameter(0,                              "Mean recursion depth",  "%2.4lf",            mean_recursion_depth);
 logs_WriteParameter(0,                                           "Mean nA",  "%2.4lf",            mean_nA);
 logs_WriteParameter(0,                            "Max A rescaling factor",  "%2.4lf",            maxnA);
 logs_WriteParameter(0,                                 "Mean config. sign",  "%2.4lf",            mean_sign);
 logs_WriteParameter(0,                                  "Mean return time",  "%2.4lf",            mean_return_time);
 logs_Write(0, "\n");
 
 int action_counter_total = 0;
 for(int i=0; i<action_collection_size; i++)
  action_counter_total += action_counter[i];
 
 logs_Write(0, "\tFACTUAL PROBABILITIES OF ACTIONS (over %i calls in %i mc steps (%i%% of all mc steps)): ",  action_counter_total, nmc, (int)round(100.0*(double)(action_counter_total)/(double)nmc) );
 for(int i=0; i<action_collection_size; i++)
 {
  double act_prob = (double)(action_counter[i])/(double)action_counter_total;
  logs_Write(0, "   %30s [id = %i]: \t %2.4E \t (%02i%% of all actions, %i calls)", action_collection_name[i], i, act_prob, (int)(round(100.0*act_prob)), action_counter[i]);
 };
}

void init_pplus_tuning()
{
 ppl_trial_min           = 0.0;
 ppl_trial_max           = 1.00;
 n_ppl_trials            = 10;
 
 dppl = (ppl_trial_max - ppl_trial_min)/(double)(n_ppl_trials - 1);
 
 SAFE_MALLOC(acc_estimates, double, n_ppl_trials);
 for(int ippl=0; ippl<n_ppl_trials; ippl++)
  acc_estimates[ippl] = 0.0;
 pplus_tuning_ndata = 0; 
}

void collect_pplus_tuning_data()
{
 for(int ippl=0; ippl<n_ppl_trials; ippl++)
 {
  double ppl_trial = ppl_trial_min + (double)ippl*dppl;
  acc_estimates[ippl]  += MIN(ppl_trial, nA[ns]*(1.0 - ppl_trial));   
 };
 pplus_tuning_ndata ++;
}

void tune_pplus()
{
 double max_acceptance = 0.0;
 double ppl_max        = 0.0;  
 for(int ippl=0; ippl<n_ppl_trials; ippl++)
 {
  double ppl_trial = ppl_trial_min + (double)ippl*dppl;
  double acceptance = 2.0*acc_estimates[ippl]/(double)pplus_tuning_ndata + (1.0 - ppl_trial)*(1.0 - anA/(double)nmc);
  if(acceptance>max_acceptance)
  {
   max_acceptance = acceptance;
   ppl_max        = ppl_trial;
  };
  acc_estimates[ippl] = 0.0;
 };
 logs_Write(0, "\n Auto-tuning of p_plus at step %04i, acceptance so far %2.4lf: ", nmc, aac/(double)nmc);
 logs_WriteParameter(0, "Optimal value of p_plus", "%2.4lf (estimate over %i MC steps), expected acceptance %2.4lf", ppl_max, pplus_tuning_ndata, max_acceptance);
 if(ppl_max>0.0 && fabs((p_plus-ppl_max)/ppl_max)<0.02)
 {
  p_plus_tuning = 0;
  logs_WriteWarning("Automatic tuning of p_plus switched off, as the optimal value is already very close");
 }
 else
 {
  p_plus = ppl_max;
  ppl_trial_min = MAX(ppl_max - dppl, 0.0);
  ppl_trial_max = fabs(ppl_max + dppl);
  dppl = (ppl_trial_max - ppl_trial_min)/(double)(n_ppl_trials - 1);
  logs_WriteParameter(0, "New range of trial values of p_plus", "%2.4lf ... %2.4lf (in steps of %2.4E)", ppl_trial_min, ppl_trial_max, dppl);
 }; 
 logs_Write(0, "\n");
 pplus_tuning_ndata = 0; 
}

void free_pplus_tuning()
{
 SAFE_FREE(acc_estimates);
}

