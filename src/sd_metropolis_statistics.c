#include "sd_metropolis_statistics.h"

double                   maxnA    = 0; //Max. value of nA
double                     anA    = 0; //Expectation value of nA[x] - it is necessary in order to restore the correct weight
double                     dnA    = 0; //Expectation value of its square - necessary to estimate the error!
double                   msign    = 0.0; //Expectation value of the sign - also necessary for recovering the correct normalization
int                        aac    = 0;
double                     ans    = 0.0;
int                        nmc    = 0;
int*            action_counter    = NULL;

//Variables for tuning of p_plus
double  ppl_trial_min           = 0.00;
double  ppl_trial_max           = 1.00;
double  dppl                    = 0.111111;
int     n_ppl_trials            = 10;
double* acc_estimates           = NULL;
int     pplus_tuning_ndata      = 0;

//If control_max_ampl_sum=1 and max_ampl_sum is set to some nonzero value, 
//an error message is generated everytime nA exceeds max_ampl_sum, 
//this is useful for debugging
int       control_max_ampl_sum     = 0; 
double            max_ampl_sum     = 0.0;
double        max_ampl_sum_tol     = 0.0;

//These are the variables which are only set by process_mc_stat
double acceptance_rate      = 0.0;
double mean_recursion_depth = 0.0;
double mean_nA              = 0.0;
double err_nA               = 0.0;
double mean_sign            = 0.0;

void init_metropolis_statistics()
{
 maxnA = 0.0;
 aac   = 0;
 ans   = 0.0;
 anA   = 0.0;
 dnA   = 0.0;
 msign = 0.0;
 ans   = 0.0;
 nmc   = 0;
 SAFE_MALLOC_IF_NULL(action_counter, int, action_collection_size);
 for(int i=0; i<action_collection_size; i++)
  action_counter[i] = 0;
  
 max_ampl_sum = f_max_ampl_sum(); 
 if(max_ampl_sum>0.0)
 {
  control_max_ampl_sum = 1;
  max_ampl_sum_tol     = 0.001*max_ampl_sum;                   
 }; 
 //These are the variables which are only set by process_mc_stat
 acceptance_rate      = 0.0;
 mean_recursion_depth = 0.0;
 mean_nA              = 0.0;
 err_nA               = 0.0;
 mean_sign            = 0.0;
 
 if(p_plus_tuning)
  init_pplus_tuning();
}

void free_metropolis_statistics()
{
 if(p_plus_tuning)
  free_pplus_tuning();    
}

void gather_mc_stat()
{
 if(control_max_ampl_sum && (nA[ns] > (max_ampl_sum + max_ampl_sum_tol)))
  logs_WriteError("nA[%i] = %2.4E > max_ampl_sum = %2.4E", ns, nA[ns], max_ampl_sum);
 ans   += (double)(ns+1);
 anA   += nA[ns];
 dnA   += SQR(nA[ns]);
 maxnA  = MAX(nA[ns], maxnA);
 msign += asign[ns];
 
 if(p_plus_tuning)
  collect_pplus_tuning_data();
 if(p_plus_tuning && nmc>0 && nmc%p_plus_tuning_interval==0)
  tune_pplus();  
  
 nmc   ++;
}

void process_mc_stat(const char* prefix, int save_to_files)
{
 //Summarizing the post-run properties of the MC process
 acceptance_rate      = (double)aac/(double)nmc;
 mean_recursion_depth =         ans/(double)nmc;
 mean_nA              =         anA/(double)nmc;
 err_nA               = sqrt((dnA/(double)nmc - SQR(mean_nA))/(double)(nmc-1));
 mean_sign            =       msign/(double)nmc;
 
 logs_Write(0, "\nSTATISTICS ON THE MC PROCESS (over %i steps): ",       nmc);
 logs_WriteParameter(0,         "Acceptance rate",    "%2.4lf",            acceptance_rate);
 logs_WriteParameter(0,    "Mean recursion depth",    "%2.4lf",            mean_recursion_depth);
 logs_WriteParameter(0, "Mean A rescaling factor",    "%2.4lf +/- %2.4lf", mean_nA, err_nA);
 logs_WriteParameter(0,  "Max A rescaling factor",    "%2.4lf",            maxnA);
 logs_WriteParameter(0,       "Mean config. sign",    "%2.4lf",            mean_sign);
 logs_Write(0, "\n");
 
 int action_counter_total = 0;
 for(int i=0; i<action_collection_size; i++)
  action_counter_total += action_counter[i];
 
 char* act_stat_str = NULL;
 sprintf_append(&act_stat_str, "%s ", prefix);
 logs_Write(0, "\tFACTUAL PROBABILITIES OF ACTIONS (over %i calls in %i mc steps (%i%% of all mc steps)): ",  action_counter_total, nmc, (int)round(100.0*(double)(action_counter_total)/(double)nmc) );
 for(int i=0; i<action_collection_size; i++)
 {
  double act_prob = (double)(action_counter[i])/(double)action_counter_total;
  logs_Write(0, "   %30s [id = %i]: \t %2.4E \t (%02i%% of all actions, %i calls)", action_collection_name[i], i, act_prob, (int)(round(100.0*act_prob)), action_counter[i]);
  sprintf_append(&act_stat_str, "%16s %2.4E ", action_collection_name[i], act_prob);
 };
 sprintf_append(&act_stat_str, "\n");
 if(save_to_files && action_stat_file!=NULL)
  safe_append_str_to_file(action_stat_file, act_stat_str, io_sleep_time, io_write_attempts);
 SAFE_FREE(act_stat_str); 
   
 logs_Write(0, "\n");
 //Saving the statistical characteristics of the MC process
 if(save_to_files && mc_stat_file!=NULL)
 {
  int res = safe_append_to_file(mc_stat_file, io_sleep_time, io_write_attempts, "%s %2.4E %2.4E %2.4E %2.4E %2.4E %2.4E\n", prefix, acceptance_rate, mean_recursion_depth, mean_nA, err_nA, maxnA, p_plus);
  if(res!=0)
   logs_WriteError("safe_append_to_file %s failed with code %i", mc_stat_file, res);
 }; 
}

void print_max_amplitudes()
{
 int i, adata = -1;
 double ampl, ampl_sum = 0.0;
 logs_Write(0, "\t MAXIMAL AMPLITUDES OF ELEMENTARY ACTIONS");
 for(i=0; i<action_collection_size; i++)
 {
  ampl = (action_collection_amplitude[i])(&adata);
  logs_Write(0, "   %20s [action_id = %i]:\t %+2.4E", action_collection_name[i], i, ampl);
  ampl_sum += fabs(ampl);
 };
 logs_Write(0, " TOTAL: %2.4E\n", ampl_sum);
}

//Max. sum of all amplitudes - for the NAIVE tuning of cc and NN
double f_max_ampl_sum()
{
 int adata = -1, iaction;
 double ampl_sum = 0.0;
 for(iaction=0; iaction<action_collection_size; iaction++)
  ampl_sum += fabs((action_collection_amplitude[iaction])(&adata) );
 return ampl_sum;
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

