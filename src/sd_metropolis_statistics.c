#include "sd_metropolis_statistics.h"

double                   maxnA    = 0; //Max. value of nA
double                     anA    = 0; //Expectation value of nA[x] - it is necessary in order to restore the correct weight
double                     dnA    = 0; //Expectation value of its square - necessary to estimate the error!
double                   msign    = 0.0; //Expectation value of the sign - also necessary for recovering the correct normalization
int                        aac    = 0;
double                     ans    = 0.0;
int                        nmc    = 0;
int*            action_counter    = NULL;

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
 SAFE_MALLOC(action_counter, int, action_collection_size);
 for(int i=0; i<action_collection_size; i++)
  action_counter[i] = 0;
 //These are the variables which are only set by process_mc_stat
 acceptance_rate      = 0.0;
 mean_recursion_depth = 0.0;
 mean_nA              = 0.0;
 err_nA               = 0.0;
 mean_sign            = 0.0;
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
 nmc   ++;
}

void process_mc_stat(const char* prefix)
{
 //Summarizing the post-run properties of the MC process
 acceptance_rate      = (double)aac/(double)nmc;
 mean_recursion_depth =         ans/(double)nmc;
 mean_nA              =         anA/(double)nmc;
 err_nA               = sqrt((dnA/(double)nmc - SQR(mean_nA))/(double)(nmc-1));
 mean_sign            =       msign/(double)nmc;
 
 logs_Write(0, "\nSTATISTICS ON THE MC PROCESS (over %i steps): ",       nmc);
 logs_WriteParameter(         "Acceptance rate",    "%2.4lf",            acceptance_rate);
 logs_WriteParameter(    "Mean recursion depth",    "%2.4lf",            mean_recursion_depth);
 logs_WriteParameter( "Mean A rescaling factor",    "%2.4lf +/- %2.4lf", mean_nA, err_nA);
 logs_WriteParameter(  "Max A rescaling factor",    "%2.4lf",            maxnA);
 logs_WriteParameter(       "Mean config. sign",    "%2.4lf",            mean_sign);
 logs_Write(0, "\n");
 
 int action_counter_total = 0;
 for(int i=0; i<action_collection_size; i++)
  action_counter_total += action_counter[i];
 
 char* act_stat_str = NULL;
 sprintf_append(&act_stat_str, "%s ", prefix);
 logs_Write(0, "\tFACTUAL PROBABILITIES OF ACTIONS (over %i calls in %i mc steps): ",  action_counter_total, nmc);
 for(int i=0; i<action_collection_size; i++)
 {
  double act_prob = (double)(action_counter[i])/(double)action_counter_total;
  logs_Write(0, "   %30s [id = %i]: \t %2.4E \t (%02i%% of all actions, %i calls)", action_collection_name[i], i, act_prob, (int)(round(100.0*act_prob)), action_counter[i]);
  sprintf_append(&act_stat_str, "%16s %2.4E ", action_collection_name[i], act_prob);
 };
 sprintf_append(&act_stat_str, "\n");
 if(action_stat_file!=NULL)
  safe_append_str_to_file(action_stat_file, act_stat_str);
 SAFE_FREE(act_stat_str); 
   
 logs_Write(0, "\n");
 //Saving the statistical characteristics of the MC process
 if(mc_stat_file!=NULL)
 {
  int res = safe_append_to_file(mc_stat_file, "%s %2.4E %2.4E %2.4E %2.4E %2.4E %2.4E\n", prefix, acceptance_rate, mean_recursion_depth, mean_nA, err_nA, maxnA, mean_sign);
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

