#include "largeN_QFT_parameters.h"

//Parameters of a generic SD equations for practically any large-N QFT
double   lambda                   = 0.05;     //tHooft coupling constant
double   meff_sq                  = 0.0;      //Square of the effective mass
double   cc                       = 1.0;      //Rescaling of observables according to the number of fields in the correlator
double   NN                       = 1.0;      //Overall rescaling of observables
int      DIM                      = 1;        //Space-time dimensionality
int      LT                       = 2;        //Temporal size of the system
int      LS                       = 2;        //Spatial size of the system 
//Parameters of the statistical analysis of the MC data
int      max_stack_nel            = 10000;    //Maximal number of elements in the stack characterizing the system state
int      max_history_nel          = 10000;    //Maximal number of elements in the stack containing the history of momenta contractions
int      max_correlator_order     = 5;        //Maximal correlator order to trace
int      min_observables_order    = 0;        //Minimal order of observables which are included into statistics in some formal expansion (e.g. SC/WC expansion)
int      max_observables_order    = INT_MAX;  //Correspondingly, maximal order
//Output files
char*    observables_file         = NULL;     //File for the expectation values of the correlators
//Parameter tuning parameters
int      param_auto_tuning        = 1;        //Automatic tuning of transition amplitudes so that nAs are minimized
double   param_tuning_accuracy    = 0.000001; //Accuracy of parameter auto-tuning
int      param_tuning_max_iter    = 1000;     //Max. allowed number of iterations in param auto-tuning
//In debug mode, we can also check the stack consistency at every step
int      check_stack              = 1;


void print_largeN_QFT_parameters()
{
 logs_Write(0, "\tPARAMETERS OF A GENERIC SIMULATION OF A LARGE-N QFT");
 logs_WriteParameter(0,                                   "DIM",       "%i",      DIM);
 logs_WriteParameter(0,                                    "LT",       "%i",      LT);
 logs_WriteParameter(0,                                    "LS",       "%i",      LS);
 logs_WriteParameter(0,                                "lambda",    "%2.4E",      lambda);
 logs_WriteParameter(0,          "Square of the effective mass",    "%2.4E",      meff_sq);
 logs_WriteParameter(0,                                    "cc", "%2.4E %s",      cc, (param_auto_tuning? "(Automatically tuned)" : ""));
 logs_WriteParameter(0,                                    "NN", "%2.4E %s",      NN, (param_auto_tuning? "(Automatically tuned)" : ""));
 if(param_auto_tuning)
 {
 logs_WriteParameter(0,     "Accuracy of parameter auto-tuning",    "%2.4E",      param_tuning_accuracy);
 logs_WriteParameter(0,        "Max. iterations of auto-tuning",       "%i",      param_tuning_max_iter);
 };
 logs_WriteParameter(0,               "Max.correlator to trace",       "%i",      max_correlator_order);
 logs_WriteParameter(0,                 "Min.observables order",       "%i",      min_observables_order);
 logs_WriteParameter(0,                 "Max.observables order",       "%i",      max_observables_order);
 logs_WriteParameter(0,                         "max_stack_nel",       "%i",      max_stack_nel);
 logs_WriteParameter(0,                       "max_history_nel",       "%i",      max_history_nel);
 logs_WriteParameter(0, "Check stack consistency at every step",       "%s",      (check_stack? "YES" : "NO"));
 if(observables_file!=NULL)
 logs_WriteParameter(0,                      "Observables file",       "%s",      observables_file);
}

void largeN_QFT_prefix(char* prefix) //Prints lambda, cc, NN, LT, LS to prefix
{
 sprintf(prefix, "%2.4E %2.4E %2.4E %2.4E %i %i ", lambda, meff_sq, cc, NN, LT, LS);    
}

//The functions below implement auto-check of the minimization 
void cc_NN_vicinity(double epsilon, double* data)
{
 double cc0 = cc;
 double NN0 = NN;
 
 cc = cc0*(1.0 - epsilon);
 data[0] = f_max_ampl_sum();
 cc = cc0*(1.0 + epsilon);
 data[1] = f_max_ampl_sum();
 cc = cc0;
 
 NN = NN0*(1.0 - epsilon);
 data[2] = f_max_ampl_sum();
 NN = NN0*(1.0 + epsilon);
 data[3] = f_max_ampl_sum();
 NN = NN0; 
}

int  check_cc_NN_minimum(double tol)
{
 double vdata[4], S0;
 int res = 1, i;
 S0 = f_max_ampl_sum();
 cc_NN_vicinity(tol, vdata);
 for(i=0; i<4; i++)
  res = res && (S0<vdata[i]); 
 if(res) 
  logs_Write(0, "Max. amplitude sum has minimal value %2.4E at cc = %2.4E, NN = %2.4E", S0, cc, NN);
 else
  logs_WriteWarning("Max. amplitude sum value %2.4E at cc = %2.4E, NN = %2.4E seems not to be the minimum...", S0, cc, NN); 
 return res; 
}

int  find_cc_NN_minimum(double tol, double* min_val)
{
 double vdata[4], S0, my_min_val = 0.0;
 int res = 0, i, min_i, step_count;
 double epsilon = 0.1;

 while(epsilon>MIN(tol,0.01))
 {
  logs_Write(1, "Searching for the minimal values of cc and NN with tolerance %2.4E, starting with cc = %2.4E, NN = %2.4E", epsilon, cc, NN);
  res = 0;
  step_count = 0;
  while(!res && step_count<param_tuning_max_iter)
  {
   S0 = f_max_ampl_sum();
   cc_NN_vicinity(epsilon, vdata);
   my_min_val = vdata[0];
   min_i      = 0;
   for(i=1; i<4; i++)
    if(vdata[i]<my_min_val)
    {
     min_i      = i;
     my_min_val = vdata[i];
    };
   res = (my_min_val>S0);
   
   if(!res && min_i==0)
    cc *= (1.0 - epsilon);
   if(!res && min_i==1)
    cc *= (1.0 + epsilon);
   if(!res && min_i==2)
    NN *= (1.0 - epsilon);
   if(!res && min_i==3)
    NN *= (1.0 + epsilon);
   step_count ++;
   logs_Write(2, "Auto-minimizer step %03i: cc = %2.4E, NN = %2.4E, S = %2.4E", step_count, cc, NN, my_min_val);
  };
  logs_Write(1, "Optimal values with precision %2.4E: cc = %2.4E, NN = %2.4E", epsilon, cc, NN);
  epsilon *= 0.1;
 };
 
 logs_Write(0, "Minimal values of cc and NN (with precision %2.4E): \t cc = %2.4E, NN = %2.4E, max. amplitude sum = %2.4E", MIN(tol, 0.01), cc, NN, my_min_val);
 if(min_val!=NULL)
  *min_val = my_min_val;

 return res;
}

void cc_NN_vicinity_metropolis(double epsilon, int tune_mc_steps, double* v, double* dv, double* ccs, double* NNs, double* ms)
{
 double cc0 = cc;
 double NN0 = NN;
 double* sparam[2]  = { &cc,  &NN};
 double* sparam0[2] = {&cc0, &NN0};
 
 for(int i=0; i<4; i++)
 {
  (*(sparam[i/2])) = (*(sparam0[i/2]))*(1.0 + (double)(2*(i%2) - 1)*epsilon);
  init_metropolis();
  for(int imc=0; imc<tune_mc_steps; imc++)
   metropolis_step();
  //And now we should get the estimate of nA and its error
  //We should be able to tell the minimal value within the error 
  process_mc_stat("", 0);
  ccs[i] = cc;
  NNs[i] = NN;
    v[i] = mean_nA;
   dv[i] = err_nA;
   ms[i] = mean_sign;
 };
 
 cc = cc0;
 NN = NN0;
 init_metropolis();
}

//Tunes cc and NN using a real MC process in such a way that <nA> is minimized...
double   tune_cc_NN_minimum(double tol, int tune_mc_steps)
{
 double v0, dv0, v[4], dv[4], ccs[4], NNs[4], ms[4], ret_val = 0.0;
 //First find the naive values of cc and NN - just to start with...
 
 logs_Write(-1, "Start auto-tuning of cc and NN with cc = %2.4E, NN = %2.4E\n", cc, NN);
 
 init_metropolis();
 for(int imc=0; imc<tune_mc_steps; imc++)
  metropolis_step();
 process_mc_stat("", 0);
 v0  = mean_nA;
 dv0 =  err_nA;
 
 int step_count = 0, res = 0;
 while(!res && step_count<param_tuning_max_iter)
 {
  cc_NN_vicinity_metropolis(0.05, tune_mc_steps, v, dv, ccs, NNs, ms);
  double minv = v[0];
  int    mini =    0;
  logs_Write(-1, " Search step %i of %i: ", step_count, param_tuning_max_iter);
  for(int i=0; i<4; i++)
  {
   logs_Write(-1, "\t At cc = %2.4E,\t NN = %2.4E,\t <nA> = %2.4E +/- %2.4E,\t <sign> = %+2.4E", ccs[i], NNs[i], v[i], dv[i], ms[i]);
   if(v[i]<minv)
   {
    minv = v[i];
    mini = i; 
   };
  }; 
  if(minv<v0)
  {
   cc = ccs[mini];
   NN = NNs[mini];
   v0 = minv;
   dv0 = dv[mini];
   ret_val = minv;
   logs_Write(-1, "  => Changed to cc = %2.4E, NN = %2.4E, <nA> = %2.4E, <sign> = %2.4E\n", cc, NN, minv, ms[mini]);
  }
  else
  {
   logs_Write(-1, " => Stopping the search at cc = %2.4E, NN = %2.4E, <nA> = %2.4E +/- %2.4E\n", cc, NN, v0, dv0);
   res = 1; //Stop the search
  };
  step_count ++; 
 };
 
 init_metropolis(); 
 return ret_val;
}
