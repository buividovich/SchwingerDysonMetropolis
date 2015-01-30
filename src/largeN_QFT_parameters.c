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
//Output files
char*    observables_file         = NULL;     //File for the expectation values of the correlators
int      param_auto_tuning        = 1;        //Automatic tuning of transition amplitudes so that nAs are minimized
double   param_tuning_accuracy    = 0.000001; //Accuracy of parameter auto-tuning

void print_largeN_QFT_parameters()
{
 logs_Write(0, "\tPARAMETERS OF A GENERIC SIMULATION OF A LARGE-N QFT");
 logs_WriteParameter(                                   "DIM",       "%i",      DIM);
 logs_WriteParameter(                                    "LT",       "%i",      LT);
 logs_WriteParameter(                                    "LS",       "%i",      LS);
 logs_WriteParameter(                                "lambda",    "%2.4E",      lambda);
 if(fabs(meff_sq)>0.0)
 logs_WriteParameter(          "Square of the effective mass",    "%2.4E",      meff_sq);
 logs_WriteParameter(                                    "cc", "%2.4E %s",      cc, (param_auto_tuning? "(Automatically tuned)" : ""));
 logs_WriteParameter(                                    "NN", "%2.4E %s",      NN, (param_auto_tuning? "(Automatically tuned)" : ""));
 if(param_auto_tuning)
  logs_WriteParameter(    "Accuracy of parameter auto-tuning",    "%2.4E",      param_tuning_accuracy);
 logs_WriteParameter(               "Max.correlator to trace",       "%i",      max_correlator_order);
 logs_WriteParameter(                         "max_stack_nel",       "%i",      max_stack_nel);
 logs_WriteParameter(                       "max_history_nel",       "%i",      max_history_nel);
 if(observables_file!=NULL)
  logs_WriteParameter(                     "Observables file",       "%s",      observables_file);
}

void largeN_QFT_prefix(char* prefix) //Prints lambda, cc, NN, LT, LS to prefix
{
 sprintf(prefix, "%2.4E %2.4E %2.4E %2.4E %i %i ", lambda, meff_sq, cc, NN, LT, LS);    
}

//The functions below implement auto-check of the minimization 

void cc_NN_vicinity(t_amplitude_sum S, double epsilon, double* data)
{
 double cc0 = cc;
 double NN0 = NN;
 
 cc = cc0*(1.0 - epsilon);
 data[0] = S();
 cc = cc0*(1.0 + epsilon);
 data[1] = S();
 cc = cc0;
 
 NN = NN0*(1.0 - epsilon);
 data[2] = S();
 NN = NN0*(1.0 + epsilon);
 data[3] = S();
 NN = NN0; 
}

int  check_cc_NN_minimum(t_amplitude_sum S, double tol)
{
 double vdata[4], S0;
 int res = 1, i;
 S0 = S();
 cc_NN_vicinity(S, tol, vdata);
 for(i=0; i<4; i++)
  res = res && (S0<vdata[i]); 
 if(res) 
  logs_Write(0, "Max. amplitude sum has minimal value %2.4E at cc = %2.4E, NN = %2.4E", S0, cc, NN);
 else
  logs_WriteWarning("Max. amplitude sum value %2.4E at cc = %2.4E, NN = %2.4E seems not to be the minimum...", S0, cc, NN); 
 return res; 
}

int  find_cc_NN_minimum(t_amplitude_sum S, double tol)
{
 double vdata[4], S0;
 int res = 0, i, min_i, step_count;
 double epsilon = 0.1, min_val;
 
 while(epsilon>MIN(tol,0.01))
 {
  logs_Write(1, "Searching for the minimal values of cc and NN with tolerance %2.4E, starting with cc = %2.4E, NN = %2.4E", epsilon, cc, NN);
  res = 0;
  step_count = 0;
  while(!res)
  {
   S0 = S();
   cc_NN_vicinity(S, epsilon, vdata);
   min_val = vdata[0];
   min_i   = 0;
   for(i=1; i<4; i++)
    if(vdata[i]<min_val)
    {
     min_i   = i;
     min_val = vdata[i];
    };
   res = (min_val>S0);
   if(!res && min_i==0)
    cc *= (1.0 - epsilon);
   if(!res && min_i==1)
    cc *= (1.0 + epsilon);
   if(!res && min_i==2)
    NN *= (1.0 - epsilon);
   if(!res && min_i==3)
    NN *= (1.0 + epsilon);
   step_count ++;
   logs_Write(2, "Auto-minimizer step %i: cc = %2.4E, NN = %2.4E", step_count, cc, NN);
  };
  logs_Write(1, "Optimal values with precision %2.4E: cc = %2.4E, NN = %2.4E", epsilon, cc, NN);
  epsilon *= 0.1;
 };
 
 logs_Write(0, "Minimal values of cc and NN (with precision %2.4E): \t cc = %2.4E, NN = %2.4E", MIN(tol, 0.01), cc, NN);
  
 return res;
}

