#include "largeN_QFT_parameters.h"

//Parameters of a generic SD equations for practically any large-N QFT
double   lambda                   = 0.05;    //tHooft coupling constant
double   cc                       = 1.0;     //Rescaling of observables according to the number of fields in the correlator
double   NN                       = 1.0;     //Overall rescaling of observables
int      DIM                      = 1;       //Space-time dimensionality
int      LT                       = 2;       //Temporal size of the system
int      LS                       = 2;       //Spatial size of the system 
//Parameters of the statistical analysis of the MC data
int      max_stack_nel            = 10000;   //Maximal number of elements in the stack characterizing the system state
int      max_history_nel          = 10000;   //Maximal number of elements in the stack containing the history of momenta contractions
int      max_correlator_order     = 5;       //Maximal correlator order to trace
//Output files
char*    observables_file         = NULL;    //File for the expectation values of the correlators
int      param_auto_tuning        = 1;       //Automatic tuning of transition amplitudes so that nAs are minimized
double   param_tuning_accuracy    = 0.00001; //Accuracy of parameter auto-tuning

void print_largeN_QFT_parameters()
{
 logs_Write(0, "\tPARAMETERS OF A GENERIC SIMULATION OF A LARGE-N QFT");
 logs_WriteParameter(                                   "DIM",       "%i",      DIM);
 logs_WriteParameter(                                    "LT",       "%i",      LT);
 logs_WriteParameter(                                    "LS",       "%i",      LS);
 logs_WriteParameter(                                "lambda",    "%2.4E",      lambda);
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
 sprintf(prefix, "%2.4E %2.4E %2.4E %i %i ", lambda, cc, NN, LT, LS);    
}
