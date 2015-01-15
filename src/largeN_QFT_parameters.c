#include "largeN_QFT_parameters.h"

//Parameters of a generic SD equations for practically any large-N QFT
double   lambda                   = 0.05;    //tHooft coupling constant
double   cc                       = 1.0;     //Rescaling of observables according to the number of fields in the correlator
double   NN                       = 1.0;     //Overall rescaling of observables 
//Parameters of the statistical analysis of the MC data
int      max_correlator_order     = 100;     //Maximal correlator order to trace
//Output files
char*    observables_file         = NULL;    //File for the expectation values of the correlators
int      param_auto_tuning        = 0;       //Automatic tuning of transition amplitudes so that nAs are minimized

void print_largeN_QFT_parameters()
{
 logs_Write(0, "\tPARAMETERS OF A GENERIC SIMULATION OF A LARGE-N QFT");
 logs_WriteParameter(                                "lambda",    "%2.4E",      lambda);
 logs_WriteParameter(                                    "cc", "%2.4E %s",      cc, (param_auto_tuning? "(Automatically tuned)" : ""));
 logs_WriteParameter(                                    "NN", "%2.4E %s",      NN, (param_auto_tuning? "(Automatically tuned)" : ""));
 logs_WriteParameter(               "Max.correlator to trace",       "%i",      max_correlator_order);
 if(observables_file!=NULL)
  logs_WriteParameter(                     "Observables file",       "%s",      observables_file);
}
