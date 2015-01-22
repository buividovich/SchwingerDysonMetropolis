#include "sd_metropolis_statistics.h"

double                   maxnA    = 0; //Max. value of nA
double                     anA    = 0; //Expectation value of nA[x] - it is necessary in order to restore the correct weight
double                     dnA    = 0; //Expectation value of its square - necessary to estimate the error!
double                   msign    = 0.0; //Expectation value of the sign - also necessary for recovering the correct normalization
int                        aac    = 0;
double                     ans    = 0.0;
int                        nmc    = 0;

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
 //These are the variables which are only set by process_mc_stat
 acceptance_rate      = 0.0;
 mean_recursion_depth = 0.0;
 mean_nA              = 0.0;
 err_nA               = 0.0;
 mean_sign            = 0.0;
}

void gather_mc_stat(int accepted)
{
 ans   += (double)(ns+1);
 anA   += nA[ns];
 dnA   += SQR(nA[ns]);
 maxnA  = MAX(nA[ns], maxnA);
 msign += asign[ns]; 
 aac   += accepted;
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
 //Saving the statistical characteristics of the MC process
 if(mc_stat_file!=NULL)
 {
  FILE *f = fopen(mc_stat_file, "a");
  fprintf(f, "%s %2.4E %2.4E %2.4E %2.4E %2.4E %2.4E\n", prefix, acceptance_rate, mean_recursion_depth, mean_nA, err_nA, maxnA, mean_sign);
  fclose(f);
 }; 
}

