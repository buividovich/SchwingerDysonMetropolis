#include "largeN_QFT_parameters.h"

//Parameters of a generic SD equations for practically any large-N QFT
double   lambda                   = 0.05;     //tHooft coupling constant
double   meff_sq                  = 0.0;      //Square of the effective mass
double   cc                       = 1.0;      //Rescaling of observables according to the number of fields in the correlator
double   NN                       = 1.0;      //Overall rescaling of observables
double   genus_A                  = 1.0;      //Constant A in recursion for cc[g]
double   genus_nu                 = 1.1;      //Constant nu in recursion for cc[g]
double   genus_B                  = 1.0;      //Constant B in recursion for NN[g]
double   genus_mu                 = 1.1;      //Constant mu in recursion for NN[g]
double   genus_f_exponent         = 1.0;      //This is an exponent before the gamma function in f[g] 
int      DIM                      = 1;        //Space-time dimensionality
int      LT                       = 2;        //Temporal size of the system
int      LS                       = 2;        //Spatial size of the system 
//Parameters of the statistical analysis of the MC data
int      max_stack_nel            = 10000;    //Maximal number of elements in the stack characterizing the system state
int      max_history_nel          = 10000;    //Maximal number of elements in the stack containing the history of momenta contractions
int      max_correlator_order     = 2;        //Maximal correlator order to trace
int      min_observables_order    = 0;        //Minimal order of observables which are included into statistics in some formal expansion (e.g. SC/WC expansion)
int      max_observables_order    = INT_MAX;  //Correspondingly, maximal order      
int      max_genus                = 0;        //Max order of 1/N^2 expansion
//Output files
char*    observables_file         = NULL;     //File for the expectation values of the correlators
char*    stack_stat_file          = NULL;     //File for histograms characterizing the stack usage
//Parameter tuning parameters
int      param_auto_tuning        = 1;        //Automatic tuning of transition amplitudes so that nAs are minimized
double   param_tuning_accuracy    = 0.000001; //Accuracy of parameter auto-tuning
int      param_tuning_max_iter    = 1000;     //Max. allowed number of iterations in param auto-tuning
//In debug mode, we can also check the stack consistency at every step
int      check_stack              = 1;
//Calculable parameters - to be filled by initialization procedures
double*  cc_genus                 = NULL;
double*  NN_genus                 = NULL;
double*   f_genus                 = NULL;
//Variable for genus counting
int         genus                 = 0;


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
 if(max_genus>0)
 {
  logs_WriteParameter(0,               "Max. observables genus",       "%i",      max_genus);
  logs_WriteParameter(0,      "Constant A in genus reweighting",    "%2.4E",      genus_A);
  logs_WriteParameter(0,     "Constant nu in genus reweighting",   "%2.2lf",      genus_nu);
  logs_WriteParameter(0,      "Constant B in genus reweighting",    "%2.4E",      genus_B);
  logs_WriteParameter(0,     "Constant mu in genus reweighting",   "%2.2lf",      genus_mu);
  logs_WriteParameter(0,      "f exponent in genus reweighting",    "%2.4E",      genus_f_exponent);
 };
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
 if(stack_stat_file!=NULL)
 logs_WriteParameter(0,                 "Stack statistics file",       "%s",      stack_stat_file);
}

void largeN_QFT_prefix(char* prefix) //Prints lambda, cc, NN, LT, LS to prefix
{
 sprintf(prefix, "%2.4E %2.4E %2.4E %2.4E %i %i ", lambda, meff_sq, cc, NN, LT, LS);    
}

double f_max_ampl_sum_genus()
{
 genus = 0;
 double res = f_max_ampl_sum();
 for(genus=1; genus<=max_genus; genus++)
  res = MAX(res, f_max_ampl_sum());
 return res; 
}

//The functions below implement auto-check of the minimization 
void param_vicinity(double epsilon, double** params, double* data, int nparams)
{
 for(int ipar=0; ipar<nparams; ipar++)
  for(int idir=0; idir<2; idir++)
  {
   double tmp = (params[ipar])[0];
   (params[ipar])[0] = tmp*(1.0 + epsilon*(double)(2*idir - 1));
   init_genus_constants(3);
   data[2*ipar + idir] = f_max_ampl_sum_genus();
   (params[ipar])[0] = tmp;
  };
 init_genus_constants(3); 
}

int  check_param_minimum(double epsilon, double** params, int nparams)
{
 DECLARE_AND_MALLOC(vdata, double, 2*nparams);
 max_ampl_sum = f_max_ampl_sum_genus();
 param_vicinity(epsilon, params, vdata, nparams);
 
 int res = 1;
 for(int i=0; i<2*nparams; i++)
  res = res && (max_ampl_sum<vdata[i]); 
  
 if(res) 
  logs_Write(0, "Max. amplitude sum has minimal value %2.4E at the current set of parameters", max_ampl_sum);
 else
  logs_WriteWarning("Max. amplitude sum value %2.4E seems not to be the minimum...", max_ampl_sum);
 SAFE_FREE(vdata);
 return res;
}

int  find_param_minimum(double tol, double** params, double* min_val, int nparams)
{
 DECLARE_AND_MALLOC(vdata, double, 2*nparams);
     
 double S0, my_min_val = 0.0;
 int res = 0, i, min_i, step_count;
 double epsilon = 0.1;

 //TODO: better way of parameter output?
 while(epsilon>MIN(tol,0.01))
 {
  logs_Write(1, "Searching for the minimal values of parameters with tolerance %2.4E", epsilon);
  res = 0;
  step_count = 0;
  init_genus_constants(3);
  while(!res && step_count<param_tuning_max_iter)
  {
   S0 = f_max_ampl_sum_genus();
   param_vicinity(epsilon, params, vdata, nparams);
   my_min_val = vdata[0];
   min_i      = 0;
   for(i=1; i<2*nparams; i++)
    if(vdata[i]<my_min_val)
    {
     min_i      = i;
     my_min_val = vdata[i];
    };
   res = (my_min_val>S0);
   
   int iparam = min_i/2;
   int idir   = min_i%2;
   
   if(!res)
   {
    (params[iparam])[0] *= (1.0 + epsilon*(double)(2*idir - 1));
    init_genus_constants(3);
   }; 
      
   step_count ++;
   logs_Write(2, "Auto-minimizer step %03i: cc = %2.4E, NN = %2.4E, S = %2.4E", step_count, cc, NN, my_min_val);
  };
  logs_Write(1, "Optimal values with precision %2.4E: cc = %2.4E, NN = %2.4E", epsilon, cc, NN);
  epsilon *= 0.1;
 };
 
 init_genus_constants(3);

 logs_Write(0, "Minimal values of cc and NN (with precision %2.4E): \t cc = %2.4E, NN = %2.4E, max. amplitude sum = %2.4E", MIN(tol, 0.01), cc, NN, my_min_val);
 if(min_val!=NULL)
  *min_val = my_min_val;

 SAFE_FREE(vdata);
 return res;
}

double n2an_sup(double a)
{
 RETURN_IF_FALSE(a<1.0, 0.0);
 double rp = 0.0;
 double rn = 3.0*a;
 int n = 2;
 while(rn>=rp)
 {
  rp = rn;
  rn = 0.5*(double)((n+1)*(n+2))*pow(a, (double)n);
  n++;            
 };
 return rp;  
}


void init_genus_constants(int noise_level)
{
 SAFE_MALLOC_IF_NULL(cc_genus, double, max_genus+2);
 SAFE_MALLOC_IF_NULL(NN_genus, double, max_genus+2);
 SAFE_MALLOC_IF_NULL( f_genus, double, max_genus+2);

 cc_genus[0] = cc;
 NN_genus[0] = NN;
  f_genus[0] = 1.0;
  
 logs_Write(noise_level, "\t g=%02i,\t cc = %2.4E, f = %2.4E, sup = %2.4E", 0, cc_genus[0], f_genus[0], 0.0);
 
 if(max_genus>0)
 {
  logs_Write(noise_level, "\nReweighting coefficients for the genus expansion");
  for(int g=1; g<=max_genus; g++)
  {
   cc_genus[g] = cc_genus[g-1]*(1.0 + genus_A/pow((double)g, genus_nu));
   NN_genus[g] = NN_genus[g-1]*(1.0 + genus_B/pow((double)g, genus_mu));
  
   double sup = n2an_sup(cc_genus[g-1]/cc_genus[g]);
  
   f_genus[g] = genus_f_exponent*f_genus[g-1]*sup;
   logs_Write(noise_level, "\t g=%02i,\t cc = %2.4E, f = %2.4E, sup = %2.4E", g, cc_genus[g], f_genus[g], sup);
  };
  logs_Write(noise_level, " ");
 };
}

void free_genus_constants()
{
 SAFE_FREE(cc_genus);
 SAFE_FREE(NN_genus);
 SAFE_FREE( f_genus);
}

void free_largeN_QFT_parameters()
{
 SAFE_FREE(observables_file);
 SAFE_FREE(stack_stat_file);
}
