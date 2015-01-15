#include "hmm_parameters.h"

//Parameters of SD equations
double   lambda                   = 0.05;    //Coupling constant in the Hermitian matrix model
double   cc                       = 1.0;     //Rescaling of observables
double   NN                       = 1.0;     //Rescaling of observables
int      param_auto_tuning        = 0;       //Automatic tuning of NN and c so that the transition probabilities are minimized
//Parameters of the statistical analysis of the MC data
int      prod_mc_steps            = 10000;   //Number of MC steps for production runs
int      therm_mc_steps           = 0;       //Number of MC steps for thermalization
int      mc_interval              = 1;       //Interval between the successive presumably uncorrelated measurements
int      maxg                     = 100;     //Maximal correlator order to trace
//Output files
char*    observables_file         = NULL;    //File for the expectation values of the correlators
char*    mc_stat_file             = NULL;    //File for the quantities characterizing the MC process itself
char*    ns_history_file          = NULL;    //File for saving the MC history of ns

static struct option long_options[] =
{
 METROPOLIS_LONG_OPTIONS,
 {                   "lambda",  required_argument,                       NULL, 'A'},
 {                       "cc",  required_argument,                       NULL, 'B'},
 {                       "NN",  required_argument,                       NULL, 'C'},
 {            "prod_mc_steps",  required_argument,                       NULL, 'D'},
 {           "therm_mc_steps",  required_argument,                       NULL, 'E'},
 {              "mc_interval",  required_argument,                       NULL, 'F'},
 {                     "maxg",  required_argument,                       NULL, 'G'},
 {         "observables-file",  required_argument,                       NULL, 'H'},
 {             "mc-stat-file",  required_argument,                       NULL, 'I'},
 {          "ns-history-file",  required_argument,                       NULL, 'J'},
 {        "param-auto-tuning",        no_argument,         &param_auto_tuning,   1},
 {                          0,                  0,                       NULL,   0}
};

static char short_option_list[] = "A:B:C:D:E:F:G:H:I:J:";

int parse_command_line_options(int argc, char **argv)
{
 int gc, option_index;
 //Compose the joint list of all short options
 int short_option_list_length = strlen(metropolis_short_option_list) + strlen(short_option_list) + 3;
 char* short_option_list = NULL;
 SAFE_MALLOC(short_option_list, char, short_option_list_length);
 sprintf(short_option_list, "%s%s", metropolis_short_option_list, short_option_list); 
 
 while(1)
 {
  gc = getopt_long(argc, argv, short_option_list, long_options, &option_index);
  if(gc==-1)
   break;
  switch(gc)
  {
   PARSE_METROPOLIS_OPTIONS;
   case   0:
   break;
   case 'A':            //                 "lambda",
    SAFE_SSCANF_BREAK(optarg, "%lf", lambda);
   break; 
   case 'B':            //                 "cc",
    SAFE_SSCANF_BREAK(optarg, "%lf", cc);
   break;
   case 'C':            //                 "NN",
    SAFE_SSCANF_BREAK(optarg, "%lf", NN);
   break;
   case 'D':            //                 "prod_mc_steps",
    SAFE_SSCANF_BREAK(optarg, "%i", prod_mc_steps);
    ASSERT(prod_mc_steps<0);
   break;
   case 'E':            //                 "therm_mc_steps",
    SAFE_SSCANF_BREAK(optarg, "%i", therm_mc_steps);
    ASSERT(therm_mc_steps<0);
   break;
   case 'F':            //                 "mc_interval",
    SAFE_SSCANF_BREAK(optarg, "%i", mc_interval);
    ASSERT(mc_interval<1);
   break;
   case 'G':            //                 "maxg",
    SAFE_SSCANF_BREAK(optarg, "%i", maxg);
    ASSERT(maxg<0);
   break;
   case 'H':            //                 "observables_file",
    COPY_FILE_NAME(optarg, observables_file);
   break;
   case 'I':            //                 "mc_stat_file",
    COPY_FILE_NAME(optarg, mc_stat_file);
   break;
   case 'J':            //                 "ns_history_file",
    COPY_FILE_NAME(optarg, ns_history_file);
   break;
   case '?':
	printhelp();
   break;	
   default:
    printhelp();       
   break;
  }; 
 }; 
 return 0;
}

void init_parameters()
{
 if(param_auto_tuning)
 {
  cc = 2.0/sqrt(fabs(lambda));
  NN = cc;
 };    
}

void print_parameters()
{
 print_metropolis_parameters();
 logs_WriteParameter(                                "lambda",    "%2.4E",      lambda);
 logs_WriteParameter(                                    "cc", "%2.4E %s",      cc, (param_auto_tuning? "(Automatically tuned)" : ""));
 logs_WriteParameter(                                    "NN", "%2.4E %s",      NN, (param_auto_tuning? "(Automatically tuned)" : ""));
 logs_WriteParameter(     "Number of MC steps for production",       "%i",      prod_mc_steps);
 logs_WriteParameter( "Number of MC steps for thermalization",       "%i",      therm_mc_steps);
 logs_WriteParameter(         "Interval between measurements",       "%i",      mc_interval);
 logs_WriteParameter(               "Max.correlator to trace",       "%i",      maxg);
 if(observables_file!=NULL)
  logs_WriteParameter(                     "Observables file",       "%s",      observables_file);
 if(mc_stat_file!=NULL)
  logs_WriteParameter(                   "MC statistics file",       "%s",      mc_stat_file);
 if(ns_history_file!=NULL)
  logs_WriteParameter(         "Recursion depth history file",       "%s",      ns_history_file);
}

void printhelp()
{
 logs_WriteError("Wrong options!");
 exit(EXIT_FAILURE);
}


