#include "parameters.h"

int      maxn                     = 10000;   //Follow the expansion up to this order
double   lambda                   = 0.05;    //Coupling constant in the Hermitian matrix model
double   cc                       = 1.0;     //Rescaling of observables
double   NN                       = 1.0;     //Rescaling of observables
double   pplus                    = 0.5;     //Probability of the "forward" step
int      nmc                      = 100000;  //Number of MC steps
int      maxg                     = 100;     //Maximal correlator order to trace
char*    observables_file         = NULL;    //File for the expectation values of the correlators
char*    mc_stat_file             = NULL;    //File for the quantities characterizing the MC process itself
char*    ns_hist_file             = NULL;    //File for saving the histogram of the recursion depth
char*    ns_history_file          = NULL;    //File for saving the MC history of ns
int      param_auto_tuning        = 0;       //Automatic tuning of NN and c so that the transition probabilities are minimized

static struct option long_options[] =
{
 {                     "maxn",  required_argument,                       NULL, 'A'},
 {                   "lambda",  required_argument,                       NULL, 'B'},
 {                       "cc",  required_argument,                       NULL, 'C'},
 {                       "NN",  required_argument,                       NULL, 'D'},
 {                    "pplus",  required_argument,                       NULL, 'E'},
 {                      "nmc",  required_argument,                       NULL, 'F'},
 {                     "maxg",  required_argument,                       NULL, 'G'},
 {         "observables-file",  required_argument,                       NULL, 'H'},
 {             "mc-stat-file",  required_argument,                       NULL, 'I'},
 {             "ns-hist-file",  required_argument,                       NULL, 'J'},
 {          "ns-history-file",  required_argument,                       NULL, 'K'},
 {        "param-auto-tuning",        no_argument,         &param_auto_tuning,   1},
 {                          0,                  0,                       NULL,   0}
};

static char short_option_list[] = "A:B:C:D:E:F:G:H:I:J:K:";

int parse_command_line_options(int argc, char **argv)
{
 int gc, option_index;
 //Compose the joint list of all short options
 
 while(1)
 {
  gc = getopt_long(argc, argv, short_option_list, long_options, &option_index);
  if(gc==-1)
   break;
  switch(gc)
  {
   case   0:
   break;
   case 'A':            //                 "maxn",
    SAFE_SSCANF_BREAK(optarg,  "%i", maxn);
    ASSERT(maxn<=1);
   break;
   case 'B':            //                 "lambda",
    SAFE_SSCANF_BREAK(optarg, "%lf", lambda);
   break; 
   case 'C':            //                 "cc",
    SAFE_SSCANF_BREAK(optarg, "%lf", cc);
   break;
   case 'D':            //                 "NN",
    SAFE_SSCANF_BREAK(optarg, "%lf", NN);
   break;
   case 'E':            //                 "pplus",
    SAFE_SSCANF_BREAK(optarg, "%lf", pplus);
   break;
   case 'F':            //                 "nmc",
    SAFE_SSCANF_BREAK(optarg, "%i", nmc);
    ASSERT(nmc<=1);
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
   case 'J':            //                 "ns_hist_file",
    COPY_FILE_NAME(optarg, ns_hist_file);
   break;
   case 'K':            //                 "ns_history_file",
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
 int memsize = (maxn*maxn)/2*sizeof(int);
 logs_Write(0, " SD-METROPOLIS SIMULATION OF A HERMITIAN MATRIX MODEL ");
 logs_Write(0, "  max. iteration depth...................%i",          maxn);
 logs_Write(0, "  sizeof(int)............................%i",        sizeof(int));
 logs_Write(0, "  Memory size required for the stack.....%i Mb",      memsize/(1024*1024));
 logs_Write(0, "  lambda.................................%2.4E",     lambda);
 logs_Write(0, "  cc.....................................%2.4E %s",      cc, (param_auto_tuning? "(Automatically tuned)" : ""));
 logs_Write(0, "  NN.....................................%2.4E %s",      NN, (param_auto_tuning? "(Automatically tuned)" : ""));
 logs_Write(0, "  pplus..................................%2.4lf",     pplus);
 logs_Write(0, "  Number of MC steps.....................%i",           nmc);
 logs_Write(0, "  Max.correlator to trace................%i",          maxg);
 if(observables_file!=NULL)
 logs_Write(0, "  Observables file.......................%s",     observables_file);
 if(mc_stat_file!=NULL)
 logs_Write(0, "  MC statistics file.....................%s",     mc_stat_file);
 if(ns_hist_file!=NULL)
 logs_Write(0, "  Recursion depth histogram file.........%s",     ns_hist_file);
 if(ns_history_file!=NULL)
 logs_Write(0, "  Recursion depth history file...........%s",     ns_history_file);
}

void printhelp()
{
 logs_WriteError("Wrong options!");
 exit(EXIT_FAILURE);
}



