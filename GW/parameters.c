#include "parameters.h"

const char * const       mode_names[] = {"Strong-coupling expansion", "Weak-coupling expansion", "Matrix Lagrange multiplier"};
const char * const short_mode_names[] = {"sc", "wc", "lm"};

int      mode                     = 0;       //0 is for the strong-coupling expansion
int      maxn                     = 10000;   //Follow the expansion up to this order
double   lambda                   = 0.05;    //Coupling constant in the Hermitian matrix model
double   cc                       = 1.0;     //Rescaling of observables
double   NN                       = 1.0;     //Rescaling of observables
double   pplus                    = 0.5;     //Probability of the "forward" step
int      nmc                      = 100000;  //Number of MC steps
int      maxg                     = 20;      //Maximal correlator order to trace
char*    observables_file         = NULL;    //File for the expectation values of the correlators
char*    mc_stat_file             = NULL;    //File for the quantities characterizing the MC process itself
char*    ns_hist_file             = NULL;    //File for saving the histogram of the recursion depth
char*    ns_history_file          = NULL;    //File for saving the MC history of ns
int      param_auto_tuning        = 0;       //Automatic tuning of NN and c so that the transition probabilities are minimized
double   param_tuning_accuracy    = 0.00001; //Accuracy of parameter tuning 

double   alpha_wc                 = 0.0; //Parameter alpha of the weak-coupling expansion

static struct option long_options[] =
{
 {                     "mode",  required_argument,                       NULL, 'A'},
 {                     "maxn",  required_argument,                       NULL, 'B'},
 {                   "lambda",  required_argument,                       NULL, 'C'},
 {                       "cc",  required_argument,                       NULL, 'D'},
 {                       "NN",  required_argument,                       NULL, 'E'},
 {                    "pplus",  required_argument,                       NULL, 'F'},
 {                      "nmc",  required_argument,                       NULL, 'G'},
 {                     "maxg",  required_argument,                       NULL, 'H'},
 {         "observables-file",  required_argument,                       NULL, 'I'},
 {             "mc-stat-file",  required_argument,                       NULL, 'J'},
 {             "ns-hist-file",  required_argument,                       NULL, 'K'},
 {          "ns-history-file",  required_argument,                       NULL, 'L'},
 {        "param-auto-tuning",        no_argument,         &param_auto_tuning,   1},
 {                          0,                  0,                       NULL,   0}
};

static char short_option_list[] = "A:B:C:D:E:F:G:H:I:J:K:L:";

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
   case 'A':            //                 "mode",
    SAFE_SSCANF_BREAK(optarg,  "%i", mode);
    ASSERT(mode<0 || mode>2);
   break;
   case 'B':            //                 "maxn",
    SAFE_SSCANF_BREAK(optarg,  "%i", maxn);
    ASSERT(maxn<=1);
   break;
   case 'C':            //                 "lambda",
    SAFE_SSCANF_BREAK(optarg, "%lf", lambda);
   break; 
   case 'D':            //                 "cc",
    SAFE_SSCANF_BREAK(optarg, "%lf", cc);
   break;
   case 'E':            //                 "NN",
    SAFE_SSCANF_BREAK(optarg, "%lf", NN);
   break;
   case 'F':            //                 "pplus",
    SAFE_SSCANF_BREAK(optarg, "%lf", pplus);
   break;
   case 'G':            //                 "nmc",
    SAFE_SSCANF_BREAK(optarg, "%i", nmc);
    ASSERT(nmc<=1);
   break;
   case 'H':            //                 "maxg",
    SAFE_SSCANF_BREAK(optarg, "%i", maxg);
    ASSERT(maxg<0);
   break;
   case 'I':            //                 "observables_file",
    COPY_FILE_NAME(optarg, observables_file);
   break;
   case 'J':            //                 "mc_stat_file",
    COPY_FILE_NAME(optarg, mc_stat_file);
   break;
   case 'K':            //                 "ns_hist_file",
    COPY_FILE_NAME(optarg, ns_hist_file);
   break;
   case 'L':            //                 "ns_history_file",
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
 alpha_wc = fabs(lambda)/(2.0*fabs(lambda) + 8.0);
 if(param_auto_tuning)
 {
  if(mode==0)
   tune_parameters_sc();
  if(mode==1)
   tune_parameters_wc(); 
 }; 
}

const char ansgreen[] = "\x1b[1;32m";

void print_parameters()
{
 int memsize = (maxn*maxn)/2*sizeof(int);
 logs_Write(0, " SD-METROPOLIS SIMULATION OF THE GROSS-WITTEN MATRIX MODEL ");
 logs_Write(0, "  Mode...................................%s%s",        ansgreen,     mode_names[mode]);
 logs_Write(0, "  max. iteration depth...................%s%i",        ansgreen,     maxn);
 logs_Write(0, "  sizeof(int)............................%s%i",        ansgreen,      sizeof(int));
 logs_Write(0, "  Memory size required for the stack.....%s%i Mb",     ansgreen,   memsize/(1024*1024));
 logs_Write(0, "  lambda.................................%s%2.4E",     ansgreen,   lambda);
 logs_Write(0, "  alpha_wc...............................%s%2.4E",     ansgreen,   alpha_wc);
 logs_Write(0, "  cc.....................................%s%2.4E %s",  ansgreen,       cc, (param_auto_tuning? "(Automatically tuned)" : ""));
 logs_Write(0, "  NN.....................................%s%2.4E %s",  ansgreen,       NN, (param_auto_tuning? "(Automatically tuned)" : ""));
 
 if(mode==0)
 logs_Write(0, "  ptot...................................%s%2.4E %s",  ansgreen, ptot_sc(cc),(param_auto_tuning? "(Automatically tuned)" : ""));
 if(mode==1)
 logs_Write(0, "  ptot...................................%s%2.4E %s",  ansgreen, ptot_wc(cc),(param_auto_tuning? "(Automatically tuned)" : ""));
 
 logs_Write(0, "  pplus..................................%s%2.4lf",    ansgreen, pplus);
 logs_Write(0, "  Number of MC steps.....................%s%i",        ansgreen,   nmc);
 logs_Write(0, "  Max.correlator to trace................%s%i",        ansgreen,  maxg);
 if(observables_file!=NULL)
 logs_Write(0, "  Observables file.......................%s%s",        ansgreen, observables_file);
 if(mc_stat_file!=NULL)
 logs_Write(0, "  MC statistics file.....................%s%s",        ansgreen, mc_stat_file);
 if(ns_hist_file!=NULL)
 logs_Write(0, "  Recursion depth histogram file.........%s%s",        ansgreen, ns_hist_file);
 if(ns_history_file!=NULL)
 logs_Write(0, "  Recursion depth history file...........%s%s",        ansgreen, ns_history_file);
}

void printhelp()
{
 logs_WriteError("Wrong options!");
 exit(EXIT_FAILURE);
}

double ptot_sc(double c)
{
 return 2.0/sqrt(fabs(lambda)*c) + 1.0/(fabs(lambda)*c) + c/fabs(lambda);
}

double dptot_sc(double c)
{
 return -1.0/sqrt(fabs(lambda)*c*c*c) - 1.0/(fabs(lambda)*c*c) + 1.0/fabs(lambda);
}

int tune_parameters_sc()
{
 double c1 = 0.5, c2 = 20.0, cm, dpm;
 double dp1 = dptot_sc(c1);
 double dp2 = dptot_sc(c2);
 logs_Write(0, "Tuning parameters cc and NN in the strong-coupling regime...");
 while(fabs(c1-c2)>param_tuning_accuracy || fabs(dp1-dp2)>param_tuning_accuracy)
 {
  logs_Write(1, "c1 = %2.4lf, c2 = %2.4lf, dp1 = %2.4lf, dp2 = %2.4lf", c1, c2, dp1, dp2);
  if(dp1*dp2>0)
  {
   logs_WriteError("dp1*dp2>0 in tune_parameters...");
   return -1;       
  };
  cm = 0.5*(c1 + c2);
  dpm = dptot_sc(cm);
  if(dpm*dp1<0)
  {
   c2  = cm;
   dp2 = dpm;
  }
  else
  {
   c1  = cm;
   dp1 = dpm;
  };
 };
 cc = 0.5*(c1 + c2);
 NN = 1.0/sqrt(fabs(lambda)*cc);
 double ptot_check = 1.0/(fabs(lambda)*NN*cc) + NN + 1.0/(fabs(lambda)*cc) + cc/fabs(lambda);
 logs_Write(1, " => cc = %2.4lf, NN = %2.4lf, ptot = %2.4lf\n", cc, NN, ptot_check);
 return 0;
}

double ptot_wc(double c)
{
 return 4.0/c + 2.0*c*SQR(alpha_wc)/(lambda*(1.0 - c*alpha_wc))*(4.0 + lambda + 4.0/(1.0 - c*alpha_wc));
}

double dptot_wc(double c)
{
 return -4.0/SQR(c) + 2.0*SQR(alpha_wc)/(lambda*SQR(1.0 - alpha_wc*c))*(8.0/(1.0 - c*alpha_wc) + lambda);
}

int tune_parameters_wc()
{
 double c1 = 0.1, c2 = 20.0, cm, dpm;
 double dp1 = dptot_wc(c1);
 double dp2 = dptot_wc(c2);
 logs_Write(0, "Tuning parameters cc and NN in the weak-coupling regime...");
 while(fabs(c1-c2)>param_tuning_accuracy || fabs(dp1-dp2)>param_tuning_accuracy)
 {
  logs_Write(1, "c1 = %2.4lf, c2 = %2.4lf, dp1 = %2.4lf, dp2 = %2.4lf", c1, c2, dp1, dp2);
  if(dp1*dp2>0)
  {
   logs_WriteError("dp1*dp2>0 in tune_parameters...");
   return -1;       
  };
  cm = 0.5*(c1 + c2);
  dpm = dptot_wc(cm);
  if(dpm*dp1<0)
  {
   c2  = cm;
   dp2 = dpm;
  }
  else
  {
   c1  = cm;
   dp1 = dpm;
  };
 };
 cc = 0.5*(c1 + c2);
 NN = 1.0;
 double ptot_check = 1.0/(NN*cc) + NN/cc + 2.0/cc + 2.0*cc*SQR(alpha_wc)/(lambda*(1.0 - cc*alpha_wc))*(4.0 + lambda + 4.0/(1.0 - cc*alpha_wc));
 logs_Write(1, " => cc = %2.4lf, NN = %2.4lf, ptot = %2.4lf\n", cc, NN, ptot_check);
 return 0;
}



