#include "parameters.h"

static struct option long_options[] =
{
 METROPOLIS_LONG_OPTIONS,
 LARGEN_QFT_LONG_OPTIONS,
 {                       0,                  0,                       NULL,   0}
};

static char my_short_option_list[] = "";

//Parameters which are calculated from the user input
double   alpha_wc                 = 0.0; //Parameter alpha of the weak-coupling expansion

int parse_command_line_options(int argc, char **argv)
{
 int gc, option_index;
 //Compose the joint list of all short options
 int short_option_list_length = 3;
 short_option_list_length += strlen(   metropolis_short_option_list);
 short_option_list_length += strlen(   largeN_QFT_short_option_list);
 short_option_list_length += strlen(           my_short_option_list);
 char* short_option_list = NULL;
 SAFE_MALLOC(short_option_list, char, short_option_list_length);
 sprintf(short_option_list, "%s",     my_short_option_list);
 strcat(short_option_list,    metropolis_short_option_list);
 strcat(short_option_list,    largeN_QFT_short_option_list);

 while(1)
 {
  gc = getopt_long(argc, argv, short_option_list, long_options, &option_index);
  if(gc==-1)
   break;
  switch(gc)
  {
   PARSE_METROPOLIS_OPTIONS;
   PARSE_LARGEN_QFT_OPTIONS;
   case   0:
   break;
   case '?':
	printhelp();
   break;	
   default:
    printhelp();       
   break;
  }; 
 }; 
 
 SAFE_FREE(short_option_list);
 return 0;
}

void init_parameters()
{
 alpha_wc = fabs(lambda)/(2.0*fabs(lambda) + 8.0);
 if(param_auto_tuning)
 {
  cc = 1.0; NN = 1.0;
  find_cc_NN_minimum(param_tuning_accuracy, &max_ampl_sum);
  if(max_ampl_sum>0.0)
  {
   control_max_ampl_sum = 1;
   max_ampl_sum_tol     = 0.001*max_ampl_sum;                   
  };
 }; 
 check_cc_NN_minimum(0.05);    
}

void print_parameters()
{
 logs_Write(0, "\n");
 print_metropolis_parameters();
 print_largeN_QFT_parameters();
}

double ptot(double c)
{
 return 4.0/c + 2.0*c*SQR(alpha_wc)/(lambda*(1.0 - c*alpha_wc))*(4.0 + lambda + 4.0/(1.0 - c*alpha_wc));
}

double dptot(double c)
{
 return -4.0/SQR(c) + 2.0*SQR(alpha_wc)/(lambda*SQR(1.0 - alpha_wc*c))*(8.0/(1.0 - c*alpha_wc) + lambda);
}

int tune_parameters()
{
 double c1, c2, cm, dpm, dp1, dp2;
 //Empirical initial guesses for c1 and c2
 c1 = 3.0/sqrt(lambda);
 c2 = 5.0/sqrt(lambda);
 dp1 = dptot(c1);
 dp2 = dptot(c2);
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
  dpm = dptot(cm);
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




