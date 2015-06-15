#include "parameters.h"

static struct option long_options[] =
{
 METROPOLIS_LONG_OPTIONS, //0-9, Z, Y, X
 LARGEN_QFT_LONG_OPTIONS, //A-P
 {             0,                  0,    NULL,   0}
};

static char my_short_option_list[] = "";

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
 double* my_params[2]     = {&cc, &NN};
 int     n_tunable_params = 2;
 if(param_auto_tuning)
 {
  find_param_minimum(param_tuning_accuracy, my_params, &max_ampl_sum, n_tunable_params);  
  if(max_ampl_sum>0.0)
  {
   control_max_ampl_sum = 1;
   max_ampl_sum_tol     = 0.001*max_ampl_sum;                   
  };
 };
 init_genus_constants(1);
 check_param_minimum(0.05, my_params, n_tunable_params);
 genus = 0; 
}

void free_parameters()
{
 free_genus_constants();
 free_largeN_QFT_parameters();
 free_metropolis_parameters();
}

void print_parameters()
{
 logs_Write(0, "\n");
 print_metropolis_parameters();
 print_largeN_QFT_parameters();
}



