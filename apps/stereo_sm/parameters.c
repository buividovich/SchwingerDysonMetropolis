#include "parameters.h"

double alpha           = 0.0;
int    max_alpha_order = 50;

static struct option long_options[] =
{
 METROPOLIS_LONG_OPTIONS,
 LARGEN_QFT_LONG_OPTIONS,
 {                 "alpha",  required_argument,                       NULL, 'a'},
 {       "max-alpha-order",  required_argument,                       NULL, 'b'},
 {                       0,                  0,                       NULL,   0}
};

static char my_short_option_list[] = "a:b:";

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
   case   'a':
    SAFE_SSCANF_BREAK(optarg, "%lf", alpha);
   break;
   case   'b':
    SAFE_SSCANF_BREAK(optarg,  "%i", max_alpha_order);
   break;
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
 init_lat_propagator(&P, 1, 0.25*lambda);
 
 double* my_params[2] = {&cc, &NN};
 if(param_auto_tuning)
 {
  cc = 1.0; NN = 1.0;
  find_param_minimum(param_tuning_accuracy, my_params, &max_ampl_sum, 2);
 };
 check_param_minimum(param_tuning_accuracy, my_params, 2);
}

void free_parameters()
{
 free_genus_constants();
 free_largeN_QFT_parameters();
}

void print_parameters()
{
 int mu;
 logs_Write(0, "\n");
 print_metropolis_parameters();
 print_largeN_QFT_parameters( 1, 1);
 logs_Write(0, "\tPARAMETERS OF LATTICE PROPAGATOR");
 logs_WriteParameter(0,            "alpha", "%2.4E", alpha);
 logs_WriteParameter(0, "Max. alpha order",    "%i", max_alpha_order);
 logs_WriteParameter(0,            "Sigma", "%2.4E", P.sigma);
 logs_WriteParameter(0,     "Mass squared", "%2.4E", P.mass_sq);
 logs_WriteParameter(0,          "lat_dim",    "%i", lat_dim);
 logs_WriteParameter(0,          "lat_vol",    "%i", lat_vol);
 for(mu=0; mu<lat_dim; mu++)
  logs_WriteParameter(0, "lat_size",    "[%i]: %i", mu, lat_size[mu]);
 print_max_amplitudes(); 
}
