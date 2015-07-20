#include "parameters.h"

int    model   = 0;   //0 is for free Gaussian, 1 is for LM with inverses, 2 is for LM with positive powers
double epsilon = 0.0;
//Calculable and model-dependent parameters
double   create_amplitudes[NMODELS] = {0.0, 0.0, 0.0};
double increase_amplitudes[NMODELS] = {0.0, 0.0, 0.0};
double     join_amplitudes[NMODELS] = {0.0, 0.0, 0.0};


static struct option long_options[] =
{
 METROPOLIS_LONG_OPTIONS, //0-9, Z, Y, X
 LARGEN_QFT_LONG_OPTIONS, //A-V
 {               "epsilon",  required_argument,                       NULL, 'a'},
 {                 "model",  required_argument,                       NULL, 'b'},
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
   case 'a':
    SAFE_SSCANF_BREAK(optarg, "%lf", epsilon);
   break;
   case 'b':
    SAFE_SSCANF_BREAK(optarg, "%i", model);
    ASSERT(model<0 || model>=NMODELS);
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
 double* my_params[2] = {&cc, &NN};
 if(param_auto_tuning)
 {
  cc = 1.0; NN = 1.0;
  find_param_minimum(param_tuning_accuracy, my_params, &max_ampl_sum, 2);
 };
 check_param_minimum(param_tuning_accuracy, my_params, 2);
 //Setting the model-dependent amplitudes
   create_amplitudes[FREE_GAUSS]  = 1.0/(NN*cc);
   create_amplitudes[LM_INVERSE]  = 1.0/(NN*cc);
   create_amplitudes[LM_POSITIVE] = 2.0/(NN*cc);
 
 increase_amplitudes[FREE_GAUSS]  = 2.0/cc;
 increase_amplitudes[LM_INVERSE]  = 1.0/cc;
 increase_amplitudes[LM_POSITIVE] = 3.0/cc;
 
     join_amplitudes[FREE_GAUSS]  = NN/cc;
     join_amplitudes[LM_INVERSE]  = NN;
     join_amplitudes[LM_POSITIVE] = NN/cc;
}

void free_parameters()
{
 free_genus_constants();
 free_largeN_QFT_parameters();
}

void print_parameters()
{
 logs_Write(0, "\n");
 logs_Write(0, "\tPARAMETERS OF THE ALGORITHM");
 logs_WriteParameter(0, "Model",    "%i %s", model, model_names[model]);
 logs_WriteParameter(0, "epsilon", "%2.4lf", epsilon);
 print_metropolis_parameters();
 print_largeN_QFT_parameters();
 print_max_amplitudes(); 
}
