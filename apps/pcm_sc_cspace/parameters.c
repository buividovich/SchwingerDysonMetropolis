#include "parameters.h"

double meff_sq              =   1.0; //Effective mass term - for resummations 
//Some switches
int    save_sampling_hist   =   0; //Whether to save the histogram of sampling in sectors of different n and order
int    save_correlators     =   0; //Whether to save correlators <gx gy>
int    resummation          =   0; //Whether to re-sum the kinetic "random walk" term
//Useful calculable parameters
double beta                 =   0.0; // == 1/\lambda, strong-coupling expansion parameter
double sigma                =   0.0; // Self-energy, is trivial without resummation
double mass2                =   0.0;   //Squared mass determining the pole of the bare propagator
double coord_factor         =   1.0;

static struct option long_options[] =
{
 METROPOLIS_LONG_OPTIONS,
 LARGEN_QFT_LONG_OPTIONS,
 {               "meff-sq",  required_argument,                 NULL, 'a'},
 {    "save-sampling-hist",        no_argument,  &save_sampling_hist,   1},
 {      "save-correlators",        no_argument,    &save_correlators,   1},
 {           "resummation",        no_argument,         &resummation,   1},
 {                       0,                  0,                 NULL,   0}
};

static char my_short_option_list[] = "a:";

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
    SAFE_SSCANF_BREAK(optarg, "%lf", meff_sq);
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
 beta = 1.0/lambda;
 /*if(!resummation)
  meff_sq = -2.0*(double)(DIM);*/

 mass2 = lambda + meff_sq;
 
 if(DIM==1 && LT==2)
  coord_factor = 2.0;
 if(DIM==2 && ((LT==2 && LS>2 ) || (LT>2 && LS==2)))
  coord_factor = 2.0;
 if(DIM==2 && LT==2 && LS==2)
  coord_factor = 4.0;

 init_lattice();

 if(strlen(data_dir)==0)
  sprintf(data_dir, "./data/pcm_sc_cspace/");
 if(strlen(suffix)==0)
  largeN_QFT_suffix(suffix);

 if(resummation)
  snprintf(&(suffix[strlen(suffix)]), 512, "_rsm");

 if(strlen(label)>0)
  snprintf(&(suffix[strlen(suffix)]), 512, "_%s", label);
}

void print_parameters()
{
 logs_Write(0, "\n");
 print_metropolis_parameters();
 print_largeN_QFT_parameters();
 logs_WriteParameter(0,                    "beta=1/lambda",   "%2.4lf", beta);
 if(resummation)
 {
  logs_WriteParameter(0,              "Effective mass term",   "%2.4lf", meff_sq);
  logs_WriteParameter(0, "Pole mass in the bare propagator",   "%2.4lf", mass2);
 }; 
 
 logs_WriteParameter(0,      "Resumming the kinetic term?",       "%s", (resummation?        "YES" : "NO"));
 logs_WriteParameter(0,  "Saving the sampling histograms?",       "%s", (save_sampling_hist? "YES" : "NO"));
 logs_WriteParameter(0,          "Saving the correlators?",       "%s", (save_correlators?   "YES" : "NO"));
  
 logs_Write(0, "\n"); 
}
