#include "parameters.h"

//Some switches
int    save_sampling_hist   =   0;  //Whether to save the histogram of sampling in sectors of different n and order
int    save_correlators     =   0;  //Whether to save correlators <gx gy>
int    save_metropolis_stat =   0;  //Whether to save the data on the performance of Metropolis algorithm

//Useful calculable parameters
double    stereo_alpha    = 0.0;  // -\lambda/8, expansion parameter...

static struct option long_options[] =
{
 METROPOLIS_LONG_OPTIONS,
 LARGEN_QFT_LONG_OPTIONS,
 {    "save-sampling-hist",        no_argument,    &save_sampling_hist,   1},
 {      "save-correlators",        no_argument,      &save_correlators,   1},
 {  "save-metropolis-stat",        no_argument,  &save_metropolis_stat,   1},
 {                       0,                  0,                   NULL,   0}
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
 init_lat_propagator(&P, 0.25*lambda, 1.0);
 
 stereo_alpha = -0.125*lambda;
 
 if(strlen(data_dir)==0)
  sprintf(data_dir, "./data/pcm_wc_mspace/");
 if(strlen(suffix)==0)
  largeN_QFT_suffix(suffix);
 
 if(strlen(label)>0)
  snprintf(&(suffix[strlen(suffix)]), 512, "_%s", label); 
}

void print_parameters()
{
 logs_Write(0, "\n");
 print_metropolis_parameters();
 print_largeN_QFT_parameters();
 logs_WriteParameter(0,  "Saving the sampling histograms?",  "%s", (save_sampling_hist?     "YES" : "NO"));
 logs_WriteParameter(0,          "Saving the correlators?",  "%s", (save_correlators?       "YES" : "NO"));
 logs_WriteParameter(0,    "Saving Metropolis statistics?",  "%s", (save_metropolis_stat?   "YES" : "NO"));
 
 logs_Write(0, "\tPARAMETERS OF LATTICE PROPAGATOR");
 logs_WriteParameter(0,            "Sigma", "%2.4E", P.sigma);
 logs_WriteParameter(0,     "Mass squared", "%2.4E", P.mass_sq);
 logs_WriteParameter(0,          "lat_dim",    "%i", lat_dim);
 logs_WriteParameter(0,          "lat_vol",    "%i", lat_vol);
 for(int mu=0; mu<lat_dim; mu++)
  logs_WriteParameter(0, "lat_size",    "[%i]: %i", mu, lat_size[mu]);
 logs_Write(0, "\n"); 
}
