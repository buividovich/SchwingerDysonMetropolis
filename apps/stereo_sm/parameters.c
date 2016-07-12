#include "parameters.h"

double    alpha           = 0.002;
int       max_alpha_order = 50;
double    average_seq_len = 6.0; //Desired number of momenta in the sequence
double    average_num_seq = 3.0; //Desired number of sequences

static struct option long_options[] =
{
 METROPOLIS_LONG_OPTIONS,
 LARGEN_QFT_LONG_OPTIONS,
 {                 "alpha", required_argument,  NULL, 'a'},
 {       "max-alpha-order", required_argument,  NULL, 'b'},
 {       "average-seq-len", required_argument,  NULL, 'c'},
 {       "average-num-seq", required_argument,  NULL, 'd'},
 {                       0,                  0,                       NULL,   0}
};

static char my_short_option_list[] = "a:b:c:d:";

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
   case   'c':
    SAFE_SSCANF_BREAK(optarg, "%lf", average_seq_len);
   break;
   case   'd':
    SAFE_SSCANF_BREAK(optarg, "%lf", average_num_seq);
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
 //More motivated tuning
 if(param_auto_tuning)
 {
  double x = 0.25 - 0.25/SQR((double)(average_seq_len));
  cc = P.sigma/x;
  double u = (1.0 - sqrt(1.0 - 4.0*x))/(2.0*x) - 1.0;
  NN = u*(double)(average_num_seq)/(double)(average_num_seq - 1);
 }; 
}

void free_parameters()
{
 free_genus_constants();
 free_largeN_QFT_parameters();
}

void print_parameters()
{
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
 for(int mu=0; mu<lat_dim; mu++)
  logs_WriteParameter(0, "lat_size",    "[%i]: %i", mu, lat_size[mu]);
 logs_Write(0, "\tPARAMETERS OF STOCHASTIC SAMPLING");
 logs_WriteParameter(0, "Desired average sequence length", "%2.2lf", average_seq_len);
 logs_WriteParameter(0, "Desired average no.  sequences",  "%2.2lf", average_num_seq); 
 print_max_amplitudes(); 
}
