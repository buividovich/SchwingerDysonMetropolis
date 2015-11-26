#include "parameters.h"

int    mmax        =  1;
int    s1          = +1;          //sigma parameters controlling the thimble structure
int    s2s3        = -1;          //ratio of s2 to s3
double mean_link   = -1.0;
double m2          =  0.1;         //Mass squared 

static struct option long_options[] =
{
 LARGEN_QFT_LONG_OPTIONS, //A-P
 {        "mmax",  required_argument,    NULL,   'a'},
 {          "s1",  required_argument,    NULL,   'b'},
 {        "s2s3",  required_argument,    NULL,   'c'},
 {   "mean-link",  required_argument,    NULL,   'd'},
 {             0,                  0,    NULL,     0}
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
   case 'a':
    SAFE_SSCANF(optarg,  "%i", mmax);  
   break;
   case 'b':
    SAFE_SSCANF(optarg,  "%i",   s1);    
   break;
   case 'c':
    SAFE_SSCANF(optarg,  "%i", s2s3);    
   break;
   case 'd':
    SAFE_SSCANF(optarg, "%lf", mean_link);    
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

void print_parameters()
{
 print_largeN_QFT_parameters();
 logs_Write(0, "\t\tPARAMETERS OF SERIES GENERATOR: ");
 logs_WriteParameter(0,   "Max. order",                 "%i", mmax);
 logs_WriteParameter(0,  "[s1, s2/s3]",       "[%+1i, %+1i]", s1, s2s3);
 logs_WriteParameter(0,    "Mean link",             "%2.4lf", mean_link);
 logs_WriteParameter(0, "Mass squared",             "%2.4lf", m2);
}

void init_parameters()
{
 if(mean_link<0.0)
 {
  if(LS==2)   
   mean_link = (lambda < 4.0? 1.0 - 0.125*lambda : 2.0/lambda);
  else
   mean_link = 1.0; 
 };
 
 m2 = s2s3*(2*DIM*(1.0 - mean_link) - s1*lambda);
}

