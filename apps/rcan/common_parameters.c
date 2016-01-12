#include "common_parameters.h"

int    mmax          =  0;          //Maximal order
int    mmin_prc      =  0;          //Minimal order for which the results are saved in memory

char*  out_file      =  NULL;
int    append_mode   =  0;
int    auto_naming   =  1;

static struct option long_options[] =
{
 LARGEN_QFT_LONG_OPTIONS, //A-P
 METROPOLIS_LONG_OPTIONS, //0-9, X - Z  
 {           "mmax",  required_argument,          NULL,   'a'},
 {       "mmin-prc",  required_argument,          NULL,   'b'},
 {       "out-file",  required_argument,          NULL,   'c'},
 { "no-auto-naming",        no_argument,  &auto_naming,     0},
 {    "append-mode",        no_argument,  &append_mode,     1},
 {                0,                  0,          NULL,     0}
};

static char my_short_option_list[] = "a:b:c:";

int parse_common_command_line_options(int argc, char **argv)
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
    SAFE_SSCANF(optarg,  "%i", mmin_prc);  
   break;
   case 'c':
    COPY_FILE_NAME(optarg, out_file);
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

void print_common_parameters()
{
 print_largeN_QFT_parameters(1, 0);
 logs_Write(0, "\tGENERIC PARAMETERS OF SERIES GENERATOR: ");
 logs_WriteParameter(0,                    "Max. order",                 "%i", mmax);
 logs_WriteParameter(0,                "sizeof(double)",           "%i bytes", sizeof(double));
 if(out_file!=NULL)
 logs_WriteParameter(0,                      "Output file",                 "%s %s %s", out_file, (auto_naming? "[Auto naming]" : ""), (append_mode? "[Append mode]" : ""));
}

void init_common_parameters()
{
 if(out_file==NULL && auto_naming)
 {
  sprintf_append(&out_file, "./data/smlm_analytics/lu_s%i_l%2.2lf.dat", LS, lambda);
 };
}

