#include "parameters.h"

static struct option long_options[] =
{
 METROPOLIS_LONG_OPTIONS,
 LARGEN_QFT_LONG_OPTIONS,
 LATTICE_STACK_LONG_OPTIONS,
 {                       0,                  0,                       NULL,   0}
};

static char my_short_option_list[] = "";

int parse_command_line_options(int argc, char **argv)
{
 int gc, option_index;
 //Compose the joint list of all short options
 int short_option_list_length = 3;
 short_option_list_length += strlen(   metropolis_short_option_list);
 short_option_list_length += strlen(   largeN_QFT_short_option_list);
 short_option_list_length += strlen(lattice_stack_short_option_list);
 short_option_list_length += strlen(           my_short_option_list);
 char* short_option_list = NULL;
 SAFE_MALLOC(short_option_list, char, short_option_list_length);
 sprintf(short_option_list, "%s",     my_short_option_list);
 strcat(short_option_list,    metropolis_short_option_list);
 strcat(short_option_list,    largeN_QFT_short_option_list);
 strcat(short_option_list, lattice_stack_short_option_list);

 while(1)
 {
  gc = getopt_long(argc, argv, short_option_list, long_options, &option_index);
  if(gc==-1)
   break;
  switch(gc)
  {
   PARSE_METROPOLIS_OPTIONS;
   PARSE_LARGEN_QFT_OPTIONS;
   PARSE_LATTICE_STACK_OPTIONS;
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
 return 0;
}

void init_parameters()
{
 if(param_auto_tuning)
 {
  cc = 2.0/sqrt(fabs(lambda));
  NN = cc;
 };    
}

void print_parameters()
{
 print_metropolis_parameters();
 print_largeN_QFT_parameters();
 print_lattice_stack_parameters();
}



