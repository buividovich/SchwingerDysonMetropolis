#include "parameters.h"

char* nphi_str = NULL;
int   ntr      = 0;    //Number of traces in the correlator
int*  ngs      = NULL; //Array with number of fields in each trace 
int   nphi     = 0;    //Total number of field operators inside each trace
int   npairs   = 0;    //Total number of Wick contractions of all the fields

static struct option long_options[] =
{
 METROPOLIS_LONG_OPTIONS, //0-9, Z, Y, X
 LARGEN_QFT_LONG_OPTIONS, //A-P
 {    "nphi-str",  required_argument,    NULL,   'a'},
 {             0,                  0,    NULL,     0}
};

static char my_short_option_list[] = "a";

int parse_nphi_str(char* s, int** ns)
{ 
 int nnums = 0;  
 int slen  = strlen(s);
 long int lnum = 0;
 char* endptr = NULL;
 char* spart  = s;
 do{
  lnum = strtol(spart, &endptr, 10);
  spart = endptr;
  if(lnum!=0)
  {
   nnums ++;
   int num = (int)lnum;
   logs_Write(3, "%i: %i", nnums, num);
   SAFE_REALLOC((*ns), int, nnums);
   (*ns)[nnums-1] = num;         
  };
 }while((endptr!=&(s[slen]) && lnum!=0));
 return nnums;
}

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
    SAFE_MALLOC(nphi_str, char, strlen(optarg) + 1);
    strcpy(nphi_str, optarg);
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
 
 ASSERT(nphi_str==NULL);
 
 ntr = parse_nphi_str(nphi_str, &ngs);
 char* my_nphi_str = NULL;
 sprintf_append(&my_nphi_str, "[%i", ngs[0]);
 nphi = ngs[0];
 for(int i=1; i<ntr; i++)
 {
  nphi += ngs[i];
  sprintf_append(&my_nphi_str, " %i", ngs[i]);
 };
 sprintf_append(&my_nphi_str, "]");
 npairs = get_num_pairings(nphi);
 
 logs_WriteParameter(0,       "Configuration of correlator", "%s", my_nphi_str);
 logs_WriteParameter(0,            "Total number of traces", "%i", ntr);
 logs_WriteParameter(0,            "Total number of fields", "%i", nphi);
 logs_WriteParameter(0, "Total number of Wick contractions", "%i", npairs);
 
 SAFE_FREE(my_nphi_str);
 SAFE_FREE(short_option_list);  
 return 0;
}
