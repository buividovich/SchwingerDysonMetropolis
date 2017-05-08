#include "parameters.h"

char data_suffix[512]          =   "";   //Whether to save the histogram of sampling in sectors of different n and order
char scan_suffix[512]          =   "";
char scan_label[512]           =   "";
int  save_summary              =   0;
char raw_data_dir[512]         =   "C:\\DATA\\pcm_wc_mspace_raw_data";
char exact_data_dir[512]       =   "G:\\LAT\\sd_metropolis\\data\\pcm_wc_mspace_exact";
int  calculate_correlators     =   0;

int number_of_clusters         =   0;
char** cluster_data_suffixes   =   NULL;

//static const char * const cluster_data_suffixes[] = {"", "_idc", "_idc1", "_itep1", "_itep2"}; //
//#define NCLUSTERS (5)

static struct option long_options[] =
{
 LARGEN_QFT_LONG_OPTIONS,
 {          "data-suffix", required_argument, NULL,                      'a'},
 {          "scan-suffix", required_argument, NULL,                      'b'},
 {           "scan-label", required_argument, NULL,                      'c'},
 {        "cluster-label", required_argument, NULL,                      'd'},
 {         "raw-data-dir", required_argument, NULL,                      'e'},
 {       "exact-data-dir", required_argument, NULL,                      'f'},
 {"calculate-correlators",       no_argument, &calculate_correlators,      1},
 {                      0,                 0, NULL,                        0}
};

static char my_short_option_list[] = "a:b:c:d:e:f:";

int parse_command_line_options(int argc, char **argv)
{
 int gc, option_index;
 //Compose the joint list of all short options
 int short_option_list_length = 3;
 short_option_list_length += strlen(   largeN_QFT_short_option_list);
 short_option_list_length += strlen(           my_short_option_list);
 char* short_option_list = NULL;
 SAFE_MALLOC(short_option_list, char, short_option_list_length);
 sprintf(short_option_list, "%s",     my_short_option_list);
 strcat(short_option_list,    largeN_QFT_short_option_list);
 
 //Preparing the data to populate the cluster labels
 cluster_data_suffixes = (char** )malloc(sizeof(char*));
 number_of_clusters = 0;
 
 cluster_data_suffixes[0] = (char* )malloc(sizeof(char)*3);
 strcpy(cluster_data_suffixes[0], "");
 number_of_clusters ++;
 
 while(1)
 {
  gc = getopt_long(argc, argv, my_short_option_list, long_options, &option_index);
  if(gc==-1)
   break;
  switch(gc)
  {
   PARSE_LARGEN_QFT_OPTIONS; 
   case 'a':
    strcpy(data_suffix, optarg);
   break;
   case 'b':
    strcpy(scan_suffix, optarg);
   break;
   case 'c':
    strcpy(scan_label, optarg);
   break;
   case 'd':
    cluster_data_suffixes = (char** )realloc(cluster_data_suffixes, (number_of_clusters+1)*sizeof(char*));
    cluster_data_suffixes[number_of_clusters] = (char* )malloc(sizeof(char)*(strlen(optarg)+4));
    sprintf(cluster_data_suffixes[number_of_clusters], "_%s", optarg);
    number_of_clusters ++;
   break;
   case 'e':
    strcpy(raw_data_dir, optarg);
   break;
   case 'f':
    strcpy(exact_data_dir, optarg);
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
 if(strlen(data_suffix)==0)
  logs_WriteErrorAndTerminate("Data suffix should not be empty!!!");
 save_summary = (strlen(scan_suffix)>0 && strlen(scan_label)>0);  
}

void print_parameters()
{
 logs_Write(0, "\n");
 logs_WriteParameter(0,                   "LS",  "%i", LS);
 logs_WriteParameter(0,            "max_order",  "%i", max_order);
 logs_WriteParameter(0, "X space correlators?",  "%s", (calculate_correlators? "YES" : "NO"));
 logs_WriteParameter(0,          "Data suffix",  "%s", data_suffix);
 logs_WriteParameter(0,          "Scan suffix",  "%s", scan_suffix);
 logs_WriteParameter(0,           "Scan label",  "%s", scan_label);
 logs_WriteParameter(0,   "Raw data directory",  "%s", raw_data_dir);
 logs_WriteParameter(0, "Exact data directory",  "%s", exact_data_dir);
 logs_WriteParameter(0,    "Saving summaries?",   "%s", (save_summary? "YES" : "NO"));
 logs_Write(0, "CLUSTER LABELS (%i in total)", number_of_clusters);
 for(int icl=0; icl<number_of_clusters; icl++)
  logs_WriteParameter(0, "Cluster",  "%i: [%s]", icl, cluster_data_suffixes[icl]);
 
 logs_Write(0, "\n");
}
