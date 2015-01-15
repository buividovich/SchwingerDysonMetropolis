#ifndef _SD_METROPOLIS_PARAMETERS_H_
#define _SD_METROPOLIS_PARAMETERS_H_

#include <clue_logs.h>

extern int      max_recursion_depth; //Max. possible depth of the sequence of states
extern double   p_plus;              //Probability of the "Forward" move
extern int      prod_mc_steps;       //Number of MC steps for production runs
extern int      therm_mc_steps;      //Number of MC steps for thermalization
extern int      mc_interval;         //Interval between the successive presumably uncorrelated measurements
extern char*    mc_stat_file;        //File for the quantities characterizing the MC process itself
extern char*    ns_history_file;     //File for saving the MC history of ns

#define METROPOLIS_LONG_OPTIONS                                            \
 {      "max_recursion_depth",  required_argument,      NULL,  '0'},       \
 {                   "p_plus",  required_argument,      NULL,  '1'},       \
 {            "prod_mc_steps",  required_argument,      NULL,  '2'},       \
 {           "therm_mc_steps",  required_argument,      NULL,  '3'},       \
 {              "mc_interval",  required_argument,      NULL,  '4'},       \
 {             "mc_stat_file",  required_argument,      NULL,  '5'},       \
 {          "ns_history_file",  required_argument,      NULL,  '6'}
 
#define PARSE_METROPOLIS_OPTIONS                                           \
 case '0':                                                                 \
  SAFE_SSCANF_BREAK(optarg,  "%i", max_recursion_depth);                   \
  ASSERT(max_recursion_depth<0);                                           \
 break;                                                                    \
 case '1':                                                                 \
  SAFE_SSCANF_BREAK(optarg, "%lf", p_plus);                                \
  ASSERT((p_plus<0.0) || (p_plus>1.0));                                    \
 break;                                                                    \
 case '2':                                                                 \
  SAFE_SSCANF_BREAK(optarg,  "%i", prod_mc_steps);                         \
  ASSERT(prod_mc_steps<0);                                                 \
 break;                                                                    \
 case '3':                                                                 \
  SAFE_SSCANF_BREAK(optarg,  "%i", therm_mc_steps);                        \
  ASSERT(therm_mc_steps<0);                                                \
 break;                                                                    \
 case '4':                                                                 \
  SAFE_SSCANF_BREAK(optarg,  "%i", mc_interval);                           \
  ASSERT(mc_interval<0);                                                   \
 break;                                                                    \
 case '5':                                                                 \
  COPY_FILE_NAME(optarg, mc_stat_file);                                    \
 break;                                                                    \
 case '6':                                                                 \
  COPY_FILE_NAME(optarg, ns_history_file);                                 \
 break;
 
static const char  metropolis_short_option_list[] = "0:1:";  

void print_metropolis_parameters();                                                                

#endif
