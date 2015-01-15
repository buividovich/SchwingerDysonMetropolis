#ifndef _SD_METROPOLIS_PARAMETERS_H_
#define _SD_METROPOLIS_PARAMETERS_H_

#include <clue_logs.h>

extern int    max_recursion_depth;
extern double p_plus;

#define METROPOLIS_LONG_OPTIONS                                            \
 {"max_recursion_depth",        required_argument,      NULL,  '0'},       \
 {             "p_plus",   	    required_argument,      NULL,  '1'}
 
#define PARSE_METROPOLIS_OPTIONS                                           \
 case '0':                                                                 \
  SAFE_SSCANF_BREAK(optarg,  "%i", max_recursion_depth);                   \
  ASSERT(max_recursion_depth>1);                                           \
 break;                                                                    \
 case '1':                                                                 \
  SAFE_SSCANF_BREAK(optarg, "%lf", p_plus);                                \
  ASSERT((p_plus<0.0) || (p_plus>1.0));                                    \
 break;  
 
static const char  metropolis_short_option_list[] = "0:1:";  

void print_metropolis_parameters();                                                                

#endif
