#ifndef _SD_METROPOLIS_PARAMETERS_H_
#define _SD_METROPOLIS_PARAMETERS_H_

#include <math.h>
#include <limits.h>

#include <clue_logs.h>

extern int      max_recursion_depth;    //Max. possible depth of the sequence of states
extern double   p_plus;                 //Probability of the "Forward" move
extern int      p_plus_tuning;          //If 1 (default), p_plus will be automatically tuned to maximize acceptance
extern int      p_plus_tuning_interval; //p_plus will be updated in this interval
extern int      number_mc_steps;        //Number of MC steps for production runs
extern int      mc_reporting_interval;  //The interval at which the current state of the MC process is reported
extern int      exit_upon_overflow;     //Whether to stop the MC process once the history or the stack overflow is detected

//0-9, X - Z  
#define METROPOLIS_LONG_OPTIONS                                                            \
 {      "max-recursion-depth",  required_argument,                      NULL,  '0'},       \
 {                   "p-plus",  required_argument,                      NULL,  '1'},       \
 {   "p-plus-tuning-interval",  required_argument,                      NULL,  '2'},       \
 {          "number-mc-steps",  required_argument,                      NULL,  '3'},       \
 {    "mc-reporting-interval",  required_argument,                      NULL,  '4'},       \
 {         "logs-noise-level",  required_argument,                      NULL,  'X'},       \
 {          "no-pplus-tuning",        no_argument,            &p_plus_tuning,    0},       \
 {           "no-ansi-colors",        no_argument,              &ansi_colors,    0},       \
 {   "print-errors-to-stderr",        no_argument,   &print_errors_to_stderr,    1},       \
 {       "exit-upon-overflow",        no_argument,       &exit_upon_overflow,    1}
 
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
  SAFE_SSCANF_BREAK(optarg,  "%i", p_plus_tuning_interval);                \
 break;                                                                    \
 case '3':                                                                 \
  SAFE_SSCANF_BREAK(optarg,  "%i", number_mc_steps);                       \
  ASSERT(number_mc_steps<0);                                               \
 break;                                                                    \
 case '4':                                                                 \
  SAFE_SSCANF_BREAK(optarg,  "%i", mc_reporting_interval);                 \
  ASSERT(mc_reporting_interval<0);                                         \
 break;                                                                    \
 case 'X':                                                                 \
  SAFE_SSCANF_BREAK(optarg,  "%i", logs_noise_level);                      \
 break;                                                                    

static const char metropolis_short_option_list[] = "0:1:2:3:4:X:";

void print_metropolis_parameters();

#endif
