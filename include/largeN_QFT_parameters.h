#ifndef _LARGEN_QFT_PARAMETERS_H_
#define _LARGEN_QFT_PARAMETERS_H_

#include <clue_logs.h>

//Parameters of a generic SD equations for practically any large-N QFT
extern double   lambda;                   //tHooft coupling constant
extern double   cc;                       //Rescaling of observables according to the number of fields in the correlator
extern double   NN;                       //Overall rescaling of observables 
//Parameters of the statistical analysis of the MC data
extern int      max_correlator_order;     //Maximal correlator order to trace
//Output files
extern char*    observables_file;         //File for the expectation values of the correlators
extern int      param_auto_tuning;        //Automatic tuning of transition amplitudes so that nAs are minimized

void print_largeN_QFT_parameters();

#define LARGEN_QFT_LONG_OPTIONS                                                      \
 {                   "lambda",  required_argument,                       NULL, 'A'}, \
 {                       "cc",  required_argument,                       NULL, 'B'}, \
 {                       "NN",  required_argument,                       NULL, 'C'}, \
 {     "max-correlator-order",  required_argument,                       NULL, 'D'}, \
 {         "observables-file",  required_argument,                       NULL, 'E'}, \
 {        "param-auto-tuning",        no_argument,         &param_auto_tuning,   1}

#define PARSE_LARGEN_QFT_OPTIONS                                         \
   case 'A':                                                             \
    SAFE_SSCANF_BREAK(optarg, "%lf", lambda);                            \
   break;                                                                \
   case 'B':                                                             \
    SAFE_SSCANF_BREAK(optarg, "%lf", cc);                                \
   break;                                                                \
   case 'C':                                                             \
    SAFE_SSCANF_BREAK(optarg, "%lf", NN);                                \
   break;                                                                \
   case 'D':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", max_correlator_order);              \
    ASSERT(max_correlator_order<0);                                      \
   break;                                                                \
   case 'E':                                                             \
    COPY_FILE_NAME(optarg, observables_file);                            \
   break;

static const char largeN_QFT_short_option_list[] = "A:B:C:D:E:";


#endif
