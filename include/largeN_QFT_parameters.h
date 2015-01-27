#ifndef _LARGEN_QFT_PARAMETERS_H_
#define _LARGEN_QFT_PARAMETERS_H_

#include <clue_logs.h>

//Parameters of a generic SD equations for practically any large-N QFT
extern double   lambda;                   //tHooft coupling constant
extern double   cc;                       //Rescaling of observables according to the number of fields in the correlator
extern double   NN;                       //Overall rescaling of observables 
extern int      DIM;                      //Space-time dimensionality
extern int      LT;                       //Temporal size of the system
extern int      LS;                       //Spatial size of the system
//Parameters of a generic stack-based MC process
extern int      max_stack_nel;            //Maximal number of elements in the stack characterizing the system state
extern int      max_history_nel;          //Maximal number of elements in the stack containing the history of momenta contractions
//Parameters of the statistical analysis of the MC data
extern int      max_correlator_order;     //Maximal correlator order to trace
//Output files
extern char*    observables_file;         //File for the expectation values of the correlators
extern int      param_auto_tuning;        //Automatic tuning of transition amplitudes so that nAs are minimized
extern double   param_tuning_accuracy;    //Accuracy of parameter tuning

void print_largeN_QFT_parameters();

void largeN_QFT_prefix(char* prefix); //Prints lambda, cc, NN, LT, LS to prefix

//int check_cc_NN();

#define LARGEN_QFT_LONG_OPTIONS                                                      \
 {                   "lambda",  required_argument,                       NULL, 'A'}, \
 {                       "cc",  required_argument,                       NULL, 'B'}, \
 {                       "NN",  required_argument,                       NULL, 'C'}, \
 {                      "DIM",  required_argument,                       NULL, 'D'}, \
 {                       "LT",  required_argument,                       NULL, 'E'}, \
 {                       "LS",  required_argument,                       NULL, 'F'}, \
 {            "max_stack_nel",  required_argument,                       NULL, 'G'}, \
 {          "max_history_nel",  required_argument,                       NULL, 'H'}, \
 {     "max-correlator-order",  required_argument,                       NULL, 'I'}, \
 {    "param-tuning-accuracy",  required_argument,                       NULL, 'J'}, \
 {         "observables-file",  required_argument,                       NULL, 'K'}, \
 {     "no-param-auto-tuning",        no_argument,         &param_auto_tuning,   0}

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
    SAFE_SSCANF_BREAK(optarg,  "%i", DIM);                               \
    ASSERT(DIM<0);                                                       \
   break;                                                                \
   case 'E':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", LT);                                \
    ASSERT(LT<0);                                                        \
   break;                                                                \
   case 'F':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", LS);                                \
    ASSERT(LS<0);                                                        \
   break;                                                                \
   case 'G':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", max_stack_nel);                     \
    ASSERT(max_stack_nel<10);                                            \
   break;                                                                \
   case 'H':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", max_history_nel);                   \
    ASSERT(max_history_nel<10);                                          \
   break;                                                                \
   case 'I':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", max_correlator_order);              \
    ASSERT(max_correlator_order<0);                                      \
   break;                                                                \
   case 'J':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%lf", param_tuning_accuracy);            \
    ASSERT(param_tuning_accuracy<0.0);                                   \
   break;                                                                \
   case 'K':                                                             \
    COPY_FILE_NAME(optarg, observables_file);                            \
   break;

static const char largeN_QFT_short_option_list[] = "A:B:C:D:E:F:G:H:I:J:K:";

#endif
