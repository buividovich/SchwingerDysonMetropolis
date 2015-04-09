#ifndef _LARGEN_QFT_PARAMETERS_H_
#define _LARGEN_QFT_PARAMETERS_H_

#include <math.h>
#include <limits.h>
#include <clue_logs.h>
#include <clue_utils.h>

#include <sd_metropolis.h>
#include <sd_metropolis_parameters.h>
#include <sd_metropolis_statistics.h>

//Parameters of a generic SD equations for practically any large-N QFT
extern double   lambda;                   //tHooft coupling constant
extern double   meff_sq;                  //Square of the effective mass
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
extern int      min_observables_order;    //Minimal order of observables which are included into statistics in some formal expansion (e.g. SC/WC expansion)
extern int      max_observables_order;    //Correspondingly, maximal order
//Output files
extern char*    observables_file;         //File for the expectation values of the correlators
extern int      param_auto_tuning;        //Automatic tuning of transition amplitudes so that nAs are minimized
extern double   param_tuning_accuracy;    //Accuracy of parameter tuning
extern int      param_tuning_max_iter;    //Max. allowed number of iterations in param auto-tuning
//In debug mode, we can also check the stack consistency at every step
extern int      check_stack; 

void print_largeN_QFT_parameters();
void largeN_QFT_prefix(char* prefix); //Prints lambda, meff_sq, cc, NN, LT, LS to prefix

void      cc_NN_vicinity(double epsilon, double* data);
int  check_cc_NN_minimum(double tol);
int   find_cc_NN_minimum(double tol, double* min_val);


//Tunes cc and NN using a real MC process in such a way that <nA> is minimized...
double   tune_cc_NN_minimum(double tol, int tune_mc_steps); 

void init_lat_size_array();
//int check_cc_NN();

#define LARGEN_QFT_LONG_OPTIONS                                                      \
 {                     "lambda",  required_argument,                       NULL, 'A'}, \
 {                    "meff-sq",  required_argument,                       NULL, 'B'}, \
 {                         "cc",  required_argument,                       NULL, 'C'}, \
 {                         "NN",  required_argument,                       NULL, 'D'}, \
 {                        "DIM",  required_argument,                       NULL, 'E'}, \
 {                         "LT",  required_argument,                       NULL, 'F'}, \
 {                         "LS",  required_argument,                       NULL, 'G'}, \
 {              "max-stack-nel",  required_argument,                       NULL, 'H'}, \
 {            "max-history-nel",  required_argument,                       NULL, 'I'}, \
 {       "max-correlator-order",  required_argument,                       NULL, 'J'}, \
 {      "min-observables-order",  required_argument,                       NULL, 'K'}, \
 {      "max-observables-order",  required_argument,                       NULL, 'L'}, \
 {      "param-tuning-accuracy",  required_argument,                       NULL, 'M'}, \
 {      "param-tuning-max-iter",  required_argument,                       NULL, 'N'}, \
 {           "observables-file",  required_argument,                       NULL, 'O'}, \
 {       "no-param-auto-tuning",        no_argument,         &param_auto_tuning,   0}, \
 {             "no-stack-check",        no_argument,               &check_stack,   0}

#define PARSE_LARGEN_QFT_OPTIONS                                         \
   case 'A':                                                             \
    SAFE_SSCANF_BREAK(optarg, "%lf", lambda);                            \
   break;                                                                \
   case 'B':                                                             \
    SAFE_SSCANF_BREAK(optarg, "%lf", meff_sq);                           \
   break;                                                                \
   case 'C':                                                             \
    SAFE_SSCANF_BREAK(optarg, "%lf", cc);                                \
    param_auto_tuning = 0;                                               \
   break;                                                                \
   case 'D':                                                             \
    SAFE_SSCANF_BREAK(optarg, "%lf", NN);                                \
    param_auto_tuning = 0;                                               \
   break;                                                                \
   case 'E':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", DIM);                               \
    ASSERT(DIM<0);                                                       \
   break;                                                                \
   case 'F':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", LT);                                \
    ASSERT(LT<0);                                                        \
   break;                                                                \
   case 'G':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", LS);                                \
    ASSERT(LS<0);                                                        \
   break;                                                                \
   case 'H':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", max_stack_nel);                     \
    ASSERT(max_stack_nel<10);                                            \
   break;                                                                \
   case 'I':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", max_history_nel);                   \
    ASSERT(max_history_nel<10);                                          \
   break;                                                                \
   case 'J':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", max_correlator_order);              \
    ASSERT(max_correlator_order<0);                                      \
   break;                                                                \
   case 'K':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", min_observables_order);             \
   break;                                                                \
   case 'L':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", max_observables_order);             \
   break;                                                                \
   case 'M':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%lf", param_tuning_accuracy);            \
    ASSERT(param_tuning_accuracy<0.0);                                   \
   break;                                                                \
   case 'N':                                                             \
    SAFE_SSCANF_BREAK(optarg,   "%i", param_tuning_max_iter);            \
    ASSERT(param_tuning_max_iter<1);                                     \
   break;                                                                \
   case 'O':                                                             \
    COPY_FILE_NAME(optarg, observables_file);                            \
   break;

static const char largeN_QFT_short_option_list[] = "A:B:C:D:E:F:G:H:I:J:K:L:M:N:O:";

#endif
