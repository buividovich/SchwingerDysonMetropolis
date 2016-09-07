#ifndef _LARGEN_QFT_PARAMETERS_H_
#define _LARGEN_QFT_PARAMETERS_H_

#include <math.h>
#include <limits.h>
#include <clue_logs.h>
#include <clue_utils.h>

//Parameters of a generic SD equations for practically any large-N QFT
extern int      DIM;                      //Dimension of space
extern int      LT;                       //Temporal size of the system
extern int      LS;                       //Spatial size of the system
extern double   lambda;                   //tHooft coupling constant
extern int      max_order;                //Max. order of the expansion
//Algorithmic parameters which control statistical sampling
extern double   alpha;                    //Rescaling of observables according to the order
extern double   cc;                       //Rescaling of observables according to the number of fields in the correlator
extern double   NN;                       //Overall rescaling of observables, also genus-dependent
//In debug mode, we can also check the stack consistency at every step
extern int      check_stack;
//Output parameters
extern char     data_dir[512];
extern char     suffix[512];

#define LARGEN_QFT_LONG_OPTIONS                                                        \
 {                        "DIM",  required_argument,                       NULL, 'A'}, \
 {                         "LT",  required_argument,                       NULL, 'B'}, \
 {                         "LS",  required_argument,                       NULL, 'C'}, \
 {                     "lambda",  required_argument,                       NULL, 'D'}, \
 {                  "max-order",  required_argument,                       NULL, 'E'}, \
 {                      "alpha",  required_argument,                       NULL, 'F'}, \
 {                         "cc",  required_argument,                       NULL, 'G'}, \
 {                         "NN",  required_argument,                       NULL, 'H'}, \
 {                   "data-dir",  required_argument,                       NULL, 'I'}, \
 {                     "suffix",  required_argument,                       NULL, 'J'}, \
 {                "check-stack",        no_argument,               &check_stack,   1}

#define PARSE_LARGEN_QFT_OPTIONS                                         \
   case 'A':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", DIM);                               \
    ASSERT(DIM<0);                                                       \
   break;                                                                \
   case 'B':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", LT);                                \
    ASSERT(LT<0);                                                        \
   break;                                                                \
   case 'C':                                                             \
    SAFE_SSCANF_BREAK(optarg,  "%i", LS);                                \
    ASSERT(LS<0);                                                        \
   break;                                                                \
   case 'D':                                                             \
    SAFE_SSCANF_BREAK(optarg, "%lf", lambda);                            \
   break;                                                                \
   case 'E':                                                             \
    SAFE_SSCANF_BREAK(optarg, "%i", max_order);                          \
   break;                                                                \
   case 'F':                                                             \
    SAFE_SSCANF_BREAK(optarg, "%lf", alpha);                             \
   break;                                                                \
   case 'G':                                                             \
    SAFE_SSCANF_BREAK(optarg, "%lf", cc);                                \
   break;                                                                \
   case 'H':                                                             \
    SAFE_SSCANF_BREAK(optarg, "%lf", NN);                                \
   break;                                                                \
   case 'I':                                                             \
    strcpy(optarg, data_dir);                                            \
   break;                                                                \
   case 'J':                                                             \
    strcpy(optarg, suffix);                                              \
   break;

static const char largeN_QFT_short_option_list[] = "A:B:C:D:E:F:G:H:I:J:";

void print_largeN_QFT_parameters();
void       largeN_QFT_suffix(char* s);

#endif
