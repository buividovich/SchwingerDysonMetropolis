#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <stdlib.h>
#include <stdio.h>
#include <getopt.h> 
#include <clue_utils.h>
#include <clue_logs.h>

extern int      maxn;                  //Follow the expansion up to this order
extern double   lambda;                //Coupling constant in the Hermitian matrix model
extern double   cc;                    //Rescaling of observables
extern double   NN;                    //Rescaling of observables
extern double   pplus;                 //Probability of the forward evolution
extern int      nmc;                   //Number of MC steps
extern int      maxg;                  //Maximal correlator order to trace
extern char*    ns_hist_file;          //File for saving the histogram of the recursion depth
extern char*    observables_file;      //File for the expectation values of the correlators
extern char*    mc_stat_file;          //File for the quantities characterizing the MC process itself
extern char*    ns_history_file;       //File for saving the MC history of NS

extern int      param_auto_tuning;     //Automatic tuning of NN and c so that the transition probabilities are minimized

int  parse_command_line_options(int argc, char **argv);
void init_parameters();
void print_parameters();
void printhelp();

#endif
