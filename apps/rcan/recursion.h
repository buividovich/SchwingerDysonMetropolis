#ifndef _RECURSION_H_
#define _RECURSION_H_

#include <clue_io.h>

#include <largeN_QFT_parameters.h>
#include "parameters.h"

typedef unsigned int uint;

void getG2(double* G2, int m);

void testG4(int m);
void get_momentum_str(uint P, int n, char** s, char* separator);
uint   total_momentum(uint P, int n);

void print_buffer_usage();

void init_recursion();
void free_recursion();

#endif
