#ifndef _RECURSION_H_
#define _RECURSION_H_

#include <clue_io.h>

#include <largeN_QFT_parameters.h>
#include "common_parameters.h"

typedef unsigned int uint;

//General functions for working with sequences of momenta 
void  get_momentum_str(uint P, int n, char** s, char* separator);
uint    total_momentum(uint P, int n);
void    unpack_momenta(uint P, int n, uint* ps);
uint      pack_momenta(int n, uint* ps);

//General functions for recursions of any kind
void    test_recursion();
void    free_recursion();

//Generalized strong-coupling expansion
void    init_generalized_sc();
void     run_generalized_sc();
void get_Gxy_generalized_sc(double* res);

#endif
