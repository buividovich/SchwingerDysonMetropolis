#ifndef _RECURSION_H_
#define _RECURSION_H_

#include <clue_io.h>

#include <largeN_QFT_parameters.h>
#include "common_parameters.h"

typedef long double  t_real;
#define ACALL(_f)  _f ## l

typedef unsigned int  uint;

//General functions for working with sequences of momenta 
void  get_momentum_str(uint P, int n, char* s, const char* separator);
void  get_momentum_str_append(uint P, int n, char** s, char* separator);
uint      pack_momenta(int n, uint* ps);

void printf_momentum_str(uint P, int n, char* prefix, char* separator, char* suffix);

void alloc_recursion_packed();
void init_lattice_constants();

//General functions for recursions of any kind
void    test_recursion();
void    free_recursion();

void    print_Gmn();

//Generalized strong-coupling expansion
void    init_generalized_sc();
void     run_generalized_sc();
void get_Gxy_generalized_sc(double* res);

//Stereographic projection
void get_Gxy_stereographic(double* Gxy, double* Gx);
void run_stereographic();
void init_stereographic();

#endif
