#ifndef _RECURSION_H_
#define _RECURSION_H_

#include <largeN_QFT_parameters.h>
#include "parameters.h"

typedef unsigned int uint;

#define  MAX_M   (32)    //Maximal orders - to define storage of pre-computed values as static array
#define  MAX_MX2 (64)    //Maximal orders - to define storage of pre-computed values as static array
#define  MAXLS   (16)    //Maximal LS  - for static arrays
#define  MAXLS2  (256)   //Maximal LS2 - for static arrays containing bare propagators etc.

double G(uint P, int n, int m);

void init_recursion();
void free_recursion();

#endif
