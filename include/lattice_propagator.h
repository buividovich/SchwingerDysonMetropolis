#ifndef _LATTICE_PROPAGATOR_H_
#define _LATTICE_PROPAGATOR_H_

#include <square_lattice.h>
#include <rand_num_generators.h>

#include "largeN_QFT_parameters.h"
#include "lattice_stack.h"

//Lattice propagator structure
typedef struct 
{
 int     allocated;     //=1 if all the memory has already been allocated
 double* probabilities;
 double  sigma; // = sum of propagator over all momenta
 double  mass_sq; //Squared mass
} t_lat_propagator;

void    init_lat_propagator(t_lat_propagator* P, int allocate, double mass_sq);
void    free_lat_propagator(t_lat_propagator* P);

double  lat_momentum_sq(int *m);
double  lat_propagator(int* m, double mass_sq);

void    rand_momentum(t_lat_propagator* P, int* m);

void      zero_momentum(int* m);
void    assign_momentum(int* m_out, int* m_in);
void    invert_momentum(int* m_out, int* m_in);

void    addto_momentum(int* m, int s1, int* m1);                  //m += s1*m1
void     addto2momenta(int* m, int s1, int* m1, int s2, int* m2); //m += s1*m1 + s2*m2

void    add2momenta(int* m, int* m1, int* m2);          //m = m1 + m2
void    add3momenta(int* m, int* m1, int* m2, int* m3); //m = m1 + m2 + m3

int     check_momentum_conservation(t_lat_stack* lat_stack, const char* stack_name);

#endif
