#include "lattice_propagator.h"

void init_lat_propagator(t_lat_propagator* P, int allocate, double mass_sq)
{
 int my_lat_sizes[26], m[26], i; //26 is the largest space-time dimensionality I can think of - the critical dimensionality of a bosonic string
 if(allocate)
 {
  //Initializing the square_lattice unit
  my_lat_sizes[0] = LT;
  for(i=1; i<DIM; i++)
   my_lat_sizes[i] = LS;
  lat_init(DIM, my_lat_sizes);
  SAFE_MALLOC(P->probabilities, double, lat_vol);
  P->allocated = 1;
 };
 //Saving propagator values for all lattice momenta - to speed up the generation of random momenta
 //At the same time, calculating sigma
 P->mass_sq = mass_sq;
 P->sigma   = 0.0;
 for(i=0; i<lat_vol; i++)
 {
  lat_idx2coords(i, m);
  P->probabilities[i] = fabs(lat_propagator(m, mass_sq));
  P->sigma += P->probabilities[i];
 };
 //Finally, normalizing the probabilities to unity
 for(i=0; i<lat_vol; i++)
  P->probabilities[i] /= P->sigma;
 //And the sigma itself should be normalized by volume
 P->sigma /= (double)lat_vol; 
 logs_Write(2, "init_lat_propagator: P->sigma = %2.4E, lat_vol = %i", P->sigma, lat_vol);
}

void free_lat_propagator(t_lat_propagator* P)
{
 lat_free();
 SAFE_FREE(P->probabilities);
 P->allocated = 0;
}

double  lat_momentum_sq(int *m)
{
 int mu;
 double res = 0.0;
 for(mu=0; mu<DIM; mu++)
  res += 4.0*SQR(sin( M_PI*(double)(m[mu])/(double)lat_size[mu] ));
 return res; 
}

double  lat_propagator(int* m, double mass_sq)
{
 return 1.0/(mass_sq + lat_momentum_sq(m)); 
}

void rand_momentum(t_lat_propagator* P, int* m)
{
 int i = rand_choice(P->probabilities, lat_vol);
 lat_idx2coords(i, m);    
}

void      zero_momentum(int* m)
{
 int mu;
 for(mu=0; mu<DIM; mu++)
  m[mu] = 0;
}

void    assign_momentum(int* m_out, int* m_in)
{
 int mu;
 for(mu=0; mu<DIM; mu++)
  m_out[mu] = m_in[mu];
}

void    invert_momentum(int* m_out, int* m_in)
{
 int mu;
 for(mu=0; mu<DIM; mu++)
  m_out[mu] = -m_in[mu];
}

void    addto_momentum(int* m, int s1, int* m1)                  //m += s1*m1
{
 int mu;
 for(mu=0; mu<DIM; mu++)
  m[mu] += s1*m1[mu];
}

void    addto2momenta(int* m, int s1, int* m1, int s2, int* m2) //m += m1 + m2
{
 int mu;
 for(mu=0; mu<DIM; mu++)
  m[mu] += s1*m1[mu] + s2*m2[mu];
}

void    add2momenta(int* m, int* m1, int* m2)
{
 int mu;
 for(mu=0; mu<DIM; mu++)
  m[mu] = m1[mu] + m2[mu];
}

void    add3momenta(int* m, int* m1, int* m2, int* m3)
{
 int mu;
 for(mu=0; mu<DIM; mu++)
  m[mu] = m1[mu] + m2[mu] + m3[mu];
}

int   check_momentum_conservation(t_lat_stack* lat_stack, const char* stack_name)
{
 int iseq, i, mu, res_seq, res = 1;
 
 DECLARE_AND_MALLOC(tm, int, lat_stack->dim);
 
 for(iseq=0; iseq<lat_stack->top; iseq++)
 {
  zero_momentum(tm);
  for(i=0; i<lat_stack->len[iseq]; i++)
   addto_momentum(tm, +1, lat_stack->stack[lat_stack->start[iseq] + i]);
  res_seq = 1; 
  for(mu=0; mu<lat_stack->dim; mu++)
   res_seq = res_seq && (tm[mu]%lat_size[mu]==0);
  if(!res_seq)
  {
   logs_WriteError("Momentum conservation violated in sequence %i of stack %s", iseq, stack_name);
   for(mu=0; mu<DIM; mu++)
    logs_Write(0, "tm[%i] = %i, lat_size[%i] = %i", mu, tm[mu], mu, lat_size[mu]);
  };
  res = res && res_seq;  
 };
 SAFE_FREE(tm);
 
 return res;
}

