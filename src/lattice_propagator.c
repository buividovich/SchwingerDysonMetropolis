#include "lattice_propagator.h"

void init_lat_propagator(t_lat_propagator* P, double mass_sq, double norm_factor)
{
 if(!lat_initialized_flag)
  init_lattice();
 //Calculating the array of probabilities in order to initialize fast random choice
 DECLARE_AND_MALLOC(probs, double, lat_vol);
 //At the same time, calculating sigma
 P->mass_sq = mass_sq;
 P->sigma   = 0.0;
 for(int i=0; i<lat_vol; i++)
 {
  int m[4];
  lat_idx2coords(i, m);
  probs[i] = fabs(lat_propagator(m, mass_sq, norm_factor));
  P->sigma += probs[i];
 };
 //Finally, normalizing the probabilities to unity
 for(int i=0; i<lat_vol; i++)
  probs[i] /= P->sigma;
 //Initializing the fast random choice algorithm
 SAFE_MALLOC(P->frc_data, t_frc_data, 1); 
 init_fast_rand_choice(probs, lat_vol, P->frc_data);
 //And the sigma itself should be normalized by volume
 P->sigma /= (double)lat_vol; 
 
 logs_Write(2, "init_lat_propagator: P->sigma = %2.4E, lat_vol = %i", P->sigma, lat_vol);
 
 P->allocated = 1; 
 SAFE_FREE(probs);
}

void free_lat_propagator(t_lat_propagator* P)
{
 lat_free();
 free_fast_rand_choice(P->frc_data);
 SAFE_FREE(P->frc_data);
 P->allocated = 0;
}

double  lat_momentum_sq(int *m)
{
 double res = 0.0;
 for(int mu=0; mu<DIM; mu++)
  res += 4.0*SQR(sin( M_PI*(double)(m[mu])/(double)lat_size[mu] ));
 return res; 
}

double  lat_propagator(int* m, double mass_sq, double norm_factor)
{
 return norm_factor/(mass_sq + lat_momentum_sq(m)); 
}

void rand_momentum(t_lat_propagator* P, int* m)
{
 //printf("A call to rand_momentum...\n");
 //fflush(stdout);
 int i = fast_rand_choice(P->frc_data);
 //printf("After fast_rand_choice...\n");
 //fflush(stdout);
 lat_idx2coords(i, m);
 //printf("After lat_idx2coords...\n");
 //fflush(stdout);
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
 for(int mu=0; mu<DIM; mu++)
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

int   check_momentum_conservation(t_lat_momentum_stack* lat_stack, const char* stack_name)
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

