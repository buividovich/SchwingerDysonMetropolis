#include "stack_statistics.h"

double* seq_len_hist     = NULL;
double* num_seq_hist     = NULL;
int     nstat            = 0;
int     my_max_hist_size = 0;

void init_stack_statistics(int max_hist_size)
{
 SAFE_MALLOC(seq_len_hist, double, max_hist_size);
 SAFE_MALLOC(num_seq_hist, double, max_hist_size);
 my_max_hist_size = max_hist_size;
 for(int i=0; i<max_hist_size; i++)
 {
  seq_len_hist[i] = 0.0;
  num_seq_hist[i] = 0.0;
 };
 nstat = 0; 
}

void free_stack_statistics()
{
 SAFE_FREE(seq_len_hist);
 SAFE_FREE(num_seq_hist);   
}

void gather_stack_statistics(t_lat_stack* X)
{
 int len = X->len[X->top-1]/2;
 if(len>0 && len<my_max_hist_size)
  seq_len_hist[len-1]    += 1.0;
 if(X->top>0 && X->top<my_max_hist_size)
  num_seq_hist[X->top-1] += 1.0;
 nstat ++;
}

void print_stack_statistics(FILE* f)
{
 if(f==stdout)
  logs_Write(0, "\nHistogram of sequence lengths (over %i instances):", nstat);

 for(int i=0; i<my_max_hist_size-1; i++)
 {
  double p1  = seq_len_hist[i]/(double)nstat;
  double dp1 = sqrt(seq_len_hist[i])/(double)nstat;
  double p2  = seq_len_hist[i+1]/(double)nstat;
  double dp2 = sqrt(seq_len_hist[i+1])/(double)nstat;
  
  if(p1>0.0 && p2>0.0 && dp1/p1<0.1 && dp2/p2<0.1)
  {
   double r  = p2/p1;
   double dr = sqrt(SQR( dp2/p1 ) + SQR( p2*dp1/(p1*p1) )); 
   fprintf(f, "%03i %2.4E %2.4E %2.4E %2.4E %2.4E %2.4E %2.4E\n", 2*i, cc, p1, dp1, r, dr, cc*r, cc*dr);
  }
  else
  {
   logs_Write(0, "Histogram of sequence lengths contains large errors at i=%i, p1=%2.4E +/- %2.4E, p2=%2.4E +/- %2.4E, terminating output", i, p1, dp1, p2, dp2);
   break; 
  }; 
 }; 

 if(f==stdout)
  logs_Write(0, " ");
}

double mean_nA_prediction(double NN_new, double cc_new, double source_norm)
{
 double source_norm_new      = cc/cc_new*source_norm;
 double normalization_factor = source_norm/(1.0 - mean_nA);
 
 double mean_cc_ratio        = 0;
 double fact                 = cc/cc_new;
 for(int i=0; i<my_max_hist_size-1; i++)
 {
  mean_cc_ratio += fact*seq_len_hist[i];
  fact *= cc/cc_new;
 }; 
 mean_cc_ratio = mean_cc_ratio/(double)nstat;
 
 double eta = normalization_factor/(1.0 + normalization_factor)*mean_cc_ratio;
 
 return 1.0 - source_norm_new*(1.0/eta - NN/NN_new);
}
