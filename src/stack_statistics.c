#include "stack_statistics.h"

t_stack_stat* init_stack_statistics(int max_hist_size)
{
 t_stack_stat* stat = (t_stack_stat* )malloc(sizeof(t_stack_stat));
 SAFE_MALLOC(stat->seq_len_hist, int, max_hist_size);
 SAFE_MALLOC(stat->num_seq_hist, int, max_hist_size);
 stat->my_max_hist_size = max_hist_size;
 for(int i=0; i<max_hist_size; i++)
 {
  stat->seq_len_hist[i] = 0;
  stat->num_seq_hist[i] = 0;
 };
 stat->nstat = 0; 
 return stat;
}

void free_stack_statistics(t_stack_stat* stat)
{
 SAFE_FREE(stat->seq_len_hist);
 SAFE_FREE(stat->num_seq_hist);   
}

void gather_stack_statistics(t_stack_stat* stat, t_lat_stack* X)
{
 int len = X->len[X->top-1]/2;
 if(len>0 && len<stat->my_max_hist_size)
  stat->seq_len_hist[len-1]    ++;
 if(X->top>0 && X->top<stat->my_max_hist_size)
  stat->num_seq_hist[X->top-1] ++;
 stat->nstat ++;
}

void print_stack_statistics(t_stack_stat* stat)
{
 double mean_seq_len = 0.0;
 for(int ilen=0; ilen<stat->my_max_hist_size; ilen++)
  mean_seq_len += (double)(ilen+1)*(double)(stat->seq_len_hist[ilen]);
 mean_seq_len /= (double)(stat->nstat);
 
 double mean_num_seq = 0.0;
 for(int imns=0; imns<stat->my_max_hist_size; imns++)
  mean_num_seq += (double)(imns+1)*(double)(stat->num_seq_hist[imns]);
 mean_num_seq /= (double)(stat->nstat);
 
 logs_Write(0, "STATISTICS ON STACK USAGE:");
 logs_WriteParameter(0, "Mean sequence length",              "%2.2lf", mean_seq_len);
 logs_WriteParameter(0, "Mean number of sequences in stack", "%2.2lf", mean_num_seq);
 logs_Write(0, "");
}

void print_stack_histogram(t_stack_stat* stat, FILE* f)
{
 if(f==stdout)
  logs_Write(0, "\nHistogram of sequence lengths (over %i instances):", stat->nstat);

 for(int i=0; i<stat->my_max_hist_size-1; i++)
 {
  double p1  =      stat->seq_len_hist[i]   /(double)(stat->nstat);
  double dp1 = sqrt(stat->seq_len_hist[i]  )/(double)(stat->nstat);
  double p2  =      stat->seq_len_hist[i+1] /(double)(stat->nstat);
  double dp2 = sqrt(stat->seq_len_hist[i+1])/(double)(stat->nstat);
  
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

double mean_nA_prediction(t_stack_stat* stat, double NN_new, double cc_new, double source_norm)
{
 double source_norm_new      = cc/cc_new*source_norm;
 double normalization_factor = source_norm/(1.0 - mean_nA);
 
 double mean_cc_ratio        = 0;
 double fact                 = cc/cc_new;
 for(int i=0; i<stat->my_max_hist_size-1; i++)
 {
  mean_cc_ratio += fact*stat->seq_len_hist[i];
  fact *= cc/cc_new;
 }; 
 mean_cc_ratio = mean_cc_ratio/(double)(stat->nstat);
 
 double eta = normalization_factor/(1.0 + normalization_factor)*mean_cc_ratio;
 
 return 1.0 - source_norm_new*(1.0/eta - NN/NN_new);
}
