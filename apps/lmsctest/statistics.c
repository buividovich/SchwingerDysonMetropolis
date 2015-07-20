#include "statistics.h"

double                     astop  = 0.0;
int*                    G_hist    = NULL;

void init_observable_stat()
{
 SAFE_MALLOC(G_hist, int, max_correlator_order);
 for(int i=0; i<max_correlator_order; i++)
  G_hist[i] = 0;
 init_stack_statistics(max_stack_nel);
}

void gather_observable_stat()
{   
 if(X.len[X.top-1]<max_correlator_order)
  G_hist[X.len[X.top-1]] ++;
 gather_stack_statistics(&X); 
}

void process_observable_stat() //It is assumed that process_mc_stat was already called!!!
{
 if(stack_stat_file!=NULL)
 {
  FILE* f = fopen(stack_stat_file, "w");
  if(f==NULL)
   logs_WriteError("Cannot open the file %s for writing", stack_stat_file);
  else
   print_stack_statistics(f);
  fclose(f);
 };
 
 if(logs_noise_level>=2)
  print_stack_statistics(stdout);
 
 double source_norm          = action_create_amplitude(NULL);
 double normalization_factor = source_norm/(1.0 - mean_nA);
 logs_Write(0, "Normalization factor of multiple-trace correlators: %2.4E\n", normalization_factor);
 normalization_factor = normalization_factor/(1.0 + mean_sign*normalization_factor);
 logs_Write(0, "Normalization factor of factorized single-trace correlators: %2.4E\n", normalization_factor);

 logs_Write(0, "Collected data for G:");
 for(int i=0; i<max_correlator_order; i++)
 {
  double rescaling_factor     = NN*pow(cc, (double)i);
  double G  =      (double)(G_hist[i]) /(double)nmc*rescaling_factor;
  double dG = sqrt((double)(G_hist[i]))/(double)nmc*rescaling_factor;
  logs_Write(0, "\t G%01i:\t %2.4E +/- %2.4E", i, G, dG);
 };

}

void free_observable_stat()
{
 free_stack_statistics();
 SAFE_FREE(G_hist);
}

