#include "statistics.h"

double                     astop  = 0.0;
int                   *G_hist[2]  = {NULL, NULL};

void init_observable_stat()
{
 int s, i;
 for(s=0; s<2; s++)
 {
  SAFE_MALLOC(G_hist[s], int, max_correlator_order);
  for(i=0; i<max_correlator_order; i++)
   G_hist[s][i] = 0;
 };
}

void gather_observable_stat()
{
 int gn = X.len[X.top-1]; 
 ASSERT(gn<=0);
 if(gn<=max_correlator_order)
  G_hist[(asign[ns]>0? 0 : 1)][gn-1]++;
}

void process_observable_stat() //It is assumed that process_mc_stat was already called!!!
{
 double source_norm          = action_create_amplitude(NULL);
 double normalization_factor = source_norm/(1.0 - mean_nA);
 logs_Write(0, "Normalization factor of multiple-trace correlators: %2.4E\n", normalization_factor);
 normalization_factor = normalization_factor/(1.0 + mean_sign*normalization_factor);
 logs_Write(0, "Normalization factor of factorized single-trace correlators: %2.4E\n", normalization_factor);

 char* observables_datastr = NULL;
 char* gtotal_datastr      = NULL;
 sprintf_append(&observables_datastr, "%2.4E %2.4E %2.4E ", lambda, cc, NN);
  
 for(int ign=0; ign<max_correlator_order; ign++)
 {
  double rescaling_factor     = NN*pow(cc,(double)(ign+1));
    
  double G  = (double)(G_hist[0][ign] - G_hist[1][ign])/(double)nmc;
  double dG = sqrt((double)(G_hist[0][ign] + G_hist[1][ign]))/(double)nmc;
  double GT  = (double)(G_hist[0][ign] + G_hist[1][ign])/(double)nmc;
    
  G  *= rescaling_factor*normalization_factor;
  GT *= rescaling_factor*normalization_factor;
  dG *= rescaling_factor*normalization_factor;
    
  //SP characterizes the strength of the sign problem
  double SP = (double)(G_hist[0][ign] - G_hist[1][ign])/(double)(G_hist[0][ign] + G_hist[1][ign]);
  
  sprintf_append(&observables_datastr, "%2.4E %2.4E %2.4E ",  G, dG, SP);
  sprintf_append(&gtotal_datastr,      "%i %2.4E %2.4E\n",    ign+1, GT, dG);
  
  logs_Write(0, "G_%i:\t %2.4E +/- %2.4E,\t sp = %+2.4E, \t GT = %2.4E +/- %2.4E", ign+1, G, dG, SP, GT, dG);
 };
 if(observables_file!=NULL)
  safe_append_str_to_file(observables_file, observables_datastr, io_sleep_time, io_write_attempts);
 if(gtotal_file!=NULL)
  safe_append_str_to_file(     gtotal_file,      gtotal_datastr, io_sleep_time, io_write_attempts);
 
 SAFE_FREE(observables_datastr);
 SAFE_FREE(gtotal_datastr); 
}

void free_observable_stat()
{
 int s;
 for(s=0; s<2; s++)
  SAFE_FREE(G_hist[s]);
}

