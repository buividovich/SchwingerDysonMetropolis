#include "hmm_statistics.h"

double                     astop  = 0.0;
int                   *G_hist[2]  = {NULL, NULL};

void init_observable_stat()
{
 int s;
 for(s=0; s<2; s++)
 {
  SAFE_MALLOC(G_hist[s], int, max_correlator_order);
  set_zero_int(G_hist[s], max_correlator_order);
 };
}

void gather_observable_stat()
{
 int gn = 0; //TODO: modify here!!!
 ASSERT(gn<=0);
 if(gn<=max_correlator_order)
  G_hist[(asign[ns]>0? 0 : 1)][gn-1]++;
}

void process_observable_stat()
{
 int ign;
 if(observables_file!=NULL)
 {
  double source_norm          = 1.0/NN;
  double normalization_factor = source_norm/(1.0 - anA/(double)prod_mc_steps);
  logs_Write(0, "Normalization factor of multiple-trace correlators: %2.4E\n", normalization_factor);
  normalization_factor = normalization_factor/(1.0 + (mean_sign/(double)nmc)*normalization_factor);
  logs_Write(0, "Normalization factor of factorized single-trace correlators: %2.4E\n", normalization_factor);
  
  FILE *ofile = fopen(observables_file, "w");
  for(ign=0; ign<max_correlator_order; ign++)
   if(G_hist[0][ign]>0 || G_hist[1][ign]>0)
   {
    double rescaling_factor     = NN*pow(cc,(double)(ign));
    
    double G  = (double)(G_hist[0][ign] - G_hist[1][ign])/(double)nmc;
    double dG = sqrt((double)(G_hist[0][ign] + G_hist[1][ign]))/(double)nmc;
    
    G  *= rescaling_factor*normalization_factor;
    dG *= rescaling_factor*normalization_factor;
    
    //SP characterizes the strength of the sign problem
    double SP = (double)(G_hist[0][ign] - G_hist[1][ign])/(double)(G_hist[0][ign] + G_hist[1][ign]);
    fprintf(ofile, "%i %2.4E %2.4E %2.4E\n", ign+1, G, dG, SP);
   }; 
  fclose(ofile);   
 };
}

void free_observable_stat()
{
 int s;
 for(s=0; s<2; s++)
  SAFE_FREE(G_hist[s]);
}

