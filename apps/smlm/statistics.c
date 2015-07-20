#include "statistics.h"

double                     astop  = 0.0;
int*                   G2_hist[2] = {NULL, NULL};
int*                    G_hist[2] = {NULL, NULL};
int                  max_G2_order = 0;

void init_observable_stat()
{
 max_G2_order = 0;

 for(int s=0; s<2; s++)
 {
  SAFE_MALLOC(G2_hist[s], int, lat_vol);
  for(int i=0; i<lat_vol; i++)
   G2_hist[s][i] = 0;
  SAFE_MALLOC(G_hist[s], int, max_correlator_order);
  for(int i=0; i<max_correlator_order; i++)
   G_hist[s][i] = 0; 
 };
 
 init_stack_statistics(max_stack_nel);
}

void gather_observable_stat()
{   
 if(X.len[X.top-1]==2)
 {
  int m = lat_coords2idx_safe(STACK_EL(X,0));
  G2_hist[(asign[ns]>0? 0 : 1)][m] ++;
 };
 int i = X.len[X.top-1]/2 - 1;
 if(i>=0 && i<max_correlator_order)
  G_hist[(asign[ns]>0? 0 : 1)][i] ++;
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

 //Saving the data for G2 
 logs_Write(0, "Collected data for G2:");
 double G2_total = 0.0, G2_total_err = 0.0;
 double G2_link  = 0.0, G2_link_err  = 0.0;
 
 double rescaling_factor     = NN*cc;

 for(int i=0; i<lat_vol; i++)
 {   
  int m[4];
  lat_idx2coords(i, m); 
  
  double G  =      (double)(G2_hist[0][i] - G2_hist[1][i]) /(double)nmc;
  double dG = sqrt((double)(G2_hist[0][i] + G2_hist[1][i]))/(double)nmc;
    
  G  *= rescaling_factor*normalization_factor;
  dG *= rescaling_factor*normalization_factor;
  
  G2_total     += G;
  G2_total_err += SQR(dG);
 
  for(int mu=0; mu<DIM; mu++)
  {
   double cf = cos(2.0*M_PI*(double)m[mu]/(double)lat_size[mu]);
   G2_link      += G*cf;
   G2_link_err  += SQR(dG*cf);  
  }; 
 };
 
 G2_total_err = sqrt(G2_total_err);
 
 G2_link      = G2_link/(double)DIM;
 G2_link_err  = sqrt(G2_link_err/(double)DIM);
 
 logs_Write(0, " G2 TOTAL: %2.4E +/- %2.4E", G2_total, G2_total_err);
 logs_Write(0, " G2  LINK: %2.4E +/- %2.4E",  G2_link,  G2_link_err);
 logs_Write(0, " max. G2 order: %i", max_G2_order);
 
 logs_Write(0, "\n\t Total G for different orders:");
 for(int i=0; i<max_correlator_order; i++)
 {
  rescaling_factor = NN*pow(cc, (double)(i + 1)); 
  double G  =      (double)(G_hist[0][i] - G_hist[1][i]) /(double)nmc;
  double dG = sqrt((double)(G_hist[0][i] + G_hist[1][i]))/(double)nmc;
  G  *= rescaling_factor*normalization_factor;
  dG *= rescaling_factor*normalization_factor;
  logs_Write(0, "\t G%02i:\t %2.4E +/- %2.4E", 2*i + 2, G, dG);
 };
 logs_Write(0, "");
 
 if(observables_file!=NULL)
 {
  FILE* f = fopen(observables_file, "a");
  if(f==NULL)
   logs_WriteError("Cannot open the file %s for writing", observables_file);
  else
   fprintf(f, "%2.4E %2.4E %2.4E %2.4E %2.4E\n", lambda, G2_link, G2_link_err, G2_total, G2_total_err);
  fclose(f);
 };
}

void free_observable_stat()
{
 free_stack_statistics();
 for(int s=0; s<2; s++)
 {
  SAFE_FREE(G2_hist[s]);
  SAFE_FREE(G_hist[s]);
 };  
}

