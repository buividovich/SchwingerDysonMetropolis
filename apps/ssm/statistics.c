#include "statistics.h"

double                     astop  = 0.0;
int*                   G2_hist[2] = {NULL, NULL};

void init_observable_stat()
{
 int s, i;
 for(s=0; s<2; s++)
 {
  SAFE_MALLOC(G2_hist[s], int, lat_vol);
  for(i=0; i<lat_vol; i++)
   G2_hist[s][i] = 0;
 };
}

void gather_observable_stat()
{
 int i;
 if(X.len[X.top-1]==2)
 {
  i = lat_coords2idx(STACK_EL(X,0));
  G2_hist[(asign[ns]>0? 0 : 1)][i] ++;
 };
}

void process_observable_stat() //It is assumed that process_mc_stat was already called!!!
{
 int i;
 FILE* ofile = NULL;
 char cstr[30];
 
 double source_norm          = action_create_amplitude(NULL);
 double normalization_factor = source_norm/(1.0 - mean_nA);
 logs_Write(0, "Normalization factor of multiple-trace correlators: %2.4E\n", normalization_factor);
 normalization_factor = normalization_factor/(1.0 + mean_sign*normalization_factor);
 logs_Write(0, "Normalization factor of factorized single-trace correlators: %2.4E\n", normalization_factor);
 
 if(observables_file!=NULL)
 { 
  ofile = fopen(observables_file, "w");
  if(ofile!=NULL)
   fprintf(ofile, "%2.4E ", lambda);
  else
   logs_WriteError("Could not open the file %s for writing", observables_file); 
 };
 
 //Saving the data for G2 
 double rescaling_factor     = NN*cc;
 
 for(i=0; i<lat_vol; i++)
 {   
  double G  =      (double)(G2_hist[0][i] - G2_hist[1][i])/(double)nmc;
  double dG = sqrt((double)(G2_hist[0][i] + G2_hist[1][i]))/(double)nmc;
    
  G  *= rescaling_factor*normalization_factor;
  dG *= rescaling_factor*normalization_factor;
   
  //SP characterizes the strength of the sign problem
  double SP = 0.0;
  if(abs(G2_hist[0][i] + G2_hist[1][i])!=0)
   SP = (double)(G2_hist[0][i] - G2_hist[1][i])/(double)(G2_hist[0][i] + G2_hist[1][i]);
   
  sprintf_coords(cstr, i);
   
  if(ofile!=NULL)
   fprintf(ofile, "%s %2.4E %2.4E %2.4E ", cstr, G, dG, SP);
     
  logs_Write(0, "G%s:\t %2.4E +/- %2.4E,\t sp = %2.4E", cstr, G, dG, SP);
 };
 
 if(ofile!=NULL)
 {
  fprintf(ofile, "\n"); 
  fclose(ofile);   
 };
}

void free_observable_stat()
{
 int s;
 for(s=0; s<2; s++)
  SAFE_FREE(G2_hist[s]);
}

