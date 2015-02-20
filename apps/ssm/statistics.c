#include "statistics.h"

double                     astop  = 0.0;
int*                   G2_hist[2] = {NULL, NULL};
int*                   G4_hist[2] = {NULL, NULL};
int                  max_G2_order = 0;
int                  max_G4_order = 0;
int                      lat_vol2 = 0;
int                      lat_vol3 = 0;
int                      lat_vol4 = 0;

void init_observable_stat()
{
 int s, i;
 lat_vol2 = lat_vol*lat_vol;
 lat_vol3 = lat_vol2*lat_vol;
 lat_vol4 = lat_vol3*lat_vol;
 
 max_G2_order = 0;
 max_G4_order = 0;

 for(s=0; s<2; s++)
 {
  SAFE_MALLOC(G2_hist[s], int, lat_vol);
  SAFE_MALLOC(G4_hist[s], int, lat_vol4);
  for(i=0; i<lat_vol; i++)
   G2_hist[s][i] = 0;
  for(i=0; i<lat_vol4; i++)
   G4_hist[s][i] = 0;
 };
}

void gather_observable_stat()
{
 if(X.len[X.top-1]==2)
  max_G2_order = MAX(O[X.top-1], max_G2_order);
 if(X.len[X.top-1]==4)
  max_G4_order = MAX(O[X.top-1], max_G4_order);
   
 if(O[X.top-1]<=max_observables_order && O[X.top-1]>=min_observables_order)
 {
  if(X.len[X.top-1]==2)
  {
   int m = lat_coords2idx_safe(STACK_EL(X,0));
   G2_hist[(asign[ns]>0? 0 : 1)][m] ++;
  };
  if(X.len[X.top-1]==4)
  {
   int m =          lat_coords2idx_safe(STACK_EL(X,0)) + 
           lat_vol *lat_coords2idx_safe(STACK_EL(X,1)) +
           lat_vol2*lat_coords2idx_safe(STACK_EL(X,2)) +
           lat_vol3*lat_coords2idx_safe(STACK_EL(X,3)); 
   G4_hist[(asign[ns]>0? 0 : 1)][m] ++;
  };
 }; 
}

void process_observable_stat() //It is assumed that process_mc_stat was already called!!!
{
 char   *outstr = NULL;
 char   cstr[20];
 double rescaling_factor, order0_err;
 
 double source_norm          = action_create_amplitude(NULL);
 double normalization_factor = source_norm/(1.0 - mean_nA);
 logs_Write(0, "Normalization factor of multiple-trace correlators: %2.4E\n", normalization_factor);
 normalization_factor = normalization_factor/(1.0 + mean_sign*normalization_factor);
 logs_Write(0, "Normalization factor of factorized single-trace correlators: %2.4E\n", normalization_factor);
 
 sprintf_append(&outstr, "%2.4E ", lambda);
 
 //Saving the data for G2 
 logs_Write(0, "Collected data for G2:");
 rescaling_factor     = NN*cc;
 order0_err = 0.0;
 for(int i=0; i<lat_vol; i++)
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
  sprintf_append(&outstr, "%s %2.4E %2.4E %2.4E ", cstr, G, dG, SP);
   
  order0_err += SQR(G - 0.5); 
     
  logs_Write(0, "\t G%s:\t %+2.4E +/- %+2.4E,\t sp = %2.4E, \t max. order = %i", cstr, G, dG, SP, max_G2_order);
 };
 sprintf_append(&outstr, "\n");
  
 logs_Write(0, "Overall error for G2: %2.4E", sqrt(order0_err));
 
 if(observables_file!=NULL)
 {
  int res = safe_append_str_to_file(observables_file, outstr);
  if(res!=0)
   logs_WriteError("safe_append_str_to_file %s failed with code %i", observables_file, res);
 }; 
 
 //Printing out the data for G4
 logs_Write(0, "\nCollected data for G4:");
 rescaling_factor     = NN*SQR(cc);
 order0_err           = 0.0;
 for(int i=0; i<lat_vol4; i++)
 {   
  double G  =      (double)(G4_hist[0][i] - G4_hist[1][i])/(double)nmc;
  double dG = sqrt((double)(G4_hist[0][i] + G4_hist[1][i]))/(double)nmc;
    
  G  *= rescaling_factor*normalization_factor;
  dG *= rescaling_factor*normalization_factor;
   
  //SP characterizes the strength of the sign problem
  double SP = 0.0;
  if(abs(G4_hist[0][i] + G4_hist[1][i])!=0)
   SP = (double)(G4_hist[0][i] - G4_hist[1][i])/(double)(G4_hist[0][i] + G4_hist[1][i]);
  
  int tmp = i; 
  int p1 = tmp%lat_vol; tmp = tmp/lat_vol;
  int q1 = tmp%lat_vol; tmp = tmp/lat_vol;
  int p2 = tmp%lat_vol; tmp = tmp/lat_vol;
  int q2 = tmp%lat_vol;    
  
  #define DELTA(_x) ( (_x)%lat_vol==0? 1.0 : 0.0 )
  double sbG4 = (DELTA(p1 + q1)*DELTA(p2 + q2) + DELTA(p1 + q2)*DELTA(p2 + q1))/(double)lat_vol2 - 1.0/(double)lat_vol3;
  #undef DELTA
    
  if((p1 + q1 + p2 + q2)%lat_vol==0)
  {
   logs_Write(0, "\tG[%i%i%i%i]:\t %+2.4E +/- %+2.4E, should be %+2.4E, \t SP = %+2.4E, \t diff = %+2.4E, \t max.order = %i", p1, q1, p2, q2, G, dG, sbG4, SP, SQR(G - sbG4), max_G4_order);
   order0_err += SQR(G - sbG4); 
  }
  else
  {
   if(fabs(G)>0.0) 
    logs_WriteError("Nonzero result for G4 at (p1 + q1 + p2 + q2)!=0!!! p1 = %i, q1 = %i, p2 = %i, q2 = %i, G = %2.4E +/- %2.4E", p1, q1, p2, q2, G, dG);
  };  
 };
 logs_Write(0, "Overall error for G4: %2.4E", order0_err);
 
 SAFE_FREE(outstr);
}

void free_observable_stat()
{
 int s;
 for(s=0; s<2; s++)
  SAFE_FREE(G2_hist[s]);
}

