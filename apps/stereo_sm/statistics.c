#include "statistics.h"

t_observable_stat* init_observable_stat()
{
 t_observable_stat* my_observable_stat = (t_observable_stat *)malloc(sizeof(t_observable_stat));

 for(int s=0; s<2; s++)
 {
  SAFE_MALLOC(my_observable_stat->G2_hist[s], int, lat_vol*max_alpha_order);
  for(int i=0; i<lat_vol*max_alpha_order; i++)
   my_observable_stat->G2_hist[s][i] = 0;
  
  my_observable_stat->nG4 = SQR(lat_vol)*SQR(lat_vol);
  SAFE_MALLOC(my_observable_stat->G4_hist[s], int, my_observable_stat->nG4*max_alpha_order);
  for(int i=0; i<my_observable_stat->nG4*max_alpha_order; i++)
   my_observable_stat->G4_hist[s][i] = 0; 
   
  SAFE_MALLOC(my_observable_stat->G_hist[s], int, max_correlator_order);
  for(int i=0; i<max_correlator_order; i++)
   my_observable_stat->G_hist[s][i] = 0; 
 };
 
 my_observable_stat->nstat = 0;
 my_observable_stat->actual_max_alpha_order = 0;
 
 return my_observable_stat;
}

void gather_observable_stat(t_observable_stat* stat)
{
 int ao = alpha_order_stack[alpha_order_stop-1];
 int si =        sign_stack[       sign_stop-1];
 if(X.len[X.top-1]==2 && ao<max_alpha_order)
 {
  int m = lat_coords2idx_safe(STACK_EL(X,0));
  stat->G2_hist[(si>0? 0 : 1)][lat_vol*ao + m] ++;
  stat->actual_max_alpha_order = MAX(ao, stat->actual_max_alpha_order);
 };
 
 if(X.len[X.top-1]==4 && ao<max_alpha_order)
 {
  int mall = 0; int factor = 1;
  for(int i=0; i<4; i++)
  {
   int m = lat_coords2idx_safe(STACK_EL(X,i));
   mall += m*factor;
   factor *= lat_vol;
  }; 
  
  stat->G4_hist[(si>0? 0 : 1)][stat->nG4*ao + mall] ++;
  stat->actual_max_alpha_order = MAX(ao, stat->actual_max_alpha_order);
 };
 
 /*if(X.len[X.top-1]==2 && ao==1)
 {
  print_action_history();
  sleep(10);                   
 };*/
 /*int i = X.len[X.top-1]/2 - 1;
 if(i>=0 && i<max_correlator_order)
  G_hist[(asign[ns]>0? 0 : 1)][i] ++;*/
 stat->nstat ++;
}

void process_observable_stat(t_observable_stat* stat) //It is assumed that process_mc_stat was already called!!!
{
 /*if(stack_stat_file!=NULL)
 {
  FILE* f = fopen(stack_stat_file, "w");
  if(f==NULL)
   logs_WriteError("Cannot open the file %s for writing", stack_stat_file);
  else
   print_stack_statistics(f);
  fclose(f);
 };*/
 
 //if(logs_noise_level>=2)
 // print_stack_statistics(stdout);
 
 double source_norm          = action_create_amplitude(NULL);
 double normalization_factor = source_norm/(1.0 - mean_nA);
 logs_Write(0, "Normalization factor of multiple-trace correlators: %2.4E\n", normalization_factor);
 //normalization_factor = normalization_factor/(1.0 + mean_sign*normalization_factor);
 normalization_factor = normalization_factor/(1.0 + normalization_factor);
 
 logs_Write(0, "Normalization factor of factorized single-trace correlators: %2.4E\n", normalization_factor);

 //Saving the data for G2 
 logs_Write(0, "Collected data for G2 (up to alpha_order = %i):", stat->actual_max_alpha_order);
  
 /*FILE* fobs = fopen(observables_file, "a");
 if(fobs==NULL)
  logs_WriteError("Cannot open the file %s for writing", observables_file);*/
 
 for(int ia=0; ia<=stat->actual_max_alpha_order; ia++)
 {
  //char* outstr = NULL;
  //sprintf_append(&outstr, "Order %03i: ", ia);
  logs_Write(0, "Order %03i: ", ia);
  
  for(int m=0; m<lat_vol; m++)
  {
   double rescaling_factor     = NN*cc;
   double aG2   =     ((double)(stat->G2_hist[0][ia*lat_vol + m])  -  (double)(stat->G2_hist[1][ia*lat_vol + m]))/(double)(stat->nstat);
   int    nG2   =     (        (stat->G2_hist[0][ia*lat_vol + m])  +          (stat->G2_hist[1][ia*lat_vol + m]));
   double eG2   = sqrt((double)nG2)/(double)(stat->nstat);
   double sG2   = (nG2>0?  (double)(stat->nstat)*aG2/nG2 : 0.0);
          aG2  *= rescaling_factor*normalization_factor*pow(alpha, -(double)ia);
          eG2  *= rescaling_factor*normalization_factor*pow(alpha, -(double)ia);
   //sprintf_append(&outstr, "G[%02i] = %+2.3E +/- %2.1E (sp = %+2.2E, tnp = % 10i)  ", aG2,  eG2,  sG2,  (int)round(nG2_link));       
   logs_Write(0, "G2[%02i] = %+2.3E +/- %2.1E (sp = %+2.2E, tnp = % 10i)", m, aG2,  eG2,  sG2, nG2);
  };
  logs_Write(0, "");
  
  /*for(int mall=0; mall<stat->nG4; mall++)
  {
   double rescaling_factor     = NN*cc*cc;
   double aG4   =     ((double)(stat->G4_hist[0][ia*stat->nG4 + mall])  -  (double)(stat->G4_hist[1][ia*stat->nG4 + mall]))/(double)(stat->nstat);
   int    nG4   =     (        (stat->G4_hist[0][ia*stat->nG4 + mall])  +          (stat->G4_hist[1][ia*stat->nG4 + mall]));
   double eG4   = sqrt((double)nG4)/(double)(stat->nstat);
   double sG4   = (nG4>0?  (double)(stat->nstat)*aG4/nG4 : 0.0);
          aG4  *= rescaling_factor*normalization_factor*pow(alpha, -(double)ia);
          eG4  *= rescaling_factor*normalization_factor*pow(alpha, -(double)ia);
   //sprintf_append(&outstr, "G[%02i] = %+2.3E +/- %2.1E (sp = %+2.2E, tnp = % 10i)  ", aG2,  eG2,  sG2,  (int)round(nG2_link)); 
   if(nG4>0)
   logs_Write(0, "G4[%02i] = %+2.3E +/- %2.1E (sp = %+2.2E, tnp = % 10i)", mall, aG4,  eG4,  sG4, nG4);
  };
  logs_Write(0, "");*/
   
 }; //End of loop over ia
  
  /*double  G2_total[2] = {0.0, 0.0};
  double dG2_total[2] = {0.0, 0.0};
  double  G2_link[2]  = {0.0, 0.0};
  double dG2_link[2]  = {0.0, 0.0};
  
  for(int s=0; s<2; s++)
  {
   for(int i=0; i<lat_vol; i++)
   {
    int m[4];
    lat_idx2coords(i, m);
     G2_total[s] +=      (double)(G2_hist[s][ia*lat_vol + i]);
    dG2_total[s] +=      (double)(G2_hist[s][ia*lat_vol + i]);
    for(int mu=0; mu<DIM; mu++)
    {
     double w = cos(2.0*M_PI*(double)m[mu]/(double)lat_size[mu]);
      G2_link[s]  +=   w*(double)(G2_hist[s][ia*lat_vol + i]);
     dG2_link[s]  += w*w*(double)(G2_hist[s][ia*lat_vol + i]); 
    }; 
   };
    G2_link[s] *= 1.0/(double)DIM;
   dG2_link[s] *= 1.0/(double)DIM;
  }; 
  
  double aG2_total  =     ( G2_total[0]  -       G2_total[1])/(double)nmc;
  double nG2_total  = fabs( G2_total[0]) + fabs( G2_total[1]);     
  double eG2_total  = sqrt(dG2_total[0]  +      dG2_total[1])/(double)nmc;
  double sG2_total  = (nG2_total>0.0? (double)nmc*aG2_total/nG2_total : 0.0);
         aG2_total *= rescaling_factor*normalization_factor*pow(1.0/alpha, (double)ia)*lambda; //The last multiplication is necessary to get <g g> correlators
         eG2_total *= rescaling_factor*normalization_factor*pow(1.0/alpha, (double)ia)*lambda;
 
  sprintf_append(&outstr, "G_total = %+2.3E +/- %2.1E (sp = %+2.2E, tnp = % 10i)  ", aG2_total, eG2_total, sG2_total, (int)round(nG2_total));
  
  double aG2_link   =     ( G2_link[0]   -       G2_link[1])/(double)nmc;
  double nG2_link   = fabs( G2_link[0])  + fabs( G2_link[1]);     
  double eG2_link   = sqrt(dG2_link[0]   +      dG2_link[1])/(double)nmc;
  double sG2_link   = (nG2_link>0.0?  (double)nmc*aG2_link/nG2_link : 0.0);
         aG2_link  *= rescaling_factor*normalization_factor*pow(1.0/alpha, (double)ia)*lambda;
         eG2_link  *= rescaling_factor*normalization_factor*pow(1.0/alpha, (double)ia)*lambda;
         
  sprintf_append(&outstr, "G_link  = %+2.3E +/- %2.1E (sp = %+2.2E, tnp = % 10i)  ", aG2_link,  eG2_link,  sG2_link,  (int)round(nG2_link));
      
  logs_Write(0, "%s", outstr);
  
  fprintf(fobs, "%02i %+2.4E %+2.4E %+2.4E %+2.4E\n", ia, aG2_total, eG2_total, aG2_link, eG2_link);
 };
 
 fclose(fobs); */
 
 
 /*   G2_total_sum     += G2_total[ia];
  G2_total_sum_err += G2_total_err[ia];
  G2_link_sum      += G2_link[ia];
  G2_link_sum_err  += G2_link_err[ia]; */
  
 /*for(int ia=0; ia<=actual_max_alpha_order; ia++)
 {
  G2_total_err[ia] = sqrt(G2_total_err[ia]);
  G2_link[ia]      = G2_link[ia]/(double)DIM;
  G2_link_err[ia]  = sqrt(G2_link_err[ia]/(double)DIM);
  logs_Write(0, "ia = %i \t G2TOTAL: %+2.4E +/- %+2.4E \t G2LINK: %+2.4E +/- %+2.4E", ia, G2_total[ia], G2_total_err[ia], G2_link[ia], G2_link_err[ia]);
 };
 G2_total_sum_err = sqrt(G2_total_sum_err);
 G2_link_sum      = G2_link_sum/(double)DIM;
 G2_link_sum_err  = sqrt(G2_link_sum_err/(double)DIM);
 logs_Write(0, "Sum over all orders: G2TOTAL = %+2.4E +/- %2.4E \t G2LINK = %+2.4E +/- %2.4E\n", G2_total_sum, G2_total_sum_err, G2_link_sum, G2_link_sum_err);*/
 
 /*logs_Write(0, "\n\t Total G for different orders:");
 for(int i=0; i<max_correlator_order; i++)
 {
  rescaling_factor = NN*pow(cc, (double)(i + 1)); 
  double G  =      (double)(G_hist[0][i] - G_hist[1][i]) /(double)nmc;
  double dG = sqrt((double)(G_hist[0][i] + G_hist[1][i]))/(double)nmc;
  G  *= rescaling_factor*normalization_factor;
  dG *= rescaling_factor*normalization_factor;
  logs_Write(0, "\t G%02i:\t %2.4E +/- %2.4E", 2*i + 2, G, dG);
 };
 logs_Write(0, "");*/
 
 /*if(observables_file!=NULL)
 {
  FILE* f = fopen(observables_file, "a");
  if(f==NULL)
   logs_WriteError("Cannot open the file %s for writing", observables_file);
  else
   fprintf(f, "%2.4E %2.4E %2.4E %2.4E %2.4E\n", lambda, G2_link, G2_link_err, G2_total, G2_total_err);
  fclose(f);
 };*/
}

void free_observable_stat(t_observable_stat* stat)
{
 for(int s=0; s<2; s++)
 {
  SAFE_FREE(stat->G2_hist[s]);
  SAFE_FREE(stat->G4_hist[s]);
  SAFE_FREE(stat->G_hist[s]);
 };
}

/* G2_total[ia]     += G;
   G2_total_err[ia] += SQR(dG);
    
   for(int mu=0; mu<DIM; mu++)
   {
    double cf = cos(2.0*M_PI*(double)m[mu]/(double)lat_size[mu]);
    G2_link[ia]      += G*cf;
    G2_link_err[ia]  += SQR(dG*cf);  
   }; */
