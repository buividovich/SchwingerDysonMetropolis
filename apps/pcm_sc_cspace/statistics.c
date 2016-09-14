#include "statistics.h"

t_observable_stat* init_observable_stat()
{
 t_observable_stat* my_observable_stat = (t_observable_stat *)malloc(sizeof(t_observable_stat));

 SAFE_MALLOC(my_observable_stat->sampling_hist, int, SQR(max_order+1));
 for(int i=0; i<SQR(max_order+1); i++)
  my_observable_stat->sampling_hist[i] = 0;

 for(int s=0; s<2; s++)
 {
  SAFE_MALLOC(my_observable_stat->G2[s],     double,    lat_vol*(max_order+1));
  for(int ao=0; ao<=max_order; ao++)
  {          
   for(int x=0; x<lat_vol; x++)
    my_observable_stat->G2[s][ao*lat_vol + x] = 0.0;
  }; 
 };

 my_observable_stat->nstat                  = 0;
 my_observable_stat->nstat_useless          = 0;
 my_observable_stat->actual_max_beta_order  = 0;
 my_observable_stat->max_history_nel        = 0;
 my_observable_stat->max_stack_nel          = 0;
 
 return my_observable_stat;
}

void gather_observable_stat(t_observable_stat* stat)
{
 if(X.top==1 && X.nel==2)
 {
  int si = (asign[ns]>0? 0 : 1);
  if(beta_order <= max_order)
  {
   double W =  pow(alpha, -(double)beta_order);
   stat->G2[si][beta_order*lat_vol + lat_distance(STACK_EL(X, 0), STACK_EL(X, 1))] += W;
  }
  else
   stat->nstat_useless ++;
 };
 
 int no = X.nel/2; 
 if(no-1<=max_order && beta_order<=max_order)
  stat->sampling_hist[(max_order+1)*(no-1) + beta_order] ++;
  
 stat->actual_max_beta_order = MAX(beta_order, stat->actual_max_beta_order); 
 stat->max_stack_nel         = MAX(X.nel,      stat->max_stack_nel);
 
 stat->nstat ++;
}

void process_observable_stat(t_observable_stat* stat) //It is assumed that process_mc_stat was already called!!!
{
 double source_norm          = action_create_amplitude(NULL);
 logs_Write(0, "Source norm:\t %+2.6E", source_norm);
 logs_Write(0, "Max.history length:\t %i", stat->max_history_nel);
 logs_Write(0, "Max.stack   length:\t %i", stat->max_stack_nel);
 logs_Write(0, "Total %i calls to gather(), %s%i useless (%2.2lf%%)", stat->nstat, ANSI_COLOR(red), stat->nstat_useless, 100.0*(double)stat->nstat_useless/(double)stat->nstat);
 
 char mean_nA_filename[512];
 sprintf(mean_nA_filename, "%s/mean_nA_%s.dat", data_dir, suffix);
 FILE* mean_nA_file = fopen(mean_nA_filename, "a");
 if(mean_nA_file==NULL)
  logs_WriteError("Could not open the file %s for writing", mean_nA_filename);
 fprintf(mean_nA_file, "%+2.8E\n", 1.0 - mean_nA);
 fclose(mean_nA_file);
  
 logs_Write(0, "Collected data for G2 (up to beta_order = %i):", stat->actual_max_beta_order);
   
 char scalar_filename[512];
 sprintf(scalar_filename, "%s/scalars_%s.dat", data_dir, suffix);
 FILE* scalar_file = fopen(scalar_filename, "a");
 if(scalar_file==NULL)
  logs_WriteError("Could not open the file %s for writing", scalar_filename);
  
 //Fixing the overall scale by the order 0 result for the unitarity 
 double order_zero_unitarity_num   = 0.0;
 order_zero_unitarity_num += (stat->G2[0][0] - stat->G2[1][0])/(double)(stat->nstat);
 
 double normalization_factor = 1.0/order_zero_unitarity_num; 

 for(int bo=0; bo<=max_order; bo++)
 {
  double sml[2] = {0.0, 0.0};
  double sun[2] = {0.0, 0.0};
  for(int s=0; s<2; s++)
  {
   sun[s] = normalization_factor*stat->G2[s][bo*lat_vol + 0]/(double)(stat->nstat);
   sml[s] = 0.0;
   for(int mu=0; mu<DIM; mu++)
   {
    sml[s] += normalization_factor*stat->G2[s][bo*lat_vol + lat_shift_fwd(0, mu)]/(double)(stat->nstat);
    sml[s] += normalization_factor*stat->G2[s][bo*lat_vol + lat_shift_bwd(0, mu)]/(double)(stat->nstat);
   };
   sml[s] /= (double)(2*DIM);
  }; 
   
  double aml  = (sml[0] - sml[1]);
  double rml  = (sml[0] + sml[1]);
  double spml = (fabs(rml)>1.0E-10? aml/rml : 1.0);
  
  double aun  = (sun[0] - sun[1]);
  double run  = (sun[0] + sun[1]);
  double spun = (fabs(run)>1.0E-10? aun/run : 1.0);
  
  if(scalar_file!=NULL)
   fprintf(scalar_file, "%i %+2.4E %+2.4E %+2.4E %+2.4E\n", bo, aun, spun, aml, spml);
  logs_Write(0, "Order %03i: link = %+2.4E (sp=%+2.2E), unitarity = %+2.4E (sp = %2.2E)", bo, aml, spml, aun, spun);
  
 }; //End of loop over ao
 
 if(scalar_file!=NULL)
  fclose(scalar_file);
  
 logs_Write(0, "\n");
 logs_Write(0, "Frequency of visits in different sectors:");
 for(int no=1; no<=max_order+1; no++)
 {
  printf("n = %02i:\t", no);
  for(int bo=0; bo<=max_order; bo++)
   if(stat->sampling_hist[(max_order+1)*(no-1) + bo]>0)
    printf(" %2.2lf%%", 100.0*(double)(stat->sampling_hist[(max_order+1)*(no-1) + bo])/(double)(stat->nstat));
  printf("\n");
 };
 printf("\n");
 
 logs_Write(0, "\n");
 logs_Write(0, "Numbers of visits in different sectors:");
 for(int no=1; no<=max_order+1; no++)
 {
  printf("n = %02i:\t", no);
  for(int ao=0; ao<=max_order; ao++)
  {
   if(no+ao-1<=max_order)
    printf(" % 8i", stat->sampling_hist[(max_order+1)*(no-1) + ao]);
  };  
  printf("\n");
 };
 printf("\n");
 
 if(save_sampling_hist)
 {
  char sampling_hist_filename[512];
  sprintf(sampling_hist_filename, "%s/sampling_hist_%s.dat", data_dir, suffix);
  FILE* sampling_hist_file = fopen(sampling_hist_filename, "w");
  if(sampling_hist_file!=NULL)
  {
   for(int no=1; no<=max_order+1; no++)
   {
    for(int ao=0; ao<max_order+1; ao++)
     fprintf(sampling_hist_file, "%2.4E ", (double)(stat->sampling_hist[(max_order+1)*(no-1) + ao])/(double)(stat->nstat));
    fprintf(sampling_hist_file, "\n");
   };
   fclose(sampling_hist_file);
  }; 
 };
 
 //for(int s=0; s<2; s++)
 // SAFE_FREE(sG2[s]);
}

void free_observable_stat(t_observable_stat* stat)
{
 SAFE_FREE(stat->sampling_hist);
 for(int s=0; s<2; s++)
  SAFE_FREE(stat->G2[s]);
}

/* FILE* correlator_file = NULL;
  
if(save_correlators)
{
 char correlator_filename[512];
 sprintf(correlator_filename, "%s/Gxy_%s_o%i.dat", data_dir, suffix, bo);
 correlator_file = fopen(correlator_filename, "a");
 if(correlator_file==NULL)
  logs_WriteError("Could not open the file %s for writing", correlator_filename);
};  */

/*  if(save_correlators && correlator_file!=NULL)
  {
   for(int pt=0; pt<LT; pt++)
   {
    double aG2  = normalization_factor*(sG2[0][pt] - sG2[1][pt]);
    double rG2  = normalization_factor*(sG2[0][pt] + sG2[1][pt]);
    double spG2 = (fabs(rG2)>1.0E-10? aG2/rG2 : 0.0);
    fprintf(correlator_file, "%i %+2.4E %+2.4E\n", pt, aG2, spG2);
   };
   fclose(correlator_file);
  }; */
