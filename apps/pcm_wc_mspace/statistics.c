#include "statistics.h"

t_observable_stat* init_observable_stat()
{
 t_observable_stat* my_observable_stat = (t_observable_stat *)malloc(sizeof(t_observable_stat));
 
 SAFE_MALLOC(my_observable_stat->sampling_hist, int, SQR(max_order+1));
 for(int i=0; i<SQR(max_order+1); i++)
  my_observable_stat->sampling_hist[i] = 0;
 
 for(int s=0; s<2; s++)
 {
  SAFE_MALLOC(my_observable_stat->Gx[s],      double,                 (max_order+1));
  SAFE_MALLOC(my_observable_stat->Gxy[s],     double, lat_size[DIM-1]*(max_order+1));
  for(int ao=0; ao<=max_order; ao++)
  {          
   my_observable_stat->Gx[s][ao] = 0.0;
   for(int ps=0; ps<lat_size[DIM-1]; ps++)
    my_observable_stat->Gxy[s][ao*lat_size[DIM-1] + ps] = 0.0;
  }; 
 };
 
 my_observable_stat->nstat                  = 0;
 my_observable_stat->nstat_useless          = 0;
 my_observable_stat->actual_max_alpha_order = 0;
 
 return my_observable_stat;
}

void gather_observable_stat(t_observable_stat* stat)
{
 if(X.top==1)
 {
  int si = (asign[ns]>0? 0 : 1);
  int no = X.len[X.top-1]/2;
 
  int eff_order = no + alpha_order - 1;
 
  if(eff_order <= max_order)
  {
   double W =  pow(cc, (double)no)*pow(alpha, -(double)alpha_order);
   stat->Gx[si][eff_order] += W;
   int Ps[4] = {0.0, 0.0, 0.0, 0.0}; 
   int sg = -1;
   for(int i=0; i<X.len[X.top-1]-1; i++)
   {
    addto_momentum(Ps, +1, STACK_EL(X, i));
    reduce_torus(&(Ps[DIM-1]), lat_size[DIM-1]);
    int sgi = (sg>0? si : 1 - si);
    stat->Gxy[sgi][eff_order*lat_size[DIM-1] + Ps[DIM-1]] += W;
    sg *= -1;
   };
  }
  else
   stat->nstat_useless ++;
  
  if(no-1<=max_order && alpha_order<=max_order)
   stat->sampling_hist[(max_order+1)*(no-1) + alpha_order] ++;
  
  stat->actual_max_alpha_order = MAX(alpha_order, stat->actual_max_alpha_order);
 }; 
 stat->nstat ++;
}

void process_observable_stat(t_observable_stat* stat) //It is assumed that process_mc_stat was already called!!!
{
 logs_Write(0, "Collected data for G2 (up to alpha_order = %i):", stat->actual_max_alpha_order);
 
 double f = stereo_alpha;
 double  sGx[2]  = { 0.0,  0.0}; //{1.0, 0.0}
 double* sGxy[2] = {NULL, NULL};
 for(int s=0; s<2; s++)
 {
  SAFE_MALLOC(sGxy[s], double, lat_size[DIM-1]);
  for(int ps=0; ps<lat_size[DIM-1]; ps++)
   sGxy[s][ps] = 0.0;
 };  
 sGxy[0][0] = 0.0; //1.0
 
 char scalar_filename[512];
 sprintf(scalar_filename, "%s/scalars_%s.dat", data_dir, suffix);
 FILE* scalar_file = fopen(scalar_filename, "a");
 if(scalar_file==NULL)
  logs_WriteError("Could not open the file %s for writing", scalar_filename);
 
 for(int ao=0; ao<=max_order; ao++)
 {
  FILE* correlator_file = NULL;
  
  if(save_correlators)
  {
   char correlator_filename[512];
   sprintf(correlator_filename, "%s/Gxy_%s_o%i.dat", data_dir, suffix, ao);
   correlator_file = fopen(correlator_filename, "ab");
   if(correlator_file==NULL)
    logs_WriteError("Could not open the file %s for binary appending", correlator_filename); 
  };  
  
  for(int s=0; s<2; s++)
  {
   double aGx = (stat->Gx[s][ao])/(double)(stat->nstat);
   
   int sg = (f>0.0? s : 1 - s);
   
   sGx[sg] += 2.0*fabs(f)*aGx;
   sGxy[sg][0] += 4.0*fabs(f)*aGx; 
   
   for(int ps=0; ps<lat_size[DIM-1]; ps++)
   {
    double aGxy  = (stat->Gxy[s][ao*lat_size[DIM-1] + ps])/(double)(stat->nstat);
    sGxy[sg][ps] += 4.0*fabs(f)*aGxy;
   };
  };

  double sml[2] = {0.0, 0.0};
  for(int s=0; s<2; s++)
   for(int ps=0; ps<lat_size[DIM-1]; ps++)
   {
    double link_phase = cos(2.0*M_PI*(double)ps/(double)lat_size[DIM-1]);
    int sg = (link_phase>0? s : 1 - s);
    sml[sg] += sGxy[s][ps]*fabs(link_phase);
   }; 
  
  double aGx  = sGx[0] - sGx[1];
  double rGx  = sGx[0] + sGx[1];
  double spGx = (fabs(rGx)>1.0E-10? aGx/rGx : 0.0);
  
  double aml  = sml[0] - sml[1];
  double rml  = sml[0] + sml[1];
  double spml = (fabs(rml)>1.0E-10? aml/rml : 0.0);
  
  if(scalar_file!=NULL)
   fprintf(scalar_file, "%i %+2.4E %+2.4E %+2.4E %+2.4E\n", ao, aGx, spGx, aml, spml);
  logs_Write(0, "Order %03i: Gx = %+2.4E (sp=%+2.2E), link = %+2.4E (sp=%+2.2E)", ao, aGx, spGx, aml, spml);
  
  if(save_correlators && correlator_file!=NULL)
  {
   for(int ps=0; ps<lat_size[DIM-1]; ps++)
   {
    double aGxy  = sGxy[0][ps] - sGxy[1][ps];
    //double rGxy  = sGxy[0][ps] + sGxy[1][ps]; //We do not need this so far...
    //double spGxy = (fabs(rGxy)>1.0E-10? aGxy/rGxy : 0.0);
    float buf    = (float)(aGxy);
    int res = fwrite(&(buf), sizeof(float), 1, correlator_file);
    if(res!=1)
     logs_WriteError("Could not write binary data to the correlator file");
   };
   fclose(correlator_file);
  }; 
 
  f *= stereo_alpha;
 }; //End of loop over ao
 
 if(scalar_file!=NULL)
  fclose(scalar_file);
  
 logs_Write(0, "\n");
 logs_Write(0, "Frequency of visits in different sectors:");
 for(int no=1; no<=max_order+1; no++)
 {
  printf("n = %02i:\t", no);
  for(int ao=0; ao<max_order+1; ao++)
   if(stat->sampling_hist[(max_order+1)*(no-1) + ao]>0)
    printf(" %2.2lf%%", 100.0*(double)(stat->sampling_hist[(max_order+1)*(no-1) + ao])/(double)(stat->nstat));
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
 
 for(int s=0; s<2; s++)
  SAFE_FREE(sGxy[s]);
}

void free_observable_stat(t_observable_stat* stat)
{
 SAFE_FREE(stat->sampling_hist);
 for(int s=0; s<2; s++)
 {
  SAFE_FREE(stat->Gx[s]);
  SAFE_FREE(stat->Gxy[s]);
 };
}
