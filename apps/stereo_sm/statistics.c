#include "statistics.h"

t_observable_stat* init_observable_stat()
{
 t_observable_stat* my_observable_stat = (t_observable_stat *)malloc(sizeof(t_observable_stat));
 
 SAFE_MALLOC(my_observable_stat->sampling_hist, int, SQR(max_alpha_order+1));
 for(int i=0; i<SQR(max_alpha_order+1); i++)
  my_observable_stat->sampling_hist[i] = 0;
 
 for(int s=0; s<2; s++)
 {
  SAFE_MALLOC(my_observable_stat->Gx[s],  double,    (max_alpha_order+1));
  SAFE_MALLOC(my_observable_stat->Gxy[s], double, LT*(max_alpha_order+1));
  for(int ao=0; ao<=max_alpha_order; ao++)
  {          
   my_observable_stat->Gx[s][ao] = 0.0;
   for(int t=0; t<LT; t++)
    my_observable_stat->Gxy[s][ao*LT + t] = 0.0;
  }; 
 };
 
 my_observable_stat->Gtotal                 = 0.0;
 
 my_observable_stat->nstat                  = 0;
 my_observable_stat->nstat_useless          = 0;
 my_observable_stat->actual_max_alpha_order = 0;
 
 return my_observable_stat;
}

void gather_observable_stat(t_observable_stat* stat)
{
 if(X.top==1)
 {
  int ao = alpha_order;
  int no = X.len[X.top-1]/2;
  int si = (asign[ns]>0? 0 : 1);
 
  int eff_order = no + ao;
 
  if(eff_order <= (max_alpha_order+1))
  {
   double W =  pow(cc, (double)no)*pow(alpha, -(double)ao);
   stat->Gx[si][eff_order-1] += W;
   stat->Gtotal += W;
   int Pt[4] = {0.0, 0.0, 0.0, 0.0}; 
   int sg = -1;
   for(int i=0; i<X.len[X.top-1]-1; i++)
   {
    addto_momentum(Pt, +1, STACK_EL(X, i));
    reduce_torus(&(Pt[0]), lat_size[0]);
    int sgi = (sg>0? si : 1 - si);
    stat->Gxy[sgi][(eff_order-1)*LT + Pt[0]] += W;
    sg *= -1;
   };
  }
  else
   stat->nstat_useless ++;
  
  if(no-1<=max_alpha_order && ao<=max_alpha_order)
   stat->sampling_hist[(max_alpha_order+1)*(no-1) + ao] ++;
  
  stat->actual_max_alpha_order = MAX(ao, stat->actual_max_alpha_order);
 }; 
 stat->nstat ++;
}

void process_observable_stat(t_observable_stat* stat) //It is assumed that process_mc_stat was already called!!!
{
 double source_norm          = action_create_amplitude(NULL);
 logs_Write(0, "Source norm:\t %+2.6E", source_norm);
 double normalization_factor = source_norm/(1.0 - mean_nA);
 //logs_Write(0, "Normalization factor of multiple-trace correlators: %2.4E\n", normalization_factor);
 //normalization_factor = normalization_factor/(1.0 + normalization_factor);
 
 logs_Write(0, "Normalization factor of factorized single-trace correlators: %2.4E\n", normalization_factor);
 logs_Write(0, "Total %i calls to gather(), %s%i useless (%2.2lf%%)", stat->nstat, ANSI_COLOR(red), stat->nstat_useless, 100.0*(double)stat->nstat_useless/(double)stat->nstat);
 
 char mean_nA_filename[512];
 sprintf(mean_nA_filename, "%s/mean_nA_%s.dat", data_dir, suffix);
 FILE* mean_nA_file = fopen(mean_nA_filename, "a");
 if(mean_nA_file==NULL)
  logs_WriteError("Could not open the file %s for writing", mean_nA_filename);
 fprintf(mean_nA_file, "%+2.8E\n", 1.0 - mean_nA);
 fclose(mean_nA_file);
  
 logs_Write(0, "Collected data for G2 (up to alpha_order = %i):", stat->actual_max_alpha_order);
 
 double f = stereo_alpha;
 double  sGx[2]  = { 0.0,  0.0}; //{1.0, 0.0}
 double* sGxy[2] = {NULL, NULL};
 for(int s=0; s<2; s++)
 {
  SAFE_MALLOC(sGxy[s], double, LT);
  for(int pt=0; pt<LT; pt++)
   sGxy[s][pt] = 0.0;
 };  
 sGxy[0][0] = 0.0; //1.0
 
 char scalar_filename[512];
 sprintf(scalar_filename, "%s/scalars_%s.dat", data_dir, suffix);
 FILE* scalar_file = fopen(scalar_filename, "a");
 if(scalar_file==NULL)
  logs_WriteError("Could not open the file %s for writing", scalar_filename);
 
 for(int ao=0; ao<=max_alpha_order; ao++)
 {
  char correlator_filename[512];
  sprintf(correlator_filename, "%s/Gxy_%s_o%i.dat", data_dir, suffix, ao);
  FILE* correlator_file = fopen(correlator_filename, "a");
  if(correlator_file==NULL)
   logs_WriteError("Could not open the file %s for writing", correlator_filename);
  
  for(int s=0; s<2; s++)
  {
   double aGx = (stat->Gx[s][ao])/(double)(stat->nstat);
   aGx *= NN; //*normalization_factor;
   
   int sg = (f>0.0? s : 1 - s);
   
   sGx[sg] += 2.0*fabs(f)*aGx;
   sGxy[sg][0] += 4.0*fabs(f)*aGx; 
   
   for(int pt=0; pt<LT; pt++)
   {
    double aGxy  = (stat->Gxy[s][ao*LT + pt])/(double)(stat->nstat);
    aGxy *= NN; //*normalization_factor;
    sGxy[sg][pt] += 4.0*fabs(f)*aGxy;
   };
  };

  double sml[2] = {0.0, 0.0};
  for(int s=0; s<2; s++)
   for(int pt=0; pt<LT; pt++)
   {
    double link_phase = cos(2.0*M_PI*(double)pt/(double)LT);
    int sg = (link_phase>0? s : 1 - s);
    sml[sg] += sGxy[s][pt]*fabs(link_phase);
   }; 
  
  double aGx  = sGx[0] - sGx[1];
  double rGx  = sGx[0] + sGx[1];
  double spGx = (fabs(rGx)>1.0E-10? aGx/rGx : 0.0);
  
  double aml  = sml[0] - sml[1];
  double rml  = sml[0] + sml[1];
  double spml = (fabs(rml)>1.0E-10? aml/rml : 0.0);
  
  if(scalar_file!=NULL)
   fprintf(scalar_file, "%i %+2.4E %+2.4E %+2.4E %+2.4E\n", ao, aGx, spGx, aml, spml);
  logs_Write(0, "Order %03i: Gx = %+2.4lf (sp=%+2.2lf), link = %+2.4lf (sp=%+2.2lf)", ao, aGx, spGx, aml, spml);
  
  if(correlator_file!=NULL)
  {
   for(int pt=0; pt<LT; pt++)
   {
    double aGxy  = sGxy[0][pt] - sGxy[1][pt];
    double rGxy  = sGxy[0][pt] + sGxy[1][pt];
    double spGxy = (fabs(rGxy)>1.0E-10? aGxy/rGxy : 0.0);
    fprintf(correlator_file, "%i %+2.4E %+2.4E\n", pt, aGxy, spGxy);
   };
   fclose(correlator_file);
  }; 
 
  f *= stereo_alpha;
 }; //End of loop over ao
 
 if(scalar_file!=NULL)
  fclose(scalar_file);
  
 logs_Write(0, "\n");
 logs_Write(0, "Frequency of visits in different sectors:");
 for(int no=1; no<=max_alpha_order+1; no++)
 {
  printf("n = %i:\t", no);
  for(int ao=0; ao<max_alpha_order+1; ao++)
   if(stat->sampling_hist[(max_alpha_order+1)*(no-1) + ao]>0)
    printf(" % 8i", stat->sampling_hist[(max_alpha_order+1)*(no-1) + ao]);
    //printf(" % 2.2lf%%", 100.0*(double)(stat->sampling_hist[(max_alpha_order+1)*(no-1) + ao])/(double)(stat->nstat));
  printf("\n");
 };
 
 stat->Gtotal *= NN*normalization_factor/(double)(stat->nstat);
 
 logs_Write(0, "Total G:\t%2.4E", stat->Gtotal);
 
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
