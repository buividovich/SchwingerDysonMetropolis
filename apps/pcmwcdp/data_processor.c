#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <clue_logs.h>
#include <fftw3.h>

#include "parameters.h"

static const char   raw_data_dir[] = "C:\\DATA\\pcm_wc_mspace_raw_data";
static const char       data_dir[] = "G:\\LAT\\sd_metropolis\\data\\pcm_wc_mspace";
static const char exact_data_dir[] = "G:\\LAT\\sd_metropolis\\data\\pcm_wc_mspace_exact";

static const char * const cluster_data_suffixes[] = {"", "_idc", "_idc1", "_itep1", "_itep2"}; //
#define NCLUSTERS (5)
#define MMAX      (12)
#define LS        (108)

//TODO: all exact + MC data should be in the .dat file
//TODO: extrapolations

int main(int argc, char *argv[])
{
 ansi_colors = 1;  print_errors_to_stderr = 0;
 logs_noise_level = 1;
 logs_Write(0, "");
 
 parse_command_line_options(argc, argv);
 init_parameters();
 print_parameters();
 
 //Reading in the exact data 
 char scalars_exact_filename[512];
 sprintf(scalars_exact_filename, "%s\\scalars_%s.exact", exact_data_dir, data_suffix);
 FILE* scalars_exact_file = fopen(scalars_exact_filename, "r");
 if(scalars_exact_file==NULL)
  logs_WriteErrorAndTerminate("The file %s not found!", scalars_exact_filename);
 
 double invorder, GxExact0, linkExact0, mpntExact0;
 fscanf(scalars_exact_file, "%lf %lf %lf %lf\n", &invorder, &GxExact0, &linkExact0, &mpntExact0);
 fclose(scalars_exact_file);
 logs_WriteParameter(0, "  Gx at order 0, exact", "%+2.4E",   GxExact0);
 logs_WriteParameter(0, "link at order 0, exact", "%+2.4E", linkExact0);
 logs_WriteParameter(0, "mpnt at order 0, exact", "%+2.4E", mpntExact0);
 logs_Write(0, "");
 
 /*********** Reading in the numerical data ****************/
 
 int num_data_points_from_cluster[NCLUSTERS], num_data_points = 0;
 
 double   Gx[MMAX],   GxCovariance[MMAX][MMAX];
 double link[MMAX], linkCovariance[MMAX][MMAX];
 double mpnt[MMAX], mpntCovariance[MMAX][MMAX];
 
 double   Gx_sign[MMAX],   Gx_sign_error[MMAX];
 double link_sign[MMAX], link_sign_error[MMAX];
 
 for(int order=0; order<MMAX; order++)
 {
  Gx[order] = 0.0; link[order] = 0.0; mpnt[order] = 0.0;
    Gx_sign[order] = 0.0;   Gx_sign_error[order] = 0.0;
  link_sign[order] = 0.0; link_sign_error[order] = 0.0;
  for(int order1=0; order1<MMAX; order1++)
  {
     GxCovariance[order][order1] = 0.0;
   linkCovariance[order][order1] = 0.0;
   mpntCovariance[order][order1] = 0.0;
  };        
 };
 
 logs_Write(0, "Loading the scalar data... ");
 for(int icluster=0; icluster<NCLUSTERS; icluster++)
 {
  char scalars_filename[512];
  sprintf(scalars_filename, "%s\\scalars_%s%s.dat", raw_data_dir, data_suffix, cluster_data_suffixes[icluster]);
  FILE* scalars_file = fopen(scalars_filename, "r");
  if(scalars_file==NULL)
   logs_WriteError("The file %s could not be opened for reading...", scalars_filename);
  if(scalars_file!=NULL)
  {
   int order_count = 0, line_count = 0, data_count = 0, error_flag = 0;
   char buf[512];
   while(fgets(buf, 512, scalars_file)!=NULL)
   {
    double aGx[MMAX], spGx[MMAX], aml[MMAX], spml[MMAX]; int ao;
    sscanf(buf, "%i %lf %lf %lf %lf\n", &ao, &(aGx[order_count]), &(spGx[order_count]), &(aml[order_count]), &(spml[order_count]));
    if(ao!=order_count)
    {
     logs_WriteError("ao=%i is not the expected order %i at line %i of file %s", ao, order_count, line_count, scalars_filename);
     error_flag = 1;
    };
    if(order_count==MMAX-1)
    {
     //Calculating the mean and the covariance matrix
     for(int order=0; order<MMAX; order++)
     {
      Gx[order] += aGx[order]; link[order] += aml[order];
      for(int order1=0; order1<MMAX; order1++)
      {
         GxCovariance[order][order1] += aGx[order]*aGx[order1];
       linkCovariance[order][order1] += aml[order]*aml[order1];
      };
        Gx_sign[order] += spGx[order];   Gx_sign_error[order] += SQR(spGx[order]);
      link_sign[order] += spml[order]; link_sign_error[order] += SQR(spml[order]); 
     }; 
     //Updating counters
     data_count ++;
     num_data_points ++;
     order_count = 0;
    }
    else
     order_count ++; 
    line_count++;
   };
   fclose(scalars_file);
   num_data_points_from_cluster[icluster] = data_count;
   if(!error_flag)
    logs_Write(1, "Done for %.70s\t\t [%7d lines, %7d data points]", scalars_filename, line_count, data_count);
  };
 };
 logs_Write(0, "");
 
 //Preparing the FFTW3 stuff
 fftw_complex fftw_in[LS], fftw_out[LS];
 fftw_plan p = fftw_plan_dft_1d(LS, fftw_in, fftw_out, FFTW_FORWARD, FFTW_ESTIMATE);
 
 double correlatorMean[MMAX][LS], correlatorError[MMAX][LS];
 for(int order=0; order<MMAX; order++)
  for(int x=0; x<LS; x++)
  {
   correlatorMean[order][x]  = 0.0;
   correlatorError[order][x] = 0.0;
  }; 
 
 logs_Write(0, "Loading the correlators... ");
 for(int icluster=0; icluster<NCLUSTERS; icluster++)
 {
  //Open files for all orders
  FILE* correlator_files[MMAX];
  char correlator_filename_mask[512];
  sprintf(correlator_filename_mask, "%s\\Gxy_%s%s_o", raw_data_dir, data_suffix, cluster_data_suffixes[icluster]);
  int error_flag = 0;
  for(int order=0; order<MMAX; order++)
  {
   char correlator_filename[512];
   sprintf(correlator_filename, "%s%i.dat", correlator_filename_mask, order);
   correlator_files[order] = fopen(correlator_filename, "rb");
   if(correlator_files[order]==NULL)
   {
    logs_WriteError("The file %s could not be opened for reading...", correlator_filename);
    error_flag = 1;
   }; 
  };
  //Now reading all orders in parallel  
  if(!error_flag)
  {
   int eof_not_reached = 1, data_count = 0;
   float float_buf[MMAX][LS];
   double correlator_xspace[MMAX][LS];
   double ampnt[MMAX];
   while(eof_not_reached)
   {
    for(int order=0; order<MMAX; order++)
     eof_not_reached = eof_not_reached && (fread(float_buf[order], sizeof(float), LS, correlator_files[order])==LS);
    if(eof_not_reached)
    {
     //FFT of data for all orders
     for(int order=0; order<MMAX; order++)
     {
      for(int x=0; x<LS; x++)
       fftw_in[x] = (double)(float_buf[order][x]) + I*0.0;
      fftw_execute(p);
      for(int x=0; x<LS; x++)
       correlator_xspace[order][x] = creal(fftw_out[x]);
      ampnt[order] = correlator_xspace[order][LS/2];
     };
     //Correlator mean and covariance
     for(int order=0; order<MMAX; order++)
     {
      mpnt[order] += ampnt[order];
      for(int order1=0; order1<MMAX; order1++)
       mpntCovariance[order][order1] += ampnt[order]*ampnt[order1];
      for(int x=0; x<LS; x++)
      {
       correlatorMean[order][x]  += correlator_xspace[order][x];
       correlatorError[order][x] += SQR(correlator_xspace[order][x]);
      }; 
     };
     //Increasing data counter
     data_count ++;
    }; 
   };
   for(int order=0; order<MMAX; order++)
    fclose(correlator_files[order]);
   if(data_count!=num_data_points_from_cluster[icluster])
    logs_WriteError("%i data points in files %s*.dat, expected %i", data_count, correlator_filename_mask, num_data_points_from_cluster[icluster]);
   else 
    logs_Write(1, "Done for %.70s*.dat\t\t [%7d data points]", correlator_filename_mask, data_count);
  };
 };
 logs_Write(0, "");
 
 fftw_destroy_plan(p);
 
 //Processing Gx, link and midpoint data
 char   GxMeanFileName[512],   GxCovarianceFileName[512];
 char linkMeanFileName[512], linkCovarianceFileName[512];
 char mpntMeanFileName[512], mpntCovarianceFileName[512];
 sprintf(  GxMeanFileName,       "%s\\Gx_%s.mean",   data_dir, data_suffix);
 sprintf(linkMeanFileName,       "%s\\link_%s.mean", data_dir, data_suffix);
 sprintf(mpntMeanFileName,       "%s\\mpnt_%s.mean", data_dir, data_suffix);
 sprintf(  GxCovarianceFileName, "%s\\Gx_%s.cov",    data_dir, data_suffix);
 sprintf(linkCovarianceFileName, "%s\\link_%s.cov",  data_dir, data_suffix);
 sprintf(mpntCovarianceFileName, "%s\\mpnt_%s.cov",  data_dir, data_suffix);
 FILE*   GxMeanFile              = fopen(  GxMeanFileName       , "w"); 
 FILE* linkMeanFile              = fopen(linkMeanFileName       , "w");
 FILE* mpntMeanFile              = fopen(mpntMeanFileName       , "w");
 FILE*   GxCovarianceFile        = fopen(  GxCovarianceFileName , "w");
 FILE* linkCovarianceFile        = fopen(linkCovarianceFileName , "w");
 FILE* mpntCovarianceFile        = fopen(mpntCovarianceFileName , "w");
 ASSERT(GxMeanFile==NULL || linkMeanFile==NULL || mpntMeanFile==NULL || GxCovarianceFile==NULL || linkCovarianceFile==NULL || mpntCovarianceFile==NULL);
 
 char SignsFileName[512];
 sprintf(SignsFileName, "%s\\signs_%s.mean", data_dir, data_suffix);
 FILE* SignsFile = fopen(SignsFileName, "w");
 ASSERT(SignsFile==NULL);
 
 double GxNFactor, linkNFactor, mpntNFactor, NFactor;
 
 for(int order=0; order<MMAX; order++)
 {
    Gx[order] /= (double)num_data_points; 
  link[order] /= (double)num_data_points;
  mpnt[order] /= (double)num_data_points;
 
  if(order==0) //Get NFactors
  {
     GxNFactor=(  GxExact0 - 1.0)/  Gx[0];
   linkNFactor=(linkExact0 - 1.0)/link[0];
   mpntNFactor=(mpntExact0 - 1.0)/mpnt[0];
   
   NFactor = 0.5*(GxNFactor + linkNFactor);
   
   logs_Write(0, "");
   logs_Write(0, "NFactors from %i data points in total", num_data_points);
   logs_WriteParameter(0, "NFactor from   Gx", "%+2.4E",   GxNFactor);
   logs_WriteParameter(0, "NFactor from link", "%+2.4E", linkNFactor);
   logs_WriteParameter(0, "NFactor from mpnt", "%+2.4E", mpntNFactor);
   logs_WriteParameter(0, "NFactor to use   ", "%+2.4E",     NFactor);
  };
  
  for(int order1=0; order1<MMAX; order1++)
  {
     GxCovariance[order][order1] = (  GxCovariance[order][order1]/(double)num_data_points -   Gx[order]*  Gx[order1])/(double)(num_data_points-1);
   linkCovariance[order][order1] = (linkCovariance[order][order1]/(double)num_data_points - link[order]*link[order1])/(double)(num_data_points-1);
   mpntCovariance[order][order1] = (mpntCovariance[order][order1]/(double)num_data_points - mpnt[order]*mpnt[order1])/(double)(num_data_points-1);
   
     GxCovariance[order][order1] *= SQR(NFactor);
   linkCovariance[order][order1] *= SQR(NFactor);
   mpntCovariance[order][order1] *= SQR(NFactor);
   
   fprintf(  GxCovarianceFile, "%+2.6E ",   GxCovariance[order][order1]);
   fprintf(linkCovarianceFile, "%+2.6E ", linkCovariance[order][order1]);
   fprintf(mpntCovarianceFile, "%+2.6E ", mpntCovariance[order][order1]);
  };
  fprintf(  GxCovarianceFile, "\n");
  fprintf(linkCovarianceFile, "\n");
  fprintf(mpntCovarianceFile, "\n");
  
    Gx[order] = 1.0 + NFactor*  Gx[order];
  link[order] = 1.0 + NFactor*link[order];
  mpnt[order] = 1.0 + NFactor*mpnt[order];
  
  fprintf(  GxMeanFile, "%2.4E %+2.4E %2.4E\n", 1.0/(double)(order+1),   Gx[order], sqrt(fabs(  GxCovariance[order][order])) );
  fprintf(linkMeanFile, "%2.4E %+2.4E %2.4E\n", 1.0/(double)(order+1), link[order], sqrt(fabs(linkCovariance[order][order])) );
  fprintf(mpntMeanFile, "%2.4E %+2.4E %2.4E\n", 1.0/(double)(order+1), mpnt[order], sqrt(fabs(mpntCovariance[order][order])) );
  
  //Processing and saving average signs
    Gx_sign[order]       /= (double)num_data_points;
  link_sign[order]       /= (double)num_data_points;
    Gx_sign_error[order] /= (double)num_data_points;
  link_sign_error[order] /= (double)num_data_points;
  
    Gx_sign_error[order] = sqrt(fabs(  Gx_sign_error[order] - SQR(  Gx_sign[order]))/(double)(num_data_points - 1));
  link_sign_error[order] = sqrt(fabs(link_sign_error[order] - SQR(link_sign[order]))/(double)(num_data_points - 1));
  fprintf(SignsFile, "%02i %+2.4E %+2.4E %+2.4E %+2.4E\n", order, Gx_sign[order], Gx_sign_error[order], link_sign[order], link_sign_error[order]);
 };
 
 fclose(   GxMeanFile      ); 
 fclose( linkMeanFile      );
 fclose( mpntMeanFile      );
 fclose(   GxCovarianceFile);
 fclose( linkCovarianceFile);
 fclose( mpntCovarianceFile);
 fclose(SignsFile);
 
 //Processing and saving the correlators data
 for(int x=0; x<LS; x++)
  for(int order=0; order<MMAX; order++)
  {
   correlatorMean[order][x]  /= (double)num_data_points;
   correlatorError[order][x]  = sqrt(fabs((correlatorError[order][x]/(double)num_data_points - SQR(correlatorMean[order][x]))/(double)(num_data_points-1)));
   //Rescaling with NFactor
   correlatorMean[order][x]   = 1.0 + NFactor*correlatorMean[order][x];
   correlatorError[order][x] *= NFactor;
  };
 
 char CorrelatorFileName[512];
 sprintf(CorrelatorFileName, "%s\\Gxy_%s.dat", data_dir, data_suffix);
 FILE* CorrelatorFile = fopen(CorrelatorFileName, "w");
 ASSERT(CorrelatorFile==NULL);
 for(int x=0; x<=LS; x++)
 {
  fprintf(CorrelatorFile, "%03i ", x);
  for(int order=0; order<MMAX; order++)
   fprintf(CorrelatorFile, "%+2.4E %+2.4E ", correlatorMean[order][x%LS], correlatorError[order][x%LS]);
  fprintf(CorrelatorFile, "\n");       
 };
 fclose(CorrelatorFile);
 
 if(!save_summary)
  return EXIT_SUCCESS; 
 
 double acceptance_rate     = 0.0, mean_recursion_depth     = 0.0, mean_stack_top     = 0.0, one_minus_mean_nA     = 0.0, maxnA     = 0.0, mean_sign     = 0.0, mean_return_time     = 0.0, max_ns     = 0.0, max_stop     = 0.0;
 double acceptance_rate_err = 0.0, mean_recursion_depth_err = 0.0, mean_stack_top_err = 0.0, one_minus_mean_nA_err = 0.0, maxnA_err = 0.0, mean_sign_err = 0.0, mean_return_time_err = 0.0, max_ns_err = 0.0, max_stop_err = 0.0;
 
 //Now reading in the mcstat files and averaging MC characteristics
 logs_Write(0, "");
 logs_Write(0, "Loading the MC stat files...");
 for(int icluster=0; icluster<NCLUSTERS; icluster++)
 {
  char MCStatFilename[512];
  sprintf(MCStatFilename, "%s\\metropolis_stat_%s%s.dat", raw_data_dir, data_suffix, cluster_data_suffixes[icluster]);
  FILE* MCStatFile = fopen(MCStatFilename, "r");
  if(MCStatFile==NULL)
   logs_WriteError("The file %s could not be opened for reading...", MCStatFilename);
  if(MCStatFile!=NULL)
  {
   int data_count = 0; char buf[512];
   while(fgets(buf, 512, MCStatFile)!=NULL)
   {
    double my_acceptance_rate, my_mean_recursion_depth, my_mean_stack_top, my_one_minus_mean_nA, my_maxnA, my_mean_sign, my_mean_return_time;
    int my_max_ns, my_max_stop;
    sscanf(buf, "%lf %lf %3d %lf %3d %lf %lf %lf %lf\n", 
     &my_acceptance_rate, 
     &my_mean_recursion_depth, 
     &my_max_ns, 
     &my_mean_stack_top, 
     &my_max_stop, 
     &my_one_minus_mean_nA, 
     &my_maxnA, 
     &my_mean_sign,
     &my_mean_return_time
    );
    
    acceptance_rate      += my_acceptance_rate;           acceptance_rate_err      += SQR(my_acceptance_rate        );
    mean_recursion_depth += my_mean_recursion_depth;      mean_recursion_depth_err += SQR(my_mean_recursion_depth   );
    max_ns               += my_max_ns;                    max_ns_err               += SQR(my_max_ns                 );
    mean_stack_top       += my_mean_stack_top;            mean_stack_top_err       += SQR(my_mean_stack_top         );
    max_stop             += my_max_stop;                  max_stop_err             += SQR(my_max_stop               );
    one_minus_mean_nA    += my_one_minus_mean_nA;         one_minus_mean_nA_err    += SQR(my_one_minus_mean_nA      );
    maxnA                += my_maxnA;                     maxnA_err                += SQR(my_maxnA                  );
    mean_sign            += my_mean_sign;                 mean_sign_err            += SQR(my_mean_sign              );
    mean_return_time     += my_mean_return_time;          mean_return_time_err     += SQR(my_mean_return_time       );
    
    data_count ++;
   };
   if(data_count!=num_data_points_from_cluster[icluster])
    logs_WriteError("%i data points in file %s, expected %i", data_count, MCStatFilename, num_data_points_from_cluster[icluster]);
   else 
    logs_Write(1, "Done for %.70s\t\t [%7d data points]", MCStatFilename, data_count);
   fclose(MCStatFile); 
  };
 };
 
 //Processing MC stat
 #define STATPROC(val_)                                                                \
 val_           /= (double)num_data_points;                                            \
 (val_ ## _err) /= (double)num_data_points;                                            \
 (val_ ## _err)  = sqrt(fabs((val_ ## _err) - SQR(val_))/(double)(num_data_points-1)); \
 fprintf(MCStatSummaryFile, "%+2.4E %+2.4E ", val_, (val_ ## _err));
 
 char MCStatSummaryFilename[512];
 sprintf(MCStatSummaryFilename, "%s\\mc_stat_summary_%s.dat", data_dir, scan_suffix);
 FILE* MCStatSummaryFile = fopen(MCStatSummaryFilename, "a");
 if(MCStatSummaryFile!=NULL)
 {
  fprintf(MCStatSummaryFile, "%s ", scan_label); //Col 1
  STATPROC(acceptance_rate       );     //Col 2,3
  STATPROC(mean_recursion_depth  );     //Col 4,5
  STATPROC(max_ns                );     //Col 6,7
  STATPROC(mean_stack_top        );     //Col 8,9
  STATPROC(max_stop              );     //Col 10,11
  STATPROC(one_minus_mean_nA     );     //Col 12,13
  STATPROC(maxnA                 );     //Col 14,15
  STATPROC(mean_sign             );     //Col 16,17
  STATPROC(mean_return_time      );     //Col 18,19
  fprintf(MCStatSummaryFile, "\n");
  fclose(MCStatSummaryFile);
 }
 else
  logs_WriteError("Could not open the file %s for appending...", MCStatSummaryFilename);
 
 //Summaries of observables
 char GxSummaryFilename[512], linkSummaryFilename[512], mpntSummaryFilename[512];
 sprintf(  GxSummaryFilename, "%s\\Gx_summary_%s.dat",   data_dir, scan_suffix);
 sprintf(linkSummaryFilename, "%s\\link_summary_%s.dat", data_dir, scan_suffix);
 sprintf(mpntSummaryFilename, "%s\\mpnt_summary_%s.dat", data_dir, scan_suffix);
 
 FILE*   GxSummaryFile = fopen(  GxSummaryFilename, "a");
 FILE* linkSummaryFile = fopen(linkSummaryFilename, "a");
 FILE* mpntSummaryFile = fopen(mpntSummaryFilename, "a");
 
 if(GxSummaryFile!=NULL && linkSummaryFile!=NULL && mpntSummaryFile!=NULL)
 {
  fprintf(  GxSummaryFile, "%s ", scan_label);
  fprintf(linkSummaryFile, "%s ", scan_label);
  fprintf(mpntSummaryFile, "%s ", scan_label);
  for(int order=0; order<MMAX; order++)
  {
   fprintf(  GxSummaryFile, "%+2.4E %2.4E ",   Gx[order], sqrt(fabs(  GxCovariance[order][order])) );
   fprintf(linkSummaryFile, "%+2.4E %2.4E ", link[order], sqrt(fabs(linkCovariance[order][order])) );
   fprintf(mpntSummaryFile, "%+2.4E %2.4E ", mpnt[order], sqrt(fabs(mpntCovariance[order][order])) );       
  };
  fprintf(  GxSummaryFile, "\n");
  fprintf(linkSummaryFile, "\n");
  fprintf(mpntSummaryFile, "\n");
  fclose(  GxSummaryFile);
  fclose(linkSummaryFile);
  fclose(mpntSummaryFile);
 }
 else
  logs_WriteError("One of the files %s, %s or %s could not be opened...", GxSummaryFilename, linkSummaryFilename, mpntSummaryFilename);
 
 logs_Write(0, "");
 return EXIT_SUCCESS;   
}

/*  if(data_count<4)
    {
    printf(">>>%+2.4E %+2.4E %i %+2.4E %i %+2.4E %+2.4E %+2.4E %+2.4E\n",
     my_acceptance_rate, 
     my_mean_recursion_depth, 
     my_max_ns, 
     my_mean_stack_top, 
     my_max_stop, 
     my_one_minus_mean_nA, 
     my_maxnA, 
     my_mean_sign,
     my_mean_return_time);
    };*/
