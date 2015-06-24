#include "detailed_statistics.h"

int my_genus = 1;

int h15 = 0;
int h24 = 0;
int h33 = 0;
int h42 = 0;
int h51 = 0;

void gather_detailed_statistics()
{
 if(X.top==2 && genus==my_genus && X.nel==6)
 {
  if(X.len[0]==1 && X.len[1]==5)
   h15 ++;
  if(X.len[0]==2 && X.len[1]==4)
   h24 ++;
  if(X.len[0]==3 && X.len[1]==3)
   h33 ++;
  if(X.len[0]==4 && X.len[1]==2)
   h42 ++;
  if(X.len[0]==5 && X.len[1]==1)
   h51 ++;   
 };
}

void print_detailed_statistics()
{
 genus = 0;
 double source_norm          = action_create_amplitude(NULL);
 double normalization_factor = source_norm/(1.0 - mean_nA);
 logs_Write(0, "Normalization factor of multiple-trace correlators: %2.4E\n", normalization_factor);
 double rescaling_factor     = f_genus[my_genus]*SQR(NN_genus[my_genus])*pow(cc_genus[my_genus], 6);
 double G15 = (double)h15/(double)nmc*normalization_factor*rescaling_factor;
 double G24 = (double)h24/(double)nmc*normalization_factor*rescaling_factor;
 double G33 = (double)h33/(double)nmc*normalization_factor*rescaling_factor;
 double G42 = (double)h42/(double)nmc*normalization_factor*rescaling_factor;
 double G51 = (double)h51/(double)nmc*normalization_factor*rescaling_factor;
 double Gtotal = G15 + G24 + G33 + G42 + G51;
 logs_Write(0, "Detailed statistics: ");
 logs_WriteParameter(0, "G15",     "%02.2lf, should be 10", G15);
 logs_WriteParameter(0, "G24",     "%02.2lf, should be 9 ", G24);
 logs_WriteParameter(0, "G33",     "%02.2lf, should be 12", G33);
 logs_WriteParameter(0, "G42",     "%02.2lf, should be 9 ", G42);
 logs_WriteParameter(0, "G51",     "%02.2lf, should be 10", G51);
 logs_WriteParameter(0, "\tTotal", "%02.2lf", Gtotal);
 
 logs_Write(0, "");
}
