#include "detailed_statistics.h"

int my_genus = 1;

int h11     = 0;
int h1111   = 0;
int h111111 = 0;

int h15 = 0;
int h24 = 0;
int h33 = 0;
int h42 = 0;
int h51 = 0;

int h121 = 0;
int h211 = 0;
int h112 = 0;

int h51f2 = 0;
int h51f4 = 0;
int h51f5 = 0;

void check_h15_action_history()
{
 int res = 0, c0 = 0, c1 = 0, c2 = 0;
 if(ns==2)
 {
  c0 = (action_history[0].action_id==6 && action_history[0].action_data_in==0); //State_initializer
  c1 = (action_history[1].action_id==2                                       ); //Evolve line
  c2 = (action_history[2].action_id==5 && action_history[0].action_data_in==0); //Split
  res = c0 && c1 && c2;
 };
 
 if(!res) 
 {
  logs_WriteError("15 state reached with %i steps, res = %i, c0 = %i, c1 = %i, c2 = %i", ns, res, c0, c1, c2);
  print_action_history();
  system("PAUSE");    
 };
}

void gather_detailed_statistics()
{
 if(X.top==2 && genus==my_genus && X.nel==6)
 {
  if(X.len[1]==1 && X.len[0]==5)
  {
   h15 ++;
   check_h15_action_history();
  }; 
  if(X.len[1]==2 && X.len[0]==4)
   h24 ++;
  if(X.len[1]==3 && X.len[0]==3)
   h33 ++;
  if(X.len[1]==4 && X.len[0]==2)
   h42 ++;
  if(X.len[1]==5 && X.len[0]==1)
  {
   h51 ++;
   /*logs_WriteWarning("51");
   print_action_history();
   system("PAUSE");*/
   if(action_history[ns].action_id==2)
    h51f2++;
   if(action_history[ns].action_id==4)
    h51f4++;
   if(action_history[ns].action_id==5)
    h51f5++;
  }; 
 };
 
 if(X.top<=6 && X.nel==X.top)
 {
  if(X.nel==2 && genus==1)
   h11 ++;
  if(X.nel==4 && genus==2)
   h1111 ++;
  if(X.nel==6 && genus==3)
   h111111 ++; 
 };
 
 if(X.top==3 && genus==my_genus && X.nel==4)
 {
  if(X.len[2]==1 && X.len[1]==2 && X.len[0]==1)
   h121 ++;
  if(X.len[2]==2 && X.len[1]==1 && X.len[0]==1)
   h211 ++;
  if(X.len[2]==1 && X.len[1]==1 && X.len[0]==2)
   h112 ++; 
 };
}

void print_detailed_statistics()
{
 genus = 0;
 double source_norm          = my_source_norm();
 double normalization_factor = source_norm/(1.0 - mean_nA);
 logs_Write(0, "Normalization factor of multiple-trace correlators: %2.4E\n", normalization_factor);
 
 double rescaling_factor     = f_genus[1]*pow(NN_genus[1], 2)*pow(cc_genus[1], 2);
 double G11     = (double)h11/(double)nmc*normalization_factor*rescaling_factor;
 
        rescaling_factor     = f_genus[2]*pow(NN_genus[2], 4)*pow(cc_genus[2], 4);
 double G1111   = (double)h1111/(double)nmc*normalization_factor*rescaling_factor;
 
        rescaling_factor     = f_genus[3]*pow(NN_genus[3], 6)*pow(cc_genus[3], 6);
 double G111111 = (double)h111111/(double)nmc*normalization_factor*rescaling_factor;
 
 
        rescaling_factor     = f_genus[my_genus]*SQR(NN_genus[my_genus])*pow(cc_genus[my_genus], 6);
 double G15 = (double)h15/(double)nmc*normalization_factor*rescaling_factor;
 double G24 = (double)h24/(double)nmc*normalization_factor*rescaling_factor;
 double G33 = (double)h33/(double)nmc*normalization_factor*rescaling_factor;
 double G42 = (double)h42/(double)nmc*normalization_factor*rescaling_factor;
 double G51 = (double)h51/(double)nmc*normalization_factor*rescaling_factor;
 double Gtotal = G15 + G24 + G33 + G42 + G51;
 
 double G51f2 = (double)h51f2/(double)nmc*normalization_factor*rescaling_factor;
 double G51f4 = (double)h51f4/(double)nmc*normalization_factor*rescaling_factor;
 double G51f5 = (double)h51f5/(double)nmc*normalization_factor*rescaling_factor;
 
        rescaling_factor     = f_genus[my_genus]*pow(NN_genus[my_genus], 3)*pow(cc_genus[my_genus], 4);
 double G121  = (double)h121/(double)nmc*normalization_factor*rescaling_factor;
 double G211  = (double)h211/(double)nmc*normalization_factor*rescaling_factor;
 double G112  = (double)h112/(double)nmc*normalization_factor*rescaling_factor;

 logs_Write(0, "Detailed statistics on tr(phi) correlators: ");
 logs_WriteParameter(0, "G11",       "%2.2lf, should be 1  (% 8i occurences)", G11,     h11);
 logs_WriteParameter(0, "G1111",     "%2.2lf, should be 3  (% 8i occurences)", G1111,   h1111);
 logs_WriteParameter(0, "G111111",   "%2.2lf, should be 15 (% 8i occurences)", G111111, h111111);
 
 logs_Write(0, "Detailed statistics: ");
 logs_WriteParameter(0, "G15",     "%2.2lf, should be 10 (% 8i occurences)", G15, h15);
 logs_WriteParameter(0, "G24",     "%2.2lf, should be 9  (% 8i occurences)", G24, h24);
 logs_WriteParameter(0, "G33",     "%2.2lf, should be 12 (% 8i occurences)", G33, h33);
 logs_WriteParameter(0, "G42",     "%2.2lf, should be 9  (% 8i occurences)", G42, h42);
 logs_WriteParameter(0, "G51",     "%2.2lf, should be 10 (% 8i occurences)", G51, h51);
 logs_WriteParameter(0, "\tTotal", "%2.2lf", Gtotal);
 
 logs_Write(0, "Detailed statistics for G51: ");
 logs_WriteParameter(0, "G51f2",     "%2.2lf, should be 6 (% 8i occurences)", G51f2, h51f2);
 logs_WriteParameter(0, "G51f4",     "%2.2lf, should be 2 (% 8i occurences)", G51f4, h51f4);
 logs_WriteParameter(0, "G51f5",     "%2.2lf, should be 2 (% 8i occurences)", G51f5, h51f5);
 logs_WriteParameter(0, "\tTotal", "%2.2lf (%i occurences)", G51f2 + G51f4 + G51f5, h51f2 + h51f4 + h51f5);
 
 logs_Write(0, "Detailed statistics for G121 and G211: ");
 logs_WriteParameter(0, "G121",     "%2.2lf, should be 1 (% 8i occurences)", G121, h121);
 logs_WriteParameter(0, "G211",     "%2.2lf, should be 1 (% 8i occurences)", G211, h211);
 logs_WriteParameter(0, "G112",     "%2.2lf, should be 1 (% 8i occurences)", G112, h112);
 
 logs_Write(0, "");
}
