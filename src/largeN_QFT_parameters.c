#include "largeN_QFT_parameters.h"

//Parameters of a generic SD equations for practically any large-N QFT
int    DIM               = 1;                      //Dimension of space
int    LT                = 2;                      //Temporal size of the system
int    LS                = 2;                      //Spatial size of the system
double lambda            = 1.31;                   //tHooft coupling constant
int    max_order         = 1;
//Algorithmic parameters which control statistical sampling
double alpha               =   1.0;                    //Rescaling of observables according to the order
double cc                  =   1.0;                    //Rescaling of observables according to the number of fields in the correlator
double NN                  =   1.0;                    //Overall rescaling of observables, also genus-dependent
int    max_stack_nel       =   2;            //Maximal number of elements in the stack characterizing the system state
int    check_stack         =   0;
//Output parameters
char   data_dir[512]       =   "";           //Directory for automatically named data files
char     suffix[512]       =   "";           //Suffix for custom labelling
char      label[512]       =   "";           //A label in addition to automatically generated suffix
 
void                   print_largeN_QFT_parameters()
{
 logs_Write(0, "\tPHYSICAL PARAMETERS OF A GENERIC LARGE-N QFT");
 logs_WriteParameter(0,                                   "DIM",       "%i",      DIM);
 logs_WriteParameter(0,                                    "LT",       "%i",      LT);
 logs_WriteParameter(0,                                    "LS",       "%i",      LS);
 logs_WriteParameter(0,                                "lambda",    "%2.4E",      lambda);
 logs_WriteParameter(0,                            "Max. order",       "%i",      max_order);
 logs_Write(0, "\tALGORITHMIC PARAMETERS OF A GENERIC SIMULATION OF LARGE-N QFT");
 logs_WriteParameter(0,                                    "cc",    "%2.4E",      cc);
 logs_WriteParameter(0,                                    "NN",    "%2.4E",      NN);
 logs_WriteParameter(0,      "Max. no of elements in the stack",       "%i",      max_stack_nel);
 logs_WriteParameter(0, "Check stack consistency at every step",       "%s",      (check_stack? "YES" : "NO"));
 logs_Write(0, "\tOUTPUT PARAMETERS");
 logs_WriteParameter(0,                        "Data directory",       "%s",      data_dir);
 logs_WriteParameter(0,                                "Suffix",       "%s",      suffix);
 if(strlen(label)>0)
 logs_WriteParameter(0,                                 "Label",       "%s",      label);
}

void largeN_QFT_suffix(char* s)
{
 sprintf(s,   "d%i_t%i_s%i_l%2.4lf", DIM, LT, LS, lambda);
}

void init_lattice()
{
 int my_lat_sizes[4];
 //Initializing the square_lattice unit
 my_lat_sizes[0] = LT;
 for(int i=1; i<DIM; i++)
  my_lat_sizes[i] = LS;
 lat_init(DIM, my_lat_sizes);
}
