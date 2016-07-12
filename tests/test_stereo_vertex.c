#include "test_stereo_vertex.h"

void test_stereo_vertex(uint64_t ntrials)
{
 int* Q[3];
 for(int i=0; i<3; i++)
  SAFE_MALLOC(Q[i], int, 1);
 logs_Write(0, "");
 logs_Write(0, "Lowest-order vertex:"); 
 for((*(Q[0]))=0; (*(Q[0]))<LS; (*(Q[0]))++)
  for((*(Q[1]))=0; (*(Q[1]))<LS; (*(Q[1]))++)
   for((*(Q[2]))=0; (*(Q[2]))<LS; (*(Q[2]))++)
   {
    int Pt[2];
    double res = vertex(Q, Pt, 3);
    logs_Write(0, "%i%i%i\t %+2.4lf [%i -> %+2.4lf]", (*(Q[0])), (*(Q[1])), (*(Q[2])), res, Pt[0], lat_propagator(Pt, P.mass_sq));
   };
 logs_Write(0, "");  
  
}

