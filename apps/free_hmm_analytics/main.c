#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include <clue_logs.h>
#include <rand_num_generators.h>

#include <unistd.h>

#include <sd_metropolis.h>
#include <free_hmm_genus_expansion.h>

#include "parameters.h"


int* genus_counter = NULL;

void my_pairing_processor(int* p, int n, int* c)
{
 int nloops = count_contractions(p, n, c);
 int genus  = nphi/2 + ntr - nloops;
 logs_Write(1, "%i loops, genus = %i\n", nloops, genus);
  
 if(genus%2!=0)
  logs_WriteError("Odd genus g = %i!!!", genus);
 if(genus/2>=max_genus)
  logs_WriteError("genus g = %i > max_genus = %i!!!", genus, max_genus); 
 genus_counter[genus/2] ++;
 
 if(logs_noise_level>=1)
 {
  char* pstr = print_pairing_as_pairs(p, n);
  logs_Write(0, "%s", pstr);
  SAFE_FREE(pstr);
 }; 
}


int main(int argc, char *argv[])
{
 ansi_colors = 1;
 print_errors_to_stderr = 0;
 logs_Write(0, "\n Experimenting with 1/N expansion for the free HMM model\n");
 
 parse_command_line_options(argc, argv);
 
 //Initialize contractions array for the single-trace correlator
 DECLARE_AND_MALLOC(contractions, int, nphi);
 init_contractions(contractions, ngs, ntr);
 print_contractions(contractions, nphi);
 
 max_genus = nphi/4 + ntr;
 SAFE_MALLOC(genus_counter, int, max_genus);
 for(int i=0; i<max_genus; i++)
  genus_counter[i] = 0;
   
 generate_pairings(my_pairing_processor, nphi, contractions); 
 
 logs_Write(0, "\nGeneration of pairings finished... \n");
 
 char* trstr = NULL;
 sprintf_append(&trstr, "<1/N tr(phi^%i) ", ngs[0]);
 for(int itr=1; itr<ntr; itr++)
  sprintf_append(&trstr, " 1/N tr(phi^%i)", ngs[itr]);
 sprintf_append(&trstr, "> = ");
 int sum_started = 0; int checksum = 0;
 for(int i=0; i<max_genus; i++)
  if(genus_counter[i]>0)
  {
   sprintf_append(&trstr, "%s%i", (sum_started? " + " : ""), genus_counter[i]);
   if(i>0)
    sprintf_append(&trstr, "/N^%i", 2*i);
   sum_started = 1; 
   checksum += genus_counter[i];
  };
 logs_Write(0, "%s", trstr);
 
 if(checksum!=npairs)
  logs_WriteError("Total sum of coefficients = %i does not match the number of pairings = %i", checksum, npairs);
  
 SAFE_FREE(contractions);
 SAFE_FREE(genus_counter);
 SAFE_FREE(trstr);
 
 return EXIT_SUCCESS;
}
