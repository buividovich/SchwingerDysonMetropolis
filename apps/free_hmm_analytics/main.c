#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <clue_logs.h>
#include <rand_num_generators.h>

#include <sd_metropolis.h>
#include <free_hmm_genus_expansion.h>

#include "parameters.h"

int main(int argc, char *argv[])
{
 ansi_colors = 1;
 print_errors_to_stderr = 0;
 logs_Write(0, "\n Experimenting with 1/N expansion for the free HMM model\n");
 
 parse_command_line_options(argc, argv);
 
 int npairings = get_num_pairings(nphi);
 logs_Write(0, "There are %i ways to group %i objects in pairs", npairings, nphi);
 DECLARE_AND_MALLOC(pairings, int*, npairings);
 for(int i=0; i<npairings; i++)
  SAFE_MALLOC(pairings[i], int, nphi);
  
 generate_pairings(nphi, pairings); 
 
 logs_Write(1, "\nGeneration of pairings finished... \n");
 
 if(logs_noise_level>=1)
 {
  for(int i=0; i<npairings; i++)
  {
   char* ps1 =  print_pairing_as_list(pairings[i], nphi);
   char* ps2 = print_pairing_as_pairs(pairings[i], nphi);
   logs_Write(1, "\t %04i \t %s => %s", i, ps1, ps2);
   SAFE_FREE(ps1);
   SAFE_FREE(ps2);
  };
  logs_Write(0, "");
 };
 
  
 for(int i=0; i<npairings; i++)
  SAFE_FREE(pairings[i]);
 SAFE_FREE(pairings);
 
 return EXIT_SUCCESS;
}
