#include <stdio.h>
#include <stdlib.h>

#include <clue_logs.h>

#include <sd_metropolis.h>

#include "hmm_parameters.h"
#include "hmm_actions.h"

int main(int argc, char *argv[])
{
 ansi_colors = 1;
 logs_Write(0, "DIAGRAMMATIC MONTE-CARLO FOR HERMITIAN PHI4 MATRIX MODEL IN THE PLANAR LIMIT");
 
 parse_command_line_options(argc, argv);
 init_parameters();
 print_parameters();
 
 init_hmm_actions();
 
 
 free_hmm_actions();
 
 return 0;
}
