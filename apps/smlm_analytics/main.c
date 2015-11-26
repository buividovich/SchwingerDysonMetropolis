#include <stdio.h>
#include <stdlib.h>

#include <largeN_QFT_parameters.h>

#include "parameters.h"
#include "recursion.h"

//TODO Plan of implementation:
// 1. Fast implementation for any LS using uint-packing, without buffering
// 2. Implement buffering
// 3. Implement fast operations for LS==2 using #defines/ inline functions

int main(int argc, char *argv[])
{
 ansi_colors = 1;
 print_errors_to_stdout = 1;
 parse_command_line_options(argc, argv);
 init_parameters();
 print_parameters();
 init_recursion();
 
 free_recursion();
 return EXIT_SUCCESS;
}
