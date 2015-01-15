#include "sd_metropolis_parameters.h"

int    max_recursion_depth = 10000;
double p_plus              = 0.5;

void print_metropolis_parameters()
{
 logs_Write(0, "\t PARAMETERS OF THE METROPOLIS ALGORITHM");
 logs_WriteParameter(       "Max. recursion depth",     "%i", max_recursion_depth);
 logs_WriteParameter("Probability of forward move", "%2.4lf", p_plus);
}

