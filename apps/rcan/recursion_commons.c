#include "recursion_commons.h" 

void printf_momentum_str(uint P, int n, char* prefix, char* separator, char* suffix)
{
 uint mP = P;
 fprintf(stdout, "%s", prefix);
 for(int i=0; i<2*n-1; i++)
 {
  fprintf(stdout, "%1u%s", mP%LS, separator);
  mP /= LS;
 };
 fprintf(stdout, "%1u%s\n", mP%LS, suffix);
 fflush(stdout);
}

void get_momentum_str(uint P, int n, char** s, char* separator)
{
 uint mP = P;
 for(int i=0; i<2*n-1; i++)
 {
  sprintf_append(s, "%1u%s", mP%LS, separator);
  mP /= LS;
 };
 sprintf_append(s, "%1u", mP%LS);
}

uint total_momentum(uint P, int n)
{
 uint mP = P, res = 0;
 for(int i=0; i<2*n; i++)
 {
  res += mP%LS;
  mP /= LS;
 };
 return (res%LS);
}

void unpack_momenta(uint P, int n, uint* ps)
{
 uint mP = P;
 for(int i=0; i<2*n; i++)
 {
  ps[i] = mP%LS;
  mP /= LS;
 };
}

uint pack_momenta(int n, uint* ps)
{
 uint P = 0, F = 1;
 for(int i=0; i<2*n; i++)
 {
  P += ps[i]*F;
  F *= LS;
 };
 return P;
}

