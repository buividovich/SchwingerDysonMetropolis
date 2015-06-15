#include "free_hmm_genus_expansion.h"

int*  pairs   = NULL; //plist has size of 2n
int   plevel  = 0;
int   nlist   = 0;
int   pcount  = 0;

char* print_pairing_as_pairs(int* p, int n)
{
 char* pstr = NULL;
 DECLARE_AND_MALLOC(pair_printed, int, n);
 for(int i=0; i<n; i++)
  pair_printed[i] = 0;
 for(int i=0; i<n; i++)
 {
  int j = p[i];
  if(!pair_printed[i] && !pair_printed[j])
  {
   sprintf_append(&pstr, "[%i %i] ", i, j);
   pair_printed[i] = 1;
   pair_printed[j] = 1;
  };
 };
 SAFE_FREE(pair_printed);
 return pstr;
}

char* print_pairing_as_list(int* p, int n)
{
 char* pstr = NULL;
 sprintf_append(&pstr, "[%i", p[0]);
 for(int i=1; i<n; i++)
  sprintf_append(&pstr, " %i", p[i]);
 sprintf_append(&pstr, "]");
 return pstr;
}

void logs_print_pairing(int msg_level, int* p, int n, char *fmt, ...)
{
 if(logs_noise_level>=msg_level)
 {
  if(ansi_colors && msg_level>=0 && msg_level<4)
   fprintf(stdout, "%s", ansi_color_codes[msg_level]);
  va_list va;
  va_start(va, fmt);
  vfprintf(stdout, fmt, va);
  va_end(va);
  
  fprintf(stdout, "[%i", p[0]);
  for(int i=1; i<n; i++)
   fprintf(stdout, " %i", p[i]);
  fprintf(stdout, "]"); 
  
  if(ansi_colors)
   fprintf(stdout, "\x1b[0m");
  fprintf(stdout, "\n");
  fflush(stdout);
 };    
}

int pairing(int** plist, int n)
{
 int res = 0;
 if(plevel==n/2)
 { //All elements are already paired, append the list to the global list of pairings
  if(pcount<nlist)
  {
   SAFE_MALLOC(plist[pcount], int, n);
   memcpy(plist[pcount], pairs, n*sizeof(int));
   pcount ++;
   logs_print_pairing(2, pairs, n, "pcount = %i, adding pairs ", pcount);   
  }
  else
  {
   logs_WriteError("pcount = %i >= nlist = %i", pcount, nlist);
   logs_print_pairing(0, pairs, n, "");
   return -1;  
  };  
 }
 else
 {
  logs_print_pairing(3, pairs, n, "plevel = %i, pcount = %i", plevel, pcount);
  int i = 0;
  while(pairs[i]!=-1)
   i++; //i is the first unpaired element in the list
  for(int j=i+1; j<n; j++)
   if(pairs[j]==-1) //j and i are the non-coinciding unpaired elements
   {
    logs_Write(3, "Pair [%i %i] found", i, j);
    pairs[i] = j;
    pairs[j] = i;
    plevel ++;
    res = pairing(plist, n);
    plevel --;
    pairs[i] = -1;
    pairs[j] = -1;
   };
 };
 return res;
}

int get_num_pairings(int n)
{
 RETURN_IF_FALSE(n%2==0, -2);
 int q = (n - 1);
 for(int i=3; i<n; i+=2)
  q *= (n - i);
 return q;
}

int generate_pairings(int n, int** plist)
{
 RETURN_IF_FALSE(n%2==0, -2);
 nlist = get_num_pairings(n);
 SAFE_MALLOC(pairs, int, n);
 
 for(int i=0; i<n; i++)
  pairs[i] = -1;
 
 plevel = 0;
 pcount = 0; 
 int res = pairing(plist, n);
 
 SAFE_FREE(pairs);
 return res;
}
