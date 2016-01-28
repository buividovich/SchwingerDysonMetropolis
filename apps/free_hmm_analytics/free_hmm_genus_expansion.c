#include "free_hmm_genus_expansion.h"

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

int get_num_pairings(int n)
{
 RETURN_IF_FALSE(n%2==0, -2);
 int q = (n - 1);
 for(int i=3; i<n; i+=2)
  q *= (n - i);
 return q;
}

int pairing(t_pairing_processor P, int* p, int n, int np, int plevel, int* pcount, int* data)
{
 int res = 0;
 if(plevel==n/2)
 { //All elements are already paired, append the list to the global list of pairings
  if(*pcount<np)
  {
   logs_print_pairing(2, p, n, "pcount = %i, processing pairs ", pcount[0]);
   P(p, n, data);
   pcount[0] ++;
  }
  else
  {
   logs_WriteError("pcount = %i exceeds the total number %i of pairings", pcount[0], np);
   logs_print_pairing(0, p, n, "");
   return -1;
  };
 }
 else
 {
  logs_print_pairing(3, p, n, "plevel = %i, pcount = %i ", plevel, pcount[0]);
  int i = 0;
  while(p[i]!=-1)
   i++; //i is the first unpaired element in the list
  for(int j=i+1; j<n; j++)
   if(p[j]==-1) //j and i are the non-coinciding unpaired elements
   {
    logs_Write(3, "Pair [%i %i] found", i, j);
    p[i] = j;
    p[j] = i;
    res = pairing(P, p, n, np, plevel+1, pcount, data);
    p[i] = -1;
    p[j] = -1;
   };
 };
 return res;
}

int generate_pairings(t_pairing_processor P, int n, int* data)
{
 RETURN_IF_FALSE(n%2==0, -2);
 int np = get_num_pairings(n);
 
 DECLARE_AND_MALLOC(p, int, n);
 for(int i=0; i<n; i++)
  p[i] = -1;
 
 int pcount = 0; 
 int res = pairing(P, p, n, np, 0, &pcount, data);
 if(pcount<np)
  logs_WriteError("pcount = %i < number of pairings = %i", pcount, np);
 
 SAFE_FREE(p);
 return res;
}

void init_contractions(int* c, int* gs, int nt)
{
 int istart = 0;
 for(int itr=0; itr<nt; itr++)
 {
  for(int i=istart; i<istart+gs[itr]-1; i++)
   c[i] = i+1;
  c[istart+gs[itr]-1] = istart;
  istart += gs[itr];
 };
}

void print_contractions(int* c, int n)
{
 char* s1 = NULL;
 char* s2 = NULL;
 logs_Write(0, "\nContraction of %i indices:", n);
 for(int i=0; i<n; i++)
 {
  sprintf_append(&s1, "% 4i", i);
  sprintf_append(&s2, "% 4i", c[i]);
 };
 logs_Write(0, "\t[%s]", s1);
 logs_Write(0, "\t[%s]", s2);
 SAFE_FREE(s1);
 SAFE_FREE(s2);
 logs_Write(0, "");
}


//p is the list of pairings of n elements in the Wick's theorem
//c is the list specifying the contraction of indices in the correlator
//return value is the number of closed loops in the contractions
int count_contractions(int* p, int n, int* c)
{
 int i, j, istart, ncon = 0, nloops = 0;
 DECLARE_AND_MALLOC(icon, int, n);
 for(i=0; i<n; i++)
  icon[i] = 0;
  
 char* pstr = print_pairing_as_pairs(p, n);
 logs_Write(2, "Starting contractions of the pairing %s", pstr);
 SAFE_FREE(pstr);
 
 while(ncon<n)
 {
  //Find the first non-contracted index
  i = 0;
  while(icon[i])
   i++;
  istart = i;
  logs_Write(2, "Currently %i indices are contracted, starting contractions with i_%i", ncon, istart);
  //Start the contractions
  do{
   j = p[i]; //Pairing from Wick theorem
   logs_Write(3, "Moving from i_%i to j_%i", i, j);
   i = c[j]; //Contraction coming from correlator structure
   logs_Write(3, "Moving from j_%i to i_%i", j, i);
   icon[i] = 1;
   ncon ++;
  }while(i!=istart);
  nloops ++;
 };
 SAFE_FREE(icon);
 return nloops;
}
