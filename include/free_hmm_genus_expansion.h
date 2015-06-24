#ifndef _FREE_HMM_GENUS_EXPANSION_H_
#define _FREE_HMM_GENUS_EXPANSION_H_

#include <stdarg.h>

#include <clue_errors.h>
#include <clue_logs.h>
#include <clue_io.h>

typedef void (*t_pairing_processor)(int* p, int n, int* data);

int get_num_pairings(int n);

//Generate all possible pairings of n objects and feed them to the user-defined function P
int generate_pairings(t_pairing_processor P, int n, int* data);

void init_contractions(int* c, int* gs, int nt);
void print_contractions(int* c, int n);

//p is the list of pairings of n elements in the Wick's theorem
//c is the list specifying the contraction of indices in the correlator
//return value is the number of closed loops in the contractions
int count_contractions(int* pairing, int n, int* c);

void  logs_print_pairing(int msg_level, int* p, int n, char *fmt, ...);
char* print_pairing_as_pairs(int* p, int n);
char* print_pairing_as_list(int* p, int n);

#endif
