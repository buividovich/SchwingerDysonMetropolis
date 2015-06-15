#ifndef _FREE_HMM_GENUS_EXPANSION_H_
#define _FREE_HMM_GENUS_EXPANSION_H_

#include <stdarg.h>

#include <clue_errors.h>
#include <clue_logs.h>
#include <clue_io.h>

int get_num_pairings(int n);
int generate_pairings(int n, int** plist);

void logs_print_pairing(int msg_level, int* p, int n, char *fmt, ...);
char* print_pairing_as_pairs(int* p, int n);
char* print_pairing_as_list(int* p, int n);

#endif
