#ifndef _HMM_ACTIONS_H_
#define _HMM_ACTIONS_H_

//Create new factorized-out line
double action_create_amplitude(void* data_in);
int    action_create_do(void** data_out, void* data_in);
int    action_create_undo(void* data_out, void* data_in);

//Create new factorized-in line
double action_evolve_line_amplitude(void* data_in);
int    action_evolve_line_do(void** data_out, void* data_in);
int    action_evolve_line_undo(void* data_out, void* data_in);

//Create new vertex
double action_evolve_vertex_amplitude(void* data_in);
int    action_evolve_vertex_do(void** data_out, void* data_in);
int    action_evolve_vertex_undo(void* data_out, void* data_in);

//Join two sets of lines
double action_join_amplitude(void* data_in);
int    action_join_do(void** data_out, void* data_in);
int    action_join_undo(void* data_out, void* data_in);


#endif
