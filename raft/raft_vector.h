#ifndef RAFT_VECTOR_H
#define RAFT_VECTOR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h> 

typedef struct {

   double * p_data;
   int size;
   int stride;

}raft_vector;

// Standard for an empty vector:
static const raft_vector RAFT_VECTOR_EMPTY = { 0, 0, 0 };

raft_vector raft_vector_create( int size );
void raft_vector_destroy( raft_vector * p_vec );
raft_vector raft_vector_scan(FILE *fp, int size);
int raft_vector_print(FILE *fp, raft_vector vec);


// Linear index for an element in the vector:
#define raft_vector_linear_index( vec, i )\
( ( i ) * ( vec ).stride )

// Reference for an element of the vector:
#define raft_vector_element( vec, i )\
( *( ( vec ).p_data + raft_vector_linear_index( vec, i ) ) )

// Reference for the number of elements of the vector:
#define raft_vector_size( vec ) (vec.size)


#ifdef __cplusplus
} //extern "C" {
#endif

#endif // #ifndef RAFT_VECTOR_H
