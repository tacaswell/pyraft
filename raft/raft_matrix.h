#ifndef RAFT_MATRIX_H
#define RAFT_MATRIX_H

#include "raft_vector.h" // For raft_vector
#include <stdio.h> // For FILE

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {

   double * p_data;  // Pointer to first element;
   int lines;        // Number of lines in matrix;
   int line_stride;  // Distance between successive elements in a line;
   int columns;      // Number of columns in matrix;
   int column_stride;// Distance between successive elements in a column;

} raft_matrix;

// Standard for an empty matrix:
static const raft_matrix RAFT_MATRIX_EMPTY = { 0, 0, 0, 0, 0 };
// Verify whether is an empty matrix:
#define raft_matrix_is_empty( mat ) ( !( ( mat ).p_data ) )

// Creates raft_matrix from FORTRAN-style matrix:
raft_matrix raft_matrix_column_major( double * p_data, int lines, int columns );
// Creates embedded raft_matrix from FORTRAN-style matrix:
raft_matrix raft_matrix_column_major_emb( double * p_data, int lines, int columns, int column_stride );
// Creates raft_matrix from C-style matrix:
raft_matrix raft_matrix_row_major( double * p_data, int lines, int columns );
// Creates embedded raft_matrix from C-style matrix:
raft_matrix raft_matrix_row_major_emb( double * p_data, int lines, int columns, int line_stride );
// Reads raft_matrix from file:
raft_matrix raft_matrix_read( FILE * p_file );
// Reads raft_matrix from file:  
raft_matrix raft_matrix_scan(FILE * p_file, int lines, int columns);
// Transpose matrix:
raft_matrix raft_matrix_transpose( raft_matrix mat );
// Flip matrix upside-down:
raft_matrix raft_matrix_flipud( raft_matrix mat );
// Create matrix:
raft_matrix raft_matrix_create( int lines, int columns );


// Destroy matrix:
void raft_matrix_destroy( raft_matrix * p_mat );

// Write raft_matrix to file. Returns 0 if successful:
int raft_matrix_write( raft_matrix mat, FILE * p_file );
// Print raft_matrix to file. Returns 0 if successful:
int raft_matrix_print(raft_matrix mat, FILE * p_file);


// Get a line from a matrix as a vector:
raft_vector raft_matrix_get_line( raft_matrix mat, int line );
// Get a column from a matrix as a vector:
raft_vector raft_matrix_get_column( raft_matrix mat, int column );

// Linear index for an element in the matrix:
#define raft_matrix_linear_index( mat, i, j )\
( ( ( i ) * ( mat ).column_stride ) +  ( ( j ) * ( mat ).line_stride ) )

// Is it possible to use a macro inside a macro? If so,
// we should implement this in terms of the above. TODO!
// Reference for an element of the matrix:
#define raft_matrix_element( mat, i, j )\
( *( ( mat ).p_data + ( ( i ) * ( mat ).column_stride ) +  ( ( j ) * ( mat ).line_stride ) ) )

// Reference for the number of lines of matrix:
#define raft_matrix_nlines( mat ) ( mat.lines )

// Reference for the number of columns of matrix:
#define raft_matrix_ncolumns( mat ) ( mat.columns )

#ifdef __cplusplus
} // extern "C" {
#endif

#endif // #ifndef RAFT_MATRIX_H
