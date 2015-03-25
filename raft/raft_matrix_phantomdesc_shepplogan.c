#include "raft_matrix_phantomdesc_shepplogan.h"

#define SL( i, j ) \
raft_matrix_element( result, i, j )

#define PIo10 0.31415926535897932384626433832795

raft_matrix raft_matrix_phantomdesc_shepplogan( void )
{
   // Initialize matrix:
   raft_matrix result = raft_matrix_create( 10, 6 );
   // Check if allocation was successful
   if ( raft_matrix_is_empty( result ) )
      return result;

   // Ellipse data:
   SL( 0, 0 ) = 1.0;   SL( 0, 1 ) = 0.69;    SL( 0, 2 ) = 0.92;
   SL( 0, 3 ) = 0.0;   SL( 0, 4 ) = 0.0;     SL( 0, 5 ) = 0.0;

   SL( 1, 0 ) = -0.8;  SL( 1, 1 ) = 0.6624;  SL( 1, 2 ) = 0.874;
   SL( 1, 3 ) = 0.0;   SL( 1, 4 ) = -0.0184; SL( 1, 5 ) = 0.0;

   SL( 2, 0 ) = -0.2;  SL( 2, 1 ) = 0.11;    SL( 2, 2 ) = 0.31;
   SL( 2, 3 ) = 0.22;  SL( 2, 4 ) =  0.0;    SL( 2, 5 ) = -PIo10;

   SL( 3, 0 ) = -0.2;  SL( 3, 1 ) = 0.16;    SL( 3, 2 ) = 0.41;
   SL( 3, 3 ) = -0.22; SL( 3, 4 ) =  0.0;    SL( 3, 5 ) = PIo10;

   SL( 4, 0 ) = 0.1;   SL( 4, 1 ) = 0.21;    SL( 4, 2 ) = 0.25;
   SL( 4, 3 ) = 0.0;   SL( 4, 4 ) =  0.35;   SL( 4, 5 ) = 0.0;

   SL( 5, 0 ) = 0.1;   SL( 5, 1 ) = 0.046;   SL( 5, 2 ) = 0.046;
   SL( 5, 3 ) = 0.0;   SL( 5, 4 ) =  0.1;    SL( 5, 5 ) = 0.0;

   SL( 6, 0 ) = 0.1;   SL( 6, 1 ) = 0.046;   SL( 6, 2 ) = 0.046;
   SL( 6, 3 ) = 0.0;   SL( 6, 4 ) = -0.1;    SL( 6, 5 ) = 0.0;

   SL( 7, 0 ) = 0.1;   SL( 7, 1 ) = 0.046;   SL( 7, 2 ) = 0.023;
   SL( 7, 3 ) = -0.08; SL( 7, 4 ) = -0.605;  SL( 7, 5 ) = 0.0;

   SL( 8, 0 ) = 0.1;   SL( 8, 1 ) = 0.023;   SL( 8, 2 ) = 0.023;
   SL( 8, 3 ) = 0.0;   SL( 8, 4 ) = -0.606;  SL( 8, 5 ) = 0.0;

   SL( 9, 0 ) = 0.1;   SL( 9, 1 ) = 0.023;   SL( 9, 2 ) = 0.046;
   SL( 9, 3 ) = 0.06;  SL( 9, 4 ) = -0.605;  SL( 9, 5 ) = 0.0;

   return result;
}
