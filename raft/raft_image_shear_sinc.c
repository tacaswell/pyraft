#include "raft_image_shear_sinc.h"
#include "raft_matrix_translate_sinc.h"

void raft_image_shear_sinc( raft_image image, double alpha, int nthreads )
{
   double h_to_v_sampling_ratio = ( ( image.data.columns - 1 ) / ( image.tl_y - image.br_y ) ) *
                                  ( ( image.br_x - image.tl_x ) / ( image.data.lines - 1 ) );

   // Horizontal displacements:
   raft_vector deltas = raft_vector_create( image.data.lines );
   double center = ( deltas.size - 1 ) / 2.0;
   int i = 0;
   for ( ; i < deltas.size; ++i )
      raft_vector_element( deltas, i ) = alpha * ( center - i ) * h_to_v_sampling_ratio;

   // Apply translations:
   raft_matrix_translate_many_sinc( deltas, image.data, nthreads );

   // Release resources:
   raft_vector_destroy( &deltas );
}
