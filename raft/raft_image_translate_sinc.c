#include "raft_image_translate_sinc.h"
#include "raft_matrix_translate_sinc.h"

void raft_image_translate_sinc( raft_image image, double delta, int nthreads )
{
   // Horizontal displacement:
   delta *= ( image.data.columns - 1 );
   delta /= ( image.br_x - image.tl_x );

   // Phony vector:
   raft_vector deltas;
   deltas.p_data = &delta;
   deltas.size = image.data.lines;
   deltas.stride = 0;

   // Apply translations:
   raft_matrix_translate_many_sinc( deltas, image.data, nthreads );
}
