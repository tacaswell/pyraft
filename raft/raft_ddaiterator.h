#ifndef RAFT_DDAITERATOR_H
#define RAFT_DDAITERATOR_H

#include "raft_image.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {

   int i; // Current line;
   int j; // Current column;
   double sigma; // Current intersection.

   // Private members:
   int slow_move; // Is this a slow direction move?
   int * p_slow_idx; // Pointer to slow index;
   int * p_fast_idx; // Pointer to fast index;
   double epsilon; // Current slow direction position in pixel;
   double epsilon_update; // Slow direction position update per fast move;
   double delta_x; // Horizontal pixel width;
   double delta_y; // Vertical pixel width;
   double delta_slow; // Pixel width in slow direction;
   int Delta_i; // Vertical index update direction;
   int Delta_j; // Horizontal index update direction;
   int Delta_slow; // Direction of slow index update;
   int Delta_fast; // Direction of fast index update;
   double sigma_tilde; // Typical intersection value;
   int lines; // Number of lines in the image;
   int columns; // Number of columns in image;
   double bbox_tl_x; // x-coordinate of top-left of bounding box;
   double bbox_tl_y; // y-coordinate of top-left of bounding box;
   double bbox_br_x; // x-coordinate of bottom-right of bounding box;
   double bbox_br_y; // y-coordinate of bottom-right of bounding box;
   double tl_x; // x-coordinate of top-left of sampling space;
   double tl_y; // y-coordinate of top-left of sampling space;
   double br_x; // x-coordinate of bottom-right of sampling space;
   double br_y; // y-coordinate of bottom-right of sampling space;
   double x0; // Value of first x-boundary to cross integration path;
   double y0; // Value of first y-boundary to cross integration path;
   int i_init; // Possible starting line;
   int j_init; // Possible starting column;
   double sin_theta; // Sine of current angle;
   double cos_theta; // Cosine of current angle.

} raft_ddaiterator;

// Creates an iterator for a given image geometry.
void raft_ddaiterator_set_image( raft_ddaiterator * p_it, raft_image image );
// Sets the ray angle for the iterator:
void raft_ddaiterator_set_theta( raft_ddaiterator * p_it, double theta );
// Sets the ray position for the iterator:
void raft_ddaiterator_set_t( raft_ddaiterator * p_it, double t );
// Update iterator:
void raft_ddaiterator_update( raft_ddaiterator * p_it );
// Has the iterator passed the end of the image?
int raft_ddaiterator_end( raft_ddaiterator const * p_it );

#ifdef __cplusplus
} // extern "C" {
#endif

#endif // #ifndef RAFT_DDAITERATOR_H
