#include "raft_ddaiterator.h"
#include <math.h>

#define SIGNAL( x )\
( ( ( x ) > 0.0 ) ? 1.0 : ( ( ( x ) < 0.0 ) ? -1.0 : 0.0 ) )
#define MAX( x, y )\
( ( ( x ) > ( y ) ) ? ( x ) : ( y ) )
#define MIN( x, y )\
( ( ( x ) < ( y ) ) ? ( x ) : ( y ) )
#define INDEX( x, x_i, delta_x, Delta_j )\
( ( Delta_j >= 0 ) ? \
  floor( ( ( x - x_i ) / delta_x ) - 0.5 ) + 1 : \
  ceil( ( ( x - x_i ) / delta_x ) - 0.5 ) \
)

void raft_ddaiterator_set_image( raft_ddaiterator * p_it, raft_image image )
{
   // Pixel widths:
   p_it->delta_x = raft_image_hsamplingdistance( image );
   p_it->delta_y = raft_image_vsamplingdistance( image );

   // Bounding-box limits:
   p_it->bbox_tl_x = image.tl_x - ( p_it->delta_x / 2.0 );
   p_it->bbox_br_x = image.br_x + ( p_it->delta_x / 2.0 );
   p_it->bbox_tl_y = image.tl_y - ( p_it->delta_y / 2.0 );
   p_it->bbox_br_y = image.br_y + ( p_it->delta_y / 2.0 );
   // Sampling limits:
   p_it->tl_x = image.tl_x;
   p_it->br_x = image.br_x;
   p_it->tl_y = image.tl_y;
   p_it->br_y = image.br_y;

   // Image discretization:
   p_it->lines = image.data.lines;
   p_it->columns  = image.data.columns;
}

inline void raft_ddaiterator_setslowfast( raft_ddaiterator * p_it,
                                          int * p_fast_idx,
                                          int * p_slow_idx,
                                          double delta_fast,
                                          double delta_slow,
                                          double small_trig,
                                          double large_trig,
                                          int Delta_fast,
                                          int Delta_slow
                                        )
{
   // Pointers to slow and fast index:
   p_it->p_fast_idx = p_fast_idx;
   p_it->p_slow_idx = p_slow_idx;

   //Dimensions:
   p_it->epsilon_update = fabs( delta_fast * small_trig / large_trig );
   p_it->delta_slow = fabs( delta_slow );

   //Pixel update directions:
   p_it->Delta_slow = Delta_slow;
   p_it->Delta_fast = Delta_fast;

   // Typical intersection value:
   p_it->sigma_tilde = fabs( delta_fast / large_trig );
}

void raft_ddaiterator_set_theta( raft_ddaiterator * p_it, double theta )
{
   // Trigonometric values:
   p_it->cos_theta = cos( theta ); p_it->sin_theta = sin( theta );
   // Auxiliary variable:
   int sign;

   // Possible limiting boundaries in the beginning of
   // intersection with bounding box:
   sign = SIGNAL( p_it->sin_theta );
   p_it->x0 = sign * MAX( sign * p_it->bbox_br_x, sign * p_it->bbox_tl_x );
   p_it->Delta_j = -SIGNAL( p_it->delta_x ) * sign;

   sign = SIGNAL( p_it->cos_theta );
   p_it->y0 = sign * MIN( sign * p_it->bbox_br_y, sign * p_it->bbox_tl_y );
   p_it->Delta_i = SIGNAL( p_it->delta_y ) * sign;

   // Possible starting pixel values:
   p_it->i_init = ( p_it->Delta_i < 0 ) ? ( p_it->lines - 1 ) : 0;
   p_it->j_init = ( p_it->Delta_j < 0 ) ? ( p_it->columns - 1 ) : 0;

   if ( fabs( p_it->delta_x * p_it->sin_theta ) <=
        fabs( p_it->delta_y * p_it->cos_theta )
      ) // Vertical is fast direction:
   {
      raft_ddaiterator_setslowfast( p_it,
                                    &( p_it->i ), &( p_it->j ),
                                    p_it->delta_y, p_it->delta_x,
                                    p_it->sin_theta, p_it->cos_theta,
                                    p_it->Delta_i, p_it->Delta_j
                                  );
   }
   else // Horizontal is slow direction:
   {
      raft_ddaiterator_setslowfast( p_it,
                                    &( p_it->j ), &( p_it->i ),
                                    p_it->delta_x, p_it->delta_y,
                                    p_it->cos_theta, p_it->sin_theta,
                                    p_it->Delta_j, p_it->Delta_i
                                  );
   }
}

void raft_ddaiterator_set_t( raft_ddaiterator * p_it, double t )
{
   // Always start with fast direction:
   p_it->slow_move = 0;

   // Parametric points in intersections with bounding box:
   double sx, sy;
   // Position in boundary:
   double xy;

   if ( p_it->sin_theta == 0.0 )
   {
      p_it->j = INDEX( t, p_it->tl_x, p_it->delta_x, p_it->Delta_j );
      p_it->i = p_it->i_init;
      p_it->epsilon = 0.0;
      p_it->epsilon_update = 0.0;
   } else if ( p_it->cos_theta == 0.0 )
   {
      p_it->i = INDEX( t, p_it->tl_y, p_it->delta_y, p_it->Delta_i );
      p_it->j = p_it->j_init;
      p_it->epsilon = 0.0;
      p_it->epsilon_update = 0.0;
   } else if ( ( sx = ( ( t * p_it->cos_theta - p_it->x0 ) / p_it->sin_theta ) ) >
               ( sy = ( ( p_it->y0 - t * p_it->sin_theta ) / p_it->cos_theta ) )
             ) // Hits bounding box at x == x0, need to determinal i:
   {
      // y-value at boundary:
      xy = ( sx * p_it->cos_theta ) + ( t * p_it->sin_theta );
      // Corresponding index:
      p_it->i = INDEX( xy, p_it->tl_y, p_it->delta_y, p_it->Delta_i );
      p_it->j = p_it->j_init;

      // Displacement from previous pixel:
      p_it->epsilon = fabs( xy - ( p_it->tl_y +
                                   ( p_it->i - ( 0.5 * p_it->Delta_i ) ) * p_it->delta_y
                                 )
                          );

      // Final corrections:
      if ( p_it->p_slow_idx == &( p_it->i ) )
      { // Computed displacement was from previous pixel:
         p_it->j -= p_it->Delta_j;
         // But we need to be at the correct pixel:
         raft_ddaiterator_update( p_it );
      }
      else
      { // Required displacement is in the other axis!
         // Need difference in the other direction:
         p_it->epsilon -= fabs( p_it->delta_y );
         // First intersection value:
         p_it->sigma = fabs( p_it->epsilon /= p_it->cos_theta );
         // Displacement required in the other axis:
         p_it->epsilon = fabs( p_it->sigma * p_it->sin_theta );
      }
   }
   else // Hits bounding box at y == y0, need to determine i:
   {
      // x-value at boundary:
      xy = ( -sy * p_it->sin_theta ) + ( t * p_it->cos_theta );
      // Corresponding index:
      p_it->j = INDEX( xy, p_it->tl_x, p_it->delta_x, p_it->Delta_j );
      p_it->i = p_it->i_init;

      // Displacement from previous pixel:
      p_it->epsilon = fabs( xy - ( p_it->tl_x +
                                   ( p_it->j - ( 0.5 * p_it->Delta_j ) ) * p_it->delta_x
                                 )
                          );

      // Final corrections:
      if ( p_it->p_slow_idx == &( p_it->j ) )
      { // Computed displacement was from previous pixel:
         p_it->i -= p_it->Delta_i;
         // But we need to be at the correct pixel:
         raft_ddaiterator_update( p_it );
      }
      else
      { // Required displacement is in the other axis!
         // Need difference in the other direction:
         p_it->epsilon -= fabs( p_it->delta_x );
         // First intersection value:
         p_it->sigma = fabs( p_it->epsilon /= p_it->sin_theta );
         // Displacement required in the other axis:
         p_it->epsilon = fabs( p_it->sigma * p_it->cos_theta );
      }
   }
}

void raft_ddaiterator_update( raft_ddaiterator * p_it )
{
   // Auxiliary variables:
   double epsilon_previous, alpha;

   if ( p_it->slow_move )
   {
      // Update slow index:
      *( p_it->p_slow_idx ) += p_it->Delta_slow;

      // Update intersection value:
      p_it->sigma = p_it->sigma_tilde - p_it->sigma;

      // Next is a fast move:
      p_it->slow_move = 0;
   }
   else
   {
      // Update fast index:
      *( p_it->p_fast_idx ) += p_it->Delta_fast;

      // Records previous epsilon value:
      epsilon_previous = p_it->epsilon;
      // Update epsilon:
      p_it->epsilon += p_it->epsilon_update;

      // Will next be a slow move?
      p_it->slow_move = ( p_it->epsilon > p_it->delta_slow );
      // If so, use fancy intersection computation.
      if ( p_it->slow_move )
      {
         // Update epsilon value:
         p_it->epsilon -= p_it->delta_slow;
         if ( p_it->epsilon )
         {
            // $alpha \leftarrow (\Delta_{slow} - \epsilon_p) / \epsilon$
            alpha = ( p_it->delta_slow - epsilon_previous ) / p_it->epsilon;
            // $\alpha \leftarrow \alpha/(1 + \alpha) = (\Delta_{slow} - \epsilon_p)/\Delta$
            alpha /= ( 1.0 + alpha );
         }
         else
            alpha = 1.0;

         // Update intersection value:
         p_it->sigma = alpha * p_it->sigma_tilde;
      }
      else
      {
         // Update intersection value:
         p_it->sigma = p_it->sigma_tilde;
      }
   }
}

int raft_ddaiterator_end( raft_ddaiterator const * p_it )
{
   return( ( p_it->i < 0 ) ||
           ( p_it->i >= p_it->lines ) ||
           ( p_it->j < 0 ) ||
           ( p_it->j >= p_it->columns )
         );
}

