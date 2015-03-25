#ifndef RAFT_IMAGE_H
#define RAFT_IMAGE_H

#include "raft_matrix.h" // For raft_matrix
#include <stdio.h> // For FILE

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {

   raft_matrix data; // Holds pixel data;
   double tl_x; // x-coordinate of top-left corner;
   double tl_y; // y-coordinate of top-left corner;
   double br_x; // x-coordinate of bottom-right corner;
   double br_y; // y-coordinate of bottom-right corner;

} raft_image;

// Standard for an empty image:
static const raft_image RAFT_IMAGE_EMPTY = { { 0, 0, 0, 0, 0 }, 0.0, 0.0, 0.0, 0.0 };
// Verify whether is an empty image:
#define raft_image_is_empty( im ) ( raft_matrix_is_empty( ( im ).data ) )

// Reads raft_image from file:
raft_image raft_image_read( FILE * p_file );
// Transpose image:
raft_image raft_image_transpose( raft_image im );
// Create image:
raft_image raft_image_create( int lines, int columns );

// Create image for phantom usage
raft_image raft_image_create_phantom( int lines, int columns);

// Create image for sinogram usage
raft_image raft_image_create_sinogram( int lines, int columns);

// Create image given viewport:
raft_image raft_image_create_withviewport( int lines, int columns,
                                           double tl_x, double tl_y,
                                           double br_x, double br_y);

// Destroy image:
void raft_image_destroy( raft_image * p_im );

// Write raft_image to file. Returns 0 if successful:
int raft_image_write( raft_image image, FILE * p_file );

// Set image corner
void raft_image_set_corner(raft_image *image, 
			   double tl_x, double tl_y,
			   double br_x, double br_y);

// Horizontal sampling distance:
#define raft_image_hsamplingdistance( mat ) \
( ( ( mat ).br_x - ( mat ).tl_x ) / ( ( mat ).data.columns - 1 ) )

// Vertical sampling distance:
#define raft_image_vsamplingdistance( mat ) \
( ( ( mat ).br_y - ( mat ).tl_y ) / ( ( mat ).data.lines - 1 ) )

// Reference for a sample of the image:
#define raft_image_sample( image, i, j ) \
( raft_matrix_element( ( image ).data, i, j ) )

#ifdef __cplusplus
} // extern "C" {
#endif

#endif // #ifndef RAFT_IMAGE_H
