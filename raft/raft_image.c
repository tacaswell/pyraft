#include "raft_image.h"

#define PI 3.1415926535897932384626433832795

raft_image raft_image_read( FILE * p_file )
{
   raft_image result;
   result.data = raft_matrix_read( p_file );

   if ( raft_image_is_empty( result ) )
   {
      result.tl_x = result.tl_y = result.br_x = result.br_y = 0.0;
      return result;
   }

   if( !fread( &( result.tl_x ), sizeof( double ), 1, p_file ) )
   {
      raft_image_destroy( &result );
      return( result );
   }
   if( !fread( &( result.tl_y ), sizeof( double ), 1, p_file ) )
   {
      raft_image_destroy( &result );
      return( result );
   }
   if( !fread( &( result.br_x ), sizeof( double ), 1, p_file ) )
   {
      raft_image_destroy( &result );
      return( result );
   }
   if( !fread( &( result.br_y ), sizeof( double ), 1, p_file ) )
   {
      raft_image_destroy( &result );
      return( result );
   }

   return result;
}

raft_image raft_image_transpose( raft_image im )
{
   raft_image result;

   result.data = raft_matrix_transpose( im.data );
   result.tl_x = im.tl_y;
   result.tl_y = im.tl_x;
   result.br_x = im.br_y;
   result.br_y = im.br_x;

   return result;
}

raft_image raft_image_create( int lines, int columns )
{
   raft_image result;
   result.data = raft_matrix_create( lines, columns );

   if ( raft_image_is_empty( result ) )
   {
      result.tl_x = result.tl_y = result.br_x = result.br_y = 0.0;
      return result;
   }

   result.tl_x = result.br_y = -1.0;
   result.br_x = result.tl_y =  1.0;

   return result;
}

raft_image raft_image_create_phantom( int lines, int columns )
{
   raft_image result;
   result.data = raft_matrix_create( lines, columns );

   if ( raft_image_is_empty( result ) )
   {
      result.tl_x = result.tl_y = result.br_x = result.br_y = 0.0;
      return result;
   }

   result.tl_x = result.br_y = -1.0;
   result.br_x = result.tl_y =  1.0;

   return result;
}

raft_image raft_image_create_sinogram( int lines, int columns )
{
   raft_image result;
   result.data = raft_matrix_create( lines, columns );

   if ( raft_image_is_empty( result ) )
   {
      result.tl_x = result.tl_y = result.br_x = result.br_y = 0.0;
      return result;
   }

   result.tl_x = 0.0;     
   result.br_y = 1.0;
   
   result.br_x = PI * ( 1.0 - 1.0 / columns );
   result.tl_y = -1.0;

   return result;
}


raft_image raft_image_create_withviewport( int lines, int columns,
                                           double tl_x, double tl_y,
                                           double br_x, double br_y)
{
   raft_image result;
   result.data = raft_matrix_create( lines, columns );

   if ( raft_image_is_empty( result ) )
   {
      result.tl_x = result.tl_y = result.br_x = result.br_y = 0.0;
      return result;
   }

   result.tl_x = tl_x;
   result.tl_y = tl_y;
   result.br_x = br_x;
   result.br_y = br_y;

   return result;
}

void raft_image_destroy( raft_image * p_im )
{
   raft_matrix_destroy( &( p_im->data ) );
   p_im->tl_x = p_im->tl_y = p_im->br_x = p_im->br_y = 0.0;
}

int raft_image_write( raft_image image, FILE * p_file )
{
   int count = raft_matrix_write( image.data, p_file );
   if ( count )
      return count;

   if( !fwrite( &( image.tl_x ), sizeof( double ), 1, p_file ) )
      return -1;
   if( !fwrite( &( image.tl_y ), sizeof( double ), 1, p_file ) )
      return -2;
   if( !fwrite( &( image.br_x ), sizeof( double ), 1, p_file ) )
      return -3;
   if( !fwrite( &( image.br_y ), sizeof( double ), 1, p_file ) )
      return -4;

   return 0;
}

void raft_image_set_corner(raft_image *image, 
			   double tl_x, 
			   double tl_y,
			   double br_x,
			   double br_y)
{
  image->tl_x = tl_x;
  image->tl_y = tl_y;
  image->br_x = br_x;
  image->br_y = br_y;
}
