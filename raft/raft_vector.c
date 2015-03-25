#include "raft_vector.h"
#include <stdlib.h>

/*=====================================================*/

/*!
 * \brief Create a raft vector.
 * \param size vector size.
 */

raft_vector raft_vector_create(int size)
{
   raft_vector result;

   result.p_data = malloc( size * sizeof( double ) );

   if ( !result.p_data )
   {
     result.size = 0;
     result.stride = 0;
     return RAFT_VECTOR_EMPTY;
   }
   else
   {
     result.size = size;
     result.stride = 1;
   }

   return result;
}

/*=====================================================*/

/*!
 * \brief Destroy a raft vector.
 * \param p_vec raft vector.
 */

void raft_vector_destroy(raft_vector * p_vec)
{
   free( p_vec->p_data );
   p_vec->p_data = 0;
   p_vec->size = 0;
   p_vec->stride = 0;
}

/*=====================================================*/

/*!
 * \brief Scan an ascii vector.
 * \param p_vec raft vector.
 * \param size vector size.
 */

raft_vector raft_vector_scan(FILE *fp,
			     int size)
{
  raft_vector result;
  int i, count;
  
  if (!fp) 
    return RAFT_VECTOR_EMPTY;
  
  // Alloc memory:
  result = raft_vector_create(size);
  
  if ( !result.p_data )
    return RAFT_VECTOR_EMPTY;
  
  // Read data:
  for (i = 0 ; i < size; i ++)
    {
      count = fscanf(fp,"%lf", &raft_vector_element(result, i));
    }
  
  return result;
}

/*=====================================================*/

/*!
 * \brief Print an ascii vector.
 * \param p_vec raft vector.
 */

int raft_vector_print(FILE *fp,
		      raft_vector vec)
{
  int i, count;

  if (!fp) 
    return 1;
  
  // Write out data:
  for ( i = 0; i < raft_vector_size(vec); ++i )
    {
      count = fprintf(fp, "%lf ", raft_vector_element(vec, i));
    }
  
  return 0;
}
