#include "raft_matrix.h"
#include <stdio.h>
#include <stdlib.h>

/*=====================================================*/

/*!
 * \brief Return a submatrix with FORTRAN standard.
 * \param p_data vector data.
 * \param lines number of lines.
 * \param columns number of columns.
 * \param line_stride stride for lines.
 */

raft_matrix raft_matrix_column_major_emb(double * p_data, 
					 int lines, 
					 int columns, 
					 int line_stride)
{
   raft_matrix result;

   result.p_data = p_data;
   result.lines = lines;
   result.line_stride = line_stride;
   result.columns = columns;
   result.column_stride = 1;

   return result;
}

/*=====================================================*/

/*!
 * \brief Return a submatrix with FORTRAN standard. 
 * \param p_data vector data.
 * \param lines number of lines.
 * \param columns number of columns.
 */

raft_matrix raft_matrix_column_major(double * p_data, 
				     int lines, 
				     int columns)
{
   return raft_matrix_column_major_emb(p_data, 
				       lines, 
				       columns, 
				       lines);
}

/*=====================================================*/

/*!
 * \brief Return a submatrix with FORTRAN standard. 
 * \param p_data vector data.
 * \param lines number of lines.
 * \param columns number of columns.
 * \param column_stride stride for columns.
 */

raft_matrix raft_matrix_row_major_emb(double * p_data, 
				      int lines, 
				      int columns, 
				      int column_stride)
{
   raft_matrix result;

   result.p_data = p_data;
   result.lines = lines;
   result.line_stride = 1;
   result.columns = columns;
   result.column_stride = column_stride;

   return result;
}

/*=====================================================*/

/*!
 * \brief Return a submatrix with FORTRAN standard.  
 * \param p_data vector data.
 * \param lines number of lines.
 * \param columns number of columns.
 */

raft_matrix raft_matrix_row_major(double * p_data, 
				  int lines, 
				  int columns)
{
  return raft_matrix_row_major_emb(p_data, 
				   lines, 
				   columns, 
				   columns);
}

/*=====================================================*/

/*!
 * \brief Read binary matrix from file
 * \param p_file file pointer. 
 */

raft_matrix raft_matrix_read(FILE * p_file)
{
   // Resulting matrix:
   raft_matrix result;

   // Read count:
   int count;

   // Read number of lines:
   if ( !fread( &(result.lines), sizeof( int ), 1, p_file ) )
      return RAFT_MATRIX_EMPTY;

   // Read number of columns:
   if ( !fread( &(result.columns), sizeof( int ), 1, p_file ) )
      return RAFT_MATRIX_EMPTY;

   // Alloc memory:
   result.p_data = malloc( result.lines * result.columns * sizeof( double ) );
   if ( !( result.p_data ) )
      return RAFT_MATRIX_EMPTY;

   // Read data:
   count = fread( result.p_data, sizeof( double ), result.lines * result.columns, p_file );
   if ( count != ( result.lines * result.columns ) )
   {
      raft_matrix_destroy( &result );
      return( result );
   }

   // Data read successfull! Write remaining fields:
   result.line_stride = 1;
   result.column_stride = result.columns;

   return result;
}

/*=====================================================*/

/*!
 * \brief Read ascii matrix from file
 * \param p_file file pointer. 
 * \param lines number of lines.
 * \param columns number of columns.
 */

raft_matrix raft_matrix_scan(FILE *p_file,
			     int lines,
			     int columns)
{
  raft_matrix result;
  int i,j;
  
  if (!p_file) 
    return RAFT_MATRIX_EMPTY;
  
  // Alloc memory:
  result = raft_matrix_create(lines, columns);
  
  if ( !result.p_data )
    return RAFT_MATRIX_EMPTY;
  
  // Read data:
  for (i = 0 ; i < lines; i ++)
    {
      for (j = 0 ; j < columns; j++ )
	{
	  fscanf(p_file,"%lf", &raft_matrix_element(result, i, j));
	}
    }
  
  return result;
}

/*=====================================================*/

/*!
 * \brief Transpose matrix
 * \param mat raft matrix. 
 */

raft_matrix raft_matrix_transpose( raft_matrix mat )
{
   raft_matrix result;

   result.p_data = mat.p_data;
   result.lines = mat.columns;
   result.line_stride = mat.column_stride;
   result.columns = mat.lines;
   result.column_stride = mat.line_stride;

   return result;
}

/*=====================================================*/

/*!
 * \brief Flips matrix upside down
 * \param mat raft matrix.
 */

raft_matrix raft_matrix_flipud( raft_matrix mat )
{
   raft_matrix result;

   result.p_data = mat.p_data + ( mat.lines - 1 ) * mat.column_stride;
   result.lines = mat.lines;
   result.line_stride = mat.line_stride;
   result.columns = mat.columns;
   result.column_stride = -mat.column_stride;

   return result;
}

/*=====================================================*/

/*!
 * \brief Create matrix
 * \param lines number of lines.
 * \param columns number of columns.
 */

raft_matrix raft_matrix_create(int lines, 
			       int columns)
{
   raft_matrix result;

   result.p_data = malloc( lines * columns * sizeof( double ) );
   if( !( result.p_data ) )
      return RAFT_MATRIX_EMPTY;

   result.lines = lines;
   result.line_stride = 1;
   result.columns = columns;
   result.column_stride = columns;

   return result;
}

/*=====================================================*/

/*!
 * \brief Destroy matrix
 * \param p_mat raft matrix.
 */

void raft_matrix_destroy(raft_matrix * p_mat)
{
   free( p_mat->p_data );
   *p_mat = RAFT_MATRIX_EMPTY;
}

/*=====================================================*/

/*!
 * \brief Write binary matrix at file
 * \param mat raft matrix.
 * \param p_file file pointer.
 */

int raft_matrix_write(raft_matrix mat, 
		      FILE * p_file)
{
   // Counters:
   int i, j, count;

   // Write number of lines:
   count = fwrite( &( mat.lines ), sizeof( int ), 1, p_file );
   if ( !count )
      return 1;

   // Write number of columns:
   count = fwrite( &( mat.columns ), sizeof( int ), 1, p_file );
   if ( !count )
      return 2;

   // Write out data:
   for ( i = 0; i < raft_matrix_nlines(mat); ++i )
     for ( j = 0; j < raft_matrix_ncolumns(mat); ++j )
      {
         count = fwrite( &raft_matrix_element( mat, i, j ), sizeof( double ), 1, p_file );
         if ( !count )
            return 3 + i * mat.columns + j;
      }

   return 0;
}

/*=====================================================*/

/*!
 * \brief Print ascii matrix at file
 * \param mat raft matrix.
 * \param p_file file pointer.
 * \param lines number of lines.
 * \param columns number of columns.
 */

int raft_matrix_print(raft_matrix mat, 
		      FILE *p_file)
{
  int i, j, count;

  if (!p_file) 
    return 1;
  
  // Write out data:
  for ( i = 0; i < raft_matrix_nlines(mat); ++i )
    {
      for ( j = 0; j < raft_matrix_ncolumns(mat); ++j )
	{
	  count = fprintf(p_file, "%lf ", raft_matrix_element(mat, i, j));
	}
      fprintf(p_file,"\n");
    }
  
  return 0;
}

/*=====================================================*/

/*!
 * \brief Extract line from matrix
 * \param mat raft matrix.
 * \param line line number.
 */

raft_vector raft_matrix_get_line(raft_matrix mat, 
				 int line)
{
   raft_vector result;
   result.p_data = mat.p_data + line * mat.column_stride;
   result.size = mat.columns;
   result.stride = mat.line_stride;
   return result;
}

/*=====================================================*/

/*!
 * \brief Extract column from matrix
 * \param mat raft matrix.
 * \param column matrix.
 */

raft_vector raft_matrix_get_column(raft_matrix mat,
				   int column)
{
   raft_vector result;

   result.p_data = mat.p_data + column * mat.line_stride;
   result.size = mat.lines;
   result.stride = mat.column_stride;

   return result;
}
