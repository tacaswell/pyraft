#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../raft/raft_image.h"

#define MAX(x,y) (((x)>(y)) ? (x):(y))
#define MIN(x,y) (((x)<(y)) ? (x):(y)) 

void getFluorFromSpec(raft_image sino,
		      raft_image table, //spec table from XFCT 
	              int s,      //slice number 
		      int x,      //number of points: x 
		      int y,      //number of points: y
		      int a,      //number of angles: a
		      int t,      //total number of columns @ file
		      int col)    //column to extract fluorescence data
{
  int k, j, i, b, total;
  double q;

  total = x*y*a;
  
  /* building slice: s=0,1,2,3,4,... */
  s++;

  for(k=0; k < a; k++)
    {
      if(k%2==0) // even
	{
	  b =  k*x*y + (s-1)*x;
	  
          i = 0;
	  for( j = b; j < b + x; j++)
	    {
	      q = raft_matrix_element(table.data, j, col); // (table[j][col]);
 
	      raft_matrix_element(sino.data,k,i) = q;

	      i++;
	    }
	}
      else
	{
	  b =  (k+1)*x*y + (s-1)*x;
	  
          i = 0;
	  for(j = b; j < b + x; j++)
	    {
	      q = raft_matrix_element(table.data, j, col); // (table[j][col]);

	      raft_matrix_element(sino.data,k,i) = q;

	      i++;
	    }	  
	}
    }

}


/*-----------------------*/

void getTransFromSpec(raft_image sino,
		      raft_image table, //spec table from XFCT 
	              int s,      //slice number 
		      int x,      //number of points: x 
		      int y,      //number of points: y
		      int a,      //number of angles: a
		      int t,      //total number of columns @ file
		      int c,	  //column to extract transmission data
		      int n)	  //column to extract normalization factor
{
  int k, j, i, b, total;
  double q, I0;
  
  total = x*y*a;
  
  /* building slice: s=0,1,2,3,4,... */
  
  s++;

  for(k=0; k < a; k++)
    {
      if(k%2==0) // even
	{
	  b =  k*x*y + (s-1)*x;
	  
          i = 0;

	  I0 = raft_matrix_element(table.data, b, n);

	  for( j = b; j < b + x; j++)
	    {
	      q = raft_matrix_element(table.data, j, c);
 
	      raft_matrix_element(sino.data,k,i) = q/I0;

	      i++;
	    }
	}
      else
	{
	  b =  (k+1)*x*y + (s-1)*x;
	  
          i = 0;

  	  I0 = raft_matrix_element(table.data, b, n);

	  for(j = b; j < b + x; j++)
	    {
	      q = raft_matrix_element(table.data, j, c);

	      raft_matrix_element(sino.data,k,i) = q/I0;

	      i++;
	    }	  
	}
    }

}

