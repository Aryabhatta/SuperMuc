/*
 * relax_jacobi.c
 *
 * Jacobi Relaxation
 *
 */
#include <stdbool.h>
#include "heat.h"
#include <omp.h>

/*
 * Residual (length of error vector)
 * between current solution and next after a Jacobi step
 */

// ****NOTE: This function is not used anymore****
double residual_jacobi( double *u, 
			unsigned sizex, unsigned sizey )
{
    unsigned i, j;
    double unew, diff, sum=0.0;
/*
    for( i=1; i<sizey-1; i++ )
    {
	for( j=1; j<sizex-1; j++ )
        {
	    unew = 0.25 * (u[ i*sizex     + (j-1) ]+  // left
			   u[ i*sizex     + (j+1) ]+  // right
			   u[ (i-1)*sizex + j     ]+  // top
			   u[ (i+1)*sizex + j     ]); // bottom

	    diff = unew - u[i*sizex + j];
	    sum += diff * diff; 
	}
    }
*/
    return sum;
}
int min( int a , int b )
{
   if(a < b)
	return a;
   else
	return b;
}

/*
 * One Jacobi iteration step
 */
double relax_jacobi( double *u, double *utmp,
		   unsigned sizex, unsigned sizey )
{
  int i, j, b = 0;
  double diff, sum = 0.0;
    
  // Blocking of loops - split into horizontal 4
  int CurIter = 1, NoBlocks = 4;
  int blockSz = (sizey/4);
  int BlockIterNo[4] = {0,0,0,0};

  int NoIter = 4;

  for( CurIter = 1 ; CurIter <= NoIter ; CurIter ++) 
  {
    for( b= CurIter-1; b >= 0 ; b--) // goes from 0 to 3
    {
      	BlockIterNo[b] += 1;

      	if( BlockIterNo[b]%2 ==1 )// odd iteration
    	{
         omp_set_num_threads(32);
         #pragma omp parallel
         {
          #pragma omp for private (diff,j) reduction(+:sum)
          
       	   for( i= b*blockSz + 1; i<=min( (b+1)*blockSz, sizey-2 ); i++ )
           { 
              for( j=1; j<sizex-1; j++ )
              {
              	utmp[i*sizex+j]= 0.25 * ( u[ i*sizex     + (j-1) ]+  // left
                	                  u[ i*sizex     + (j+1) ]+  // right
                        	          u[ (i-1)*sizex + j     ]+  // top
                                	  u[ (i+1)*sizex + j     ]); // bottom
	        diff = utmp[ i*sizex+j] - u[i*sizex+j];
                sum += diff * diff;
              }
           }
         }
        }
     
        else // even iteration
        {
          omp_set_num_threads(32);
          #pragma omp parallel
          {
          #pragma omp for private (diff,j) reduction(+:sum)

           for( i= b*blockSz + 1; i<=min( (b+1)*blockSz, sizey-2); i++ )
           {
              for( j=1; j<sizex-1; j++ )
              {
            	u[i*sizex+j]= 0.25 * (utmp[ i*sizex     + (j-1) ]+  // left
                                     utmp[ i*sizex     + (j+1) ]+  // right
                                     utmp[ (i-1)*sizex + j     ]+  // top
                                     utmp[ (i+1)*sizex + j     ]); // bottom

           	diff = u[ i*sizex+j] - utmp[i*sizex+j];
           	sum += diff * diff;
              }
	   }
          }
        }
    }
  }

  // lower triangle of Interleaving
  // for( CurIter = 2 ; CurIter <= NoIter ; CurIter ++) 
  for( CurIter = 2 ; CurIter <= NoBlocks ; CurIter ++) 
  {
     for( b=NoBlocks-1; b >= CurIter-1 ; b--) // goes from 4 to 2, then 4 to 3, then 4 to 4
     {
    	BlockIterNo[b] += 1;

    	if( BlockIterNo[b]%2 == 1 )// odd iteration
    	{
             omp_set_num_threads(32);
         #pragma omp parallel
         {
          #pragma omp for private (diff,j) reduction(+:sum)


       	   for( i= b*blockSz + 1; i<=min( (b+1)*blockSz, sizey-2); i++ )
           {
          	for( j=1; j<sizex-1; j++ )
          	{
            	   utmp[i*sizex+j]= 0.25 * (u[ i*sizex     + (j-1) ]+  // left
                                            u[ i*sizex     + (j+1) ]+  // right
                                            u[ (i-1)*sizex + j     ]+  // top
                                            u[ (i+1)*sizex + j     ]); // bottom

	           diff = utmp[ i*sizex+j] - u[i*sizex+j];
       		   sum += diff * diff;
          	}
	   }
         }
        }  
        else // even iteration
        {
             omp_set_num_threads(32);
         #pragma omp parallel
         {
          #pragma omp for private (diff,j) reduction(+:sum)


       	   for( i= b*blockSz + 1; i<=min( (b+1)*blockSz,sizey-2); i++ )
           {
          	for( j=1; j<sizex-1; j++ )
          	{
            	   u[i*sizex+j]= 0.25 * (utmp[ i*sizex     + (j-1) ]+  // left
                                         utmp[ i*sizex     + (j+1) ]+  // right
                                         utmp[ (i-1)*sizex + j     ]+  // top
                                         utmp[ (i+1)*sizex + j     ]); // bottom

	           diff = u[ i*sizex+j] - utmp[i*sizex+j];
        	   sum += diff * diff;
          	}
           }
          }
        }    
      }
  }
  
  // u has most updated value
  // no need to copy utmp to u
  return sum;
}


  //for( b= (CuIter-1)>(NoBlocks-1)?(NoBlocks-1):(CurIter-1); b >= 0 ; b--) // goes from 0 to 3
