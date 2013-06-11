/*
 * relax_jacobi.c
 *
 * Jacobi Relaxation
 *
 */
#include <stdbool.h>
#include "heat.h"

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

    return sum;
}


/*
 * One Jacobi iteration step
 */
double relax_jacobi( double **u, double **utmp,
	   unsigned sizex, unsigned sizey, int xDim, int yDim, int rank )
{
   int i, j;
   double diff, sum = 0.0;
   int xStart, xEnd, yStart, yEnd;

   if( yDim == 1 )
   {
     xStart =  ((sizey-2)/xDim) * rank + 1;
     xEnd = ((sizey-2)/xDim) * (rank+1) ;
     yStart = 1;
     yEnd = sizex-2;
   }
 
  if( xDim == 1 )
   {
     yStart =  (sizex-2)/yDim * rank + 1;
     yEnd = (sizex-2)/yDim * (rank+1);
     xStart = 1;
     xEnd = sizey-2;
   }
   
    if(xDim == 2 && yDim == 2 )
    {
      if( rank == 0 )
      {
        xStart = (sizey-2)/xDim * rank + 1;
        xEnd = (sizey-2)/xDim * (rank+1);
        yStart = (sizex-2)/yDim * rank + 1;
        yEnd = (sizex-2)/yDim * (rank+1);
      }
      else if( rank == 1)
      {
        xStart = (sizey-2)/xDim * (rank-1) + 1;
        xEnd = (sizey-2)/xDim * rank;
        yStart = (sizex-2)/yDim * rank + 1;
        yEnd = (sizex-2)/yDim * (rank+1);
      }
      else if( rank == 3 )
      {
        xStart = (sizey-2)/xDim * (rank-2) + 1;
        xEnd = (sizey-2)/xDim * (rank-1);
        yStart = (sizex-2)/yDim * (rank-2) + 1;
        yEnd = (sizex-2)/yDim * (rank-1);
      }
      else if( rank == 2)
      {
        xStart = (sizey-2)/xDim * (rank-1) + 1;
        xEnd = (sizey-2)/xDim * (rank);
        yStart = (sizex-2)/yDim * (rank-2) + 1;
        yEnd = (sizex-2)/yDim * (rank-1);
      }
    }

    // relaxation
    for( i=xStart; i<=xEnd; i++ )
    {
        for( j=yStart; j<=yEnd; j++ )
        {
            *(*utmp + (i*sizex+j))= 0.25 * (
                                     *(*u + (i*sizex+(j-1)) ) +  // left
                                     *(*u + (i*sizex+(j+1)) ) +  // right
                                     *(*u + ((i-1)*sizex+j) ) +  // top
                                     *(*u + ((i+1)*sizex+j) ));  // bottom

         diff = *(*utmp + (i*sizex+j)) - *(*u + (i*sizex+j));
         sum += diff * diff;
        }
    }

    // avoid copy operation
    double *temp = *u;
    *u = *utmp;
    *utmp = temp;

   return sum;
}

double relax_jacobiInner( double **u, double **utmp,
	   unsigned sizex, unsigned sizey, int xDim, int yDim, int rank )
{
   int i, j;
   double diff, sum = 0.0;
   int xStart, xEnd, yStart, yEnd;

   if( yDim == 1 )
   {
     xStart =  ((sizey-2)/xDim) * rank + 2;
     xEnd = ((sizey-2)/xDim) * (rank+1) - 1;
     yStart = 1;
     yEnd = sizex-2;
   }
 
  if( xDim == 1 )
   {
     yStart =  (sizex-2)/yDim * rank + 2;
     yEnd = (sizex-2)/yDim * (rank+1) - 1;
     xStart = 1;
     xEnd = sizey-2;
   }
   
    if(xDim == 2 && yDim == 2 )
    {
      if( rank == 0 )
      {
        xStart = (sizey-2)/xDim * rank + 2;
        xEnd = (sizey-2)/xDim * (rank+1)-1;
        yStart = (sizex-2)/yDim * rank + 2;
        yEnd = (sizex-2)/yDim * (rank+1) - 1;
      }
      else if( rank == 1)
      {
        xStart = (sizey-2)/xDim * (rank-1) + 2;
        xEnd = (sizey-2)/xDim * rank - 1;
        yStart = (sizex-2)/yDim * rank + 2;
        yEnd = (sizex-2)/yDim * (rank+1) - 1 ;
      }
      else if( rank == 3 )
      {
        xStart = (sizey-2)/xDim * (rank-2) + 2;
        xEnd = (sizey-2)/xDim * (rank-1) -1 ;
        yStart = (sizex-2)/yDim * (rank-2) + 2;
        yEnd = (sizex-2)/yDim * (rank-1)-1;
      }
      else if( rank == 2)
      {
        xStart = (sizey-2)/xDim * (rank-1) + 2;
        xEnd = (sizey-2)/xDim * (rank)-1;
        yStart = (sizex-2)/yDim * (rank-2) + 2;
        yEnd = (sizex-2)/yDim * (rank-1)-1;
      }
    }

    // relaxation
    for( i=xStart; i<=xEnd; i++ )
    {
        for( j=yStart; j<=yEnd; j++ )
        {
            *(*utmp + (i*sizex+j))= 0.25 * (
                                     *(*u + (i*sizex+(j-1)) ) +  // left
                                     *(*u + (i*sizex+(j+1)) ) +  // right
                                     *(*u + ((i-1)*sizex+j) ) +  // top
                                     *(*u + ((i+1)*sizex+j) ));  // bottom

         diff = *(*utmp + (i*sizex+j)) - *(*u + (i*sizex+j));
         sum += diff * diff;
        }
    }

   return sum;
}

double relax_jacobiBoundary( double **u, double **utmp,
	   unsigned sizex, unsigned sizey, int xDim, int yDim, int rank )
{
   int i, j;
   double diff, sum = 0.0;
   int xStart, xEnd, yStart, yEnd;

   if( yDim == 1 )
   {
     xStart =  ((sizey-2)/xDim) * rank + 1;
     xEnd = ((sizey-2)/xDim) * (rank+1);
     yStart = 1;
     yEnd = sizex-2;
   }
 
  if( xDim == 1 )
   {
     yStart =  (sizex-2)/yDim * rank + 1;
     yEnd = (sizex-2)/yDim * (rank+1);
     xStart = 1;
     xEnd = sizey-2;
   }

    if(xDim == 2 && yDim == 2 )
    {
      if( rank == 0 )
      {
        xStart = (sizey-2)/xDim * rank + 1;
        xEnd = (sizey-2)/xDim * (rank+1);
        yStart = (sizex-2)/yDim * rank + 1;
        yEnd = (sizex-2)/yDim * (rank+1);
      }
      else if( rank == 1)
      {
        xStart = (sizey-2)/xDim * (rank-1) + 1;
        xEnd = (sizey-2)/xDim * rank;
        yStart = (sizex-2)/yDim * rank + 1;
        yEnd = (sizex-2)/yDim * (rank+1);
      }
      else if( rank == 3 )
      {
        xStart = (sizey-2)/xDim * (rank-2) + 1;
        xEnd = (sizey-2)/xDim * (rank-1);
        yStart = (sizex-2)/yDim * (rank-2) + 1;
        yEnd = (sizex-2)/yDim * (rank-1);
      }
      else if( rank == 2)
      {
        xStart = (sizey-2)/xDim * (rank-1) + 1;
        xEnd = (sizey-2)/xDim * (rank);
        yStart = (sizex-2)/yDim * (rank-2) + 1;
        yEnd = (sizex-2)/yDim * (rank-1);
      }
    }

    // relaxation working for x boundaries
    for( i=xStart; i<=xEnd; i+= (xEnd-xStart) )
    {
        for( j=yStart; j<=yEnd; j++ )
        {
            *(*utmp + (i*sizex+j))= 0.25 * (
                                     *(*u + (i*sizex+(j-1)) ) +  // left
                                     *(*u + (i*sizex+(j+1)) ) +  // right
                                     *(*u + ((i-1)*sizex+j) ) +  // top
                                     *(*u + ((i+1)*sizex+j) ));  // bottom

         diff = *(*utmp + (i*sizex+j)) - *(*u + (i*sizex+j));
         sum += diff * diff;
        }
    }

    // relaxation working for y boundaries
    for( i=xStart; i<=xEnd; i++ )
    {
        for( j=yStart; j<=yEnd; j+=(yEnd-yStart) )
        {
            *(*utmp + (i*sizex+j))= 0.25 * (
                                     *(*u + (i*sizex+(j-1)) ) +  // left
                                     *(*u + (i*sizex+(j+1)) ) +  // right
                                     *(*u + ((i-1)*sizex+j) ) +  // top
                                     *(*u + ((i+1)*sizex+j) ));  // bottom

         diff = *(*utmp + (i*sizex+j)) - *(*u + (i*sizex+j));
         sum += diff * diff;
        }
    }

    // Swap u & utmp only after boundary relaxation is done
    double *temp = *u;
    *u = *utmp;
    *utmp = temp;

   return sum;
}
