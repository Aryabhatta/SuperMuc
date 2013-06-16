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

   // code generalised
   if( yDim == 1 )
   {
     xStart =  ((sizey-2)/xDim) * rank + 1;
     xEnd = ((sizey-2)/xDim) * (rank+1) ;
     yStart = 1;
     yEnd = sizex-2;
   }

   // code generalised 
  if( xDim == 1 )
   {
     yStart =  (sizex-2)/yDim * rank + 1;
     yEnd = (sizex-2)/yDim * (rank+1);
     xStart = 1;
     xEnd = sizey-2;
   }
   
    if(xDim == yDim ) // topology p * p
    {
        int stepX = (sizey-2)/xDim;

        xStart = stepX * (rank/yDim) + 1;
        xEnd = stepX * ( (rank/yDim) + 1 );
        
        int stepY = (sizex-2)/yDim;
        
        yStart = stepY * (rank % xDim) + 1;
        yEnd = stepY * ( (rank % xDim)+1);
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

   // generalised code
   if( yDim == 1 )
   {
     xStart =  ((sizey-2)/xDim) * rank + 2;
     xEnd = ((sizey-2)/xDim) * (rank+1) - 1;
     yStart = 1;
     yEnd = sizex-2;
   }
 
  // generalised code
  if( xDim == 1 )
   {
     yStart =  (sizex-2)/yDim * rank + 2;
     yEnd = (sizex-2)/yDim * (rank+1) - 1;
     xStart = 1;
     xEnd = sizey-2;
   }
   
    if( xDim == yDim )
    {
        int stepX = (sizey-2)/xDim;

        xStart = stepX * (rank/yDim) + 1;
        xEnd = stepX * ( (rank/yDim) + 1 );
        
        int stepY = (sizex-2)/yDim;
        
        yStart = stepY * (rank % xDim) + 1;
        yEnd = stepY * ( (rank % xDim)+1);

        // since inner
        xStart = xStart + 1;
        xEnd = xEnd - 1;
        yStart = yStart + 1;
        yEnd = yEnd - 1;
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

   // generalised code
   if( yDim == 1 )
   {
     xStart =  ((sizey-2)/xDim) * rank + 1;
     xEnd = ((sizey-2)/xDim) * (rank+1);
     yStart = 1;
     yEnd = sizex-2;
   }
 
   // generalised code
  if( xDim == 1 )
   {
     yStart =  (sizex-2)/yDim * rank + 1;
     yEnd = (sizex-2)/yDim * (rank+1);
     xStart = 1;
     xEnd = sizey-2;
   }

    if(xDim == yDim )
    {
        int stepX = (sizey-2)/xDim;

        xStart = stepX * (rank/yDim) + 1;
        xEnd = stepX * ( (rank/yDim) + 1 );
        
        int stepY = (sizex-2)/yDim;
        
        yStart = stepY * (rank % xDim) + 1;
        yEnd = stepY * ( (rank % xDim)+1);
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
