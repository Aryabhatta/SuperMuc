/*
 * relax_gauss.c
 *
 * Gauss-Seidel Relaxation
 *
 */

#include "heat.h"

/*
 * Residual (length of error vector)
 * between current solution and next after a Gauss-Seidel step
 *
 * Temporary array utmp needed to not change current solution
 *
 * Flop count in inner body is 7
 */

double residual_gauss( double *u, double *utmp,
		       unsigned sizex, unsigned sizey )
{
    unsigned i, j;
    double unew, diff, sum=0.0;

    // first row (boundary condition) into utmp
    for( j=1; j<sizex-1; j++ )
	utmp[0*sizex + j] = u[0*sizex + j];
    // first column (boundary condition) into utmp
    for( i=1; i<sizey-1; i++ )
	utmp[i*sizex + 0] = u[i*sizex + 0];

    for( i=1; i<sizey-1; i++ )
    {
	for( j=1; j<sizex-1; j++ )
	{
	    unew = 0.25 * (utmp[ i*sizex     + (j-1) ]+  // new left
			   u   [ i*sizex     + (j+1) ]+  // right  
			   utmp[ (i-1)*sizex + j     ]+  // new top
			   u   [ (i+1)*sizex + j     ]); // bottom

	    diff = unew - u[i*sizex + j];
	    sum += diff * diff; 

	    utmp[i*sizex + j] = unew;
	}
    }

    return sum;
}


/*
 * One Gauss-Seidel iteration step
 *
 * Flop count in inner body is 4
 */
double relax_gauss( double *u, 
		  unsigned sizex, unsigned sizey, int xDim, int yDim, int rank )
{
    int i, j;
    int red = 0;  
    double uold, diff, sum=0.0;

   int xStart, xEnd, yStart, yEnd;

   if( yDim == 1 )
   {
     xStart =  ((sizey-2)/xDim) * rank + 1;
     xEnd = ((sizey-2)/xDim) * (rank+1) ;
     yStart = 1;
     yEnd = sizex-2;
   }
   else if( xDim == 1 )
   {
     yStart =  (sizex-2)/yDim * rank + 1;
     yEnd = (sizex-2)/yDim * (rank+1);
     xStart = 1;
     xEnd = sizey-2;
   }
   else if(xDim == yDim ) // topology p*p
    {
        int stepX = (sizey-2)/xDim;

        xStart = stepX * (rank/yDim) + 1;
        xEnd = stepX * ( (rank/yDim) + 1 );

        int stepY = (sizex-2)/yDim;

        yStart = stepY * (rank % xDim) + 1;
        yEnd = stepY * ( (rank % xDim)+1);
    }

    // Red Black approach
    if( red == 0 )
    {
        for( i=xStart; i<=xEnd; i++ )
        {
    	    for( j=yStart; j<=yEnd; j++ )
	        {
                if( (i+j)%2 == 0)
                {
                  uold = u[i*sizex+j];
	              u[i*sizex+j]= 0.25 * (u[ i*sizex     + (j-1) ]+
            			          	    u[ i*sizex     + (j+1) ]+
			            	            u[ (i-1)*sizex + j     ]+
    				                    u[ (i+1)*sizex + j     ]);

	              diff = u[i*sizex + j] - uold;
                  sum += diff * diff; 
                }
	        }
        }
    }

    red = 1;
    
    // Red Black approach
    if( red == 1 )
    {
        for( i=xStart; i<=xEnd; i++ )
        {
    	    for( j=yStart; j<=yEnd; j++ )
	        {
                if( (i+j)%2 == 1)
                {
                  uold = u[i*sizex+j];
	              u[i*sizex+j]= 0.25 * (u[ i*sizex     + (j-1) ]+
			          	    u[ i*sizex     + (j+1) ]+
				            u[ (i-1)*sizex + j     ]+
    				        u[ (i+1)*sizex + j     ]);

	              diff = u[i*sizex + j] - uold;
                  sum += diff * diff; 
                }
	        }
        }
    }
    return sum;
}

double relax_gaussInner( double *u, 
		  unsigned sizex, unsigned sizey, int xDim, int yDim, int rank )
{
    unsigned i, j;
    int red = 0;  
    double uold, diff, sum=0.0;
   int xStart, xEnd, yStart, yEnd;

   if( yDim == 1 )
   {
     xStart =  ((sizey-2)/xDim) * rank + 2;
     xEnd = ((sizey-2)/xDim) * (rank+1) - 1;
     yStart = 1;
     yEnd = sizex-2;
   }
   else if( xDim == 1 )
   {
     yStart =  (sizex-2)/yDim * rank + 2;
     yEnd = (sizex-2)/yDim * (rank+1) - 1;
     xStart = 1;
     xEnd = sizey-2;
   }
   else if(xDim == yDim ) //topology p*p
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


    if( red == 0 )
    {
        for( i=xStart; i<=xEnd; i++ )
        {
    	    for( j=yStart; j<=yEnd; j++ )
	        {
                if( (i+j)%2 == 0)
                {
                  uold = u[i*sizex+j];
	              u[i*sizex+j]= 0.25 * (u[ i*sizex     + (j-1) ]+
			          	    u[ i*sizex     + (j+1) ]+
				            u[ (i-1)*sizex + j     ]+
    				        u[ (i+1)*sizex + j     ]);

	              diff = u[i*sizex + j] - uold;
                  sum += diff * diff; 
                }
	        }
        }
    }

    red = 1;

    if( red == 1 )
    {
        for( i=xStart; i<=xEnd; i++ )
        {
    	    for( j=yStart; j<=yEnd; j++ )
	        {
                if( (i+j)%2 == 1)
                {
                  uold = u[i*sizex+j];
	              u[i*sizex+j]= 0.25 * (u[ i*sizex     + (j-1) ]+
			          	    u[ i*sizex     + (j+1) ]+
				            u[ (i-1)*sizex + j     ]+
    				        u[ (i+1)*sizex + j     ]);

	              diff = u[i*sizex + j] - uold;
                  sum += diff * diff; 
                }
	        }
        }
    }
    return sum;
}

double relax_gaussBoundary( double *u, 
		  unsigned sizex, unsigned sizey, int xDim, int yDim, int rank )
{
    unsigned i, j;
    int red = 0;  
    double uold, diff, sum=0.0;

   int xStart, xEnd, yStart, yEnd;

   if( yDim == 1 )
   {
     xStart =  ((sizey-2)/xDim) * rank + 1;
     xEnd = ((sizey-2)/xDim) * (rank+1) ;
     yStart = 1;
     yEnd = sizex-2;
   }
   else if( xDim == 1 )
   {
     yStart =  (sizex-2)/yDim * rank + 1;
     yEnd = (sizex-2)/yDim * (rank+1);
     xStart = 1;
     xEnd = sizey-2;
   }
   else if(xDim == yDim ) // topology p*p
    {
        int stepX = (sizey-2)/xDim;

        xStart = stepX * (rank/yDim) + 1;
        xEnd = stepX * ( (rank/yDim) + 1 );

        int stepY = (sizex-2)/yDim;

        yStart = stepY * (rank % xDim) + 1;
        yEnd = stepY * ( (rank % xDim)+1);
    }

    if( red == 0 )
    {
        for( i=xStart; i<=xEnd; i+=(xEnd-xStart) )
        {
    	    for( j=yStart; j<=yEnd; j++ )
	        {
                if( (i+j)%2 == 0)
                {
                  uold = u[i*sizex+j];
	              u[i*sizex+j]= 0.25 * (u[ i*sizex     + (j-1) ]+
			          	    u[ i*sizex     + (j+1) ]+
				            u[ (i-1)*sizex + j     ]+
    				        u[ (i+1)*sizex + j     ]);

	              diff = u[i*sizex + j] - uold;
                  sum += diff * diff; 
                }
	        }
        }
        for( i=xStart; i<=xEnd; i++ )
        {
    	    for( j=yStart; j<=yEnd; j+=(yEnd-yStart) )
	        {
                if( (i+j)%2 == 0)
                {
                  uold = u[i*sizex+j];
	              u[i*sizex+j]= 0.25 * (u[ i*sizex     + (j-1) ]+
			          	    u[ i*sizex     + (j+1) ]+
				            u[ (i-1)*sizex + j     ]+
    				        u[ (i+1)*sizex + j     ]);

	              diff = u[i*sizex + j] - uold;
                  sum += diff * diff; 
                }
	        }
        }
    }

    red = 1;

    if( red == 1 )
    {
        for( i=xStart; i<=xEnd; i+=(xEnd-xStart) )
        {
    	    for( j=yStart; j<=yEnd; j++ )
	        {
                if( (i+j)%2 == 1)
                {
                  uold = u[i*sizex+j];
	              u[i*sizex+j]= 0.25 * (u[ i*sizex     + (j-1) ]+
			          	    u[ i*sizex     + (j+1) ]+
				            u[ (i-1)*sizex + j     ]+
    				        u[ (i+1)*sizex + j     ]);

	              diff = u[i*sizex + j] - uold;
                  sum += diff * diff; 
                }
	        }
        }
        for( i=xStart; i<=xEnd; i++ )
        {
    	    for( j=yStart; j<=yEnd; j+=(yEnd-yStart) )
	        {
                if( (i+j)%2 == 1)
                {
                  uold = u[i*sizex+j];
	              u[i*sizex+j]= 0.25 * (u[ i*sizex     + (j-1) ]+
			          	    u[ i*sizex     + (j+1) ]+
				            u[ (i-1)*sizex + j     ]+
    				        u[ (i+1)*sizex + j     ]);

	              diff = u[i*sizex + j] - uold;
                  sum += diff * diff; 
                }
	        }
        }
    }
    return sum;
}
