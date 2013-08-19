/*
 * heat.h
 *
 * Iterative solver for heat distribution
 */


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <stdbool.h>

#include "input.h"
#include "heat.h"
#include "timing.h"

#define BLOCK 0

void updateu( double * buf, int size, int row, double * u, int rowStart, int np);
void updateucol( double * buf, int size, int row, double * u, int colStart, int np);
void commBoundary(int xDim, int yDim, int rank, int np, double * u, MPI_Comm comm);

void commBoundaryNB(int xDim, int yDim, int rank, int np, double * u, MPI_Comm comm, double * buf1, MPI_Request * request1, double * buf2, MPI_Request * request2, double * buf3, MPI_Request * request3, double * buf4, MPI_Request * request4);
void updateBoundary( int xDim, int yDim, int rank, int np, double * u, MPI_Comm comm, double * buf1, MPI_Request * request1, double * buf2, MPI_Request * request2, double * buf3, MPI_Request * request3, double * buf4, MPI_Request * request4);

void sendU(int xDim, int yDim, int rank, int np, int NumTask, double * u );
void sendResidual(int rank, int numtask, double * residual);

void usage( char *s )
{
    fprintf(stderr, 
	    "Usage: %s <input file> [X-dim] [Y-dim] [result file] \n\n", s);
}

int main( int argc, char *argv[] )
{
    unsigned iter;
    FILE *infile, *resfile;
    char *resfilename;

    // algorithmic parameters
    algoparam_t param;
    int np, i=0,j=0;

    double runtime, flop;
    double residual;
    double totaltime=0;

    totaltime = wtime();

    // check arguments
    if( argc < 3 )
    {
	usage( argv[0] );
	return 1;
    }

    // check input file
    if( !(infile=fopen(argv[1], "r"))  ) 
    {
	fprintf(stderr, 
		"\nError: Cannot open \"%s\" for reading.\n\n", argv[1]);
      
	usage(argv[0]);
	return 1;
    }
  
    // check result file
    resfilename= (argc>=5) ? argv[4]:"heat.ppm";

    if( !(resfile=fopen(resfilename, "w")) )
    {
	    fprintf(stderr, 
		"\nError: Cannot open \"%s\" for writing.\n\n", 
		resfilename);
      
	    usage(argv[0]);
	    return 1;
    }

    int NumTask = 0, rank = 0;
    int xDim =0, yDim = 0;

    xDim = atoi(argv[2]);
    yDim = atoi(argv[3]);

    // MPI initialise
    int status = MPI_Init( &argc, &argv);

    if( status != MPI_SUCCESS )
    {
	printf("\nError in MPI\n");
	return 0;
    }
    
    MPI_Status stat;
    MPI_Comm comm;

    int dims[2];
    dims[0] = xDim;
    dims[1] = yDim;

    int periods[2]={0, 0};
    bool reorder = false;

    // Creating virtual topology cart
    MPI_Cart_create(MPI_COMM_WORLD,2,dims,periods, reorder, &comm );

    MPI_Comm_size( MPI_COMM_WORLD, &NumTask );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );

    // check input
    if( !read_input(infile, &param) )
    {
	    fprintf(stderr, "\nError: Error parsing input file.\n\n");      
	    usage(argv[0]);
	    return 1;
    }

    print_params(&param);

    // set the visualization resolution
    param.visres = param.max_res;

    param.u     = 0;
    param.uhelp = 0;
    param.uvis  = 0;

    param.act_res = param.initial_res;

    // only for non blocking communication
    #if BLOCK == 0
        MPI_Request request1, request2, request3, request4;
    #endif

    // loop over different resolutions
    while(1) {

	// free allocated memory of previous experiment
	if (param.u != 0)
	    finalize(&param);

	if( !initialize(&param) )
	{
	    fprintf(stderr, "Error in Jacobi initialization.\n\n");
	    usage(argv[0]);
	}

	// full size (param.act_res are only the inner points)
	np = param.act_res + 2;

    #if BLOCK == 0
    double buf1[(np-2)];
    double buf2[(np-2)];
    double buf3[(np-2)];
    double buf4[(np-2)];
    #endif
    
	// starting time
	runtime = wtime();
	residual = 999999999;

	iter = 0;
	while(1) 
    {
        #if BLOCK
        if( xDim != 1 && yDim != 1) // Do nt cal for 1x1 topology
    	    commBoundary( xDim, yDim, rank, np, param.u, comm);
        #endif

        #if BLOCK == 0
        commBoundaryNB( xDim, yDim, rank, np, param.u, comm, buf1, &request1, buf2, &request2, buf3, &request3, buf4, &request4);
        #endif

    	switch( param.algorithm ) 
        {

            case 0: // JACOBI
                    #if BLOCK
      	            //modifed residual calculation
	                residual = relax_jacobi(&param.u, &param.uhelp, np, np, \
                    xDim, yDim, rank);
                    #endif

                    #if BLOCK == 0
                  	//modifed residual calculation
                	residual = relax_jacobiInner(&param.u, &param.uhelp, np, np, \
                           xDim, yDim, rank);
                    #endif
        		    break;

            case 1: // GAUSS
                    #if BLOCK
        		    residual = relax_gauss(param.u, np, np, xDim, yDim, rank);
                    #endif

                    #if BLOCK == 0
        		    residual = relax_gaussInner(param.u, np, np, xDim, yDim, rank);
                    #endif
        		    break;
	    }

        #if BLOCK
	    // Synchronise
	    MPI_Barrier(MPI_COMM_WORLD);
        #endif

	    iter++;
    
        #if BLOCK == 0
        //check for boundary
        updateBoundary( xDim, yDim, rank, np, param.u, comm, buf1, &request1, buf2, &request2, buf3, &request3, buf4, &request4);
        
        if( param.algorithm == 0)
        {
            // relax boundary
            residual += relax_jacobiBoundary(&param.u, &param.uhelp, np, np,xDim, yDim, rank);
        }
        if( param.algorithm == 1 )
        {
            // relax boundary
	        residual += relax_gaussBoundary(param.u, np, np,xDim, yDim, rank);
        }
        #endif

	    // max. iteration reached ? (no limit with maxiter=0)
	    if (param.maxiter>0 && iter>=param.maxiter) break;
	}
	
    // stopping time
	runtime = wtime() - runtime;

	sendResidual(rank, NumTask, &residual);

	if( rank == 0 )
	{
        // change in no of flops/iteration due to change in the way of
	    // calculation residual
	    flop = iter * 4.0 * param.act_res * param.act_res;

		fprintf(stderr, "Resolution: %5u, ", param.act_res);
		fprintf(stderr, "Time: %04.6f ", runtime);
		fprintf(stderr, "(%3.3f GFlop => %6.2f MFlop/s, ", 
			    flop/1000000000.0,
			    flop/runtime/1000000);
	    fprintf(stderr, "residual %f, %d iterations)\n", residual, iter);

	    // for plot...
	    printf("%5d %f\n", param.act_res, flop/runtime/1000000);
	}

	if (param.act_res + param.res_step_size > param.max_res) 
	{
		break;
	}
	param.act_res += param.res_step_size;

    } // Grid loop ends
   
    sendU( xDim, yDim, rank, np, NumTask, param.u );
    
    if( rank == 0 )
    {
        coarsen( param.u, np, np,
	        param.uvis, param.visres+2, param.visres+2 );
  
        write_image( resfile, param.uvis,  
		    param.visres+2, 
		    param.visres+2 );
    }

    finalize( &param );

    totaltime = wtime() - totaltime;
    printf("\nTotal Program run time = %lf, Rank = %d\n",totaltime, rank );	
    MPI_Finalize();
    return 0;
}

void updateu( double * buf, int size, int row, double * u, int rowStart, int np)
{
	int i = 0;
	for( i=rowStart; i < (size+rowStart); i++ )
	{
		u[row*np + i] = buf[i-rowStart];
	}
}

void updateucol( double * buf, int size, int col, double * u, int colStart, int np)
{
	int i = 0;
	for( i=colStart; i < (size+colStart); i++ )
	{
		u[i*np + col] = buf[i-colStart];
	}
}

void sendResidual(int rank, int numtask, double * residual)
{
	MPI_Status stat;
    int i = 0;

	if(rank==0)
	{
		double buffer;
        for( i=1; i < numtask; i++ )
        {
		    MPI_Recv(&buffer,1, MPI_DOUBLE, i, 10, MPI_COMM_WORLD, &stat);
    		*residual += buffer;
        }
	}
    else
    {
		MPI_Send(residual,1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
    }
}

void sendU(int xDim, int yDim, int rank, int np, int NumTask, double * u )
{
    double * upart = 0;
    int i=0, j=0, k=0;
    int upartSz = (np-2)/NumTask  * (np-2);

    upart = (double *) malloc(sizeof(double) * upartSz);
    int xStart=0, yStart=0, xEnd=0, yEnd= 0;

	MPI_Status stat;

    // Topology 4 * 1
    if( yDim == 1 )
    {
      if( rank != 0 )
      {	     
        xStart =  (np-2)/xDim * rank + 1;
     	xEnd = (np-2)/xDim * (rank+1);
     	yStart = 1;
    	yEnd = np-2;  	

    	// copy
    	for( i=xStart; i<=xEnd; i++ )
    	{
           for( j=yStart; j<=yEnd; j++ )
           {
		        upart[(i-xStart)*(np-2)+(j-yStart)] = u[i*np+j];
           }
        }

	    MPI_Send(upart,upartSz, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD );
      }
     if( rank == 0 )
     {
        for( k=1; k< NumTask; k++ )
        {
            MPI_Recv(upart, upartSz, MPI_DOUBLE, k, 10, MPI_COMM_WORLD, &stat);

	        xStart =  (np-2)/xDim * k + 1; 
         	xEnd = (np-2)/xDim * (k+1);    
         	yStart = 1;
    	    yEnd = np-2;  	
    	
            // copy
    	    for( i=xStart; i<=xEnd; i++ )
    	    {
                for( j=yStart; j<=yEnd; j++ )
                {
        	        u[i*np+j] = upart[(i-xStart)*(np-2)+(j-yStart)];
                }
            }
        }
     }// else closes
   }
   else if( xDim == 1 )// Topology 1 * 4
    {
      if( rank != 0 )
      {	     
        yStart =  (np-2)/yDim * rank + 1;
     	yEnd = (np-2)/yDim * (rank+1);
     	xStart = 1;
    	xEnd = np-2;  	

    	// copy
    	for( i=xStart; i<=xEnd; i++ )
    	{
           for( j=yStart; j<=yEnd; j++ )
           {
		        upart[(i-xStart)*((np-2)/yDim)+(j-yStart)] = u[i*np+j];
           }
        }

	    MPI_Send(upart,upartSz, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD );
     }
     if( rank == 0 )
     {
        for( k=1; k< NumTask; k++ )
        {
            MPI_Recv(upart, upartSz, MPI_DOUBLE, k, 10, MPI_COMM_WORLD, &stat);
	
            yStart =  (np-2)/yDim * k + 1; 
     	    yEnd = (np-2)/yDim * (k+1);    
     	    xStart = 1;
    	    xEnd = np-2;  	
    	
            // copy
    	    for( i=xStart; i<=xEnd; i++ )
    	    {
                for( j=yStart; j<=yEnd; j++ )
                {
            	    u[i*np+j] = upart[(i-xStart)*((np-2)/yDim)+(j-yStart)];
                }
            }
        }
     }// else closes
   }
   else if( yDim == xDim ) // toppology p * p only// Topology 2 * 2
    {
        int stepX = (np-2)/xDim;
        int stepY = (np-2)/yDim;

        if( rank != 0 )
        {
            // Carve out self boundaries
            xStart = stepX * (rank/yDim) + 1 ;
            xEnd = stepX * ( (rank/yDim) + 1);
            yStart = stepY * (rank%xDim) + 1;
            yEnd = stepY * ( (rank%xDim) + 1);

    	    // copy
    	    for( i=xStart; i<=xEnd; i++ )
    	    {
                for( j=yStart; j<=yEnd; j++ )
                {
		            upart[(i-xStart)*((np-2)/yDim)+(j-yStart)] = u[i*np+j];
                }
            }

	        MPI_Send(upart,upartSz, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD );
        }
        else if( rank == 0)
        {
            for( k=1; k< NumTask; k++ )
            {
                MPI_Recv(upart, upartSz, MPI_DOUBLE, k, 10, MPI_COMM_WORLD, &stat);

                xStart = stepX * (k/yDim) + 1 ;
                xEnd = stepX * ( (k/yDim) + 1);
                yStart = stepY * (k%xDim) + 1;
                yEnd = stepY * ( (k%xDim) + 1);
            
                // copy
        	    for( i=xStart; i<=xEnd; i++ )
        	    {
                    for( j=yStart; j<=yEnd; j++ )
                    {
            	        u[i*np+j] = upart[(i-xStart)*((np-2)/yDim)+(j-yStart)];
                    }
                }
            }
        }
   }
    
    free(upart);
}

void commBoundary(int xDim, int yDim, int rank, int np, double * u, MPI_Comm comm)
{

    int coordinates[2];
    int Urank,Lrank,Drank,Rrank;

	MPI_Status stat;
    MPI_Cart_coords( comm, rank, 2, coordinates);
    MPI_Cart_shift( comm, 1, -1, &Rrank, &Lrank);
    MPI_Cart_shift( comm, 0, -1, &Drank, &Urank );

	int i=0, j=0;

    // boundary communication
    if( yDim == 1 ) // config x * 1
	{
	    // Array for MPI sendrecv
	    double * buf; 
	    buf = (double * ) malloc(sizeof(double)*(np-2));
		if( rank == 0)
		{
  		   int index = (np-2)/xDim;

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
    	        buf[i-1] = u[index*np+ i];
   		   }

		   MPI_Send( buf, (np-2), MPI_DOUBLE, Drank, 10, comm);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Drank, 10, comm, &stat);		   
		   updateu( buf, (np-2), index+1, u,1, np);
		}
        else if( rank == xDim -1 )
        {
  		   int index = (np-2)/xDim * rank + 1;

		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Urank, 10, comm, &stat);

		   updateu( buf, (np-2), index-1, u,1, np);
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			    buf[i-1] = u[index*np + i];
   		   }

		   MPI_Send( buf, (np-2), MPI_DOUBLE, Urank, 10, comm);
        }
        else // for all other ranks
        {
            // Communication with upper rank
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Urank, 10, comm, &stat);
  		   int index = (np-2)/xDim * rank + 1;

		   // Update U
		   updateu( buf, (np-2), index-1, u,1, np);

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			    buf[i-1] = u[index*np + i];
   		   }
		   MPI_Send( buf, (np-2), MPI_DOUBLE, Urank, 10, comm);

            // Communication with lower rank
		   index = (np-2)/xDim * (rank+1);

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			    buf[i-1] = u[index*np + i];
   		   }
		   MPI_Send( buf, (np-2), MPI_DOUBLE, Drank, 10, comm);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Drank, 10, comm, &stat);		
		   updateu( buf, (np-2), index+1, u,1, np);
        }
            free(buf);
    }
    else if( xDim == 1 )// for 1 * 4
	{
	    // Array for MPI sendrecv
	    double * buf; 
	    buf = (double * ) malloc(sizeof(double)*(np-2));

		if( rank == 0)
		{
  		   int index = (np-2)/yDim;
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[i*np+ index];
   		   }

		   MPI_Send( buf, (np-2), MPI_DOUBLE, Rrank, 10, comm);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Rrank, 10, comm, &stat);		   
		   updateucol( buf, (np-2), index+1, u,1, np);
		}
        else if( rank == yDim - 1)
        {
  		   int index = (np-2)/yDim * rank + 1;

		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Lrank, 10, comm, &stat);

		   updateucol( buf, (np-2), index-1, u,1, np);
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			    buf[i-1] = u[i*np + index];
   		   }

		   MPI_Send( buf, (np-2), MPI_DOUBLE, Lrank, 10, comm);
        }
        else
        {
           // Communication with left
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Lrank, 10, comm, &stat);
  		   int index = (np-2)/yDim * rank + 1;

		   // Update U
		   updateucol( buf, (np-2), index-1, u,1, np);
 
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			    buf[i-1] = u[i*np + index];
   		   }

		   MPI_Send( buf, (np-2), MPI_DOUBLE, Lrank, 10, comm);

           // Communication with right
		   index = (np-2)/yDim * (rank+1);

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			    buf[i-1] = u[i*np + index];
   		   }
		   MPI_Send( buf, (np-2), MPI_DOUBLE, Rrank, 10, comm);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Rrank, 10, comm, &stat);		
		   updateucol( buf, (np-2), index+1, u,1, np);
        }
            free(buf);
	    }
    else if( xDim == yDim ) // works only for k*k topologies
    {
        double * buf1;
        int size = (np-2)/xDim; 
        buf1 = (double*) malloc(sizeof(double)* size );

        int stepX = size;
        int stepY = size;

        // Carve out self boundaries
        int xStart = stepX * (rank/yDim) + 1 ;
        int xEnd = stepX * ( (rank/yDim) + 1);
        int yStart = stepY * (rank%xDim) + 1;
        int yEnd = stepY * ( (rank%xDim) + 1);

        if( Urank != -1 )
        {
            // Send to Urank
            for( i=yStart; i<=yEnd; i++ )
            {
                // copy row xStart, y index
                buf1[i-yStart] = u[ xStart*np+ i ];
            }        
            MPI_Send( buf1, size, MPI_DOUBLE, Urank, 10, comm);

            // Receive from Urank
            
            // receive from Urank , row (xStart-1)
            MPI_Recv( buf1, size, MPI_DOUBLE, Urank, 10, comm, &stat );
            updateu( buf1, size, (xStart-1), u,yStart, np);
        }
        if( Drank != -1 )
        {
            // Send to Drank
            // copy boundary
            for( i=yStart; i<=yEnd; i++ )
            {
                // buf 1 send to 2
                buf1[i-yStart] = u[ xEnd*np+ i ];
            }        
            MPI_Send( buf1, size, MPI_DOUBLE, Drank, 10, comm);

            // Receive from Drank, row (xEnd+1)
            MPI_Recv( buf1, size, MPI_DOUBLE, Drank, 10, comm, &stat );
            updateu( buf1, size, (xEnd+1), u,yStart,np);            
        }
        if( Rrank != -1 )
        {
            // Send to Rrank
            for( i=xStart; i<=xEnd; i++ )
            {
                buf1[i-xStart] = u[ i*np + yEnd ];
            }
            MPI_Send( buf1, size, MPI_DOUBLE, Rrank, 10, comm);

            // Recieve from Rrank, col (yEnd+1)
            MPI_Recv( buf1, size, MPI_DOUBLE, Rrank, 10, comm, &stat );
		    updateucol( buf1, size, (yEnd+1), u,xStart, np);
        }
        if( Lrank != -1 )
        {
            // Send to Lrank, col yStart
            for( i=xStart; i<=xEnd; i++ )
            {
                // buf2 send to 0
                buf1[i-xStart] = u[ i*np + yStart ];
            }
            MPI_Send( buf1, size, MPI_DOUBLE, Lrank, 10, comm);

            // Receive from Lrank, col (yStart-1)
            MPI_Recv( buf1, size, MPI_DOUBLE, Lrank, 10, comm, &stat );
		    updateucol( buf1, size, (yStart-1), u,xStart, np);
        }

        free(buf1);
    }

}
void updateBoundary( int xDim, int yDim, int rank, int np, double * u, MPI_Comm comm, double * buf1, MPI_Request * request1, double * buf2, MPI_Request * request2, double * buf3, MPI_Request * request3, double * buf4, MPI_Request * request4)
{
    MPI_Status stat;
    if( yDim == 1)
    {
        if( rank == 0)
        {
  		   int index = (np-2)/xDim;
           MPI_Wait( request1, &stat);
		   updateu( buf1, (np-2), index+1, u,1, np);
        }
        else if( rank == xDim - 1)
        {
  		   int index = (np-2)/xDim * rank + 1;
           MPI_Wait( request1, &stat);
		   updateu( buf1, (np-2), index-1, u,1, np);
        }
        else
        {
  		   int index = (np-2)/xDim * rank + 1;
           MPI_Wait( request1, &stat);
		   // Update U
		   updateu( buf1, (np-2), index-1, u,1, np);

		   index = (np-2)/xDim * (rank+1);
           MPI_Wait( request2, &stat);
		   updateu( buf2, (np-2), index+1, u,1, np); 
        }
    }
    else if( xDim == 1)
    {
        if( rank == 0)
        {
  		   int index = (np-2)/yDim;
           MPI_Wait( request1, &stat);
		   updateucol( buf1, (np-2), index+1, u,1, np);
        }
        else if( rank == yDim - 1 )
        {
  		   int index = (np-2)/yDim * rank + 1;
           MPI_Wait( request1, &stat);
		   updateucol( buf1, (np-2), index-1, u,1, np);
        }
        else
        {
  		   int index = (np-2)/yDim * rank + 1;
           MPI_Wait( request1, &stat);
		   // Update U
		   updateucol( buf1, (np-2), index-1, u,1, np);
 
		   index = (np-2)/yDim * (rank+1);
           MPI_Wait( request2, &stat);
		   updateucol( buf2, (np-2), index+1, u,1, np);
        }
    }
    else if( yDim == xDim ) // topology p * p
    {
        int coordinates[2];
        int Urank,Lrank,Drank,Rrank;

        MPI_Cart_coords( comm, rank, 2, coordinates);
        MPI_Cart_shift( comm, 1, -1, &Rrank, &Lrank);
        MPI_Cart_shift( comm, 0, -1, &Drank, &Urank );

        int size = (np-2)/xDim; 
        int stepX = size;
        int stepY = size;

        // Carve out self boundaries
        int xStart = stepX * (rank/yDim) + 1 ;
        int xEnd = stepX * ( (rank/yDim) + 1);
        int yStart = stepY * (rank%xDim) + 1;
        int yEnd = stepY * ( (rank%xDim) + 1);

        if( Urank != -1 )
        {
           MPI_Wait( request1, &stat);
           updateu( buf1, size, (xStart-1), u,yStart, np);
        }
        if( Drank != -1 )
        {
           MPI_Wait( request2, &stat);
           updateu( buf2, size, (xEnd+1), u,yStart,np);            
        }
        if( Rrank != -1 )
        {
           MPI_Wait( request3, &stat);
           updateucol( buf3, size, (yEnd+1), u,xStart, np);
        }
        if( Lrank != -1 )
        {
           MPI_Wait( request4, &stat);
           updateucol( buf4, size, (yStart-1), u,xStart, np);
        }
    }
}

void commBoundaryNB(int xDim, int yDim, int rank, int np, double * u, MPI_Comm comm, double * buf1, MPI_Request * request1, double * buf2, MPI_Request * request2, double * buf3, MPI_Request * request3, double * buf4, MPI_Request * request4)
{
    int coordinates[2];
    int Urank,Lrank,Drank,Rrank;

	MPI_Status stat;
    MPI_Request request;
    MPI_Cart_coords( comm, rank, 2, coordinates);
    MPI_Cart_shift( comm, 1, -1, &Rrank, &Lrank);
    MPI_Cart_shift( comm, 0, -1, &Drank, &Urank );

	int i=0, j=0;

    // boundary communication
    if( yDim == 1 ) // config 4 * 1
	{
	    // Array for MPI sendrecv
        double buf[(np-2)];

		if( rank == 0)
		{
  		    int index = (np-2)/xDim;
		    // Send buffer
		    for(  i=1; i< np-1; i++)
		    {
		    	buf[i-1] = u[index*np+ i];
   		    }
		    MPI_Isend( buf, (np-2), MPI_DOUBLE, Drank, 10, comm, &request);
		    MPI_Irecv( buf1, (np-2), MPI_DOUBLE, Drank, 10, comm, request1);		   
		}
        else if( rank == xDim - 1)
        {
  		    int index = (np-2)/xDim * rank + 1;

		    MPI_Irecv( buf1, (np-2), MPI_DOUBLE, Urank, 10, comm, request1);

		    // Send buffer
		    for(  i=1; i< np-1; i++)
		    {
			    buf[i-1] = u[index*np + i];
   		    }

		    MPI_Isend( buf, (np-2), MPI_DOUBLE, Urank, 10, comm,&request);
        }
        else
        {
            // Communication with upper
	        MPI_Irecv( buf1, (np-2), MPI_DOUBLE, Urank, 10, comm, request1);
  		    int index = (np-2)/xDim * rank + 1;

		    // Send buffer
		    for(  i=1; i< np-1; i++)
		    {
			    buf[i-1] = u[index*np + i];
   		    }

		    MPI_Isend( buf, (np-2), MPI_DOUBLE, Urank, 10, comm, &request);

            // Communication with lower
		    index = (np-2)/xDim * (rank+1);

		    // Send buffer
		    for(  i=1; i< np-1; i++)
		    {
		    	buf[i-1] = u[index*np + i];
   		    }
		    MPI_Isend( buf, (np-2), MPI_DOUBLE, Drank, 10, comm, &request);
		    MPI_Irecv( buf2, (np-2), MPI_DOUBLE, Drank, 10, comm, request2);		
        }
    }
	else if( xDim == 1 )// for 1 * 4
	    {
	        // Array for MPI sendrecv
            double buf[(np-2)];

		    if( rank == 0)
		    {
  		        int index = (np-2)/yDim;
		        // Send buffer
		        for(  i=1; i< np-1; i++)
		        {
			        buf[i-1] = u[i*np+ index];
   		        }

		        MPI_Isend( buf, (np-2), MPI_DOUBLE, Rrank, 10, comm, &request);
		        MPI_Irecv( buf1, (np-2), MPI_DOUBLE, Rrank, 10, comm, request1);		   
		    }
            else if( rank == yDim -1 )
            {
               int index = (np-2)/yDim * rank + 1;

    		   MPI_Irecv( buf1, (np-2), MPI_DOUBLE, Lrank, 10, comm, request1);
	    	   // Send buffer
		       for(  i=1; i< np-1; i++)
		       {
			        buf[i-1] = u[i*np + index];
   		       }

    		   MPI_Isend( buf, (np-2), MPI_DOUBLE, Lrank, 10, comm, &request);
            }
            else
            {
                // Communication with left
		        MPI_Irecv( buf1, (np-2), MPI_DOUBLE, Lrank, 10, comm, request1);
  		        int index = (np-2)/yDim * rank + 1;
 
		        // Send buffer
		        for(  i=1; i< np-1; i++)
		        {
			        buf[i-1] = u[i*np + index];
   		        }

		        MPI_Isend( buf, (np-2), MPI_DOUBLE, Lrank, 10, comm, &request);

                // Communication with right
		        index = (np-2)/yDim * (rank+1);

		        // Send buffer
		        for(  i=1; i< np-1; i++)
		        {
			        buf[i-1] = u[i*np + index];
   		        }
		        MPI_Isend( buf, (np-2), MPI_DOUBLE, Rrank, 10, comm, &request);
		        MPI_Irecv( buf2, (np-2), MPI_DOUBLE, Rrank, 10, comm, request2);
            }
	}
    else if( xDim == yDim  ) // topology p*p
    {
        int size = (np-2)/xDim; 
        int stepX = size;
        int stepY = size;

        double * buf11;
        buf11 = (double*) malloc(sizeof(double)* size );

        // Carve out self boundaries
        int xStart = stepX * (rank/yDim) + 1 ;
        int xEnd = stepX * ( (rank/yDim) + 1);
        int yStart = stepY * (rank%xDim) + 1;
        int yEnd = stepY * ( (rank%xDim) + 1);

        if( Urank != -1 )
        {
            // Send to Urank
            for( i=yStart; i<=yEnd; i++ )
            {
                // copy row xStart, y index
                buf11[i-yStart] = u[ xStart*np+ i ];
            }        
            MPI_Isend( buf11, size, MPI_DOUBLE, Urank, 10, comm, &request);

            // Receive from Urank
            
            // receive from Urank , row (xStart-1)
            MPI_Irecv( buf1, size, MPI_DOUBLE, Urank, 10, comm, request1 );
        }
        if( Drank != -1 )
        {
            // Send to Drank
            // copy boundary
            for( i=yStart; i<=yEnd; i++ )
            {
                // buf 1 send to 2
                buf11[i-yStart] = u[ xEnd*np+ i ];
            }        
            MPI_Isend( buf11, size, MPI_DOUBLE, Drank, 10, comm, &request);

            // Receive from Drank, row (xEnd+1)
            MPI_Irecv( buf2, size, MPI_DOUBLE, Drank, 10, comm, request2 );
        }
        if( Rrank != -1 )
        {
            // Send to Rrank
            for( i=xStart; i<=xEnd; i++ )
            {
                buf11[i-xStart] = u[ i*np + yEnd ];
            }
            MPI_Isend( buf11, size, MPI_DOUBLE, Rrank, 10, comm, &request);

            // Recieve from Rrank, col (yEnd+1)
            MPI_Irecv( buf3, size, MPI_DOUBLE, Rrank, 10, comm, request3 );
        }
        if( Lrank != -1 )
        {
            // Send to Lrank, col yStart
            for( i=xStart; i<=xEnd; i++ )
            {
                // buf2 send to 0
                buf11[i-xStart] = u[ i*np + yStart ];
            }
            MPI_Isend( buf11, size, MPI_DOUBLE, Lrank, 10, comm, &request);

            // Receive from Lrank, col (yStart-1)
            MPI_Irecv( buf4, size, MPI_DOUBLE, Lrank, 10, comm, request4 );
        }
        free(buf11);
    }
}
