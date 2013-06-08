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

void commBoundaryNB(int xDim, int yDim, int rank, int np, double * u, MPI_Comm comm, double * buf1, MPI_Request * request1, double * buf2, MPI_Request * request2);
void updateBoundary( int xDim, int yDim, int rank, int np, double * u, double * buf1, MPI_Request * request1, double * buf2, MPI_Request * request2);

void sendU(int xDim, int yDim, int rank, int np, int NumTask, double * u );
void sendResidual(int rank, double * residual);

void usage( char *s )
{
    fprintf(stderr, 
	    "Usage: %s <input file> [X-dim] [Y-dim] [result file] \n\n", s);
}

void foo( MPI_Comm comm )
{
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

    // check argumentsi
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

    // MPI initialise
    int NumTask = 0, rank = 0;
    int xDim =0, yDim = 0;

    xDim = atoi(argv[2]);
    yDim = atoi(argv[3]);

    int status = MPI_Init( &argc, &argv);

    if( status != MPI_SUCCESS )
    {
	printf("\nError in MPI\n");
	return 0;
    }
    
    MPI_Status stat;
    MPI_Comm comm;
    int dims[2];// = {2,2};
    dims[0] = xDim;
    dims[1] = yDim;
    int periods[2]={0, 0};
    bool reorder = false;

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

#if BLOCK == 0
        MPI_Request request1, request2;
        double * buf1, * buf2; 
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

	fprintf(stderr, "Resolution: %5u\r", param.act_res);

	// full size (param.act_res are only the inner points)
	np = param.act_res + 2;

#if BLOCK == 0
	    buf1 = (double * ) malloc(sizeof(double)*(np-2));
    	buf2 = (double * ) malloc(sizeof(double)*(np-2));
#endif
    
	// starting time
	runtime = wtime();
	residual = 999999999;

	iter = 0;
	while(1) {
#if BLOCK
	    commBoundary( xDim, yDim, rank, np, param.u, comm);
#endif

#if BLOCK == 0
//	    commBoundaryNB( xDim, yDim, rank, np, param.u, comm);
        commBoundaryNB( xDim, yDim, rank, np, param.u, comm, buf1, &request1, buf2, &request2);
#endif

	switch( param.algorithm ) {

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

		    relax_gauss(param.u, np, np);
		    residual = residual_gauss( param.u, param.uhelp, np, np);
		    break;
	    }
#if BLOCK
	    // Synchronise
	    MPI_Barrier(MPI_COMM_WORLD);
#endif


	    iter++;

	    // solution good enough ?
//	    if (residual < 0.000005){printf("\nResidual less, rank %d exits", rank);} break;

	    // max. iteration reached ? (no limit with maxiter=0)
	    if (param.maxiter>0 && iter>=param.maxiter) break;
	    
	  //  if (iter % 100 == 0)
	//	fprinitf(stderr, "residual %f, %d iterations\n", residual, iter);

#if BLOCK == 0
    //check for boundary
    updateBoundary( xDim, yDim, rank, np, param.u, buf1, &request1, buf2, &request2);
    // relax boundary
	residual += relax_jacobiBoundary(&param.u, &param.uhelp, np, np, \
                   xDim, yDim, rank);
#endif
	}

	// stopping time
	runtime = wtime() - runtime;

	sendResidual(rank, &residual);

	if( rank == 0 )
	{
	// change in no of flops/iteration due to change in the way of
	// calculation residual
	flop = iter * 4.0 * param.act_res * param.act_res;

		fprintf(stderr, "Resolution: %5u, ", param.act_res);
		fprintf(stderr, "Time: %04.3f ", runtime);
		fprintf(stderr, "(%3.3f GFlop => %6.2f MFlop/s, ", 
			flop/1000000000.0,
			flop/runtime/1000000);
	fprintf(stderr, "residual %f, %d iterations)\n", residual, iter);

	// for plot...
	printf("%5d %f\n", param.act_res, flop/runtime/1000000);
	}

	if (param.act_res + param.res_step_size > param.max_res) 
	{
	//	printf("\n Rank %d exits...@resolution=%d", rank, param.act_res);
		break;
	}
	param.act_res += param.res_step_size;
    }

#if BLOCK == 0
    free(buf1);
    free(buf2);
#endif    

    sendU( xDim, yDim, rank, np, NumTask, param.u );

    MPI_Finalize();
    
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

void sendResidual(int rank, double * residual)
{
	MPI_Status stat;
	if(rank==0)
	{
		double buffer;
		MPI_Recv(&buffer,1, MPI_DOUBLE,1, 10, MPI_COMM_WORLD, &stat);
		*residual += buffer;
		MPI_Recv(&buffer,1, MPI_DOUBLE,2, 10, MPI_COMM_WORLD, &stat);
		*residual += buffer;
		MPI_Recv(&buffer,1, MPI_DOUBLE,3, 10, MPI_COMM_WORLD, &stat);
		*residual += buffer;
	}
	if(rank==1 || rank==2 || rank ==3 )
	{
		MPI_Send(residual,1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
	}
}

void sendU(int xDim, int yDim, int rank, int np, int NumTask, double * u )
{
    double * upart = 0;
    int i=0, j=0;
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

       // Receive from 1
       MPI_Recv(upart, upartSz, MPI_DOUBLE, 1, 10, MPI_COMM_WORLD, &stat);

	xStart =  (np-2)/xDim * 1 + 1; 
     	xEnd = (np-2)/xDim * (1+1);    
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

       // Receive from 2
       MPI_Recv(upart, upartSz, MPI_DOUBLE, 2, 10, MPI_COMM_WORLD, &stat);

	xStart =  (np-2)/xDim * 2 + 1;  
     	xEnd = (np-2)/xDim * (2+1);	  // 75

    	// copy
    	for( i=xStart; i<=xEnd; i++ )
    	{
           for( j=yStart; j<=yEnd; j++ )
           {
        	u[i*np+j] = upart[(i-xStart)*(np-2)+(j-yStart)];
           }
        }

       // Receive from 3
       MPI_Recv(upart, upartSz, MPI_DOUBLE, 3, 10, MPI_COMM_WORLD, &stat);

	xStart =  (np-2)/xDim * 3 + 1; // 76
     	xEnd = (np-2)/xDim * (3+1);	 // 100	

    	// copy
    	for( i=xStart; i<=xEnd; i++ )
    	{
           for( j=yStart; j<=yEnd; j++ )
           {
        	u[i*np+j] = upart[(i-xStart)*(np-2)+(j-yStart)];
           }
        }
     }// else closes
   }
    // Topology 1 * 4
    if( xDim == 1 )
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

       // Receive from 1
       MPI_Recv(upart, upartSz, MPI_DOUBLE, 1, 10, MPI_COMM_WORLD, &stat);

	yStart =  (np-2)/yDim * 1 + 1; 
     	yEnd = (np-2)/yDim * (1+1);    
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

       // Receive from 2
       MPI_Recv(upart, upartSz, MPI_DOUBLE, 2, 10, MPI_COMM_WORLD, &stat);

	yStart =  (np-2)/yDim * 2 + 1;  
     	yEnd = (np-2)/yDim * (2+1);	  // 75

    	// copy
    	for( i=xStart; i<=xEnd; i++ )
    	{
           for( j=yStart; j<=yEnd; j++ )
           {
        	u[i*np+j] = upart[(i-xStart)*((np-2)/yDim)+(j-yStart)];
           }
        }

       // Receive from 3
       MPI_Recv(upart, upartSz, MPI_DOUBLE, 3, 10, MPI_COMM_WORLD, &stat);

	    yStart =  (np-2)/yDim * 3 + 1; // 76
     	yEnd = (np-2)/yDim * (3+1);	 // 100	

    	// copy
    	for( i=xStart; i<=xEnd; i++ )
    	{
           for( j=yStart; j<=yEnd; j++ )
           {
        	u[i*np+j] = upart[(i-xStart)*((np-2)/yDim)+(j-yStart)];
           }
        }
     }// else closes
   }

    // Topology 2 * 2
    if( yDim == 2 && xDim == 2 )
    {
      if( rank == 1 )
      {	     

        xStart = (np-2)/xDim * (rank-1) + 1; // 1
        xEnd = (np-2)/xDim * rank;          // 50
        yStart = (np-2)/yDim * rank + 1;    // 51
        yEnd = (np-2)/yDim * (rank+1);      // 100

    	// copy
    	for( i=xStart; i<=xEnd; i++ )
    	{
           for( j=yStart; j<=yEnd; j++ )
           {
    		upart[(i-xStart)*((np-2)/2)+(j-yStart)] = u[i*np+j];
           }
        }

	    MPI_Send(upart,upartSz, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD );
     }
     if( rank == 3 )
     {
    
        xStart = (np-2)/xDim * (rank-2) + 1;    // 51
        xEnd = (np-2)/xDim * (rank-1);            // 100
        yStart = (np-2)/yDim * (rank-2) + 1;    // 51
        yEnd = (np-2)/yDim * (rank-1);            // 100

    	// copy
    	for( i=xStart; i<=xEnd; i++ )
    	{
           for( j=yStart; j<=yEnd; j++ )
           {
    		upart[(i-xStart)*((np-2)/2)+(j-yStart)] = u[i*np+j];
           }
        }

	    MPI_Send(upart,upartSz, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD );
     }
     if( rank == 2 )
     {
        xStart = (np-2)/xDim * (rank-1) + 1;    // 51
        xEnd = (np-2)/xDim * (rank);          // 100
        yStart = (np-2)/yDim * (rank-2) + 1;    // 1
        yEnd = (np-2)/yDim * (rank-1);          // 50

    	// copy
    	for( i=xStart; i<=xEnd; i++ )
    	{
           for( j=yStart; j<=yEnd; j++ )
           {
    		upart[(i-xStart)*((np-2)/2)+(j-yStart)] = u[i*np+j];
           }
        }

	    MPI_Send(upart,upartSz, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD );
     }
     if( rank == 0 )
     {

       // Receive from 1
       MPI_Recv(upart, upartSz, MPI_DOUBLE, 1, 10, MPI_COMM_WORLD, &stat);

        xStart = (np-2)/xDim * (1-1) + 1;
        xEnd = (np-2)/xDim * 1;
        yStart = (np-2)/yDim * 1 + 1;
        yEnd = (np-2)/yDim * (1+1);

    	// copy
    	for( i=xStart; i<=xEnd; i++ )
    	{
           for( j=yStart; j<=yEnd; j++ )
           {
    		 u[i*np+j] = upart[(i-xStart)*((np-2)/2)+(j-yStart)] ;
           }
        }

       // Receive from 3
       MPI_Recv(upart, upartSz, MPI_DOUBLE, 3, 10, MPI_COMM_WORLD, &stat);

        xStart = (np-2)/xDim * (3-2) + 1;
        xEnd = (np-2)/xDim * (3-1);
        yStart = (np-2)/yDim * (3-2) + 1;
        yEnd = (np-2)/yDim * (3-1);

    	// copy
    	for( i=xStart; i<=xEnd; i++ )
    	{
           for( j=yStart; j<=yEnd; j++ )
           {
    	      u[i*np+j] = upart[(i-xStart)*((np-2)/2)+(j-yStart)];
           }
        }


       // Receive from 2
       MPI_Recv(upart, upartSz, MPI_DOUBLE, 2, 10, MPI_COMM_WORLD, &stat);

        xStart = (np-2)/xDim * (2-1) + 1;
        xEnd = (np-2)/xDim * (2);
        yStart = (np-2)/yDim * (2-2) + 1;
        yEnd = (np-2)/yDim * (2-1);

    	// copy
    	for( i=xStart; i<=xEnd; i++ )
    	{
           for( j=yStart; j<=yEnd; j++ )
           {
    	      u[i*np+j] = upart[(i-xStart)*((np-2)/2)+(j-yStart)];
           }
        }

     }// else closes
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
//    printf("\nRank = %d, Rrank = %d, LRank= %d, Urank=%d, Drank = %d ", rank, Rrank, Lrank, Urank, Drank);

	int i=0, j=0;

    // boundary communication
    if( yDim == 1 ) // config 4 * 1
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
//		   printf("\nRank %d sends row %d", rank, index);
		   MPI_Send( buf, (np-2), MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD, &stat);		   
//		   printf("\nRank %d receives row %d", rank, index+1);
		   updateu( buf, (np-2), index+1, u,1, np);
		}
		if( rank == 1 )
		{
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD, &stat);
  		   int index = (np-2)/xDim * rank + 1;

//		   printf("\nRank %d receives row %d", rank, index-1);
		   // Update U
		   updateu( buf, (np-2), index-1, u,1, np);
 
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[index*np + i];
   		   }

//		   printf("\nRank %d sends row %d", rank, index);
		   MPI_Send( buf, (np-2), MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD);
		
		   index = (np-2)/xDim * (rank+1);

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[index*np + i];
   		   }
//		   printf("\nRank %d sends row %d", rank, index);
		   MPI_Send( buf, (np-2), MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD, &stat);		
//		   printf("\nRank %d receives row %d", rank, index+1);
		   updateu( buf, (np-2), index+1, u,1, np);
		}
		if( rank == 2 )
		{
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD, &stat);

  		   int index = (np-2)/xDim * rank + 1;

//		   printf("\nRank %d receives row %d", rank, index-1);
		   // Update U
		   updateu( buf, (np-2), index-1, u,1, np);
 
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[index*np +i];
   		   }

//		   printf("\nRank %d sends row %d", rank, index);
		   MPI_Send( buf, (np-2), MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD);
		
		   index =  (np-2)/xDim * (rank+1);

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[index*np +i];
   		   }
//		   printf("\nRank %d sends row %d", rank, index);
		   MPI_Send( buf, (np-2), MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD, &stat);
		
//		   printf("\nRank %d receives row %d", rank, index+1);
		   updateu( buf, (np-2), index+1, u,1, np);
		}
		if( rank == 3 )
		{
  		   int index = (np-2)/xDim * rank + 1;

		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD, &stat);

//		   printf("\nRank %d receives row %d", rank, index-1);
		   updateu( buf, (np-2), index-1, u,1, np);
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[index*np + i];
   		   }

//		   printf("\nRank %d sends row %d", rank, index);
		   MPI_Send( buf, (np-2), MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD);
		}
		
    free(buf);
	    }
	
	    // for 1 * 4
	    if( xDim == 1 )
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

//		   printf("\nRank %d sends row %d", rank, index);
		   MPI_Send( buf, (np-2), MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD, &stat);		   
//		   printf("\nRank %d receives row %d", rank, index+1);
		   updateucol( buf, (np-2), index+1, u,1, np);
		}
		if( rank == 1 )
		{
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD, &stat);
  		   int index = (np-2)/yDim * rank + 1;

//		   printf("\nRank %d receives row %d", rank, index-1);
		   // Update U
		   updateucol( buf, (np-2), index-1, u,1, np);
 
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[i*np + index];
   		   }

//		   printf("\nRank %d sends row %d", rank, index);
		   MPI_Send( buf, (np-2), MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD);
		
		   index = (np-2)/yDim * (rank+1);

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[i*np + index];
   		   }
//		   printf("\nRank %d sends row %d", rank, index);
		   MPI_Send( buf, (np-2), MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD, &stat);		
//		   printf("\nRank %d receives row %d", rank, index+1);
		   updateucol( buf, (np-2), index+1, u,1, np);
		}
		if( rank == 2 )
		{
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD, &stat);

  		   int index = (np-2)/yDim * rank + 1;

//		   printf("\nRank %d receives row %d", rank, index-1);
		   // Update U
		   updateucol( buf, (np-2), index-1, u,1, np);
 
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[i*np +index];
   		   }

//		   printf("\nRank %d sends row %d", rank, index);
		   MPI_Send( buf, (np-2), MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD);
		
		   index =  (np-2)/yDim * (rank+1);

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[i*np +index];
   		   }
//		   printf("\nRank %d sends row %d", rank, index);
		   MPI_Send( buf, (np-2), MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD, &stat);
		
//		   printf("\nRank %d receives row %d", rank, index+1);
		   updateucol( buf, (np-2), index+1, u,1, np);
		}
		if( rank == 3 )
		{
  		   int index = (np-2)/yDim * rank + 1;

		   MPI_Recv( buf, (np-2), MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD, &stat);

//		   printf("\nRank %d receives row %d", rank, index-1);
		   updateucol( buf, (np-2), index-1, u,1, np);
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[i*np + index];
   		   }

//		   printf("\nRank %d sends row %d", rank, index);
		   MPI_Send( buf, (np-2), MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD);
		}
		
    free(buf);
	    }

    if( xDim == 2 && yDim == 2 )
    {
        double * buf1, * buf2;
        int size = (np-2)/2;
        buf1 = (double*) malloc(sizeof(double)* size );
        buf2 = (double*) malloc(sizeof(double)* size );

        if( rank == 0 )
        {
            int xIndex = size;
            int yIndex = size;

            for( i=0; i<size; i++ )
            {
                // buf 1 send to 2
                buf1[i] = u[ xIndex*np+ (i+1) ];
            
                // buf2 send to 1
                buf2[i] = u[ (i+1)*np + yIndex ];
            }
        
            MPI_Send( buf1, size, MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD);
            MPI_Send( buf2, size, MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD);

            // receive from 2 , row 51
            MPI_Recv( buf1, size, MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD, &stat );
            updateu( buf1, size, size+1, u,1,np);

            // receive from 1 , col 51
            MPI_Recv( buf2, size, MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD, &stat );
		    updateucol( buf2, size, size+1, u,1, np);
        }
        if( rank == 1 )
        {            
            // receive from 3 , row 51
            MPI_Recv( buf1, size, MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD, &stat );
            updateu( buf1, size, size+1, u,(size+1), np);

            // receive from 0 , col 50
            MPI_Recv( buf2, size, MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD, &stat );
		    updateucol( buf2, size, size, u,1, np);

            int xIndex = size;
            int yIndex = size + 1;

            for( i=size; i<(size+size); i++ )
            {
                // buf 1 send to 3
                buf1[i-size] = u[ xIndex*np+ (i+1) ];
            }
            
            for( i=0; i<size; i++ )
            {
                // buf2 send to 0
                buf2[i] = u[ (i+1)*np + yIndex ];
            }
        
            MPI_Send( buf1, size, MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD);
            MPI_Send( buf2, size, MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD);
        }
        if( rank == 3 )
        {
            int xIndex = size + 1;
            int yIndex = size + 1;

            for( i=size; i<(size+size); i++ )
            {
                // buf 1 send to 1
                buf1[i-size] = u[ xIndex*np+ (i+1) ];
            }

            MPI_Send( buf1, size, MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD);

            // receive from 2 , col 50
            MPI_Recv( buf2, size, MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD, &stat );
		    updateucol( buf2, size, size, u,(size+1), np);

            // receive from 1 , row 50
            MPI_Recv( buf1, size, MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD, &stat );
            updateu( buf1, size, size, u,(size+1), np);

            for( i=size; i<(size+size); i++ )
            {
                // buf2 send to 2
                buf2[i-size] = u[ (i+1)*np + yIndex ];
            }
        
            MPI_Send( buf2, size, MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD);
        }
        if( rank == 2 )
        {
            int xIndex = size + 1;
            int yIndex = size;

            // receive from 0 , row 50
            MPI_Recv( buf1, size, MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD, &stat );
            updateu( buf1, size, size, u,1, np);

            for( i=size; i<(size+size); i++ )
            {
                // buf2 send to 3
                buf2[i-size] = u[ (i+1)*np + yIndex ];
            }

            MPI_Send( buf2, size, MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD);

            for( i=0; i<size; i++ )
            {
                // buf 1 send to 0
                buf1[i] = u[ xIndex*np+ (i+1) ];
            }
        
            MPI_Send( buf1, size, MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD);

            // receive from 3 , col 51
            MPI_Recv( buf2, size, MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD, &stat );
		    updateucol( buf2, size, size+1, u, (size+1), np);
        }

        free(buf1);
        free(buf2);
    }

}
void updateBoundary( int xDim, int yDim, int rank, int np, double * u, double * buf1, MPI_Request * request1, double * buf2, MPI_Request * request2)
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
        if( rank == 1)
        {
  		   int index = (np-2)/xDim * rank + 1;
           MPI_Wait( request1, &stat);
		   // Update U
		   updateu( buf1, (np-2), index-1, u,1, np);

		   index = (np-2)/xDim * (rank+1);
           MPI_Wait( request2, &stat);
		   updateu( buf2, (np-2), index+1, u,1, np); 
        }
        if( rank == 2)
        {
  		   int index = (np-2)/xDim * rank + 1;
           MPI_Wait( request1, &stat);
		   // Update U
		   updateu( buf1, (np-2), index-1, u,1, np);

		   index =  (np-2)/xDim * (rank+1);
           MPI_Wait( request2, &stat);
		   updateu( buf2, (np-2), index+1, u,1, np); 
        }
        if( rank == 3)
        {
  		   int index = (np-2)/xDim * rank + 1;
           MPI_Wait( request1, &stat);
		   updateu( buf1, (np-2), index-1, u,1, np);
        }
        
    }
    if( xDim == 1)
    {
        if( rank == 0)
        {
  		   int index = (np-2)/yDim;
           MPI_Wait( request1, &stat);
		   updateucol( buf1, (np-2), index+1, u,1, np);
        }
        if( rank == 1)
        {
  		   int index = (np-2)/yDim * rank + 1;
           MPI_Wait( request1, &stat);
		   // Update U
		   updateucol( buf1, (np-2), index-1, u,1, np);
 
		   index = (np-2)/yDim * (rank+1);
           MPI_Wait( request2, &stat);
		   updateucol( buf2, (np-2), index+1, u,1, np);
        }
        if( rank == 2)
        {
  		   int index = (np-2)/yDim * rank + 1;
           MPI_Wait( request1, &stat);
    	   // Update U
		   updateucol( buf1, (np-2), index-1, u,1, np);
 
		   index =  (np-2)/yDim * (rank+1);
           MPI_Wait( request2, &stat);
		   updateucol( buf2, (np-2), index+1, u,1, np);
        }
        if( rank == 3)
        {
  		   int index = (np-2)/yDim * rank + 1;
           MPI_Wait( request1, &stat);
		   updateucol( buf1, (np-2), index-1, u,1, np);
        }
    }
    if( yDim ==2 && xDim == 2)
    {
        int size = (np-2)/2;
        if( rank == 0)
        {
           MPI_Wait( request1, &stat);
           updateu( buf1, size, size+1, u,1,np);

           MPI_Wait( request2, &stat);
		   updateucol( buf2, size, size+1, u,1, np);
        }
        if( rank == 1)
        {
           MPI_Wait( request1, &stat);
           updateu( buf1, size, size+1, u,(size+1), np);

           MPI_Wait( request2, &stat);
 	       updateucol( buf2, size, size, u,1, np);
        }
        if( rank == 3)
        {
           MPI_Wait( request2, &stat);
		   updateucol( buf2, size, size, u,(size+1), np);

           MPI_Wait( request1, &stat);
           updateu( buf1, size, size, u,(size+1), np);
        }
        if( rank == 2)
        {
           MPI_Wait( request1, &stat);
           updateu( buf1, size, size, u,1, np);

           MPI_Wait( request2, &stat);
           updateucol( buf2, size, size+1, u, (size+1), np);
        }
    }
}

void commBoundaryNB(int xDim, int yDim, int rank, int np, double * u, MPI_Comm comm, double * buf1, MPI_Request * request1, double * buf2, MPI_Request * request2)
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
		   MPI_Isend( buf, (np-2), MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD, &request);
		   MPI_Irecv( buf1, (np-2), MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD, request1);		   
		}
		if( rank == 1 )
		{
		   MPI_Irecv( buf1, (np-2), MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD, request1);
  		   int index = (np-2)/xDim * rank + 1;

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[index*np + i];
   		   }

		   MPI_Isend( buf, (np-2), MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD, &request);
		
		   index = (np-2)/xDim * (rank+1);

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[index*np + i];
   		   }
		   MPI_Isend( buf, (np-2), MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD, &request);
		   MPI_Irecv( buf2, (np-2), MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD, request2);		
		}
		if( rank == 2 )
		{
		   MPI_Irecv( buf1, (np-2), MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD, request1);

  		   int index = (np-2)/xDim * rank + 1;

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[index*np +i];
   		   }

		   MPI_Isend( buf, (np-2), MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD, &request);
		
		   index =  (np-2)/xDim * (rank+1);

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[index*np +i];
   		   }
		   MPI_Isend( buf, (np-2), MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD, &request);
		   MPI_Irecv( buf2, (np-2), MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD, request2);
		
		}
		if( rank == 3 )
		{
  		   int index = (np-2)/xDim * rank + 1;

		   MPI_Irecv( buf1, (np-2), MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD, request1);

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[index*np + i];
   		   }

		   MPI_Isend( buf, (np-2), MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD,&request);
		}
		
    free(buf);
	    }
	
	    // for 1 * 4
	    if( xDim == 1 )
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

		   MPI_Isend( buf, (np-2), MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD, &request);
		   MPI_Irecv( buf1, (np-2), MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD, request1);		   
		}
		if( rank == 1 )
		{
		   MPI_Irecv( buf1, (np-2), MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD, request1);
  		   int index = (np-2)/yDim * rank + 1;
 
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[i*np + index];
   		   }

		   MPI_Isend( buf, (np-2), MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD, &request);
		
		   index = (np-2)/yDim * (rank+1);

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[i*np + index];
   		   }
		   MPI_Isend( buf, (np-2), MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD, &request);
		   MPI_Irecv( buf2, (np-2), MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD, request2);		
		}
		if( rank == 2 )
		{
		   MPI_Irecv( buf1, (np-2), MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD, request1);

  		   int index = (np-2)/yDim * rank + 1;

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[i*np +index];
   		   }

		   MPI_Isend( buf, (np-2), MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD, &request);
		
		   index =  (np-2)/yDim * (rank+1);

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[i*np +index];
   		   }
		   MPI_Isend( buf, (np-2), MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD, &request);
		   MPI_Irecv( buf2, (np-2), MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD, request2);
		}
		if( rank == 3 )
		{
  		   int index = (np-2)/yDim * rank + 1;

		   MPI_Irecv( buf1, (np-2), MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD, request1);
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i-1] = u[i*np + index];
   		   }

		   MPI_Isend( buf, (np-2), MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD, &request);
		}
		
    free(buf);
	    }

    if( xDim == 2 && yDim == 2 )
    {
        double * buf11, * buf22;
        int size = (np-2)/2;
        buf11 = (double*) malloc(sizeof(double)* size );
        buf22 = (double*) malloc(sizeof(double)* size );

        if( rank == 0 )
        {
            int xIndex = size;
            int yIndex = size;

            for( i=0; i<size; i++ )
            {
                // buf 1 send to 2
                buf11[i] = u[ xIndex*np+ (i+1) ];
            
                // buf2 send to 1
                buf22[i] = u[ (i+1)*np + yIndex ];
            }
        
            MPI_Isend( buf11, size, MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD, &request);
            MPI_Isend( buf22, size, MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD, &request);

            // receive from 2 , row 51
            MPI_Irecv( buf1, size, MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD, request1 );

            // receive from 1 , col 51
            MPI_Irecv( buf2, size, MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD, request2 );
        }
        if( rank == 1 )
        {            
            // receive from 3 , row 51
            MPI_Irecv( buf1, size, MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD, request1 );

            // receive from 0 , col 50
            MPI_Irecv( buf2, size, MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD, request2 );

            int xIndex = size;
            int yIndex = size + 1;

            for( i=size; i<(size+size); i++ )
            {
                // buf 1 send to 3
                buf11[i-size] = u[ xIndex*np+ (i+1) ];
            }
            
            for( i=0; i<size; i++ )
            {
                // buf2 send to 0
                buf22[i] = u[ (i+1)*np + yIndex ];
            }
        
            MPI_Isend( buf11, size, MPI_DOUBLE, Drank, 10, MPI_COMM_WORLD, &request);
            MPI_Isend( buf22, size, MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD, &request);
        }
        if( rank == 3 )
        {
            int xIndex = size + 1;
            int yIndex = size + 1;

            for( i=size; i<(size+size); i++ )
            {
                // buf 1 send to 1
                buf11[i-size] = u[ xIndex*np+ (i+1) ];
            }

            MPI_Isend( buf11, size, MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD, &request);

            // receive from 2 , col 50
            MPI_Irecv( buf2, size, MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD, request2 );

            // receive from 1 , row 50
            MPI_Irecv( buf1, size, MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD, request1 );

            for( i=size; i<(size+size); i++ )
            {
                // buf2 send to 2
                buf22[i-size] = u[ (i+1)*np + yIndex ];
            }
        
            MPI_Isend( buf22, size, MPI_DOUBLE, Lrank, 10, MPI_COMM_WORLD, &request);
        }
        if( rank == 2 )
        {
            int xIndex = size + 1;
            int yIndex = size;

            // receive from 0 , row 50
            MPI_Irecv( buf1, size, MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD, request1 );

            for( i=size; i<(size+size); i++ )
            {
                // buf2 send to 3
                buf22[i-size] = u[ (i+1)*np + yIndex ];
            }

            MPI_Isend( buf22, size, MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD, &request);

            for( i=0; i<size; i++ )
            {
                // buf 1 send to 0
                buf11[i] = u[ xIndex*np+ (i+1) ];
            }
        
            MPI_Isend( buf11, size, MPI_DOUBLE, Urank, 10, MPI_COMM_WORLD, &request);

            // receive from 3 , col 51
            MPI_Irecv( buf2, size, MPI_DOUBLE, Rrank, 10, MPI_COMM_WORLD, request2 );
        }

        free(buf11);
        free(buf22);
    }

}
