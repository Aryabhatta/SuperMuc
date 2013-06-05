/*
 * heat.h
 *
 * Iterative solver for heat distribution
 */


#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "input.h"
#include "heat.h"
#include "timing.h"


void updateu( double * buf, int size, int row, double * u);
void updateucol( double * buf, int size, int row, double * u);
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
    int np, i=0;

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

    MPI_Comm_size( MPI_COMM_WORLD, &NumTask );

    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Status stat;

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

    // Array for MPI sendrecv
    double * buf; 

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

	buf = (double * ) malloc(sizeof(double)*(np-2));
    
	// starting time
	runtime = wtime();
	residual = 999999999;

	iter = 0;
	while(1) {

	switch( param.algorithm ) {

        case 0: // JACOBI

      	//modifed residual calculation
	residual = relax_jacobi(&param.u, &param.uhelp, np, np, \
                   xDim, yDim, rank);

		    break;

	case 1: // GAUSS

		    relax_gauss(param.u, np, np);
		    residual = residual_gauss( param.u, param.uhelp, np, np);
		    break;
	    }
	    // Synchronise
	    MPI_Barrier(MPI_COMM_WORLD);

	    // boundary communication
	    if( yDim == 1 )
	    {
		if( rank == 0)
		{
  		   int index = (np/xDim)-1;
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i] = param.u[index*np+ i];
   		   }

		   MPI_Send( buf, (np-2), MPI_DOUBLE, 1, 10, MPI_COMM_WORLD);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, 1, 10, MPI_COMM_WORLD, &stat);
		   
		   updateu( buf, (np-2), index+1, param.u);
		}
		if( rank == 1 )
		{
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, 0, 10, MPI_COMM_WORLD, &stat);
  		   int index = ((np/xDim)-1 * rank) +1;

		   // Update U
		   updateu( buf, (np-2), index-1, param.u);
 
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i] = param.u[index*np + i];
   		   }

		   MPI_Send( buf, (np-2), MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
		
		   index = ((np/xDim)-1 * (rank+1));

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i] = param.u[index*np + i];
   		   }
		   MPI_Send( buf, (np-2), MPI_DOUBLE, 2, 10, MPI_COMM_WORLD);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, 2, 10, MPI_COMM_WORLD, &stat);
		
		   updateu( buf, (np-2), index+1, param.u);
		}
		if( rank == 2 )
		{
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, 1, 10, MPI_COMM_WORLD, &stat);
  		   int index = ((np/xDim)-1 * rank) +1;

		   // Update U
		   updateu( buf, (np-2), index-1, param.u);
 
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i] = param.u[index*np +i];
   		   }

		   MPI_Send( buf, (np-2), MPI_DOUBLE, 1, 10, MPI_COMM_WORLD);
		
		   index = ((np/xDim)-1 * (rank+1));

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i] = param.u[index*np +i];
   		   }
		   MPI_Send( buf, (np-2), MPI_DOUBLE, 3, 10, MPI_COMM_WORLD);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, 3, 10, MPI_COMM_WORLD, &stat);
		
		   updateu( buf, (np-2), index+1, param.u);
		}
		if( rank == 3 )
		{
  		   int index = (((np/xDim)-1) * rank)+1;
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, 2, 10, MPI_COMM_WORLD, &stat);

		   updateu( buf, (np-2), index-1, param.u);
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i] = param.u[index*np + i];
   		   }

		   MPI_Send( buf, (np-2), MPI_DOUBLE, 2, 10, MPI_COMM_WORLD);
		}
		
	    }
	
	    // for 1 * 4
	    if( xDim == 1 )
	    {
		if( rank == 0)
		{
  		   int index = (np/yDim)-1;

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i] = param.u[i*np+ index];
   		   }

		   MPI_Send( buf, (np-2), MPI_DOUBLE, 1, 10, MPI_COMM_WORLD);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, 1, 10, MPI_COMM_WORLD, &stat);
		   
		   updateucol( buf, (np-2), index+1, param.u);
		}
		if( rank == 1 )
		{
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, 0, 10, MPI_COMM_WORLD, &stat);
  		   int index = ((np/yDim)-1 * rank) +1;

		   // Update U
		   updateucol( buf, (np-2), index-1, param.u);
 
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i] = param.u[i*np + index];
   		   }

		   MPI_Send( buf, (np-2), MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
		
		   index = ((np/yDim)-1 * (rank+1));

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i] = param.u[i*np + index];
   		   }
		   MPI_Send( buf, (np-2), MPI_DOUBLE, 2, 10, MPI_COMM_WORLD);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, 2, 10, MPI_COMM_WORLD, &stat);
		
		   updateucol( buf, (np-2), index+1, param.u);
		}
		if( rank == 2 )
		{
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, 1, 10, MPI_COMM_WORLD, &stat);
  		   int index = ((np/yDim)-1 * rank) +1;

		   // Update U
		   updateucol( buf, (np-2), index-1, param.u);
 
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i] = param.u[i*np +index];
   		   }

		   MPI_Send( buf, (np-2), MPI_DOUBLE, 1, 10, MPI_COMM_WORLD);
		
		   index = ((np/yDim)-1 * (rank+1));

		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i] = param.u[i*np +index];
   		   }
		   MPI_Send( buf, (np-2), MPI_DOUBLE, 3, 10, MPI_COMM_WORLD);
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, 3, 10, MPI_COMM_WORLD, &stat);
		
		   updateucol( buf, (np-2), index+1, param.u);
		}
		if( rank == 3 )
		{
  		   int index = (((np/yDim)-1) * rank)+1;
		   MPI_Recv( buf, (np-2), MPI_DOUBLE, 2, 10, MPI_COMM_WORLD, &stat);

		   updateucol( buf, (np-2), index-1, param.u);
		   // Send buffer
		   for(  i=1; i< np-1; i++)
		   {
			buf[i] = param.u[i*np + index];
   		   }

		   MPI_Send( buf, (np-2), MPI_DOUBLE, 2, 10, MPI_COMM_WORLD);
		}		
	    }

	    iter++;

	    // solution good enough ?
//	    if (residual < 0.000005){printf("\nResidual less, rank %d exits", rank);} break;

	    // max. iteration reached ? (no limit with maxiter=0)
	    if (param.maxiter>0 && iter>=param.maxiter) break;
	    
	    if (iter % 100 == 0)
		fprintf(stderr, "residual %f, %d iterations\n", residual, iter);
	}

	// Flop count after <i> iterations
	//flop = iter * 11.0 * param.act_res * param.act_res;

	// change in no of flops/iteration due to change in the way of
	// calculation residual
	flop = iter * 4.0 * param.act_res/xDim * param.act_res/yDim;

	// stopping time
	runtime = wtime() - runtime;
	if(rank==0)
	{
		double buffer;
		MPI_Recv(&buffer,1, MPI_DOUBLE,1, 10, MPI_COMM_WORLD, &stat);
		residual += buffer;
		MPI_Recv(&buffer,1, MPI_DOUBLE,2, 10, MPI_COMM_WORLD, &stat);
		residual += buffer;
		MPI_Recv(&buffer,1, MPI_DOUBLE,3, 10, MPI_COMM_WORLD, &stat);
		residual += buffer;
	}
	if(rank==1 || rank==2 || rank ==3 )
	{
		MPI_Send(&residual,1, MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
	}

	if(rank==0)
	{
		fprintf(stderr, "Resolution: %5u, ", param.act_res);
		fprintf(stderr, "Time: %04.3f ", runtime);
		fprintf(stderr, "(%3.3f GFlop => %6.2f MFlop/s, ", 
			flop/1000000000.0,
			flop/runtime/1000000);
	fprintf(stderr, "residual %f, %d iterations)\n", residual, iter);

	// for plot...
	printf("%5d %f\n", param.act_res, flop/runtime/1000000);
	}
//	printf("Residual=%lf  ", residual);

	if (param.act_res + param.res_step_size > param.max_res) 
	{
		printf("\n Rank %d exits...", rank);
		break;
	}
	param.act_res += param.res_step_size;
    }

    MPI_Finalize();

    coarsen( param.u, np, np,
	     param.uvis, param.visres+2, param.visres+2 );
  
    write_image( resfile, param.uvis,  
		 param.visres+2, 
		 param.visres+2 );

    finalize( &param );

    totaltime = wtime() - totaltime;
    printf("\nTotal Program run time = %lf, Rank = %d\n",totaltime, rank );	
    free(buf);
    return 0;
}

void updateu( double * buf, int size, int row, double * u)
{
	int i = 0;
	int np = size + 2;
	for( i=1; i<= size; i++ )
	{
		u[row*np + i] = buf[i-1];
	}
}

void updateucol( double * buf, int size, int col, double * u)
{
	int i = 0;
	int np = size + 2;
	for( i=1; i<= size; i++ )
	{
		u[i*np + col] = buf[i-1];
	}
}
