/*
 * heat.h
 *
 * Iterative solver for heat distribution
 */


#include <stdio.h>
#include <stdlib.h>
#include "input.h"
#include "heat.h"
#include "timing.h"
#include <omp.h>

void usage( char *s )
{
    fprintf(stderr, 
	    "Usage: %s <input file> [result file]\n\n", s);
}


int main( int argc, char *argv[] )
{
    unsigned iter;
    FILE *infile, *resfile;
    char *resfilename;
    int i =0 , NoIter = 0;

    // algorithmic parameters
    algoparam_t param;
    int np;

    double runtime, flop;
    double residual;
    double totaltime=0;

    totaltime = wtime();

    // check argumentsi
    if( argc < 2 )
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
    resfilename= (argc>=3) ? argv[2]:"heat.ppm";

    if( !(resfile=fopen(resfilename, "w")) )
    {
	fprintf(stderr, 
		"\nError: Cannot open \"%s\" for writing.\n\n", 
		resfilename);
      
	usage(argv[0]);
	return 1;
    }

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
   // printf("the number of threads %d",omp_get_num_threads());

     while(1)
     {
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
    
	// starting time
	runtime = wtime();
	residual = 999999999;

        iter = 0;
	while(1) {

	    switch( param.algorithm ) {

		case 0: // JACOBI

      		    //modifed residual calculation
		    residual = relax_jacobi(param.u, param.uhelp, np, np);
		   // residual = residual_jacobi( param.u, param.uhelp, np, np);
		    break;

		case 1: // GAUSS

		    relax_gauss(param.u, np, np);
		    residual = residual_gauss( param.u, param.uhelp, np, np);
		    break;
	    }
            iter ++;

	    // solution good enough ?
	    if (residual < 0.000005) break;

            if (param.maxiter>0 && iter>=param.maxiter) break;

            if (iter % 100 ==0)
		fprintf(stderr,"residual %f, %d iterations\n", residual, iter);
	}
//       }

	// Flop count after <i> iterations
	//flop = iter * 11.0 * param.act_res * param.act_res;

	// change in no of flops/iteration due to change in the way of
	// calculation residual & 4 iterations in single relax_jacobi call
	flop = iter * 16.0 * param.act_res * param.act_res;

	// stopping time
	runtime = wtime() - runtime;

	fprintf(stderr, "Resolution: %5u, ", param.act_res);
	fprintf(stderr, "Time: %04.3f ", runtime);
	fprintf(stderr, "(%3.3f GFlop => %6.2f MFlop/s, ", 
		flop/1000000000.0,
		flop/runtime/1000000);
	fprintf(stderr, "residual %f, %d iterations)\n", residual, iter);
	printf("%5d %f\n", param.act_res, flop/runtime/1000000);
        if(param.act_res + param.res_step_size > param.max_res) break;
        param.act_res += param.res_step_size;
	
      }

    coarsen( param.u, np, np,
	     param.uvis, param.visres+2, param.visres+2 );
  
    write_image( resfile, param.uvis,  
		 param.visres+2, 
		 param.visres+2 );

    finalize( &param );

    totaltime = wtime() - totaltime;
    printf("\nTotal Program run time = %lf\n",totaltime );	
    return 0;
}
