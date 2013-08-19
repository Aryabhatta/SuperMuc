#include <stdio.h>
#include <math.h>
#include "mpi.h"

// Multiple iterations to take average values
#define ITER 1000000
int main(int argc, char * argv[] )
{

int NumTask = 0;
int rank = 0;

// Destination & source taken through arguments
int dest = atoi(argv[2]);
int source = atoi(argv[1]);
long i=0,j=0;

int status = MPI_Init(&argc, &argv);

// out array for source, in for destination
char *out, *in;
long bufsize=1;
double time1, time2;
double nettime=0, bandw=0;

MPI_Status stat;

if( status != MPI_SUCCESS )
{
   printf("\nError starting MPI Program.. terminating\n");
   MPI_Abort(MPI_COMM_WORLD, status);
}

MPI_Comm_size(MPI_COMM_WORLD, &NumTask);

MPI_Comm_rank(MPI_COMM_WORLD, &rank);

if( rank == 0 )
{
  printf("\n");
  printf("\n time(s)  bandwidth(MB/s), bufsize(bytes)" );
}

for( j=0; j<21; j++)
{
  bufsize = pow(2,j);
  out = (char*)malloc(sizeof(char)*bufsize);
  in = (char*)malloc(sizeof(char)*bufsize);
 
  // Initialising buffer
  for(i=0;i<bufsize;i++)
  {
    out[i] = '1';
    in[i] = '2';
  }

  // Synchronize processes before send receive
  MPI_Barrier(MPI_COMM_WORLD);

  if( rank == source )
  {
        time1 = MPI_Wtime();
	for( i=0; i<ITER; i++)
	{
          status = MPI_Send(out, bufsize, MPI_CHAR, dest, 5, MPI_COMM_WORLD);
          status = MPI_Recv(out,bufsize,MPI_CHAR, dest, 5,MPI_COMM_WORLD,&stat);
	}
	time2 = MPI_Wtime();

	// subtracting latency from nettime	
	nettime = (((time2-time1)/ITER) - 0.00000355);
	if( nettime < 0 )
		nettime = fabs(nettime);

	// calculate bandwidth, nettime divided by 2 (it is roundtrip netime)
	bandw = (2*bufsize / (nettime*1024*1024));
	printf("\n %.10lf \t %.10lf \t %ld",nettime, bandw, bufsize );
  }

  if( rank == dest )
  {
     for( i=0; i<ITER; i++)
     {
      status = MPI_Recv(in,bufsize,MPI_CHAR, source, 5,MPI_COMM_WORLD,&stat);
      status = MPI_Send(in, bufsize, MPI_CHAR, source, 5, MPI_COMM_WORLD);
     }
  }
  MPI_Barrier(MPI_COMM_WORLD);
 
  // deallocate memory 
  if(j>=0)
  {
	free(out);
	free(in);
  }
}

MPI_Finalize();
return 0;
}
