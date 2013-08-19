#include <stdio.h>
#include "mpi.h"

int main(int argc, char * argv[] )
{

int NumTask = 0;
int rank = 0;

int status = MPI_Init(&argc, &argv);

if( status != MPI_SUCCESS )
{
   printf("\nError starting MPI Program.. terminating\n");
   MPI_Abort(MPI_COMM_WORLD, status);
}

MPI_Comm_size(MPI_COMM_WORLD, &NumTask);

MPI_Comm_rank(MPI_COMM_WORLD, &rank);

if( rank == 5 )
{
	printf("\nHello World! Processor %d\n", rank);
}

MPI_Finalize();
return 0;
}
