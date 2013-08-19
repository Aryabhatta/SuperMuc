#include <stdio.h>
#include <math.h>
#include "mpi.h"

int main(int argc, char * argv[] )
{
int np = 0; // #process
int rank = 0;

// arraysize taken through argument
int ArraySz = atoi(argv[1]);

int status = MPI_Init(&argc, &argv);

MPI_Comm_size(MPI_COMM_WORLD, &np);

int localsum = 0; // localsum @ core
int RecvBuf = 0;  // Receive buffer

// size of array to be initialized by each process
int localSz;

// ArraySize exact multiple of # process
if( ArraySz % np == 0 )
{ 
   localSz = ArraySz / np;
}
else if( ArraySz < np )// ArraySize less than # process
{
   localSz = 0;
}
else // ArraySz not exact multiple, load balancing through ceil
{
   localSz = ceil((double)ArraySz/(double)np);
}
int A[localSz];

int dest;
int source;
int elem;
int i=0,j=0;

double time1, time2;
double nettime=0;

MPI_Status stat;

if( status != MPI_SUCCESS )
{
   printf("\nError starting MPI Program.. terminating\n");
   MPI_Abort(MPI_COMM_WORLD, status);
}

MPI_Comm_rank(MPI_COMM_WORLD, &rank);

if(rank==0)
{
   time1 = MPI_Wtime();
}

if( ArraySz < np ) // only localsum set to 'i' (till arraySz else 0)
{
  for( i=0; i < np; i++)
  {
   if( rank==i  && rank < ArraySz)
   {
      localsum = rank;
   }
  }
}
else if (ArraySz %np ==0)
{
  for( i=0; i < np; i++)
  {
   if( rank==i )
   {
	for( j=0; j < localSz;j++ )
	{
	  A[j] = rank * localSz + j;
	  localsum += A[j];
	}
   }
  }
}
else
{
  for( i=0; i < np; i++)
  {
   if( rank==i )
   {
	for( j=0; j < localSz;j++ )
	{
          elem = rank * localSz + j;
          if(elem < ArraySz)
          { 	  
		A[j] = elem;
          }
	  else //after valid arraySz, buffer=0
	  {
		A[j] = 0;
	  }
	   localsum += A[j];
	}
   }
  }
}

// let every process initialise localsum before starting reduction
MPI_Barrier( MPI_COMM_WORLD );

/* 
// Print
for( i=0;i<np;i++ )
  if( rank == i )
  {
    printf("\nrank=%d", rank);
   for(j=0;j<localSz;j++)
    printf("\tA[%d]=%d",j,A[j]);  
   printf("\tlocalsum=%d\n", localsum);
  }
MPI_Barrier( MPI_COMM_WORLD );
*/

  // REDUCTION TREE
  for( i=0; i<log(np); i++ )
  {
      
   if( (rank % (int)pow(2,i) == 0) && ((rank/(i+1)) %2 == 0 ) )
   {
     source = rank + pow(2,i);
     MPI_Recv( &RecvBuf, 1, MPI_INT,source ,10, MPI_COMM_WORLD, &stat);
     localsum += RecvBuf;
   }
   
   if ( (rank % (int)pow(2,i) == 0 ) && ((rank/(i+1)) %2==1) )
   {
     dest = rank - pow(2,i);
     MPI_Send(&localsum, 1 , MPI_INT, dest , 10, MPI_COMM_WORLD);
   }
  }

MPI_Barrier( MPI_COMM_WORLD );

if( rank == 0)
{
  printf("\nSum = %d\n", localsum);
  time2 = MPI_Wtime();
  printf("\nTime taken by program:%.10lf\n", (time2-time1));
}

MPI_Finalize();
return 0;
}
