#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <unistd.h>

#define N 4

int size ;

/* matrix multiplication */
void
matmult ( double  *A, double *B, double *C, int n)
{
  int i,j,k ;

  for (i=0 ; i<n ; i++)
    for (j=0 ; j<n ; j++)
      for (k=0 ; k<n ; k++) {
	C[i*n +j] += A[i*n + k] * B[k*n + j] ;
      }
}

/* matrix initialization with random numbers */
void
init(double *M, int size){
  int i,j;
  for(i=0;i<size;i++)
    for(j=0;j<size;j++)
      M[i*size+j] = rand()%100;

}

/* displays a matrix */
void
print(double *M, int size){
  int i,j;
  for(i=0;i<size;i++){
    for(j=0;j<size;j++)
      printf(" %2f ",M[i*size+j]);
    printf("\n");
  } 
}

/* Fox's algorithm */
int
main (int argc, char** argv)
{
  int p, sp,k ;
  int coords[2], pos[2], dims[2]={0,0}, periods[2]={0,0};
  int i,j,subrankrow,subrankcol,n,rank;
  double *myB, *myA, *A, *B, *C, *tmp, *MA, *MB ;
  double timer;
  MPI_Status status;
  MPI_Datatype block;
  MPI_Comm proc_grid, proc_row, proc_col;
  MPI_Request request_1, request_2;

  MPI_Init(&argc, &argv);               /* Initialize MPI               */
  timer =MPI_Wtime();
  MPI_Comm_size(MPI_COMM_WORLD, &p); /* Get the number of processors */

  sp = sqrt (p) ;
  if(p != sp*sp){
    printf("p must be a perfect square\n");
    return EXIT_FAILURE;
  }
  if(argc != 2){
    printf("you must specify the size of the matrices\n");
    return EXIT_FAILURE;
  }
  
  size = atoi(argv[1]);
  if(size % sp){
    printf("must be divisible par srqt(p)\n");
    return EXIT_FAILURE;
  }

  n = size/sp;

  myA = malloc (n*n*sizeof(double)) ;
  myB = malloc(n*n*sizeof(double));
  A = malloc (n*n*sizeof(double)) ;
  B = malloc (n*n*sizeof(double)) ;
  C = calloc (n*n, sizeof(double)) ;

  /* Create a virtual 2D-grid topology */
  MPI_Dims_create(p, 2, dims);
  MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &proc_grid);
  MPI_Comm_rank(proc_grid, &rank);    /* Note: use proc_grid */
  /* gets the coords of the processor in the grid */
  MPI_Cart_coords(proc_grid, rank, 2, coords); 

  /* Create a communicator for each row */
  MPI_Comm_split(proc_grid, coords[0], coords[1], &proc_row);
  MPI_Comm_rank(proc_row, &subrankrow); 

  /* Create a communicator for each column */ 
  MPI_Comm_split(proc_grid, coords[1], coords[0], &proc_col);
  MPI_Comm_rank(proc_col, &subrankcol);

  if (rank == 0){

      /* Create a datatype to represent a block  */
    MPI_Type_vector(n, n, size, MPI_DOUBLE, &block);
    MPI_Type_commit(&block);

      srand(time(NULL));
      MA = malloc(size*size*sizeof(double));
      MB = malloc(size*size*sizeof(double));

      init(MA, size);
      //printf("Matrix A:\n");
      //print(MA,size);

      init(MB, size);
      // printf("Matrix B:\n");
      //print(MB,size);
      // the matrices are created, we send each block to its "owner" using non blocking communication
      for (i = 0; i < sp; i++ )
	for (j = 0; j < sp; j++ ) {
	      int target ;
	      pos[0] = i ;
	      pos[1] = j ;

	      MPI_Cart_rank(proc_grid, pos, &target) ;
	      MPI_Isend(&MA[(i*size*n)+(j*n)], 1, block, target , 1, proc_grid,&request_1);
	      MPI_Isend(&MB[(i*size*n)+(j*n)], 1, block, target , 2, proc_grid,&request_2);
	}
    }
  MPI_Recv(myA, n*n, MPI_DOUBLE, 0, 1,  proc_grid, &status );
  MPI_Recv(myB, n*n, MPI_DOUBLE, 0, 2,  proc_grid, &status );
  /* everyone has its blocks, the algorithm can start */
  for (k = 0; k < sp; k++)
    {
      int m ;
      m = (coords[0] + k) % sp ;
      if (m == subrankrow)
	{
	  MPI_Bcast (myA, n*n, MPI_DOUBLE, m, proc_row) ;
	  tmp = myA;
	}
      else
	{
	  MPI_Bcast (A, n*n, MPI_DOUBLE, m, proc_row) ;
	  tmp = A;
	}
      //before the computation, let's start the communication
      MPI_Isend(myB, n*n, MPI_DOUBLE, ((subrankcol+sp-1)%sp), subrankcol+sp*k, proc_col,&request_1);
      MPI_Irecv(B, n*n, MPI_DOUBLE, ((subrankcol+1)%sp), ((subrankcol+1)%sp)+sp*k, proc_col, &request_2);
      matmult (tmp, myB, C, n) ;
      MPI_Wait(&request_1, &status);
      MPI_Wait(&request_2, &status);
      // now, we work using the newly received B
      tmp = myB;
      myB = B;
      B = tmp;		
    }
  //printf("results\n");
  // gathering the results
  MPI_Isend(C, n*n, MPI_DOUBLE, 0, 1, proc_grid, &request_1);
  if (rank == 0)
    {
      double *MC;
      MC = malloc(size*size*sizeof(double));

      for(i=0; i < p; i++){
	MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, proc_grid, &status);
	MPI_Cart_coords(proc_grid, status.MPI_SOURCE, 2, coords);
	MPI_Recv(&MC[(coords[0]*size*n)+(coords[1]*n)], 1, block, status.MPI_SOURCE, 1, proc_grid, &status); 
      }
      /* displays the result */
      //printf("Result:\n");
      //print(MC, size);
      printf("\nTime needed: %lf\n",MPI_Wtime()-timer);
      free(MA); free(MB); free(MC);
    }

  MPI_Wait(&request_1, &status);
  MPI_Finalize() ;
  free(myA); free(myB); free(A); free(B); free(C);
	
}
