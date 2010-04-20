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

	int
main (char **argv, int argc)
{
	int p, sp,k ;
	int coords[2], pos[2], reorder=1, dims[2]={0,0}, periods[2]={0,0};
	int myid,i,j,subrankrow,subrankcol,rank,n;
	       double *myB, *myA, *A, *B, *C, *tmp, *AA ;
	MPI_Status status;
	MPI_Datatype block;
	MPI_Comm proc_grid, proc_row, proc_col;
	MPI_Request send, recv;

	MPI_Init(&argc, &argv);               /* Initialize MPI               */
	double timer1 =MPI_Wtime();
	MPI_Comm_size(MPI_COMM_WORLD, &p); /* Get the number of processors */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */

	sp = sqrt (p) ;
	assert (p == sp*sp) ;
	//size = atoi (argv[1]) ;
	size = N ;
	assert (0 == size % sp);

	n = size/sp;

	myA = malloc (n*n*sizeof(double)) ;
	myB = malloc(n*n*sizeof(double));
	A = malloc (n*n*sizeof(double)) ;
	B = malloc (n*n*sizeof(double)) ;
	C = calloc (n*n, sizeof(double)) ;

	/* Create a virtual 2D-grid topology */
	MPI_Dims_create(p, 2, dims);
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &proc_grid);
	MPI_Comm_rank(proc_grid, &myid);    /* Note: use proc_grid */
	/* gets the coords of the processor in the grid */
	MPI_Cart_coords(proc_grid, myid, 2, coords); 

	/* Create a communicator for each row */
	MPI_Comm_split(proc_grid, coords[0], coords[1], &proc_row);
	MPI_Comm_rank(proc_row, &subrankrow); 

	/* Create a communicator for each column */ 
	MPI_Comm_split(proc_grid, coords[1], coords[0], &proc_col);
	MPI_Comm_rank(proc_col, &subrankcol); 

	/* Create a datatype to represent a block  */
	MPI_Type_vector(n, n, size, MPI_DOUBLE, &block);
	MPI_Type_commit(&block);

	if (!coords[0] && !coords[1])
	{
	  double MA[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16} , MB[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16} ;
	  // the matrices are created, we send each block to its "owner"
		for (i = 0; i < sp; i++ )
			for (j = 0; j < sp; j++ ) {
				if ( !i && 	!j )
				{
					int i,j ;
					for(i=0; i < n; i++ )
						for(j=0; j < n; j++ )
						{
							myA[i*n+j] = MA[i*size + j] ;
							myB[i*n+j] = MB[i*size + j] ;
						}						
				}
				else
				{
					int rank ;
					pos[0] = i ;
					pos[1] = j ;

					MPI_Cart_rank(proc_grid, pos, &rank) ;
					MPI_Send(&MA[(i*size*n)+(j*n)], 1, block, rank , 111, proc_grid);
					MPI_Send(&MB[(i*size*n)+(j*n)], 1, block, rank , 111, proc_grid);
				}
			}
	}
	else
	{

		int rank ;
		pos[0] = 0 ;
		pos[1] = 0 ;

		MPI_Cart_rank(proc_grid, pos, &rank) ;
		MPI_Recv(myA, n*n, MPI_DOUBLE, rank, 111,  proc_grid, &status );
		MPI_Recv(myB, n*n, MPI_DOUBLE, rank, 111,  proc_grid, &status );
	}

		/* everyone has its blocks, the algorithm can start */
	for (k = 0; k < sp; k++)
	{
		int m ;
		m = (coords[0] + k) % sp ;
		if (m == subrankrow)
		{
			MPI_Bcast (myA, n*n, MPI_DOUBLE, m, proc_row) ;
			//memcpy (A, myA, n*n*sizeof(double)) ;
			AA = myA;
		}
		else
		{
			MPI_Bcast (A, n*n, MPI_DOUBLE, m, proc_row) ;
			AA = A;
		}
		//before the computation, let's start the communication
		MPI_Isend(myB, n*n, MPI_DOUBLE, ((subrankcol+sp-1)%sp), subrankcol+sp*k, proc_col,&send);
		MPI_Irecv(B, n*n, MPI_DOUBLE, ((subrankcol+1)%sp), ((subrankcol+1)%sp)+sp*k, proc_col, &recv);
		matmult (AA, myB, C, n) ;
		MPI_Wait(&send, &status);
		MPI_Wait(&recv, &status);
		// now, we work using the newly received B
		tmp = myB;
		myB = B;
		B = tmp;		
	}

	// building the result
	if (!coords[0] && !coords[1])
	{
	  double* MC = calloc (size*size, sizeof(double)) ;
		for (i = 0; i < sp; i++ )
			for (j = 0; j < sp; j++ ) {
				if ( !i && 	!j )
				{
					int i,j ;
					for(i=0; i < n; i++ )
						for(j=0; j < n; j++ )
						{
						  MC[(i*size)+ j] = C[i*n + j] ;
						}
				}
				else
				{
					int rank ;
					pos[0] = i ;
					pos[1] = j ;

					MPI_Cart_rank(proc_grid, pos, &rank) ;
					MPI_Recv(&MC[(i*size*n)+(j*n)], 1, block, rank, 111,  proc_grid, &status );
				}
			}
	      	// displays the result matrix
		if (!coords[0] && !coords[1]) {
		  for (i=0;i<size;i++){
		    for(j=0;j<size;j++)
		      printf("%f ",MC[i*size + j]);
		    printf("\n");
		  }
		  // free(MA); free(MB); free(MC);
		  free(MC);
		}
	}
	else
	{
		int rank ;
		pos[0] = 0 ;
		pos[1] = 0 ;

		MPI_Cart_rank(proc_grid, pos, &rank) ;
		MPI_Send(C, n*n, MPI_DOUBLE, rank, 111,  proc_grid);
	}
	
printf("time needed: %lf\n",MPI_Wtime() - timer1);

	MPI_Finalize() ;
	free(myA); free(myB); free(A); free(B); free(C);
	
}
