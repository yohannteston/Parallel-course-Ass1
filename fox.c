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
double *MC;

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
	       double MA[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16} , MB[16] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16} ;
	double *myA, *A, *B, *C ;
	MPI_Status status;
	MPI_Datatype block;
	MPI_Comm proc_grid, proc_row, proc_col;

	MPI_Init(&argc, &argv);               /* Initialize MPI               */
	MPI_Comm_size(MPI_COMM_WORLD, &p); /* Get the number of processors */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); /* Get my number                */

	sp = sqrt (p) ;
	assert (p == sp*sp) ;
	//size = atoi (argv[1]) ;
	size = N ;
	assert (0 == size % sp);

	n = size/sp;

	myA = malloc (n*n*sizeof(double)) ;
	A = malloc (n*n*sizeof(double)) ;
	B = malloc (n*n*sizeof(double)) ;
	C = calloc (n*n, sizeof(double)) ;

	/* Create a virtual 2D-grid topology */
	MPI_Dims_create(p, 2, dims);
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &proc_grid);
	MPI_Comm_rank(proc_grid, &myid);    /* Note: use proc_grid */

	/* Create a communicator for each row */
	MPI_Cart_coords(proc_grid, myid, 2, coords);
	MPI_Comm_split(proc_grid, coords[0], coords[1], &proc_row);
	MPI_Comm_rank(proc_row, &subrankrow); 

	/* Create a communicator for each column */ 
	MPI_Comm_split(proc_grid, coords[1], coords[0], &proc_col);
	MPI_Comm_rank(proc_col, &subrankcol); 

	/*   */
	MPI_Type_vector(n, n, size, MPI_DOUBLE, &block);
	MPI_Type_commit(&block);

	if (!coords[0] && !coords[1])
	{
		MC = calloc (size*size, sizeof(double)) ;

		for (i = 0; i < sp; i++ )
			for (j = 0; j < sp; j++ ) {
				if ( !i && 	!j )
				{
					int i,j ;
					for(i=0; i < n; i++ )
						for(j=0; j < n; j++ )
						{
							myA[i*n+j] = MA[i*size + j] ;
							B[i*n+j] = MB[i*size + j] ;
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
		MPI_Recv(B, n*n, MPI_DOUBLE, rank, 111,  proc_grid, &status );
	}

		/* everyone has its blocks, the algorithm can start */
	for (k = 0; k < sp; k++)
	{
		int m ;
		m = (coords[0] + k) % sp ;
		if (m == subrankrow)
		{
			MPI_Bcast (myA, n*n, MPI_DOUBLE, m, proc_row) ;
			memcpy (A, myA, n*n*sizeof(double)) ;
		}
		else
		{
			MPI_Bcast (A, n*n, MPI_DOUBLE, m, proc_row) ;
		}
		matmult (A, B, C, n) ;
		
		MPI_Sendrecv_replace (B, n*n, MPI_DOUBLE, ((subrankcol+sp-1)%sp), subrankcol+sp*k, ((subrankcol+1)%sp), ((subrankcol+1)%sp)+sp*k, proc_col, &status) ;
	}

	if (!coords[0] && !coords[1])
	{

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
	}
	else
	{
		int rank ;
		pos[0] = 0 ;
		pos[1] = 0 ;

		MPI_Cart_rank(proc_grid, pos, &rank) ;
		MPI_Send(C, n*n, MPI_DOUBLE, rank, 111,  proc_grid);
	}

	if (!coords[0] && !coords[1])
	  for (i=0;i<size;i++){
	    for(j=0;j<size;j++)
			printf("%f ",MC[i*size + j]);
		  printf("\n");
	  }
	

	MPI_Finalize() ;
}
