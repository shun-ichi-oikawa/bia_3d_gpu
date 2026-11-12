#include <stdio.h>
#include <cutil_inline.h>
#include <cublas.h>
#define BLOCKS	256

void   SetTimer(unsigned int *t) {	 cutCreateTimer(t); cutStartTimer(*t); }
double EndTimer(unsigned int *t)
{	cutStopTimer(*t);	double ms =	cutGetTimerValue(*t); cutDeleteTimer(*t);
	return ms;
}

__global__
void inv_device ( const int N, double *x, double *y, double *z, double *p )
{	unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if( tid < N - 1 )
	{	double sum = 0; double xi = x[tid]; double yi = y[tid]; double zi = z[tid];
		for ( int j = tid + 1; j < N; j++ )
		{	double xj = x[j]; double yj = y[j]; double zj = z[j];
//			sum += 1 / sqrt( ( xi - xj ) * ( xi - xj ) + ( yi - yj ) * ( yi - yj ) + ( zi - zj ) * ( zi - zj ) );
			sum +=    rsqrt( ( xi - xj ) * ( xi - xj ) + ( yi - yj ) * ( yi - yj ) + ( zi - zj ) * ( zi - zj ) );
		}	p[tid] = sum;
	}	__syncthreads();
}
//
__global__ void sum_device
	( const int N, // const double eps_th,
	  double *x, double *y, double *z, double *u, double *v, double *w,
	  double *p, double *c
	)
{	unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
	double sum = 0;
	if( tid ==  0 ) { for ( int j = 0; j < N  ; j++ ) { sum += x[j]; } }
	if( tid ==  1 ) { for ( int j = 0; j < N  ; j++ ) { sum += y[j]; } }
	if( tid ==  2 ) { for ( int j = 0; j < N  ; j++ ) { sum += z[j]; } }
	if( tid ==  3 ) { for ( int j = 0; j < N  ; j++ ) { sum += u[j]; } }
	if( tid ==  4 ) { for ( int j = 0; j < N  ; j++ ) { sum += v[j]; } }
	if( tid ==  5 ) { for ( int j = 0; j < N  ; j++ ) { sum += w[j]; } }
	if( tid ==  6 ) { for ( int j = 0; j < N  ; j++ ) { sum += y[j]*w[j] - z[j]*v[j]; } }
	if( tid ==  7 ) { for ( int j = 0; j < N  ; j++ ) { sum += z[j]*u[j] - x[j]*w[j]; } }
	if( tid ==  8 ) { for ( int j = 0; j < N  ; j++ ) { sum += x[j]*v[j] - y[j]*u[j]; } }
	if( tid ==  9 ) { for ( int j = 0; j < N  ; j++ ) { sum += u[j]*u[j] + v[j]*v[j] + w[j]*w[j]; } }
//	if( tid == 10 ) { for ( int j = 0; j < N-1; j++ ) { sum += p[j]; } }
	if( tid == 10 ) { for ( int j = N-2; j >=0; j-- ) { sum += p[j]; } }		// more accurate than above
	c[tid] = sum;
}

extern "C" void invariance_gpu_
	( const int N, const double eps_th,
	  double *X , double *Y , double *Z , double *U , double *V , double *W,
	  double *KE, double *PE, double *Lx, double *Ly, double *Lz,
	  double *Ua, double *Va, double *Wa, double *Xa, double *Ya, double *Za
	)
{	int GRIDS = N / BLOCKS; if ( N % BLOCKS != 0 ) GRIDS += 1;
	printf("N_p, BLOCKS, GRIDS = %d, %d, %d\n", N, BLOCKS, GRIDS);
	unsigned int timer = 0; SetTimer(&timer);
	double *x; cudaMalloc( (void **) &x, sizeof(double) * N );
	double *y; cudaMalloc( (void **) &y, sizeof(double) * N );
	double *z; cudaMalloc( (void **) &z, sizeof(double) * N );
	double *u; cudaMalloc( (void **) &u, sizeof(double) * N );
	double *v; cudaMalloc( (void **) &v, sizeof(double) * N );
	double *w; cudaMalloc( (void **) &w, sizeof(double) * N );
	cudaMemcpy( x, X, sizeof(double) * N, cudaMemcpyHostToDevice);
	cudaMemcpy( y, Y, sizeof(double) * N, cudaMemcpyHostToDevice);
	cudaMemcpy( z, Z, sizeof(double) * N, cudaMemcpyHostToDevice);
	cudaMemcpy( u, U, sizeof(double) * N, cudaMemcpyHostToDevice);
	cudaMemcpy( v, V, sizeof(double) * N, cudaMemcpyHostToDevice);
	cudaMemcpy( w, W, sizeof(double) * N, cudaMemcpyHostToDevice);
	double *P =            (double*) malloc( sizeof(double) * ( N - 1 ) );
	double *p; cudaMalloc( (void **)     &p, sizeof(double) * ( N - 1 ) );
	//
	inv_device <<< GRIDS, BLOCKS>>> ( N, x, y, z, p );
	//
	cudaMemcpy( P, p, sizeof(double) * ( N - 1 ), cudaMemcpyDeviceToHost);
	//
	SetTimer(&timer);
	double *C =            (double*) malloc( sizeof(double) * 11 );
	double *c; cudaMalloc( (void **) &c,     sizeof(double) * 11 );
	//
	sum_device <<<1,11>>> ( N, x, y, z, u, v, w, p, c );
	//
	cudaMemcpy( C, c, sizeof(double) * 11, cudaMemcpyDeviceToHost);
	*Xa = C[0]; *Ya = C[1]; *Za = C[2];
	*Ua = C[3]; *Va = C[4]; *Wa = C[5];
	*Lx = C[6]; *Ly = C[7]; *Lz = C[8];
	*KE = C[9]; *PE = eps_th * C[10];
	//
	/*
	printf("CUDA  PE =%25.16e   \n", *PE);
	double *ONE = (double*) malloc( sizeof(double) ); ONE[0] = 1;
	double *one; cudaMalloc( (void **) &one, sizeof(double) );
	cudaMemcpy( one, ONE, sizeof(double), cudaMemcpyHostToDevice);
	printf("Dasum PE =%25.16e\n", eps_th * cublasDasum( N - 1,         p, 1 ) );
	printf("Ddot  PE =%25.16e\n", eps_th * cublasDdot ( N - 1, one, 0, p, 1 ) );
	printf("Ddot  Xa =%25.16e\n",          cublasDdot ( N    , one, 0, x, 1 ) );
	printf("Ddot  Ya =%25.16e\n",          cublasDdot ( N    , one, 0, y, 1 ) );
	printf("Ddot  Za =%25.16e\n",          cublasDdot ( N    , one, 0, z, 1 ) );
	printf("Ddot  Ua =%25.16e\n",          cublasDdot ( N    , one, 0, u, 1 ) );
	printf("Ddot  Va =%25.16e\n",          cublasDdot ( N    , one, 0, v, 1 ) );
	printf("Ddot  Wa =%25.16e\n",          cublasDdot ( N    , one, 0, w, 1 ) );
	cudaFree(one); free(ONE);
	//*/
	cudaFree( x ); cudaFree( y ); cudaFree( z );
	cudaFree( u ); cudaFree( v ); cudaFree( w );
	cudaFree( p ); cudaFree( c );
	printf("計算時間 = %8.5f sec\n", EndTimer(&timer) / 1000);
	return;
}
