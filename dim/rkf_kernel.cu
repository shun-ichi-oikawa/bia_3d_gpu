#include<stdio.h>
#include<stdlib.h>
#include<algorithm>
#include<limits>
#include<cutil_inline.h>
#define BLOCKS		128
//#define BLOCKS		192
//#define BLOCKS		256
//#define BLOCKS		512
#define DEBUG

#define GPU

__device__ int    NP,NEQ;
__device__ double EPS_TH;
int np; double eps;

#ifdef GPU
const double AL2 = 1e0 / 4, AL3 = 3e0 / 8, AL4 = 12e0 / 13, AL5 = 1e0, AL6 = 1e0 / 2;
#endif

__global__
void SetConst ( const int n_p, const double eps_th)
{	unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if ( tid == 0 ) { NP = n_p; NEQ = 6 * n_p; EPS_TH = eps_th;}
	__syncthreads();
}

__global__
void func( const double x, const double *y, double *f)
{	unsigned int ip  = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int TID =                           threadIdx.x;
	__shared__ double XI[BLOCKS], YI[BLOCKS], ZI[BLOCKS], AX[BLOCKS], AY[BLOCKS], AZ[BLOCKS];
	if( ip < NP )
	{	int IP  = 6 * ip;
		f[IP+0] = y[IP+3]; f[IP+1] = y[IP+4]; f[IP+2] = y[IP+5];
		XI[TID] = y[IP+0]; YI[TID] = y[IP+1]; ZI[TID] = y[IP+2];
		AX[TID] =      0 ; AY[TID] =      0 ; AZ[TID] =      0 ;
		for( int jp = 0; jp < NP; jp++)
		{	int JP = 6 * jp;
			double xij = XI[TID] - y[JP+0];
			double yij = YI[TID] - y[JP+1];
			double zij = ZI[TID] - y[JP+2];
			double rsq = xij * xij + yij * yij + zij * zij + 1e-99;
			double eee = EPS_TH / ( rsq * sqrt(rsq) );
			AX[TID] += xij * eee;
			AY[TID] += yij * eee;
			AZ[TID] += zij * eee;
		}
		f[IP+3] = AX[TID]; f[IP+4] = AY[TID]; f[IP+5] = AZ[TID];
	}
	__syncthreads();
}

void moon( const double x, const double *y, double *f)
{	for ( int ip = 0; ip < np; ip++)
	{	int IP  = 6 * ip;
		f[IP] = y[IP+3]; f[IP+1] = y[IP+4]; f[IP+2] = y[IP+5];
		double Ex_i = 0, Ey_i = 0, Ez_i = 0;
		for( int jp = 0; jp < np; jp++)
		{	int JP = 6 * jp;
			double xij = y[IP+0] - y[JP+0];
			double yij = y[IP+1] - y[JP+1];
			double zij = y[IP+2] - y[JP+2];
			double rsq = xij * xij + yij * yij + zij * zij + 1e-99;
			double eee = eps / ( rsq * sqrt(rsq) );
			Ex_i += xij * eee;
			Ey_i += yij * eee;
			Ez_i += zij * eee;
		}
		f[IP+3] = Ex_i; f[IP+4] = Ey_i; f[IP+5] = Ez_i;
	}
}

__global__
void fehl_1 ( const double H, const double *y0, double *wk)
{	unsigned int ip = blockIdx.x * blockDim.x + threadIdx.x;
	const double B21 =    1e0 /    4;
	if( ip < NP )
	{	const int BASE0 = 6 * ip, BASE6 = BASE0 + 6 * NEQ;
		for( int i = 0; i < 6; i++)
		{	wk[BASE6+i] = y0[BASE0+i] + H * B21 * wk[BASE0+i];}
	}	__syncthreads();
}
__global__
void fehl_2 ( const double H, const double *y0, double *wk)
{	unsigned int ip = blockIdx.x * blockDim.x + threadIdx.x;
	const double B31 =    3e0 /   32, B32 =     9e0 /    32;
	if( ip < NP )
	{	const int BASE0 = 6 * ip, BASE1 = BASE0 + NEQ, BASE7 = BASE0 + 7 * NEQ;
		for( int i = 0; i < 6; i++)
		{	wk[BASE7+i] = y0[BASE0+i] + H * ( B31 * wk[BASE0+i] + B32 * wk[BASE1+i]);}
	}
	__syncthreads();
}
__global__
void fehl_3 ( const double H, const double *y0, double *wk)
{	unsigned int ip = blockIdx.x * blockDim.x + threadIdx.x;
	const double B41 = 1932e0 / 2197, B42 = -7200e0 /  2197, B43 =  7296e0 /  2197;
	if( ip < NP )
	{	const int BASE0 = 6 * ip, BASE1 = BASE0 + NEQ, BASE2 = BASE1 + NEQ, BASE8 = BASE0 + 8 * NEQ;
		for( int i = 0; i < 6; i++)
		{	wk[BASE8+i] = y0[BASE0+i] + H * ( B41 * wk[BASE0+i] + B42 * wk[BASE1+i] + B43 * wk[BASE2+i] );}
	}	__syncthreads();
}
__global__
void fehl_4 ( const double H, const double *y0, double *wk)
{	unsigned int ip = blockIdx.x * blockDim.x + threadIdx.x;
	const double B51 =  439e0 /  216, B52 =    -8e0        , B53 =  3680e0 /   513, B54 = -845e0 / 4104;
	if( ip < NP )
	{	const int BASE0 = 6 * ip, BASE1 = BASE0 + NEQ, BASE2 = BASE1 + NEQ, BASE3 = BASE2 + NEQ, BASE9 = BASE0 + 9 * NEQ;
		for( int i = 0; i < 6; i++)
		{	wk[BASE9+i] = y0[BASE0+i] + H * ( B51 * wk[BASE0+i] + B52 * wk[BASE1+i] + B53 * wk[BASE2+i] + B54 * wk[BASE3+i]);}
	}	__syncthreads();
}
__global__
void fehl_5 ( const double H, const double *y0, double *wk)
{	unsigned int ip = blockIdx.x * blockDim.x + threadIdx.x;
	const double B61 =   -8e0 /   27, B62 =     2e0        , B63 = -3544e0 /  2565, B64 = 1859e0 / 4104, B65 = -11e0 / 40;
	if( ip < NP )
	{	const int BASE0 = 6 * ip, BASE1 = BASE0 + NEQ, BASE2 = BASE1 + NEQ, BASE3 = BASE2 + NEQ, BASE4 = BASE3 + NEQ, BASE10 = BASE0 + 10 * NEQ;
		for( int i = 0; i < 6; i++)
		{	wk[BASE10+i]= y0[BASE0+i] + H * ( B61 * wk[BASE0+i] + B62 * wk[BASE1+i] + B63 * wk[BASE2+i] + B64 * wk[BASE3+i] + B65 * wk[BASE4+i]);}
	}	__syncthreads();
}
__global__
void fehl_6 ( const double H, const double *y0, double *wk, double *yn)
{	unsigned int ip = blockIdx.x * blockDim.x + threadIdx.x;
	const double GA1 =   16e0 /  135, GA3 =  6656e0 / 12825, GA4 = 28561e0 / 56430, GA5 =   -9e0 /   50, GA6 =   2e0 / 55;
	if( ip < NP )
	{	const int BASE0 = 6 * ip, BASE1 = BASE0 + NEQ, BASE2 = BASE1 + NEQ, BASE3 = BASE2 + NEQ, BASE4 = BASE3 + NEQ, BASE5 = BASE4 + NEQ;
		for( int i = 0; i < 6; i++)
		{	yn[BASE0+i] = y0[BASE0+i] + H * ( GA1 * wk[BASE0+i]                     + GA3 * wk[BASE2+i] + GA4 * wk[BASE3+i] + GA5 * wk[BASE4+i] + GA6 * wk[BASE5+i]);}
	}	__syncthreads();
}
__global__
void fehl_e ( const double H, double *er, const double *wk)
{	unsigned int ip = blockIdx.x * blockDim.x + threadIdx.x;
	const double DA1 =    1e0 /  360, DA3 =  -128e0 /  4275, DA4 = -2197e0 / 75240, DA5 =    1e0 /   50, DA6 =   2e0 / 55;
	if( ip < NP )
	{	const int BASE0 = 6 * ip, BASE1 = BASE0 + NEQ, BASE2 = BASE1 + NEQ, BASE3 = BASE2 + NEQ, BASE4 = BASE3 + NEQ, BASE5 = BASE4 + NEQ;
		for( int i = 0; i < 6; i++)
		{	er[BASE0+i] =               H * ( DA1 * wk[BASE0+i]                     + DA3 * wk[BASE2+i] + DA4 * wk[BASE3+i] + DA5 * wk[BASE4+i] + DA6 * wk[BASE5+i]);}
	}	__syncthreads();
}

__global__
void calc_esterr( const double *yln, const double *yhn, double *esterr)
{	unsigned int ip = blockIdx.x * blockDim.x + threadIdx.x;
	if( ip < NP )
	{	const int BASE = 6 * ip;
		for( int i = 0; i < 6; i++)
		{	int I = BASE + i;
//			esterr[I] = ( yln[I] - yhn[I] ) / 31;
			esterr[I] = ( yln[I] - yhn[I] ) * 0.032258064516129032258064516129032;
		}
	}
	__syncthreads();
}

void fehl(const int LINDEX, const int NEQ, const double X, const double H, double *Y0, double *YN, double *ERR, double *WK)
{	const double B21 =    1e0 /    4;
	const double B31 =    3e0 /   32, B32 =     9e0 /    32;
	const double B41 = 1932e0 / 2197, B42 = -7200e0 /  2197, B43 =  7296e0 /  2197;
	const double B51 =  439e0 /  216, B52 =    -8e0        , B53 =  3680e0 /   513, B54 = -845e0 / 4104;
	const double B61 =   -8e0 /   27, B62 =     2e0        , B63 = -3544e0 /  2565, B64 = 1859e0 / 4104, B65 = -11e0 / 40;
	const double GA1 =   16e0 /  135, GA3 =  6656e0 / 12825, GA4 = 28561e0 / 56430, GA5 =   -9e0 /   50, GA6 =   2e0 / 55;
	const double DA1 =    1e0 /  360, DA3 =  -128e0 /  4275, DA4 = -2197e0 / 75240, DA5 =    1e0 /   50, DA6 =   2e0 / 55;
	const double AL2 =    1e0 /    4, AL3 =     3e0 /     8, AL4 =    12e0 /    13, AL5 =    1e0       , AL6 =   1e0 /  2;
	moon(X      ,        Y0, WK      ); for(int i = 0; i < NEQ; i++){ WK[ 6*NEQ+i] = Y0[i] + H *   B21 * WK[i];}
	moon(X+AL2*H, WK+ 6*NEQ, WK+  NEQ); for(int i = 0; i < NEQ; i++){ WK[ 7*NEQ+i] = Y0[i] + H * ( B31 * WK[i] + B32 * WK[NEQ+i]);}
	moon(X+AL3*H, WK+ 7*NEQ, WK+2*NEQ); for(int i = 0; i < NEQ; i++){ WK[ 8*NEQ+i] = Y0[i] + H * ( B41 * WK[i] + B42 * WK[NEQ+i] + B43 * WK[2*NEQ+i]);	}
	moon(X+AL4*H, WK+ 8*NEQ, WK+3*NEQ); for(int i = 0; i < NEQ; i++){ WK[ 9*NEQ+i] = Y0[i] + H * ( B51 * WK[i] + B52 * WK[NEQ+i] + B53 * WK[2*NEQ+i] + B54 * WK[3*NEQ+i]);}
	moon(X+AL5*H, WK+ 9*NEQ, WK+4*NEQ); for(int i = 0; i < NEQ; i++){ WK[10*NEQ+i] = Y0[i] + H * ( B61 * WK[i] + B62 * WK[NEQ+i] + B63 * WK[2*NEQ+i] + B64 * WK[3*NEQ+i] + B65 * WK[4*NEQ+i]);}
	moon(X+AL6*H, WK+10*NEQ, WK+5*NEQ); for(int i = 0; i < NEQ; i++){ YN[       i] = Y0[i] + H * ( GA1 * WK[i]                   + GA3 * WK[2*NEQ+i] + GA4 * WK[3*NEQ+i] + GA5 * WK[4*NEQ+i] + GA6 * WK[5*NEQ+i]); }
	if ( LINDEX == 1 )					 { for(int i = 0; i < NEQ; i++){ERR[       i] =         H * ( DA1 * WK[i]                   + DA3 * WK[2*NEQ+i] + DA4 * WK[3*NEQ+i] + DA5 * WK[4*NEQ+i] + DA6 * WK[5*NEQ+i]);}}
}

void flstep(const int NEQ, const double X, const double H, double *Y0, double *YN, double *ERR, double *WK)
{	int LINDEX = 1;
	fehl( LINDEX, NEQ, X, H, Y0, YN, ERR, WK);
}

void fhstep(const int NEQ, const double X, const double H, double *Y0, double *YN, double *ERR, double *WK)
{	double X1 = X, H1 = H / 2; int LINDEX = 0;
	fehl( LINDEX, NEQ, X1, H1, Y0, YN, ERR, WK); X1 = X1 + H1; for(int i = 0; i < NEQ; i++){ Y0[i] = YN[i]; }
	fehl( LINDEX, NEQ, X1, H1, Y0, YN, ERR, WK);
}

void   SetTimer(unsigned int *t);
double EndTimer(unsigned int *t);

//         call RKF_C ( NEQ, eps_th,X0,XE,Y0,INIT,RELERR,ABSERR,YN,ESTERR,NDE,IER,WK0,WK1,WK2,WK3,WK4)
extern "C" void rkf_c_
//	(	const int *neq, const double *eps_th, double *X0, double *XE, double *YH0, int *INIT,
//	(	const int _neq, const double _eps_th, double *X0, double *XE, double *YH0, int *INIT,
	(	const int  neq, const double  eps_th, double *X0, double *XE, double *YH0, int *INIT,
		const double *RELERR, const double *ABSERR, double *YHN,
		double *ESE, int *NDE, int *IER, double *H, double *YL0, double *YLN, double *ERR, double *WRK
	)
{	const double EPSMIN = 1e-15;
	const int    ITEMAX = 10000;
//	np = NP = *neq / 6; NEQ = *neq; eps = *eps_th;
//	np = NP = _neq / 6; NEQ = _neq; eps = _eps_th;
	np = NP =  neq / 6; NEQ =  neq; eps =  eps_th;
	*NDE  = 0;
	#ifdef DEBUG
		#ifdef GPU
			char str[20];sprintf(str,"../OUT/%07d.DIM",np);
			FILE* fp; fp = fopen(str,"w");
		#else
			FILE* fp; fp = fopen("xyzuvw.cpu","w");
		#endif
		fprintf(fp, "%24s%25s%25s%25s%25s%25s%25s\n", "time", "x", "y", "z", "u", "v", "w");
		fprintf(fp, "%24.16e%25.16e%25.16e%25.16e%25.16e%25.16e%25.16e\n",
				  *X0, YH0[0], YH0[1], YH0[2], YH0[3], YH0[4], YH0[5] );
//		printf(       "%9.6f%20.16f%20.16f%20.16f%20.16f%20.16f%20.16f\n",
//				  *X0, YH0[0], YH0[1], YH0[2], YH0[3], YH0[4], YH0[5] );
	#endif

	#ifdef GPU
	//cudaSetDevice(1);
	int GRIDS = NP / BLOCKS; if( NP % BLOCKS != 0 ) GRIDS += 1;
	const int SIZED = NEQ * sizeof(double);
	const int SIZEW = 11 * SIZED;
	double *yh0; cutilSafeCall(cudaMalloc( (void**) &yh0, SIZED ) );
	double *yhn; cutilSafeCall(cudaMalloc( (void**) &yhn, SIZED ) );
//	double *ese; cutilSafeCall(cudaMalloc( (void**) &ese, SIZED ) );	// esterr
	double *yl0; cutilSafeCall(cudaMalloc( (void**) &yl0, SIZED ) );	// work1
	double *yln; cutilSafeCall(cudaMalloc( (void**) &yln, SIZED ) );	// work2
	double *err; cutilSafeCall(cudaMalloc( (void**) &err, SIZED ) );	// work3
	double *wrk;               cudaMalloc( (void**) &wrk, SIZEW )  ;	// work4
	printf("NDE=%d NEQ=%d NP=%d GRIDS=%d\n", *NDE, NEQ, NP, GRIDS);
	#endif
	SetConst <<< 1, 1 >>> ( NP, eps);

	if( *INIT == 1 )
	{	for(int i = 0; i < NEQ; i++) YL0[i] = YH0[i];
		#ifndef GPU
		moon	( *X0, YL0, WRK);
		#else
		cutilSafeCall( cudaMemcpy	( yh0, YH0, SIZED, cudaMemcpyHostToDevice));
		cutilSafeCall( cudaMemcpy	( yl0, YL0, SIZED, cudaMemcpyHostToDevice));
		func	<<< GRIDS, BLOCKS >>>( *X0, yl0,  wrk); cudaThreadSynchronize() ;	// funcでは初めのNEQ個(SIZEDバイト)に
		               cudaMemcpy	( WRK, wrk, SIZED, cudaMemcpyDeviceToHost) ; // 時間微分を入れる
//		               cudaMemcpy	( WRK, wrk, SIZEW, cudaMemcpyDeviceToHost) ;
		cudaThreadSynchronize();
		#endif
//		printf("NP=%d\n",NP);	printf("NEQ=%d\n",NEQ);		printf("EPS_TH=%e\n",eps);
//		for(int i=0;i<12;i++) printf("WRK[%2d] = %e YL0[%2d] = %e\n", i, WRK[i], i, YL0[i]);
		double TOLDY = std::numeric_limits<double>::max();
		bool changed = false;
		for( int i = 0; i < NEQ; i++)
		{	if( fabs( WRK[i] ) > 0.0 )
			{	double TOL = *RELERR * fabs( YL0[i] ) + *ABSERR;
				TOLDY = min( TOLDY, fabs( TOL / WRK[i] ) );
				changed = true;
			}
		}
		double HMIN = EPSMIN * ( *XE - *X0 );
		if(!changed)	{ *H = HMIN; }
		else				{ *H = min ( ( *XE - *X0 ) / 100e0, pow( TOLDY, 0.2 ) ); }
		*INIT = 2;
		*NDE  = 1;
	}
	// main iteration
	double X = *X0;
	for( int ITER = 0; ITER < ITEMAX; ITER++)
	{	double h = *H;
		#ifndef GPU
		flstep ( NEQ, X, h, YL0, YLN, ERR, WRK);
		#else
		{	int LINDEX = 1;
			//	FEHL 		<<< GRIDS, BLOCKS >>> ( LINDEX, NEQ, X, *H, yl0, yln, err, wrk);
			{	func 		<<< GRIDS, BLOCKS >>> ( X      , yl0       , wrk      ); cudaThreadSynchronize();
				fehl_1 	<<< GRIDS, BLOCKS >>> ( h      , yl0       , wrk      ); cudaThreadSynchronize();
				func 		<<< GRIDS, BLOCKS >>> ( X+h*AL2, wrk+NEQ* 6, wrk+NEQ*1); cudaThreadSynchronize();
				fehl_2 	<<< GRIDS, BLOCKS >>> ( h      , yl0       , wrk      ); cudaThreadSynchronize();
				func 		<<< GRIDS, BLOCKS >>> ( X+h*AL3, wrk+NEQ* 7, wrk+NEQ*2); cudaThreadSynchronize();
				fehl_3 	<<< GRIDS, BLOCKS >>> ( h      , yl0       , wrk      ); cudaThreadSynchronize();
				func 		<<< GRIDS, BLOCKS >>> ( X+h*AL4, wrk+NEQ* 8, wrk+NEQ*3); cudaThreadSynchronize();
				fehl_4 	<<< GRIDS, BLOCKS >>> ( h      , yl0       , wrk      ); cudaThreadSynchronize();
				func 		<<< GRIDS, BLOCKS >>> ( X+h*AL5, wrk+NEQ* 9, wrk+NEQ*4); cudaThreadSynchronize();
				fehl_5 	<<< GRIDS, BLOCKS >>> ( h      , yl0       , wrk      ); cudaThreadSynchronize();
				func 		<<< GRIDS, BLOCKS >>> ( X+h*AL6, wrk+NEQ*10, wrk+NEQ*5); cudaThreadSynchronize();
				fehl_6 	<<< GRIDS, BLOCKS >>> ( h      , yl0       , wrk,yln  ); cudaThreadSynchronize();
				if ( LINDEX == 1 ) fehl_e <<< GRIDS, BLOCKS >>> ( h, err, wrk ); cudaThreadSynchronize();
			}
		}
		cudaMemcpy( YLN, yln, SIZED, cudaMemcpyDeviceToHost);
		cudaMemcpy( ERR, err, SIZED, cudaMemcpyDeviceToHost);
		cudaMemcpy( WRK, wrk, SIZEW, cudaMemcpyDeviceToHost);
		#endif
		*NDE += 6;
		double ERRET = 0;
		for( int i = 0; i < NEQ; i++)
		{	double ET = *RELERR * ( fabs( YL0[i] ) + fabs( YLN[i] ) ) / 2 + *ABSERR;
			if( ET == 0 )
			{	printf("(SUBR -RKF) TROUBLE(INAPPROPRIATE ERROR TOLERANCE)\n");
				*IER = 20000;
				cudaFree(yh0); cudaFree(yhn); cudaFree(yl0); cudaFree(yln); cudaFree(err);cudaFree(wrk);
				exit(0);
			}
			else
			{	ERRET = max( ERRET, fabs( ERR[i] ) / ET );}
		}
		if( ERRET > 1 )
		{	if( ERRET >= 59049 )	{ *H *= 0.1                   ;}
			else							{ *H *= 0.9 / pow( ERRET, 0.2);}
			if( *H <= EPSMIN )
			{	printf("(SUBR -RKF) TROUBLE(TOO SMALL STEP SIZE) H %25.16e ITER %d\n", *H, ITER);
				*IER = 30000;
				cudaFree(yh0); cudaFree(yhn); cudaFree(yl0); cudaFree(yln); cudaFree(err);cudaFree(wrk);
				exit(0);
			}
		}
		else if( X + *H < *XE )
		{	h = *H;
			#ifndef GPU
			fhstep( NEQ, X, h, YH0, YHN, ERR, WRK);
			///*
			#else
			// FHSTEP(const int NEQ, const double X, const double H, double *Y0, double *YN, double *ERR, double *WK)
			{	double X1 = X, H1 = *H / 2;
				{	//	FEHL <<< GRIDS, BLOCKS >>> ( LINDEX, NEQ, X, *H, yh0, yhn, err, wrk);
					func		<<< GRIDS, BLOCKS >>> ( X1       , yh0       , wrk      ); cudaThreadSynchronize();
					fehl_1	<<< GRIDS, BLOCKS >>> ( H1, yh0, wrk)     ; cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL2, wrk+NEQ* 6, wrk+NEQ  ); cudaThreadSynchronize();
					fehl_2	<<< GRIDS, BLOCKS >>> ( H1, yh0, wrk)     ; cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL3, wrk+NEQ* 7, wrk+NEQ*2); cudaThreadSynchronize();
					fehl_3	<<< GRIDS, BLOCKS >>> ( H1, yh0, wrk)     ; cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL4, wrk+NEQ* 8, wrk+NEQ*3); cudaThreadSynchronize();
					fehl_4	<<< GRIDS, BLOCKS >>> ( H1, yh0, wrk)     ; cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL5, wrk+NEQ* 9, wrk+NEQ*4); cudaThreadSynchronize();
					fehl_5	<<< GRIDS, BLOCKS >>> ( H1, yh0, wrk)     ; cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL6, wrk+NEQ*10, wrk+NEQ*5); cudaThreadSynchronize();
					fehl_6	<<< GRIDS, BLOCKS >>> ( H1, yh0, wrk, yhn); cudaThreadSynchronize();
				}
				cudaMemcpy( yh0, yhn, SIZED, cudaMemcpyDeviceToDevice);
				cudaMemcpy( YHN, yhn, SIZED, cudaMemcpyDeviceToHost);
				cudaMemcpy( yh0, YHN, SIZED, cudaMemcpyHostToDevice);
				cudaMemcpy( WRK, wrk, SIZEW, cudaMemcpyDeviceToHost);
				X1 += H1;
				for( int i = 0; i < NEQ; i++){ YH0[i] = YHN[i];}
				cudaMemcpy( yh0, yhn, SIZED, cudaMemcpyDeviceToDevice);
				{	//	FEHL <<< GRIDS, BLOCKS >>> ( LINDEX, NEQ, X, *H, yh0, yhn, err, wrk);
					func		<<< GRIDS, BLOCKS >>> ( X1       , yh0       , wrk      ); cudaThreadSynchronize();
					fehl_1	<<< GRIDS, BLOCKS >>> ( H1,   yh0, wrk                  ); cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL2, wrk+NEQ* 6, wrk+NEQ  ); cudaThreadSynchronize();
					fehl_2	<<< GRIDS, BLOCKS >>> ( H1,   yh0, wrk                  ); cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL3, wrk+NEQ* 7, wrk+NEQ*2); cudaThreadSynchronize();
					fehl_3	<<< GRIDS, BLOCKS >>> ( H1,   yh0, wrk                  ); cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL4, wrk+NEQ* 8, wrk+NEQ*3); cudaThreadSynchronize();
					fehl_4	<<< GRIDS, BLOCKS >>> ( H1,   yh0, wrk                  ); cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL5, wrk+NEQ* 9, wrk+NEQ*4); cudaThreadSynchronize();
					fehl_5	<<< GRIDS, BLOCKS >>> ( H1,   yh0, wrk                  ); cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL6, wrk+NEQ*10, wrk+NEQ*5); cudaThreadSynchronize();
					fehl_6	<<< GRIDS, BLOCKS >>> ( H1,   yh0, wrk       , yhn      ); cudaThreadSynchronize();
				}
				cudaMemcpy( YH0, yh0, SIZED, cudaMemcpyDeviceToHost);
				cudaMemcpy( YHN, yhn, SIZED, cudaMemcpyDeviceToHost);
				cudaMemcpy( WRK, wrk, SIZEW, cudaMemcpyDeviceToHost);
			}
			#endif
			//*/
			*NDE += 12;
			X += *H; *X0 = X;
			for( int i = 0; i < NEQ; i++) { YL0[i] = YLN[i]; YH0[i] = YHN[i]; }
/*
			#ifdef GPU
			cudaMemcpy( yl0, yln, SIZED, cudaMemcpyDeviceToDevice );
			cudaMemcpy( yh0, yhn, SIZED, cudaMemcpyDeviceToDevice );
			#endif
*/
			#ifdef DEBUG
			fprintf(fp, "%24.16e%25.16e%25.16e%25.16e%25.16e%25.16e%25.16e\n",
					 *X0, YH0[0], YH0[1], YH0[2], YH0[3], YH0[4], YH0[5] );
			if (ITER % 10 == 0 )
			{	if ( ITER == 0 ) printf("%s ", "                 ");
				printf("%s ", "*"); fflush(stdout);
				if (ITER % 100 == 0 ) printf("X = %8.6f ITER=%d\n", *X0, ITER); fflush(stdout);
			}
			#endif
			if( ERRET <= 1.889568e-4 )	{ *H *= 5.0;}
			else									{ *H *= 0.9 / pow( ERRET, 0.2);}
		}
		else	// i.e. ( X + *H > *XE )
		{	h = *H = *XE - X;
			#ifndef GPU
			flstep( NEQ, X, h, YL0, YLN, ERR, WRK);
			//*
			#else
			{	int LINDEX = 1;
				//	FEHL <<< GRIDS, BLOCKS >>> ( LINDEX, NEQ, X, *H, yl0, yln, err, wrk);
				func		<<< GRIDS, BLOCKS >>> ( X      , yl0       , wrk      );
				fehl_1	<<< GRIDS, BLOCKS >>> ( h      , yl0       , wrk      );
				func		<<< GRIDS, BLOCKS >>> ( X+h*AL2, wrk+NEQ* 6, wrk+NEQ  );
				fehl_2	<<< GRIDS, BLOCKS >>> ( h      , yl0       , wrk      );
				func		<<< GRIDS, BLOCKS >>> ( X+h*AL3, wrk+NEQ* 7, wrk+NEQ*2);
				fehl_3	<<< GRIDS, BLOCKS >>> ( h      , yl0       , wrk      );
				func		<<< GRIDS, BLOCKS >>> ( X+h*AL4, wrk+NEQ* 8, wrk+NEQ*3);
				fehl_4	<<< GRIDS, BLOCKS >>> ( h      , yl0       , wrk      );
				func		<<< GRIDS, BLOCKS >>> ( X+h*AL5, wrk+NEQ* 9, wrk+NEQ*4);
				fehl_5	<<< GRIDS, BLOCKS >>> ( h      , yl0       , wrk      );
				func		<<< GRIDS, BLOCKS >>> ( X+h*AL6, wrk+NEQ*10, wrk+NEQ*5);
				fehl_6	<<< GRIDS, BLOCKS >>> ( h      , yl0       , wrk,  yln);
				if ( LINDEX == 1 )	fehl_e<<< GRIDS, BLOCKS >>>   (	h, err,  wrk);
			}
			cudaMemcpy( YHN, yhn, SIZED, cudaMemcpyDeviceToHost);
			cudaMemcpy( ERR, err, SIZED, cudaMemcpyDeviceToHost);
			cudaMemcpy( WRK, wrk, SIZEW, cudaMemcpyDeviceToHost);
			#endif
			//*/
			#ifndef GPU
			fhstep( NEQ, X, h, YH0, YHN, ERR, WRK);
			///*
			#else
			// FHSTEP(const int NEQ, const double X, const double H, double *Y0, double *YN, double *ERR, double *WK)
			{	double X1 = X, H1 = *H / 2;
				{	//	FEHL	<<< GRIDS, BLOCKS >>> ( LINDEX, NEQ, X, *H, yh0, yhn, err, wrk);
					func 		<<< GRIDS, BLOCKS >>> ( X1, yh0,   wrk                  )	; cudaThreadSynchronize();
					fehl_1	<<< GRIDS, BLOCKS >>> ( H1, yh0,   wrk                  )	; cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL2, wrk+NEQ* 6, wrk+NEQ  )	; cudaThreadSynchronize();
					fehl_2	<<< GRIDS, BLOCKS >>> ( H1, yh0,   wrk                  )	; cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL3, wrk+NEQ* 7, wrk+NEQ*2)	; cudaThreadSynchronize();
					fehl_3	<<< GRIDS, BLOCKS >>> ( H1, yh0,   wrk                  )	; cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL4, wrk+NEQ* 8, wrk+NEQ*3)	; cudaThreadSynchronize();
					fehl_4	<<< GRIDS, BLOCKS >>> ( H1, yh0,   wrk                  )	; cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL5, wrk+NEQ* 9, wrk+NEQ*4)	; cudaThreadSynchronize();
					fehl_5	<<< GRIDS, BLOCKS >>> ( H1, yh0,   wrk                  )	; cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL6, wrk+NEQ*10, wrk+NEQ*5)	; cudaThreadSynchronize();
					fehl_6	<<< GRIDS, BLOCKS >>> ( H1, yh0,   wrk       , yhn      )	; cudaThreadSynchronize();
				}
				cudaMemcpy( yh0, yhn, SIZED, cudaMemcpyDeviceToDevice);
				cudaMemcpy( YHN, yhn, SIZED, cudaMemcpyDeviceToHost);
				cudaMemcpy( yh0, YHN, SIZED, cudaMemcpyHostToDevice);
				cudaMemcpy( WRK, wrk, SIZEW, cudaMemcpyDeviceToHost);
				X1 += H1;
				for( int i = 0; i < NEQ; i++){ YH0[i] = YHN[i];}
				cudaMemcpy( yh0, yhn, SIZED, cudaMemcpyDeviceToDevice);
				{	//	FEHL <<< GRIDS, BLOCKS >>> ( LINDEX, NEQ, X, *H, yh0, yhn, err, wrk);
					func		<<< GRIDS, BLOCKS >>> ( X1,   yh0, wrk                  )	; cudaThreadSynchronize();
					fehl_1	<<< GRIDS, BLOCKS >>> ( H1,   yh0, wrk                  )	; cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL2, wrk+NEQ* 6, wrk+NEQ  )	; cudaThreadSynchronize();
					fehl_2	<<< GRIDS, BLOCKS >>> ( H1,   yh0, wrk                  )	; cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL3, wrk+NEQ* 7, wrk+NEQ*2)	; cudaThreadSynchronize();
					fehl_3	<<< GRIDS, BLOCKS >>> ( H1,   yh0, wrk                  )	; cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL4, wrk+NEQ* 8, wrk+NEQ*3)	; cudaThreadSynchronize();
					fehl_4	<<< GRIDS, BLOCKS >>> ( H1,   yh0, wrk                  )	; cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL5, wrk+NEQ* 9, wrk+NEQ*4)	; cudaThreadSynchronize();
					fehl_5	<<< GRIDS, BLOCKS >>> ( H1,   yh0, wrk                  )	; cudaThreadSynchronize();
					func		<<< GRIDS, BLOCKS >>> ( X1+H1*AL6, wrk+NEQ*10, wrk+NEQ*5)	; cudaThreadSynchronize();
					fehl_6	<<< GRIDS, BLOCKS >>> ( H1,   yh0, wrk       , yhn      )	; cudaThreadSynchronize();
				}
				cudaMemcpy( YH0, yh0, SIZED, cudaMemcpyDeviceToHost);
				cudaMemcpy( YHN, yhn, SIZED, cudaMemcpyDeviceToHost);
				cudaMemcpy( WRK, wrk, SIZEW, cudaMemcpyDeviceToHost);
			}
			#endif
			//*/
			*NDE += 18;
			//
			for(int i = 0; i < NEQ; i++)	ESE[i] = ( YLN[i] - YHN[i] ) / 31;
			//
			*X0 = *XE;
			for( int i = 0; i < NEQ; i++)
			{	YL0[i] = YLN[i];
				YH0[i] = YHN[i];
			}
			#ifdef GPU
//			cudaMemcpy( yl0, yln, SIZED, cudaMemcpyDeviceToDevice);
//			cudaMemcpy( yh0, yhn, SIZED, cudaMemcpyDeviceToDevice);
			cudaFree(yh0); cudaFree(yhn); cudaFree(yl0); cudaFree(yln); cudaFree(err);cudaFree(wrk);
			#endif
			#ifdef DEBUG
			fprintf(fp,"%24.16e%25.16e%25.16e%25.16e%25.16e%25.16e%25.16e\n",
					 *X0, YH0[0], YH0[1], YH0[2], YH0[3], YH0[4], YH0[5] );
			#endif
//			printf(    "\n%9.6f%20.16f%20.16f%20.16f%20.16f%20.16f%20.16f\n",
//					 *X0, YH0[0], YH0[1], YH0[2], YH0[3], YH0[4], YH0[5] );
			printf("\nNP=%d ITER=%d NDE=%d\n", NP, ITER, *NDE);
			return;
		}
		#ifdef GPU
		cudaMemcpy( yhn, YHN, SIZED, cudaMemcpyHostToDevice);
		cudaMemcpy( yh0, YH0, SIZED, cudaMemcpyHostToDevice);
		cudaMemcpy( yl0, YL0, SIZED, cudaMemcpyHostToDevice);
		cudaMemcpy( yln, YLN, SIZED, cudaMemcpyHostToDevice);
		cudaMemcpy( wrk, WRK, SIZEW, cudaMemcpyHostToDevice);
		#endif
	}	printf("(SUBR -RKF) TROUBLE(TOO MANY ITERATIONS/ ITER =%d)\n", ITEMAX);
	*IER = 10000;
	cudaFree(yh0); cudaFree(yhn); cudaFree(yl0); cudaFree(yln); cudaFree(err);cudaFree(wrk);
	exit(0);
}
