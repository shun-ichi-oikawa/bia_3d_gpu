#pragma once

__device__	__constant__ int    ITEMAX   = 100000;
__device__	__constant__ double MINLIMIT = 1e-300, EPSMIN = 1e-15;
__device__	__constant__ double RELERR = 0, ABSERR = 1e-16,
	AL2d =    1e0 /    4, AL3d =     3e0 /     8, AL4d =    12e0 /    13, AL5d = 1e0, AL6d = 1e0 / 2,
	B21d =    1e0 /    4, B31d =     3e0 /    32, B32d =     9e0 /    32,
	B41d = 1932e0 / 2197, B42d = -7200e0 /  2197, B43d =  7296e0 /  2197,
	B51d =  439e0 /  216, B52d =    -8e0        , B53d =  3680e0 /   513, B54d = -845e0 / 4104,
	B61d =   -8e0 /   27, B62d =     2e0        , B63d = -3544e0 /  2565, B64d = 1859e0 / 4104, B65d = -11e0 / 40,
	GA1d =   16e0 /  135, GA3d =  6656e0 / 12825, GA4d = 28561e0 / 56430, GA5d =   -9e0 /   50, GA6d =   2e0 / 55,
	DA1d =    1e0 /  360, DA3d =  -128e0 /  4275, DA4d = -2197e0 / 75240, DA5d =    1e0 /   50, DA6d =   2e0 / 55;

__device__	//	dY/dX = F(X,Y)
void FUNC	( const double X, const double *Y, double *F )
{	double rsq = Y[0] * Y[0] + Y[1] * Y[1];
	double aaa = eps_th / ( rsq * sqrt(rsq) );
	F[0] = Y[2]      ;	F[1] = Y[3]      ;
	F[2] = Y[0] * aaa;	F[3] = Y[1] * aaa;
}

__device__ void FEHL( const int LINDEX, const double X, const double H, double *Y0, double *YN, double *ER, double *WK )
{	FUNC(X       ,Y0   ,WK   );for(int i=0;i<4;i++){WK[24+i]=Y0[i]+H*(B21d*WK[i]                                                                    );}
	FUNC(X+AL2d*H,WK+24,WK+ 4);for(int i=0;i<4;i++){WK[28+i]=Y0[i]+H*(B31d*WK[i]+B32d*WK[4+i]                                                       );}
	FUNC(X+AL3d*H,WK+28,WK+ 8);for(int i=0;i<4;i++){WK[32+i]=Y0[i]+H*(B41d*WK[i]+B42d*WK[4+i]+B43d*WK[8+i]                                          );}
	FUNC(X+AL4d*H,WK+32,WK+12);for(int i=0;i<4;i++){WK[36+i]=Y0[i]+H*(B51d*WK[i]+B52d*WK[4+i]+B53d*WK[8+i]+B54d*WK[12+i]                            );}
	FUNC(X+AL5d*H,WK+36,WK+16);for(int i=0;i<4;i++){WK[40+i]=Y0[i]+H*(B61d*WK[i]+B62d*WK[4+i]+B63d*WK[8+i]+B64d*WK[12+i]+B65d*WK[16+i]              );}
	FUNC(X+AL6d*H,WK+40,WK+20);for(int i=0;i<4;i++){YN[   i]=Y0[i]+H*(GA1d*WK[i]             +GA3d*WK[8+i]+GA4d*WK[12+i]+GA5d*WK[16+i]+GA6d*WK[20+i]);}
	if  ( LINDEX == 1        ) for(int i=0;i<4;i++){ER[   i]=      H*(DA1d*WK[i]             +DA3d*WK[8+i]+DA4d*WK[12+i]+DA5d*WK[16+i]+DA6d*WK[20+i]);}
}

__device__
void FLSTEP ( const double X, const double H, double *Y0, double *YN, double *ER, double *WK)
{	int LINDEX = 1;
	FEHL( LINDEX, X, H, Y0, YN, ER, WK);
}

__device__
void FHSTEP ( const double X, const double H, double *Y0, double *YN, double *ER, double *WK)
{	double X1 = X, H1 = H / 2; int LINDEX = 0;
	FEHL( LINDEX, X1, H1, Y0, YN, ER, WK); X1 += H1; for( int i = 0; i < 4; i++ ) Y0[i] = YN[i];
	FEHL( LINDEX, X1, H1, Y0, YN, ER, WK);
}

__device__
//void rkf_device( double *YH0, double *YHN )
void rkf_device( const double *YHZ, double *YHN )
{	double YH0[4] = { YHZ[0], YHZ[1], YHZ[2], YHZ[3]};
	double YL0[4], YLN[4], ERR[4], WK[44];
	//ier = 0;
	for(int i = 0; i < 4; i++) YL0[i] = YH0[i];
	double X = 0;
	FUNC( X, YL0, WK );	//int nde += 1;
	double TOLDY = 1e+75;
	bool changed = false;
	for( int i = 0; i < 4; i++)
	{	if( fabs( WK[i] ) > MINLIMIT )
		{	double TOL = RELERR * fabs( YL0[i] ) + ABSERR;
			TOLDY = min( TOLDY, fabs( TOL / WK[i] ) );
			changed = true;
		}
	}
	double H;
	double HMIN = EPSMIN * DT;
	if( !changed )	{ H =                            HMIN ; }
	else					{ H = min( DT / 100, pow(TOLDY, 0.2) ); }
	//	main iteration
	for(int ITER = 0; ITER < ITEMAX; ITER++)
	{	FLSTEP( X, H, YL0, YLN, ERR, WK);	//nde += 12;
		double ERRET = 0;
		for(int i = 0; i < 4; i++)
		{	double ET = RELERR * ( fabs(YL0[i]) + fabs(YLN[i] ) ) / 2 + ABSERR;
			// if ( ET  == 0 ) 	{ ier = 20000; break; }
			if( ET < MINLIMIT )	{ ITER  = ITEMAX; /*err = 1;*/ }
			else							{ ERRET = max(ERRET, fabs(ERR[i]) / ET);}
		}
		if	  ( ERRET >      1 )	// not successful, reduce timestep H
		{	if( ERRET >= 59049 ) { H *= 0.1                  ;}
			else							{ H *= 0.9 / pow(ERRET, 0.2);}
			if( H <= EPSMIN )
			{	ITER = ITEMAX;
				// ier = 30000; break;
			}
		}
		else if ( X + H < DT )	// successful, next integration phase
		{	FHSTEP( X, H, YH0, YHN, ERR, WK); //nde += 12;
			X += H;
			for(int i = 0; i < 4; i++) { YL0[i] = YLN[i]; YH0[i] = YHN[i];	}
			if( ERRET <= 1.889568e-4 )	 { H *= 5.0                  ;}
			else									 { H *= 0.9 / pow(ERRET, 0.2);}
		}
		else								// successful, last integration phase
		{	H = DT - X;
			FLSTEP( X, H, YL0, YLN, ERR, WK);
			FHSTEP( X, H, YH0, YHN, ERR, WK); //nde += 18;
			//for( int i = 0; i < 4; i++ ) esterr[i] = ( YLN[i] - YHN[i] ) / 31;
			for( int i = 0; i < 4; i++ ) { YL0[i] = YLN[i]; YH0[i] = YHN[i]; }
			ITER = ITEMAX;
			//ier = 0; break;
		}
	}
	__syncthreads();
}
