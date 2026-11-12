#include <stdio.h>
#include <cutil_inline.h>
#include <cublas.h>

#define     	BLOCKS	  128	// N = 15,626までOK
//#define     	BLOCKS	64
//#define     	BLOCKS	  30	// GTX295/285では倍精度演算器は30個
#define     	epsilon 2.22e-16
//#define     	BLOCKS	  4
//#define     	BLOCKS	512
//#define DEBUG
//void   SetTimer(unsigned int *t);
//double EndTimer(unsigned int *t);
//
//-----------------------------------------------------------------------------------------------
//グローバル変数 on device
__device__ int N, ITR;
__device__ double eps_th;
__device__ double DT;
__device__ __shared__ double TIME_END_[BLOCKS];
__device__ __shared__ double PARAMS[BLOCKS][5]; //theta_0, sn_0, cs_0, tn_0, g0_sn0; // hyperbola
//
__device__ __constant__ double pi = 3.14159265358979323846264338327950288419716939937510582097494;

#include "bia-rkf.cu"

//__device__ double r_min, g_max, r_min_by_g_max;	// parabola	標準横長放物線⊃
//
__device__ inline void vector ( const double *Y, double *dot, double *crs )
{	*dot = Y[0] * Y[2] + Y[1] * Y[3];  // r dot g	{r} = ( Y[0], Y[1] ), {g} = ( Y[2], Y[3] )
	*crs = Y[0] * Y[3] - Y[1] * Y[2];  // r cross g の ζ成分：
	return;
}
//-----------------------------------------------------------------------------------------------

__device__
double solution_sn_p_cu(double theta, double theta_0)
{	// s_pを求める関数
	return sin( ( theta_0 + theta ) / 2 );
}
__device__
double solution_sn_m_cu(double theta, double theta_0)
{	// s_mを求める関数
	return sin( ( theta_0 - theta ) / 2 );
}
//-----------------------------------------------------------------------------------------------
__device__		// 双曲線
double find_theta_end_cu( double theta, double time_end_, double PARAMS[5], double *error_)
{	double sn_p = solution_sn_p_cu(theta, PARAMS[0]);
	double sn_m = solution_sn_m_cu(theta, PARAMS[0]);
	double cs = cos(theta); double sn = sin(theta);
	double D_inv = 1 / ( cs - PARAMS[2] );
	double dfdth_ = PARAMS[3] * ( 1 - PARAMS[2] * cs ) * D_inv * D_inv + PARAMS[1] / ( 2 * sn_p * sn_m );
	*error_       = PARAMS[3] * sn * D_inv - time_end_ + log( sn_p / sn_m );
	return theta - *error_ / dfdth_;
}
//-----------------------------------------------------------------------------------------------
__device__		// 双曲線	time_hyperbola
double func_time_theta_cu( double theta, double com_dt, double PARAMS[5], double b_0, double g_0)
{	double bbb = 0;
	double sn_p = solution_sn_p_cu( theta, PARAMS[0]);
	double sn_m = solution_sn_m_cu( theta, PARAMS[0]);
	double cs   = cos(theta); double sn = sin(theta);
	double aaa  = PARAMS[3] * sn / ( cs - PARAMS[2] );
	if( sn_p * sn_m < 0 )	bbb = 0;
	else 							bbb = log(sn_p / sn_m);
	if( theta == 0 ) bbb = 0;
	return  b_0 * ( aaa + bbb ) / g_0 - com_dt;
}
//-----------------------------------------------------------------------------------------------
__device__
double solution_energy_cu( double x, double y, double u, double v)
{  //エネルギ
	double E = u * u + v * v + 2 * eps_th / sqrt( x * x + y * y );
	return E;
}
//-----------------------------------------------------------------------------------------------
__device__
double angular_momentum_cu( double x, double y, double u, double v)
{  //角運動量
	return x * v - y * u;
}
//-----------------------------------------------------------------------------------------------
__device__ void hyperbola ( double *YZ, double *YN )
{	//二次元の変化量
	__shared__ double PARAMS[BLOCKS][5]; // 0 = theta_0, 1 = sn_0, 2 = cs_0, 3 = tn_0, 4 = g0_sn0;
	__shared__ double TIME_END_[BLOCKS];
	double W_ini = solution_energy_cu( YZ[0], YZ[1], YZ[2], YZ[3] );	//	ξ-η系
	//
	double X = YZ[0]; double Y = YZ[1]; double U = YZ[2]; double V = YZ[3];	// とりあえず成分の初期位置と速度
	double r_g =       X * U + Y * V  ;		// 内積
	double L_ini  =    X * V - Y * U  ;		// 角運動量
	double r   = sqrt( X * X + Y * Y );		// 距離
	//
	double g_0 =   sqrt ( W_ini );
	double b_0 = eps_th / W_ini  ;
	double b   = fabs( L_ini ) / g_0;
	double theta_0 = PARAMS[threadIdx.x][0] = atan2( fabs(L_ini), b_0 * g_0 );		// [ 0, pi/2]
	double    sn_0 = PARAMS[threadIdx.x][1] = fabs(L_ini) / sqrt( b_0 * b_0 * g_0 * g_0 + L_ini * L_ini );
	double    cs_0 = PARAMS[threadIdx.x][2] = b_0 * g_0   / sqrt( b_0 * b_0 * g_0 * g_0 + L_ini * L_ini );

	double    tn_0 = PARAMS[threadIdx.x][3] = sn_0 / cs_0;
//	double    tn_0 = PARAMS[threadIdx.x][3] = fabs( L_ini ) / ( b_0 * g_0 );



	double theta_ini;
	double cs; double sn;
	if ( fabs( cs_0 + b * sn_0 / r ) > 1 )
	{	theta_ini = 1e-30;
		cs = cos(theta_ini);
	}
	else
	{	theta_ini = acos( cs_0 + b * sn_0 / r );	//とりあえず正で定義
		cs = cs_0 + b * sn_0 / r;
	}
	if( L_ini * r_g < 0 ) theta_ini = - theta_ini;	//これで正しい符号
	sn = sin(theta_ini);

	//if( L_ini < 0 ) printf("/n/n/nNegative L_ini =%g\n\n\n", L_ini) ;

	/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
	double x = r * cs;						//新座標系での初期位置が求まった
	double y = r * sn;
	double g0_sn0;
	if ( L_ini > 0 )
	{										//これで正しい符号
		g0_sn0 = PARAMS[threadIdx.x][4] =   g_0 / sn_0;
	}
	else
	{	g0_sn0 = PARAMS[threadIdx.x][4] = - g_0 / sn_0;
	}
	//
	double u = g0_sn0 * (     cs_0 * sn );	//新座標系での初期速度が求まった
	double v = g0_sn0 * ( 1 - cs_0 * cs );
//	double alpha = atan2( x * Y - y * X, x * X + y * Y );
//	double beta  = atan2( u * V - v * U, u * U + v * V );
//	double cs_a = cos(alpha);
//	double sn_a = sin(alpha);
	double cs_a =( x * X + y * Y ) / sqrt( ( X * X + Y * Y ) * ( x * x + y * y ) );  //空間座標の回転角
	double sn_a =( x * Y - y * X ) / sqrt( ( X * X + Y * Y ) * ( x * x + y * y ) );  //速度座標の回転角
//-----------------------------------------------------------------------------------------------------
	double time_ini = (   tn_0 * sn / ( cs - cs_0 )
						+ log( solution_sn_p_cu(theta_ini, theta_0) / solution_sn_m_cu(theta_ini, theta_0) )
					  ) * b_0 / g_0;
//-----------------------------------------------------------------------------------------------------
//
	double theta_end;
	double    error_;
//
	{	double time_end;
		if ( L_ini > 0 ) {
			time_end = time_ini + DT;// double delta_t =   dt;	//下 y < 0 から上 y > 0 へ。時間変化は正
		} else {
			time_end = time_ini - DT;// double delta_t = - dt;	//上 y > 0 から下 y < 0 へ。時間変化は負
		}
		TIME_END_[threadIdx.x] = time_end * g_0 / b_0;
	}
	{	double            theta_ini_end, error_ini_;
		double theta_gss, theta_gss_end, error_gss_;
		//
		theta_gss = 2 / pi * theta_0 * atan( TIME_END_[threadIdx.x] / tn_0 );
		theta_ini_end = find_theta_end_cu(theta_ini, TIME_END_[threadIdx.x], &PARAMS[threadIdx.x][0], &error_ini_ );
		theta_gss_end = find_theta_end_cu(theta_gss, TIME_END_[threadIdx.x], &PARAMS[threadIdx.x][0], &error_gss_ );
		if ( fabs( error_gss_ ) > 10 * fabs( error_ini_ ) )
		{	//無次元化しないと正確に error_ini  = - Δt = - 1
			theta_end = theta_ini_end;
			error_    = error_ini_;
		}
		else
		{	// 多くの場合、こちら/
			theta_end = theta_gss_end;
			error_    = error_gss_;
		}
	}
//-----------------------------------------------------------------------------------------------
   int ir  = 0;
   //最小のエラーを選択するニュートン法
	double theta_now = 0;
	double err_stock = error_;
	double theta_stock = 0;
	double error_restore_[20];
	int num = 0;
	while( ir < ITR )
	{	theta_now = theta_end;
//		if ( solution_sn_p_cu(theta_now, theta_0) * solution_sn_m_cu(theta_now, theta_0) < 0 )
//		{
//			//対数の中身が負のため、二分法によって解を求める
//			bisectional_method(theta_ini, time_ini, L_ini, b_0, g_0, &error_);
//			break;
//		}
		ir++;
		theta_end = find_theta_end_cu(theta_now, TIME_END_[threadIdx.x], &PARAMS[threadIdx.x][0], &error_ );
		num = ir - 1;
		while( num >= 1 )
		{	if(error_restore_[num - 1] == error_) goto loop;
			num--;
		}
		error_restore_[ir - 1] = error_;
		if( abs(error_) < abs(err_stock) )	//もしエラーの小さい角度がでたら、エラーと角度を更新
		{	err_stock   = error_;				//エラーの小さい法を保存
			theta_stock = theta_end;			//そのときの角度も保存
		}
		if ( abs(err_stock) < epsilon ) break;
	}
	loop:
	theta_end = theta_stock;       //最終的にループ中で最もエラーが小さかった角度を採用
	//-------------------------------------------------------------------------------/
	//DT後の粒子の情報 ( x, y, u, v)
	double cs_end = cos(theta_end);
	double sn_end = sin(theta_end);
	double  u_end = g0_sn0 *       cs_0 * sn_end  ;
	double  v_end = g0_sn0 * ( 1 - cs_0 * cs_end );
	double du_end = u_end - u;
	double dv_end = v_end - v;
	double  r_end = b * sn_0 / ( cs_end - cs_0 );
	double  x_end = r_end * cs_end;
	double  y_end = r_end * sn_end;
	//-------------------------------------------------------------------------------/
	double r1end  =  x_end * cs_a -  y_end * sn_a;   //  X_end    実験室系への変換
	double r2end  =  x_end * sn_a +  y_end * cs_a;   //  Y_end //
	double dg1end = du_end * cs_a - dv_end * sn_a;   // dU_end //
	double dg2end = du_end * sn_a + dv_end * cs_a;   // dV_end //
	double g1end  = YZ[2] + dg1end;
	double g2end  = YZ[3] + dg2end;
	YN[0]  = r1end; YN[1] =  r2end;
	YN[2]  = g1end; YN[3] =  g2end;
}
//
//-----------------------------------------------------------------------------------------------
//
__device__ void newton_ellipse
	(	double dt, double xi_ini, double ecc, double omega_a, double time_end, double tau_end, 
		double *xi_end, double *error )
{	double xi;
	int ir = 0;
	double error_ini = tau_end - xi_ini + ecc * sin(xi_ini) ;
	double delta_ini = error_ini *  ( 1 - ecc * cos(xi_ini));
	if ( error_ini == 0 )
	{	*xi_end = xi_ini;
		return;
	}
	double xi_2 = xi_ini + omega_a * dt;  // initial guess, exact when ecc = 0
	if ( cos(tau_end) < 0.9 ) { xi_2 = tau_end;}
	else
	{	if(sin(tau_end) >= 0 )	{ xi_2 = pi / 2 ;}
		else								{ xi_2 =-pi / 2 ;}
	}
//	double  xi_gss = xi_2;
//	double tau_gss = xi_gss;
//	double error_tm2 = - ( xi_2  - ecc * sin(xi_2) - omega_a * time_end );
	double error_tm2 =  tau_end - xi_2 + ecc * sin(xi_2);
	double delta_xi2 = error_tm2 / ( 1 - ecc * cos(xi_2) );
	if ( abs(error_ini) < abs(error_tm2) ) {xi = xi_2 + delta_ini;}
	else                                   {xi = xi_2 + delta_xi2;}
	/*!!!!!!!自分で予想!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
	double xi_1 = xi ;
	/*!!!!!!!自分で予想!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
	//	while ( xi_1 != xi_2 && ir <= ITR ){
	double err_stock   = 100;
	double theta_stock = 0;
	double error_time  = 0;
	double delta_xi1   = 0;
	while ( ir <= ITR )
	{	ir++;
		xi_1 = xi;
		error_time  = tau_end - xi_1 + ecc * sin(xi_1)  ;
		delta_xi1 = error_time / ( 1 - ecc * cos(xi_1) );
		xi = xi_1 + delta_xi1;
		if(abs(error_time) < abs(err_stock))	//もしエラーの小さい角度がでたら、エラーと角度を更新
		{	err_stock   = error_time;				//エラーの小さい法を保存
			theta_stock = xi;							//そのときの角度も保存
		}
		if ( abs(err_stock) < 2e-16 ) break;
	}
	*xi_end = theta_stock;       //最終的にループ中で最もエラーが小さかった角度を採用
	return;
}

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

__device__ void ellipse ( double *YZ, double *YN )
{	//二次元の変化量
	double error; double xi_end;
	double X_ini = YZ[0]; double Y_ini = YZ[1]; double U_ini = YZ[2]; double V_ini = YZ[3];
//	double dt = DT;
	double g2_ini =       U_ini*U_ini + V_ini*V_ini  ;	// 相対速度の自乗
	double r_ini  = sqrt( X_ini*X_ini + Y_ini*Y_ini );	// 相対距離
	double W      = g2_ini + 2 * eps_th / r_ini;      	// W = eps_th / a = - g_a**2
//
	double L       = X_ini * V_ini - Y_ini * U_ini;   	// 角運動量
	double r_g_ini = X_ini * U_ini + Y_ini * V_ini;   	// r cdot g
//
//	double g_a = sign( sqrt( - W ), L );              	// 反時計回りが正：角運動量と同符号
	double g_a;
	if( L >= 0 )	{ g_a =   sqrt( - W );}               	// 反時計回りが正：角運動量と同符号
	else      		{ g_a = - sqrt( - W );}               	// 反時計回りが正：角運動量と同符号
	double a  =  eps_th / W;                       		 	// 長半径：正で定義
	double b  =  L / g_a;                             	// 短半径：正で定義
	double omega_a = g_a / a;                         	// g_a / a
	double xi_ecc = asin( L / ( a * g_a ) ); double ecc = cos(xi_ecc);
//	double ecc = sqrt( L * L - a * ,a * g_a * g_a ) / L;
//  -----------------------------------------------------------------------------------
	double xi_ini = 0;
	if ( r_g_ini == 0 )
	{	// 初期条件が近日点または遠日点でない場合
		if			( r_ini <  a ) { xi_ini = 0; }						// 近日点
//		else if ( r_ini >  a ) { xi_ini = sign ( pi, L );}	// 遠日点
		else
		{	if		( r_ini == a )
			{	xi_ini = atan2( Y_ini, X_ini );
				ecc    =      0;
				xi_ecc = pi / 2;
			}
			else
			{	if ( L >= 0 )	{ xi_ini =   pi; }
				else				{ xi_ini = - pi; }
			}
		}                  // 遠日点
	}
   else               	{ xi_ini = atan2(r_g_ini / g_a, a - r_ini); } // 円運動
//	-------------------------------------------------------------------------------------------
//	theta = atan2( b * sin(xi_ini), a * ( cos(xi_ini) - ecc ) )                                
	double alpha = atan2( Y_ini, X_ini ) - atan2( b * sin(xi_ini), a * ( cos(xi_ini) - ecc ) );
//  ------------------------------------------------------------------------------------------
	double  tau_ini = xi_ini - ecc * sin(xi_ini);
	double time_ini = tau_ini / omega_a;
	double time_end = time_ini + DT;
	double  tau_end = omega_a * time_end;
	if ( DT == 0 )
	{	error = 0;
		YN = YZ;
		return;
	}
	//ここでニュートン法を使ってxi_endをもとめる
	newton_ellipse( DT, xi_ini, ecc, omega_a, time_end, tau_end, &xi_end, &error );
	//
	double cs_end = cos( xi_end );
	double sn_end = sin( xi_end );
	double xx_end =   a * ( cs_end - ecc );
	double yy_end =   b *   sn_end;
	double rr_end =   a * ( 1 - ecc * cs_end );
	double uu_end = - a * g_a * sn_end / rr_end;
	double vv_end =   b * g_a * cs_end / rr_end;
	//
	double cs_alp = cos(alpha);
	double sn_alp = sin(alpha);
	YN[0] =   xx_end * cs_alp - yy_end * sn_alp;
	YN[1] =   xx_end * sn_alp + yy_end * cs_alp;
	YN[2] =   uu_end * cs_alp - vv_end * sn_alp;
	YN[3] =   uu_end * sn_alp + vv_end * cs_alp;
	return;
//	return omega_a;
}
//-----------------------------------------------------------------------------------------------
//
__device__ inline double time_parabola	( const double theta, const double r_min_by_g_max )  // t = t(θ)
{	double cs, sn; sincos( theta, &sn, &cs );
	return ( 2 + cs ) * sn / ( 1 + cs ) / ( 1 + cs ) * 4 * r_min_by_g_max / 3;
//	return ( 2 + cs ) * sn / ( 1 + cs ) / ( 1 + cs ) * 4 * r_min /  g_max / 3;
}

__device__ inline double dt_dtheta 		( const double theta, const double r_min_by_g_max )
{	double tmp = 1 / cos( theta / 2 );
	tmp = tmp * tmp;	// sec^2 beta
	tmp = tmp * tmp;	// sec^4 beta
	return r_min_by_g_max * tmp;
//	return r_min  / g_max * tmp;
}

__device__ void parabola	// (xi,eta)系での位置・速度を横長楕円の(x,y)系に変換して時間DTだけ進める
	(	const	double *YZ,		// ( xi, eta, g_xi, g_eta)_ini
				double *YN			// ( xi, eta, g_xi, g_eta)_end
	)
{	double r_ini = sqrt( YZ[0] * YZ[0] + YZ[1] * YZ[1] );	// |r|
	double r_dot_g, r_crs_g; vector ( &YZ[0], &r_dot_g, &r_crs_g ); //{r}と{g}の内積・外積の計算
//	
//	最短距離 r_min at (x,y) = ( r_min, 0) で最速 g_max: ( g_x, g_y ) = ( 0, g_max )
	const double g_max = - 2 * eps_th / r_crs_g; // g_max * L = - eps > 0
	const double r_min =      r_crs_g /   g_max; // L = r_min * g_max > 0
//	const double r_min = -    r_crs_g * r_crs_g / ( 2 * eps_th );
	const double r_min_by_g_max = r_min / g_max;

	double x_ini = 2 * r_min - r_ini; double tmp = 2 * r_min * sqrt( r_ini / r_min - 1 );
	double y_ini = ( r_dot_g > 0 )  ? tmp : - tmp;	// y_ini = y_ini と r_dot_g は同符号
//	double x_ini =           2 * r_min -       r_ini                         ;
//	double y_ini = copysign( 2 * r_min * sqrt( r_ini / r_min - 1 ), r_dot_g );	// OK?

	double theta_ini = atan2( y_ini, x_ini );
	double  time_ini = time_parabola ( theta_ini, r_min_by_g_max );
	double  time_end = time_ini + DT;
	double  time_end_= 3 * g_max / ( 4 * r_min ) * time_end;	// x^3 + 3 x + 2 time_end_ = 0
	double Fp = cbrt( sqrt( time_end_ * time_end_ + 1 ) + time_end_ );
	double Fm = cbrt( sqrt( time_end_ * time_end_ + 1 ) - time_end_ );
	double theta_gss = 2 * atan( Fp - Fm );	// Newton法の初期値をCardanoの公式で与える

/*
	double theta_end;
	int ierr = newton_parabola ( time_end, theta_gss, &theta_end ); printf("ierr = %d\n",ierr);
*/
	//	begin Newton法
	const int ITR = 20;
	double theta[ITR], err_t[ITR];
	double theta_end, err_t_min;
	theta[0] = theta_gss;
	err_t[0] = time_parabola(theta[0], r_min_by_g_max ) - time_end;
	//printf("\n%2d%24.16e%24.16e%2d\n", 0, theta_gss, err_t[0], 0 );
	int i_min = 0; err_t_min = fabs(err_t[0]);
	for ( int i = 1; i < ITR; i++ )
	{	if( err_t[i-1] == 0 ) { i_min = i-1; break; }
		theta[i] = theta[i-1] - err_t[i-1] / dt_dtheta( theta[i-1], r_min_by_g_max );
		err_t[i] = time_parabola( theta[i], r_min_by_g_max ) - time_end;
		if( theta[i] == theta[i-1] )
		{	i_min = i; break; }
		else
		{	if(err_t_min >	fabs( err_t[i] ) )
			{	i_min = i;
				err_t_min =	fabs( err_t[i] ) ;
			}
			for ( int j = 0; j < i; j++)
			{	if(err_t[i] == err_t[j] )
				{	err_t_min = fabs(err_t[i_min]);
					i_min = j;
					break;	// 二重ループを抜けられないので、後でもう一度breakしていることに注意。
				}
			}
		}
		//printf("%2d%24.16e%24.16e%2d\n", i, theta[i], err_t[i], i_min );
		if( i != i_min && err_t[i] == err_t[i_min] ) break;
	}	theta_end = theta[i_min] ;
	//	end Newton法
//	*/

	double beta_end = theta_end / 2;
	double sn_beta_end, cs_beta_end	; sincos(  beta_end,	&sn_beta_end,	&cs_beta_end );
	double r_end =   r_min / ( cs_beta_end * cs_beta_end );
	double g_end =   g_max *   cs_beta_end;
	double sn_theta_end, cs_theta_end; sincos( theta_end,	&sn_theta_end,	&cs_theta_end);
	double x_end =   r_end * cs_theta_end;
	double y_end =   r_end * sn_theta_end;
	double u_end = - g_end * sn_beta_end;
	double v_end =   g_end * cs_beta_end;
	//
	double gamma = atan2( x_ini * YZ[1] - y_ini * YZ[0], x_ini * YZ[0] + y_ini * YZ[1] );
	double sn_gamma, cs_gamma; sincos( gamma, &sn_gamma, &cs_gamma);
	//
	YN[0] = x_end * cs_gamma - y_end * sn_gamma;  // r_xi_end
	YN[1] = x_end * sn_gamma + y_end * cs_gamma;  // r_eta_end
	YN[2] = u_end * cs_gamma - v_end * sn_gamma;  // g_xi_end
	YN[3] = u_end * sn_gamma + v_end * cs_gamma;  // g_eta_end
	//
	//*error = time_end - time_parabola(theta_end);
	//
	return;
}
//-----------------------------------------------------------------------------------------------
__device__	inline //	三次元(X,Y,Z)から二次元(ξ,η)への変換、および単位ベクトルの計算
void xyz_to_xi_eta_zeta( double *R, double *G, double *e_xi, double *e_eta, double *YZ)
{	double gm = sqrt( G[0] * G[0] + G[1] * G[1] + G[2] * G[2] );
	e_eta [0] = G[0] / gm;
	e_eta [1] = G[1] / gm;
	e_eta [2] = G[2] / gm;
	{	double L[3];
		L[0] = R[1] * G[2] - R[2] * G[1];
		L[1] = R[2] * G[0] - R[0] * G[2];
		L[2] = R[0] * G[1] - R[1] * G[0];
		double Lm = sqrt( L[0] * L[0] + L[1] * L[1] + L[2] * L[2] );
		double e_zeta[3];
		e_zeta[0] = L[0] / Lm;
		e_zeta[1] = L[1] / Lm;
		e_zeta[2] = L[2] / Lm;
		e_xi  [0] = e_eta[1] * e_zeta[2] - e_eta[2] * e_zeta[1];
		e_xi  [1] = e_eta[2] * e_zeta[0] - e_eta[0] * e_zeta[2];
		e_xi  [2] = e_eta[0] * e_zeta[1] - e_eta[1] * e_zeta[0];
	}
	YZ[0] = R[0] * e_xi [0] + R[1] * e_xi [1] + R[2] * e_xi [2]; YZ[2] =  0;
	YZ[1] = R[0] * e_eta[0] + R[1] * e_eta[1] + R[2] * e_eta[2]; YZ[3] = gm;
}

//-----------------------------------------------------------------------------------------------

__global__			//並列にニュートン法の開始: メインカーネルコード
void BIA_Device ( const int METHOD, const int i, double *RXYZUVW, double *dXYZUVW)	// 粒子 i 固定
{	unsigned int j = blockIdx.x * blockDim.x + threadIdx.x;											// 粒子 j に関して並列
	if( i < j && j < N )
	{	double e_xi[3], e_eta[3], YZ[4], YN[4];
		//	----------------------------------------------------------------------------------------------------
		//	( xi, eta, zeta )系への変換 ( R[0:2], G[0:2]) → ( YZ[0:1], 0, YZ[2:3], 0 ) // YZ: 相互作用前
		xyz_to_xi_eta_zeta( &RXYZUVW[6*j], &RXYZUVW[6*j+3], &e_xi[0], &e_eta[0], &YZ[0]);
		//	----------------------------------------------------------------------------------------------------
		if ( METHOD == 0 )			rkf_device	( &YZ[0], &YN[0] );
		else
		{	double W = sqrt( YZ[0] * YZ[0] + YZ[1] * YZ[1] ) * ( YZ[2] * YZ[2] + YZ[3] * YZ[3] ) + 2 * eps_th;
			if 		( W >  0 )		hyperbola	( &YZ[0], &YN[0] );	//	双曲線
			else if	( W == 0 )		parabola	( &YZ[0], &YN[0] );	//	放物線
			else	/*	( W <  0 )	*/	ellipse		( &YZ[0], &YN[0] );	//	楕　円
		}	double dr_xi, dr_eta, dg_xi, dg_eta;
		dr_xi = YN[0] - YZ[0]; dr_eta = YN[1] - YZ[1];
		dg_xi = YN[2] - YZ[2]; dg_eta = YN[3] - YZ[3];
		dXYZUVW[6*j  ] = ( dr_xi * e_xi[0] + dr_eta * e_eta[0] - RXYZUVW[6*j+3] * DT ) / 2; // δXij
		dXYZUVW[6*j+1] = ( dr_xi * e_xi[1] + dr_eta * e_eta[1] - RXYZUVW[6*j+4] * DT ) / 2; // δYij
		dXYZUVW[6*j+2] = ( dr_xi * e_xi[2] + dr_eta * e_eta[2] - RXYZUVW[6*j+5] * DT ) / 2; // δZij
		dXYZUVW[6*j+3] = ( dg_xi * e_xi[0] + dg_eta * e_eta[0]                       ) / 2; // ΔUij
		dXYZUVW[6*j+4] = ( dg_xi * e_xi[1] + dg_eta * e_eta[1]                       ) / 2; // ΔVij
		dXYZUVW[6*j+5] = ( dg_xi * e_xi[2] + dg_eta * e_eta[2]                       ) / 2; // ΔWij
	}
	__syncthreads();
}
//
//-----------------------------------------------------------------------------------------------------
__global__ void SetConst	( const int N_p, const int ITR_, const double eps, const double dt)
{	unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if ( tid == 0 )	{ N = N_p; ITR = ITR_; DT = dt; eps_th = eps;}
}
//-----------------------------------------------------------------------------------------------------
extern "C" void find_dxyuv_kernel_
	( const int METHOD, const int N, const int I, const double EPS, const double DT, const int ITR,
	  double *dX, double *dY, double *dZ, double *dU, double *dV, double *dW,
	  double *RXYZUVW, double *dXYZUVW
	)
{	int i = I - 1;  // I in FORTRAN i in C/CUDA
	double *RXYZUVW_cu; cudaMalloc( (void **) &RXYZUVW_cu, 6 * N * sizeof( double) );
	double *dXYZUVW_cu; cudaMalloc( (void **) &dXYZUVW_cu, 6 * N * sizeof( double) );
	cudaMemcpy( RXYZUVW_cu, RXYZUVW, 6 * N * sizeof(double), cudaMemcpyHostToDevice);
	//
	int GRIDS = ( N % BLOCKS == 0 ) ? N / BLOCKS : N / BLOCKS + 1;
	if ( I == 1 ) printf("Started GPU with BLOCKS=%d, GRIDS =%d\n", BLOCKS, GRIDS);
	SetConst	<<<     1,      1 >>>	( N, ITR, EPS, DT);
	BIA_Device	<<< GRIDS, BLOCKS >>>	( METHOD, i, RXYZUVW_cu, dXYZUVW_cu);
	cudaMemcpy( dXYZUVW, dXYZUVW_cu, 6 * N * sizeof(double), cudaMemcpyDeviceToHost);
	for( int j = i + 1; j < N; j++)
	{	dX[i] += dXYZUVW[6*j  ]; dY[i] += dXYZUVW[6*j+1]; dZ[i] += dXYZUVW[6*j+2];
		dX[j] -= dXYZUVW[6*j  ]; dY[j] -= dXYZUVW[6*j+1]; dZ[j] -= dXYZUVW[6*j+2];
		dU[i] += dXYZUVW[6*j+3]; dV[i] += dXYZUVW[6*j+4]; dW[i] += dXYZUVW[6*j+5];
		dU[j] -= dXYZUVW[6*j+3]; dV[j] -= dXYZUVW[6*j+4]; dW[j] -= dXYZUVW[6*j+5];
	}
	/*
	if ( i == 0 ) printf("copysign(2,-1)=%g\n", copysign(2.0,-1.0));
	if ( I == N )
	{	double *ONE = (double*) malloc( sizeof(double) ); ONE[0] = 1;
		double *one; cudaMalloc( (void **) &one, sizeof(double) );
		cudaMemcpy( one, ONE, sizeof(double), cudaMemcpyHostToDevice);
		printf("\nDdot dUa[0] =%25.16e\n", cublasDdot ( 6 * N, one, 0, &dXYZUVW_cu[3], 6 ) );
		printf(  "CUDA dUa[0] =%25.16e\n",                                      dU[0]      );
		cudaFree(one); free(ONE);
	}
	*/
	//メモリ開放
	cudaFree  ( RXYZUVW_cu );
	cudaFree  ( dXYZUVW_cu );
}

