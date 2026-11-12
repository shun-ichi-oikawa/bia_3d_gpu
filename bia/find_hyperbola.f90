subroutine find_hyperbola( eps_th, dt, YZ, YN, error_end, ERR, ITR, IR, IER )
   IMPLICIT REAL*8( A-H, K, L, O-Z )	! YZ, YNはいずれも相対座標・速度 eps_thは運動方程式の二倍
   real(8), dimension(4)       :: YN, YZ
   real(8), dimension(ITR - 1) :: error_restore_
   common / com_time_theta / b_0, g_0, theta_0, tn_0, cs_0, sn_0, time_end
   data pi / 3.14159265358979323846264338327950288419716939937510582097494d0 /

   X_ini = YZ(1); Y_ini = YZ(2); U_ini = YZ(3); V_ini = YZ(4)  ! (X,Y,U,V)系での相対座標・速度
   r_g = X_ini * U_ini + Y_ini * V_ini       ! r_ cdot  g_
   Lz  = X_ini * V_ini - Y_ini * U_ini       ! r_ times g_
   r   = sqrt( X_ini**2 + Y_ini**2 )
   g2  =       U_ini**2 + V_ini**2  	! g^2
   W   = g2 + 2 * eps_th / r     		! energy
!  W   = g2 + 2 * dprod(eps_th, 1/r)	! energy
   if ( W <= 0 ) then
      print *, 'W = ', W
      stop
   else
      g_0 =    sqrt( W )
      b_0 = eps_th / W  
      b   = abs( Lz ) / g_0
   end if
!
   theta_0 = atan2 ( abs(Lz), b_0 * g_0 ) 		! [ 0, pi/2]
!  theta_0 = atan2 ( abs(Lz), dprod(b_0,g_0) ) 	! [ 0, pi/2] ここだけ丁寧に計算
!  theta_0 = theta_0_q
!  cs_0 = cos(theta_0); sn_0 = sin(theta_0);
   cs_0 = b_0 * g_0 / sqrt( b_0**2 * g_0**2 + Lz**2 );
   sn_0 = abs(Lz)   / sqrt( b_0**2 * g_0**2 + Lz**2 );
   tn_0 = abs(Lz)   /     ( b_0    * g_0            ) 
!  cs_0 = eps_th / sqrt(eps_th**2+Lz**2*W); sn_0 = abs(Lz) / sqrt((eps_th**2 / W)+Lz**2);
!  tn_0 = abs(Lz)*sqrt(W) / eps_th
!  tn_0 = sn_0 / cs_0
!  cs_0 = sin( pi / 2 - theta_0 ); sn_0 = sin(theta_0); tn_0 = sn_0 / cs_0
   if ( abs( cs_0 + b * sn_0 / r ) > 1 ) then
      theta_ini = 1d-30
   else
      theta_ini = acos( cs_0 +        b * sn_0 / r )   ! とりあえず正で定義
!     theta_ini = acos( cs_0 + dprod( b , sn_0 / r ) )	! とりあえず正で定義
   end if
   if ( Lz * r_g < 0 ) then                  ! これで正しい符号
      theta_ini = - theta_ini
   end if
!
!	cs = cos(theta_ini); sn = sin(theta_ini)
	call sincos( theta_ini, sn, cs )	! sincos.c 利用
   x = r * cs; !  x = r * cs_0 + b * sn_0 より
   y = r * sn
   if ( Lz > 0 ) then                  ! これで正しい符号
      g0_sn0 =   g_0 / sn_0
   else
      g0_sn0 = - g_0 / sn_0
   end if
   u = g0_sn0 * (     cs_0 * sn )
   v = g0_sn0 * ( 1 - cs_0 * cs )
!
   alpha = atan2 ( x * Y_ini - y * X_ini, x * X_ini + y * Y_ini )
   cs_a = ( x * X_ini + y * Y_ini ) / sqrt( ( X_ini**2 + Y_ini**2 ) * ( x**2 + y**2 ) )
   sn_a = ( x * Y_ini - y * X_ini ) / sqrt( ( X_ini**2 + Y_ini**2 ) * ( x**2 + y**2 ) )
!  beta  = atan2( u * V_ini - v * U_ini, u * U_ini + v * V_ini ); sn_a = sin(alpha)
!  -------------------------------------------------------------------------------------
   time_ini = (    tn_0 * sn / ( cs - cs_0 )                                &
                 + log (   sin( ( theta_0 + theta_ini ) / 2 )               &
                         / sin( ( theta_0 - theta_ini ) / 2 ) ) ) * b_0 / g_0
!  -------------------------------------------------------------------------------------
!
   if ( Lz > 0 ) then; time_end = time_ini + dt; delta_t =   dt   ! 下y<0から上y>0へ。時間変化は正
                 else; time_end = time_ini - dt; delta_t = - dt   ! 上y>0から下y<0へ。時間変化は負
   end if
   time_end_ = time_end * g_0 / b_0
   theta_gss = 2 / pi * theta_0 * atan( time_end_ / tn_0 )	! atan で近似
!
   IR = -1
   call find_theta( theta_ini, theta_ini_end, dfdth_ini_, error_ini_ );!print *, 'ini: ', theta_ini/pi, theta_ini_end/pi, error_ini_ * b_0 / g_0
   IR =  0
   call find_theta( theta_gss, theta_gss_end, dfdth_gss_, error_gss_ );!print *, 'gss: ', theta_gss/pi, theta_gss_end/pi, error_gss_ * b_0 / g_0
   error_gss = error_gss_ * b_0 / g_0
   error_ini = error_ini_ * b_0 / g_0
   if ( abs(error_gss_) > 10 * abs(error_ini_) ) then	! 無次元化しないと正確に	error_ini  = - Δt = - 1
!     theta_now = theta_ini                        		! 無次元化量で          	error_ini_ = - b_0 / g_0
      theta_end = theta_ini_end
      error_    = error_ini_
!     dfdth_    = dfdth_ini_
   else                                   ! 多くの場合、こちら
!     theta_now = theta_gss
      theta_end = theta_gss_end
      error_    = error_gss_
!     dfdth_    = dfdth_gss_
   end if
!  IR = 0; theta_now = theta_end; call find_theta( theta_now, theta_end, dfdth_, error_ )
!
!  -----------------------------------------------------------------------------------------
   err_stock = error_
   loop : do while( IR < ITR )
      theta_now = theta_end
      sn_p = sin( ( theta_0 + theta_now ) / 2 )
      sn_m = sin( ( theta_0 - theta_now ) / 2 )
      if ( sn_p * sn_m < 0 ) then
         IR = -1
         return
      end if
      call find_theta( theta_now, theta_end, dfdth_, error_ ); IR = IR + 1
		 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!エラーがでたら確認
      num = IR - 1
      do while (num >= 1)
        if(error_restore_(num) == error_)then
      		exit loop
        end if
        num = num - 1
      end do
      error_restore_(IR) = error_
		 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(abs(error_) < abs(err_stock)) then
         err_stock = error_* b_0 / g_0
         theta_stock = theta_end
        if ( abs(err_stock) < 2d-16 ) exit
      end if
   end do loop
   theta_end = theta_stock
   error_end = error_ * b_0 / g_0
!  pause
!
!  -----------------------------------------------------------------------------------------
!
!  cs_end = cos(theta_end); sn_end = sin(theta_end)
	call sincos( theta_end, sn_end, cs_end )	! sincos.c 利用
   u_end = g0_sn0 *       cs_0 * sn_end	; du_end = u_end - u
   v_end = g0_sn0 * ( 1 - cs_0 * cs_end )	; dv_end = v_end - v
   r_end = b * sn_0 / ( cs_end - cs_0 )
!
!  -------------------------------------------------------------------------------------------------------
!
   r1end =  r_end * ( cs_end * cs_a - sn_end * sn_a ) !  X_end
   r2end =  r_end * ( cs_end * sn_a + sn_end * cs_a ) !  Y_end

   dg1end = du_end * cs_a - dv_end * sn_a ! dU_end
   dg2end = du_end * sn_a + dv_end * cs_a ! dV_end
!
   g1end = YZ(3) + dg1end
   g2end = YZ(4) + dg2end
   YN(1) = r1end; YN(2) = r2end
   YN(3) = g1end; YN(4) = g2end
!
   return
!
!  -------------------------------------------------------------------------------------
   contains
!
   subroutine find_theta( theta_now, theta_end, dfdth_, error_)
      implicit real(8) (a-h,k-z)
      common /comij/ i,j
      sn_p = sin( ( theta_0 + theta_now ) / 2 )
      sn_m = sin( ( theta_0 - theta_now ) / 2 )
!     cs_now = cos(theta_now); sn_now = sin(theta_now)
	   call sincos( theta_now, sn_now, cs_now )	! sincos.c 利用
      D_inv  =  1 / ( cs_now - cs_0 )
      error_ =  tn_0 * sn_now * D_inv - time_end_ + log( sn_p / sn_m )
      if ( error_ == 0 ) then
         theta_end = theta_now
      else
         dfdth_ =  tn_0 * ( 1 - cs_0 * cs_now ) * D_inv**2 + sn_0 / ( 2 * sn_p * sn_m )
         theta_end = theta_now - error_ / dfdth_
      end if
!     print '(3i3,2es24.16,2es11.3)', i,j,IR, theta_now / pi, theta_end / pi, error_ / dfdth_, error_ * b_0 / g_0
   end subroutine find_theta
end subroutine find_hyperbola
