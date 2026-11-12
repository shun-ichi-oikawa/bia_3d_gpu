subroutine find_parabola( eps_th, DT, YZ, YN, error_time )
   implicit real(8) ( a - h, o - z )      !放物線軌道の初期値を(ξ,η)系から( x, y)系へ変換する。
   real(8), dimension(2) ::  r,  g        ! (ξ,η)系、角運動量Lは+ζ軸方向 L = r_min * g_max
   real(8), dimension(4) :: YZ, YN        ! (ξ,η)系の位置＆速度、角運動量Lは+ z軸方向
   integer, parameter :: NI = 99
   real(8), dimension(NI) :: theta_Newton, error_Newton
!
   r_ini = sqrt( YZ(1)**2 + YZ(2)**2 )    ! |r|
   g_ini =                  YZ(4)         ! |g| = g(2) のみ、g(1) = g_xi = 0
   call VECTOR( YZ, r_dot_g, r_crs_g )    ! r と g の内積・外積の計算
!
!  alpha = atan2( r_crs_g, r_dot_g )      ! ベクトルrからgへの角度
   tmp = 1 / sqrt( r_crs_g**2 + r_dot_g**2 )
   sn_alpha = r_crs_g * tmp      ! sin_alpha
   cs_alpha = r_dot_g * tmp      ! cos_alpha
!
!  最短距離 r_min at (x,y) = ( r_min, 0) で最速 g_max: ( g_x, g_y ) = ( 0, g_max )
   g_max = - 2 * eps_th / r_crs_g   ! g_max * L = - eps > 0
   r_min =   r_crs_g    / g_max     ! L = r_min * g_max > 0
!  r_min = - r_crs_g**2 / (2 * eps_th)

   x_ini =       2 * r_min - r_ini                                !
   y_ini = sign( 2 * r_min * sqrt( r_ini / r_min - 1 ), r_dot_g ) ! y_ini と r_dot_g は同符号

   tmp   =   g_ini / r_ini  ! ベクトルr=(x,y) を角αだけ回転し、長さを g_ini / r_ini 倍する。
   u_ini = ( x_ini * cs_alpha - y_ini * sn_alpha ) * tmp ! u_ini = u_ini
   v_ini = ( x_ini * sn_alpha + y_ini * cs_alpha ) * tmp ! v_ini = v_ini

   theta_ini = atan2( y_ini, x_ini )
   time_ini  = time_parabola(theta_ini)
!
!
!  Cardanoによる三次方程式+Newton法の初期値
   time_end  = time_ini + DT
   time_end_ = 3 * g_max / ( 4 * r_min ) * time_end
!  theta_gss = 2 * atan( CBRT ( time_end_ ) )
   theta_gss = 2 * atan(    F ( time_end_ ) )
   theta     = theta_gss
   err_tm_min = 1d99; i_min = 1
   print '(a2,3a24,a4)', 'i', 'θ', 'Δθ', 'error in time', 'min'
   NEWTON: do i = 1, NI
      err_tm = time_parabola(theta) - time_end
      error_Newton(i) = err_tm
      dtheta = - err_tm / dt_dtheta(theta)
      if ( abs(err_tm) < err_tm_min ) then
         err_tm_min = abs(err_tm)
         thetamin = theta
         i_min = i
      end if

      theta = theta + dtheta
      theta_Newton(i) = theta
      print '(i3,es23.16,2es24.16,i4)', i, theta, dtheta, - err_tm, i_min
      if (         dtheta == 0     ) exit
      if ( theta - dtheta == theta ) exit

      do j = 1, i - 1
         if( error_Newton(i) == error_Newton(j) ) then
               err_tm = error_Newton(i_min)
               theta  = theta_Newton(i_min)
               exit NEWTON
         end if
      end do
   end do NEWTON

   print *, i_min

   if( i== NI + 1) theta = thetamin	! 収束せず
   theta_end = theta
   beta_end  = theta_end / 2
   print '(a,2 a 25  )', '    ','time'   ,'theta'
   print '(a,2es25.16)', 'ini:', time_ini, theta_ini
   print '(a,2es25.16)', 'gss:', time_parabola(theta_gss), theta_gss
   print '(a,2es25.16)', 'end:', time_end, theta_end
   print '(a, es25.16)', 'end:', time_parabola(theta_end)

!	cs_beta_end = cos( beta_end ); sn_beta_end = sin( beta_end )
	call sincos( beta_end, sn_beta_end, cs_beta_end)	! sincos.c 利用
   r_end =   r_min /   cs_beta_end**2
   g_end =   g_max *   cs_beta_end
   x_end =   r_end * cos(theta_end)
   y_end =   r_end * sin(theta_end)
   u_end = - g_end *   sn_beta_end
   v_end =   g_end *   cs_beta_end
!
!  gamma = atan2( x_ini * YZ(2) - y_ini * YZ(1), x_ini * YZ(1) + y_ini * YZ(2) )
!  cs_gamma = cos(gamma); sn_gamma = sin(gamma)
!
	tmp = 1 / sqrt( ( x_ini * YZ(2) - y_ini * YZ(1) )**2 + ( x_ini * YZ(1) + y_ini * YZ(2) )**2 )
   sn_gamma = ( x_ini * YZ(2) - y_ini * YZ(1) ) * tmp
   cs_gamma = ( x_ini * YZ(1) + y_ini * YZ(2) ) * tmp

   YN(1) =   x_end * cs_gamma - y_end * sn_gamma   ! r_xi_end
   YN(2) =   x_end * sn_gamma + y_end * cs_gamma   ! r_eta_end
   YN(3) =   u_end * cs_gamma - v_end * sn_gamma   ! g_xi_end
   YN(4) =   u_end * sn_gamma + v_end * cs_gamma   ! g_eta_end
   error_time = time_end - time_parabola(theta_end)

   print '(4x,4a18   )',        'x'   ,'y'   ,'u'   ,'v'
   print '(a,4es18.10)', 'ini:', x_ini, y_ini, u_ini, v_ini
   print '(a,4es18.10)', 'end:', x_end, y_end, u_end, v_end
   print *, 'L_ini =', x_ini * v_ini - y_ini * u_ini, u_ini**2 + v_ini**2 + 2 * eps_th / sqrt( x_ini**2 + y_ini**2 )
   print *, 'L_end =', x_end * v_end - y_end * u_end, u_end**2 + v_end**2 + 2 * eps_th / sqrt( x_end**2 + y_end**2 )
   return

   !----------------------------------------------------------
   contains
   real(8) function F(time_)  ! t = t(θ)
      implicit real(8) (a-h,o-z)
      FP = CBRT( sqrt( 1 + time_**2 ) + time_ )
      FM = CBRT( sqrt( 1 + time_**2 ) - time_ )
      F  = FP - FM
   end function F
   real(8) function CBRT(a)
      implicit real(8) (a-z)
      real( 8), parameter :: one = 1, pow = one / 3
      CBRT = a**pow
   end function CBRT
   real(8) function time_parabola(theta)  ! t = t(θ)
      implicit real(8) (a-h,o-z)
		call sincos( theta, sn, cs )	! sincos.c 利用
      time_parabola = ( 2 + cs ) * sn / ( 1 + cs )**2 * 4 * r_min / ( 3 * g_max )
   end function time_parabola
!
   real(8) function dt_dtheta(theta)
      implicit real(8) (a-h,o-z)
      dt_dtheta = r_min / ( g_max * cos(theta/2)**4 )
   end function dt_dtheta
!
   real(8) function energy( r, g )
      real(8) r(2), g(2)
      energy = sum( g**2 ) + 2 * eps_th / sqrt( sum( r**2 ) )
   end  function energy
!
end subroutine find_parabola
!
subroutine VECTOR( Y, DOT, CRS )
   real(8) DOT, CRS
   real(8), dimension(4) ::   Y(4)  ! 二次元ベクトル in (ξ,η)
   DOT = Y(1) * Y(3) + Y(2) * Y(4)  ! r dot g
   CRS = Y(1) * Y(4) - Y(2) * Y(3)  ! r cross g の ζ成分：
   return
end subroutine VECTOR
