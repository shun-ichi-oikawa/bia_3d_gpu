subroutine FIND_ELLIPSE( eps_th, dt, YZ, YN, error, ERR, ITR, IR, IER )
   implicit real(8) ( a - h, K - L, o - z )
   real(8), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494d0
   real(8), dimension(4) ::  YZ, YN
!
   X_ini = YZ(1); Y_ini = YZ(2); U_ini = YZ(3); V_ini = YZ(4)
!
   g2_ini =       U_ini**2 + V_ini**2        ! 相対速度の二乗
   r_ini  = sqrt( X_ini**2 + Y_ini**2 )      ! 相対距離
   W      = g2_ini + 2 * eps_th / r_ini      ! W = eps_th / a = - g_a**2
   if ( W >= 0 ) then                        ! 楕円ではない：双曲線または放物線運動
!     print '(a,3f11.6)', 'W, K, U = ', W, g2_ini, 2 * eps_th / r_ini
      stop 'Not an ellipse!'
   end if
!
   L       = X_ini * V_ini - Y_ini * U_ini   ! 角運動量
   r_g_ini = X_ini * U_ini + Y_ini * V_ini   ! r cdot g
!
   g_a = sign( sqrt( - W ), L )              ! 反時計回りが正：角運動量と同符号
   a  =  eps_th / W                          ! 長半径：正で定義
   b  =  L / g_a                             ! 短半径：正で定義
   omega_a = g_a / a                         ! g_a / a
!  xi_ecc = asin( L / ( a * g_a ) ); ecc = cos(xi_ecc)
   xi_ecc = asin( L / ( a * g_a ) ); ecc = sqrt( (a * g_a)**2 - L**2 ) / ( a * g_a )
!  -----------------------------------------------------------------------------------
   if ( r_g_ini == 0 )   then          ! 初期条件が近日点または遠日点でない場合
   	if ( r_ini < a ) then          ! 近日点または遠日点
      	xi_ini = 0                       ! 近日点
		else
     		if ( r_ini == a ) then
     			xi_ini = atan2( Y_ini, X_ini ); ecc = 0; xi_ecc = pi / 2
	   	else
				xi_ini = sign ( pi, L )          ! 遠日点
	   	end if   
      end if
  	else                                ! 円運動
		xi_ini = atan2(r_g_ini / g_a, a - r_ini)
!     print *, 'xi_ini       =',xi_ini
!     print '(2i6)',i, j
	end if
!  -----------------------------------------------------------------------------------
!  theta = atan2( b * sin(xi_ini), a * ( cos(xi_ini) - ecc ) )
   alpha = atan2( Y_ini, X_ini ) - atan2( b * sin(xi_ini), a * ( cos(xi_ini) - ecc ) )
!  -----------------------------------------------------------------------------------
   tau_ini  = xi_ini - ecc * sin(xi_ini)
   time_ini = tau_ini / omega_a 
   time_end = time_ini + dt
   tau_end  = omega_a * time_end
!  print *, 'theta_lab =', atan2(          Y_ini ,            X_ini          ) / pi
!  print *, 'theta_x-y =', atan2( b * sin(xi_ini), a * ( cos(xi_ini) - ecc ) ) / pi
! 	print *, 'g_a       =', g_a
! 	print *, 'a         =', a
! 	print *, 'b         =', b
! 	print *, 'xi_ecc    =', xi_ecc
! 	print *, 'alpha /pi =', alpha / pi
! 	print *, 'ξ_ini/pi =', xi_ini/ pi
! 	print *, 'W         =', W
! 	print '(2i6)',i, j
! 	print *, 'omega_a   =', omega_a
! 	print *, 'tau_ini/pi=', tau_ini/ pi
! 	print *, 'tau_end/pi=', tau_end/ pi
!  print *, 'ecc       =', ecc
!  print *, 'x_ini     =', a * ( cos(xi_ini) - ecc )
!  print *, 'y_ini     =', b *   sin(xi_ini)
!  print *, 'r_ini     =', sqrt( a**2 * ( cos(xi_ini) - ecc )**2 + b**2 * sin(xi_ini)**2 )
!  print *, 'r_ini     =', a * ( 1 - ecc * cos(xi_ini) )
!  print *, 'R_ini     =', r_ini
!  print *, 'L         =', L
!  print *, 'omega_a dt=', omega_a * dt
!  print *, 'time_ini  =', time_ini
!  print *, 'time_end  =', time_end
!  print *, 'period    =', 2 * pi / omega_a
   if ( dt == 0 ) then
      error = 0
      YN = YZ
      return
   end if
!
   call Newton( xi_ini, xi_end, error )   ! ここでニュートン法を使ってxi_endをもとめる
!
! 	print *, 'xi_end    =', xi_end
   cs_end = cos( xi_end )
   sn_end = sin( xi_end )
   xx_end =   a * ( cs_end - ecc )
   yy_end =   b *   sn_end
   rr_end =   a * ( 1 - ecc * cs_end )
   uu_end = - a * g_a * sn_end / rr_end
   vv_end =   b * g_a * cs_end / rr_end
!
   cs_alp = cos(alpha)
   sn_alp = sin(alpha)
   YN(1) =   xx_end * cs_alp - yy_end * sn_alp
   YN(2) =   xx_end * sn_alp + yy_end * cs_alp
   YN(3) =   uu_end * cs_alp - vv_end * sn_alp
   YN(4) =   uu_end * sn_alp + vv_end * cs_alp
! 	print '(          3x, 4a24)',                'X'   ,        'Y'   ,        'U   ',        'V'
! 	print '(a,es24.16,3es25.16)', 'ini:',         X_ini,         Y_ini,         U_ini,         V_ini
! 	print '(a,es24.16,3es25.16)', 'end:', YN(1)        , YN(2)        , YN(3)        , YN(4)
!  print '(a,es23.16,3es24.16)', 'Δ :', YN(1) - X_ini, YN(2) - Y_ini, YN(3) - U_ini, YN(3) - V_ini
   g2_end =       uu_end**2 + vv_end**2         ! 相対速度の自乗
    r_end = sqrt( xx_end**2 + yy_end**2 )    ! 相対距離
    W_end = g2_end + 2 * eps_th /  r_end     ! W = eps_th / a = - g_a**2
    L_end = xx_end *vv_end - yy_end * uu_end ! 角運動量
   return
   contains
   subroutine Newton( xi_ini, xi_end, error )
      implicit real*8 ( a-h, j-z )
      ir = 0
      error_ini = tau_end - xi_ini + ecc * sin(xi_ini)
      delta_ini = error_ini * (1 - ecc * cos(xi_ini))
      xi_ = xi_ini + delta_ini
         if ( dt == 0 ) then
         xi_end = xi_ini
         return
      end if
      xi_2 = xi_ini + omega_a * dt  ! initial guess, exact when ecc = 0
      if(cos(tau_end) < 0.9) then
      	 xi_2 = tau_end
      else
      	 xi_2 = sign(pi/2, sin(tau_end))
      end if
      xi_gss = xi_2; tau_gss = xi_gss
      error_tm2 = tau_end - xi_2 + ecc * sin(xi_2)
      delta_xi2 = error_tm2 / ( 1 - ecc * cos(xi_2) )
     	if ( abs(error_ini) < abs(error_tm2) ) then
         xi = xi_2 + delta_ini
     	else
        xi = xi_2 + delta_xi2
     	end if
! 		print *, 'xi    =', xi
!     print '(a3,4 a24   )', 'ir', 'xi_  / pi', 'xi / pi', 'delta_xi / pi', 'error in time / pi'
!     print '(i3,4es24.16)',  ir ,  xi_2 / pi ,  xi / pi ,  delta_xi2/ pi ,  error_tm2     / pi
      do while ( xi_1 /= xi_2 .and. ir <= ITR )
         ir = ir + 1
         xi_1 = xi
         error_tm1 = tau_end - xi_1 + ecc * sin(xi_1)
         delta_xi1 = error_tm1 / ( 1 - ecc * cos(xi_1) )
         xi = xi_1 + delta_xi1
!        print '(i3,4es24.16)', ir, xi_1/ pi, xi / pi, delta_xi1, error_tm1
         if ( error_tm1 == 0 .or. error_tm1 == error_tm2 .or. error_tm1 + error_tm2 == 0 ) then
            error = error_tm1;  exit
         end if
         ir = ir + 1
         xi_2 = xi
         error_tm2 = tau_end - xi_2 + ecc * sin(xi_2)
         delta_xi2 = error_tm2 / ( 1 - ecc * cos(xi_2) )
         xi = xi_2 + delta_xi2
!        print '(i3,4es24.16)', ir, xi_2/ pi, xi / pi, delta_xi2, error_tm2
         if ( error_tm2 == 0 .or. error_tm1 == error_tm2 .or. error_tm1 + error_tm2 == 0 ) then
            error = error_tm2 ; exit
         end if
      end do
      xi_end = xi
!     time_end = ( xi_end - ecc * sin(xi_end) ) / omega_a - time_ini
!     print *, 'time_ini / pi =',   time_ini / pi
!     print *, 'time_end / pi =',   time_end / pi 
   end subroutine Newton
end subroutine FIND_ELLIPSE