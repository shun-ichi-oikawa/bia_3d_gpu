subroutine find_dXYUV4      ( eps_th, dt, YZ, YN, theta_end, error_end, ERR, ITR, IR, IER )
   IMPLICIT REAL*8( A-H, K, L, O-Z )! YZ, YNはいずれも相対座標・速度 eps_thは運動方程式の二倍
   real(8), dimension(4) :: YN, YZ
   common / com_time_theta / b_0, g_0, theta_0, tn_0, cs_0, sn_0, com_dt
   data pi / 3.14159265358979323846264338327950288419716939937510582097494d0 /
   external func_time_theta
   X_ini = YZ(1); Y_ini = YZ(2); U_ini = YZ(3); V_ini = YZ(4)	! (X,Y,U,V)系での相対座標・速度
   r_g = X_ini * U_ini + Y_ini * V_ini			! r_ cdot  g_
   Lz  = X_ini * V_ini - Y_ini * U_ini			! r_ times g_
   PHI = atan2( Lz, r_g )						! 位置ベクトルに対して速度ベクトルのなす角度
!  print *, 'r_g      = ', r_g
!  print *, 'Lz       = ', Lz
!  print *, 'PHI / pi = ', PHI / pi
   r   = sqrt( X_ini**2 + Y_ini**2 )
   g2  =       U_ini**2 + V_ini**2	! g^2
   W   = g2 + 2 * eps_th / r		! energy
   g_0 = sqrt(W)
   b_0 = eps_th / W
   b   = abs(Lz) / g_0
   theta_0 = atan2 ( abs(Lz), b_0 * g_0 )	! [ 0, pi/2]
   chi = pi - 2 * theta_0
   cs_0 = cos(theta_0)
   sn_0 = sin(theta_0); tn_0 = sn_0 / cs_0
   if ( abs( cs_0 + b * sn_0 / r ) > 1 ) then
      theta = 1d-30
   else
      theta = acos( cs_0 + b * sn_0 / r )	! θをとりあえず正で定義する
   end if
   if ( Lz * r_g < 0 ) then					! θは、これで正しい符号となる
      theta = - theta
   end if
!
   cs = cos(theta); sn = sin(theta)
   x = r * cs;	!  x = r * cs_0 + b * sn_0	! x > 0
   y = r * sn								!
   u = g_0 / sn_0 * (     cs_0 * sn )
   v = g_0 / sn_0 * ( 1 - cs_0 * cs )
   alpha = atan2( x * Y_ini - y * X_ini, x * X_ini + y * Y_ini )
   beta  = atan2( u * V_ini - v * U_ini, u * U_ini + v * V_ini )
!
   time = ( y / sn_0 + b_0 * log (   sin( ( theta_0 + theta ) / 2 )         &
	                               / sin( ( theta_0 - theta ) / 2 ) ) ) / g_0
!  print *, 'time_ini     =', time
!  -------------------------------------------------------------------------------------------------
!goto 9
   com_dt = time                ; error_NOW1 = func_time_theta(theta)
   com_dt = time + error_NOW1   ; error_NOW2 = func_time_theta(theta)
   time   = time + error_NOW1
9 continue
!  -------------------------------------------------------------------------------------------------
   if ( Lz < 0 ) time = - time	! 時間を正方向に進める（時間軸に対して角度θは負方向に進む）
!  -------------------------------------------------------------------------------------------------
   com_dt = time + dt
   theta_L = theta; if ( Lz < 0 ) theta_L = - theta_L
   if ( time * com_dt < 0 .and. time + com_dt > 0 ) theta_L = - theta_L
   if ( com_dt < 0 ) then
      theta_R = -1d-16
   else
      theta_R =   theta_0 * ( 1 - 100 * epsilon(1d0) )
   end if
!  print *, 'LZ =', Lz, theta_0 / pi
!  print *, 'end of error_NOW'
   FL = func_time_theta(theta_L)
   FR = func_time_theta(theta_R)
   if ( FL * FR > 0 ) then
!     theta_R = - theta_R
      theta_x =   theta_R
      theta_R =   theta_L
      theta_L = - theta_x

!     theta_x =   theta_R
!     theta_R =   theta_L
!     theta_L =   theta_x
!     pause
   end if
   FL = func_time_theta(theta_L)
   FR = func_time_theta(theta_R)
!  print *, theta_L / pi, FL
!  print *, theta_R / pi, FR
!  print *, ''
!
    call NLBIT_( theta_end, func_time_theta, theta_L, theta_R, ERR, ITR,IR,IER)
!  print *, 'theta_end =', theta_end 
    error_end = func_time_theta(theta_end)
!
!  ------------------------------------------
!
23   if ( Lz < 0 ) theta_end = - theta_end	! 時間を正方向に進める（時間軸に対して角度θは負方向に進む）
!
   cs_end = cos(theta_end)
   sn_end = sin(theta_end)

    g_end = g_0 / sn_0
    u_end = g_end *       cs_0 * sn_end
    v_end = g_end * ( 1 - cs_0 * cs_end )

   du_end = u_end - u
   dv_end = v_end - v

    r_end = b * sn_0 / ( cs_end - cs_0 )
    x_end = r_end * cs_end
    y_end = r_end * sn_end
!
!  -------------------------------------------------------------------------------------------------------------------------
!
    r1end =  x_end * cos(alpha) -  y_end * sin(alpha)
    r2end =  x_end * sin(alpha) +  y_end * cos(alpha)
   dg1end = du_end * cos(alpha) - dv_end * sin(alpha)
   dg2end = du_end * sin(alpha) + dv_end * cos(alpha)
   if ( Lz < 0 ) then
      dg1end = - dg1end; dg2end = - dg2end
   end if
!
   g1end = YZ(3) + dg1end
   g2end = YZ(4) + dg2end
   YN(1) = r1end; YN(2) = r2end
   YN(3) = g1end; YN(4) = g2end

  return

!  -------------------------------------------------------------------------------------------------------------------------
!  print *, 'dg_ell  =', 2 * g_0 / tan( theta_0 ) * sin( theta_ell )
!  print *, 'dg_∞   =', 2 * g_0                  * sin( chi / 2   )
   print '(3x,4a24)',  'x', 'y', 'u', 'v'
   print '(a,es23.16,3es24.16)', 'ini ', x,         y,         u,         v
!  print '(a,es23.16,3es24.16)', 'est ', x + u * dt, y + v * dt
   print '(a,es23.16,3es24.16)', 'end ', x_end    , y_end    , u_end    , v_end
   print '(a,es23.16,3es24.16)', 'Δ  ', x_end - x, y_end - y, u_end - u, v_end - v
   print *, ''
   print *, 'impact param b =', b
   print *, 'alpha     / pi =', alpha     / pi
   print *, 'beta      / pi =', beta      / pi
   print *, 'chi       / pi =', chi       / pi
   print *, 'theta_0   / pi =', theta_0   / pi
   print *, 'theta_ini / pi =', theta     / pi
   print *, 'theta_end / pi =', theta_end / pi
   print *, 'theta_Δ  / pi =',(theta_end-theta) / pi
!  print *, 'Lzini = ',                                   x * v - y * u
!  print *, 'Lzend = ', x_end * v_end - y_end * u_end                   , &
   print *, 'r dot g ', r_g
   print *, 'Lzini = ', YZ(1) * YZ(4) - YZ(2) * YZ(3)
   print *, 'Lzend = ', YN(1) * YN(4) - YN(2) * YN(3)                   , &
            'error = ', YN(1) * YN(4) - YN(2) * YN(3) - ( YZ(1) * YZ(4) - YZ(2) * YZ(3) )
   print *, 'W_ini = ', u    **2 + v    **2               + 2 * eps_th / sqrt(x    **2 + y    **2)
   print *, 'W_end = ', u_end**2 + v_end**2               + 2 * eps_th / sqrt(x_end**2 + y_end**2), &
            'error = ', u_end**2 - u**2 + v_end**2 - v**2 + 2 * eps_th / sqrt(x_end**2 + y_end**2)  &
                                                          - 2 * eps_th / sqrt(x    **2 + y    **2)
   print *, '---------------------------------------------------------------------------'
   print *, ''
!  stop
!
end subroutine find_dXYUV4
!
real(8) function func_time_theta ( theta )
   IMPLICIT REAL*8( A-H, O-Z )
   common / com_time_theta / b_0, g_0, theta_0, tan_0, cos_0, sin_0, time_to_find
   data pi / 3.14159265358979323846264338327950288419716939937510582097494d0 /
!
!  print *, 'theta =', theta / pi
   if ( theta == 0 ) then
      F = 1
      func_time_theta = - time_to_find
      return
   end if
   aaa = tan_0 * sin(theta) / ( cos(theta) - cos_0 )
   theta_p = ( theta_0 + theta ) / 2
   theta_m = ( theta_0 - theta ) / 2
   F = sin ( theta_p ) / sin( theta_m )
   if ( F < 0 ) then
      print *, 'theta / pi = ', theta / pi; print *, 'F =', F
!     pause
   end if
   bbb = log(F)
   func_time_theta = b_0 * ( aaa + bbb ) / g_0 - time_to_find
!  print *, aaa, bbb, aaa+bbb
!  print '(es35.26,2es12.3)', theta / pi, func_time_theta
end function func_time_theta
!
!--------------------------------------------------------------------------------
!
real(8) function func_time_theta2 ( d_theta )
   IMPLICIT REAL*8( A-H, O-Z )
   common / com_time_theta / b_0, g_0, theta_0, tan_0, cos_0, sin_0, time_to_find
   data pi / 3.14159265358979323846264338327950288419716939937510582097494d0 /
!
   theta = theta_0 + d_theta
   if ( d_theta < 0 ) then
!     F = 1
!     func_time_theta = - time_to_find
!     pause 'STOP theta = theta_0'
 	  STOP
   else
      aaa = tan_0 * ( tan_0 * cos(d_theta) - sin(d_theta) ) / ( tan_0 * sin(d_theta) - 1 + cos_0 )
      theta_p = theta_0 - d_theta / 2
      theta_m =           d_theta / 2
      F = sin ( theta_p ) / sin( theta_m )
      if ( F < 0 ) then
         print *, 'theta / pi = ', theta / pi; print *, 'F =', F
!        pause
      end if
   end if
   bbb = log(F)
   func_time_theta2 = b_0 * ( aaa + bbb ) / g_0 - time_to_find
!  print '(es25.16,2es12.3)', theta / pi, func_time_theta
end function func_time_theta2
!
SUBROUTINE NLBIT_(X, FUNC, XL, XR, EPS, ITR, IR, IER)
!***********************************************************************
!  SOLUTION OF A NON-LINEAR EQUATION BY BISECTION METHOD AND REGULA    *
!     FALSI METHOD.                                                    *
!  PARAMETERS                                                          *
!    (1) X: VALUE OF THE SOLUTION                                      *
!    (2) FUNCT: THE FUNCTION OF THE NON-LINEAR EQUATION                *
!    (3) XL: LEFT HAND SIDE VALUE OF THE SOLUTION                      *
!    (4) XR: RIGHT HAND SIDE VALUE OF THE SOLUTION                     *
!    (5) EPS: TOLERANCE FOR CONVERGENCE                                *
!    (6) ITEND: MAXIMUM NUMBER OF REPITITION                           *
!    (7) IER: ERROR CODE                                               *
!  COPYRIGHT  T. OGUNI   JUNE 30 1989   VERSION 1.0                    *
!***********************************************************************
       IMPLICIT REAL*8(A-H,O-Z)
!     data pi/3.14159265358979323846264338327950288419716939937510582097494d0 /
      data pi/3.14159265358979323846264338327950288419716939937510582d0/
!
      IER = 0
!     EPS1 = 1.0D-10
      EPS1 = 1.0d-12
      X = XL; FL = FUNC(X); IF ( FL == 0 ) RETURN
      X = XR; FR = FUNC(X); IF ( FR == 0 ) RETURN
      IF ( FR * FL > 0 ) THEN; IER = 2
         X = ( XL * XR ) / 2
         WRITE(*,*) '(SUBR. NLBIT) NO CALUCULATION FOR F(XL)･F(XR) > 0.'
         print *, FL, ' at theta_L / pi = ', XL / pi
         print *, FR, ' at theta_R / pi = ', XR / pi
!         RETURN
!        XL = - XL
         XR = - XL
!     pause
      ENDIF
!
      dx = 100
!
      DO 250 N = 1, ITR; IR = N; dx0 = dx
         IER = 0
         X  = ( XL + XR ) / 2 ; F  = FUNC(X)      ; IF ( F == 0 ) RETURN
         IF ( F * FL .LT. 0 ) THEN
            XR = X
            FR = F
         ELSE
            XL = X
            FL = F
         END IF
!     print *, XL/pi,XR/pi
!	print *, FL,FR
!     if (abs(FL) < abs(FR)) then
         DF = FL - FR; DX = ( XR - XL ) * FL / DF
         X  = XL + DX; F  = FUNC(X); IF ( F == 0 ) RETURN
!     else
!        DF = FL - FR; DX = ( XR - XL ) * FR / DF
!        X  = XR + DX; F  = FUNC(X); IF ( F == 0 ) RETURN
!     end if
         IF ( FL * F < 0 ) THEN
            XR = X; IER = 88; if ( XR == xrr ) return	! 収束
            FR = F
         ELSE
            XL = X; IER = 99; if ( XL == xll ) return	! 収束
            FL = F
         END IF
!
         SX = XR - XL
         IF ( ABS( X) .GT. 1.0 ) SX = SX / X
         IF ( ABS(SX) .LE. EPS ) THEN
            X = XR; IF ( ABS(FR) .LE. EPS1) RETURN
            X = XL; IF ( ABS(FL) .LE. EPS1) RETURN
         ENDIF
	   xrr = XR; xll = XL
  250 CONTINUE
!
      IER = 1
!     WRITE(*,*) '(SUBR. NLBIT) REQUIRED ACCRACY IS NOT OBTAINED.'
	if ( abs(FL) < abs(FR) ) then
          x = XL; F = FL
      else
          x = XR; F = FR
	end if
      if( abs(F) < EPS ) then; IER = 3
      else                   ; IER = 1
!      	print *, IR, X/pi, F
	end if
!	print *, dx, dx0, dx+dx0
!     print *, XL- (XL+DX); pause
      RETURN
END SUBROUTINE NLBIT_