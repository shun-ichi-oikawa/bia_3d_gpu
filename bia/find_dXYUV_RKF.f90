subroutine find_dXYUV_RKF( eps_th_, xe, YZ, YN, error, IER )
   IMPLICIT REAL*8( A - H, O - Z )  ! YZ, YNはいずれも相対座標・速度 eps_thは運動方程式の二倍
   integer, parameter         :: Neq = 4
   real(8), parameter         :: relerr = 0, abserr = 1d-16    ! 1d-10以上厳しくしても意味が無い
   real(8), dimension(Neq   ) :: Y0, YZ, YN, esterr, WK1, WK2, WK3
   real(8), dimension(Neq,11) ::                               WK4
   external moon_2D
!
   Y0 = YZ; x0 = 0; init = 1
   CALL RKF(Neq,moon_2D,x0,xe,Y0,INIT,RELERR,ABSERR,YN,ESTERR,NDE,IER,WK0,WK1,WK2,WK3,WK4)
   error = maxval(esterr)
!
end subroutine find_dXYUV_RKF
!
subroutine moon_2D( X, Y, F ) ! (xi,eta)での二次元相対運動
   implicit real*8 ( a - z )
   real(8), dimension(4) :: Y, F
   common / const / eps_th
   rsq  = Y(1)**2 + Y(2)**2
   aaa  = eps_th / ( rsq * sqrt(rsq) )
   F(1) = Y(3); F(3) = Y(1) * aaa
   F(2) = Y(4); F(4) = Y(2) * aaa
end subroutine moon_2D
!
subroutine rkf(neq,func,x0,xe,yh0,init,relerr,abserr,yhn,esterr,nde,ier,h,yl0,yln,er,wk)
   implicit real*8(a-h,o-z)
   external func
   dimension yh0(neq),yhn(neq),yl0(neq),yln(neq),esterr(neq),er(neq),wk(neq,11)
   data itemax / 100000 /
   data epsmin / 1d-15 /
!
   if ( init == 1 ) then
!     -------------- initialization (first call)--------------------------
      yl0(:) = yh0(:)
!     -------- set initial step size ---------
      call func(x0,yl0,wk(1,1))
      toldy = huge(toldy)
      do i = 1, neq
         if ( wk(i,1) /= 0 ) then
            tol = relerr * abs( yl0(i) ) + abserr
            toldy = min( toldy, abs( tol / wk(i,1) ) )
         end if
      end do
      hmin = epsmin * ( xe - x0 )
      if ( toldy == huge(toldy) ) then
         h = hmin
      else
         h = min( ( xe - x0 ) / 100, toldy**0.2d0 )
      end if
!     ----------------------------------------
      init = 2
      nde = 1
   else
      nde = 0
   end if
!  ------------------------------------------------------------------
   ier = 0
   x = x0
!
!  begin main iteration
   do iter = 1, itemax
      call flstep( neq, func, x, h, yl0, yln, er, wk ); nde = nde + 6
!     ----------- compute error tolerances ---------------------
      erret = 0
      do i = 1, neq
         et = relerr / 2 * ( abs( yl0(i) ) +  abs( yln(i) ) ) + abserr
         if ( et == 0 ) then
            go to 20000
         else
            erret = max( erret, abs(er(i)) / et )
         end if
      end do
      if ( erret > 1 ) then
!        ----------- unsuccessful step ------------------------------
         if ( erret >= 59049 ) then
            h = 0.1d0 * h
         else
            h = 0.9d0 * h / erret**0.2d+0
         end if
         if ( h <= epsmin ) go to 30000
      else if ( x + h < xe ) then
!        ------------ successful step ( x + h < end point)--
         call fhstep( neq, func, x, h, yh0, yhn, er, wk ); nde = nde + 12
         x = x + h
         x0 = x
         yl0(:) = yln(:)
         yh0(:) = yhn(:)
!
!        write(1,'(5es24.16)') x, yhn(1:4)
!
!        ------- choose next step -------------------------------------
         if ( erret <= 1.889568d-4 ) then; h = 5 * h
                                     else; h = 0.9d0 * h / erret**0.2d0
         end if
      else
!		 ------------ successful step ( x + h > = end point)
         h = xe - x
         call flstep( neq, func, x, h, yl0, yln, er, wk )
         call fhstep( neq, func, x, h, yh0, yhn, er, wk ); nde = nde + 18
!      	 -------- estimate global error -----
         esterr(:) = ( yln(:) - yhn(:) ) / 31
!			 ------------------------------------------
         x0 = xe
         yl0(:) = yln(:)
         yh0(:) = yhn(:)
         return
      end if
   end do
!  end main iteration
10000 continue; ier = 10000; write(*,10001) x; return
20000 continue; ier = 20000; write(*,20001) x; return
30000 continue; ier = 30000; write(*,30001) x; return
!
10001 format('(subr.-rkf) trouble(too many iterations))'         ,' at x = ', es15.7)
20001 format('(subr.-rkf) trouble(inappropriate error tolerance)',' at x = ', es15.7)
30001 format('(subr.-rkf) trouble(too small step size)'          ,' at x = ', es15.7)
end subroutine rkf
!
!-----------------------------------------------------------------------------------------------------------------------------
!
subroutine flstep( neq, func, x, h, y0, yn, er, wk )
   implicit real*8 ( a - h, o - z )
   external func
   dimension y0(neq), yn(neq), er(neq), wk(neq,11)
   lindex = 1
   call fehl( lindex, neq, func, x, h, y0, yn, er, &
              wk(1,1),wk(1,2),wk(1,3),wk(1,4),wk(1,5),wk(1,6),wk(1,7),wk(1,8),wk(1,9),wk(1,10),wk(1,11))
   return
end subroutine flstep
!
!-----------------------------------------------------------------------------------------------------------------------------
!
subroutine fhstep( neq, func, x, h, y0, yn, er, wk )
   implicit real*8(a-h,o-z)
   external func
   real(8) y0(neq), yn(neq), er(neq), wk(neq,11)
   x1 = x
   h1 = h / 2
   lindex = 0
   call fehl( lindex, neq, func, x1, h1, y0, yn, er, &
              wk(1,1),wk(1,2),wk(1,3),wk(1,4),wk(1,5),wk(1,6),wk(1,7),wk(1,8),wk(1,9),wk(1,10),wk(1,11))
   x1 = x1 + h1
   y0 = yn
   call fehl( lindex, neq, func, x1, h1, y0, yn, er, &
              wk(1,1),wk(1,2),wk(1,3),wk(1,4),wk(1,5),wk(1,6),wk(1,7),wk(1,8),wk(1,9),wk(1,10),wk(1,11))
   return
end subroutine fhstep
!
!-----------------------------------------------------------------------------------------------------------------------------
!
subroutine fehl( lindex, neq, func, x ,h, y0, yn, er, ak1, ak2, ak3, ak4, ak5, ak6, w2, w3, w4, w5, w6 )
   implicit real*8 ( a - h, o - z )
   real(8), dimension (neq) :: y0, yn, ak1, ak2, ak3, ak4, ak5, ak6, er, w2, w3, w4, w5, w6
   parameter( one = 1, two = 2, thr = 3, twl = 12 )
   parameter( al2 =    one /    4, al3 = thr /  8, al4 = twl / 13, al5 = one, al6 = one / 2 )
   parameter( b21 =    one /    4, &
              b31 =    thr /   32, b32 =     9d0 /   32, &
              b41 = 1932d0 / 2197, b42 = -7200d0 / 2197, b43 =  7296d0 /  2197, &
              b51 =  439d0 /  216, b52 =      -8       , b53 =  3680d0 /   513, b54 = -845d0 / 4104, &
              b61 =   -8d0 /   27, b62 =       2       , b63 = -3544d0 /  2565, b64 = 1859d0 / 4104, b65 = -11d0 / 40 )
   parameter( ga1 =   16d0 /  135, ga3 = 6656d0 / 12825, ga4 = 28561d0 / 56430, ga5 =   -9d0 /   50, ga6 =   two / 55 )
   parameter( da1 =    one /  360, da3 = -128d0 /  4275, da4 = -2197d0 / 75240, da5 =    one /   50, da6 =   two / 55 )
   call func( x          , y0, ak1 ); w2 = y0 + h *   b21 * ak1
   call func( x + al2 * h, w2, ak2 ); w3 = y0 + h * ( b31 * ak1 + b32 * ak2 )
   call func( x + al3 * h, w3, ak3 ); w4 = y0 + h * ( b41 * ak1 + b42 * ak2 + b43 * ak3 )
   call func( x + al4 * h, w4, ak4 ); w5 = y0 + h * ( b51 * ak1 + b52 * ak2 + b53 * ak3 + b54 * ak4 )
   call func( x + al5 * h, w5, ak5 ); w6 = y0 + h * ( b61 * ak1 + b62 * ak2 + b63 * ak3 + b64 * ak4 + b65 * ak5 )
   call func( x + al6 * h, w6, ak6 ); yn = y0 + h * ( ga1 * ak1             + ga3 * ak3 + ga4 * ak4 + ga5 * ak5 + ga6 * ak6 )
   if ( lindex == 1 )                 er =      h * ( da1 * ak1             + da3 * ak3 + da4 * ak4 + da5 * ak5 + da6 * ak6 )
end subroutine fehl
