subroutine Invariance_CPU( N_p, eps_th, X, Y, Z, U, V, W, KE, PE, Lx, Ly, Lz, Ua, Va, Wa, Xa, Ya, Za )
   IMPLICIT REAL*8( A - H, K - L, O - Z )
   real(08), dimension(N_p) ::  X, Y, Z, U, V, W


!

   Xa = sum(  X ); Ua = sum(  U ); Lx = sum(  Y *  W -  Z *  V )
   Ya = sum(  Y ); Va = sum(  V ); Ly = sum(  Z *  U -  X *  W )
   Za = sum(  Z ); Wa = sum(  W ); Lz = sum(  X *  V -  Y *  U )
   KE = sum(  U *  U +  V *  V +  W *  W )
!
   PE = 0
   do i = 1, N_p - 1   ; Xi = X(i); Yi = Y(i); Zi = Z(i)
      do j = i + 1, N_p; Xj = X(j); Yj = Y(j); Zj = Z(j)
         PE = PE + 1 / sqrt( ( Xi - Xj )**2 + ( Yi - Yj )**2 + ( Zi - Zj )**2 )
      end do
   end do
   PE = eps_th * PE  ! 同種コードにおいてeps_th は、運動方程式の二倍（ 2 b_0_th ）なので、ここでは係数2が不要
!
!	PE = 0
!	do i = 1, N_p - 1; tmp = 0	;	Xi = X(i); Yi = Y(i); Zi = Z(i)
!		do j = N_p, i + 1, -1	;	Xj = X(j); Yj = Y(j); Zj = Z(j)
!			tmp = tmp + 1 / sqrt( ( Xi - Xj )**2 + ( Yi - Yj )**2 + ( Zi - Zj )**2 )
!		end do
!		PE = PE + tmp
!	end do
!	PE = eps_th * PE  ! eps_th は、運動方程式の二倍（ 2 b_0_th ）なので、ここでは係数2が不要
!
end subroutine Invariance_CPU
