subroutine Invariance( NP, eps_th, Y, KE, PE, Lx, Ly, Lz, Ua, Va, Wa, Xa, Ya, Za )
	IMPLICIT REAL*8( A - H, K - M, O - Z )
	real(8), dimension(6*NP) :: Y
	real(8), dimension(0:NP-1) :: DUM, ONE
	real(8) aaa(1) / 1 /
	external cublas_ddot ! GPU上で実行されるblas関数
	NEQ = 6 * NP
!
	Xa = sum( Y(1:NEQ:6) ); Ua = sum( Y(4:NEQ:6) )
	Ya = sum( Y(2:NEQ:6) ); Va = sum( Y(5:NEQ:6) )
	Za = sum( Y(3:NEQ:6) ); Wa = sum( Y(6:NEQ:6) )
	Lx = sum( Y(2:NEQ:6) * Y(6:NEQ:6) - Y(3:NEQ:6) * Y(5:NEQ:6) )
	Ly = sum( Y(3:NEQ:6) * Y(4:NEQ:6) - Y(1:NEQ:6) * Y(6:NEQ:6) )
	Lz = sum( Y(1:NEQ:6) * Y(5:NEQ:6) - Y(2:NEQ:6) * Y(4:NEQ:6) )
	KE = sum( Y(4:NEQ:6)**2 + Y(5:NEQ:6)**2 + Y(6:NEQ:6)**2 )
!
	PE = 0; DUM = 0
	do i = 0, NP - 2     ; Xi = Y(6*i+1); Yi = Y(6*i+2); Zi = Y(6*i+3); tmp = 0
		do j = i + 1, NP-1	; Xj = Y(6*j+1); Yj = Y(6*j+2); Zj = Y(6*j+3)
!		do j = NP-1,i+1,-1	; Xj = Y(6*j+1); Yj = Y(6*j+2); Zj = Y(6*j+3)
			tmp = tmp + 1 / sqrt( ( Xi - Xj )**2 + ( Yi - Yj )**2 + ( Zi - Zj )**2 )
		end do
		PE = PE + tmp
!		exit
!		DUM(i) = sum( 1 / sqrt( ( Xi - Y(6*i+6:6*NP:6) )**2 + ( Yi - Y(6*i+7:6*NP:6) )**2 + ( Zi - Y(6*i+8:6*NP:6) )**2 ) )
		DUM(i) = tmp
	end do
	PE = eps_th * PE  ! eps_th は、運動方程式の二倍（ 2 b_0_th ）なので、ここでは係数2が不要
	return
	print '(a,es24.16)', 'PE_CPU =', PE; ONE = 1
	print '(a,es24.16)', 'PE_GPU =', cublas_ddot  ( NP, one(0), 1, DUM(0), 1 ) * eps_th ! 内積
!	return
	print '(a,es24.16)', 'Ua_CPU =', Ua
	print '(a,es24.16)', 'Ua_GPU =', cublas_ddot  ( NP, one, 1, Y(4), 6 )   ! 内積
	print '(a,es24.16)', 'Va_CPU =', Va
	print '(a,es24.16)', 'Va_GPU =', cublas_ddot  ( NP, one, 1, Y(5), 6 )   ! 内積
	print '(a,es24.16)', 'Wa_CPU =', Wa
	print '(a,es24.16)', 'Wa_GPU =', cublas_ddot  ( NP, one, 1, Y(6), 6 )   ! 内積
	print '()'
!
end subroutine Invariance
