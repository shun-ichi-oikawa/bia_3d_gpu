PROGRAM DIM_3D
	IMPLICIT REAL*8 (A-H,K-M,O-Z)
	integer, parameter :: NL = 5, N_calc_EXP = 0, N_time = 1
	real(8), parameter :: time_end = 1
	real(8), parameter :: RELERR = 0, ABSERR = 1d-16, TH = time_end / N_time
	integer, parameter :: N_P = ( 2 * NL + 1 )**3 + 1, NEQ = 6 * N_P, N_calc = 10**N_calc_EXP
	real(8) YZ(NEQ), Y0(NEQ), YN(NEQ), ESTERR(NEQ), WK1(NEQ), WK2(NEQ), WK3(NEQ), WK4(NEQ,11)
!	common / E_  / eps_th; data eps_th /   3d-7 /
	common / E_  / eps_th; data eps_th / - 3d-7 /
	common / N_P_ZZ / NPP; data    NPP /   N_P  /
	character cN_p * 7, cN_calc_EXP * 1
!
!	print *, 'N_p =', N_p, 'Neq =', Neq
	write(cN_P  ,'(i7.7)') N_P; write(cN_calc_EXP,'(i1.1)') N_calc_EXP
	open (1, file = '../OUT/' // cN_P // '_10^' // cN_calc_EXP // '.bbb', form='unformatted', action='read' )
!  open (1, file = '../OUT/' // cN_P // '_10^' // cN_calc_EXP // '.ccc', form='unformatted', action='read' )
	open (2, file = '../OUT/' // cN_P // '_10^' // cN_calc_EXP // '.DIM', form='unformatted' )
  	open (3, file = '../OUT/DIM_10^' // cN_calc_EXP // '.' // cN_P )
	call cpu_time(ttt0)
	do j_calc = 0, N_calc - 1
		read(1) YZ
!		----------------------------------------------------------------------------------------------------
!		原点に静止している粒子'1'に対して相対運動が放物線軌道となる粒子'2'の位置と速度
!		----------------------------------------------------------------------------------------------------
		Y0(:) = YZ(:)
		INIT = 1; X0 = 0
		do i = 1, N_time; XE = TH * i
!!			call RKF   (     NEQ ,     MOON_3D, X0,XE,Y0,INIT,RELERR,ABSERR,YN,ESTERR,NDE,IER,WK0,WK1,WK2,WK3,WK4)
!!			call RKF_CU(     NEQ ,              X0,XE,Y0,INIT,RELERR,ABSERR,YN,ESTERR,NDE,IER,WK0,WK1,WK2,WK3,WK4)
			call RKF_C (%VAL(NEQ),%VAL(eps_th), X0,XE,Y0,INIT,RELERR,ABSERR,YN,ESTERR,NDE,IER,WK0,WK1,WK2,WK3,WK4)
			write(3,'(es23.16,6es24.16)') XE, YN(1:6)
		end do; write(2) XE, YN(:)
!		write(3,'(f12.6,6es24.16)') time_end, YN(1:6)
		if ( mod( j_calc + 1, max( N_calc / 100, 100 ) ) == 0 ) then; call cpu_time(ttt_)
			print '(i9,a,i9,2(f9.3,a))', j_calc+1, '/', N_calc, ttt_/60, ' min =', ttt_/3600, ' hours'
		end if
	end do
	call cpu_time(ttt_); print '(a,f12.3,a)', 'CPU =', ttt_ - ttt0, ' sec'
	call Invariance( N_p, 2 * eps_th, YZ, KEZ, PEZ, LxZ, LyZ, LzZ, UaZ, VaZ, WaZ, XaZ, YaZ, ZaZ )
	call Invariance( N_p, 2 * eps_th, YN, KEN, PEN, LxN, LyN, LzN, UaN, VaN, WaN, XaN, YaN, ZaN )
	print '(4x,6a24)',                'Ua',      'Va',      'Wa'
	print '(a4,6es24.16)', 'ini:',     UaZ  ,       VaZ,     WaZ
	print '(a4,6es24.16)', 'fin:', UaN      ,   VaN    , WaN
	print '(a4,6es24.16)', 'erA:', UaN-UaZ  ,   VaN-VaZ, WaN-WaZ
	print '(a4,6es24.16)', 'erR:', UaN/UaZ-1, VaN/VaZ-1, WaN/WaZ-1
	print '(4x,6 a24   )',            'Lx'  ,    'Ly'  ,    'Lz'
	print '(a4,6es24.16)', 'ini:',     LxZ  ,     LyZ  ,     LzZ
	print '(a4,6es24.16)', 'fin:', LxN      , LyN      , LzN
	print '(a4,6es24.16)', 'erA:', LxN-LxZ  , LyN-LyZ  , LzN-LzZ
	print '(a4,6es24.16)', 'erR:', LxN/LxZ-1, LyN/LyZ-1, LzN/LzZ-1
	print '(4x,6 a24   )',            'E = K + U'
	print '(a4,6es24.16)', 'ini:',            KEZ+PEZ
	print '(a4,6es24.16)', 'fin:',  KEN+PEN
	print '(a4,6es24.16)', 'erA:', (KEN-KEZ)+(PEN-PEZ)
	print '(a4,6es24.16)', 'erR:', (KEN+PEN)/(KEZ+PEZ)-1
END PROGRAM DIM_3D
include 'Invariance.f90'
