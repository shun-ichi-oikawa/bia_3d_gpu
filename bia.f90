program BIA_3D
!#define DEBUG
	implicit real(8) (a-h, k-m, o-z)
!  character*3 METHOD / 'RKF' /
   character*3 METHOD / 'ALG' /
!  character*3 PLATFM / 'CPU' /
   character*3 PLATFM / 'GPU' /
	integer, parameter :: N_calc_EXP = 0
!	integer, parameter :: N_p =      2, NL =  0, N_end = N_p  ! N_end = N_p のとき、すべての粒子を解く
!	integer, parameter :: N_p =     28, NL =  1, N_end = N_p  ! N_end = N_p のとき、すべての粒子を解く
!	integer, parameter :: N_p =    126, NL =  2, N_end = N_p  ! N_end = N_p のとき、すべての粒子を解く
!	integer, parameter :: N_p =    344, NL =  3, N_end = N_p  ! N_end = N_p のとき、すべての粒子を解く
!	integer, parameter :: N_p =    730, NL =  4, N_end = N_p  ! N_end = 1 のとき一番目の粒子のみを解く
	integer, parameter :: N_p =   1332, NL =  5, N_end = N_p
!	integer, parameter :: N_p =   2198, NL =  6, N_end = N_p
!	integer, parameter :: N_p =   3376, NL =  7, N_end = N_p
!	integer, parameter :: N_p =   4914, NL =  8, N_end = N_p
!	integer, parameter :: N_p =   9262, NL = 10, N_end = N_p
!	integer, parameter :: N_p =  15626, NL = 12, N_end = N_p
!	integer, parameter :: N_p =  19684, NL = 13, N_end = N_p
!	integer, parameter :: N_p =  24390, NL = 14, N_end = N_p
!	integer, parameter :: N_p =  29792, NL = 15, N_end = N_p
!	integer, parameter :: N_p =  29792, NL = 15, N_end =   1
! 	integer, parameter :: N_p =  68922, NL = 20, N_end = N_p
!	integer, parameter :: N_p = 132652, NL = 25, N_end = N_p
!	integer, parameter :: N_p = 226982, NL = 30, N_end = N_p
!	integer, parameter :: N_p =1030302, NL = 50, N_end = N_p
!	integer, parameter :: N_p =1030302, NL = 50, N_end =   1
	integer, parameter :: N_calc = 10**N_calc_EXP
!
	real(8), parameter :: time_end = 1
!	real(8), parameter :: eps_th =   2 * 3d-7 ! 多体RKFの二倍（換算質量の逆数）
	real(8), parameter :: eps_th = - 2 * 3d-7 ! 多体RKFの二倍（換算質量の逆数）
	real(8), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494d+00
	real(8), dimension(3) :: r, g, L, e_xi, e_eta, e_zeta ! 実験室系(X,Y,Z)
	real(8), dimension(4) :: YN, YZ                       ! 相対系(xi, eta) zeta成分なし
	real(8), dimension(N_p) ::  X , Y , Z , U , V , W , dX, dY, dZ, dU, dV, dW
	real(8), dimension(N_p) ::  XF, YF, ZF, UF, VF, WF, RX, RY, RZ, RU, RV, RW
	real(8), dimension(6*N_p) :: RXYZUVW, dXYZUVW

	common / com_time_theta / b_0, g_0, theta_0, tn_0, cs_0, sn_0, tm_end
	common / comij / i, j
	common / const / eps_th_; data eps_th_ / eps_th /
	character cN_p * 7, cN_calc_EXP * 1
!
!	real(16) dum(N_P) ! required ifort or gcc4.6+
!
	ERR = 1d-12; ITR = 18
	print *, 'N_p =', N_p   !, 'Neq =', Neq
!	print *, 'eps_th =', eps_th
	write(cN_P  ,'(i7.7)') N_P; write(cN_calc_EXP,'(i1.1)') N_calc_EXP
	open (1, file = '../OUT/' // cN_P // '_10^' // cN_calc_EXP // '.bbb', form='unformatted', action='read' )
	open (2, file = '../OUT/' // cN_P // '_10^' // cN_calc_EXP // '.DIM', form='unformatted', action='read' )
	open (3, file = '../OUT/' // cN_P // '_10^' // cN_calc_EXP // '.BIA', form='unformatted' )
	open (4, file = '../OUT/' //      'BIA_10^' // cN_calc_EXP // '.' // cN_P )
!
	call cpu_time(ttt0)
	do i_calc = 0, N_calc - 1
!		--	begin 初期条件作成・読み込み
		read (1)      ( X (i_p), Y (i_p), Z (i_p), U (i_p), V (i_p), W (i_p), i_p = 1, N_p )
!		X(2) = 2 * eps_th; Y(2) = 0; Z(2) = 0
!		U(2) = 0         ; V(2) = 1; W(2) = 0
!		X(2) = X(2) * 1000000; V(2) = V(2) / 1000
		call cpu_time(tta)
		call Invariance_GPU ( %VAL(N_p), %VAL(eps_th), X, Y, Z, U, V, W, &
				KEGini, PEGini, LxGini, LyGini, LzGini, UaGini, VaGini, WaGini, XaGini, YaGini, ZaGini )
		call cpu_time(ttb); print '(a,f12.3,a)', ' GPU ', ttb-tta, ' sec'
		call Invariance_CPU (      N_p ,      eps_th , X, Y, Z, U, V, W, &
				KE_ini, PE_ini, Lx_ini, Ly_ini, Lz_ini, Ua_ini, Va_ini, Wa_ini, Xa_ini, Ya_ini, Za_ini )
		call cpu_time(ttc)
		print '(a,5 a24   )', 'method', 'Xa', 'Ya', 'Za'
		print '(a,5es24.16)', 'CUDA: ',        XaGini  ,        YaGini  ,        ZaGini
		print '(a,5es24.16)', 'FORT: ', Xa_ini         , Ya_ini         , Za_ini
		print '(a,5es24.16)', 'Rdif: ', Xa_ini/XaGini-1, Ya_ini/YaGini-1, Za_ini/ZaGini-1
		print '(a,5 a24   )', 'method','KE', 'PE', 'KE + PE'
		print '(a,5es24.16)', 'CUDA: ',         KEGini ,        PEGini  ,                 KEGini+PEGini
		print '(a,5es24.16)', 'FORT: ', KE_ini         , PE_ini         , KE_ini+PE_ini
		print '(a,5es24.16)', 'Rdif: ', KE_ini/KEGini-1, PE_ini/PEGini-1,(KE_ini+PE_ini)/(KEGini+PEGini)-1
		print '(a,f12.3,a)', ' CPU ', ttc-ttb, ' sec'
		call Invariance_QPU (      N_p ,      eps_th , X, Y, Z, U, V, W, &
				KEQini, PEQini, LxQini, LyQini, LzQini, UaQini, VaQini, WaQini, XaQini, YaQini, ZaQini )
		call cpu_time(ttd)
		print '(a,5 a24   )', 'method', 'Xa', 'Ya', 'Za'
		print '(a,5es24.16)', 'CUDA: ',        XaGini  ,        YaGini  ,        ZaGini
		print '(a,5es24.16)', 'FORQ: ', XaQini         , YaQini         , ZaQini
		print '(a,5es24.16)', 'diff: ', XaQini/XaGini-1, YaQini/YaGini-1, ZaQini/ZaGini-1
		print '(a,5 a24   )', 'method','KE', 'PE', 'KE + PE'
		print '(a,5es24.16)', 'CUDA: ',        KEGini  ,        PEGini  ,                  KEGini+PEGini
		print '(a,5es24.16)', 'FORQ: ', KEQini         , PEQini         ,  KEQini+PEQini
		print '(a,5es24.16)', 'diff: ', KEQini/KEGini-1, PEQini/PEGini-1, (KEQini+PEQini)/(KEGini+PEGini)-1
		print '(a,f12.3,a)', ' QPU ', ttd-ttc, ' sec'
!	  	STOP
		print '()'
		print '(a,5 a24   )', 'method', 'Xa', 'Ya', 'Za'
		print '(a,5es24.16)', 'DBLE: ', Xa_ini        ,  Ya_ini         ,  Z__ini
		print '(a,5es24.16)', 'QUAD: ',        XaQini ,         YaQini  ,         ZaQini
		print '(a,5es24.16)', 'diff: ', Xa_ini/XaQini-1, Ya_ini/YaQini-1,  Za_ini/ZaQini-1
		print '(a,5 a24   )', 'method', 'Ua', 'Va', 'Wa'
		print '(a,5es24.16)', 'DBLE: ', Ua_ini        ,  Va_ini         ,  W__ini
		print '(a,5es24.16)', 'QUAD: ',        UaQini ,         VaQini  ,         WaQini
		print '(a,5es24.16)', 'diff: ', Ua_ini/UaQini-1, Va_ini/VaQini-1,  Wa_ini/WaQini-1
		print '(a,5 a24   )', 'method','KE', 'PE', 'KE + PE'
		print '(a,5es24.16)', 'DBLE: ', KE_ini         , PE_ini         ,  KE_ini+PE_ini
		print '(a,5es24.16)', 'QUAD: ',        KEQini  ,        PEQini  ,                  KEQini+PEQini
		print '(a,5es24.16)', 'diff: ', KE_ini/KEQini-1, PE_ini/PEQini-1, (KE_ini+PE_ini)/(KEQini+PEQini)-1
!
		print '(a,5 a24   )', 'method', 'Xa', 'Ya', 'Za'
		print '(a,5es24.16)', 'CUDA: ', XaGini         , YaGini         ,  ZaGini
		print '(a,5es24.16)', 'QUAD: ',        XaQini  ,         YaQini ,         ZaQini
		print '(a,5es24.16)', 'diff: ', XaGini/XaQini-1, YaGini/YaQini-1,  ZaGini/ZaQini-1
		print '(a,5 a24   )', 'method', 'Ua', 'Va', 'Wa'
		print '(a,5es24.16)', 'CUDA: ', UaGini         , VaGini         ,  WaGini
		print '(a,5es24.16)', 'QUAD: ',        UaQini  ,        VaQini  ,         WaQini
		print '(a,5es24.16)', 'diff: ', UaGini/UaQini-1, VaGini/VaQini-1,  WaGini/WaQini-1
		print '(a,5 a24   )', 'method','KE', 'PE', 'KE + PE'
		print '(a,5es24.16)', 'CUDA: ', KEGini         , PEGini         ,  KEGini+PEGini
		print '(a,5es24.16)', 'QUAD: ',        KEQini  ,        PEQini  ,                  KEQini+PEQini
		print '(a,5es24.16)', 'diff: ', KEGini/KEQini-1, PEGini/PEQini-1, (KEGini+PEGini)/(KEQini+PEQini)-1
!
		if ( N_end == N_p ) then                           ! 読み込み：バイナリデータ
			read (2) ttt, ( XF(i_p), YF(i_p), ZF(i_p), UF(i_p), VF(i_p), WF(i_p), i_p = 1, N_p )
		end if
!		-- end 初期条件・読み込み
		dX = 0; dY = 0; dZ = 0; dU = 0; dV = 0; dW = 0
		call cpu_time(ttB1)
		do i = 1, N_end;
			if( PLATFM == 'GPU' ) then
				RX(:) = X(i) - X(:); RY(:) = Y(i) - Y(:); RZ(:) = Z(i) - Z(:)  !  相対位置と速度の配列作成
				RU(:) = U(i) - U(:); RV(:) = V(i) - V(:); RW(:) = W(i) - W(:)
				RXYZUVW(1::6) = RX(:); RXYZUVW(2::6) = RY(:); RXYZUVW(3::6) = RZ(:)
				RXYZUVW(4::6) = RU(:); RXYZUVW(5::6) = RV(:); RXYZUVW(6::6) = RW(:)
				if(METHOD=='RKF') then
					call find_dxyuv_kernel( %VAL(0), %VAL(N_p), %VAL(I), %VAL(eps_th), %VAL(time_end), %VAL(ITR), &
														dX, dY, dZ, dU, dV, dW, RXYZUVW, dXYZUVW)
				else
					call find_dxyuv_kernel( %VAL(1), %VAL(N_p), %VAL(I), %VAL(eps_th), %VAL(time_end), %VAL(ITR), &
														dX, dY, dZ, dU, dV, dW, RXYZUVW, dXYZUVW)
				end if
			else        ! 'CPU'
				do j = i + 1, N_p
					r(1) = X(i) - X(j); r(2) = Y(i) - Y(j); r(3) = Z(i) - Z(j)
					g(1) = U(i) - U(j); g(2) = V(i) - V(j); g(3) = W(i) - W(j)
!  				---------------------------------------------------------------------------------
!  				( xi, eta, zeta )系への変換 (r(1:3), g(1:3)) → ( YZ(1:2), 0, YZ(3:4), 0)
					call xyz_to_xi_eta_zeta( r, g, e_xi, e_eta, e_zeta, YZ ) ! YZ: 相互作用前
!  				----------------------------------------------------------------------------------
					if		( METHOD == 'RKF' ) then      	! Runge-Kutta-Fehlberg
						call find_dXYUV_RKF( eps_th, time_end, YZ, YN, error, IER ) ! 二次元プログラム流用
						if	( abs(error) > ERR ) then
						   write(88,'(i9,2i6,i4,es12.4)') i_calc, i, j, IER, error
						   print    '(i9,2i6,i4,es12.4)', i_calc, i, j, IER, error
						end if
					else	! METHOD == 'ALG'            	! Exact with Newton-Raphson
						r_ij = sqrt( r(1)**2 + r(2)**2 + r(3)**2 )                     ! 相対距離
						E_ij =       g(1)**2 + g(2)**2 + g(3)**2 + 2 * eps_th / r_ij   ! 全エネルギ
						if      (   E_ij >  0 ) then  ! 双曲線	hyperbola
							call find_hyperbola ( eps_th, time_end, YZ, YN, error, ERR,     ITR, IR, IER )
						else if (   E_ij == 0 ) then  ! 放物線	parabola
							call find_parabola  ( eps_th, time_end, YZ, YN, error                        )
						else if (   E_ij <  0 ) then  ! 楕　円	ellipse
							call find_ellipse   ( eps_th, time_end, YZ, YN, error, ERR,     ITR, IR, IER )
						else
							stop                       ! ここには絶対に来ることはない。
						end if
						if ( abs(error) > ERR ) then
							write(99,'(i9,2i6,2i4,2es12.4)') i_calc, i, j, IER, IR, theta_0/pi, error; IER = 0
							print    '(i9,2i6,2i4,2es12.4)', i_calc, i, j, IER, IR, theta_0/pi, error
						end if
					end if
!  				------------------------------------------------------------------------------
					dr_xi =   YN(1) -   YZ(1) ; dr_eta = YN(2) - YZ(2)
					dX_ij =   dr_xi * e_xi(1) + dr_eta * e_eta(1)
					dY_ij =   dr_xi * e_xi(2) + dr_eta * e_eta(2)
					dZ_ij =   dr_xi * e_xi(3) + dr_eta * e_eta(3)
!
					dg_xi =   YN(3) -   YZ(3) ; dg_eta = YN(4) - YZ(4)
					dUijh = ( dg_xi * e_xi(1) + dg_eta * e_eta(1) ) / 2
					dVijh = ( dg_xi * e_xi(2) + dg_eta * e_eta(2) ) / 2
					dWijh = ( dg_xi * e_xi(3) + dg_eta * e_eta(3) ) / 2
!
					dXXXh = ( dX_ij - g(1) * time_end ) / 2   ! 相対位置変化-X
					dYYYh = ( dY_ij - g(2) * time_end ) / 2   ! 相対位置変化-Y
					dZZZh = ( dZ_ij - g(3) * time_end ) / 2   ! 相対位置変化-Z
!
					dX(i)  = dX(i) + dXXXh; dY(i) = dY(i) + dYYYh; dZ(i) = dZ(i) + dZZZh
					dX(j)  = dX(j) - dXXXh; dY(j) = dY(j) - dYYYh; dZ(j) = dZ(j) - dZZZh
!
					dU(i)  = dU(i) + dUijh; dV(i) = dV(i) + dVijh; dW(i) = dW(i) + dWijh
					dU(j)  = dU(j) - dUijh; dV(j) = dV(j) - dVijh; dW(j) = dW(j) - dWijh
				end do ! loop-粒子j
			end if
			dX(i) = dX(i) + U(i) * time_end
			dY(i) = dY(i) + V(i) * time_end
			dZ(i) = dZ(i) + W(i) * time_end
			if ( i == 1 ) then   ! 一個の粒子だけのCPU時間
				call cpu_time(tttt); CPU_1_ = tttt - ttB1
				print '(a,f12.3,a)', 'BIA(1) ' // PLATFM // ' with ' // METHOD // ' =', CPU_1_, ' sec'
			end if
			if ( mod(i,1000) == 0 ) print '(a,$)', '* '
		end do; print *, '' ! loop-粒子i
		X = X + dX; Y = Y + dY; Z = Z + dZ
		U = U + dU; V = V + dV; W = W + dW
!		write(3                     ) time_end, ( X(ii), Y(ii), Z(ii), U(ii), V(ii), W(ii), ii = 1, N_p )
		write(3                     ) time_end, (dX(ii),dY(ii),dZ(ii),dU(ii),dV(ii),dW(ii), ii = 1, N_p )
		write(4,'(es23.16,6es24.16)') time_end, ( X(ii), Y(ii), Z(ii), U(ii), V(ii), W(ii), ii = 1,   1 )
!		if ( mod( i_calc + 1, max( N_calc / 100, 10 ) ) == 0 ) then; call cpu_time(ttt1)
!			print '(i9,a,f9.0,a,f9.6,a)',	i_calc + 1, PLATFM // ' for BIA-' // METHOD, &
!															ttt1-ttt0, ' sec, BIA(1) =', CPU_1_, ' sec'
!		end if
		call cpu_time(ttBN); print '(a,f16.3,a,f9.6)', PLATFM // ' for BIA-' // METHOD, ttBN-ttB1, ' sec'
	end do   ! loop-i_calc
	call cpu_time(ttt2); print '(a,f12.3,a)', PLATFM // ' with ' // METHOD // ' time =', ttt2 - ttt0, ' sec'
!stop
!
	print '(a     24,2a22   )',       'X(1)'       ,'Y(1)'       ,'Z(1)'
	print '(a,f21.15,2f22.15)', 'DIM',        XF(1),        YF(1),        ZF(1)
	print '(a,f21.15,2f22.15)', 'BIA', X(1)        , Y(1)        , Z(1)
	print '(a,f21.15,2f22.15)', 'dif', X(1) - XF(1), Y(1) - YF(1), Z(1) - ZF(1)
	print '(a     24,2a22   )',       'U(1)'       ,'V(1)'       ,'W(1)'
	print '(a,f21.15,2f22.15)', 'DIM',        UF(1),        VF(1),        WF(1)
	print '(a,f21.15,2f22.15)', 'BIA', U(1)        , V(1)        , W(1)
	print '(a,f21.15,2f22.15)', 'dif', U(1) - UF(1), V(1) - VF(1), W(1) - WF(1)
	call cpu_time(ttt3); print '(a,f12.3,a)', 'CPU time =', ttt3 - ttt0, ' sec'
	if ( N_end /= N_p ) stop
!
!  ---------------------------------------------------------------------------------------------------
!
	call Invariance_CPU( N_p, eps_th, XF, YF, ZF, UF, VF, WF, KE_fin, PE_fin, Lx_fin, Ly_fin, Lz_fin, &
!	call Invariance_QPU( N_p, eps_th, XF, YF, ZF, UF, VF, WF, KE_fin, PE_fin, Lx_fin, Ly_fin, Lz_fin, &
	                          Ua_fin, Va_fin, Wa_fin, Xa_fin, Ya_fin, Za_fin )
	call Invariance_CPU( N_p, eps_th, X , Y , Z , U , V , W , KE_end, PE_end, Lx_end, Ly_end, Lz_end, &
!	call Invariance_QPU( N_p, eps_th, X , Y , Z , U , V , W , KE_end, PE_end, Lx_end, Ly_end, Lz_end, &
	                          Ua_end, Va_end, Wa_end, Xa_end, Ya_end, Za_end )
!  ---------------------------------------------------------------------------------------------------
	print '(a)', '------------------------------------------------------------------------------'
	print '(a6,3 a24   )', '    ', 'sum(U)', 'sum(V)', 'sum(W)'
	print '(a6,3es24.16)', 'ini:', Ua_ini, Va_ini, Wa_ini
	print '(a6,3es24.16)', 'fin:', Ua_fin, Va_fin, Wa_fin
	print '(a6,3es24.16)', 'end:', Ua_end, Va_end, Wa_end
!  print *, ''
	print '(a6,3es24.16)', 'err_f:', Ua_fin - Ua_ini, Va_fin - Va_ini, Wa_fin - Wa_ini
	print '(a6,3es24.16)', 'err_e:', Ua_end - Ua_ini, Va_end - Va_ini, Wa_end - Wa_ini
	print '(a6,3es24.16)', 'err_e:',     sum(dU)    ,     sum(dV)    ,     sum(dW)
!  print *, ''
	print '(a6,4 a24   )', '    ', 'sum(Lx)'  , 'sum(Ly)'  , 'sum(Lz)'
	print '(a6,4es24.16)', 'ini:',      Lx_ini,      Ly_ini,      Lz_ini
	print '(a6,4es24.16)', 'fin:',      Lx_fin,      Ly_fin,      Lz_fin
	print '(a6,4es24.16)', 'end:',      Lx_end,      Ly_end,      Lz_end
!  print *, ''
	print '(a6,4es24.16)', 'err_f:', Lx_fin - Lx_ini, Ly_fin - Ly_ini, Lz_fin - Lz_ini
	print '(a6,4es24.16)', 'err_e:', Lx_end - Lx_ini, Ly_end - Ly_ini, Lz_end - Lz_ini
	print '(a6,4 a24   )', '    ','KE + PE'       ,'KE'   ,'PE'
	print '(a6,4es24.16)', 'ini:', KE_ini + PE_ini, KE_ini, PE_ini
	print '(a6,4es24.16)', 'fin:', KE_fin + PE_fin, KE_fin, PE_fin
	print '(a6,4es24.16)', 'end:', KE_end + PE_end, KE_end, PE_end
	print '(a6,4es24.16)', 'err_f:', PE_fin - PE_ini + ( KE_fin - KE_ini )
	print '(a6,4es24.16)', 'err_e:', PE_end - PE_ini + ( KE_end - KE_ini )
	print '(a)', '------------------------------------------------------------------------------'
	print '(2x,3(a13,i9))','X at i =',    maxloc(abs(X-XF)) , 'Y at i =',maxloc(abs(Y-YF)) , 'Z at i =',maxloc(abs(Z-ZF))
	print '(a,f21.15,2f22.15)', 'DIM', XF(maxloc(abs(X-XF))),         YF(maxloc(abs(Y-YF))),         ZF(maxloc(abs(Z-ZF)))
	print '(a,f21.15,2f22.15)', 'BIA', X (maxloc(abs(X-XF))),         Y (maxloc(abs(Y-YF))),         Z (maxloc(abs(Z-ZF)))
	print '(a,f21.15,2f22.15)', 'dif',    maxval(abs(X-XF)) ,            maxval(abs(Y-YF)) ,            maxval(abs(Z-ZF))
	print '(a,f21.15,2f22.15)', 'std',    sqrt( sum((X-XF)**2) / N_p ), &
												     		 sqrt( sum((Y-YF)**2) / N_p ), &
													      sqrt( sum((Z-ZF)**2) / N_p )
	print *, ''
	print '(2x,3(a13,i9))','U at i =',    maxloc(abs(U-UF)) , 'V at i =',maxloc(abs(V-VF)) , 'W at i =',maxloc(abs(W-WF))
	print '(a,f21.15,2f22.15)', 'DIM', UF(maxloc(abs(U-UF))),         VF(maxloc(abs(V-VF))),         WF(maxloc(abs(W-WF)))
	print '(a,f21.15,2f22.15)', 'BIA', U (maxloc(abs(U-UF))),         V (maxloc(abs(V-VF))),         W (maxloc(abs(W-WF)))
	print '(a,f21.15,2f22.15)', 'dif',	  maxval(abs(U-UF)) ,            maxval(abs(V-VF)) ,            maxval(abs(W-WF))
	print '(a,f21.15,2f22.15)', 'std',	  sqrt( sum((U-UF)**2) / N_p ), &
													  sqrt( sum((V-VF)**2) / N_p ), &
													  sqrt( sum((W-WF)**2) / N_p )
	print *, ''
	print '(a,es9.2)', 'std_r =', sqrt( sum((X-XF)**2 + (Y-YF)**2 + (Z-ZF)**2 ) / N_p )
	print '(a,es9.2)', 'std_v =', sqrt( sum((U-UF)**2 + (V-VF)**2 + (W-WF)**2 ) / N_p )
	print '(/a,i6)', 'N_p =', N_p
	call cpu_time(ttt4); print '(a,f9.3,a)', 'BIA ' // PLATFM // ' with ' // METHOD // ' time =', ttt4 - ttt0, ' sec'
end program BIA_3D
!include 'xyz_to_xi_eta_zeta.f90'
subroutine xyz_to_xi_eta_zeta( r, g, e_xi, e_eta, e_zeta, Y )
	implicit real*8 (a-h,o-z)
	real(8), dimension(3) :: r, g, L, e_xi, e_eta, e_zeta ! 実験室系(x,y,z)
	real(8), dimension(4) :: Y                            ! 相対系(xi, eta) zeta成分なし
	gsq  = g(1)**2 + g(2)**2 + g(3)**2; gm = sqrt( gsq )
!	r_g  = r(1) * g(1) + r(2) * g(2) + r(3) * g(3)
!	L(1) = r(2) * g(3) - r(3) * g(2)
!	L(2) = r(3) * g(1) - r(1) * g(3)
!	L(3) = r(1) * g(2) - r(2) * g(1)
	call vector( r, g, L, r_g )

	e_eta (:) = g(:) / gm
	e_zeta(:) = L(:) / sqrt( sum( L**2 ) )
!	e_xi(1) = e_eta(2) * e_zeta(3) - e_eta(3) * e_zeta(2)
!	e_xi(2) = e_eta(3) * e_zeta(1) - e_eta(1) * e_zeta(3)
!	e_xi(3) = e_eta(1) * e_zeta(2) - e_eta(2) * e_zeta(1)
!
	call vector( e_eta, e_zeta, e_xi, zero )
!
	Y(1) = r(1) * e_xi (1) + r(2) * e_xi (2) + r(3) * e_xi (3); Y(3) =  0
	Y(2) = r(1) * e_eta(1) + r(2) * e_eta(2) + r(3) * e_eta(3); Y(4) = gm
!
	contains
	subroutine vector( a, b, crs, dot )
		real(8)  a(3),  b(3),crs(3),  dot
		crs(1) = a(2) * b(3) - a(3) * b(2)
		crs(2) = a(3) * b(1) - a(1) * b(3)
		crs(3) = a(1) * b(2) - a(2) * b(1)
		dot    = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)
	end subroutine vector
end subroutine xyz_to_xi_eta_zeta
!
include 'Invariance_CPU.f90'
include 'Invariance_QPU.f90'
include 'find_dXYUV_RKF.f90'
include 'find_hyperbola.f90'
include 'find_parabola.f90'
include 'find_ellipse.f90'
include 'find_dXYUV4.f90'
