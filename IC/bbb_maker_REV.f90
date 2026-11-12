PROGRAM bbb_maker
   IMPLICIT REAL(8) (A-H,O-Z)
!  integer, PARAMETER :: N_L =  0, N_calc_EXP = 0
!  integer, PARAMETER :: N_L = 10, N_calc_EXP = 0
   integer, PARAMETER :: N_L = 15, N_calc_EXP = 0
!  integer, PARAMETER :: N_L = 12, N_calc_EXP = 0
!  integer, PARAMETER :: N_L = 20, N_calc_EXP = 0
!  integer, PARAMETER :: N_L = 25, N_calc_EXP = 0
!  integer, PARAMETER :: N_L = 30, N_calc_EXP = 0
!  integer, PARAMETER :: N_L = 50, N_calc_EXP = 0
   integer, PARAMETER :: N_P = ( 2 * N_L + 1 )**3 + 1, NEQ = 6 * N_P, N_calc = 10**N_calc_EXP
   real(8), DIMENSION (NEQ) :: YZ, DD
   real(8), dimension (-N_L:N_L,-N_L:N_L,-N_L:N_L) :: x, y, z, u, v, w  ! 場の粒子のみ、テスト粒子を含まない！
   character cN_p * 7, cN_calc_EXP * 1
   print *, 'N_P =', N_P
   write(cN_P  ,'(i7.7)') N_P; write(cN_calc_EXP,'(i1.1)') N_calc_EXP
   open (1, file = '../OUT/' // cN_P // '_10^' // cN_calc_EXP // '.bbb', form='unformatted' ) ! 原点近くから順番に保存
!  call random_seed()
   do j_calc  = 0, N_calc - 1; call random_number(DD); DD(:) = ( 2 * DD(:) - 1 ) / 2
      YZ(1:6) = 0 ! 原点に静止   　! 新・旧粒子番号=1
      do iz = - N_L, N_L; do iy = - N_L, N_L; do ix = - N_L, N_L
         i = ( 2 * N_L + 1 ) * ( ( 2 * N_L + 1 ) * ( iz + N_L ) + iy + N_L ) + ix + N_L + 1  ! 旧粒子番号: 2--N_p
         i_= 6 * i
         YZ( 1 + i_ ) = ix + DD( 1 + i_ ); YZ( 2 + i_ ) = iy + DD( 2 + i_ ); YZ( 3 + i_ ) = iz + DD( 3 + i_ )
         YZ( 4 + i_ ) =      DD( 4 + i_ ); YZ( 5 + i_ ) =      DD( 5 + i_ ); YZ( 6 + i_ ) =      DD( 6 + i_ )
         x (ix,iy,iz) = ix + DD( 1 + i_ ); y (ix,iy,iz) = iy + DD( 2 + i_ ); z (ix,iy,iz) = iz + DD( 3 + i_ )
         u (ix,iy,iz) =      DD( 4 + i_ ); v (ix,iy,iz) =      DD( 5 + i_ ); w (ix,iy,iz) =      DD( 6 + i_ )
      end do; end do; end do
      YZ(7:12)= DD(7:12)   ! ix=iy=0なので、位置もOK
      j_p = 2
!
      do i_L = 1, N_L
!
         do iz = - i_L, i_L, 2 * i_L
            do iy =   - i_L, i_L
               do ix =   - i_L, i_L
                  j_p_ = 6 * j_p
                  YZ(j_p_+1) = DD(j_p_+1) + ix; YZ(j_p_+2) = DD(j_p_+2) + iy; YZ(j_p_+3) = DD(j_p_+3) + iz
                  YZ(j_p_+4) = DD(j_p_+4)     ; YZ(j_p_+5) = DD(j_p_+5)     ; YZ(j_p_+6) = DD(j_p_+6)
                  j_p = j_p + 1
!                 print '(i8,3i4)', j_p, ix, iy, iz
               end do
            end do
         end do
!
         do iy = - i_L, i_L, 2 * i_L
            do iz = 1 - i_L, i_L - 1
               do ix =   - i_L, i_L
                  j_p_ = 6 * j_p
                  YZ(j_p_+1) = DD(j_p_+1) + ix; YZ(j_p_+2) = DD(j_p_+2) + iy; YZ(j_p_+3) = DD(j_p_+3) + iz
                  YZ(j_p_+4) = DD(j_p_+4)     ; YZ(j_p_+5) = DD(j_p_+5)     ; YZ(j_p_+6) = DD(j_p_+6)
                  j_p = j_p + 1
!                 print '(i8,3i4)', j_p, ix, iy, iz
               end do
            end do
         end do
!
         do ix = - i_L, i_L, 2 * i_L
            do iz = 1 - i_L, i_L - 1
               do iy = 1 - i_L, i_L - 1
                  j_p_ = 6 * j_p
                  YZ(j_p_+1) = DD(j_p_+1) + ix; YZ(j_p_+2) = DD(j_p_+2) + iy; YZ(j_p_+3) = DD(j_p_+3) + iz
                  YZ(j_p_+4) = DD(j_p_+4)     ; YZ(j_p_+5) = DD(j_p_+5)     ; YZ(j_p_+6) = DD(j_p_+6)
                  j_p = j_p + 1
!                 print '(i8,3i4)', j_p, ix, iy, iz
               end do
            end do
         end do
!
!        print *, i_L, j_p
      end do
      write(1) YZ
!     if ( mod( i_calc + 1, max( N_calc / 100, 100  )  ) == 0 ) print *, i_calc, N_P, j_p
   end do
   call cpu_time(ttt); print '(a,f9.3,a)', 'CPU =', ttt, ' sec'
END PROGRAM bbb_maker
