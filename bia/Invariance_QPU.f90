subroutine Invariance_QPU( N_p, eps_th, X, Y, Z, U, V, W, KE, PE, Lx, Ly, Lz, Ua, Va, Wa, Xa, Ya, Za )
   IMPLICIT REAL*8( A - H, K - M, O - Z )
   real(08), dimension(N_p) ::  X, Y, Z, U, V, W
   real(16), dimension(N_p) :: QX,QY,QZ,QU,QV,QW
   real(16) QPE, QXi, QXj, QYi, QYj, QZi, QZj
!
   QX = X; QY = Y; QZ = Z; QU = U; QV = V; QW = W
   Xa = sum( QX ); Ua = sum( QU ); Lx = sum( QY * QW - QZ * QV )
   Ya = sum( QY ); Va = sum( QV ); Ly = sum( QZ * QU - QX * QW )
   Za = sum( QZ ); Wa = sum( QW ); Lz = sum( QX * QV - QY * QU )
   KE = sum( QU * QU + QV * QV +QW * QW )
!
   QPE = 0
   do i = 1, N_p - 1   ; QXi = QX(i); QYi = QY(i); QZi = QZ(i)
      do j = i + 1, N_p; QXj = QX(j); QYj = QY(j); QZj = QZ(j)
         QPE = QPE + 1 / sqrt( ( QXi - QXj )**2 + ( QYi - QYj )**2 + (QZi - QZj )**2 )
      end do
   end do
   PE = eps_th * QPE  ! eps_th は、運動方程式の二倍（ 2 b_0_th）なので、ここでは係数2が不要
!
end subroutine Invariance_QPU
