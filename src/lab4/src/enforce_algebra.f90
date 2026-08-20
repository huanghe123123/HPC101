
!-----------------------------------------------------------------------------
!
! remove the trace of Aij
! trace-free Aij and enforce the determinant of bssn metric to one
!-----------------------------------------------------------------------------

  subroutine enforce_ag(ex,  dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                             Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz)
  implicit none

!~~~~~~> Input parameters:

  integer,                              intent(in)    :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: Axx,Axy,Axz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: Ayy,Ayz,Azz

!~~~~~~~> Local variable:
  
  real*8, dimension(ex(1),ex(2),ex(3)) :: trA,detg
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz 
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz,gupyy,gupyz,gupzz
  real*8, parameter :: F1o3 = 1.D0 / 3.D0, ONE = 1.D0, TWO = 2.D0

!~~~~~~>

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE

  detg =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
          gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / detg
  gupxy = - ( gxy * gzz - gyz * gxz ) / detg
  gupxz =   ( gxy * gyz - gyy * gxz ) / detg
  gupyy =   ( gxx * gzz - gxz * gxz ) / detg
  gupyz = - ( gxx * gyz - gxy * gxz ) / detg
  gupzz =   ( gxx * gyy - gxy * gxy ) / detg

  trA =         gupxx * Axx + gupyy * Ayy + gupzz * Azz &
       + TWO * (gupxy * Axy + gupxz * Axz + gupyz * Ayz)

  Axx = Axx - F1o3 * gxx * trA
  Axy = Axy - F1o3 * gxy * trA
  Axz = Axz - F1o3 * gxz * trA
  Ayy = Ayy - F1o3 * gyy * trA
  Ayz = Ayz - F1o3 * gyz * trA
  Azz = Azz - F1o3 * gzz * trA

  detg = ONE / ( detg ** F1o3 ) 
  
  gxx = gxx * detg
  gxy = gxy * detg
  gxz = gxz * detg
  gyy = gyy * detg
  gyz = gyz * detg
  gzz = gzz * detg

  dxx = gxx - ONE
  dyy = gyy - ONE
  dzz = gzz - ONE

  return

  end subroutine enforce_ag

!----------------------------------------------------------------------------------  
! swap the turn of a and g
!----------------------------------------------------------------------------------
  subroutine enforce_ga(ex,  dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                             Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz)
  implicit none

!~~~~~~> Input parameters:

  integer,                              intent(in)    :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: Axx,Axy,Axz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: Ayy,Ayz,Azz

!~~~~~~~> Local variable:
  
  real*8, dimension(ex(1),ex(2),ex(3)) :: trA
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz 
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz,gupyy,gupyz,gupzz
  real*8, parameter :: F1o3 = 1.D0 / 3.D0, ONE = 1.D0, TWO = 2.D0

!~~~~~~>

  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
! for g
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz

  gupzz = ONE / ( gupzz ** F1o3 ) 
  
  gxx = gxx * gupzz
  gxy = gxy * gupzz
  gxz = gxz * gupzz
  gyy = gyy * gupzz
  gyz = gyz * gupzz
  gzz = gzz * gupzz

  dxx = gxx - ONE
  dyy = gyy - ONE
  dzz = gzz - ONE
! for A  

  gupxx =   ( gyy * gzz - gyz * gyz )
  gupxy = - ( gxy * gzz - gyz * gxz )
  gupxz =   ( gxy * gyz - gyy * gxz )
  gupyy =   ( gxx * gzz - gxz * gxz )
  gupyz = - ( gxx * gyz - gxy * gxz )
  gupzz =   ( gxx * gyy - gxy * gxy )

  trA =         gupxx * Axx + gupyy * Ayy + gupzz * Azz &
       + TWO * (gupxy * Axy + gupxz * Axz + gupyz * Ayz)

  Axx = Axx - F1o3 * gxx * trA
  Axy = Axy - F1o3 * gxy * trA
  Axz = Axz - F1o3 * gxz * trA
  Ayy = Ayy - F1o3 * gyy * trA
  Ayz = Ayz - F1o3 * gyz * trA
  Azz = Azz - F1o3 * gzz * trA

  return

  end subroutine enforce_ga

! Pointwise form of enforce_ga.  The array-syntax version above materializes
! several full-block temporaries (gxx/gyy/gzz, inverse-metric cofactors and
! trA).  This variant preserves the same algebraic sequence at each point
! while keeping those values in scalar temporaries.
  subroutine enforce_ga_pointwise(ex,  dxx,  gxy,  gxz,  dyy,  gyz,  dzz, &
                                  Axx,  Axy,  Axz,  Ayy,  Ayz,  Azz)
  implicit none

  integer, intent(in) :: ex(1:3)
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)), intent(inout) :: Axx,Axy,Axz,Ayy,Ayz,Azz

  integer :: i,j,k
  real*8 :: gxxs,gyys,gzzs,detg,scale
  real*8 :: gupxxs,gupxys,gupxzs,gupyys,gupyzs,gupzzs,trAs
  real*8, parameter :: F1o3 = 1.D0 / 3.D0, ONE = 1.D0, TWO = 2.D0

  do k = 1, ex(3)
    do j = 1, ex(2)
      do i = 1, ex(1)
        gxxs = dxx(i,j,k) + ONE
        gyys = dyy(i,j,k) + ONE
        gzzs = dzz(i,j,k) + ONE

        detg = gxxs * gyys * gzzs + gxy(i,j,k) * gyz(i,j,k) * gxz(i,j,k) + &
               gxz(i,j,k) * gxy(i,j,k) * gyz(i,j,k) - &
               gxz(i,j,k) * gyys * gxz(i,j,k) - &
               gxy(i,j,k) * gxy(i,j,k) * gzzs - &
               gxxs * gyz(i,j,k) * gyz(i,j,k)
        scale = ONE / (detg ** F1o3)

        gxxs = gxxs * scale
        gxy(i,j,k) = gxy(i,j,k) * scale
        gxz(i,j,k) = gxz(i,j,k) * scale
        gyys = gyys * scale
        gyz(i,j,k) = gyz(i,j,k) * scale
        gzzs = gzzs * scale

        dxx(i,j,k) = gxxs - ONE
        dyy(i,j,k) = gyys - ONE
        dzz(i,j,k) = gzzs - ONE

        gupxxs = gyys * gzzs - gyz(i,j,k) * gyz(i,j,k)
        gupxys = -(gxy(i,j,k) * gzzs - gyz(i,j,k) * gxz(i,j,k))
        gupxzs =  (gxy(i,j,k) * gyz(i,j,k) - gyys * gxz(i,j,k))
        gupyys =  gxxs * gzzs - gxz(i,j,k) * gxz(i,j,k)
        gupyzs = -(gxxs * gyz(i,j,k) - gxy(i,j,k) * gxz(i,j,k))
        gupzzs =  gxxs * gyys - gxy(i,j,k) * gxy(i,j,k)

        trAs = gupxxs * Axx(i,j,k) + gupyys * Ayy(i,j,k) + gupzzs * Azz(i,j,k) + &
               TWO * (gupxys * Axy(i,j,k) + gupxzs * Axz(i,j,k) + gupyzs * Ayz(i,j,k))

        Axx(i,j,k) = Axx(i,j,k) - F1o3 * gxxs * trAs
        Axy(i,j,k) = Axy(i,j,k) - F1o3 * gxy(i,j,k) * trAs
        Axz(i,j,k) = Axz(i,j,k) - F1o3 * gxz(i,j,k) * trAs
        Ayy(i,j,k) = Ayy(i,j,k) - F1o3 * gyys * trAs
        Ayz(i,j,k) = Ayz(i,j,k) - F1o3 * gyz(i,j,k) * trAs
        Azz(i,j,k) = Azz(i,j,k) - F1o3 * gzzs * trAs
      end do
    end do
  end do

  return
  end subroutine enforce_ga_pointwise
