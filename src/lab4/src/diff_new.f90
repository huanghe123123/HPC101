

#include "macrodef.fh"

#ifdef AMSS_RHS_HALO
#define FH(i,j,k) halo_get(i,j,k)
#else
#define FH(i,j,k) fh(i,j,k)
#endif


! ---------------------------------------------------------------------------
! Reusable stencil workspace buffer (optimization: avoid per-call malloc).
! fh is declared automatic in each subroutine with lower bound -1; gfortran
! heap-allocates it on every call. This module-level target buffer replaces
! them, lazily (re)allocated to the largest ex seen. Single-threaded (OMP=1).
! ---------------------------------------------------------------------------
module diff_buffers
  implicit none
  save
  real*8, target, allocatable :: fh_buf(:,:,:)
  integer :: fh_ex(3) = 0
contains
  subroutine fh_ensure(ex)
    integer, intent(in) :: ex(3)
    if (.not. allocated(fh_buf) .or. any(ex > fh_ex)) then
      if (allocated(fh_buf)) deallocate(fh_buf)
      allocate(fh_buf(-1:ex(1), -1:ex(2), -1:ex(3)))
      fh_ex = ex
    end if
  end subroutine fh_ensure
end module diff_buffers

#ifdef AMSS_RHS_HALO
! ---------------------------------------------------------------------------
! Equatorial-symmetry halo view.  MPI ghost cells are already present in the
! original field.  Only the negative-z planes needed by the stencil are
! materialized here; positive indices read the original field directly.
! Symmetry>1 keeps the legacy full halo as a correctness fallback.
! ---------------------------------------------------------------------------
module rhs_halo
  implicit none
  save
  real*8, pointer :: source(:,:,:) => null()
  real*8, allocatable, target :: halo2(:,:,:), halo3(:,:,:)
  real*8, allocatable, target :: legacy(:,:,:)
  integer :: halo_ex(3) = 0
  integer :: halo_ord = 0
  logical :: halo_active = .false.
  logical :: legacy_active = .false.
  real*8 :: halo_soa(3) = 1.d0
  external :: symmetry_bd
contains
  subroutine halo_bind(f, ex, SoA, symmetry, ord, kmin)
    real*8, target, intent(in) :: f(:,:,:)
    integer, intent(in) :: ex(3), symmetry, ord, kmin
    real*8, intent(in) :: SoA(3)
    integer :: previous_ord
    previous_ord = halo_ord
    source => f
    halo_soa = SoA
    halo_ord = ord
    halo_active = .false.
    legacy_active = .false.

    if (symmetry > 1) then
      if (.not. allocated(legacy) .or. any(halo_ex /= ex) .or. previous_ord /= ord) then
        if (allocated(legacy)) deallocate(legacy)
        allocate(legacy(-ord+1:ex(1), -ord+1:ex(2), -ord+1:ex(3)))
      end if
      call symmetry_bd(ord, ex, f, legacy, SoA)
      legacy_active = .true.
      halo_ex = ex
      return
    end if

    if (any(halo_ex /= ex)) then
      if (allocated(halo2)) deallocate(halo2)
      if (allocated(halo3)) deallocate(halo3)
      halo_ex = ex
    end if
    if (ord == 2) then
      if (.not. allocated(halo2)) allocate(halo2(ex(1), ex(2), 2))
    else
      if (.not. allocated(halo3)) allocate(halo3(ex(1), ex(2), 3))
    end if

    if (kmin < 1 .and. symmetry > 0) then
      halo_active = .true.
      if (ord == 2) then
        halo2(:,:,1) = f(:,:,1) * SoA(3)
        halo2(:,:,2) = f(:,:,2) * SoA(3)
      else
        halo3(:,:,1) = f(:,:,1) * SoA(3)
        halo3(:,:,2) = f(:,:,2) * SoA(3)
        halo3(:,:,3) = f(:,:,3) * SoA(3)
      end if
    end if
  end subroutine halo_bind

  real*8 function halo_get(i, j, k) result(value)
    integer, intent(in) :: i, j, k
    if (legacy_active) then
      value = legacy(i,j,k)
    else if (k >= 1 .or. .not. halo_active) then
      value = source(i,j,k)
    else if (halo_ord == 2) then
      value = halo2(i,j,1-k)
    else
      value = halo3(i,j,1-k)
    end if
  end function halo_get
end module rhs_halo
#endif

#ifdef AMSS_RHS_POINTWISE
! ---------------------------------------------------------------------------
! Point-local derivative helpers used by the POINTWISE RHS path.  The
! production input uses equatorial symmetry, so negative z indices are
! reflected directly and MPI ghost values remain ordinary positive indices.
! ---------------------------------------------------------------------------
module point_derivs
  implicit none
contains
  real*8 function point_value(f, i, j, k, SoA, symmetry) result(v)
    real*8, intent(in) :: f(:,:,:), SoA(3)
    integer, intent(in) :: i, j, k, symmetry
    v = point_value3(f, i, j, k, SoA(1), SoA(2), SoA(3), symmetry)
  end function point_value

  ! Scalar-sign form used by the preprocessor-expanded pointwise stencils.
  ! Keeping the signs scalar avoids constructing a temporary (/s1,s2,s3/)
  ! for every neighboring value in the fused loop.
  real*8 function point_value3(f, i, j, k, s1, s2, s3, symmetry) result(v)
    real*8, intent(in) :: f(:,:,:), s1, s2, s3
    integer, intent(in) :: i, j, k, symmetry
    integer :: ii, jj, kk
    real*8 :: sign
    ii = i
    jj = j
    kk = k
    sign = 1.d0
    if (ii < 1) then
      ii = 1-ii
      sign = sign*s1
    end if
    if (jj < 1) then
      jj = 1-jj
      sign = sign*s2
    end if
    if (kk < 1) then
      kk = 1-kk
      sign = sign*s3
    end if
    v = sign*f(ii,jj,kk)
  end function point_value3

  subroutine point_d1(ex, X, Y, Z, f, i, j, k, SoA, symmetry, fx, fy, fz)
    integer, intent(in) :: ex(3), i, j, k, symmetry
    real*8, intent(in) :: X(:), Y(:), Z(:), f(:,:,:), SoA(3)
    real*8, intent(out) :: fx, fy, fz
    integer :: imin, jmin, kmin
    real*8 :: dx, dy, dz, d12x, d12y, d12z, d2x, d2y, d2z
    dx = X(2)-X(1); dy = Y(2)-Y(1); dz = Z(2)-Z(1)
    d12x = 1.d0/(12.d0*dx); d12y = 1.d0/(12.d0*dy); d12z = 1.d0/(12.d0*dz)
    d2x = 0.5d0/dx; d2y = 0.5d0/dy; d2z = 0.5d0/dz
    imin = 1; jmin = 1; kmin = 1
    if (symmetry > 0 .and. abs(Z(1)) < dz) kmin = -1
    if (symmetry > 1 .and. abs(X(1)) < dx) imin = -1
    if (symmetry > 1 .and. abs(Y(1)) < dy) jmin = -1
    fx = 0.d0; fy = 0.d0; fz = 0.d0
    if (i+2 <= ex(1) .and. i-2 >= imin .and. j+2 <= ex(2) .and. j-2 >= jmin .and. k+2 <= ex(3) .and. k-2 >= kmin) then
      fx = d12x*(point_value(f,i-2,j,k,SoA,symmetry)-8.d0*point_value(f,i-1,j,k,SoA,symmetry)+8.d0*point_value(f,i+1,j,k,SoA,symmetry)-point_value(f,i+2,j,k,SoA,symmetry))
      fy = d12y*(point_value(f,i,j-2,k,SoA,symmetry)-8.d0*point_value(f,i,j-1,k,SoA,symmetry)+8.d0*point_value(f,i,j+1,k,SoA,symmetry)-point_value(f,i,j+2,k,SoA,symmetry))
      fz = d12z*(point_value(f,i,j,k-2,SoA,symmetry)-8.d0*point_value(f,i,j,k-1,SoA,symmetry)+8.d0*point_value(f,i,j,k+1,SoA,symmetry)-point_value(f,i,j,k+2,SoA,symmetry))
    else if (i+1 <= ex(1) .and. i-1 >= imin .and. j+1 <= ex(2) .and. j-1 >= jmin .and. k+1 <= ex(3) .and. k-1 >= kmin) then
      fx = d2x*(-point_value(f,i-1,j,k,SoA,symmetry)+point_value(f,i+1,j,k,SoA,symmetry))
      fy = d2y*(-point_value(f,i,j-1,k,SoA,symmetry)+point_value(f,i,j+1,k,SoA,symmetry))
      fz = d2z*(-point_value(f,i,j,k-1,SoA,symmetry)+point_value(f,i,j,k+1,SoA,symmetry))
    end if
  end subroutine point_d1

  subroutine point_d2(ex, X, Y, Z, f, i, j, k, SoA, symmetry, fxx, fxy, fxz, fyy, fyz, fzz)
    integer, intent(in) :: ex(3), i, j, k, symmetry
    real*8, intent(in) :: X(:), Y(:), Z(:), f(:,:,:), SoA(3)
    real*8, intent(out) :: fxx, fxy, fxz, fyy, fyz, fzz
    integer :: imin, jmin, kmin
    real*8 :: dx, dy, dz, sxx, syy, szz, fxxc, fxyc, fxzc, fyyc, fyzc, fzzc
    dx = X(2)-X(1); dy = Y(2)-Y(1); dz = Z(2)-Z(1)
    sxx = 1.d0/(dx*dx); syy = 1.d0/(dy*dy); szz = 1.d0/(dz*dz)
    fxxc = 1.d0/(12.d0*dx*dx); fyyc = 1.d0/(12.d0*dy*dy); fzzc = 1.d0/(12.d0*dz*dz)
    fxyc = 1.d0/(144.d0*dx*dy); fxzc = 1.d0/(144.d0*dx*dz); fyzc = 1.d0/(144.d0*dy*dz)
    imin = 1; jmin = 1; kmin = 1
    if (symmetry > 0 .and. abs(Z(1)) < dz) kmin = -1
    if (symmetry > 1 .and. abs(X(1)) < dx) imin = -1
    if (symmetry > 1 .and. abs(Y(1)) < dy) jmin = -1
    fxx = 0.d0; fxy = 0.d0; fxz = 0.d0; fyy = 0.d0; fyz = 0.d0; fzz = 0.d0
    if (i+2 <= ex(1) .and. i-2 >= imin .and. j+2 <= ex(2) .and. j-2 >= jmin .and. k+2 <= ex(3) .and. k-2 >= kmin) then
      fxx = fxxc*(-point_value(f,i-2,j,k,SoA,symmetry)+16.d0*point_value(f,i-1,j,k,SoA,symmetry)-30.d0*point_value(f,i,j,k,SoA,symmetry)+16.d0*point_value(f,i+1,j,k,SoA,symmetry)-point_value(f,i+2,j,k,SoA,symmetry))
      fyy = fyyc*(-point_value(f,i,j-2,k,SoA,symmetry)+16.d0*point_value(f,i,j-1,k,SoA,symmetry)-30.d0*point_value(f,i,j,k,SoA,symmetry)+16.d0*point_value(f,i,j+1,k,SoA,symmetry)-point_value(f,i,j+2,k,SoA,symmetry))
      fzz = fzzc*(-point_value(f,i,j,k-2,SoA,symmetry)+16.d0*point_value(f,i,j,k-1,SoA,symmetry)-30.d0*point_value(f,i,j,k,SoA,symmetry)+16.d0*point_value(f,i,j,k+1,SoA,symmetry)-point_value(f,i,j,k+2,SoA,symmetry))
      fxy = fxyc*((point_value(f,i-2,j-2,k,SoA,symmetry)-8.d0*point_value(f,i-1,j-2,k,SoA,symmetry)+8.d0*point_value(f,i+1,j-2,k,SoA,symmetry)-point_value(f,i+2,j-2,k,SoA,symmetry))-8.d0*(point_value(f,i-2,j-1,k,SoA,symmetry)-8.d0*point_value(f,i-1,j-1,k,SoA,symmetry)+8.d0*point_value(f,i+1,j-1,k,SoA,symmetry)-point_value(f,i+2,j-1,k,SoA,symmetry))+8.d0*(point_value(f,i-2,j+1,k,SoA,symmetry)-8.d0*point_value(f,i-1,j+1,k,SoA,symmetry)+8.d0*point_value(f,i+1,j+1,k,SoA,symmetry)-point_value(f,i+2,j+1,k,SoA,symmetry))-(point_value(f,i-2,j+2,k,SoA,symmetry)-8.d0*point_value(f,i-1,j+2,k,SoA,symmetry)+8.d0*point_value(f,i+1,j+2,k,SoA,symmetry)-point_value(f,i+2,j+2,k,SoA,symmetry)))
      fxz = fxzc*((point_value(f,i-2,j,k-2,SoA,symmetry)-8.d0*point_value(f,i-1,j,k-2,SoA,symmetry)+8.d0*point_value(f,i+1,j,k-2,SoA,symmetry)-point_value(f,i+2,j,k-2,SoA,symmetry))-8.d0*(point_value(f,i-2,j,k-1,SoA,symmetry)-8.d0*point_value(f,i-1,j,k-1,SoA,symmetry)+8.d0*point_value(f,i+1,j,k-1,SoA,symmetry)-point_value(f,i+2,j,k-1,SoA,symmetry))+8.d0*(point_value(f,i-2,j,k+1,SoA,symmetry)-8.d0*point_value(f,i-1,j,k+1,SoA,symmetry)+8.d0*point_value(f,i+1,j,k+1,SoA,symmetry)-point_value(f,i+2,j,k+1,SoA,symmetry))-(point_value(f,i-2,j,k+2,SoA,symmetry)-8.d0*point_value(f,i-1,j,k+2,SoA,symmetry)+8.d0*point_value(f,i+1,j,k+2,SoA,symmetry)-point_value(f,i+2,j,k+2,SoA,symmetry)))
      fyz = fyzc*((point_value(f,i,j-2,k-2,SoA,symmetry)-8.d0*point_value(f,i,j-1,k-2,SoA,symmetry)+8.d0*point_value(f,i,j+1,k-2,SoA,symmetry)-point_value(f,i,j+2,k-2,SoA,symmetry))-8.d0*(point_value(f,i,j-2,k-1,SoA,symmetry)-8.d0*point_value(f,i,j-1,k-1,SoA,symmetry)+8.d0*point_value(f,i,j+1,k-1,SoA,symmetry)-point_value(f,i,j+2,k-1,SoA,symmetry))+8.d0*(point_value(f,i,j-2,k+1,SoA,symmetry)-8.d0*point_value(f,i,j-1,k+1,SoA,symmetry)+8.d0*point_value(f,i,j+1,k+1,SoA,symmetry)-point_value(f,i,j+2,k+1,SoA,symmetry))-(point_value(f,i,j-2,k+2,SoA,symmetry)-8.d0*point_value(f,i,j-1,k+2,SoA,symmetry)+8.d0*point_value(f,i,j+1,k+2,SoA,symmetry)-point_value(f,i,j+2,k+2,SoA,symmetry)))
    else if (i+1 <= ex(1) .and. i-1 >= imin .and. j+1 <= ex(2) .and. j-1 >= jmin .and. k+1 <= ex(3) .and. k-1 >= kmin) then
      fxx = sxx*(point_value(f,i-1,j,k,SoA,symmetry)-2.d0*point_value(f,i,j,k,SoA,symmetry)+point_value(f,i+1,j,k,SoA,symmetry))
      fyy = syy*(point_value(f,i,j-1,k,SoA,symmetry)-2.d0*point_value(f,i,j,k,SoA,symmetry)+point_value(f,i,j+1,k,SoA,symmetry))
      fzz = szz*(point_value(f,i,j,k-1,SoA,symmetry)-2.d0*point_value(f,i,j,k,SoA,symmetry)+point_value(f,i,j,k+1,SoA,symmetry))
      fxy = 0.25d0/(dx*dy)*(point_value(f,i-1,j-1,k,SoA,symmetry)-point_value(f,i+1,j-1,k,SoA,symmetry)-point_value(f,i-1,j+1,k,SoA,symmetry)+point_value(f,i+1,j+1,k,SoA,symmetry))
      fxz = 0.25d0/(dx*dz)*(point_value(f,i-1,j,k-1,SoA,symmetry)-point_value(f,i+1,j,k-1,SoA,symmetry)-point_value(f,i-1,j,k+1,SoA,symmetry)+point_value(f,i+1,j,k+1,SoA,symmetry))
      fyz = 0.25d0/(dy*dz)*(point_value(f,i,j-1,k-1,SoA,symmetry)-point_value(f,i,j+1,k-1,SoA,symmetry)-point_value(f,i,j-1,k+1,SoA,symmetry)+point_value(f,i,j+1,k+1,SoA,symmetry))
    end if
  end subroutine point_d2
end module point_derivs
#endif


! we need only distinguish different finite difference order
! Vertex or Cell is distinguished in routine symmetry_bd which locates in
! file "fmisc.f90"

! fourth order code

!-----------------------------------------------------------------------------------------------------------------
!
! General first derivatives of 4_th oder accurate
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
!
!-----------------------------------------------------------------------------------------------------------------

  subroutine fderivs(ex,f,fx,fy,fz,X,Y,Z,SYM1,SYM2,SYM3,symmetry,onoff)
  use diff_buffers
#ifdef AMSS_RHS_HALO
  use rhs_halo
#endif
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff
  real*8,  dimension(ex(1),ex(2),ex(3)), target, intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fx,fy,fz
  real*8,                                intent(in) :: X(ex(1)),Y(ex(2)),Z(ex(3))
  real*8,                                intent(in ):: SYM1,SYM2,SYM3

!~~~~~~ other variables

  real*8 :: dX,dY,dZ
  real*8, pointer :: fh(:,:,:)
  real*8, dimension(3) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k,klo,khi
  real*8 :: d12dx,d12dy,d12dz,d2dx,d2dy,d2dz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -1
  if(Symmetry > EQ_SYMM .and. dabs(X(1)) < dX) imin = -1
  if(Symmetry > EQ_SYMM .and. dabs(Y(1)) < dY) jmin = -1

  SoA(1) = SYM1
  SoA(2) = SYM2
  SoA(3) = SYM3

#ifdef AMSS_RHS_HALO
  call halo_bind(f, ex, SoA, symmetry, 2, kmin)
#else
  call fh_ensure(ex)
  fh(-1:ex(1),-1:ex(2),-1:ex(3)) => fh_buf
  call symmetry_bd(2,ex,f,fh,SoA)
#endif

  d12dx = ONE/F12/dX
  d12dy = ONE/F12/dY
  d12dz = ONE/F12/dZ

  d2dx = ONE/TWO/dX
  d2dy = ONE/TWO/dY
  d2dz = ONE/TWO/dZ

  fx = ZEO
  fy = ZEO
  fz = ZEO

#ifdef AMSS_RHS_HALO
  if (legacy_active) then
    klo = 1
    khi = ex(3)-1
#define FVAL(i,j,k) halo_get(i,j,k)
#include "fderivs_loop.fh"
#undef FVAL
  else if (halo_active) then
    klo = 1
    khi = min(2, ex(3)-1)
#define FVAL(i,j,k) halo_get(i,j,k)
#include "fderivs_loop.fh"
#undef FVAL
    klo = 3
    khi = ex(3)-1
#define FVAL(i,j,k) f(i,j,k)
#include "fderivs_loop.fh"
#undef FVAL
  else
    klo = 1
    khi = ex(3)-1
#define FVAL(i,j,k) f(i,j,k)
#include "fderivs_loop.fh"
#undef FVAL
  end if
#else
  klo = 1
  khi = ex(3)-1
#define FVAL(i,j,k) fh(i,j,k)
#include "fderivs_loop.fh"
#undef FVAL
#endif

  return

  end subroutine fderivs
!-----------------------------------------------------------------------------
!
! single derivatives dx
!
!-----------------------------------------------------------------------------
  subroutine fdx(ex,f,fx,X,SYM1,symmetry,onoff)
  use diff_buffers
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fx
  real*8,                                intent(in ):: X(ex(1)),SYM1

!~~~~~~ other variables

  real*8 :: dX
  real*8, pointer :: fh(:,:,:)
  real*8, dimension(3) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d12dx,d2dx
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dX = X(2)-X(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > EQ_SYMM .and. dabs(X(1)) < dX) imin = -1

  SoA(1) = SYM1
! no use  
  SoA(2) = SYM1
  SoA(3) = SYM1

  call fh_ensure(ex)
  fh(-1:ex(1),-1:ex(2),-1:ex(3)) => fh_buf
  call symmetry_bd(2,ex,f,fh,SoA)

  d12dx = ONE/F12/dX

  d2dx = ONE/TWO/dX

  fx = ZEO

  do k=1,ex(3)-1
  do j=1,ex(2)-1
  do i=1,ex(1)-1
! x direction   
        if(i+2 <= imax .and. i-2 >= imin)then
!
!              f(i-2) - 8 f(i-1) + 8 f(i+1) - f(i+2)
!  fx(i) = ---------------------------------------------
!                             12 dx
      fx(i,j,k)=d12dx*(fh(i-2,j,k)-EIT*fh(i-1,j,k)+EIT*fh(i+1,j,k)-fh(i+2,j,k))

    elseif(i+1 <= imax .and. i-1 >= imin)then
!
!              - f(i-1) + f(i+1)
!  fx(i) = --------------------------------
!                     2 dx
      fx(i,j,k)=d2dx*(-fh(i-1,j,k)+fh(i+1,j,k))

! set imax and imin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fdx
!-----------------------------------------------------------------------------
!
! single derivatives dy
!
!-----------------------------------------------------------------------------
  subroutine fdy(ex,f,fy,Y,SYM2,symmetry,onoff)
  use diff_buffers
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fy
  real*8,                                intent(in ):: Y(ex(2)),SYM2

!~~~~~~ other variables

  real*8 :: dY
  real*8, pointer :: fh(:,:,:)
  real*8, dimension(3) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d12dy,d2dy
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dY = Y(2)-Y(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > EQ_SYMM .and. dabs(Y(1)) < dY) jmin = -1

  SoA(1) = SYM2
  SoA(2) = SYM2
  SoA(3) = SYM2

  call fh_ensure(ex)
  fh(-1:ex(1),-1:ex(2),-1:ex(3)) => fh_buf
  call symmetry_bd(2,ex,f,fh,SoA)

  d12dy = ONE/F12/dY

  d2dy = ONE/TWO/dY

  fy = ZEO

  do k=1,ex(3)-1
  do j=1,ex(2)-1
  do i=1,ex(1)-1
! y direction   
        if(j+2 <= jmax .and. j-2 >= jmin)then

      fy(i,j,k)=d12dy*(fh(i,j-2,k)-EIT*fh(i,j-1,k)+EIT*fh(i,j+1,k)-fh(i,j+2,k))

    elseif(j+1 <= jmax .and. j-1 >= jmin)then

     fy(i,j,k)=d2dy*(-fh(i,j-1,k)+fh(i,j+1,k))

! set jmax and jmin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fdy
!-----------------------------------------------------------------------------
!
! single derivatives dz
!
!-----------------------------------------------------------------------------
  subroutine fdz(ex,f,fz,Z,SYM3,symmetry,onoff)
  use diff_buffers
  implicit none

  integer,                               intent(in ):: ex(1:3),symmetry,onoff
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(in ):: f
  real*8,  dimension(ex(1),ex(2),ex(3)), intent(out):: fz
  real*8,                                intent(in ):: Z(ex(3)),SYM3

!~~~~~~ other variables
  
  real*8 :: dZ
  real*8, pointer :: fh(:,:,:)
  real*8, dimension(3) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8 :: d12dz,d2dz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F60=6.d1
  real*8,  parameter :: TWO=2.d0,EIT=8.d0
  real*8,  parameter ::  F9=9.d0,F45=4.5d1,F12=1.2d1

  dZ = Z(2)-Z(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -1

  SoA(1) = SYM3
  SoA(2) = SYM3
  SoA(3) = SYM3

  call fh_ensure(ex)
  fh(-1:ex(1),-1:ex(2),-1:ex(3)) => fh_buf
  call symmetry_bd(2,ex,f,fh,SoA)

  d12dz = ONE/F12/dZ

  d2dz = ONE/TWO/dZ

  fz = ZEO

  do k=1,ex(3)-1
  do j=1,ex(2)-1
  do i=1,ex(1)-1
! z direction   
        if(k+2 <= kmax .and. k-2 >= kmin)then

      fz(i,j,k)=d12dz*(fh(i,j,k-2)-EIT*fh(i,j,k-1)+EIT*fh(i,j,k+1)-fh(i,j,k+2))

    elseif(k+1 <= kmax .and. k-1 >= kmin)then

      fz(i,j,k)=d2dz*(-fh(i,j,k-1)+fh(i,j,k+1))

! set kmax and kmin 0
    endif

  enddo
  enddo
  enddo

  return

  end subroutine fdz
!-----------------------------------------------------------------------------------------------------------------
!
! General second derivatives of 4_th oder accurate
!
!               - f(i-2) + 16 f(i-1) - 30 f(i) + 16 f(i+1) - f(i+2)
!  fxx(i) = ----------------------------------------------------------
!                                  12 dx^2 
!
!             -   ( - f(i+2,j+2) + 8 f(i+1,j+2) - 8 f(i-1,j+2) + f(i-2,j+2) )
!             + 8 ( - f(i+2,j+1) + 8 f(i+1,j+1) - 8 f(i-1,j+1) + f(i-2,j+1) )
!             - 8 ( - f(i+2,j-1) + 8 f(i+1,j-1) - 8 f(i-1,j-1) + f(i-2,j-1) )
!             +   ( - f(i+2,j-2) + 8 f(i+1,j-2) - 8 f(i-1,j-2) + f(i-2,j-2) )
!  fxy(i,j) = ----------------------------------------------------------------
!                                  144 dx dy
!
!-----------------------------------------------------------------------------------------------------------------
  subroutine fdderivs(ex,f,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z, &
                      SYM1,SYM2,SYM3,symmetry,onoff)
  use diff_buffers
#ifdef AMSS_RHS_HALO
  use rhs_halo
#endif
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry,onoff
  real*8, dimension(ex(1),ex(2),ex(3)),target,intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxx,fxy,fxz,fyy,fyz,fzz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8, pointer :: fh(:,:,:)
  real*8, dimension(3) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k,klo,khi
  real*8  :: Sdxdx,Sdydy,Sdzdz,Fdxdx,Fdydy,Fdzdz
  real*8  :: Sdxdy,Sdxdz,Sdydz,Fdxdy,Fdxdz,Fdydz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -1
  if(Symmetry > EQ_SYMM .and. dabs(X(1)) < dX) imin = -1
  if(Symmetry > EQ_SYMM .and. dabs(Y(1)) < dY) jmin = -1

  SoA(1) = SYM1
  SoA(2) = SYM2
  SoA(3) = SYM3

#ifdef AMSS_RHS_HALO
  call halo_bind(f, ex, SoA, symmetry, 2, kmin)
#else
  call fh_ensure(ex)
  fh(-1:ex(1),-1:ex(2),-1:ex(3)) => fh_buf
  call symmetry_bd(2,ex,f,fh,SoA)
#endif

  Sdxdx =  ONE /( dX * dX )
  Sdydy =  ONE /( dY * dY )
  Sdzdz =  ONE /( dZ * dZ )

  Fdxdx = F1o12 /( dX * dX )
  Fdydy = F1o12 /( dY * dY )
  Fdzdz = F1o12 /( dZ * dZ )

  Sdxdy = F1o4 /( dX * dY )
  Sdxdz = F1o4 /( dX * dZ )
  Sdydz = F1o4 /( dY * dZ )

  Fdxdy = F1o144 /( dX * dY )
  Fdxdz = F1o144 /( dX * dZ )
  Fdydz = F1o144 /( dY * dZ )

  fxx = ZEO
  fyy = ZEO
  fzz = ZEO
  fxy = ZEO
  fxz = ZEO
  fyz = ZEO

#ifdef AMSS_RHS_HALO
  if (legacy_active) then
    klo = 1
    khi = ex(3)-1
#define FVAL(i,j,k) halo_get(i,j,k)
#include "fdderivs_loop.fh"
#undef FVAL
  else if (halo_active) then
    klo = 1
    khi = min(2, ex(3)-1)
#define FVAL(i,j,k) halo_get(i,j,k)
#include "fdderivs_loop.fh"
#undef FVAL
    klo = 3
    khi = ex(3)-1
#define FVAL(i,j,k) f(i,j,k)
#include "fdderivs_loop.fh"
#undef FVAL
  else
    klo = 1
    khi = ex(3)-1
#define FVAL(i,j,k) f(i,j,k)
#include "fdderivs_loop.fh"
#undef FVAL
  end if
#else
  klo = 1
  khi = ex(3)-1
#define FVAL(i,j,k) fh(i,j,k)
#include "fdderivs_loop.fh"
#undef FVAL
#endif

  return

  end subroutine fdderivs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! only for compute_ricci.f90 usage
!-----------------------------------------------------------------------------
  subroutine fddxx(ex,f,fxx,X,Y,Z,SYM1,SYM2,SYM3,symmetry)
  use diff_buffers
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxx
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8, pointer :: fh(:,:,:)
  real*8, dimension(3) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdx,Fdxdx
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -1
  if(Symmetry > EQ_SYMM .and. dabs(X(1)) < dX) imin = -1
  if(Symmetry > EQ_SYMM .and. dabs(Y(1)) < dY) jmin = -1

  SoA(1) = SYM1
  SoA(2) = SYM2
  SoA(3) = SYM3

  call fh_ensure(ex)
  fh(-1:ex(1),-1:ex(2),-1:ex(3)) => fh_buf
  call symmetry_bd(2,ex,f,fh,SoA)

  Sdxdx =  ONE /( dX * dX )

  Fdxdx = F1o12 /( dX * dX )

  fxx = ZEO

  do k=1,ex(3)-1
  do j=1,ex(2)-1
  do i=1,ex(1)-1
!~~~~~~ fxx
        if(i+2 <= imax .and. i-2 >= imin)then
   fxx(i,j,k) = Fdxdx*(-fh(i-2,j,k)+F16*fh(i-1,j,k)-F30*fh(i,j,k) &
                       -fh(i+2,j,k)+F16*fh(i+1,j,k)              )
   elseif(i+1 <= imax .and. i-1 >= imin)then
   fxx(i,j,k) = Sdxdx*(fh(i-1,j,k)-TWO*fh(i,j,k) &
                      +fh(i+1,j,k)              )
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddxx

  subroutine fddyy(ex,f,fyy,X,Y,Z,SYM1,SYM2,SYM3,symmetry)
  use diff_buffers
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fyy
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8, pointer :: fh(:,:,:)
  real*8, dimension(3) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdydy,Fdydy
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -1
  if(Symmetry > EQ_SYMM .and. dabs(X(1)) < dX) imin = -1
  if(Symmetry > EQ_SYMM .and. dabs(Y(1)) < dY) jmin = -1

  SoA(1) = SYM1
  SoA(2) = SYM2
  SoA(3) = SYM3

  call fh_ensure(ex)
  fh(-1:ex(1),-1:ex(2),-1:ex(3)) => fh_buf
  call symmetry_bd(2,ex,f,fh,SoA)

  Sdydy =  ONE /( dY * dY )

  Fdydy = F1o12 /( dY * dY )

  fyy = ZEO

  do k=1,ex(3)-1
  do j=1,ex(2)-1
  do i=1,ex(1)-1
!~~~~~~ fyy
        if(j+2 <= jmax .and. j-2 >= jmin)then

   fyy(i,j,k) = Fdydy*(-fh(i,j-2,k)+F16*fh(i,j-1,k)-F30*fh(i,j,k) &
                       -fh(i,j+2,k)+F16*fh(i,j+1,k)              )
   elseif(j+1 <= jmax .and. j-1 >= jmin)then

   fyy(i,j,k) = Sdydy*(fh(i,j-1,k)-TWO*fh(i,j,k) &
                      +fh(i,j+1,k)              )
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddyy

  subroutine fddzz(ex,f,fzz,X,Y,Z,SYM1,SYM2,SYM3,symmetry)
  use diff_buffers
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fzz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8, pointer :: fh(:,:,:)
  real*8, dimension(3) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdzdz,Fdzdz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -1
  if(Symmetry > EQ_SYMM .and. dabs(X(1)) < dX) imin = -1
  if(Symmetry > EQ_SYMM .and. dabs(Y(1)) < dY) jmin = -1

  SoA(1) = SYM1
  SoA(2) = SYM2
  SoA(3) = SYM3

  call fh_ensure(ex)
  fh(-1:ex(1),-1:ex(2),-1:ex(3)) => fh_buf
  call symmetry_bd(2,ex,f,fh,SoA)

  Sdzdz =  ONE /( dZ * dZ )

  Fdzdz = F1o12 /( dZ * dZ )

  fzz = ZEO

  do k=1,ex(3)-1
  do j=1,ex(2)-1
  do i=1,ex(1)-1
!~~~~~~ fzz
        if(k+2 <= kmax .and. k-2 >= kmin)then

   fzz(i,j,k) = Fdzdz*(-fh(i,j,k-2)+F16*fh(i,j,k-1)-F30*fh(i,j,k) &
                       -fh(i,j,k+2)+F16*fh(i,j,k+1)              )
   elseif(k+1 <= kmax .and. k-1 >= kmin)then

   fzz(i,j,k) = Sdzdz*(fh(i,j,k-1)-TWO*fh(i,j,k) &
                      +fh(i,j,k+1)              )
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddzz

  subroutine fddxy(ex,f,fxy,X,Y,Z,SYM1,SYM2,SYM3,symmetry)
  use diff_buffers
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxy
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8, pointer :: fh(:,:,:)
  real*8, dimension(3) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdy,Fdxdy
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -1
  if(Symmetry > EQ_SYMM .and. dabs(X(1)) < dX) imin = -1
  if(Symmetry > EQ_SYMM .and. dabs(Y(1)) < dY) jmin = -1

  SoA(1) = SYM1
  SoA(2) = SYM2
  SoA(3) = SYM3

  call fh_ensure(ex)
  fh(-1:ex(1),-1:ex(2),-1:ex(3)) => fh_buf
  call symmetry_bd(2,ex,f,fh,SoA)

  Sdxdy = F1o4 /( dX * dY )

  Fdxdy = F1o144 /( dX * dY )

  fxy = ZEO

  do k=1,ex(3)-1
  do j=1,ex(2)-1
  do i=1,ex(1)-1
!~~~~~~ fxy
       if(i+2 <= imax .and. i-2 >= imin .and. j+2 <= jmax .and. j-2 >= jmin)then

   fxy(i,j,k) = Fdxdy*(     (fh(i-2,j-2,k)-F8*fh(i-1,j-2,k)+F8*fh(i+1,j-2,k)-fh(i+2,j-2,k))  &
                       -F8 *(fh(i-2,j-1,k)-F8*fh(i-1,j-1,k)+F8*fh(i+1,j-1,k)-fh(i+2,j-1,k))  &
                       +F8 *(fh(i-2,j+1,k)-F8*fh(i-1,j+1,k)+F8*fh(i+1,j+1,k)-fh(i+2,j+1,k))  &
                       -    (fh(i-2,j+2,k)-F8*fh(i-1,j+2,k)+F8*fh(i+1,j+2,k)-fh(i+2,j+2,k)))
   elseif(i+1 <= imax .and. i-1 >= imin .and. j+1 <= jmax .and. j-1 >= jmin)then

   fxy(i,j,k) = Sdxdy*(fh(i-1,j-1,k)-fh(i+1,j-1,k)-fh(i-1,j+1,k)+fh(i+1,j+1,k))
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddxy

  subroutine fddxz(ex,f,fxz,X,Y,Z,SYM1,SYM2,SYM3,symmetry)
  use diff_buffers
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fxz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8, pointer :: fh(:,:,:)
  real*8, dimension(3) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdxdz,Fdxdz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -1
  if(Symmetry > EQ_SYMM .and. dabs(X(1)) < dX) imin = -1
  if(Symmetry > EQ_SYMM .and. dabs(Y(1)) < dY) jmin = -1

  SoA(1) = SYM1
  SoA(2) = SYM2
  SoA(3) = SYM3

  call fh_ensure(ex)
  fh(-1:ex(1),-1:ex(2),-1:ex(3)) => fh_buf
  call symmetry_bd(2,ex,f,fh,SoA)

  Sdxdz = F1o4 /( dX * dZ )

  Fdxdz = F1o144 /( dX * dZ )

  fxz = ZEO

  do k=1,ex(3)-1
  do j=1,ex(2)-1
  do i=1,ex(1)-1
!~~~~~~ fxz
       if(i+2 <= imax .and. i-2 >= imin .and. k+2 <= kmax .and. k-2 >= kmin)then
   fxz(i,j,k) = Fdxdz*(     (fh(i-2,j,k-2)-F8*fh(i-1,j,k-2)+F8*fh(i+1,j,k-2)-fh(i+2,j,k-2))  &
                       -F8 *(fh(i-2,j,k-1)-F8*fh(i-1,j,k-1)+F8*fh(i+1,j,k-1)-fh(i+2,j,k-1))  &
                       +F8 *(fh(i-2,j,k+1)-F8*fh(i-1,j,k+1)+F8*fh(i+1,j,k+1)-fh(i+2,j,k+1))  &
                       -    (fh(i-2,j,k+2)-F8*fh(i-1,j,k+2)+F8*fh(i+1,j,k+2)-fh(i+2,j,k+2)))
   elseif(i+1 <= imax .and. i-1 >= imin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fxz(i,j,k) = Sdxdz*(fh(i-1,j,k-1)-fh(i+1,j,k-1)-fh(i-1,j,k+1)+fh(i+1,j,k+1))
   endif

   enddo
   enddo
   enddo

  return

  end subroutine fddxz

  subroutine fddyz(ex,f,fyz,X,Y,Z,SYM1,SYM2,SYM3,symmetry)
  use diff_buffers
  implicit none

  integer,                             intent(in ):: ex(1:3),symmetry
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ):: f
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out):: fyz
  real*8,                              intent(in ):: X(ex(1)),Y(ex(2)),Z(ex(3)),SYM1,SYM2,SYM3

!~~~~~~ other variables
  real*8 :: dX,dY,dZ
  real*8, pointer :: fh(:,:,:)
  real*8, dimension(3) :: SoA
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k
  real*8  :: Sdydz,Fdydz
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2
  real*8, parameter :: ZEO=0.d0, ONE=1.d0, TWO=2.d0, F1o4=2.5d-1, F9=9.d0,  F45=4.5d1
  real*8, parameter :: F8=8.d0, F16=1.6d1, F30=3.d1, F27=2.7d1, F270=2.7d2, F490=4.9d2
  real*8, parameter :: F1o6=ONE/6.d0, F1o12=ONE/1.2d1, F1o144=ONE/1.44d2
  real*8, parameter :: F1o180=ONE/1.8d2,F1o3600=ONE/3.6d3

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -1
  if(Symmetry > EQ_SYMM .and. dabs(X(1)) < dX) imin = -1
  if(Symmetry > EQ_SYMM .and. dabs(Y(1)) < dY) jmin = -1

  SoA(1) = SYM1
  SoA(2) = SYM2
  SoA(3) = SYM3

  call fh_ensure(ex)
  fh(-1:ex(1),-1:ex(2),-1:ex(3)) => fh_buf
  call symmetry_bd(2,ex,f,fh,SoA)

  Sdydz = F1o4 /( dY * dZ )

  Fdydz = F1o144 /( dY * dZ )

  fyz = ZEO

  do k=1,ex(3)-1
  do j=1,ex(2)-1
  do i=1,ex(1)-1
!~~~~~~ fyz
       if(j+2 <= jmax .and. j-2 >= jmin .and. k+2 <= kmax .and. k-2 >= kmin)then
   fyz(i,j,k) = Fdydz*(     (fh(i,j-2,k-2)-F8*fh(i,j-1,k-2)+F8*fh(i,j+1,k-2)-fh(i,j+2,k-2))  &
                       -F8 *(fh(i,j-2,k-1)-F8*fh(i,j-1,k-1)+F8*fh(i,j+1,k-1)-fh(i,j+2,k-1))  &
                       +F8 *(fh(i,j-2,k+1)-F8*fh(i,j-1,k+1)+F8*fh(i,j+1,k+1)-fh(i,j+2,k+1))  &
                       -    (fh(i,j-2,k+2)-F8*fh(i,j-1,k+2)+F8*fh(i,j+1,k+2)-fh(i,j+2,k+2)))
   elseif(j+1 <= jmax .and. j-1 >= jmin .and. k+1 <= kmax .and. k-1 >= kmin)then
   fyz(i,j,k) = Sdydz*(fh(i,j-1,k-1)-fh(i,j+1,k-1)-fh(i,j-1,k+1)+fh(i,j+1,k+1))
   endif 

   enddo
   enddo
   enddo

  return

  end subroutine fddyz
