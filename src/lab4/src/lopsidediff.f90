
! Compute advection terms in right hand sides of field equations

#include "macrodef.fh"

#ifdef AMSS_RHS_HALO
#define FH(i,j,k) halo_get(i,j,k)
#else
#define FH(i,j,k) fh(i,j,k)
#endif

! we need only distinguish different finite difference order
! Vertex or Cell is distinguished in routine symmetry_bd which locates in
! file "fmisc.f90"

! fourth order code

!-----------------------------------------------------------------------------
!
! Compute advection terms in right hand sides of field equations
!         v
! D f = ------[ - 3f    - 10f  + 18f    - 6f     + f     ]
!  i     12dx       i-v      i      i+v     i+2v    i+3v
!
! where
!
!        i
!      |B |
! v = -----
!        i
!       B
!
!-----------------------------------------------------------------------------

subroutine lopsided(ex,X,Y,Z,f,f_rhs,Sfx,Sfy,Sfz,Symmetry,SoA)
#ifdef AMSS_RHS_HALO
  use rhs_halo
#endif
  implicit none

!~~~~~~> Input parameters:

  integer, intent(in)  :: ex(1:3),Symmetry
  real*8,  intent(in)  :: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, target,dimension(ex(1),ex(2),ex(3)),intent(in)   :: f
  real*8,dimension(ex(1),ex(2),ex(3)),intent(in)   :: Sfx,Sfy,Sfz

  real*8,dimension(ex(1),ex(2),ex(3)),intent(inout):: f_rhs
  real*8,dimension(3),intent(in) ::SoA

!~~~~~~> local variables:
! note index -2,-1,0, so we have 3 extra points
  real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3))   :: fh
  integer :: imin,jmin,kmin,imax,jmax,kmax,i,j,k,klo,khi
  real*8 :: dX,dY,dZ
  real*8 :: d12dx,d12dy,d12dz,d2dx,d2dy,d2dz
  real*8,  parameter :: ZEO=0.d0,ONE=1.d0, F3=3.d0
  real*8,  parameter :: TWO=2.d0,F6=6.0d0,F18=1.8d1
  real*8,  parameter :: F12=1.2d1, F10=1.d1,EIT=8.d0
  integer, parameter :: NO_SYMM = 0, EQ_SYMM = 1, OCTANT = 2

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)

  d12dx = ONE/F12/dX
  d12dy = ONE/F12/dY
  d12dz = ONE/F12/dZ

  d2dx = ONE/TWO/dX
  d2dy = ONE/TWO/dY
  d2dz = ONE/TWO/dZ

  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1
  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -2
  if(Symmetry > EQ_SYMM .and. dabs(X(1)) < dX) imin = -2
  if(Symmetry > EQ_SYMM .and. dabs(Y(1)) < dY) jmin = -2

#ifdef AMSS_RHS_HALO
  call halo_bind(f, ex, SoA, Symmetry, 3, kmin)
#else
  call symmetry_bd(3,ex,f,fh,SoA)
#endif

! upper bound set ex-1 only for efficiency, 
! the loop body will set ex 0 also
#ifdef AMSS_RHS_HALO
  if (legacy_active) then
    klo = 1
    khi = ex(3)-1
#define FVAL(i,j,k) halo_get(i,j,k)
    !$omp parallel do
#include "lopsided_loop.fh"
#undef FVAL
  else if (halo_active) then
    klo = 1
    khi = min(3, ex(3)-1)
#define FVAL(i,j,k) halo_get(i,j,k)
    !$omp parallel do
#include "lopsided_loop.fh"
#undef FVAL
    klo = 4
    khi = ex(3)-1
#define FVAL(i,j,k) f(i,j,k)
    !$omp parallel do
#include "lopsided_loop.fh"
#undef FVAL
  else
    klo = 1
    khi = ex(3)-1
#define FVAL(i,j,k) f(i,j,k)
    !$omp parallel do
#include "lopsided_loop.fh"
#undef FVAL
  end if
#else
  klo = 1
  khi = ex(3)-1
#define FVAL(i,j,k) fh(i,j,k)
    !$omp parallel do
#include "lopsided_loop.fh"
#undef FVAL
#endif

  return

  end subroutine lopsided

