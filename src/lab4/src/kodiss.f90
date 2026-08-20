

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

!---------------------------------------------------------------------------------------------
!usual type Kreiss-Oliger type numerical dissipation
!We support cell center only
! Note the notation D_+ and D_- [P240 of B. Gustafsson, H.-O. Kreiss, and J. Oliger, Time
! Dependent Problems and Difference Methods (Wiley, New York, 1995).]
! D_+ = (f(i+1) - f(i))/h
! D_- = (f(i) - f(i-1))/h
! then we have D_+D_- = D_-D_+
!              D_+^3D_-^3 = (D_+D_-)^3 =
!    f(i-3) - 6 f(i-2) + 15 f(i-1) - 20 f(i) + 15 f(i+1) - 6 f(i+2) + f(i+3)
! -----------------------------------------------------------------------------
!                                    dx^6
! this is for 4th order accurate finite difference scheme
!---------------------------------------------------------------------------------------------
subroutine kodis(ex,X,Y,Z,f,f_rhs,SoA,Symmetry,eps)

#ifdef AMSS_RHS_HALO
use rhs_halo
#endif

implicit none
! argument variables
integer,intent(in) :: Symmetry
integer,dimension(3),intent(in)::ex
real*8, dimension(1:3), intent(in) :: SoA
double precision,intent(in),dimension(ex(1))::X
double precision,intent(in),dimension(ex(2))::Y
double precision,intent(in),dimension(ex(3))::Z
double precision,target,intent(in),dimension(ex(1),ex(2),ex(3))::f
double precision,intent(inout),dimension(ex(1),ex(2),ex(3))::f_rhs
real*8,intent(in) :: eps
! local variables
real*8,dimension(-2:ex(1),-2:ex(2),-2:ex(3))   :: fh
integer :: imin,jmin,kmin,imax,jmax,kmax,klo,khi
integer :: i,j,k
real*8  :: dX,dY,dZ
real*8, parameter :: ONE=1.d0,SIX=6.d0,FIT=1.5d1,TWT=2.d1
real*8,parameter::cof=6.4d1   ! 2^6
integer, parameter :: NO_SYMM=0, OCTANT=2

!rhs_i = rhs_i + eps/dx/cof*(f_i-3 - 6*f_i-2 + 15*f_i-1 - 20*f_i + 15*f_i+1 - 6*f_i+2 + f_i+3)

  dX = X(2)-X(1)
  dY = Y(2)-Y(1)
  dZ = Z(2)-Z(1)
  
  imax = ex(1)
  jmax = ex(2)
  kmax = ex(3)

  imin = 1
  jmin = 1
  kmin = 1

  if(Symmetry > NO_SYMM .and. dabs(Z(1)) < dZ) kmin = -2
  if(Symmetry == OCTANT .and. dabs(X(1)) < dX) imin = -2
  if(Symmetry == OCTANT .and. dabs(Y(1)) < dY) jmin = -2

#ifdef AMSS_RHS_HALO
  call halo_bind(f, ex, SoA, Symmetry, 3, kmin)
#else
  call symmetry_bd(3,ex,f,fh,SoA)
#endif

#ifdef AMSS_RHS_HALO
  if (legacy_active) then
    klo = 1
    khi = ex(3)
#define FVAL(i,j,k) halo_get(i,j,k)
#include "kodis_loop.fh"
#undef FVAL
  else if (halo_active) then
    klo = 1
    khi = min(3, ex(3))
#define FVAL(i,j,k) halo_get(i,j,k)
#include "kodis_loop.fh"
#undef FVAL
    klo = 4
    khi = ex(3)
#define FVAL(i,j,k) f(i,j,k)
#include "kodis_loop.fh"
#undef FVAL
  else
    klo = 1
    khi = ex(3)
#define FVAL(i,j,k) f(i,j,k)
#include "kodis_loop.fh"
#undef FVAL
  end if
#else
  klo = 1
  khi = ex(3)
#define FVAL(i,j,k) fh(i,j,k)
#include "kodis_loop.fh"
#undef FVAL
#endif

  return

  end subroutine kodis
