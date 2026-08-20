
#include "macrodef.fh"

! ---------------------------------------------------------------------------
! Batched kodis / lopsided stencils (optimization candidate, #18).
!
! The 24 lopsided + 24 kodis calls at the end of compute_rhs each perform a
! full-block scan with their own halo/mirror setup and per-point stencil-order
! guard.  Across the 24 calls the shift arrays (betax/betay/betaz, used by
! lopsided) are IDENTICAL and the guard depends only on (i,j,k), not on the
! field.  Batching B fields into one sweep amortizes: 24/B sweeps instead of
! 24, one guard evaluation per point, one shift read per point (lopsided), and
! one mirror setup per chunk.
!
! Correctness: the per-point arithmetic is bit-identical to kodis_loop.fh /
! lopsided_loop.fh (the active "bam" #else branch).  Equatorial symmetry
! (Symmetry==1) needs only z-reflection, applied inline by zref() matching
! symmetry_bd's funcc(:,:,-i) = funcc(:,:,i+1)*SoA(3).  Octant (Symmetry==2)
! also needs x/y reflection and is NOT handled here -- the caller falls back
! to the original per-variable routines.
!
! Gating: compile with -DAMSS_BATCH_STENCIL (CMake AMSS_BATCH_STENCIL=ON) and
! set env AMSS_BATCH_STENCIL_B to 2, 4, or 8 at run time.  With the flag unset
! or B unset, bssn_rhs uses the original 24+24 calls verbatim (#13 baseline).
! ---------------------------------------------------------------------------

module batch_stencils
  implicit none
  save

  public :: field_ptr, bs_enabled, bs_size, kodis_batch, lopsided_batch

  type :: field_ptr
    real*8, pointer :: p(:,:,:) => null()
  end type field_ptr

  logical :: bs_initialized = .false.
  integer :: bs_B = 0            ! batch size (2/4/8); 0 = disabled

contains

  subroutine bs_init()
    integer :: n, stat
    character(len=16) :: val
    if (bs_initialized) return
    call get_environment_variable("AMSS_BATCH_STENCIL_B", val, status=stat)
    read(val, *, iostat=stat) n
    if (stat == 0 .and. (n == 2 .or. n == 4 .or. n == 8)) bs_B = n
    bs_initialized = .true.
  end subroutine bs_init

  logical function bs_enabled()
    call bs_init()
    bs_enabled = (bs_B > 0)
  end function bs_enabled

  integer function bs_size()
    call bs_init()
    bs_size = bs_B
  end function bs_size

  ! Reflected value of field f at z-index kk (kk<1 -> s3*f(i,j,1-kk)).
  ! Matches symmetry_bd(:,:,-i) =(:,:,i+1)*SoA(3): for kk<1, 1-kk>=1.
  pure function zref(f, i, j, kk, s3) result(v)
    real*8, intent(in) :: f(:,:,:)
    integer, intent(in) :: i, j, kk
    real*8, intent(in) :: s3
    real*8 :: v
    if (kk < 1) then
      v = s3 * f(i, j, 1-kk)
    else
      v = f(i, j, kk)
    end if
  end function zref

  ! =====================================================================
  ! kodis_batch: Kreiss-Oliger 6th-order dissipation for nv fields.
  !   fp(1..nv)  : input fields f
  !   frp(1..nv) : RHS fields f_rhs (accumulated in place)
  !   s3(1..nv)  : z-parity sign (SoA(3)) of each field
  ! Arithmetic matches kodis_loop.fh (bam #else):
  !   f_rhs += eps/cof * ( (7pt_x)/dX + (7pt_y)/dY + (7pt_z)/dZ )
  !   guard: i-3>=imin .and. i+3<=imax .and. j-3>=jmin .and. j+3<=jmax
  !          .and. k-3>=kmin .and. k+3<=kmax
  ! =====================================================================
  subroutine kodis_batch(ex, X, Y, Z, nv, fp, frp, s3, Symmetry, eps)
    integer, intent(in) :: ex(3), nv, Symmetry
    real*8, intent(in) :: X(ex(1)), Y(ex(2)), Z(ex(3)), eps, s3(nv)
    type(field_ptr), intent(in) :: fp(nv), frp(nv)
    integer, parameter :: NO_SYMM = 0
    real*8, parameter :: SIX = 6.d0, FIT = 1.5d1, TWT = 2.d1, cof = 6.4d1
    integer :: imin, jmin, kmin, imax, jmax, kmax, i, j, k, b, cs, ce, Bsz, bb
    real*8 :: dX, dY, dZ, e, c, dxs, dys, dzs, sb
    logical :: reflect

    imax = ex(1); jmax = ex(2); kmax = ex(3)
    imin = 1; jmin = 1; kmin = 1
    dX = X(2)-X(1); dY = Y(2)-Y(1); dZ = Z(2)-Z(1)
    if (Symmetry > NO_SYMM .and. abs(Z(1)) < dZ) kmin = -2
    reflect = (kmin < 1)
    e = eps / cof

    Bsz = bs_B
    if (Bsz < 1) Bsz = nv
    if (Bsz > nv) Bsz = nv

#ifdef AMSS_BATCH_STENCIL_SINGLE_REGION
    ! The legacy layout starts one OpenMP team for every field batch.  With
    ! B=2 this is twelve fork/join regions per call.  Keep the arithmetic and
    ! field order unchanged, but distribute (k, field-batch) with one team.
    !$omp parallel do collapse(2) schedule(static) &
    !$omp& private(i,j,k,b,bb,c,dxs,dys,dzs,sb,cs,ce)
    do k = 1, kmax
      do cs = 1, nv, Bsz
        ce = min(cs + Bsz - 1, nv)
        do j = 1, ex(2)
          do i = 1, ex(1)
            if (i-3 >= imin .and. i+3 <= imax .and. &
                j-3 >= jmin .and. j+3 <= jmax .and. &
                k-3 >= kmin .and. k+3 <= kmax) then
              do b = cs, ce
                bb = b
                sb = s3(b)
                c = fp(b)%p(i,j,k)
                dxs = (fp(b)%p(i-3,j,k)+fp(b)%p(i+3,j,k)) &
                    - SIX*(fp(b)%p(i-2,j,k)+fp(b)%p(i+2,j,k)) &
                    + FIT*(fp(b)%p(i-1,j,k)+fp(b)%p(i+1,j,k)) &
                    - TWT*c
                dys = (fp(b)%p(i,j-3,k)+fp(b)%p(i,j+3,k)) &
                    - SIX*(fp(b)%p(i,j-2,k)+fp(b)%p(i,j+2,k)) &
                    + FIT*(fp(b)%p(i,j-1,k)+fp(b)%p(i,j+1,k)) &
                    - TWT*c
                if (reflect) then
                  dzs = (zref(fp(b)%p,i,j,k-3,sb)+fp(b)%p(i,j,k+3)) &
                      - SIX*(zref(fp(b)%p,i,j,k-2,sb)+fp(b)%p(i,j,k+2)) &
                      + FIT*(zref(fp(b)%p,i,j,k-1,sb)+fp(b)%p(i,j,k+1)) &
                      - TWT*c
                else
                  dzs = (fp(b)%p(i,j,k-3)+fp(b)%p(i,j,k+3)) &
                      - SIX*(fp(b)%p(i,j,k-2)+fp(b)%p(i,j,k+2)) &
                      + FIT*(fp(b)%p(i,j,k-1)+fp(b)%p(i,j,k+1)) &
                      - TWT*c
                end if
                frp(b)%p(i,j,k) = frp(b)%p(i,j,k) + e*(dxs/dX + dys/dY + dzs/dZ)
              end do
            end if
          end do
        end do
      end do
    end do
    !$omp end parallel do
#else
    do cs = 1, nv, Bsz
      ce = min(cs + Bsz - 1, nv)
      !$omp parallel do private(i,j,k,b,bb,c,dxs,dys,dzs,sb)
      do k = 1, kmax
        do j = 1, ex(2)
          do i = 1, ex(1)
            if (i-3 >= imin .and. i+3 <= imax .and. &
                j-3 >= jmin .and. j+3 <= jmax .and. &
                k-3 >= kmin .and. k+3 <= kmax) then
              do b = cs, ce
                bb = b
                sb = s3(b)
                c = fp(b)%p(i,j,k)
                dxs = (fp(b)%p(i-3,j,k)+fp(b)%p(i+3,j,k)) &
                    - SIX*(fp(b)%p(i-2,j,k)+fp(b)%p(i+2,j,k)) &
                    + FIT*(fp(b)%p(i-1,j,k)+fp(b)%p(i+1,j,k)) &
                    - TWT*c
                dys = (fp(b)%p(i,j-3,k)+fp(b)%p(i,j+3,k)) &
                    - SIX*(fp(b)%p(i,j-2,k)+fp(b)%p(i,j+2,k)) &
                    + FIT*(fp(b)%p(i,j-1,k)+fp(b)%p(i,j+1,k)) &
                    - TWT*c
                if (reflect) then
                  dzs = (zref(fp(b)%p,i,j,k-3,sb)+fp(b)%p(i,j,k+3)) &
                      - SIX*(zref(fp(b)%p,i,j,k-2,sb)+fp(b)%p(i,j,k+2)) &
                      + FIT*(zref(fp(b)%p,i,j,k-1,sb)+fp(b)%p(i,j,k+1)) &
                      - TWT*c
                else
                  dzs = (fp(b)%p(i,j,k-3)+fp(b)%p(i,j,k+3)) &
                      - SIX*(fp(b)%p(i,j,k-2)+fp(b)%p(i,j,k+2)) &
                      + FIT*(fp(b)%p(i,j,k-1)+fp(b)%p(i,j,k+1)) &
                      - TWT*c
                end if
                frp(b)%p(i,j,k) = frp(b)%p(i,j,k) + e*(dxs/dX + dys/dY + dzs/dZ)
              end do
            end if
          end do
        end do
      end do
      !$omp end parallel do
    end do
#endif
  end subroutine kodis_batch

  ! =====================================================================
  ! lopsided_batch: upwind advection for nv fields with a SHARED shift.
  !   Sfx/Sfy/Sfz = betax/betay/betaz (identical for all 24 calls)
  ! Per-axis branch (bam #else code) depends only on Sf sign and (i/j/k) vs
  ! imin/imax, so it is selected once per point and applied to all B fields.
  ! Loop bounds: i=1..ex(1)-1, j=1..ex(2)-1, k=1..ex(3)-1 (matches lopsided).
  ! =====================================================================
  subroutine lopsided_batch(ex, X, Y, Z, nv, fp, frp, s3, &
                            Sfx, Sfy, Sfz, Symmetry, kfp, eps, fuse_tail)
    integer, intent(in) :: ex(3), nv, Symmetry
    real*8, intent(in) :: X(ex(1)), Y(ex(2)), Z(ex(3)), s3(nv)
    real*8, intent(in) :: Sfx(ex(1),ex(2),ex(3))
    real*8, intent(in) :: Sfy(ex(1),ex(2),ex(3))
    real*8, intent(in) :: Sfz(ex(1),ex(2),ex(3))
    type(field_ptr), intent(in) :: fp(nv), frp(nv)
    type(field_ptr), intent(in), optional :: kfp(nv)
    real*8, intent(in), optional :: eps
    logical, intent(in), optional :: fuse_tail
    integer, parameter :: NO_SYMM = 0
    real*8, parameter :: ZEO = 0.d0, F3 = 3.d0, F10 = 1.d1, F18 = 1.8d1, F6 = 6.d0
    real*8, parameter :: F12 = 1.2d1, EIT = 8.d0
    integer :: imin, jmin, kmin, imax, jmax, kmax, i, j, k, b, cs, ce, Bsz
    integer :: xb, yb, zb
    real*8 :: dX, dY, dZ, sx, sy, sz, d12x, d12y, d12z, r, sb
    real*8 :: e, c, dxs, dys, dzs
    logical :: do_fuse_tail
    logical :: reflect

    imax = ex(1); jmax = ex(2); kmax = ex(3)
    imin = 1; jmin = 1; kmin = 1
    dX = X(2)-X(1); dY = Y(2)-Y(1); dZ = Z(2)-Z(1)
    if (Symmetry > NO_SYMM .and. abs(Z(1)) < dZ) kmin = -2
    ! Symmetry>1 (octant) needs x/y reflection -- caller excludes it.
    reflect = (kmin < 1)
    d12x = 1.d0/F12/dX; d12y = 1.d0/F12/dY; d12z = 1.d0/F12/dZ
    do_fuse_tail = .false.
    if (present(fuse_tail)) do_fuse_tail = fuse_tail
    e = 0.d0
    if (present(eps)) e = eps/6.4d1

    Bsz = bs_B
    if (Bsz < 1) Bsz = nv
    if (Bsz > nv) Bsz = nv

#ifdef AMSS_BATCH_STENCIL_SINGLE_REGION
    ! Keep one team alive across all B-sized field groups.  Each thread walks
    ! the same serial cs loop and shares the point loop with an omp do; the
    ! barrier at the end of every omp do preserves the original group order.
    !$omp parallel private(i,j,k,b,xb,yb,zb,sx,sy,sz,r,sb,cs,ce)
    do cs = 1, nv, Bsz
      ce = min(cs + Bsz - 1, nv)
      !$omp do schedule(static)
#else
    do cs = 1, nv, Bsz
      ce = min(cs + Bsz - 1, nv)
      !$omp parallel do private(i,j,k,b,xb,yb,zb,sx,sy,sz,r,sb)
#endif
      do k = 1, ex(3)-1
        do j = 1, ex(2)-1
          do i = 1, ex(1)-1
            sx = Sfx(i,j,k); sy = Sfy(i,j,k); sz = Sfz(i,j,k)
            ! ---- x direction branch (shared across fields) ----
            xb = 0
            if (sx > ZEO) then
              if (i+3 <= imax) then
                xb = 1
              elseif (i+2 <= imax) then
                xb = 2
              elseif (i+1 <= imax) then
                xb = 3
              end if
            elseif (sx < ZEO) then
              if (i-3 >= imin) then
                xb = 4
              elseif (i-2 >= imin) then
                xb = 5
              elseif (i-1 >= imin) then
                xb = 6
              end if
            end if
            ! ---- y direction branch (shared) ----
            yb = 0
            if (sy > ZEO) then
              if (j+3 <= jmax) then
                yb = 1
              elseif (j+2 <= jmax) then
                yb = 2
              elseif (j+1 <= jmax) then
                yb = 3
              end if
            elseif (sy < ZEO) then
              if (j-3 >= jmin) then
                yb = 4
              elseif (j-2 >= jmin) then
                yb = 5
              elseif (j-1 >= jmin) then
                yb = 6
              end if
            end if
            ! ---- z direction branch (shared) ----
            zb = 0
            if (sz > ZEO) then
              if (k+3 <= kmax) then
                zb = 1
              elseif (k+2 <= kmax) then
                zb = 2
              elseif (k+1 <= kmax) then
                zb = 3
              end if
            elseif (sz < ZEO) then
              if (k-3 >= kmin) then
                zb = 4
              elseif (k-2 >= kmin) then
                zb = 5
              elseif (k-1 >= kmin) then
                zb = 6
              end if
            end if

            if (xb /= 0 .or. yb /= 0 .or. zb /= 0 .or. &
                (do_fuse_tail .and. i-3 >= imin .and. i+3 <= imax .and. &
                 j-3 >= jmin .and. j+3 <= jmax .and. &
                 k-3 >= kmin .and. k+3 <= kmax)) then
              do b = cs, ce
                sb = s3(b)
                r = frp(b)%p(i,j,k)
                ! ---- x ----  A={1,6} B={2,5} C={3,4}
                select case (xb)
                case (1, 6)
                  r = r + sx*d12x*(-F3*fp(b)%p(i-1,j,k) - F10*fp(b)%p(i,j,k) &
                       + F18*fp(b)%p(i+1,j,k) - F6*fp(b)%p(i+2,j,k) + fp(b)%p(i+3,j,k))
                case (2, 5)
                  r = r + sx*d12x*(fp(b)%p(i-2,j,k) - EIT*fp(b)%p(i-1,j,k) &
                       + EIT*fp(b)%p(i+1,j,k) - fp(b)%p(i+2,j,k))
                case (3, 4)
                  r = r - sx*d12x*(-F3*fp(b)%p(i+1,j,k) - F10*fp(b)%p(i,j,k) &
                       + F18*fp(b)%p(i-1,j,k) - F6*fp(b)%p(i-2,j,k) + fp(b)%p(i-3,j,k))
                end select
                ! ---- y ----
                select case (yb)
                case (1, 6)
                  r = r + sy*d12y*(-F3*fp(b)%p(i,j-1,k) - F10*fp(b)%p(i,j,k) &
                       + F18*fp(b)%p(i,j+1,k) - F6*fp(b)%p(i,j+2,k) + fp(b)%p(i,j+3,k))
                case (2, 5)
                  r = r + sy*d12y*(fp(b)%p(i,j-2,k) - EIT*fp(b)%p(i,j-1,k) &
                       + EIT*fp(b)%p(i,j+1,k) - fp(b)%p(i,j+2,k))
                case (3, 4)
                  r = r - sy*d12y*(-F3*fp(b)%p(i,j+1,k) - F10*fp(b)%p(i,j,k) &
                       + F18*fp(b)%p(i,j-1,k) - F6*fp(b)%p(i,j-2,k) + fp(b)%p(i,j-3,k))
                end select
                ! ---- z (reflected for kk<1) ----
                if (reflect) then
                  select case (zb)
                  case (1, 6)
                    r = r + sz*d12z*(-F3*zref(fp(b)%p,i,j,k-1,sb) - F10*zref(fp(b)%p,i,j,k,sb) &
                         + F18*zref(fp(b)%p,i,j,k+1,sb) - F6*zref(fp(b)%p,i,j,k+2,sb) &
                         + zref(fp(b)%p,i,j,k+3,sb))
                  case (2, 5)
                    r = r + sz*d12z*(zref(fp(b)%p,i,j,k-2,sb) - EIT*zref(fp(b)%p,i,j,k-1,sb) &
                         + EIT*zref(fp(b)%p,i,j,k+1,sb) - zref(fp(b)%p,i,j,k+2,sb))
                  case (3, 4)
                    r = r - sz*d12z*(-F3*zref(fp(b)%p,i,j,k+1,sb) - F10*zref(fp(b)%p,i,j,k,sb) &
                         + F18*zref(fp(b)%p,i,j,k-1,sb) - F6*zref(fp(b)%p,i,j,k-2,sb) &
                         + zref(fp(b)%p,i,j,k-3,sb))
                  end select
                else
                  select case (zb)
                  case (1, 6)
                    r = r + sz*d12z*(-F3*fp(b)%p(i,j,k-1) - F10*fp(b)%p(i,j,k) &
                         + F18*fp(b)%p(i,j,k+1) - F6*fp(b)%p(i,j,k+2) + fp(b)%p(i,j,k+3))
                  case (2, 5)
                    r = r + sz*d12z*(fp(b)%p(i,j,k-2) - EIT*fp(b)%p(i,j,k-1) &
                         + EIT*fp(b)%p(i,j,k+1) - fp(b)%p(i,j,k+2))
                  case (3, 4)
                    r = r - sz*d12z*(-F3*fp(b)%p(i,j,k+1) - F10*fp(b)%p(i,j,k) &
                         + F18*fp(b)%p(i,j,k-1) - F6*fp(b)%p(i,j,k-2) + fp(b)%p(i,j,k-3))
                  end select
                end if
                ! Optional B=2 tail fusion: preserve the complete
                ! lopsided accumulation order, then apply KO to the same
                ! RHS value before writing it once.  kfp is the raw field
                ! source (diagonal metric fields differ from lfp by +1).
                if (do_fuse_tail .and. present(kfp) .and. &
                    i-3 >= imin .and. i+3 <= imax .and. &
                    j-3 >= jmin .and. j+3 <= jmax .and. &
                    k-3 >= kmin .and. k+3 <= kmax) then
                  c = kfp(b)%p(i,j,k)
                  dxs = (kfp(b)%p(i-3,j,k)+kfp(b)%p(i+3,j,k)) &
                      - 6.d0*(kfp(b)%p(i-2,j,k)+kfp(b)%p(i+2,j,k)) &
                      + 15.d0*(kfp(b)%p(i-1,j,k)+kfp(b)%p(i+1,j,k)) &
                      - 20.d0*c
                  dys = (kfp(b)%p(i,j-3,k)+kfp(b)%p(i,j+3,k)) &
                      - 6.d0*(kfp(b)%p(i,j-2,k)+kfp(b)%p(i,j+2,k)) &
                      + 15.d0*(kfp(b)%p(i,j-1,k)+kfp(b)%p(i,j+1,k)) &
                      - 20.d0*c
                  if (reflect) then
                    dzs = (zref(kfp(b)%p,i,j,k-3,sb)+kfp(b)%p(i,j,k+3)) &
                        - 6.d0*(zref(kfp(b)%p,i,j,k-2,sb)+kfp(b)%p(i,j,k+2)) &
                        + 15.d0*(zref(kfp(b)%p,i,j,k-1,sb)+kfp(b)%p(i,j,k+1)) &
                        - 20.d0*c
                  else
                    dzs = (kfp(b)%p(i,j,k-3)+kfp(b)%p(i,j,k+3)) &
                        - 6.d0*(kfp(b)%p(i,j,k-2)+kfp(b)%p(i,j,k+2)) &
                        + 15.d0*(kfp(b)%p(i,j,k-1)+kfp(b)%p(i,j,k+1)) &
                        - 20.d0*c
                  end if
                  r = r + e*(dxs/dX + dys/dY + dzs/dZ)
                end if
                frp(b)%p(i,j,k) = r
              end do
            end if
          end do
        end do
      end do
#ifdef AMSS_BATCH_STENCIL_SINGLE_REGION
      !$omp end do
#else
      !$omp end parallel do
#endif
    end do
#ifdef AMSS_BATCH_STENCIL_SINGLE_REGION
    !$omp end parallel
#endif
  end subroutine lopsided_batch

end module batch_stencils
