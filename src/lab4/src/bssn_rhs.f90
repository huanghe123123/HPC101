

#include "macrodef.fh"

! Workspace for the legacy/constraint branch of the fused RHS.  The
! POINTWISE evolution branch does not need any of these 111 full 3-D arrays;
! keeping them in a module-level pool prevents gfortran from allocating and
! freeing a large automatic array for every RHS call.
module rhs_legacy_buffers
  implicit none
  save
  integer, parameter :: RHS_LEGACY_SLOTS = 111
  real*8, allocatable, target :: legacy_ws(:,:,:,:)
  integer :: legacy_ex(3) = 0
contains
  subroutine rhs_legacy_ensure(ex)
    integer, intent(in) :: ex(3)
    ! Pointer aliases in the caller must have exactly the current block
    ! extent: whole-array constraint expressions use their declared shape.
    ! Reuse is still effective for the repeated same-shape calls within a
    ! level, while avoiding a larger previous block leaking into this call.
    if (.not. allocated(legacy_ws) .or. any(ex /= legacy_ex)) then
      if (allocated(legacy_ws)) deallocate(legacy_ws)
      allocate(legacy_ws(ex(1), ex(2), ex(3), RHS_LEGACY_SLOTS))
      legacy_ex = ex
    end if
  end subroutine rhs_legacy_ensure
end module rhs_legacy_buffers

! This checkout is the vacuum BSSN benchmark.  Keep the full argument list
! and field storage for ABI compatibility, but let the optimized CPU build
! compile all stress-energy element loads into exact zero constants.
#ifdef AMSS_VACUUM_BSSN
#define rho(i,j,k)  0.d0
#define Sx(i,j,k)   0.d0
#define Sy(i,j,k)   0.d0
#define Sz(i,j,k)   0.d0
#define Sxx(i,j,k)  0.d0
#define Sxy(i,j,k)  0.d0
#define Sxz(i,j,k)  0.d0
#define Syy(i,j,k)  0.d0
#define Syz(i,j,k)  0.d0
#define Szz(i,j,k)  0.d0
#endif

  function compute_rhs_bssn(ex, T,X, Y, Z,                                     &
               chi    ,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               ham_Res, movx_Res, movy_Res, movz_Res,                          &
                        Gmx_Res, Gmy_Res, Gmz_Res,                             &
               Symmetry,Lev,eps,co)  result(gont)
! calculate constraint violation when co=0               
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev,co
  real*8, intent(in ):: T
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: dtSfx,  dtSfy,  dtSfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out),target :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
! when out, physical second kind of connection  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzyy, Gamzyz, Gamzzz
! when out, physical Ricci tensor  
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8,intent(in) :: eps
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: ham_Res, movx_Res, movy_Res, movz_Res
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gmx_Res, Gmy_Res, Gmz_Res
!  gont = 0: success; gont = 1: something wrong
  integer::gont

!~~~~~~> Other variables:

  real*8, dimension(ex(1),ex(2),ex(3)) :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: chix,chiy,chiz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx,gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy,gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz,gyyz,gyzz,gzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Lapx,Lapy,Lapz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betaxx,betaxy,betaxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betayx,betayy,betayz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betazx,betazy,betazz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxx,Gamxy,Gamxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyx,Gamyy,Gamyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzx,Gamzy,Gamzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Kx,Ky,Kz,div_beta,S
  real*8, dimension(ex(1),ex(2),ex(3)) :: f,fxx,fxy,fxz,fyy,fyz,fzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxa,Gamya,Gamza,alpn1,chin1
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz

  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8            :: dX, dY, dZ, PI
  real*8, parameter :: ZEO = 0.d0,ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8, parameter :: EIGHT = 8.D0, HALF = 0.5D0, THR = 3.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  double precision,parameter::FF = 0.75d0,eta=2.d0
  real*8, parameter :: F1o3 = 1.D0/3.D0, F2o3 = 2.D0/3.D0,F3o2=1.5d0, F1o6 = 1.D0/6.D0
  real*8, parameter :: F16=1.6d1,F8=8.d0



#ifdef AMSS_RHS_NAN_CHECK
!!! sanity check
  dX = sum(chi)+sum(trK)+sum(dxx)+sum(gxy)+sum(gxz)+sum(dyy)+sum(gyz)+sum(dzz) &
      +sum(Axx)+sum(Axy)+sum(Axz)+sum(Ayy)+sum(Ayz)+sum(Azz)                   &
      +sum(Gamx)+sum(Gamy)+sum(Gamz)                                           &
      +sum(Lap)+sum(betax)+sum(betay)+sum(betaz)
  if(dX.ne.dX) then
     if(sum(chi).ne.sum(chi))write(*,*)"bssn.f90: find NaN in chi"
     if(sum(trK).ne.sum(trK))write(*,*)"bssn.f90: find NaN in trk"
     if(sum(dxx).ne.sum(dxx))write(*,*)"bssn.f90: find NaN in dxx"
     if(sum(gxy).ne.sum(gxy))write(*,*)"bssn.f90: find NaN in gxy"
     if(sum(gxz).ne.sum(gxz))write(*,*)"bssn.f90: find NaN in gxz"
     if(sum(dyy).ne.sum(dyy))write(*,*)"bssn.f90: find NaN in dyy"
     if(sum(gyz).ne.sum(gyz))write(*,*)"bssn.f90: find NaN in gyz"
     if(sum(dzz).ne.sum(dzz))write(*,*)"bssn.f90: find NaN in dzz"
     if(sum(Axx).ne.sum(Axx))write(*,*)"bssn.f90: find NaN in Axx"
     if(sum(Axy).ne.sum(Axy))write(*,*)"bssn.f90: find NaN in Axy"
     if(sum(Axz).ne.sum(Axz))write(*,*)"bssn.f90: find NaN in Axz"
     if(sum(Ayy).ne.sum(Ayy))write(*,*)"bssn.f90: find NaN in Ayy"
     if(sum(Ayz).ne.sum(Ayz))write(*,*)"bssn.f90: find NaN in Ayz"
     if(sum(Azz).ne.sum(Azz))write(*,*)"bssn.f90: find NaN in Azz"
     if(sum(Gamx).ne.sum(Gamx))write(*,*)"bssn.f90: find NaN in Gamx"
     if(sum(Gamy).ne.sum(Gamy))write(*,*)"bssn.f90: find NaN in Gamy"
     if(sum(Gamz).ne.sum(Gamz))write(*,*)"bssn.f90: find NaN in Gamz"
     if(sum(Lap).ne.sum(Lap))write(*,*)"bssn.f90: find NaN in Lap"
     if(sum(betax).ne.sum(betax))write(*,*)"bssn.f90: find NaN in betax"
     if(sum(betay).ne.sum(betay))write(*,*)"bssn.f90: find NaN in betay"
     if(sum(betaz).ne.sum(betaz))write(*,*)"bssn.f90: find NaN in betaz"
     gont = 1
     return
  endif
#endif

  PI = dacos(-ONE)

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)

  alpn1 = Lap + ONE
  chin1 = chi + ONE
  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE

  call fderivs(ex,betax,betaxx,betaxy,betaxz,X,Y,Z,ANTI, SYM, SYM,Symmetry,Lev)
  call fderivs(ex,betay,betayx,betayy,betayz,X,Y,Z, SYM,ANTI, SYM,Symmetry,Lev)
  call fderivs(ex,betaz,betazx,betazy,betazz,X,Y,Z, SYM, SYM,ANTI,Symmetry,Lev)
  
  div_beta = betaxx + betayy + betazz
 
  call fderivs(ex,chi,chix,chiy,chiz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

  chi_rhs = F2o3 *chin1*( alpn1 * trK - div_beta ) !rhs for chi

  call fderivs(ex,dxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,Lev)
  call fderivs(ex,dyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,Lev)
  call fderivs(ex,dzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)

  gxx_rhs = - TWO * alpn1 * Axx    -  F2o3 * gxx * div_beta          + &
              TWO *(  gxx * betaxx +   gxy * betayx +   gxz * betazx)

  gyy_rhs = - TWO * alpn1 * Ayy    -  F2o3 * gyy * div_beta          + &
              TWO *(  gxy * betaxy +   gyy * betayy +   gyz * betazy)

  gzz_rhs = - TWO * alpn1 * Azz    -  F2o3 * gzz * div_beta          + &
              TWO *(  gxz * betaxz +   gyz * betayz +   gzz * betazz)

  gxy_rhs = - TWO * alpn1 * Axy    +  F1o3 * gxy    * div_beta       + &
                      gxx * betaxy                  +   gxz * betazy + &
                                       gyy * betayx +   gyz * betazx   &
                                                    -   gxy * betazz

  gyz_rhs = - TWO * alpn1 * Ayz    +  F1o3 * gyz    * div_beta       + &
                      gxy * betaxz +   gyy * betayz                  + &
                      gxz * betaxy                  +   gzz * betazy   &
                                                    -   gyz * betaxx
 
  gxz_rhs = - TWO * alpn1 * Axz    +  F1o3 * gxz    * div_beta       + &
                      gxx * betaxz +   gxy * betayz                  + &
                                       gyz * betayx +   gzz * betazx   &
                                                    -   gxz * betayy     !rhs for gij

! invert tilted metric
  gupzz =  gxx * gyy * gzz + gxy * gyz * gxz + gxz * gxy * gyz - &
           gxz * gyy * gxz - gxy * gxy * gzz - gxx * gyz * gyz
  gupxx =   ( gyy * gzz - gyz * gyz ) / gupzz
  gupxy = - ( gxy * gzz - gyz * gxz ) / gupzz
  gupxz =   ( gxy * gyz - gyy * gxz ) / gupzz
  gupyy =   ( gxx * gzz - gxz * gxz ) / gupzz
  gupyz = - ( gxx * gyz - gxy * gxz ) / gupzz
  gupzz =   ( gxx * gyy - gxy * gxy ) / gupzz

  if(co == 0)then
! Gam^i_Res = Gam^i + gup^ij_,j
  Gmx_Res = Gamx - (gupxx*(gupxx*gxxx+gupxy*gxyx+gupxz*gxzx)&
                   +gupxy*(gupxx*gxyx+gupxy*gyyx+gupxz*gyzx)&
                   +gupxz*(gupxx*gxzx+gupxy*gyzx+gupxz*gzzx)&
                   +gupxx*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)&
                   +gupxy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)&
                   +gupxz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)&
                   +gupxx*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                   +gupxy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                   +gupxz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  Gmy_Res = Gamy - (gupxx*(gupxy*gxxx+gupyy*gxyx+gupyz*gxzx)&
                   +gupxy*(gupxy*gxyx+gupyy*gyyx+gupyz*gyzx)&
                   +gupxz*(gupxy*gxzx+gupyy*gyzx+gupyz*gzzx)&
                   +gupxy*(gupxy*gxxy+gupyy*gxyy+gupyz*gxzy)&
                   +gupyy*(gupxy*gxyy+gupyy*gyyy+gupyz*gyzy)&
                   +gupyz*(gupxy*gxzy+gupyy*gyzy+gupyz*gzzy)&
                   +gupxy*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                   +gupyy*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                   +gupyz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  Gmz_Res = Gamz - (gupxx*(gupxz*gxxx+gupyz*gxyx+gupzz*gxzx)&
                   +gupxy*(gupxz*gxyx+gupyz*gyyx+gupzz*gyzx)&
                   +gupxz*(gupxz*gxzx+gupyz*gyzx+gupzz*gzzx)&
                   +gupxy*(gupxz*gxxy+gupyz*gxyy+gupzz*gxzy)&
                   +gupyy*(gupxz*gxyy+gupyz*gyyy+gupzz*gyzy)&
                   +gupyz*(gupxz*gxzy+gupyz*gyzy+gupzz*gzzy)&
                   +gupxz*(gupxz*gxxz+gupyz*gxyz+gupzz*gxzz)&
                   +gupyz*(gupxz*gxyz+gupyz*gyyz+gupzz*gyzz)&
                   +gupzz*(gupxz*gxzz+gupyz*gyzz+gupzz*gzzz))
  endif

! second kind of connection
  Gamxxx =HALF*( gupxx*gxxx + gupxy*(TWO*gxyx - gxxy ) + gupxz*(TWO*gxzx - gxxz ))
  Gamyxx =HALF*( gupxy*gxxx + gupyy*(TWO*gxyx - gxxy ) + gupyz*(TWO*gxzx - gxxz ))
  Gamzxx =HALF*( gupxz*gxxx + gupyz*(TWO*gxyx - gxxy ) + gupzz*(TWO*gxzx - gxxz ))
 
  Gamxyy =HALF*( gupxx*(TWO*gxyy - gyyx ) + gupxy*gyyy + gupxz*(TWO*gyzy - gyyz ))
  Gamyyy =HALF*( gupxy*(TWO*gxyy - gyyx ) + gupyy*gyyy + gupyz*(TWO*gyzy - gyyz ))
  Gamzyy =HALF*( gupxz*(TWO*gxyy - gyyx ) + gupyz*gyyy + gupzz*(TWO*gyzy - gyyz ))

  Gamxzz =HALF*( gupxx*(TWO*gxzz - gzzx ) + gupxy*(TWO*gyzz - gzzy ) + gupxz*gzzz)
  Gamyzz =HALF*( gupxy*(TWO*gxzz - gzzx ) + gupyy*(TWO*gyzz - gzzy ) + gupyz*gzzz)
  Gamzzz =HALF*( gupxz*(TWO*gxzz - gzzx ) + gupyz*(TWO*gyzz - gzzy ) + gupzz*gzzz)

  Gamxxy =HALF*( gupxx*gxxy + gupxy*gyyx + gupxz*( gxzy + gyzx - gxyz ) )
  Gamyxy =HALF*( gupxy*gxxy + gupyy*gyyx + gupyz*( gxzy + gyzx - gxyz ) )
  Gamzxy =HALF*( gupxz*gxxy + gupyz*gyyx + gupzz*( gxzy + gyzx - gxyz ) )

  Gamxxz =HALF*( gupxx*gxxz + gupxy*( gxyz + gyzx - gxzy ) + gupxz*gzzx )
  Gamyxz =HALF*( gupxy*gxxz + gupyy*( gxyz + gyzx - gxzy ) + gupyz*gzzx )
  Gamzxz =HALF*( gupxz*gxxz + gupyz*( gxyz + gyzx - gxzy ) + gupzz*gzzx )

  Gamxyz =HALF*( gupxx*( gxyz + gxzy - gyzx ) + gupxy*gyyz + gupxz*gzzy )
  Gamyyz =HALF*( gupxy*( gxyz + gxzy - gyzx ) + gupyy*gyyz + gupyz*gzzy )
  Gamzyz =HALF*( gupxz*( gxyz + gxzy - gyzx ) + gupyz*gyyz + gupzz*gzzy )
! Raise indices of \tilde A_{ij} and store in R_ij

  Rxx =    gupxx * gupxx * Axx + gupxy * gupxy * Ayy + gupxz * gupxz * Azz + &
      TWO*(gupxx * gupxy * Axy + gupxx * gupxz * Axz + gupxy * gupxz * Ayz)

  Ryy =    gupxy * gupxy * Axx + gupyy * gupyy * Ayy + gupyz * gupyz * Azz + &
      TWO*(gupxy * gupyy * Axy + gupxy * gupyz * Axz + gupyy * gupyz * Ayz)

  Rzz =    gupxz * gupxz * Axx + gupyz * gupyz * Ayy + gupzz * gupzz * Azz + &
      TWO*(gupxz * gupyz * Axy + gupxz * gupzz * Axz + gupyz * gupzz * Ayz)

  Rxy =    gupxx * gupxy * Axx + gupxy * gupyy * Ayy + gupxz * gupyz * Azz + &
          (gupxx * gupyy       + gupxy * gupxy)* Axy                       + &
          (gupxx * gupyz       + gupxz * gupxy)* Axz                       + &
          (gupxy * gupyz       + gupxz * gupyy)* Ayz

  Rxz =    gupxx * gupxz * Axx + gupxy * gupyz * Ayy + gupxz * gupzz * Azz + &
          (gupxx * gupyz       + gupxy * gupxz)* Axy                       + &
          (gupxx * gupzz       + gupxz * gupxz)* Axz                       + &
          (gupxy * gupzz       + gupxz * gupyz)* Ayz

  Ryz =    gupxy * gupxz * Axx + gupyy * gupyz * Ayy + gupyz * gupzz * Azz + &
          (gupxy * gupyz       + gupyy * gupxz)* Axy                       + &
          (gupxy * gupzz       + gupyz * gupxz)* Axz                       + &
          (gupyy * gupzz       + gupyz * gupyz)* Ayz

! Right hand side for Gam^i without shift terms...
  call fderivs(ex,Lap,Lapx,Lapy,Lapz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  call fderivs(ex,trK,Kx,Ky,Kz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

   Gamx_rhs = - TWO * (   Lapx * Rxx +   Lapy * Rxy +   Lapz * Rxz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Rxx +   chiy * Rxy +   chiz * Rxz ) - &
              gupxx * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupxy * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupxz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamxxx * Rxx + Gamxyy * Ryy + Gamxzz * Rzz   + &
                TWO * ( Gamxxy * Rxy + Gamxxz * Rxz + Gamxyz * Ryz ) )

   Gamy_rhs = - TWO * (   Lapx * Rxy +   Lapy * Ryy +   Lapz * Ryz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Rxy +  chiy * Ryy +    chiz * Ryz ) - &
              gupxy * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupyy * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupyz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamyxx * Rxx + Gamyyy * Ryy + Gamyzz * Rzz   + &
                TWO * ( Gamyxy * Rxy + Gamyxz * Rxz + Gamyyz * Ryz ) )

   Gamz_rhs = - TWO * (   Lapx * Rxz +   Lapy * Ryz +   Lapz * Rzz ) + &
        TWO * alpn1 * (                                                &
        -F3o2/chin1 * (   chix * Rxz +  chiy * Ryz +    chiz * Rzz ) - &
              gupxz * (   F2o3 * Kx  +  EIGHT * PI * Sx            ) - &
              gupyz * (   F2o3 * Ky  +  EIGHT * PI * Sy            ) - &
              gupzz * (   F2o3 * Kz  +  EIGHT * PI * Sz            ) + &
                        Gamzxx * Rxx + Gamzyy * Ryy + Gamzzz * Rzz   + &
                TWO * ( Gamzxy * Rxy + Gamzxz * Rxz + Gamzyz * Ryz ) )

  call fdderivs(ex,betax,gxxx,gxyx,gxzx,gyyx,gyzx,gzzx,&
                X,Y,Z,ANTI,SYM, SYM ,Symmetry,Lev)
  call fdderivs(ex,betay,gxxy,gxyy,gxzy,gyyy,gyzy,gzzy,&
                X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
  call fdderivs(ex,betaz,gxxz,gxyz,gxzz,gyyz,gyzz,gzzz,&
                X,Y,Z,SYM ,SYM, ANTI,Symmetry,Lev)

  fxx = gxxx + gxyy + gxzz
  fxy = gxyx + gyyy + gyzz
  fxz = gxzx + gyzy + gzzz

  Gamxa =       gupxx * Gamxxx + gupyy * Gamxyy + gupzz * Gamxzz + &
          TWO*( gupxy * Gamxxy + gupxz * Gamxxz + gupyz * Gamxyz )
  Gamya =       gupxx * Gamyxx + gupyy * Gamyyy + gupzz * Gamyzz + &
          TWO*( gupxy * Gamyxy + gupxz * Gamyxz + gupyz * Gamyyz )
  Gamza =       gupxx * Gamzxx + gupyy * Gamzyy + gupzz * Gamzzz + &
          TWO*( gupxy * Gamzxy + gupxz * Gamzxz + gupyz * Gamzyz )

  call fderivs(ex,Gamx,Gamxx,Gamxy,Gamxz,X,Y,Z,ANTI,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamy,Gamyx,Gamyy,Gamyz,X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamz,Gamzx,Gamzy,Gamzz,X,Y,Z,SYM ,SYM ,ANTI,Symmetry,Lev)

  Gamx_rhs =               Gamx_rhs +  F2o3 *  Gamxa * div_beta        - &
                     Gamxa * betaxx - Gamya * betaxy - Gamza * betaxz  + &
             F1o3 * (gupxx * fxx    + gupxy * fxy    + gupxz * fxz    ) + &
                     gupxx * gxxx   + gupyy * gyyx   + gupzz * gzzx    + &
              TWO * (gupxy * gxyx   + gupxz * gxzx   + gupyz * gyzx  )

  Gamy_rhs =               Gamy_rhs +  F2o3 *  Gamya * div_beta        - &
                     Gamxa * betayx - Gamya * betayy - Gamza * betayz  + &
             F1o3 * (gupxy * fxx    + gupyy * fxy    + gupyz * fxz    ) + &
                     gupxx * gxxy   + gupyy * gyyy   + gupzz * gzzy    + &
              TWO * (gupxy * gxyy   + gupxz * gxzy   + gupyz * gyzy  )

  Gamz_rhs =               Gamz_rhs +  F2o3 *  Gamza * div_beta        - &
                     Gamxa * betazx - Gamya * betazy - Gamza * betazz  + &
             F1o3 * (gupxz * fxx    + gupyz * fxy    + gupzz * fxz    ) + &
                     gupxx * gxxz   + gupyy * gyyz   + gupzz * gzzz    + &
              TWO * (gupxy * gxyz   + gupxz * gxzz   + gupyz * gyzz  )    !rhs for Gam^i

!first kind of connection stored in gij,k
  gxxx = gxx * Gamxxx + gxy * Gamyxx + gxz * Gamzxx
  gxyx = gxx * Gamxxy + gxy * Gamyxy + gxz * Gamzxy
  gxzx = gxx * Gamxxz + gxy * Gamyxz + gxz * Gamzxz
  gyyx = gxx * Gamxyy + gxy * Gamyyy + gxz * Gamzyy
  gyzx = gxx * Gamxyz + gxy * Gamyyz + gxz * Gamzyz
  gzzx = gxx * Gamxzz + gxy * Gamyzz + gxz * Gamzzz

  gxxy = gxy * Gamxxx + gyy * Gamyxx + gyz * Gamzxx
  gxyy = gxy * Gamxxy + gyy * Gamyxy + gyz * Gamzxy
  gxzy = gxy * Gamxxz + gyy * Gamyxz + gyz * Gamzxz
  gyyy = gxy * Gamxyy + gyy * Gamyyy + gyz * Gamzyy
  gyzy = gxy * Gamxyz + gyy * Gamyyz + gyz * Gamzyz
  gzzy = gxy * Gamxzz + gyy * Gamyzz + gyz * Gamzzz

  gxxz = gxz * Gamxxx + gyz * Gamyxx + gzz * Gamzxx
  gxyz = gxz * Gamxxy + gyz * Gamyxy + gzz * Gamzxy
  gxzz = gxz * Gamxxz + gyz * Gamyxz + gzz * Gamzxz
  gyyz = gxz * Gamxyy + gyz * Gamyyy + gzz * Gamzyy
  gyzz = gxz * Gamxyz + gyz * Gamyyz + gzz * Gamzyz
  gzzz = gxz * Gamxzz + gyz * Gamyzz + gzz * Gamzzz

!compute Ricci tensor for tilted metric
   call fdderivs(ex,dxx,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
   Rxx =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,dyy,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
   Ryy =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,dzz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
   Rzz =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,gxy,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,ANTI, ANTI,SYM ,symmetry,Lev)
   Rxy =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,gxz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,ANTI ,SYM ,ANTI,symmetry,Lev)
   Rxz =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

   call fdderivs(ex,gyz,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM ,ANTI ,ANTI,symmetry,Lev)
   Ryz =   gupxx * fxx + gupyy * fyy + gupzz * fzz + &
         ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) * TWO

  Rxx =     - HALF * Rxx                                   + &
               gxx * Gamxx+ gxy * Gamyx   +    gxz * Gamzx + &
             Gamxa * gxxx +  Gamya * gxyx +  Gamza * gxzx  + &
   gupxx *(                                                  &
       TWO*(Gamxxx * gxxx + Gamyxx * gxyx + Gamzxx * gxzx) + &
            Gamxxx * gxxx + Gamyxx * gxxy + Gamzxx * gxxz )+ &
   gupxy *(                                                  &
       TWO*(Gamxxx * gxyx + Gamyxx * gyyx + Gamzxx * gyzx  + &
            Gamxxy * gxxx + Gamyxy * gxyx + Gamzxy * gxzx) + &
            Gamxxy * gxxx + Gamyxy * gxxy + Gamzxy * gxxz  + &
            Gamxxx * gxyx + Gamyxx * gxyy + Gamzxx * gxyz )+ &
   gupxz *(                                                  &
       TWO*(Gamxxx * gxzx + Gamyxx * gyzx + Gamzxx * gzzx  + &
            Gamxxz * gxxx + Gamyxz * gxyx + Gamzxz * gxzx) + &
            Gamxxz * gxxx + Gamyxz * gxxy + Gamzxz * gxxz  + &
            Gamxxx * gxzx + Gamyxx * gxzy + Gamzxx * gxzz )+ &
   gupyy *(                                                  &
       TWO*(Gamxxy * gxyx + Gamyxy * gyyx + Gamzxy * gyzx) + &
            Gamxxy * gxyx + Gamyxy * gxyy + Gamzxy * gxyz )+ &
   gupyz *(                                                  &
       TWO*(Gamxxy * gxzx + Gamyxy * gyzx + Gamzxy * gzzx  + &
            Gamxxz * gxyx + Gamyxz * gyyx + Gamzxz * gyzx) + &
            Gamxxz * gxyx + Gamyxz * gxyy + Gamzxz * gxyz  + &
            Gamxxy * gxzx + Gamyxy * gxzy + Gamzxy * gxzz )+ &
   gupzz *(                                                  &
       TWO*(Gamxxz * gxzx + Gamyxz * gyzx + Gamzxz * gzzx) + &
            Gamxxz * gxzx + Gamyxz * gxzy + Gamzxz * gxzz )

  Ryy =     - HALF * Ryy                                   + &
               gxy * Gamxy+  gyy * Gamyy  +  gyz * Gamzy   + &
             Gamxa * gxyy +  Gamya * gyyy +  Gamza * gyzy  + &
   gupxx *(                                                  &
       TWO*(Gamxxy * gxxy + Gamyxy * gxyy + Gamzxy * gxzy) + &
            Gamxxy * gxyx + Gamyxy * gxyy + Gamzxy * gxyz )+ &
   gupxy *(                                                  &
       TWO*(Gamxxy * gxyy + Gamyxy * gyyy + Gamzxy * gyzy  + &
            Gamxyy * gxxy + Gamyyy * gxyy + Gamzyy * gxzy) + &
            Gamxyy * gxyx + Gamyyy * gxyy + Gamzyy * gxyz  + &
            Gamxxy * gyyx + Gamyxy * gyyy + Gamzxy * gyyz )+ &
   gupxz *(                                                  &
       TWO*(Gamxxy * gxzy + Gamyxy * gyzy + Gamzxy * gzzy  + &
            Gamxyz * gxxy + Gamyyz * gxyy + Gamzyz * gxzy) + &
            Gamxyz * gxyx + Gamyyz * gxyy + Gamzyz * gxyz  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupyy *(                                                  &
       TWO*(Gamxyy * gxyy + Gamyyy * gyyy + Gamzyy * gyzy) + &
            Gamxyy * gyyx + Gamyyy * gyyy + Gamzyy * gyyz )+ &
   gupyz *(                                                  &
       TWO*(Gamxyy * gxzy + Gamyyy * gyzy + Gamzyy * gzzy  + &
            Gamxyz * gxyy + Gamyyz * gyyy + Gamzyz * gyzy) + &
            Gamxyz * gyyx + Gamyyz * gyyy + Gamzyz * gyyz  + &
            Gamxyy * gyzx + Gamyyy * gyzy + Gamzyy * gyzz )+ &
   gupzz *(                                                  &
       TWO*(Gamxyz * gxzy + Gamyyz * gyzy + Gamzyz * gzzy) + &
            Gamxyz * gyzx + Gamyyz * gyzy + Gamzyz * gyzz )

  Rzz =     - HALF * Rzz                                   + &
               gxz * Gamxz+ gyz * Gamyz  +    gzz * Gamzz  + &
             Gamxa * gxzz +  Gamya * gyzz +  Gamza * gzzz  + &
   gupxx *(                                                  &
       TWO*(Gamxxz * gxxz + Gamyxz * gxyz + Gamzxz * gxzz) + &
            Gamxxz * gxzx + Gamyxz * gxzy + Gamzxz * gxzz )+ &
   gupxy *(                                                  &
       TWO*(Gamxxz * gxyz + Gamyxz * gyyz + Gamzxz * gyzz  + &
            Gamxyz * gxxz + Gamyyz * gxyz + Gamzyz * gxzz) + &
            Gamxyz * gxzx + Gamyyz * gxzy + Gamzyz * gxzz  + &
            Gamxxz * gyzx + Gamyxz * gyzy + Gamzxz * gyzz )+ &
   gupxz *(                                                  &
       TWO*(Gamxxz * gxzz + Gamyxz * gyzz + Gamzxz * gzzz  + &
            Gamxzz * gxxz + Gamyzz * gxyz + Gamzzz * gxzz) + &
            Gamxzz * gxzx + Gamyzz * gxzy + Gamzzz * gxzz  + &
            Gamxxz * gzzx + Gamyxz * gzzy + Gamzxz * gzzz )+ &
   gupyy *(                                                  &
       TWO*(Gamxyz * gxyz + Gamyyz * gyyz + Gamzyz * gyzz) + &
            Gamxyz * gyzx + Gamyyz * gyzy + Gamzyz * gyzz )+ &
   gupyz *(                                                  &
       TWO*(Gamxyz * gxzz + Gamyyz * gyzz + Gamzyz * gzzz  + &
            Gamxzz * gxyz + Gamyzz * gyyz + Gamzzz * gyzz) + &
            Gamxzz * gyzx + Gamyzz * gyzy + Gamzzz * gyzz  + &
            Gamxyz * gzzx + Gamyyz * gzzy + Gamzyz * gzzz )+ &
   gupzz *(                                                  &
       TWO*(Gamxzz * gxzz + Gamyzz * gyzz + Gamzzz * gzzz) + &
            Gamxzz * gzzx + Gamyzz * gzzy + Gamzzz * gzzz )

  Rxy = HALF*(     - Rxy                                   + &
               gxx * Gamxy +    gxy * Gamyy + gxz * Gamzy  + &
               gxy * Gamxx +    gyy * Gamyx + gyz * Gamzx  + &
             Gamxa * gxyx +  Gamya * gyyx +  Gamza * gyzx  + &
             Gamxa * gxxy +  Gamya * gxyy +  Gamza * gxzy )+ &
   gupxx *(                                                  &
            Gamxxx * gxxy + Gamyxx * gxyy + Gamzxx * gxzy  + &
            Gamxxy * gxxx + Gamyxy * gxyx + Gamzxy * gxzx  + &
            Gamxxx * gxyx + Gamyxx * gxyy + Gamzxx * gxyz )+ &
   gupxy *(                                                  &
            Gamxxx * gxyy + Gamyxx * gyyy + Gamzxx * gyzy  + &
            Gamxxy * gxyx + Gamyxy * gyyx + Gamzxy * gyzx  + &
            Gamxxy * gxyx + Gamyxy * gxyy + Gamzxy * gxyz  + &
            Gamxxy * gxxy + Gamyxy * gxyy + Gamzxy * gxzy  + &
            Gamxyy * gxxx + Gamyyy * gxyx + Gamzyy * gxzx  + &
            Gamxxx * gyyx + Gamyxx * gyyy + Gamzxx * gyyz )+ &
   gupxz *(                                                  &
            Gamxxx * gxzy + Gamyxx * gyzy + Gamzxx * gzzy  + &
            Gamxxy * gxzx + Gamyxy * gyzx + Gamzxy * gzzx  + &
            Gamxxz * gxyx + Gamyxz * gxyy + Gamzxz * gxyz  + &
            Gamxxz * gxxy + Gamyxz * gxyy + Gamzxz * gxzy  + &
            Gamxyz * gxxx + Gamyyz * gxyx + Gamzyz * gxzx  + &
            Gamxxx * gyzx + Gamyxx * gyzy + Gamzxx * gyzz )+ &
   gupyy *(                                                  &
            Gamxxy * gxyy + Gamyxy * gyyy + Gamzxy * gyzy  + &
            Gamxyy * gxyx + Gamyyy * gyyx + Gamzyy * gyzx  + &
            Gamxxy * gyyx + Gamyxy * gyyy + Gamzxy * gyyz )+ &
   gupyz *(                                                  &
            Gamxxy * gxzy + Gamyxy * gyzy + Gamzxy * gzzy  + &
            Gamxyy * gxzx + Gamyyy * gyzx + Gamzyy * gzzx  + &
            Gamxxz * gyyx + Gamyxz * gyyy + Gamzxz * gyyz  + &
            Gamxxz * gxyy + Gamyxz * gyyy + Gamzxz * gyzy  + &
            Gamxyz * gxyx + Gamyyz * gyyx + Gamzyz * gyzx  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupzz *(                                                  &
            Gamxxz * gxzy + Gamyxz * gyzy + Gamzxz * gzzy  + &
            Gamxyz * gxzx + Gamyyz * gyzx + Gamzyz * gzzx  + &
            Gamxxz * gyzx + Gamyxz * gyzy + Gamzxz * gyzz )

  Rxz = HALF*(     - Rxz                                   + &
               gxx * Gamxz +  gxy * Gamyz + gxz * Gamzz    + &
               gxz * Gamxx +  gyz * Gamyx + gzz * Gamzx    + &
             Gamxa * gxzx +  Gamya * gyzx +  Gamza * gzzx  + &
             Gamxa * gxxz +  Gamya * gxyz +  Gamza * gxzz )+ &
   gupxx *(                                                  &
            Gamxxx * gxxz + Gamyxx * gxyz + Gamzxx * gxzz  + &
            Gamxxz * gxxx + Gamyxz * gxyx + Gamzxz * gxzx  + &
            Gamxxx * gxzx + Gamyxx * gxzy + Gamzxx * gxzz )+ &
   gupxy *(                                                  &
            Gamxxx * gxyz + Gamyxx * gyyz + Gamzxx * gyzz  + &
            Gamxxz * gxyx + Gamyxz * gyyx + Gamzxz * gyzx  + &
            Gamxxy * gxzx + Gamyxy * gxzy + Gamzxy * gxzz  + &
            Gamxxy * gxxz + Gamyxy * gxyz + Gamzxy * gxzz  + &
            Gamxyz * gxxx + Gamyyz * gxyx + Gamzyz * gxzx  + &
            Gamxxx * gyzx + Gamyxx * gyzy + Gamzxx * gyzz )+ &
   gupxz *(                                                  &
            Gamxxx * gxzz + Gamyxx * gyzz + Gamzxx * gzzz  + &
            Gamxxz * gxzx + Gamyxz * gyzx + Gamzxz * gzzx  + &
            Gamxxz * gxzx + Gamyxz * gxzy + Gamzxz * gxzz  + &
            Gamxxz * gxxz + Gamyxz * gxyz + Gamzxz * gxzz  + &
            Gamxzz * gxxx + Gamyzz * gxyx + Gamzzz * gxzx  + &
            Gamxxx * gzzx + Gamyxx * gzzy + Gamzxx * gzzz )+ &
   gupyy *(                                                  &
            Gamxxy * gxyz + Gamyxy * gyyz + Gamzxy * gyzz  + &
            Gamxyz * gxyx + Gamyyz * gyyx + Gamzyz * gyzx  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupyz *(                                                  &
            Gamxxy * gxzz + Gamyxy * gyzz + Gamzxy * gzzz  + &
            Gamxyz * gxzx + Gamyyz * gyzx + Gamzyz * gzzx  + &
            Gamxxz * gyzx + Gamyxz * gyzy + Gamzxz * gyzz  + &
            Gamxxz * gxyz + Gamyxz * gyyz + Gamzxz * gyzz  + &
            Gamxzz * gxyx + Gamyzz * gyyx + Gamzzz * gyzx  + &
            Gamxxy * gzzx + Gamyxy * gzzy + Gamzxy * gzzz )+ &
   gupzz *(                                                  &
            Gamxxz * gxzz + Gamyxz * gyzz + Gamzxz * gzzz  + &
            Gamxzz * gxzx + Gamyzz * gyzx + Gamzzz * gzzx  + &
            Gamxxz * gzzx + Gamyxz * gzzy + Gamzxz * gzzz )

  Ryz = HALF*(     - Ryz                                   + &
               gxy * Gamxz + gyy * Gamyz + gyz * Gamzz     + &
               gxz * Gamxy + gyz * Gamyy + gzz * Gamzy     + &
             Gamxa * gxzy +  Gamya * gyzy +  Gamza * gzzy  + &
             Gamxa * gxyz +  Gamya * gyyz +  Gamza * gyzz )+ &
   gupxx *(                                                  &
            Gamxxy * gxxz + Gamyxy * gxyz + Gamzxy * gxzz  + &
            Gamxxz * gxxy + Gamyxz * gxyy + Gamzxz * gxzy  + &
            Gamxxy * gxzx + Gamyxy * gxzy + Gamzxy * gxzz )+ &
   gupxy *(                                                  &
            Gamxxy * gxyz + Gamyxy * gyyz + Gamzxy * gyzz  + &
            Gamxxz * gxyy + Gamyxz * gyyy + Gamzxz * gyzy  + &
            Gamxyy * gxzx + Gamyyy * gxzy + Gamzyy * gxzz  + &
            Gamxyy * gxxz + Gamyyy * gxyz + Gamzyy * gxzz  + &
            Gamxyz * gxxy + Gamyyz * gxyy + Gamzyz * gxzy  + &
            Gamxxy * gyzx + Gamyxy * gyzy + Gamzxy * gyzz )+ &
   gupxz *(                                                  &
            Gamxxy * gxzz + Gamyxy * gyzz + Gamzxy * gzzz  + &
            Gamxxz * gxzy + Gamyxz * gyzy + Gamzxz * gzzy  + &
            Gamxyz * gxzx + Gamyyz * gxzy + Gamzyz * gxzz  + &
            Gamxyz * gxxz + Gamyyz * gxyz + Gamzyz * gxzz  + &
            Gamxzz * gxxy + Gamyzz * gxyy + Gamzzz * gxzy  + &
            Gamxxy * gzzx + Gamyxy * gzzy + Gamzxy * gzzz )+ &
   gupyy *(                                                  &
            Gamxyy * gxyz + Gamyyy * gyyz + Gamzyy * gyzz  + &
            Gamxyz * gxyy + Gamyyz * gyyy + Gamzyz * gyzy  + &
            Gamxyy * gyzx + Gamyyy * gyzy + Gamzyy * gyzz )+ &
   gupyz *(                                                  &
            Gamxyy * gxzz + Gamyyy * gyzz + Gamzyy * gzzz  + &
            Gamxyz * gxzy + Gamyyz * gyzy + Gamzyz * gzzy  + &
            Gamxyz * gyzx + Gamyyz * gyzy + Gamzyz * gyzz  + &
            Gamxyz * gxyz + Gamyyz * gyyz + Gamzyz * gyzz  + &
            Gamxzz * gxyy + Gamyzz * gyyy + Gamzzz * gyzy  + &
            Gamxyy * gzzx + Gamyyy * gzzy + Gamzyy * gzzz )+ &
   gupzz *(                                                  &
            Gamxyz * gxzz + Gamyyz * gyzz + Gamzyz * gzzz  + &
            Gamxzz * gxzy + Gamyzz * gyzy + Gamzzz * gzzy  + &
            Gamxyz * gzzx + Gamyyz * gzzy + Gamzyz * gzzz )
!covariant second derivative of chi respect to tilted metric
  call fdderivs(ex,chi,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)

  fxx = fxx - Gamxxx * chix - Gamyxx * chiy - Gamzxx * chiz
  fxy = fxy - Gamxxy * chix - Gamyxy * chiy - Gamzxy * chiz
  fxz = fxz - Gamxxz * chix - Gamyxz * chiy - Gamzxz * chiz
  fyy = fyy - Gamxyy * chix - Gamyyy * chiy - Gamzyy * chiz
  fyz = fyz - Gamxyz * chix - Gamyyz * chiy - Gamzyz * chiz
  fzz = fzz - Gamxzz * chix - Gamyzz * chiy - Gamzzz * chiz
! Store D^l D_l chi - 3/(2*chi) D^l chi D_l chi in f

  f =        gupxx * ( fxx - F3o2/chin1 * chix * chix ) + &
             gupyy * ( fyy - F3o2/chin1 * chiy * chiy ) + &
             gupzz * ( fzz - F3o2/chin1 * chiz * chiz ) + &
       TWO * gupxy * ( fxy - F3o2/chin1 * chix * chiy ) + &
       TWO * gupxz * ( fxz - F3o2/chin1 * chix * chiz ) + &
       TWO * gupyz * ( fyz - F3o2/chin1 * chiy * chiz ) 
! Add chi part to Ricci tensor:

  Rxx = Rxx + (fxx - chix*chix/chin1/TWO + gxx * f)/chin1/TWO
  Ryy = Ryy + (fyy - chiy*chiy/chin1/TWO + gyy * f)/chin1/TWO
  Rzz = Rzz + (fzz - chiz*chiz/chin1/TWO + gzz * f)/chin1/TWO
  Rxy = Rxy + (fxy - chix*chiy/chin1/TWO + gxy * f)/chin1/TWO
  Rxz = Rxz + (fxz - chix*chiz/chin1/TWO + gxz * f)/chin1/TWO
  Ryz = Ryz + (fyz - chiy*chiz/chin1/TWO + gyz * f)/chin1/TWO

! covariant second derivatives of the lapse respect to physical metric
  call fdderivs(ex,Lap,fxx,fxy,fxz,fyy,fyz,fzz,X,Y,Z, &
                SYM,SYM,SYM,symmetry,Lev)

  gxxx = (gupxx * chix + gupxy * chiy + gupxz * chiz)/chin1
  gxxy = (gupxy * chix + gupyy * chiy + gupyz * chiz)/chin1
  gxxz = (gupxz * chix + gupyz * chiy + gupzz * chiz)/chin1
! now get physical second kind of connection
  Gamxxx = Gamxxx - ( (chix + chix)/chin1 - gxx * gxxx )*HALF
  Gamyxx = Gamyxx - (                     - gxx * gxxy )*HALF
  Gamzxx = Gamzxx - (                     - gxx * gxxz )*HALF
  Gamxyy = Gamxyy - (                     - gyy * gxxx )*HALF
  Gamyyy = Gamyyy - ( (chiy + chiy)/chin1 - gyy * gxxy )*HALF
  Gamzyy = Gamzyy - (                     - gyy * gxxz )*HALF
  Gamxzz = Gamxzz - (                     - gzz * gxxx )*HALF
  Gamyzz = Gamyzz - (                     - gzz * gxxy )*HALF
  Gamzzz = Gamzzz - ( (chiz + chiz)/chin1 - gzz * gxxz )*HALF
  Gamxxy = Gamxxy - (  chiy        /chin1 - gxy * gxxx )*HALF
  Gamyxy = Gamyxy - (         chix /chin1 - gxy * gxxy )*HALF
  Gamzxy = Gamzxy - (                     - gxy * gxxz )*HALF
  Gamxxz = Gamxxz - (  chiz        /chin1 - gxz * gxxx )*HALF
  Gamyxz = Gamyxz - (                     - gxz * gxxy )*HALF
  Gamzxz = Gamzxz - (         chix /chin1 - gxz * gxxz )*HALF
  Gamxyz = Gamxyz - (                     - gyz * gxxx )*HALF
  Gamyyz = Gamyyz - (  chiz        /chin1 - gyz * gxxy )*HALF
  Gamzyz = Gamzyz - (         chiy /chin1 - gyz * gxxz )*HALF

  fxx = fxx - Gamxxx*Lapx - Gamyxx*Lapy - Gamzxx*Lapz
  fyy = fyy - Gamxyy*Lapx - Gamyyy*Lapy - Gamzyy*Lapz
  fzz = fzz - Gamxzz*Lapx - Gamyzz*Lapy - Gamzzz*Lapz
  fxy = fxy - Gamxxy*Lapx - Gamyxy*Lapy - Gamzxy*Lapz
  fxz = fxz - Gamxxz*Lapx - Gamyxz*Lapy - Gamzxz*Lapz
  fyz = fyz - Gamxyz*Lapx - Gamyyz*Lapy - Gamzyz*Lapz

! store D^i D_i Lap in trK_rhs upto chi
  trK_rhs =    gupxx * fxx + gupyy * fyy + gupzz * fzz + &
        TWO* ( gupxy * fxy + gupxz * fxz + gupyz * fyz )
#if 1        
!! follow bam code
  S =  chin1 * ( gupxx * Sxx + gupyy * Syy + gupzz * Szz + &
     TWO * ( gupxy * Sxy + gupxz * Sxz + gupyz * Syz ) )
  f = F2o3 * trK * trK -(&
       gupxx * ( &
       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz) ) + &
       gupyy * ( &
       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz) ) + &
       gupzz * ( &
       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz) ) + &
       TWO * ( &
       gupxy * ( &
       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
       gupxy * (Axx * Ayy + Axy * Axy) + &
       gupxz * (Axx * Ayz + Axz * Axy) + &
       gupyz * (Axy * Ayz + Axz * Ayy) ) + &
       gupxz * ( &
       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
       gupxy * (Axx * Ayz + Axy * Axz) + &
       gupxz * (Axx * Azz + Axz * Axz) + &
       gupyz * (Axy * Azz + Axz * Ayz) ) + &
       gupyz * ( &
       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
       gupxy * (Axy * Ayz + Ayy * Axz) + &
       gupxz * (Axy * Azz + Ayz * Axz) + &
       gupyz * (Ayy * Azz + Ayz * Ayz) ) )) -1.6d1*PI*rho + EIGHT * PI * S
  f = - F1o3 *(  gupxx * fxx + gupyy * fyy + gupzz * fzz + &
        TWO* ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) + alpn1/chin1*f)
  
  fxx = alpn1 * (Rxx - EIGHT * PI * Sxx) - fxx
  fxy = alpn1 * (Rxy - EIGHT * PI * Sxy) - fxy
  fxz = alpn1 * (Rxz - EIGHT * PI * Sxz) - fxz
  fyy = alpn1 * (Ryy - EIGHT * PI * Syy) - fyy
  fyz = alpn1 * (Ryz - EIGHT * PI * Syz) - fyz
  fzz = alpn1 * (Rzz - EIGHT * PI * Szz) - fzz
#else        
! Add lapse and S_ij parts to Ricci tensor:

  fxx = alpn1 * (Rxx - EIGHT * PI * Sxx) - fxx
  fxy = alpn1 * (Rxy - EIGHT * PI * Sxy) - fxy
  fxz = alpn1 * (Rxz - EIGHT * PI * Sxz) - fxz
  fyy = alpn1 * (Ryy - EIGHT * PI * Syy) - fyy
  fyz = alpn1 * (Ryz - EIGHT * PI * Syz) - fyz
  fzz = alpn1 * (Rzz - EIGHT * PI * Szz) - fzz

! Compute trace-free part (note: chi^-1 and chi cancel!):

  f = F1o3 *(  gupxx * fxx + gupyy * fyy + gupzz * fzz + &
        TWO* ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) )
#endif

  Axx_rhs = fxx - gxx * f
  Ayy_rhs = fyy - gyy * f
  Azz_rhs = fzz - gzz * f
  Axy_rhs = fxy - gxy * f
  Axz_rhs = fxz - gxz * f
  Ayz_rhs = fyz - gyz * f

! Now: store A_il A^l_j into fij:

  fxx =       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz)
  fyy =       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz)
  fzz =       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz)
  fxy =       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
              gupxy *(Axx * Ayy + Axy * Axy)                            + &
              gupxz *(Axx * Ayz + Axz * Axy)                            + &
              gupyz *(Axy * Ayz + Axz * Ayy)
  fxz =       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
              gupxy *(Axx * Ayz + Axy * Axz)                            + &
              gupxz *(Axx * Azz + Axz * Axz)                            + &
              gupyz *(Axy * Azz + Axz * Ayz)
  fyz =       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
              gupxy *(Axy * Ayz + Ayy * Axz)                            + &
              gupxz *(Axy * Azz + Ayz * Axz)                            + &
              gupyz *(Ayy * Azz + Ayz * Ayz)

  f = chin1
! store D^i D_i Lap in trK_rhs
  trK_rhs = f*trK_rhs
          
  Axx_rhs =           f * Axx_rhs+ alpn1 * (trK * Axx - TWO * fxx)  + &
           TWO * (  Axx * betaxx +   Axy * betayx +   Axz * betazx )- &
             F2o3 * Axx * div_beta

  Ayy_rhs =           f * Ayy_rhs+ alpn1 * (trK * Ayy - TWO * fyy)  + &
           TWO * (  Axy * betaxy +   Ayy * betayy +   Ayz * betazy )- &
             F2o3 * Ayy * div_beta

  Azz_rhs =           f * Azz_rhs+ alpn1 * (trK * Azz - TWO * fzz)  + &
           TWO * (  Axz * betaxz +   Ayz * betayz +   Azz * betazz )- &
             F2o3 * Azz * div_beta

  Axy_rhs =           f * Axy_rhs+ alpn1 *( trK * Axy  - TWO * fxy )+ &
                    Axx * betaxy                  +   Axz * betazy  + &
                                     Ayy * betayx +   Ayz * betazx  + &
             F1o3 * Axy * div_beta                -   Axy * betazz

  Ayz_rhs =           f * Ayz_rhs+ alpn1 *( trK * Ayz  - TWO * fyz )+ &
                    Axy * betaxz +   Ayy * betayz                   + &
                    Axz * betaxy                  +   Azz * betazy  + &
             F1o3 * Ayz * div_beta                -   Ayz * betaxx
 
  Axz_rhs =           f * Axz_rhs+ alpn1 *( trK * Axz  - TWO * fxz )+ &
                    Axx * betaxz +   Axy * betayz                   + &
                                     Ayz * betayx +   Azz * betazx  + &
             F1o3 * Axz * div_beta                -   Axz * betayy      !rhs for Aij

! Compute trace of S_ij

  S =  f * ( gupxx * Sxx + gupyy * Syy + gupzz * Szz + &
     TWO * ( gupxy * Sxy + gupxz * Sxz + gupyz * Syz ) )

  trK_rhs = - trK_rhs + alpn1 *( F1o3 * trK * trK         + &
                gupxx * fxx + gupyy * fyy + gupzz * fzz   + &
        TWO * ( gupxy * fxy + gupxz * fxz + gupyz * fyz ) + &
       FOUR * PI * ( rho + S ))                                !rhs for trK
  
!!!! gauge variable part

  Lap_rhs = -TWO*alpn1*trK
  betax_rhs = FF*dtSfx
  betay_rhs = FF*dtSfy
  betaz_rhs = FF*dtSfz

  dtSfx_rhs = Gamx_rhs - eta*dtSfx
  dtSfy_rhs = Gamy_rhs - eta*dtSfy
  dtSfz_rhs = Gamz_rhs - eta*dtSfz

  SSS(1)=SYM
  SSS(2)=SYM
  SSS(3)=SYM

  AAS(1)=ANTI
  AAS(2)=ANTI
  AAS(3)=SYM

  ASA(1)=ANTI
  ASA(2)=SYM
  ASA(3)=ANTI

  SAA(1)=SYM
  SAA(2)=ANTI
  SAA(3)=ANTI

  ASS(1)=ANTI
  ASS(2)=SYM
  ASS(3)=SYM

  SAS(1)=SYM
  SAS(2)=ANTI
  SAS(3)=SYM

  SSA(1)=SYM
  SSA(2)=SYM
  SSA(3)=ANTI

!!!!!!!!!advection term part

  call lopsided(ex,X,Y,Z,gxx,gxx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gxy,gxy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,gxz,gxz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,gyy,gyy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gyz,gyz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,gzz,gzz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Axx,Axx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Axy,Axy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,Axz,Axz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,Ayy,Ayy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Ayz,Ayz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,Azz,Azz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,chi,chi_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,trK,trK_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Gamx,Gamx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,Gamy,Gamy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,Gamz,Gamz_rhs,betax,betay,betaz,Symmetry,SSA)
!!
  call lopsided(ex,X,Y,Z,Lap,Lap_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,betax,betax_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,betay,betay_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,betaz,betaz_rhs,betax,betay,betaz,Symmetry,SSA)

  call lopsided(ex,X,Y,Z,dtSfx,dtSfx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,dtSfy,dtSfy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,dtSfz,dtSfz_rhs,betax,betay,betaz,Symmetry,SSA)

  if(eps>0)then 
! usual Kreiss-Oliger dissipation      
  call kodis(ex,X,Y,Z,chi,chi_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,trK,trK_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dxx,gxx_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxy,gxy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxz,gxz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dyy,gyy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gyz,gyz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dzz,gzz_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axx,Axx_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axy,Axy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axz,Axz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayy,Ayy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayz,Ayz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Azz,Azz_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamx,Gamx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamy,Gamy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamz,Gamz_rhs,SSA,Symmetry,eps)

#if 1 
!! bam does not apply dissipation on gauge variables
  call kodis(ex,X,Y,Z,Lap,Lap_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betax,betax_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betay,betay_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betaz,betaz_rhs,SSA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfx,dtSfx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfy,dtSfy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfz,dtSfz_rhs,SSA,Symmetry,eps)
#endif

  endif

  if(co == 0)then
! ham_Res = trR + 2/3 * K^2 - A_ij * A^ij - 16 * PI * rho
! here trR is respect to physical metric
  ham_Res =   gupxx * Rxx + gupyy * Ryy + gupzz * Rzz + &
        TWO* ( gupxy * Rxy + gupxz * Rxz + gupyz * Ryz )

  ham_Res = chin1*ham_Res + F2o3 * trK * trK -(&
       gupxx * ( &
       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz) ) + &
       gupyy * ( &
       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz) ) + &
       gupzz * ( &
       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz) ) + &
       TWO * ( &
       gupxy * ( &
       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
       gupxy * (Axx * Ayy + Axy * Axy) + &
       gupxz * (Axx * Ayz + Axz * Axy) + &
       gupyz * (Axy * Ayz + Axz * Ayy) ) + &
       gupxz * ( &
       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
       gupxy * (Axx * Ayz + Axy * Axz) + &
       gupxz * (Axx * Azz + Axz * Axz) + &
       gupyz * (Axy * Azz + Axz * Ayz) ) + &
       gupyz * ( &
       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
       gupxy * (Axy * Ayz + Ayy * Axz) + &
       gupxz * (Axy * Azz + Ayz * Axz) + &
       gupyz * (Ayy * Azz + Ayz * Ayz) ) ))- F16 * PI * rho

! mov_Res_j = gupkj*(-1/chi d_k chi*A_ij + D_k A_ij) - 2/3 d_j trK - 8 PI s_j where D respect to physical metric
! store D_i A_jk - 1/chi d_i chi*A_jk in gjk_i
  call fderivs(ex,Axx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Axy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0)
  call fderivs(ex,Axz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0)
  call fderivs(ex,Ayy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Ayz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0)
  call fderivs(ex,Azz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)

  gxxx = gxxx - (  Gamxxx * Axx + Gamyxx * Axy + Gamzxx * Axz &
                 + Gamxxx * Axx + Gamyxx * Axy + Gamzxx * Axz) - chix*Axx/chin1
  gxyx = gxyx - (  Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz &
                 + Gamxxx * Axy + Gamyxx * Ayy + Gamzxx * Ayz) - chix*Axy/chin1
  gxzx = gxzx - (  Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz &
                 + Gamxxx * Axz + Gamyxx * Ayz + Gamzxx * Azz) - chix*Axz/chin1
  gyyx = gyyx - (  Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz &
                 + Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz) - chix*Ayy/chin1
  gyzx = gyzx - (  Gamxxz * Axy + Gamyxz * Ayy + Gamzxz * Ayz &
                 + Gamxxy * Axz + Gamyxy * Ayz + Gamzxy * Azz) - chix*Ayz/chin1
  gzzx = gzzx - (  Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz &
                 + Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz) - chix*Azz/chin1
  gxxy = gxxy - (  Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz &
                 + Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz) - chiy*Axx/chin1
  gxyy = gxyy - (  Gamxyy * Axx + Gamyyy * Axy + Gamzyy * Axz &
                 + Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz) - chiy*Axy/chin1
  gxzy = gxzy - (  Gamxyz * Axx + Gamyyz * Axy + Gamzyz * Axz &
                 + Gamxxy * Axz + Gamyxy * Ayz + Gamzxy * Azz) - chiy*Axz/chin1
  gyyy = gyyy - (  Gamxyy * Axy + Gamyyy * Ayy + Gamzyy * Ayz &
                 + Gamxyy * Axy + Gamyyy * Ayy + Gamzyy * Ayz) - chiy*Ayy/chin1
  gyzy = gyzy - (  Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz &
                 + Gamxyy * Axz + Gamyyy * Ayz + Gamzyy * Azz) - chiy*Ayz/chin1
  gzzy = gzzy - (  Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz &
                 + Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz) - chiy*Azz/chin1
  gxxz = gxxz - (  Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz &
                 + Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz) - chiz*Axx/chin1
  gxyz = gxyz - (  Gamxyz * Axx + Gamyyz * Axy + Gamzyz * Axz &
                 + Gamxxz * Axy + Gamyxz * Ayy + Gamzxz * Ayz) - chiz*Axy/chin1
  gxzz = gxzz - (  Gamxzz * Axx + Gamyzz * Axy + Gamzzz * Axz &
                 + Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz) - chiz*Axz/chin1
  gyyz = gyyz - (  Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz &
                 + Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz) - chiz*Ayy/chin1
  gyzz = gyzz - (  Gamxzz * Axy + Gamyzz * Ayy + Gamzzz * Ayz &
                 + Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz) - chiz*Ayz/chin1
  gzzz = gzzz - (  Gamxzz * Axz + Gamyzz * Ayz + Gamzzz * Azz &
                 + Gamxzz * Axz + Gamyzz * Ayz + Gamzzz * Azz) - chiz*Azz/chin1
movx_Res = gupxx*gxxx + gupyy*gxyy + gupzz*gxzz &
          +gupxy*gxyx + gupxz*gxzx + gupyz*gxzy &
          +gupxy*gxxy + gupxz*gxxz + gupyz*gxyz
movy_Res = gupxx*gxyx + gupyy*gyyy + gupzz*gyzz &
          +gupxy*gyyx + gupxz*gyzx + gupyz*gyzy &
          +gupxy*gxyy + gupxz*gxyz + gupyz*gyyz
movz_Res = gupxx*gxzx + gupyy*gyzy + gupzz*gzzz &
          +gupxy*gyzx + gupxz*gzzx + gupyz*gzzy &
          +gupxy*gxzy + gupxz*gxzz + gupyz*gyzz

movx_Res = movx_Res - F2o3*Kx - F8*PI*sx
movy_Res = movy_Res - F2o3*Ky - F8*PI*sy
movz_Res = movz_Res - F2o3*Kz - F8*PI*sz
  endif


  gont = 0

  return

  end function compute_rhs_bssn
!=====================================================================
! compute_rhs_bssn_fused : single-point fused version of compute_rhs_bssn
! (bssn_rhs.f90 lines 5-980).
!
! The assembly part of the original (lines 134-817) is evaluated inside one
! explicit k/j/i loop with scalar temporaries (register fusion); the 43+
! temporary 3D arrays of the original are eliminated. Every stencil call
! (fderivs/fdderivs/lopsided/kodis) is unchanged and runs on full arrays
! outside the loop.
!
! Generation renaming (each stencil call keeps its exact input->output
! mapping; fderivs/fdderivs zero their outputs first, so all grid points,
! incl. boundaries, are fully defined per call):
!   gxxx..gzzz  : 1st generation, fderivs(dxx,gxy,gxz,dyy,gyz,dzz)  -> used by
!                  Gmx/Gmy/Gmz_Res (co==0) and Gamxxx..Gamzzz conformal conn.
!   bxxx..bzzz  : 2nd generation, fdderivs(betax,betay,betaz)       -> used by
!                  fxx,fxy,fxz and Gamx/Gamy/Gamz_rhs += terms
!   gxxx..gzzz  : re-written (whole array) inside the final if(co==0) block by
!                  fderivs(Axx..Azz) -- verbatim, after lopsided/kodis
!   fxx_dxx..fzz_dxx ... fzz_lap : fdderivs(dxx,dyy,dzz,gxy,gxz,gyz,chi,Lap)
!                  one dedicated 6-array set per call (Rxx..Rzz, chi part,
!                  Lap part)
! Loop-internal scalar generations of gxxx..gzzz:
!   gxxx_t..gzzz_t : 3rd generation, first-kind connection (lines 344-363)
!   gxxx2_s,gxxy2_s,gxxz2_s : chi-connection vector (lines 620-622)
!
! Arrays stored from the loop (read by the final if(co==0) block and
! by lopsided): gupxx..gupzz, Rxx..Rzz, Gamxxx..Gamzzz, gxx,gyy,gzz, chin1.
! With AMSS_RHS_SKIP_CONSTRAINT_STORES, the first three groups are written
! only for co==0; their scalar values are still used by the co==1 RHS.
!=====================================================================

  function compute_rhs_bssn_fused(ex, T,X, Y, Z,                                     &
               chi    ,   trK    ,                                             &
               dxx    ,   gxy    ,   gxz    ,   dyy    ,   gyz    ,   dzz,     &
               Axx    ,   Axy    ,   Axz    ,   Ayy    ,   Ayz    ,   Azz,     &
               Gamx   ,  Gamy    ,  Gamz    ,                                  &
               Lap    ,  betax   ,  betay   ,  betaz   ,                       &
               dtSfx  ,  dtSfy   ,  dtSfz   ,                                  &
               chi_rhs,   trK_rhs,                                             &
               gxx_rhs,   gxy_rhs,   gxz_rhs,   gyy_rhs,   gyz_rhs,   gzz_rhs, &
               Axx_rhs,   Axy_rhs,   Axz_rhs,   Ayy_rhs,   Ayz_rhs,   Azz_rhs, &
               Gamx_rhs,  Gamy_rhs,  Gamz_rhs,                                 &
               Lap_rhs,  betax_rhs,  betay_rhs,  betaz_rhs,                    &
               dtSfx_rhs,  dtSfy_rhs,  dtSfz_rhs,                              &
               rho,Sx,Sy,Sz,Sxx,Sxy,Sxz,Syy,Syz,Szz,                           &
               Gamxxx,Gamxxy,Gamxxz,Gamxyy,Gamxyz,Gamxzz,                      &
               Gamyxx,Gamyxy,Gamyxz,Gamyyy,Gamyyz,Gamyzz,                      &
               Gamzxx,Gamzxy,Gamzxz,Gamzyy,Gamzyz,Gamzzz,                      &
               Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,                                        &
               ham_Res, movx_Res, movy_Res, movz_Res,                          &
                        Gmx_Res, Gmy_Res, Gmz_Res,                             &
               Symmetry,Lev,eps,co)  result(gont)
! calculate constraint violation when co=0
#ifdef AMSS_RHS_POINTWISE
  use point_derivs
#endif
#ifdef AMSS_RHS_OMP_ASSEMBLY
! omp_get_max_threads() in the parallel-do if() clause needs an explicit
! declaration under implicit none; omp_lib provides it (linked via -fopenmp).
  use omp_lib
#endif
#ifdef AMSS_RHS_WORKSPACE_POOL
  use rhs_legacy_buffers, only: legacy_ws, rhs_legacy_ensure
#endif
#ifdef AMSS_BATCH_STENCIL
  use batch_stencils, only: field_ptr, bs_enabled, bs_size, kodis_batch, lopsided_batch
#endif
  implicit none

!~~~~~~> Input parameters:

  integer,intent(in ):: ex(1:3), Symmetry,Lev,co
  real*8, intent(in ):: T
  real*8, intent(in ):: X(1:ex(1)),Y(1:ex(2)),Z(1:ex(3))
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout),target :: chi,dxx,dyy,dzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ),target :: trK
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ),target :: gxy,gxz,gyz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ),target :: Axx,Axy,Axz,Ayy,Ayz,Azz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ),target :: Gamx,Gamy,Gamz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout),target :: Lap, betax, betay, betaz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ),target :: dtSfx,  dtSfy,  dtSfz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out),target :: chi_rhs,trK_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out),target :: gxx_rhs,gxy_rhs,gxz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out),target :: gyy_rhs,gyz_rhs,gzz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out),target :: Axx_rhs,Axy_rhs,Axz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out),target :: Ayy_rhs,Ayz_rhs,Azz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out),target :: Gamx_rhs,Gamy_rhs,Gamz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out),target :: Lap_rhs, betax_rhs, betay_rhs, betaz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out),target :: dtSfx_rhs,dtSfy_rhs,dtSfz_rhs
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: rho,Sx,Sy,Sz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(in ) :: Sxx,Sxy,Sxz,Syy,Syz,Szz
! when out, physical second kind of connection
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxxx, Gamxxy, Gamxxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamxyy, Gamxyz, Gamxzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyxx, Gamyxy, Gamyxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamyyy, Gamyyz, Gamyzz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzxx, Gamzxy, Gamzxz
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Gamzyy, Gamzyz, Gamzzz
! when out, physical Ricci tensor
  real*8, dimension(ex(1),ex(2),ex(3)),intent(out) :: Rxx,Rxy,Rxz,Ryy,Ryz,Rzz
  real*8,intent(in) :: eps
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: ham_Res, movx_Res, movy_Res, movz_Res
  real*8, dimension(ex(1),ex(2),ex(3)),intent(inout) :: Gmx_Res, Gmy_Res, Gmz_Res
!  gont = 0: success; gont = 1: something wrong
  integer::gont
#ifdef AMSS_RHS_POINTWISE
  logical :: pointwise_mode
#endif
#ifdef AMSS_BATCH_STENCIL
#ifdef AMSS_RHS_RAW_DIAG_LOPSIDED
  logical :: raw_diag_lopsided
#endif
#endif

!~~~~~~> Other variables:

! whole-array temporaries that must survive the loop:
!   gxx,gyy,gzz : read by lopsided after the loop
!   chin1       : read by final if(co==0) block
!   gupxx..gupzz: stored from the loop, read by final if(co==0) block
!   Rxx..Rzz, Gamxxx..Gamzzz : output arrays, final values stored from loop
  real*8, dimension(ex(1),ex(2),ex(3)),target :: gxx,gyy,gzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: chin1
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupxx,gupxy,gupxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gupyy,gupyz,gupzz

#ifdef AMSS_RHS_WORKSPACE_POOL
! stencil outputs (all written before the fused loop).  These are pointers
! into rhs_legacy_buffers; they remain unassociated on the POINTWISE hot
! path, where the preprocessor maps all of their uses to point-local scalars.
  real*8, pointer :: chix(:,:,:),chiy(:,:,:),chiz(:,:,:)
  real*8, pointer :: gxxx(:,:,:),gxyx(:,:,:),gxzx(:,:,:),gyyx(:,:,:),gyzx(:,:,:),gzzx(:,:,:)
  real*8, pointer :: gxxy(:,:,:),gxyy(:,:,:),gxzy(:,:,:),gyyy(:,:,:),gyzy(:,:,:),gzzy(:,:,:)
  real*8, pointer :: gxxz(:,:,:),gxyz(:,:,:),gxzz(:,:,:),gyyz(:,:,:),gyzz(:,:,:),gzzz(:,:,:)
  real*8, pointer :: Lapx(:,:,:),Lapy(:,:,:),Lapz(:,:,:)
  real*8, pointer :: betaxx(:,:,:),betaxy(:,:,:),betaxz(:,:,:)
  real*8, pointer :: betayx(:,:,:),betayy(:,:,:),betayz(:,:,:)
  real*8, pointer :: betazx(:,:,:),betazy(:,:,:),betazz(:,:,:)
  real*8, pointer :: Gamxx(:,:,:),Gamxy(:,:,:),Gamxz(:,:,:)
  real*8, pointer :: Gamyx(:,:,:),Gamyy(:,:,:),Gamyz(:,:,:)
  real*8, pointer :: Gamzx(:,:,:),Gamzy(:,:,:),Gamzz(:,:,:)
  real*8, pointer :: Kx(:,:,:),Ky(:,:,:),Kz(:,:,:)

! 2nd generation of gxxx..gzzz: fdderivs(betax/betay/betaz) outputs
! (original wrote them into gxxx..gzzz after the Gamxxx..Gamzzz block;
!  both generations are read inside the loop, so the names are split)
  real*8, pointer :: bxxx(:,:,:),bxyx(:,:,:),bxzx(:,:,:),byyx(:,:,:),byzx(:,:,:),bzzx(:,:,:)
  real*8, pointer :: bxxy(:,:,:),bxyy(:,:,:),bxzy(:,:,:),byyy(:,:,:),byzy(:,:,:),bzzy(:,:,:)
  real*8, pointer :: bxxz(:,:,:),bxyz(:,:,:),bxzz(:,:,:),byyz(:,:,:),byzz(:,:,:),bzzz(:,:,:)

! fdderivs outputs: one dedicated 6-array set per field
  real*8, pointer :: fxx_dxx(:,:,:),fxy_dxx(:,:,:),fxz_dxx(:,:,:),fyy_dxx(:,:,:),fyz_dxx(:,:,:),fzz_dxx(:,:,:)
  real*8, pointer :: fxx_dyy(:,:,:),fxy_dyy(:,:,:),fxz_dyy(:,:,:),fyy_dyy(:,:,:),fyz_dyy(:,:,:),fzz_dyy(:,:,:)
  real*8, pointer :: fxx_dzz(:,:,:),fxy_dzz(:,:,:),fxz_dzz(:,:,:),fyy_dzz(:,:,:),fyz_dzz(:,:,:),fzz_dzz(:,:,:)
  real*8, pointer :: fxx_gxy(:,:,:),fxy_gxy(:,:,:),fxz_gxy(:,:,:),fyy_gxy(:,:,:),fyz_gxy(:,:,:),fzz_gxy(:,:,:)
  real*8, pointer :: fxx_gxz(:,:,:),fxy_gxz(:,:,:),fxz_gxz(:,:,:),fyy_gxz(:,:,:),fyz_gxz(:,:,:),fzz_gxz(:,:,:)
  real*8, pointer :: fxx_gyz(:,:,:),fxy_gyz(:,:,:),fxz_gyz(:,:,:),fyy_gyz(:,:,:),fyz_gyz(:,:,:),fzz_gyz(:,:,:)
  real*8, pointer :: fxx_chi(:,:,:),fxy_chi(:,:,:),fxz_chi(:,:,:),fyy_chi(:,:,:),fyz_chi(:,:,:),fzz_chi(:,:,:)
  real*8, pointer :: fxx_lap(:,:,:),fxy_lap(:,:,:),fxz_lap(:,:,:),fyy_lap(:,:,:),fyz_lap(:,:,:),fzz_lap(:,:,:)
#else
! Original per-call arrays retained for the default/reference build.
  real*8, dimension(ex(1),ex(2),ex(3)) :: chix,chiy,chiz
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxx,gxyx,gxzx,gyyx,gyzx,gzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxy,gxyy,gxzy,gyyy,gyzy,gzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: gxxz,gxyz,gxzz,gyyz,gyzz,gzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Lapx,Lapy,Lapz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betaxx,betaxy,betaxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betayx,betayy,betayz
  real*8, dimension(ex(1),ex(2),ex(3)) :: betazx,betazy,betazz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamxx,Gamxy,Gamxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamyx,Gamyy,Gamyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Gamzx,Gamzy,Gamzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: Kx,Ky,Kz
  real*8, dimension(ex(1),ex(2),ex(3)) :: bxxx,bxyx,bxzx,byyx,byzx,bzzx
  real*8, dimension(ex(1),ex(2),ex(3)) :: bxxy,bxyy,bxzy,byyy,byzy,bzzy
  real*8, dimension(ex(1),ex(2),ex(3)) :: bxxz,bxyz,bxzz,byyz,byzz,bzzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: fxx_dxx,fxy_dxx,fxz_dxx,fyy_dxx,fyz_dxx,fzz_dxx
  real*8, dimension(ex(1),ex(2),ex(3)) :: fxx_dyy,fxy_dyy,fxz_dyy,fyy_dyy,fyz_dyy,fzz_dyy
  real*8, dimension(ex(1),ex(2),ex(3)) :: fxx_dzz,fxy_dzz,fxz_dzz,fyy_dzz,fyz_dzz,fzz_dzz
  real*8, dimension(ex(1),ex(2),ex(3)) :: fxx_gxy,fxy_gxy,fxz_gxy,fyy_gxy,fyz_gxy,fzz_gxy
  real*8, dimension(ex(1),ex(2),ex(3)) :: fxx_gxz,fxy_gxz,fxz_gxz,fyy_gxz,fyz_gxz,fzz_gxz
  real*8, dimension(ex(1),ex(2),ex(3)) :: fxx_gyz,fxy_gyz,fxz_gyz,fyy_gyz,fyz_gyz,fzz_gyz
  real*8, dimension(ex(1),ex(2),ex(3)) :: fxx_chi,fxy_chi,fxz_chi,fyy_chi,fyz_chi,fzz_chi
  real*8, dimension(ex(1),ex(2),ex(3)) :: fxx_lap,fxy_lap,fxz_lap,fyy_lap,fyz_lap,fzz_lap
#endif

! scalar temporaries of the fused loop (replace the original 3D temporaries
! alpn1, div_beta, gupxx..gupzz(loop use), Gamxxx..Gamzzz, Rxx..Rzz,
! Gamx_rhs, Gamy_rhs, Gamz_rhs, Axx_rhs..Azz_rhs, fxx..fzz, f, S, Gamxa..Gamza)
  real*8 :: alpn1_s, chin1_s, gxx_s, gyy_s, gzz_s, div_beta_s
  real*8 :: gupdet_s, gupxx_s, gupxy_s, gupxz_s, gupyy_s, gupyz_s, gupzz_s
  real*8 :: Gamxxx_s, Gamxxy_s, Gamxxz_s, Gamxyy_s, Gamxyz_s, Gamxzz_s
  real*8 :: Gamyxx_s, Gamyxy_s, Gamyxz_s, Gamyyy_s, Gamyyz_s, Gamyzz_s
  real*8 :: Gamzxx_s, Gamzxy_s, Gamzxz_s, Gamzyy_s, Gamzyz_s, Gamzzz_s
  real*8 :: Rxx_s, Rxy_s, Rxz_s, Ryy_s, Ryz_s, Rzz_s
  real*8 :: Gamx_rhs_s, Gamy_rhs_s, Gamz_rhs_s, trK_rhs_s
  real*8 :: Axx_rhs_s, Axy_rhs_s, Axz_rhs_s, Ayy_rhs_s, Ayz_rhs_s, Azz_rhs_s
  real*8 :: fxx_s, fxy_s, fxz_s, fyy_s, fyz_s, fzz_s
  real*8 :: f_s, S_s
  real*8 :: Gamxa_s, Gamya_s, Gamza_s
! 3rd generation of gxxx..gzzz: first-kind connection (original lines 344-363)
  real*8 :: gxxx_t, gxyx_t, gxzx_t, gyyx_t, gyzx_t, gzzx_t
  real*8 :: gxxy_t, gxyy_t, gxzy_t, gyyy_t, gyzy_t, gzzy_t
  real*8 :: gxxz_t, gxyz_t, gxzz_t, gyyz_t, gyzz_t, gzzz_t
! chi-connection vector (original lines 620-622)
  real*8 :: gxxx2_s, gxxy2_s, gxxz2_s

  real*8,dimension(3) ::SSS,AAS,ASA,SAA,ASS,SAS,SSA
  real*8            :: dX, dY, dZ, PI
  integer           :: i,j,k
#ifdef AMSS_RHS_TILED
  integer           :: jb,kb,jend,kend
#endif
#ifdef AMSS_RHS_POINTWISE
  integer           :: pw_imin, pw_jmin, pw_kmin, pw_order
  real*8            :: pw_d12x, pw_d12y, pw_d12z
  real*8            :: pw_d2x, pw_d2y, pw_d2z
  real*8            :: pw_sxx, pw_syy, pw_szz
  real*8            :: pw_fxxc, pw_fxyc, pw_fxzc, pw_fyyc, pw_fyzc, pw_fzzc
! Pointwise second-derivative results stay in registers/stack scalars.  The
! LEGACY/HALO paths still use the full-array fdderivs outputs above.
  real*8 :: pw_bxxx,pw_bxyx,pw_bxzx,pw_byyx,pw_byzx,pw_bzzx
  real*8 :: pw_bxxy,pw_bxyy,pw_bxzy,pw_byyy,pw_byzy,pw_bzzy
  real*8 :: pw_bxxz,pw_bxyz,pw_bxzz,pw_byyz,pw_byzz,pw_bzzz
  real*8 :: pw_fxx_dxx,pw_fxy_dxx,pw_fxz_dxx,pw_fyy_dxx,pw_fyz_dxx,pw_fzz_dxx
  real*8 :: pw_fxx_dyy,pw_fxy_dyy,pw_fxz_dyy,pw_fyy_dyy,pw_fyz_dyy,pw_fzz_dyy
  real*8 :: pw_fxx_dzz,pw_fxy_dzz,pw_fxz_dzz,pw_fyy_dzz,pw_fyz_dzz,pw_fzz_dzz
  real*8 :: pw_fxx_gxy,pw_fxy_gxy,pw_fxz_gxy,pw_fyy_gxy,pw_fyz_gxy,pw_fzz_gxy
  real*8 :: pw_fxx_gxz,pw_fxy_gxz,pw_fxz_gxz,pw_fyy_gxz,pw_fyz_gxz,pw_fzz_gxz
  real*8 :: pw_fxx_gyz,pw_fxy_gyz,pw_fxz_gyz,pw_fyy_gyz,pw_fyz_gyz,pw_fzz_gyz
  real*8 :: pw_fxx_chi,pw_fxy_chi,pw_fxz_chi,pw_fyy_chi,pw_fyz_chi,pw_fzz_chi
  real*8 :: pw_fxx_lap,pw_fxy_lap,pw_fxz_lap,pw_fyy_lap,pw_fyz_lap,pw_fzz_lap
! Pointwise first-derivative results.  They are consumed by the fused
! assembly at the same (i,j,k), so materializing them as full 3-D arrays is
! unnecessary in the POINTWISE path.
  real*8 :: pw_betaxx,pw_betaxy,pw_betaxz
  real*8 :: pw_betayx,pw_betayy,pw_betayz
  real*8 :: pw_betazx,pw_betazy,pw_betazz
  real*8 :: pw_chix,pw_chiy,pw_chiz
  real*8 :: pw_gxxx,pw_gxxy,pw_gxxz
  real*8 :: pw_gxyx,pw_gxyy,pw_gxyz
  real*8 :: pw_gxzx,pw_gxzy,pw_gxzz
  real*8 :: pw_gyyx,pw_gyyy,pw_gyyz
  real*8 :: pw_gyzx,pw_gyzy,pw_gyzz
  real*8 :: pw_gzzx,pw_gzzy,pw_gzzz
  real*8 :: pw_lapx,pw_lapy,pw_lapz
  real*8 :: pw_kx,pw_ky,pw_kz
  real*8 :: pw_gamxx,pw_gamxy,pw_gamxz
  real*8 :: pw_gamyx,pw_gamyy,pw_gamyz
  real*8 :: pw_gamzx,pw_gamzy,pw_gamzz
#endif
  real*8, parameter :: ZEO = 0.d0,ONE = 1.D0, TWO = 2.D0, FOUR = 4.D0
  real*8, parameter :: EIGHT = 8.D0, HALF = 0.5D0, THR = 3.d0
  real*8, parameter :: SYM = 1.D0, ANTI= - 1.D0
  double precision,parameter::FF = 0.75d0,eta=2.d0
  real*8, parameter :: F1o3 = 1.D0/3.D0, F2o3 = 2.D0/3.D0,F3o2=1.5d0, F1o6 = 1.D0/6.D0
  real*8, parameter :: F16=1.6d1,F8=8.d0



#ifdef AMSS_RHS_NAN_CHECK
!!! sanity check
  dX = sum(chi)+sum(trK)+sum(dxx)+sum(gxy)+sum(gxz)+sum(dyy)+sum(gyz)+sum(dzz) &
      +sum(Axx)+sum(Axy)+sum(Axz)+sum(Ayy)+sum(Ayz)+sum(Azz)                   &
      +sum(Gamx)+sum(Gamy)+sum(Gamz)                                           &
      +sum(Lap)+sum(betax)+sum(betay)+sum(betaz)
  if(dX.ne.dX) then
     if(sum(chi).ne.sum(chi))write(*,*)"bssn.f90: find NaN in chi"
     if(sum(trK).ne.sum(trK))write(*,*)"bssn.f90: find NaN in trk"
     if(sum(dxx).ne.sum(dxx))write(*,*)"bssn.f90: find NaN in dxx"
     if(sum(gxy).ne.sum(gxy))write(*,*)"bssn.f90: find NaN in gxy"
     if(sum(gxz).ne.sum(gxz))write(*,*)"bssn.f90: find NaN in gxz"
     if(sum(dyy).ne.sum(dyy))write(*,*)"bssn.f90: find NaN in dyy"
     if(sum(gyz).ne.sum(gyz))write(*,*)"bssn.f90: find NaN in gyz"
     if(sum(dzz).ne.sum(dzz))write(*,*)"bssn.f90: find NaN in dzz"
     if(sum(Axx).ne.sum(Axx))write(*,*)"bssn.f90: find NaN in Axx"
     if(sum(Axy).ne.sum(Axy))write(*,*)"bssn.f90: find NaN in Axy"
     if(sum(Axz).ne.sum(Axz))write(*,*)"bssn.f90: find NaN in Axz"
     if(sum(Ayy).ne.sum(Ayy))write(*,*)"bssn.f90: find NaN in Ayy"
     if(sum(Ayz).ne.sum(Ayz))write(*,*)"bssn.f90: find NaN in Ayz"
     if(sum(Azz).ne.sum(Azz))write(*,*)"bssn.f90: find NaN in Azz"
     if(sum(Gamx).ne.sum(Gamx))write(*,*)"bssn.f90: find NaN in Gamx"
     if(sum(Gamy).ne.sum(Gamy))write(*,*)"bssn.f90: find NaN in Gamy"
     if(sum(Gamz).ne.sum(Gamz))write(*,*)"bssn.f90: find NaN in Gamz"
     if(sum(Lap).ne.sum(Lap))write(*,*)"bssn.f90: find NaN in Lap"
     if(sum(betax).ne.sum(betax))write(*,*)"bssn.f90: find NaN in betax"
     if(sum(betay).ne.sum(betay))write(*,*)"bssn.f90: find NaN in betay"
     if(sum(betaz).ne.sum(betaz))write(*,*)"bssn.f90: find NaN in betaz"
     gont = 1
     return
  endif
#endif

  PI = dacos(-ONE)

  dX = X(2) - X(1)
  dY = Y(2) - Y(1)
  dZ = Z(2) - Z(1)
#ifdef AMSS_RHS_POINTWISE
  pointwise_mode = (co == 1 .and. Symmetry <= 1)
  pw_d12x = 1.d0/(12.d0*dX)
  pw_d12y = 1.d0/(12.d0*dY)
  pw_d12z = 1.d0/(12.d0*dZ)
  pw_d2x = 0.5d0/dX
  pw_d2y = 0.5d0/dY
  pw_d2z = 0.5d0/dZ
  pw_sxx = 1.d0/(dX*dX)
  pw_syy = 1.d0/(dY*dY)
  pw_szz = 1.d0/(dZ*dZ)
  pw_fxxc = 1.d0/(12.d0*dX*dX)
  pw_fyyc = 1.d0/(12.d0*dY*dY)
  pw_fzzc = 1.d0/(12.d0*dZ*dZ)
  pw_fxyc = 1.d0/(144.d0*dX*dY)
  pw_fxzc = 1.d0/(144.d0*dX*dZ)
  pw_fyzc = 1.d0/(144.d0*dY*dZ)
  pw_imin = 1
  pw_jmin = 1
  pw_kmin = 1
  if (Symmetry > 0 .and. dabs(Z(1)) < dZ) pw_kmin = -1
  if (Symmetry > 1 .and. dabs(X(1)) < dX) pw_imin = -1
  if (Symmetry > 1 .and. dabs(Y(1)) < dY) pw_jmin = -1
#endif

#ifdef AMSS_RHS_WORKSPACE_POOL
! The legacy branch is cold under POINTWISE.  Bind its full-array stencil
! outputs only when that branch is actually selected; this removes 111
! automatic-array allocations from every POINTWISE RHS call.
#ifdef AMSS_RHS_POINTWISE
  if (.not. pointwise_mode) then
#else
  if (.true.) then
#endif
    call rhs_legacy_ensure(ex)
    chix  => legacy_ws(:,:,:,1);   chiy  => legacy_ws(:,:,:,2);   chiz  => legacy_ws(:,:,:,3)
    gxxx  => legacy_ws(:,:,:,4);   gxyx  => legacy_ws(:,:,:,5);   gxzx  => legacy_ws(:,:,:,6)
    gyyx  => legacy_ws(:,:,:,7);   gyzx  => legacy_ws(:,:,:,8);   gzzx  => legacy_ws(:,:,:,9)
    gxxy  => legacy_ws(:,:,:,10);  gxyy  => legacy_ws(:,:,:,11);  gxzy  => legacy_ws(:,:,:,12)
    gyyy  => legacy_ws(:,:,:,13);  gyzy  => legacy_ws(:,:,:,14);  gzzy  => legacy_ws(:,:,:,15)
    gxxz  => legacy_ws(:,:,:,16);  gxyz  => legacy_ws(:,:,:,17);  gxzz  => legacy_ws(:,:,:,18)
    gyyz  => legacy_ws(:,:,:,19);  gyzz  => legacy_ws(:,:,:,20);  gzzz  => legacy_ws(:,:,:,21)
    Lapx  => legacy_ws(:,:,:,22);  Lapy  => legacy_ws(:,:,:,23);  Lapz  => legacy_ws(:,:,:,24)
    betaxx=> legacy_ws(:,:,:,25);  betaxy=> legacy_ws(:,:,:,26);  betaxz=> legacy_ws(:,:,:,27)
    betayx=> legacy_ws(:,:,:,28);  betayy=> legacy_ws(:,:,:,29);  betayz=> legacy_ws(:,:,:,30)
    betazx=> legacy_ws(:,:,:,31);  betazy=> legacy_ws(:,:,:,32);  betazz=> legacy_ws(:,:,:,33)
    Gamxx => legacy_ws(:,:,:,34);  Gamxy => legacy_ws(:,:,:,35);  Gamxz => legacy_ws(:,:,:,36)
    Gamyx => legacy_ws(:,:,:,37);  Gamyy => legacy_ws(:,:,:,38);  Gamyz => legacy_ws(:,:,:,39)
    Gamzx => legacy_ws(:,:,:,40);  Gamzy => legacy_ws(:,:,:,41);  Gamzz => legacy_ws(:,:,:,42)
    Kx    => legacy_ws(:,:,:,43);  Ky    => legacy_ws(:,:,:,44);  Kz    => legacy_ws(:,:,:,45)
    bxxx  => legacy_ws(:,:,:,46);  bxyx  => legacy_ws(:,:,:,47);  bxzx  => legacy_ws(:,:,:,48)
    byyx  => legacy_ws(:,:,:,49);  byzx  => legacy_ws(:,:,:,50);  bzzx  => legacy_ws(:,:,:,51)
    bxxy  => legacy_ws(:,:,:,52);  bxyy  => legacy_ws(:,:,:,53);  bxzy  => legacy_ws(:,:,:,54)
    byyy  => legacy_ws(:,:,:,55);  byzy  => legacy_ws(:,:,:,56);  bzzy  => legacy_ws(:,:,:,57)
    bxxz  => legacy_ws(:,:,:,58);  bxyz  => legacy_ws(:,:,:,59);  bxzz  => legacy_ws(:,:,:,60)
    byyz  => legacy_ws(:,:,:,61);  byzz  => legacy_ws(:,:,:,62);  bzzz  => legacy_ws(:,:,:,63)
    fxx_dxx=>legacy_ws(:,:,:,64);  fxy_dxx=>legacy_ws(:,:,:,65);  fxz_dxx=>legacy_ws(:,:,:,66)
    fyy_dxx=>legacy_ws(:,:,:,67);  fyz_dxx=>legacy_ws(:,:,:,68);  fzz_dxx=>legacy_ws(:,:,:,69)
    fxx_dyy=>legacy_ws(:,:,:,70);  fxy_dyy=>legacy_ws(:,:,:,71);  fxz_dyy=>legacy_ws(:,:,:,72)
    fyy_dyy=>legacy_ws(:,:,:,73);  fyz_dyy=>legacy_ws(:,:,:,74);  fzz_dyy=>legacy_ws(:,:,:,75)
    fxx_dzz=>legacy_ws(:,:,:,76);  fxy_dzz=>legacy_ws(:,:,:,77);  fxz_dzz=>legacy_ws(:,:,:,78)
    fyy_dzz=>legacy_ws(:,:,:,79);  fyz_dzz=>legacy_ws(:,:,:,80);  fzz_dzz=>legacy_ws(:,:,:,81)
    fxx_gxy=>legacy_ws(:,:,:,82);  fxy_gxy=>legacy_ws(:,:,:,83);  fxz_gxy=>legacy_ws(:,:,:,84)
    fyy_gxy=>legacy_ws(:,:,:,85);  fyz_gxy=>legacy_ws(:,:,:,86);  fzz_gxy=>legacy_ws(:,:,:,87)
    fxx_gxz=>legacy_ws(:,:,:,88);  fxy_gxz=>legacy_ws(:,:,:,89);  fxz_gxz=>legacy_ws(:,:,:,90)
    fyy_gxz=>legacy_ws(:,:,:,91);  fyz_gxz=>legacy_ws(:,:,:,92);  fzz_gxz=>legacy_ws(:,:,:,93)
    fxx_gyz=>legacy_ws(:,:,:,94);  fxy_gyz=>legacy_ws(:,:,:,95);  fxz_gyz=>legacy_ws(:,:,:,96)
    fyy_gyz=>legacy_ws(:,:,:,97);  fyz_gyz=>legacy_ws(:,:,:,98);  fzz_gyz=>legacy_ws(:,:,:,99)
    fxx_chi=>legacy_ws(:,:,:,100); fxy_chi=>legacy_ws(:,:,:,101); fxz_chi=>legacy_ws(:,:,:,102)
    fyy_chi=>legacy_ws(:,:,:,103); fyz_chi=>legacy_ws(:,:,:,104); fzz_chi=>legacy_ws(:,:,:,105)
    fxx_lap=>legacy_ws(:,:,:,106); fxy_lap=>legacy_ws(:,:,:,107); fxz_lap=>legacy_ws(:,:,:,108)
    fyy_lap=>legacy_ws(:,:,:,109); fyz_lap=>legacy_ws(:,:,:,110); fzz_lap=>legacy_ws(:,:,:,111)
#ifdef AMSS_RHS_POINTWISE
  end if
#endif
#endif

! whole-array inits needed after the loop (lopsided reads gxx,gyy,gzz;
! final if(co==0) block reads chin1). alpn1 is loop-scalar only in fused.
#ifdef AMSS_BATCH_STENCIL
#ifdef AMSS_RHS_RAW_DIAG_LOPSIDED
  raw_diag_lopsided = .false.
  if (bs_enabled() .and. Symmetry <= 1) raw_diag_lopsided = .true.
#endif
#endif
#ifdef AMSS_RHS_SKIP_CHIN1_SCAN
  if (co == 0) then
    chin1 = chi + ONE
  end if
#else
  chin1 = chi + ONE
#endif
#ifdef AMSS_BATCH_STENCIL
#ifdef AMSS_RHS_RAW_DIAG_LOPSIDED
  if (.not. raw_diag_lopsided) then
    gxx = dxx + ONE
    gyy = dyy + ONE
    gzz = dzz + ONE
  end if
#else
  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
#endif
#else
  gxx = dxx + ONE
  gyy = dyy + ONE
  gzz = dzz + ONE
#endif

!~~~~~~> all stencil calls hoisted before the fused loop (original order,
!         same input->output mapping as the original)

#ifdef AMSS_RHS_POINTWISE
  if (.not. pointwise_mode) then
#endif
  call fderivs(ex,betax,betaxx,betaxy,betaxz,X,Y,Z,ANTI, SYM, SYM,Symmetry,Lev)
  call fderivs(ex,betay,betayx,betayy,betayz,X,Y,Z, SYM,ANTI, SYM,Symmetry,Lev)
  call fderivs(ex,betaz,betazx,betazy,betazz,X,Y,Z, SYM, SYM,ANTI,Symmetry,Lev)

  call fderivs(ex,chi,chix,chiy,chiz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

  call fderivs(ex,dxx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gxy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,gxz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,Lev)
  call fderivs(ex,dyy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,gyz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,Lev)
  call fderivs(ex,dzz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,Lev)

  call fderivs(ex,Lap,Lapx,Lapy,Lapz,X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)
  call fderivs(ex,trK,Kx,Ky,Kz,X,Y,Z,SYM,SYM,SYM,symmetry,Lev)

  call fdderivs(ex,betax,bxxx,bxyx,bxzx,byyx,byzx,bzzx,&
                X,Y,Z,ANTI,SYM, SYM ,Symmetry,Lev)
  call fdderivs(ex,betay,bxxy,bxyy,bxzy,byyy,byzy,bzzy,&
                X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
  call fdderivs(ex,betaz,bxxz,bxyz,bxzz,byyz,byzz,bzzz,&
                X,Y,Z,SYM ,SYM, ANTI,Symmetry,Lev)

  call fderivs(ex,Gamx,Gamxx,Gamxy,Gamxz,X,Y,Z,ANTI,SYM ,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamy,Gamyx,Gamyy,Gamyz,X,Y,Z,SYM ,ANTI,SYM ,Symmetry,Lev)
  call fderivs(ex,Gamz,Gamzx,Gamzy,Gamzz,X,Y,Z,SYM ,SYM ,ANTI,Symmetry,Lev)

  call fdderivs(ex,dxx,fxx_dxx,fxy_dxx,fxz_dxx,fyy_dxx,fyz_dxx,fzz_dxx,&
                X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
  call fdderivs(ex,dyy,fxx_dyy,fxy_dyy,fxz_dyy,fyy_dyy,fyz_dyy,fzz_dyy,&
                X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
  call fdderivs(ex,dzz,fxx_dzz,fxy_dzz,fxz_dzz,fyy_dzz,fyz_dzz,fzz_dzz,&
                X,Y,Z,SYM ,SYM ,SYM ,symmetry,Lev)
  call fdderivs(ex,gxy,fxx_gxy,fxy_gxy,fxz_gxy,fyy_gxy,fyz_gxy,fzz_gxy,&
                X,Y,Z,ANTI, ANTI,SYM ,symmetry,Lev)
  call fdderivs(ex,gxz,fxx_gxz,fxy_gxz,fxz_gxz,fyy_gxz,fyz_gxz,fzz_gxz,&
                X,Y,Z,ANTI ,SYM ,ANTI,symmetry,Lev)
  call fdderivs(ex,gyz,fxx_gyz,fxy_gyz,fxz_gyz,fyy_gyz,fyz_gyz,fzz_gyz,&
                X,Y,Z,SYM ,ANTI ,ANTI,symmetry,Lev)

  call fdderivs(ex,chi,fxx_chi,fxy_chi,fxz_chi,fyy_chi,fyz_chi,fzz_chi,&
                X,Y,Z,SYM,SYM,SYM,Symmetry,Lev)

  call fdderivs(ex,Lap,fxx_lap,fxy_lap,fxz_lap,fyy_lap,fyz_lap,fzz_lap,&
                X,Y,Z,SYM,SYM,SYM,symmetry,Lev)
#ifdef AMSS_RHS_POINTWISE
  end if
#endif

!~~~~~~> fused single-point assembly loop (original lines 134-817)

#ifdef AMSS_RHS_TILED
  do jb = 1, ex(2), AMSS_RHS_TILE_J_VALUE
    jend = min(jb + AMSS_RHS_TILE_J_VALUE - 1, ex(2))
    do kb = 1, ex(3), AMSS_RHS_TILE_K_VALUE
      kend = min(kb + AMSS_RHS_TILE_K_VALUE - 1, ex(3))
      do k = kb, kend
        do j = jb, jend
#else
#ifdef AMSS_RHS_OMP_ASSEMBLY
! Parallelize the POINTWISE co==1 evolution hot path over k/j.  The if() clause
! keeps the loop serial when OMP=1 (no team overhead) and when pointwise_mode is
! false (co==0 constraint / octant), where the body takes the legacy full-array
! branches that read shared hoisted stencils.  Only the pw_*/_s per-point
! scalars are private; arrays and loop-invariant coefficients stay shared.
!$omp parallel do collapse(2) schedule(static) &
!$omp if(omp_get_max_threads() > 1 .and. pointwise_mode) &
!$omp private(alpn1_s, chin1_s, gxx_s, gyy_s, gzz_s, div_beta_s, &
!$omp& gupdet_s, gupxx_s, gupxy_s, gupxz_s, gupyy_s, gupyz_s, &
!$omp& gupzz_s, Gamxxx_s, Gamxxy_s, Gamxxz_s, Gamxyy_s, Gamxyz_s, &
!$omp& Gamxzz_s, Gamyxx_s, Gamyxy_s, Gamyxz_s, Gamyyy_s, Gamyyz_s, &
!$omp& Gamyzz_s, Gamzxx_s, Gamzxy_s, Gamzxz_s, Gamzyy_s, Gamzyz_s, &
!$omp& Gamzzz_s, Rxx_s, Rxy_s, Rxz_s, Ryy_s, Ryz_s, &
!$omp& Rzz_s, Gamx_rhs_s, Gamy_rhs_s, Gamz_rhs_s, trK_rhs_s, Axx_rhs_s, &
!$omp& Axy_rhs_s, Axz_rhs_s, Ayy_rhs_s, Ayz_rhs_s, Azz_rhs_s, fxx_s, &
!$omp& fxy_s, fxz_s, fyy_s, fyz_s, fzz_s, f_s, &
!$omp& S_s, Gamxa_s, Gamya_s, Gamza_s, gxxx_t, gxyx_t, &
!$omp& gxzx_t, gyyx_t, gyzx_t, gzzx_t, gxxy_t, gxyy_t, &
!$omp& gxzy_t, gyyy_t, gyzy_t, gzzy_t, gxxz_t, gxyz_t, &
!$omp& gxzz_t, gyyz_t, gyzz_t, gzzz_t, gxxx2_s, gxxy2_s, &
!$omp& gxxz2_s, i, j, k, pw_order, pw_bxxx, &
!$omp& pw_bxyx, pw_bxzx, pw_byyx, pw_byzx, pw_bzzx, pw_bxxy, &
!$omp& pw_bxyy, pw_bxzy, pw_byyy, pw_byzy, pw_bzzy, pw_bxxz, &
!$omp& pw_bxyz, pw_bxzz, pw_byyz, pw_byzz, pw_bzzz, pw_fxx_dxx, &
!$omp& pw_fxy_dxx, pw_fxz_dxx, pw_fyy_dxx, pw_fyz_dxx, pw_fzz_dxx, pw_fxx_dyy, &
!$omp& pw_fxy_dyy, pw_fxz_dyy, pw_fyy_dyy, pw_fyz_dyy, pw_fzz_dyy, pw_fxx_dzz, &
!$omp& pw_fxy_dzz, pw_fxz_dzz, pw_fyy_dzz, pw_fyz_dzz, pw_fzz_dzz, pw_fxx_gxy, &
!$omp& pw_fxy_gxy, pw_fxz_gxy, pw_fyy_gxy, pw_fyz_gxy, pw_fzz_gxy, pw_fxx_gxz, &
!$omp& pw_fxy_gxz, pw_fxz_gxz, pw_fyy_gxz, pw_fyz_gxz, pw_fzz_gxz, pw_fxx_gyz, &
!$omp& pw_fxy_gyz, pw_fxz_gyz, pw_fyy_gyz, pw_fyz_gyz, pw_fzz_gyz, pw_fxx_chi, &
!$omp& pw_fxy_chi, pw_fxz_chi, pw_fyy_chi, pw_fyz_chi, pw_fzz_chi, pw_fxx_lap, &
!$omp& pw_fxy_lap, pw_fxz_lap, pw_fyy_lap, pw_fyz_lap, pw_fzz_lap, pw_betaxx, &
!$omp& pw_betaxy, pw_betaxz, pw_betayx, pw_betayy, pw_betayz, pw_betazx, &
!$omp& pw_betazy, pw_betazz, pw_chix, pw_chiy, pw_chiz, pw_gxxx, &
!$omp& pw_gxxy, pw_gxxz, pw_gxyx, pw_gxyy, pw_gxyz, pw_gxzx, &
!$omp& pw_gxzy, pw_gxzz, pw_gyyx, pw_gyyy, pw_gyyz, pw_gyzx, &
!$omp& pw_gyzy, pw_gyzz, pw_gzzx, pw_gzzy, pw_gzzz, pw_lapx, &
!$omp& pw_lapy, pw_lapz, pw_kx, pw_ky, pw_kz, pw_gamxx, &
!$omp& pw_gamxy, pw_gamxz, pw_gamyx, pw_gamyy, pw_gamyz, pw_gamzx, &
!$omp& pw_gamzy, pw_gamzz)
#endif
  do k = 1, ex(3)
    do j = 1, ex(2)
#endif
#ifdef AMSS_RHS_BULK_SIMD
!$omp simd
#endif
      do i = 1, ex(1)

#ifdef AMSS_RHS_POINTWISE
        if (pointwise_mode) then
          ! Only the first two reflected z planes need the general helper.
          ! The bulk (k>=3) path below uses direct array references after
          ! preprocessing, so no pointwise derivative call remains in the
          ! cache-hot volume.
          if (pw_kmin < 1 .and. k < 3) then
            call point_d1(ex,X,Y,Z,betax,i,j,k,(/ANTI,SYM,SYM/),Symmetry,pw_betaxx,pw_betaxy,pw_betaxz)
            call point_d1(ex,X,Y,Z,betay,i,j,k,(/SYM,ANTI,SYM/),Symmetry,pw_betayx,pw_betayy,pw_betayz)
            call point_d1(ex,X,Y,Z,betaz,i,j,k,(/SYM,SYM,ANTI/),Symmetry,pw_betazx,pw_betazy,pw_betazz)
            call point_d1(ex,X,Y,Z,chi,i,j,k,(/SYM,SYM,SYM/),Symmetry,pw_chix,pw_chiy,pw_chiz)
            call point_d1(ex,X,Y,Z,dxx,i,j,k,(/SYM,SYM,SYM/),Symmetry,pw_gxxx,pw_gxxy,pw_gxxz)
            call point_d1(ex,X,Y,Z,gxy,i,j,k,(/ANTI,ANTI,SYM/),Symmetry,pw_gxyx,pw_gxyy,pw_gxyz)
            call point_d1(ex,X,Y,Z,gxz,i,j,k,(/ANTI,SYM,ANTI/),Symmetry,pw_gxzx,pw_gxzy,pw_gxzz)
            call point_d1(ex,X,Y,Z,dyy,i,j,k,(/SYM,SYM,SYM/),Symmetry,pw_gyyx,pw_gyyy,pw_gyyz)
            call point_d1(ex,X,Y,Z,gyz,i,j,k,(/SYM,ANTI,ANTI/),Symmetry,pw_gyzx,pw_gyzy,pw_gyzz)
            call point_d1(ex,X,Y,Z,dzz,i,j,k,(/SYM,SYM,SYM/),Symmetry,pw_gzzx,pw_gzzy,pw_gzzz)
            call point_d1(ex,X,Y,Z,Lap,i,j,k,(/SYM,SYM,SYM/),Symmetry,pw_lapx,pw_lapy,pw_lapz)
            call point_d1(ex,X,Y,Z,trK,i,j,k,(/SYM,SYM,SYM/),Symmetry,pw_kx,pw_ky,pw_kz)
            call point_d1(ex,X,Y,Z,Gamx,i,j,k,(/ANTI,SYM,SYM/),Symmetry,pw_gamxx,pw_gamxy,pw_gamxz)
            call point_d1(ex,X,Y,Z,Gamy,i,j,k,(/SYM,ANTI,SYM/),Symmetry,pw_gamyx,pw_gamyy,pw_gamyz)
            call point_d1(ex,X,Y,Z,Gamz,i,j,k,(/SYM,SYM,ANTI/),Symmetry,pw_gamzx,pw_gamzy,pw_gamzz)
            call point_d2(ex,X,Y,Z,betax,i,j,k,(/ANTI,SYM,SYM/),Symmetry,pw_bxxx,pw_bxyx,pw_bxzx,pw_byyx,pw_byzx,pw_bzzx)
            call point_d2(ex,X,Y,Z,betay,i,j,k,(/SYM,ANTI,SYM/),Symmetry,pw_bxxy,pw_bxyy,pw_bxzy,pw_byyy,pw_byzy,pw_bzzy)
            call point_d2(ex,X,Y,Z,betaz,i,j,k,(/SYM,SYM,ANTI/),Symmetry,pw_bxxz,pw_bxyz,pw_bxzz,pw_byyz,pw_byzz,pw_bzzz)
            call point_d2(ex,X,Y,Z,dxx,i,j,k,(/SYM,SYM,SYM/),Symmetry,pw_fxx_dxx,pw_fxy_dxx,pw_fxz_dxx,pw_fyy_dxx,pw_fyz_dxx,pw_fzz_dxx)
            call point_d2(ex,X,Y,Z,dyy,i,j,k,(/SYM,SYM,SYM/),Symmetry,pw_fxx_dyy,pw_fxy_dyy,pw_fxz_dyy,pw_fyy_dyy,pw_fyz_dyy,pw_fzz_dyy)
            call point_d2(ex,X,Y,Z,dzz,i,j,k,(/SYM,SYM,SYM/),Symmetry,pw_fxx_dzz,pw_fxy_dzz,pw_fxz_dzz,pw_fyy_dzz,pw_fyz_dzz,pw_fzz_dzz)
            call point_d2(ex,X,Y,Z,gxy,i,j,k,(/ANTI,ANTI,SYM/),Symmetry,pw_fxx_gxy,pw_fxy_gxy,pw_fxz_gxy,pw_fyy_gxy,pw_fyz_gxy,pw_fzz_gxy)
            call point_d2(ex,X,Y,Z,gxz,i,j,k,(/ANTI,SYM,ANTI/),Symmetry,pw_fxx_gxz,pw_fxy_gxz,pw_fxz_gxz,pw_fyy_gxz,pw_fyz_gxz,pw_fzz_gxz)
            call point_d2(ex,X,Y,Z,gyz,i,j,k,(/SYM,ANTI,ANTI/),Symmetry,pw_fxx_gyz,pw_fxy_gyz,pw_fxz_gyz,pw_fyy_gyz,pw_fyz_gyz,pw_fzz_gyz)
            call point_d2(ex,X,Y,Z,chi,i,j,k,(/SYM,SYM,SYM/),Symmetry,pw_fxx_chi,pw_fxy_chi,pw_fxz_chi,pw_fyy_chi,pw_fyz_chi,pw_fzz_chi)
            call point_d2(ex,X,Y,Z,Lap,i,j,k,(/SYM,SYM,SYM/),Symmetry,pw_fxx_lap,pw_fxy_lap,pw_fxz_lap,pw_fyy_lap,pw_fyz_lap,pw_fzz_lap)
          else
            ! The stencil order depends only on this point's location.  The
            ! field-independent bounds test is deliberately done once here,
            ! rather than once in each of the 15 first- and 11 second-
            ! derivative expansions below.
            pw_order = 0
            if (i+2 <= ex(1) .and. i-2 >= pw_imin .and. &
                j+2 <= ex(2) .and. j-2 >= pw_jmin .and. &
                k+2 <= ex(3) .and. k-2 >= pw_kmin) then
              pw_order = 4
            else if (i+1 <= ex(1) .and. i-1 >= pw_imin .and. &
                     j+1 <= ex(2) .and. j-1 >= pw_jmin .and. &
                     k+1 <= ex(3) .and. k-1 >= pw_kmin) then
              pw_order = 2
            end if
#define PW_VALUE(ii,jj,kk) PW_FIELD(ii,jj,kk)
#define PW_FIELD betax
#define PW_S1 ANTI
#define PW_S2 SYM
#define PW_S3 SYM
#define PW_FX pw_betaxx
#define PW_FY pw_betaxy
#define PW_FZ pw_betaxz
#include "point_d1_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FX
#undef PW_FY
#undef PW_FZ

#define PW_FIELD betay
#define PW_S1 SYM
#define PW_S2 ANTI
#define PW_S3 SYM
#define PW_FX pw_betayx
#define PW_FY pw_betayy
#define PW_FZ pw_betayz
#include "point_d1_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FX
#undef PW_FY
#undef PW_FZ

#define PW_FIELD betaz
#define PW_S1 SYM
#define PW_S2 SYM
#define PW_S3 ANTI
#define PW_FX pw_betazx
#define PW_FY pw_betazy
#define PW_FZ pw_betazz
#include "point_d1_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FX
#undef PW_FY
#undef PW_FZ

#define PW_FIELD chi
#define PW_S1 SYM
#define PW_S2 SYM
#define PW_S3 SYM
#define PW_FX pw_chix
#define PW_FY pw_chiy
#define PW_FZ pw_chiz
#include "point_d1_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FX
#undef PW_FY
#undef PW_FZ

#define PW_FIELD dxx
#define PW_S1 SYM
#define PW_S2 SYM
#define PW_S3 SYM
#define PW_FX pw_gxxx
#define PW_FY pw_gxxy
#define PW_FZ pw_gxxz
#include "point_d1_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FX
#undef PW_FY
#undef PW_FZ

#define PW_FIELD gxy
#define PW_S1 ANTI
#define PW_S2 ANTI
#define PW_S3 SYM
#define PW_FX pw_gxyx
#define PW_FY pw_gxyy
#define PW_FZ pw_gxyz
#include "point_d1_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FX
#undef PW_FY
#undef PW_FZ

#define PW_FIELD gxz
#define PW_S1 ANTI
#define PW_S2 SYM
#define PW_S3 ANTI
#define PW_FX pw_gxzx
#define PW_FY pw_gxzy
#define PW_FZ pw_gxzz
#include "point_d1_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FX
#undef PW_FY
#undef PW_FZ

#define PW_FIELD dyy
#define PW_S1 SYM
#define PW_S2 SYM
#define PW_S3 SYM
#define PW_FX pw_gyyx
#define PW_FY pw_gyyy
#define PW_FZ pw_gyyz
#include "point_d1_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FX
#undef PW_FY
#undef PW_FZ

#define PW_FIELD gyz
#define PW_S1 SYM
#define PW_S2 ANTI
#define PW_S3 ANTI
#define PW_FX pw_gyzx
#define PW_FY pw_gyzy
#define PW_FZ pw_gyzz
#include "point_d1_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FX
#undef PW_FY
#undef PW_FZ

#define PW_FIELD dzz
#define PW_S1 SYM
#define PW_S2 SYM
#define PW_S3 SYM
#define PW_FX pw_gzzx
#define PW_FY pw_gzzy
#define PW_FZ pw_gzzz
#include "point_d1_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FX
#undef PW_FY
#undef PW_FZ

#define PW_FIELD Lap
#define PW_S1 SYM
#define PW_S2 SYM
#define PW_S3 SYM
#define PW_FX pw_lapx
#define PW_FY pw_lapy
#define PW_FZ pw_lapz
#include "point_d1_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FX
#undef PW_FY
#undef PW_FZ

#define PW_FIELD trK
#define PW_S1 SYM
#define PW_S2 SYM
#define PW_S3 SYM
#define PW_FX pw_kx
#define PW_FY pw_ky
#define PW_FZ pw_kz
#include "point_d1_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FX
#undef PW_FY
#undef PW_FZ

#define PW_FIELD Gamx
#define PW_S1 ANTI
#define PW_S2 SYM
#define PW_S3 SYM
#define PW_FX pw_gamxx
#define PW_FY pw_gamxy
#define PW_FZ pw_gamxz
#include "point_d1_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FX
#undef PW_FY
#undef PW_FZ

#define PW_FIELD Gamy
#define PW_S1 SYM
#define PW_S2 ANTI
#define PW_S3 SYM
#define PW_FX pw_gamyx
#define PW_FY pw_gamyy
#define PW_FZ pw_gamyz
#include "point_d1_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FX
#undef PW_FY
#undef PW_FZ

#define PW_FIELD Gamz
#define PW_S1 SYM
#define PW_S2 SYM
#define PW_S3 ANTI
#define PW_FX pw_gamzx
#define PW_FY pw_gamzy
#define PW_FZ pw_gamzz
#include "point_d1_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FX
#undef PW_FY
#undef PW_FZ

#define PW_FIELD betax
#define PW_S1 ANTI
#define PW_S2 SYM
#define PW_S3 SYM
#define PW_FXX pw_bxxx
#define PW_FXY pw_bxyx
#define PW_FXZ pw_bxzx
#define PW_FYY pw_byyx
#define PW_FYZ pw_byzx
#define PW_FZZ pw_bzzx
#include "point_d2_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FXX
#undef PW_FXY
#undef PW_FXZ
#undef PW_FYY
#undef PW_FYZ
#undef PW_FZZ

#define PW_FIELD betay
#define PW_S1 SYM
#define PW_S2 ANTI
#define PW_S3 SYM
#define PW_FXX pw_bxxy
#define PW_FXY pw_bxyy
#define PW_FXZ pw_bxzy
#define PW_FYY pw_byyy
#define PW_FYZ pw_byzy
#define PW_FZZ pw_bzzy
#include "point_d2_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FXX
#undef PW_FXY
#undef PW_FXZ
#undef PW_FYY
#undef PW_FYZ
#undef PW_FZZ

#define PW_FIELD betaz
#define PW_S1 SYM
#define PW_S2 SYM
#define PW_S3 ANTI
#define PW_FXX pw_bxxz
#define PW_FXY pw_bxyz
#define PW_FXZ pw_bxzz
#define PW_FYY pw_byyz
#define PW_FYZ pw_byzz
#define PW_FZZ pw_bzzz
#include "point_d2_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FXX
#undef PW_FXY
#undef PW_FXZ
#undef PW_FYY
#undef PW_FYZ
#undef PW_FZZ

#define PW_FIELD dxx
#define PW_S1 SYM
#define PW_S2 SYM
#define PW_S3 SYM
#define PW_FXX pw_fxx_dxx
#define PW_FXY pw_fxy_dxx
#define PW_FXZ pw_fxz_dxx
#define PW_FYY pw_fyy_dxx
#define PW_FYZ pw_fyz_dxx
#define PW_FZZ pw_fzz_dxx
#include "point_d2_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FXX
#undef PW_FXY
#undef PW_FXZ
#undef PW_FYY
#undef PW_FYZ
#undef PW_FZZ

#define PW_FIELD dyy
#define PW_S1 SYM
#define PW_S2 SYM
#define PW_S3 SYM
#define PW_FXX pw_fxx_dyy
#define PW_FXY pw_fxy_dyy
#define PW_FXZ pw_fxz_dyy
#define PW_FYY pw_fyy_dyy
#define PW_FYZ pw_fyz_dyy
#define PW_FZZ pw_fzz_dyy
#include "point_d2_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FXX
#undef PW_FXY
#undef PW_FXZ
#undef PW_FYY
#undef PW_FYZ
#undef PW_FZZ

#define PW_FIELD dzz
#define PW_S1 SYM
#define PW_S2 SYM
#define PW_S3 SYM
#define PW_FXX pw_fxx_dzz
#define PW_FXY pw_fxy_dzz
#define PW_FXZ pw_fxz_dzz
#define PW_FYY pw_fyy_dzz
#define PW_FYZ pw_fyz_dzz
#define PW_FZZ pw_fzz_dzz
#include "point_d2_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FXX
#undef PW_FXY
#undef PW_FXZ
#undef PW_FYY
#undef PW_FYZ
#undef PW_FZZ

#define PW_FIELD gxy
#define PW_S1 ANTI
#define PW_S2 ANTI
#define PW_S3 SYM
#define PW_FXX pw_fxx_gxy
#define PW_FXY pw_fxy_gxy
#define PW_FXZ pw_fxz_gxy
#define PW_FYY pw_fyy_gxy
#define PW_FYZ pw_fyz_gxy
#define PW_FZZ pw_fzz_gxy
#include "point_d2_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FXX
#undef PW_FXY
#undef PW_FXZ
#undef PW_FYY
#undef PW_FYZ
#undef PW_FZZ

#define PW_FIELD gxz
#define PW_S1 ANTI
#define PW_S2 SYM
#define PW_S3 ANTI
#define PW_FXX pw_fxx_gxz
#define PW_FXY pw_fxy_gxz
#define PW_FXZ pw_fxz_gxz
#define PW_FYY pw_fyy_gxz
#define PW_FYZ pw_fyz_gxz
#define PW_FZZ pw_fzz_gxz
#include "point_d2_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FXX
#undef PW_FXY
#undef PW_FXZ
#undef PW_FYY
#undef PW_FYZ
#undef PW_FZZ

#define PW_FIELD gyz
#define PW_S1 SYM
#define PW_S2 ANTI
#define PW_S3 ANTI
#define PW_FXX pw_fxx_gyz
#define PW_FXY pw_fxy_gyz
#define PW_FXZ pw_fxz_gyz
#define PW_FYY pw_fyy_gyz
#define PW_FYZ pw_fyz_gyz
#define PW_FZZ pw_fzz_gyz
#include "point_d2_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FXX
#undef PW_FXY
#undef PW_FXZ
#undef PW_FYY
#undef PW_FYZ
#undef PW_FZZ

#define PW_FIELD chi
#define PW_S1 SYM
#define PW_S2 SYM
#define PW_S3 SYM
#define PW_FXX pw_fxx_chi
#define PW_FXY pw_fxy_chi
#define PW_FXZ pw_fxz_chi
#define PW_FYY pw_fyy_chi
#define PW_FYZ pw_fyz_chi
#define PW_FZZ pw_fzz_chi
#include "point_d2_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FXX
#undef PW_FXY
#undef PW_FXZ
#undef PW_FYY
#undef PW_FYZ
#undef PW_FZZ

#define PW_FIELD Lap
#define PW_S1 SYM
#define PW_S2 SYM
#define PW_S3 SYM
#define PW_FXX pw_fxx_lap
#define PW_FXY pw_fxy_lap
#define PW_FXZ pw_fxz_lap
#define PW_FYY pw_fyy_lap
#define PW_FYZ pw_fyz_lap
#define PW_FZZ pw_fzz_lap
#include "point_d2_inline.fh"
#undef PW_FIELD
#undef PW_S1
#undef PW_S2
#undef PW_S3
#undef PW_FXX
#undef PW_FXY
#undef PW_FXZ
#undef PW_FYY
#undef PW_FYZ
#undef PW_FZZ
           ! End of the preprocessor-expanded bulk stencil path.
#undef PW_VALUE
          end if
        end if
#endif

#ifdef AMSS_RHS_POINTWISE
        ! The legacy branch still needs the precomputed arrays, but the
        ! POINTWISE branch already populated these same values as scalars.
        if (.not. pointwise_mode) then
          pw_betaxx = betaxx(i,j,k); pw_betaxy = betaxy(i,j,k); pw_betaxz = betaxz(i,j,k)
          pw_betayx = betayx(i,j,k); pw_betayy = betayy(i,j,k); pw_betayz = betayz(i,j,k)
          pw_betazx = betazx(i,j,k); pw_betazy = betazy(i,j,k); pw_betazz = betazz(i,j,k)
          pw_chix = chix(i,j,k); pw_chiy = chiy(i,j,k); pw_chiz = chiz(i,j,k)
          pw_gxxx = gxxx(i,j,k); pw_gxxy = gxxy(i,j,k); pw_gxxz = gxxz(i,j,k)
          pw_gxyx = gxyx(i,j,k); pw_gxyy = gxyy(i,j,k); pw_gxyz = gxyz(i,j,k)
          pw_gxzx = gxzx(i,j,k); pw_gxzy = gxzy(i,j,k); pw_gxzz = gxzz(i,j,k)
          pw_gyyx = gyyx(i,j,k); pw_gyyy = gyyy(i,j,k); pw_gyyz = gyyz(i,j,k)
          pw_gyzx = gyzx(i,j,k); pw_gyzy = gyzy(i,j,k); pw_gyzz = gyzz(i,j,k)
          pw_gzzx = gzzx(i,j,k); pw_gzzy = gzzy(i,j,k); pw_gzzz = gzzz(i,j,k)
          pw_lapx = Lapx(i,j,k); pw_lapy = Lapy(i,j,k); pw_lapz = Lapz(i,j,k)
          pw_kx = Kx(i,j,k); pw_ky = Ky(i,j,k); pw_kz = Kz(i,j,k)
          pw_gamxx = Gamxx(i,j,k); pw_gamxy = Gamxy(i,j,k); pw_gamxz = Gamxz(i,j,k)
          pw_gamyx = Gamyx(i,j,k); pw_gamyy = Gamyy(i,j,k); pw_gamyz = Gamyz(i,j,k)
          pw_gamzx = Gamzx(i,j,k); pw_gamzy = Gamzy(i,j,k); pw_gamzz = Gamzz(i,j,k)
        end if
#define betaxx(ii,jj,kk) pw_betaxx
#define betaxy(ii,jj,kk) pw_betaxy
#define betaxz(ii,jj,kk) pw_betaxz
#define betayx(ii,jj,kk) pw_betayx
#define betayy(ii,jj,kk) pw_betayy
#define betayz(ii,jj,kk) pw_betayz
#define betazx(ii,jj,kk) pw_betazx
#define betazy(ii,jj,kk) pw_betazy
#define betazz(ii,jj,kk) pw_betazz
#define chix(ii,jj,kk) pw_chix
#define chiy(ii,jj,kk) pw_chiy
#define chiz(ii,jj,kk) pw_chiz
#define gxxx(ii,jj,kk) pw_gxxx
#define gxxy(ii,jj,kk) pw_gxxy
#define gxxz(ii,jj,kk) pw_gxxz
#define gxyx(ii,jj,kk) pw_gxyx
#define gxyy(ii,jj,kk) pw_gxyy
#define gxyz(ii,jj,kk) pw_gxyz
#define gxzx(ii,jj,kk) pw_gxzx
#define gxzy(ii,jj,kk) pw_gxzy
#define gxzz(ii,jj,kk) pw_gxzz
#define gyyx(ii,jj,kk) pw_gyyx
#define gyyy(ii,jj,kk) pw_gyyy
#define gyyz(ii,jj,kk) pw_gyyz
#define gyzx(ii,jj,kk) pw_gyzx
#define gyzy(ii,jj,kk) pw_gyzy
#define gyzz(ii,jj,kk) pw_gyzz
#define gzzx(ii,jj,kk) pw_gzzx
#define gzzy(ii,jj,kk) pw_gzzy
#define gzzz(ii,jj,kk) pw_gzzz
#define Lapx(ii,jj,kk) pw_lapx
#define Lapy(ii,jj,kk) pw_lapy
#define Lapz(ii,jj,kk) pw_lapz
#define Kx(ii,jj,kk) pw_kx
#define Ky(ii,jj,kk) pw_ky
#define Kz(ii,jj,kk) pw_kz
#define Gamxx(ii,jj,kk) pw_gamxx
#define Gamxy(ii,jj,kk) pw_gamxy
#define Gamxz(ii,jj,kk) pw_gamxz
#define Gamyx(ii,jj,kk) pw_gamyx
#define Gamyy(ii,jj,kk) pw_gamyy
#define Gamyz(ii,jj,kk) pw_gamyz
#define Gamzx(ii,jj,kk) pw_gamzx
#define Gamzy(ii,jj,kk) pw_gamzy
#define Gamzz(ii,jj,kk) pw_gamzz
#endif

        alpn1_s = Lap(i,j,k) + ONE
        chin1_s = chi(i,j,k) + ONE
        gxx_s = dxx(i,j,k) + ONE
        gyy_s = dyy(i,j,k) + ONE
        gzz_s = dzz(i,j,k) + ONE
        div_beta_s = betaxx(i,j,k) + betayy(i,j,k) + betazz(i,j,k)

        chi_rhs(i,j,k) = F2o3 * chin1_s * ( alpn1_s * trK(i,j,k) - div_beta_s ) !rhs for chi

#ifdef AMSS_RHS_REORDER
        ! Phase 1 (lifetime scheduling): shallow gauge RHS have no downstream
        ! consumer inside the loop.  Store them now so alpn1_s/dtSf reads die
        ! immediately.  trK(i,j,k) and dtSfx/y/z are loop-invariant inputs.
        ! dtSf_rhs (Gamma-driver) is NOT here -- it needs Gamx_rhs_s and is
        ! stored right after the Gamma_rhs store-back below.
        Lap_rhs(i,j,k) = -TWO*alpn1_s*trK(i,j,k)
        betax_rhs(i,j,k) = FF*dtSfx(i,j,k)
        betay_rhs(i,j,k) = FF*dtSfy(i,j,k)
        betaz_rhs(i,j,k) = FF*dtSfz(i,j,k)
#endif

        gxx_rhs(i,j,k) = - TWO * alpn1_s * Axx(i,j,k) - F2o3 * gxx_s * div_beta_s + &
              TWO *(  gxx_s * betaxx(i,j,k) +   gxy(i,j,k) * betayx(i,j,k) +   gxz(i,j,k) * betazx(i,j,k))

        gyy_rhs(i,j,k) = - TWO * alpn1_s * Ayy(i,j,k) - F2o3 * gyy_s * div_beta_s + &
              TWO *(  gxy(i,j,k) * betaxy(i,j,k) +   gyy_s * betayy(i,j,k) +   gyz(i,j,k) * betazy(i,j,k))

        gzz_rhs(i,j,k) = - TWO * alpn1_s * Azz(i,j,k) - F2o3 * gzz_s * div_beta_s + &
              TWO *(  gxz(i,j,k) * betaxz(i,j,k) +   gyz(i,j,k) * betayz(i,j,k) +   gzz_s * betazz(i,j,k))

        gxy_rhs(i,j,k) = - TWO * alpn1_s * Axy(i,j,k) + F1o3 * gxy(i,j,k) * div_beta_s + &
                      gxx_s * betaxy(i,j,k)                  +   gxz(i,j,k) * betazy(i,j,k) + &
                                       gyy_s * betayx(i,j,k) +   gyz(i,j,k) * betazx(i,j,k)   &
                                                    -   gxy(i,j,k) * betazz(i,j,k)

        gyz_rhs(i,j,k) = - TWO * alpn1_s * Ayz(i,j,k) + F1o3 * gyz(i,j,k) * div_beta_s + &
                      gxy(i,j,k) * betaxz(i,j,k) +   gyy_s * betayz(i,j,k)                  + &
                      gxz(i,j,k) * betaxy(i,j,k)                  +   gzz_s * betazy(i,j,k)   &
                                                    -   gyz(i,j,k) * betaxx(i,j,k)

        gxz_rhs(i,j,k) = - TWO * alpn1_s * Axz(i,j,k) + F1o3 * gxz(i,j,k) * div_beta_s + &
                      gxx_s * betaxz(i,j,k) +   gxy(i,j,k) * betayz(i,j,k)                  + &
                                       gyz(i,j,k) * betayx(i,j,k) +   gzz_s * betazx(i,j,k)   &
                                                    -   gxz(i,j,k) * betayy(i,j,k)     !rhs for gij

! invert tilted metric
        gupdet_s =  gxx_s * gyy_s * gzz_s + gxy(i,j,k) * gyz(i,j,k) * gxz(i,j,k) + gxz(i,j,k) * gxy(i,j,k) * gyz(i,j,k) - &
             gxz(i,j,k) * gyy_s * gxz(i,j,k) - gxy(i,j,k) * gxy(i,j,k) * gzz_s - gxx_s * gyz(i,j,k) * gyz(i,j,k)
        gupxx_s =   ( gyy_s * gzz_s - gyz(i,j,k) * gyz(i,j,k) ) / gupdet_s
        gupxy_s = - ( gxy(i,j,k) * gzz_s - gyz(i,j,k) * gxz(i,j,k) ) / gupdet_s
        gupxz_s =   ( gxy(i,j,k) * gyz(i,j,k) - gyy_s * gxz(i,j,k) ) / gupdet_s
        gupyy_s =   ( gxx_s * gzz_s - gxz(i,j,k) * gxz(i,j,k) ) / gupdet_s
        gupyz_s = - ( gxx_s * gyz(i,j,k) - gxy(i,j,k) * gxz(i,j,k) ) / gupdet_s
        gupzz_s =   ( gxx_s * gyy_s - gxy(i,j,k) * gxy(i,j,k) ) / gupdet_s

        if(co == 0)then
! Gam^i_Res = Gam^i + gup^ij_,j  (reads 1st generation gxxx..gzzz)
          Gmx_Res(i,j,k) = Gamx(i,j,k) - (gupxx_s*(gupxx_s*gxxx(i,j,k)+gupxy_s*gxyx(i,j,k)+gupxz_s*gxzx(i,j,k))&
                     +gupxy_s*(gupxx_s*gxyx(i,j,k)+gupxy_s*gyyx(i,j,k)+gupxz_s*gyzx(i,j,k))&
                     +gupxz_s*(gupxx_s*gxzx(i,j,k)+gupxy_s*gyzx(i,j,k)+gupxz_s*gzzx(i,j,k))&
                     +gupxx_s*(gupxy_s*gxxy(i,j,k)+gupyy_s*gxyy(i,j,k)+gupyz_s*gxzy(i,j,k))&
                     +gupxy_s*(gupxy_s*gxyy(i,j,k)+gupyy_s*gyyy(i,j,k)+gupyz_s*gyzy(i,j,k))&
                     +gupxz_s*(gupxy_s*gxzy(i,j,k)+gupyy_s*gyzy(i,j,k)+gupyz_s*gzzy(i,j,k))&
                     +gupxx_s*(gupxz_s*gxxz(i,j,k)+gupyz_s*gxyz(i,j,k)+gupzz_s*gxzz(i,j,k))&
                     +gupxy_s*(gupxz_s*gxyz(i,j,k)+gupyz_s*gyyz(i,j,k)+gupzz_s*gyzz(i,j,k))&
                     +gupxz_s*(gupxz_s*gxzz(i,j,k)+gupyz_s*gyzz(i,j,k)+gupzz_s*gzzz(i,j,k)))
          Gmy_Res(i,j,k) = Gamy(i,j,k) - (gupxx_s*(gupxy_s*gxxx(i,j,k)+gupyy_s*gxyx(i,j,k)+gupyz_s*gxzx(i,j,k))&
                     +gupxy_s*(gupxy_s*gxyx(i,j,k)+gupyy_s*gyyx(i,j,k)+gupyz_s*gyzx(i,j,k))&
                     +gupxz_s*(gupxy_s*gxzx(i,j,k)+gupyy_s*gyzx(i,j,k)+gupyz_s*gzzx(i,j,k))&
                     +gupxy_s*(gupxy_s*gxxy(i,j,k)+gupyy_s*gxyy(i,j,k)+gupyz_s*gxzy(i,j,k))&
                     +gupyy_s*(gupxy_s*gxyy(i,j,k)+gupyy_s*gyyy(i,j,k)+gupyz_s*gyzy(i,j,k))&
                     +gupyz_s*(gupxy_s*gxzy(i,j,k)+gupyy_s*gyzy(i,j,k)+gupyz_s*gzzy(i,j,k))&
                     +gupxy_s*(gupxz_s*gxxz(i,j,k)+gupyz_s*gxyz(i,j,k)+gupzz_s*gxzz(i,j,k))&
                     +gupyy_s*(gupxz_s*gxyz(i,j,k)+gupyz_s*gyyz(i,j,k)+gupzz_s*gyzz(i,j,k))&
                     +gupyz_s*(gupxz_s*gxzz(i,j,k)+gupyz_s*gyzz(i,j,k)+gupzz_s*gzzz(i,j,k)))
          Gmz_Res(i,j,k) = Gamz(i,j,k) - (gupxx_s*(gupxz_s*gxxx(i,j,k)+gupyz_s*gxyx(i,j,k)+gupzz_s*gxzx(i,j,k))&
                     +gupxy_s*(gupxz_s*gxyx(i,j,k)+gupyz_s*gyyx(i,j,k)+gupzz_s*gyzx(i,j,k))&
                     +gupxz_s*(gupxz_s*gxzx(i,j,k)+gupyz_s*gyzx(i,j,k)+gupzz_s*gzzx(i,j,k))&
                     +gupxy_s*(gupxz_s*gxxy(i,j,k)+gupyz_s*gxyy(i,j,k)+gupzz_s*gxzy(i,j,k))&
                     +gupyy_s*(gupxz_s*gxyy(i,j,k)+gupyz_s*gyyy(i,j,k)+gupzz_s*gyzy(i,j,k))&
                     +gupyz_s*(gupxz_s*gxzy(i,j,k)+gupyz_s*gyzy(i,j,k)+gupzz_s*gzzy(i,j,k))&
                     +gupxz_s*(gupxz_s*gxxz(i,j,k)+gupyz_s*gxyz(i,j,k)+gupzz_s*gxzz(i,j,k))&
                     +gupyz_s*(gupxz_s*gxyz(i,j,k)+gupyz_s*gyyz(i,j,k)+gupzz_s*gyzz(i,j,k))&
                     +gupzz_s*(gupxz_s*gxzz(i,j,k)+gupyz_s*gyzz(i,j,k)+gupzz_s*gzzz(i,j,k)))
        endif

! second kind of connection (conformal, reads 1st generation gxxx..gzzz)
        Gamxxx_s =HALF*( gupxx_s*gxxx(i,j,k) + gupxy_s*(TWO*gxyx(i,j,k) - gxxy(i,j,k) ) + gupxz_s*(TWO*gxzx(i,j,k) - gxxz(i,j,k) ))
        Gamyxx_s =HALF*( gupxy_s*gxxx(i,j,k) + gupyy_s*(TWO*gxyx(i,j,k) - gxxy(i,j,k) ) + gupyz_s*(TWO*gxzx(i,j,k) - gxxz(i,j,k) ))
        Gamzxx_s =HALF*( gupxz_s*gxxx(i,j,k) + gupyz_s*(TWO*gxyx(i,j,k) - gxxy(i,j,k) ) + gupzz_s*(TWO*gxzx(i,j,k) - gxxz(i,j,k) ))

        Gamxyy_s =HALF*( gupxx_s*(TWO*gxyy(i,j,k) - gyyx(i,j,k) ) + gupxy_s*gyyy(i,j,k) + gupxz_s*(TWO*gyzy(i,j,k) - gyyz(i,j,k) ))
        Gamyyy_s =HALF*( gupxy_s*(TWO*gxyy(i,j,k) - gyyx(i,j,k) ) + gupyy_s*gyyy(i,j,k) + gupyz_s*(TWO*gyzy(i,j,k) - gyyz(i,j,k) ))
        Gamzyy_s =HALF*( gupxz_s*(TWO*gxyy(i,j,k) - gyyx(i,j,k) ) + gupyz_s*gyyy(i,j,k) + gupzz_s*(TWO*gyzy(i,j,k) - gyyz(i,j,k) ))

        Gamxzz_s =HALF*( gupxx_s*(TWO*gxzz(i,j,k) - gzzx(i,j,k) ) + gupxy_s*(TWO*gyzz(i,j,k) - gzzy(i,j,k) ) + gupxz_s*gzzz(i,j,k))
        Gamyzz_s =HALF*( gupxy_s*(TWO*gxzz(i,j,k) - gzzx(i,j,k) ) + gupyy_s*(TWO*gyzz(i,j,k) - gzzy(i,j,k) ) + gupyz_s*gzzz(i,j,k))
        Gamzzz_s =HALF*( gupxz_s*(TWO*gxzz(i,j,k) - gzzx(i,j,k) ) + gupyz_s*(TWO*gyzz(i,j,k) - gzzy(i,j,k) ) + gupzz_s*gzzz(i,j,k))

        Gamxxy_s =HALF*( gupxx_s*gxxy(i,j,k) + gupxy_s*gyyx(i,j,k) + gupxz_s*( gxzy(i,j,k) + gyzx(i,j,k) - gxyz(i,j,k) ) )
        Gamyxy_s =HALF*( gupxy_s*gxxy(i,j,k) + gupyy_s*gyyx(i,j,k) + gupyz_s*( gxzy(i,j,k) + gyzx(i,j,k) - gxyz(i,j,k) ) )
        Gamzxy_s =HALF*( gupxz_s*gxxy(i,j,k) + gupyz_s*gyyx(i,j,k) + gupzz_s*( gxzy(i,j,k) + gyzx(i,j,k) - gxyz(i,j,k) ) )

        Gamxxz_s =HALF*( gupxx_s*gxxz(i,j,k) + gupxy_s*( gxyz(i,j,k) + gyzx(i,j,k) - gxzy(i,j,k) ) + gupxz_s*gzzx(i,j,k) )
        Gamyxz_s =HALF*( gupxy_s*gxxz(i,j,k) + gupyy_s*( gxyz(i,j,k) + gyzx(i,j,k) - gxzy(i,j,k) ) + gupyz_s*gzzx(i,j,k) )
        Gamzxz_s =HALF*( gupxz_s*gxxz(i,j,k) + gupyz_s*( gxyz(i,j,k) + gyzx(i,j,k) - gxzy(i,j,k) ) + gupzz_s*gzzx(i,j,k) )

        Gamxyz_s =HALF*( gupxx_s*( gxyz(i,j,k) + gxzy(i,j,k) - gyzx(i,j,k) ) + gupxy_s*gyyz(i,j,k) + gupxz_s*gzzy(i,j,k) )
        Gamyyz_s =HALF*( gupxy_s*( gxyz(i,j,k) + gxzy(i,j,k) - gyzx(i,j,k) ) + gupyy_s*gyyz(i,j,k) + gupyz_s*gzzy(i,j,k) )
        Gamzyz_s =HALF*( gupxz_s*( gxyz(i,j,k) + gxzy(i,j,k) - gyzx(i,j,k) ) + gupyz_s*gyyz(i,j,k) + gupzz_s*gzzy(i,j,k) )
! Raise indices of \tilde A_{ij} and store in R_ij

        Rxx_s =    gupxx_s * gupxx_s * Axx(i,j,k) + gupxy_s * gupxy_s * Ayy(i,j,k) + gupxz_s * gupxz_s * Azz(i,j,k) + &
            TWO*(gupxx_s * gupxy_s * Axy(i,j,k) + gupxx_s * gupxz_s * Axz(i,j,k) + gupxy_s * gupxz_s * Ayz(i,j,k))

        Ryy_s =    gupxy_s * gupxy_s * Axx(i,j,k) + gupyy_s * gupyy_s * Ayy(i,j,k) + gupyz_s * gupyz_s * Azz(i,j,k) + &
            TWO*(gupxy_s * gupyy_s * Axy(i,j,k) + gupxy_s * gupyz_s * Axz(i,j,k) + gupyy_s * gupyz_s * Ayz(i,j,k))

        Rzz_s =    gupxz_s * gupxz_s * Axx(i,j,k) + gupyz_s * gupyz_s * Ayy(i,j,k) + gupzz_s * gupzz_s * Azz(i,j,k) + &
            TWO*(gupxz_s * gupyz_s * Axy(i,j,k) + gupxz_s * gupzz_s * Axz(i,j,k) + gupyz_s * gupzz_s * Ayz(i,j,k))

        Rxy_s =    gupxx_s * gupxy_s * Axx(i,j,k) + gupxy_s * gupyy_s * Ayy(i,j,k) + gupxz_s * gupyz_s * Azz(i,j,k) + &
              (gupxx_s * gupyy_s       + gupxy_s * gupxy_s)* Axy(i,j,k)                       + &
              (gupxx_s * gupyz_s       + gupxz_s * gupxy_s)* Axz(i,j,k)                       + &
              (gupxy_s * gupyz_s       + gupxz_s * gupyy_s)* Ayz(i,j,k)

        Rxz_s =    gupxx_s * gupxz_s * Axx(i,j,k) + gupxy_s * gupyz_s * Ayy(i,j,k) + gupxz_s * gupzz_s * Azz(i,j,k) + &
              (gupxx_s * gupyz_s       + gupxy_s * gupxz_s)* Axy(i,j,k)                       + &
              (gupxx_s * gupzz_s       + gupxz_s * gupxz_s)* Axz(i,j,k)                       + &
              (gupxy_s * gupzz_s       + gupxz_s * gupyz_s)* Ayz(i,j,k)

        Ryz_s =    gupxy_s * gupxz_s * Axx(i,j,k) + gupyy_s * gupyz_s * Ayy(i,j,k) + gupyz_s * gupzz_s * Azz(i,j,k) + &
              (gupxy_s * gupyz_s       + gupyy_s * gupxz_s)* Axy(i,j,k)                       + &
              (gupxy_s * gupzz_s       + gupyz_s * gupxz_s)* Axz(i,j,k)                       + &
              (gupyy_s * gupzz_s       + gupyz_s * gupyz_s)* Ayz(i,j,k)

! Right hand side for Gam^i without shift terms...

        Gamx_rhs_s = - TWO * (   Lapx(i,j,k) * Rxx_s +   Lapy(i,j,k) * Rxy_s +   Lapz(i,j,k) * Rxz_s ) + &
             TWO * alpn1_s * (                                                &
             -F3o2/chin1_s * (   chix(i,j,k) * Rxx_s +   chiy(i,j,k) * Rxy_s +   chiz(i,j,k) * Rxz_s ) - &
                   gupxx_s * (   F2o3 * Kx(i,j,k)  +  EIGHT * PI * Sx(i,j,k)            ) - &
                   gupxy_s * (   F2o3 * Ky(i,j,k)  +  EIGHT * PI * Sy(i,j,k)            ) - &
                   gupxz_s * (   F2o3 * Kz(i,j,k)  +  EIGHT * PI * Sz(i,j,k)            ) + &
                             Gamxxx_s * Rxx_s + Gamxyy_s * Ryy_s + Gamxzz_s * Rzz_s   + &
                   TWO * ( Gamxxy_s * Rxy_s + Gamxxz_s * Rxz_s + Gamxyz_s * Ryz_s ) )

        Gamy_rhs_s = - TWO * (   Lapx(i,j,k) * Rxy_s +   Lapy(i,j,k) * Ryy_s +   Lapz(i,j,k) * Ryz_s ) + &
             TWO * alpn1_s * (                                                &
             -F3o2/chin1_s * (   chix(i,j,k) * Rxy_s +  chiy(i,j,k) * Ryy_s +    chiz(i,j,k) * Ryz_s ) - &
                   gupxy_s * (   F2o3 * Kx(i,j,k)  +  EIGHT * PI * Sx(i,j,k)            ) - &
                   gupyy_s * (   F2o3 * Ky(i,j,k)  +  EIGHT * PI * Sy(i,j,k)            ) - &
                   gupyz_s * (   F2o3 * Kz(i,j,k)  +  EIGHT * PI * Sz(i,j,k)            ) + &
                             Gamyxx_s * Rxx_s + Gamyyy_s * Ryy_s + Gamyzz_s * Rzz_s   + &
                   TWO * ( Gamyxy_s * Rxy_s + Gamyxz_s * Rxz_s + Gamyyz_s * Ryz_s ) )

        Gamz_rhs_s = - TWO * (   Lapx(i,j,k) * Rxz_s +   Lapy(i,j,k) * Ryz_s +   Lapz(i,j,k) * Rzz_s ) + &
             TWO * alpn1_s * (                                                &
             -F3o2/chin1_s * (   chix(i,j,k) * Rxz_s +  chiy(i,j,k) * Ryz_s +    chiz(i,j,k) * Rzz_s ) - &
                   gupxz_s * (   F2o3 * Kx(i,j,k)  +  EIGHT * PI * Sx(i,j,k)            ) - &
                   gupyz_s * (   F2o3 * Ky(i,j,k)  +  EIGHT * PI * Sy(i,j,k)            ) - &
                   gupzz_s * (   F2o3 * Kz(i,j,k)  +  EIGHT * PI * Sz(i,j,k)            ) + &
                             Gamzxx_s * Rxx_s + Gamzyy_s * Ryy_s + Gamzzz_s * Rzz_s   + &
                   TWO * ( Gamzxy_s * Rxy_s + Gamzxz_s * Rxz_s + Gamzyz_s * Ryz_s ) )

! fxx,fxy,fxz from 2nd generation (fdderivs of shift vector)
#ifdef AMSS_RHS_POINTWISE
        if (pointwise_mode) then
          fxx_s = pw_bxxx + pw_bxyy + pw_bxzz
          fxy_s = pw_bxyx + pw_byyy + pw_byzz
          fxz_s = pw_bxzx + pw_byzy + pw_bzzz
        else
#endif
        fxx_s = bxxx(i,j,k) + bxyy(i,j,k) + bxzz(i,j,k)
        fxy_s = bxyx(i,j,k) + byyy(i,j,k) + byzz(i,j,k)
        fxz_s = bxzx(i,j,k) + byzy(i,j,k) + bzzz(i,j,k)
#ifdef AMSS_RHS_POINTWISE
        end if
#endif

        Gamxa_s =       gupxx_s * Gamxxx_s + gupyy_s * Gamxyy_s + gupzz_s * Gamxzz_s + &
              TWO*( gupxy_s * Gamxxy_s + gupxz_s * Gamxxz_s + gupyz_s * Gamxyz_s )
        Gamya_s =       gupxx_s * Gamyxx_s + gupyy_s * Gamyyy_s + gupzz_s * Gamyzz_s + &
              TWO*( gupxy_s * Gamyxy_s + gupxz_s * Gamyxz_s + gupyz_s * Gamyyz_s )
        Gamza_s =       gupxx_s * Gamzxx_s + gupyy_s * Gamzyy_s + gupzz_s * Gamzzz_s + &
              TWO*( gupxy_s * Gamzxy_s + gupxz_s * Gamzxz_s + gupyz_s * Gamzyz_s )

#ifdef AMSS_RHS_POINTWISE
        if (pointwise_mode) then
          Gamx_rhs_s =               Gamx_rhs_s +  F2o3 *  Gamxa_s * div_beta_s        - &
                        Gamxa_s * betaxx(i,j,k) - Gamya_s * betaxy(i,j,k) - Gamza_s * betaxz(i,j,k)  + &
                F1o3 * (gupxx_s * fxx_s    + gupxy_s * fxy_s    + gupxz_s * fxz_s    ) + &
                        gupxx_s * pw_bxxx   + gupyy_s * pw_byyx   + gupzz_s * pw_bzzx    + &
                 TWO * (gupxy_s * pw_bxyx   + gupxz_s * pw_bxzx   + gupyz_s * pw_byzx  )
          Gamy_rhs_s =               Gamy_rhs_s +  F2o3 *  Gamya_s * div_beta_s        - &
                        Gamxa_s * betayx(i,j,k) - Gamya_s * betayy(i,j,k) - Gamza_s * betayz(i,j,k)  + &
                F1o3 * (gupxy_s * fxx_s    + gupyy_s * fxy_s    + gupyz_s * fxz_s    ) + &
                        gupxx_s * pw_bxxy   + gupyy_s * pw_byyy   + gupzz_s * pw_bzzy    + &
                 TWO * (gupxy_s * pw_bxyy   + gupxz_s * pw_bxzy   + gupyz_s * pw_byzy  )
          Gamz_rhs_s =               Gamz_rhs_s +  F2o3 *  Gamza_s * div_beta_s        - &
                        Gamxa_s * betazx(i,j,k) - Gamya_s * betazy(i,j,k) - Gamza_s * betazz(i,j,k)  + &
                F1o3 * (gupxz_s * fxx_s    + gupyz_s * fxy_s    + gupzz_s * fxz_s    ) + &
                        gupxx_s * pw_bxxz   + gupyy_s * pw_byyz   + gupzz_s * pw_bzzz    + &
                 TWO * (gupxy_s * pw_bxyz   + gupxz_s * pw_bxzz   + gupyz_s * pw_byzz  )    !rhs for Gam^i
#ifdef AMSS_RHS_POINTWISE
        else
#endif
          Gamx_rhs_s =               Gamx_rhs_s +  F2o3 *  Gamxa_s * div_beta_s        - &
                        Gamxa_s * betaxx(i,j,k) - Gamya_s * betaxy(i,j,k) - Gamza_s * betaxz(i,j,k)  + &
                F1o3 * (gupxx_s * fxx_s    + gupxy_s * fxy_s    + gupxz_s * fxz_s    ) + &
                        gupxx_s * bxxx(i,j,k)   + gupyy_s * byyx(i,j,k)   + gupzz_s * bzzx(i,j,k)    + &
                 TWO * (gupxy_s * bxyx(i,j,k)   + gupxz_s * bxzx(i,j,k)   + gupyz_s * byzx(i,j,k)  )
          Gamy_rhs_s =               Gamy_rhs_s +  F2o3 *  Gamya_s * div_beta_s        - &
                        Gamxa_s * betayx(i,j,k) - Gamya_s * betayy(i,j,k) - Gamza_s * betayz(i,j,k)  + &
                F1o3 * (gupxy_s * fxx_s    + gupyy_s * fxy_s    + gupyz_s * fxz_s    ) + &
                        gupxx_s * bxxy(i,j,k)   + gupyy_s * byyy(i,j,k)   + gupzz_s * bzzy(i,j,k)    + &
                 TWO * (gupxy_s * bxyy(i,j,k)   + gupxz_s * bxzy(i,j,k)   + gupyz_s * byzy(i,j,k)  )
          Gamz_rhs_s =               Gamz_rhs_s +  F2o3 *  Gamza_s * div_beta_s        - &
                        Gamxa_s * betazx(i,j,k) - Gamya_s * betazy(i,j,k) - Gamza_s * betazz(i,j,k)  + &
                F1o3 * (gupxz_s * fxx_s    + gupyz_s * fxy_s    + gupzz_s * fxz_s    ) + &
                        gupxx_s * bxxz(i,j,k)   + gupyy_s * byyz(i,j,k)   + gupzz_s * bzzz(i,j,k)    + &
                 TWO * (gupxy_s * bxyz(i,j,k)   + gupxz_s * bxzz(i,j,k)   + gupyz_s * byzz(i,j,k)  )    !rhs for Gam^i
#ifdef AMSS_RHS_POINTWISE
        end if
#endif
#endif

!first kind of connection stored in gij,k  (3rd generation, scalar)
#ifdef AMSS_RHS_REORDER
        ! Phase 4 (lifetime scheduling): Gamma_rhs is now complete (shift terms
        ! applied above).  Store it to the array and finish the Gamma-driver
        ! (dtSf_rhs) immediately, so Gamx_rhs_s / Gamy_rhs_s / Gamz_rhs_s and
        ! the raised-Gamma + shift-Hessian intermediates (Gamxa_s/Gamya_s/
        ! Gamza_s, fxx_s/fxy_s/fxz_s) can be released before the Ricci block.
        ! NOTE: Gamxa_s/Gamya_s/Gamza_s are ALSO read by Ricci (the Γ^i·Γ_{ijk}
        ! term at 2201+), so they survive past this point -- only the three
        ! Gamma_rhs scalars and the shift-Hessian fxx_s/fxy_s/fxz_s die here.
        Gamx_rhs(i,j,k) = Gamx_rhs_s
        Gamy_rhs(i,j,k) = Gamy_rhs_s
        Gamz_rhs(i,j,k) = Gamz_rhs_s
        dtSfx_rhs(i,j,k) = Gamx_rhs_s - eta*dtSfx(i,j,k)
        dtSfy_rhs(i,j,k) = Gamy_rhs_s - eta*dtSfy(i,j,k)
        dtSfz_rhs(i,j,k) = Gamz_rhs_s - eta*dtSfz(i,j,k)
#endif

        gxxx_t = gxx_s * Gamxxx_s + gxy(i,j,k) * Gamyxx_s + gxz(i,j,k) * Gamzxx_s
        gxyx_t = gxx_s * Gamxxy_s + gxy(i,j,k) * Gamyxy_s + gxz(i,j,k) * Gamzxy_s
        gxzx_t = gxx_s * Gamxxz_s + gxy(i,j,k) * Gamyxz_s + gxz(i,j,k) * Gamzxz_s
        gyyx_t = gxx_s * Gamxyy_s + gxy(i,j,k) * Gamyyy_s + gxz(i,j,k) * Gamzyy_s
        gyzx_t = gxx_s * Gamxyz_s + gxy(i,j,k) * Gamyyz_s + gxz(i,j,k) * Gamzyz_s
        gzzx_t = gxx_s * Gamxzz_s + gxy(i,j,k) * Gamyzz_s + gxz(i,j,k) * Gamzzz_s

        gxxy_t = gxy(i,j,k) * Gamxxx_s + gyy_s * Gamyxx_s + gyz(i,j,k) * Gamzxx_s
        gxyy_t = gxy(i,j,k) * Gamxxy_s + gyy_s * Gamyxy_s + gyz(i,j,k) * Gamzxy_s
        gxzy_t = gxy(i,j,k) * Gamxxz_s + gyy_s * Gamyxz_s + gyz(i,j,k) * Gamzxz_s
        gyyy_t = gxy(i,j,k) * Gamxyy_s + gyy_s * Gamyyy_s + gyz(i,j,k) * Gamzyy_s
        gyzy_t = gxy(i,j,k) * Gamxyz_s + gyy_s * Gamyyz_s + gyz(i,j,k) * Gamzyz_s
        gzzy_t = gxy(i,j,k) * Gamxzz_s + gyy_s * Gamyzz_s + gyz(i,j,k) * Gamzzz_s

        gxxz_t = gxz(i,j,k) * Gamxxx_s + gyz(i,j,k) * Gamyxx_s + gzz_s * Gamzxx_s
        gxyz_t = gxz(i,j,k) * Gamxxy_s + gyz(i,j,k) * Gamyxy_s + gzz_s * Gamzxy_s
        gxzz_t = gxz(i,j,k) * Gamxxz_s + gyz(i,j,k) * Gamyxz_s + gzz_s * Gamzxz_s
        gyyz_t = gxz(i,j,k) * Gamxyy_s + gyz(i,j,k) * Gamyyy_s + gzz_s * Gamzyy_s
        gyzz_t = gxz(i,j,k) * Gamxyz_s + gyz(i,j,k) * Gamyyz_s + gzz_s * Gamzyz_s
        gzzz_t = gxz(i,j,k) * Gamxzz_s + gyz(i,j,k) * Gamyzz_s + gzz_s * Gamzzz_s

!compute Ricci tensor for tilted metric
#ifdef AMSS_RHS_POINTWISE
        if (pointwise_mode) then
          Rxx_s =   gupxx_s * pw_fxx_dxx + gupyy_s * pw_fyy_dxx + gupzz_s * pw_fzz_dxx + &
                ( gupxy_s * pw_fxy_dxx + gupxz_s * pw_fxz_dxx + gupyz_s * pw_fyz_dxx ) * TWO
          Ryy_s =   gupxx_s * pw_fxx_dyy + gupyy_s * pw_fyy_dyy + gupzz_s * pw_fzz_dyy + &
                ( gupxy_s * pw_fxy_dyy + gupxz_s * pw_fxz_dyy + gupyz_s * pw_fyz_dyy ) * TWO
          Rzz_s =   gupxx_s * pw_fxx_dzz + gupyy_s * pw_fyy_dzz + gupzz_s * pw_fzz_dzz + &
                ( gupxy_s * pw_fxy_dzz + gupxz_s * pw_fxz_dzz + gupyz_s * pw_fyz_dzz ) * TWO
          Rxy_s =   gupxx_s * pw_fxx_gxy + gupyy_s * pw_fyy_gxy + gupzz_s * pw_fzz_gxy + &
                ( gupxy_s * pw_fxy_gxy + gupxz_s * pw_fxz_gxy + gupyz_s * pw_fyz_gxy ) * TWO
          Rxz_s =   gupxx_s * pw_fxx_gxz + gupyy_s * pw_fyy_gxz + gupzz_s * pw_fzz_gxz + &
                ( gupxy_s * pw_fxy_gxz + gupxz_s * pw_fxz_gxz + gupyz_s * pw_fyz_gxz ) * TWO
          Ryz_s =   gupxx_s * pw_fxx_gyz + gupyy_s * pw_fyy_gyz + gupzz_s * pw_fzz_gyz + &
                ( gupxy_s * pw_fxy_gyz + gupxz_s * pw_fxz_gyz + gupyz_s * pw_fyz_gyz ) * TWO
#ifdef AMSS_RHS_POINTWISE
        else
#endif
          Rxx_s =   gupxx_s * fxx_dxx(i,j,k) + gupyy_s * fyy_dxx(i,j,k) + gupzz_s * fzz_dxx(i,j,k) + &
                ( gupxy_s * fxy_dxx(i,j,k) + gupxz_s * fxz_dxx(i,j,k) + gupyz_s * fyz_dxx(i,j,k) ) * TWO
          Ryy_s =   gupxx_s * fxx_dyy(i,j,k) + gupyy_s * fyy_dyy(i,j,k) + gupzz_s * fzz_dyy(i,j,k) + &
                ( gupxy_s * fxy_dyy(i,j,k) + gupxz_s * fxz_dyy(i,j,k) + gupyz_s * fyz_dyy(i,j,k) ) * TWO
          Rzz_s =   gupxx_s * fxx_dzz(i,j,k) + gupyy_s * fyy_dzz(i,j,k) + gupzz_s * fzz_dzz(i,j,k) + &
                ( gupxy_s * fxy_dzz(i,j,k) + gupxz_s * fxz_dzz(i,j,k) + gupyz_s * fyz_dzz(i,j,k) ) * TWO
          Rxy_s =   gupxx_s * fxx_gxy(i,j,k) + gupyy_s * fyy_gxy(i,j,k) + gupzz_s * fzz_gxy(i,j,k) + &
                ( gupxy_s * fxy_gxy(i,j,k) + gupxz_s * fxz_gxy(i,j,k) + gupyz_s * fyz_gxy(i,j,k) ) * TWO
          Rxz_s =   gupxx_s * fxx_gxz(i,j,k) + gupyy_s * fyy_gxz(i,j,k) + gupzz_s * fzz_gxz(i,j,k) + &
                ( gupxy_s * fxy_gxz(i,j,k) + gupxz_s * fxz_gxz(i,j,k) + gupyz_s * fyz_gxz(i,j,k) ) * TWO
          Ryz_s =   gupxx_s * fxx_gyz(i,j,k) + gupyy_s * fyy_gyz(i,j,k) + gupzz_s * fzz_gyz(i,j,k) + &
                ( gupxy_s * fxy_gyz(i,j,k) + gupxz_s * fxz_gyz(i,j,k) + gupyz_s * fyz_gyz(i,j,k) ) * TWO
#ifdef AMSS_RHS_POINTWISE
        end if
#endif
#endif

        Rxx_s =     - HALF * Rxx_s                                   + &
                     gxx_s * Gamxx(i,j,k)+ gxy(i,j,k) * Gamyx(i,j,k)   +    gxz(i,j,k) * Gamzx(i,j,k) + &
                   Gamxa_s * gxxx_t +  Gamya_s * gxyx_t +  Gamza_s * gxzx_t  + &
           gupxx_s *(                                                  &
               TWO*(Gamxxx_s * gxxx_t + Gamyxx_s * gxyx_t + Gamzxx_s * gxzx_t) + &
                    Gamxxx_s * gxxx_t + Gamyxx_s * gxxy_t + Gamzxx_s * gxxz_t )+ &
           gupxy_s *(                                                  &
               TWO*(Gamxxx_s * gxyx_t + Gamyxx_s * gyyx_t + Gamzxx_s * gyzx_t  + &
                    Gamxxy_s * gxxx_t + Gamyxy_s * gxyx_t + Gamzxy_s * gxzx_t) + &
                    Gamxxy_s * gxxx_t + Gamyxy_s * gxxy_t + Gamzxy_s * gxxz_t  + &
                    Gamxxx_s * gxyx_t + Gamyxx_s * gxyy_t + Gamzxx_s * gxyz_t )+ &
           gupxz_s *(                                                  &
               TWO*(Gamxxx_s * gxzx_t + Gamyxx_s * gyzx_t + Gamzxx_s * gzzx_t  + &
                    Gamxxz_s * gxxx_t + Gamyxz_s * gxyx_t + Gamzxz_s * gxzx_t) + &
                    Gamxxz_s * gxxx_t + Gamyxz_s * gxxy_t + Gamzxz_s * gxxz_t  + &
                    Gamxxx_s * gxzx_t + Gamyxx_s * gxzy_t + Gamzxx_s * gxzz_t )+ &
           gupyy_s *(                                                  &
               TWO*(Gamxxy_s * gxyx_t + Gamyxy_s * gyyx_t + Gamzxy_s * gyzx_t) + &
                    Gamxxy_s * gxyx_t + Gamyxy_s * gxyy_t + Gamzxy_s * gxyz_t )+ &
           gupyz_s *(                                                  &
               TWO*(Gamxxy_s * gxzx_t + Gamyxy_s * gyzx_t + Gamzxy_s * gzzx_t  + &
                    Gamxxz_s * gxyx_t + Gamyxz_s * gyyx_t + Gamzxz_s * gyzx_t) + &
                    Gamxxz_s * gxyx_t + Gamyxz_s * gxyy_t + Gamzxz_s * gxyz_t  + &
                    Gamxxy_s * gxzx_t + Gamyxy_s * gxzy_t + Gamzxy_s * gxzz_t )+ &
           gupzz_s *(                                                  &
               TWO*(Gamxxz_s * gxzx_t + Gamyxz_s * gyzx_t + Gamzxz_s * gzzx_t) + &
                    Gamxxz_s * gxzx_t + Gamyxz_s * gxzy_t + Gamzxz_s * gxzz_t )

        Ryy_s =     - HALF * Ryy_s                                   + &
                     gxy(i,j,k) * Gamxy(i,j,k)+  gyy_s * Gamyy(i,j,k)  +  gyz(i,j,k) * Gamzy(i,j,k)   + &
                   Gamxa_s * gxyy_t +  Gamya_s * gyyy_t +  Gamza_s * gyzy_t  + &
           gupxx_s *(                                                  &
               TWO*(Gamxxy_s * gxxy_t + Gamyxy_s * gxyy_t + Gamzxy_s * gxzy_t) + &
                    Gamxxy_s * gxyx_t + Gamyxy_s * gxyy_t + Gamzxy_s * gxyz_t )+ &
           gupxy_s *(                                                  &
               TWO*(Gamxxy_s * gxyy_t + Gamyxy_s * gyyy_t + Gamzxy_s * gyzy_t  + &
                    Gamxyy_s * gxxy_t + Gamyyy_s * gxyy_t + Gamzyy_s * gxzy_t) + &
                    Gamxyy_s * gxyx_t + Gamyyy_s * gxyy_t + Gamzyy_s * gxyz_t  + &
                    Gamxxy_s * gyyx_t + Gamyxy_s * gyyy_t + Gamzxy_s * gyyz_t )+ &
           gupxz_s *(                                                  &
               TWO*(Gamxxy_s * gxzy_t + Gamyxy_s * gyzy_t + Gamzxy_s * gzzy_t  + &
                    Gamxyz_s * gxxy_t + Gamyyz_s * gxyy_t + Gamzyz_s * gxzy_t) + &
                    Gamxyz_s * gxyx_t + Gamyyz_s * gxyy_t + Gamzyz_s * gxyz_t  + &
                    Gamxxy_s * gyzx_t + Gamyxy_s * gyzy_t + Gamzxy_s * gyzz_t )+ &
           gupyy_s *(                                                  &
               TWO*(Gamxyy_s * gxyy_t + Gamyyy_s * gyyy_t + Gamzyy_s * gyzy_t) + &
                    Gamxyy_s * gyyx_t + Gamyyy_s * gyyy_t + Gamzyy_s * gyyz_t )+ &
           gupyz_s *(                                                  &
               TWO*(Gamxyy_s * gxzy_t + Gamyyy_s * gyzy_t + Gamzyy_s * gzzy_t  + &
                    Gamxyz_s * gxyy_t + Gamyyz_s * gyyy_t + Gamzyz_s * gyzy_t) + &
                    Gamxyz_s * gyyx_t + Gamyyz_s * gyyy_t + Gamzyz_s * gyyz_t  + &
                    Gamxyy_s * gyzx_t + Gamyyy_s * gyzy_t + Gamzyy_s * gyzz_t )+ &
           gupzz_s *(                                                  &
               TWO*(Gamxyz_s * gxzy_t + Gamyyz_s * gyzy_t + Gamzyz_s * gzzy_t) + &
                    Gamxyz_s * gyzx_t + Gamyyz_s * gyzy_t + Gamzyz_s * gyzz_t )

        Rzz_s =     - HALF * Rzz_s                                   + &
                     gxz(i,j,k) * Gamxz(i,j,k)+ gyz(i,j,k) * Gamyz(i,j,k)  +    gzz_s * Gamzz(i,j,k)  + &
                   Gamxa_s * gxzz_t +  Gamya_s * gyzz_t +  Gamza_s * gzzz_t  + &
           gupxx_s *(                                                  &
               TWO*(Gamxxz_s * gxxz_t + Gamyxz_s * gxyz_t + Gamzxz_s * gxzz_t) + &
                    Gamxxz_s * gxzx_t + Gamyxz_s * gxzy_t + Gamzxz_s * gxzz_t )+ &
           gupxy_s *(                                                  &
               TWO*(Gamxxz_s * gxyz_t + Gamyxz_s * gyyz_t + Gamzxz_s * gyzz_t  + &
                    Gamxyz_s * gxxz_t + Gamyyz_s * gxyz_t + Gamzyz_s * gxzz_t) + &
                    Gamxyz_s * gxzx_t + Gamyyz_s * gxzy_t + Gamzyz_s * gxzz_t  + &
                    Gamxxz_s * gyzx_t + Gamyxz_s * gyzy_t + Gamzxz_s * gyzz_t )+ &
           gupxz_s *(                                                  &
               TWO*(Gamxxz_s * gxzz_t + Gamyxz_s * gyzz_t + Gamzxz_s * gzzz_t  + &
                    Gamxzz_s * gxxz_t + Gamyzz_s * gxyz_t + Gamzzz_s * gxzz_t) + &
                    Gamxzz_s * gxzx_t + Gamyzz_s * gxzy_t + Gamzzz_s * gxzz_t  + &
                    Gamxxz_s * gzzx_t + Gamyxz_s * gzzy_t + Gamzxz_s * gzzz_t )+ &
           gupyy_s *(                                                  &
               TWO*(Gamxyz_s * gxyz_t + Gamyyz_s * gyyz_t + Gamzyz_s * gyzz_t) + &
                    Gamxyz_s * gyzx_t + Gamyyz_s * gyzy_t + Gamzyz_s * gyzz_t )+ &
           gupyz_s *(                                                  &
               TWO*(Gamxyz_s * gxzz_t + Gamyyz_s * gyzz_t + Gamzyz_s * gzzz_t  + &
                    Gamxzz_s * gxyz_t + Gamyzz_s * gyyz_t + Gamzzz_s * gyzz_t) + &
                    Gamxzz_s * gyzx_t + Gamyzz_s * gyzy_t + Gamzzz_s * gyzz_t  + &
                    Gamxyz_s * gzzx_t + Gamyyz_s * gzzy_t + Gamzyz_s * gzzz_t )+ &
           gupzz_s *(                                                  &
               TWO*(Gamxzz_s * gxzz_t + Gamyzz_s * gyzz_t + Gamzzz_s * gzzz_t) + &
                    Gamxzz_s * gzzx_t + Gamyzz_s * gzzy_t + Gamzzz_s * gzzz_t )

        Rxy_s = HALF*(     - Rxy_s                                   + &
                     gxx_s * Gamxy(i,j,k) +    gxy(i,j,k) * Gamyy(i,j,k) + gxz(i,j,k) * Gamzy(i,j,k)  + &
                     gxy(i,j,k) * Gamxx(i,j,k) +    gyy_s * Gamyx(i,j,k) + gyz(i,j,k) * Gamzx(i,j,k)  + &
                   Gamxa_s * gxyx_t +  Gamya_s * gyyx_t +  Gamza_s * gyzx_t  + &
                   Gamxa_s * gxxy_t +  Gamya_s * gxyy_t +  Gamza_s * gxzy_t )+ &
           gupxx_s *(                                                  &
                    Gamxxx_s * gxxy_t + Gamyxx_s * gxyy_t + Gamzxx_s * gxzy_t  + &
                    Gamxxy_s * gxxx_t + Gamyxy_s * gxyx_t + Gamzxy_s * gxzx_t  + &
                    Gamxxx_s * gxyx_t + Gamyxx_s * gxyy_t + Gamzxx_s * gxyz_t )+ &
           gupxy_s *(                                                  &
                    Gamxxx_s * gxyy_t + Gamyxx_s * gyyy_t + Gamzxx_s * gyzy_t  + &
                    Gamxxy_s * gxyx_t + Gamyxy_s * gyyx_t + Gamzxy_s * gyzx_t  + &
                    Gamxxy_s * gxyx_t + Gamyxy_s * gxyy_t + Gamzxy_s * gxyz_t  + &
                    Gamxxy_s * gxxy_t + Gamyxy_s * gxyy_t + Gamzxy_s * gxzy_t  + &
                    Gamxyy_s * gxxx_t + Gamyyy_s * gxyx_t + Gamzyy_s * gxzx_t  + &
                    Gamxxx_s * gyyx_t + Gamyxx_s * gyyy_t + Gamzxx_s * gyyz_t )+ &
           gupxz_s *(                                                  &
                    Gamxxx_s * gxzy_t + Gamyxx_s * gyzy_t + Gamzxx_s * gzzy_t  + &
                    Gamxxy_s * gxzx_t + Gamyxy_s * gyzx_t + Gamzxy_s * gzzx_t  + &
                    Gamxxz_s * gxyx_t + Gamyxz_s * gxyy_t + Gamzxz_s * gxyz_t  + &
                    Gamxxz_s * gxxy_t + Gamyxz_s * gxyy_t + Gamzxz_s * gxzy_t  + &
                    Gamxyz_s * gxxx_t + Gamyyz_s * gxyx_t + Gamzyz_s * gxzx_t  + &
                    Gamxxx_s * gyzx_t + Gamyxx_s * gyzy_t + Gamzxx_s * gyzz_t )+ &
           gupyy_s *(                                                  &
                    Gamxxy_s * gxyy_t + Gamyxy_s * gyyy_t + Gamzxy_s * gyzy_t  + &
                    Gamxyy_s * gxyx_t + Gamyyy_s * gyyx_t + Gamzyy_s * gyzx_t  + &
                    Gamxxy_s * gyyx_t + Gamyxy_s * gyyy_t + Gamzxy_s * gyyz_t )+ &
           gupyz_s *(                                                  &
                    Gamxxy_s * gxzy_t + Gamyxy_s * gyzy_t + Gamzxy_s * gzzy_t  + &
                    Gamxyy_s * gxzx_t + Gamyyy_s * gyzx_t + Gamzyy_s * gzzx_t  + &
                    Gamxxz_s * gyyx_t + Gamyxz_s * gyyy_t + Gamzxz_s * gyyz_t  + &
                    Gamxxz_s * gxyy_t + Gamyxz_s * gyyy_t + Gamzxz_s * gyzy_t  + &
                    Gamxyz_s * gxyx_t + Gamyyz_s * gyyx_t + Gamzyz_s * gyzx_t  + &
                    Gamxxy_s * gyzx_t + Gamyxy_s * gyzy_t + Gamzxy_s * gyzz_t )+ &
           gupzz_s *(                                                  &
                    Gamxxz_s * gxzy_t + Gamyxz_s * gyzy_t + Gamzxz_s * gzzy_t  + &
                    Gamxyz_s * gxzx_t + Gamyyz_s * gyzx_t + Gamzyz_s * gzzx_t  + &
                    Gamxxz_s * gyzx_t + Gamyxz_s * gyzy_t + Gamzxz_s * gyzz_t )

        Rxz_s = HALF*(     - Rxz_s                                   + &
                     gxx_s * Gamxz(i,j,k) +  gxy(i,j,k) * Gamyz(i,j,k) + gxz(i,j,k) * Gamzz(i,j,k)    + &
                     gxz(i,j,k) * Gamxx(i,j,k) +  gyz(i,j,k) * Gamyx(i,j,k) + gzz_s * Gamzx(i,j,k)    + &
                   Gamxa_s * gxzx_t +  Gamya_s * gyzx_t +  Gamza_s * gzzx_t  + &
                   Gamxa_s * gxxz_t +  Gamya_s * gxyz_t +  Gamza_s * gxzz_t )+ &
           gupxx_s *(                                                  &
                    Gamxxx_s * gxxz_t + Gamyxx_s * gxyz_t + Gamzxx_s * gxzz_t  + &
                    Gamxxz_s * gxxx_t + Gamyxz_s * gxyx_t + Gamzxz_s * gxzx_t  + &
                    Gamxxx_s * gxzx_t + Gamyxx_s * gxzy_t + Gamzxx_s * gxzz_t )+ &
           gupxy_s *(                                                  &
                    Gamxxx_s * gxyz_t + Gamyxx_s * gyyz_t + Gamzxx_s * gyzz_t  + &
                    Gamxxz_s * gxyx_t + Gamyxz_s * gyyx_t + Gamzxz_s * gyzx_t  + &
                    Gamxxy_s * gxzx_t + Gamyxy_s * gxzy_t + Gamzxy_s * gxzz_t  + &
                    Gamxxy_s * gxxz_t + Gamyxy_s * gxyz_t + Gamzxy_s * gxzz_t  + &
                    Gamxyz_s * gxxx_t + Gamyyz_s * gxyx_t + Gamzyz_s * gxzx_t  + &
                    Gamxxx_s * gyzx_t + Gamyxx_s * gyzy_t + Gamzxx_s * gyzz_t )+ &
           gupxz_s *(                                                  &
                    Gamxxx_s * gxzz_t + Gamyxx_s * gyzz_t + Gamzxx_s * gzzz_t  + &
                    Gamxxz_s * gxzx_t + Gamyxz_s * gyzx_t + Gamzxz_s * gzzx_t  + &
                    Gamxxz_s * gxzx_t + Gamyxz_s * gxzy_t + Gamzxz_s * gxzz_t  + &
                    Gamxxz_s * gxxz_t + Gamyxz_s * gxyz_t + Gamzxz_s * gxzz_t  + &
                    Gamxzz_s * gxxx_t + Gamyzz_s * gxyx_t + Gamzzz_s * gxzx_t  + &
                    Gamxxx_s * gzzx_t + Gamyxx_s * gzzy_t + Gamzxx_s * gzzz_t )+ &
           gupyy_s *(                                                  &
                    Gamxxy_s * gxyz_t + Gamyxy_s * gyyz_t + Gamzxy_s * gyzz_t  + &
                    Gamxyz_s * gxyx_t + Gamyyz_s * gyyx_t + Gamzyz_s * gyzx_t  + &
                    Gamxxy_s * gyzx_t + Gamyxy_s * gyzy_t + Gamzxy_s * gyzz_t )+ &
           gupyz_s *(                                                  &
                    Gamxxy_s * gxzz_t + Gamyxy_s * gyzz_t + Gamzxy_s * gzzz_t  + &
                    Gamxyz_s * gxzx_t + Gamyyz_s * gyzx_t + Gamzyz_s * gzzx_t  + &
                    Gamxxz_s * gyzx_t + Gamyxz_s * gyzy_t + Gamzxz_s * gyzz_t  + &
                    Gamxxz_s * gxyz_t + Gamyxz_s * gyyz_t + Gamzxz_s * gyzz_t  + &
                    Gamxzz_s * gxyx_t + Gamyzz_s * gyyx_t + Gamzzz_s * gyzx_t  + &
                    Gamxxy_s * gzzx_t + Gamyxy_s * gzzy_t + Gamzxy_s * gzzz_t )+ &
           gupzz_s *(                                                  &
                    Gamxxz_s * gxzz_t + Gamyxz_s * gyzz_t + Gamzxz_s * gzzz_t  + &
                    Gamxzz_s * gxzx_t + Gamyzz_s * gyzx_t + Gamzzz_s * gzzx_t  + &
                    Gamxxz_s * gzzx_t + Gamyxz_s * gzzy_t + Gamzxz_s * gzzz_t )

        Ryz_s = HALF*(     - Ryz_s                                   + &
                     gxy(i,j,k) * Gamxz(i,j,k) + gyy_s * Gamyz(i,j,k) + gyz(i,j,k) * Gamzz(i,j,k)     + &
                     gxz(i,j,k) * Gamxy(i,j,k) + gyz(i,j,k) * Gamyy(i,j,k) + gzz_s * Gamzy(i,j,k)     + &
                   Gamxa_s * gxzy_t +  Gamya_s * gyzy_t +  Gamza_s * gzzy_t  + &
                   Gamxa_s * gxyz_t +  Gamya_s * gyyz_t +  Gamza_s * gyzz_t )+ &
           gupxx_s *(                                                  &
                    Gamxxy_s * gxxz_t + Gamyxy_s * gxyz_t + Gamzxy_s * gxzz_t  + &
                    Gamxxz_s * gxxy_t + Gamyxz_s * gxyy_t + Gamzxz_s * gxzy_t  + &
                    Gamxxy_s * gxzx_t + Gamyxy_s * gxzy_t + Gamzxy_s * gxzz_t )+ &
           gupxy_s *(                                                  &
                    Gamxxy_s * gxyz_t + Gamyxy_s * gyyz_t + Gamzxy_s * gyzz_t  + &
                    Gamxxz_s * gxyy_t + Gamyxz_s * gyyy_t + Gamzxz_s * gyzy_t  + &
                    Gamxyy_s * gxzx_t + Gamyyy_s * gxzy_t + Gamzyy_s * gxzz_t  + &
                    Gamxyy_s * gxxz_t + Gamyyy_s * gxyz_t + Gamzyy_s * gxzz_t  + &
                    Gamxyz_s * gxxy_t + Gamyyz_s * gxyy_t + Gamzyz_s * gxzy_t  + &
                    Gamxxy_s * gyzx_t + Gamyxy_s * gyzy_t + Gamzxy_s * gyzz_t )+ &
           gupxz_s *(                                                  &
                    Gamxxy_s * gxzz_t + Gamyxy_s * gyzz_t + Gamzxy_s * gzzz_t  + &
                    Gamxxz_s * gxzy_t + Gamyxz_s * gyzy_t + Gamzxz_s * gzzy_t  + &
                    Gamxyz_s * gxzx_t + Gamyyz_s * gxzy_t + Gamzyz_s * gxzz_t  + &
                    Gamxyz_s * gxxz_t + Gamyyz_s * gxyz_t + Gamzyz_s * gxzz_t  + &
                    Gamxzz_s * gxxy_t + Gamyzz_s * gxyy_t + Gamzzz_s * gxzy_t  + &
                    Gamxxy_s * gzzx_t + Gamyxy_s * gzzy_t + Gamzxy_s * gzzz_t )+ &
           gupyy_s *(                                                  &
                    Gamxyy_s * gxyz_t + Gamyyy_s * gyyz_t + Gamzyy_s * gyzz_t  + &
                    Gamxyz_s * gxyy_t + Gamyyz_s * gyyy_t + Gamzyz_s * gyzy_t  + &
                    Gamxyy_s * gyzx_t + Gamyyy_s * gyzy_t + Gamzyy_s * gyzz_t )+ &
           gupyz_s *(                                                  &
                    Gamxyy_s * gxzz_t + Gamyyy_s * gyzz_t + Gamzyy_s * gzzz_t  + &
                    Gamxyz_s * gxzy_t + Gamyyz_s * gyzy_t + Gamzyz_s * gzzy_t  + &
                    Gamxyz_s * gyzx_t + Gamyyz_s * gyzy_t + Gamzyz_s * gyzz_t  + &
                    Gamxyz_s * gxyz_t + Gamyyz_s * gyyz_t + Gamzyz_s * gyzz_t  + &
                    Gamxzz_s * gxyy_t + Gamyzz_s * gyyy_t + Gamzzz_s * gyzy_t  + &
                    Gamxyy_s * gzzx_t + Gamyyy_s * gzzy_t + Gamzyy_s * gzzz_t )+ &
           gupzz_s *(                                                  &
                    Gamxyz_s * gxzz_t + Gamyyz_s * gyzz_t + Gamzyz_s * gzzz_t  + &
                    Gamxzz_s * gxzy_t + Gamyzz_s * gyzy_t + Gamzzz_s * gzzy_t  + &
                    Gamxyz_s * gzzx_t + Gamyyz_s * gzzy_t + Gamzyz_s * gzzz_t )
!covariant second derivative of chi respect to tilted metric
#ifdef AMSS_RHS_POINTWISE
        if (pointwise_mode) then
          fxx_s = pw_fxx_chi - Gamxxx_s * chix(i,j,k) - Gamyxx_s * chiy(i,j,k) - Gamzxx_s * chiz(i,j,k)
          fxy_s = pw_fxy_chi - Gamxxy_s * chix(i,j,k) - Gamyxy_s * chiy(i,j,k) - Gamzxy_s * chiz(i,j,k)
          fxz_s = pw_fxz_chi - Gamxxz_s * chix(i,j,k) - Gamyxz_s * chiy(i,j,k) - Gamzxz_s * chiz(i,j,k)
          fyy_s = pw_fyy_chi - Gamxyy_s * chix(i,j,k) - Gamyyy_s * chiy(i,j,k) - Gamzyy_s * chiz(i,j,k)
          fyz_s = pw_fyz_chi - Gamxyz_s * chix(i,j,k) - Gamyyz_s * chiy(i,j,k) - Gamzyz_s * chiz(i,j,k)
          fzz_s = pw_fzz_chi - Gamxzz_s * chix(i,j,k) - Gamyzz_s * chiy(i,j,k) - Gamzzz_s * chiz(i,j,k)
#ifdef AMSS_RHS_POINTWISE
        else
#endif
          fxx_s = fxx_chi(i,j,k) - Gamxxx_s * chix(i,j,k) - Gamyxx_s * chiy(i,j,k) - Gamzxx_s * chiz(i,j,k)
          fxy_s = fxy_chi(i,j,k) - Gamxxy_s * chix(i,j,k) - Gamyxy_s * chiy(i,j,k) - Gamzxy_s * chiz(i,j,k)
          fxz_s = fxz_chi(i,j,k) - Gamxxz_s * chix(i,j,k) - Gamyxz_s * chiy(i,j,k) - Gamzxz_s * chiz(i,j,k)
          fyy_s = fyy_chi(i,j,k) - Gamxyy_s * chix(i,j,k) - Gamyyy_s * chiy(i,j,k) - Gamzyy_s * chiz(i,j,k)
          fyz_s = fyz_chi(i,j,k) - Gamxyz_s * chix(i,j,k) - Gamyyz_s * chiy(i,j,k) - Gamzyz_s * chiz(i,j,k)
          fzz_s = fzz_chi(i,j,k) - Gamxzz_s * chix(i,j,k) - Gamyzz_s * chiy(i,j,k) - Gamzzz_s * chiz(i,j,k)
#ifdef AMSS_RHS_POINTWISE
        end if
#endif
#endif
! Store D^l D_l chi - 3/(2*chi) D^l chi D_l chi in f

        f_s =        gupxx_s * ( fxx_s - F3o2/chin1_s * chix(i,j,k) * chix(i,j,k) ) + &
                     gupyy_s * ( fyy_s - F3o2/chin1_s * chiy(i,j,k) * chiy(i,j,k) ) + &
                     gupzz_s * ( fzz_s - F3o2/chin1_s * chiz(i,j,k) * chiz(i,j,k) ) + &
               TWO * gupxy_s * ( fxy_s - F3o2/chin1_s * chix(i,j,k) * chiy(i,j,k) ) + &
               TWO * gupxz_s * ( fxz_s - F3o2/chin1_s * chix(i,j,k) * chiz(i,j,k) ) + &
               TWO * gupyz_s * ( fyz_s - F3o2/chin1_s * chiy(i,j,k) * chiz(i,j,k) )
! Add chi part to Ricci tensor:

        Rxx_s = Rxx_s + (fxx_s - chix(i,j,k)*chix(i,j,k)/chin1_s/TWO + gxx_s * f_s)/chin1_s/TWO
        Ryy_s = Ryy_s + (fyy_s - chiy(i,j,k)*chiy(i,j,k)/chin1_s/TWO + gyy_s * f_s)/chin1_s/TWO
        Rzz_s = Rzz_s + (fzz_s - chiz(i,j,k)*chiz(i,j,k)/chin1_s/TWO + gzz_s * f_s)/chin1_s/TWO
        Rxy_s = Rxy_s + (fxy_s - chix(i,j,k)*chiy(i,j,k)/chin1_s/TWO + gxy(i,j,k) * f_s)/chin1_s/TWO
        Rxz_s = Rxz_s + (fxz_s - chix(i,j,k)*chiz(i,j,k)/chin1_s/TWO + gxz(i,j,k) * f_s)/chin1_s/TWO
        Ryz_s = Ryz_s + (fyz_s - chiy(i,j,k)*chiz(i,j,k)/chin1_s/TWO + gyz(i,j,k) * f_s)/chin1_s/TWO

        gxxx2_s = (gupxx_s * chix(i,j,k) + gupxy_s * chiy(i,j,k) + gupxz_s * chiz(i,j,k))/chin1_s
        gxxy2_s = (gupxy_s * chix(i,j,k) + gupyy_s * chiy(i,j,k) + gupyz_s * chiz(i,j,k))/chin1_s
        gxxz2_s = (gupxz_s * chix(i,j,k) + gupyz_s * chiy(i,j,k) + gupzz_s * chiz(i,j,k))/chin1_s
! now get physical second kind of connection
        Gamxxx_s = Gamxxx_s - ( (chix(i,j,k) + chix(i,j,k))/chin1_s - gxx_s * gxxx2_s )*HALF
        Gamyxx_s = Gamyxx_s - (                     - gxx_s * gxxy2_s )*HALF
        Gamzxx_s = Gamzxx_s - (                     - gxx_s * gxxz2_s )*HALF
        Gamxyy_s = Gamxyy_s - (                     - gyy_s * gxxx2_s )*HALF
        Gamyyy_s = Gamyyy_s - ( (chiy(i,j,k) + chiy(i,j,k))/chin1_s - gyy_s * gxxy2_s )*HALF
        Gamzyy_s = Gamzyy_s - (                     - gyy_s * gxxz2_s )*HALF
        Gamxzz_s = Gamxzz_s - (                     - gzz_s * gxxx2_s )*HALF
        Gamyzz_s = Gamyzz_s - (                     - gzz_s * gxxy2_s )*HALF
        Gamzzz_s = Gamzzz_s - ( (chiz(i,j,k) + chiz(i,j,k))/chin1_s - gzz_s * gxxz2_s )*HALF
        Gamxxy_s = Gamxxy_s - (  chiy(i,j,k)        /chin1_s - gxy(i,j,k) * gxxx2_s )*HALF
        Gamyxy_s = Gamyxy_s - (         chix(i,j,k) /chin1_s - gxy(i,j,k) * gxxy2_s )*HALF
        Gamzxy_s = Gamzxy_s - (                     - gxy(i,j,k) * gxxz2_s )*HALF
        Gamxxz_s = Gamxxz_s - (  chiz(i,j,k)        /chin1_s - gxz(i,j,k) * gxxx2_s )*HALF
        Gamyxz_s = Gamyxz_s - (                     - gxz(i,j,k) * gxxy2_s )*HALF
        Gamzxz_s = Gamzxz_s - (         chix(i,j,k) /chin1_s - gxz(i,j,k) * gxxz2_s )*HALF
        Gamxyz_s = Gamxyz_s - (                     - gyz(i,j,k) * gxxx2_s )*HALF
        Gamyyz_s = Gamyyz_s - (  chiz(i,j,k)        /chin1_s - gyz(i,j,k) * gxxy2_s )*HALF
        Gamzyz_s = Gamzyz_s - (         chiy(i,j,k) /chin1_s - gyz(i,j,k) * gxxz2_s )*HALF

#ifdef AMSS_RHS_POINTWISE
        if (pointwise_mode) then
          fxx_s = pw_fxx_lap - Gamxxx_s*Lapx(i,j,k) - Gamyxx_s*Lapy(i,j,k) - Gamzxx_s*Lapz(i,j,k)
          fyy_s = pw_fyy_lap - Gamxyy_s*Lapx(i,j,k) - Gamyyy_s*Lapy(i,j,k) - Gamzyy_s*Lapz(i,j,k)
          fzz_s = pw_fzz_lap - Gamxzz_s*Lapx(i,j,k) - Gamyzz_s*Lapy(i,j,k) - Gamzzz_s*Lapz(i,j,k)
          fxy_s = pw_fxy_lap - Gamxxy_s*Lapx(i,j,k) - Gamyxy_s*Lapy(i,j,k) - Gamzxy_s*Lapz(i,j,k)
          fxz_s = pw_fxz_lap - Gamxxz_s*Lapx(i,j,k) - Gamyxz_s*Lapy(i,j,k) - Gamzxz_s*Lapz(i,j,k)
          fyz_s = pw_fyz_lap - Gamxyz_s*Lapx(i,j,k) - Gamyyz_s*Lapy(i,j,k) - Gamzyz_s*Lapz(i,j,k)
#ifdef AMSS_RHS_POINTWISE
        else
#endif
          fxx_s = fxx_lap(i,j,k) - Gamxxx_s*Lapx(i,j,k) - Gamyxx_s*Lapy(i,j,k) - Gamzxx_s*Lapz(i,j,k)
          fyy_s = fyy_lap(i,j,k) - Gamxyy_s*Lapx(i,j,k) - Gamyyy_s*Lapy(i,j,k) - Gamzyy_s*Lapz(i,j,k)
          fzz_s = fzz_lap(i,j,k) - Gamxzz_s*Lapx(i,j,k) - Gamyzz_s*Lapy(i,j,k) - Gamzzz_s*Lapz(i,j,k)
          fxy_s = fxy_lap(i,j,k) - Gamxxy_s*Lapx(i,j,k) - Gamyxy_s*Lapy(i,j,k) - Gamzxy_s*Lapz(i,j,k)
          fxz_s = fxz_lap(i,j,k) - Gamxxz_s*Lapx(i,j,k) - Gamyxz_s*Lapy(i,j,k) - Gamzxz_s*Lapz(i,j,k)
          fyz_s = fyz_lap(i,j,k) - Gamxyz_s*Lapx(i,j,k) - Gamyyz_s*Lapy(i,j,k) - Gamzyz_s*Lapz(i,j,k)
#ifdef AMSS_RHS_POINTWISE
        end if
#endif
#endif

! store D^i D_i Lap in trK_rhs upto chi
        trK_rhs_s =    gupxx_s * fxx_s + gupyy_s * fyy_s + gupzz_s * fzz_s + &
              TWO* ( gupxy_s * fxy_s + gupxz_s * fxz_s + gupyz_s * fyz_s )
#if 1
!! follow bam code
        S_s =  chin1_s * ( gupxx_s * Sxx(i,j,k) + gupyy_s * Syy(i,j,k) + gupzz_s * Szz(i,j,k) + &
           TWO * ( gupxy_s * Sxy(i,j,k) + gupxz_s * Sxz(i,j,k) + gupyz_s * Syz(i,j,k) ) )
        f_s = F2o3 * trK(i,j,k) * trK(i,j,k) -(&
             gupxx_s * ( &
             gupxx_s * Axx(i,j,k) * Axx(i,j,k) + gupyy_s * Axy(i,j,k) * Axy(i,j,k) + gupzz_s * Axz(i,j,k) * Axz(i,j,k) + &
             TWO * (gupxy_s * Axx(i,j,k) * Axy(i,j,k) + gupxz_s * Axx(i,j,k) * Axz(i,j,k) + gupyz_s * Axy(i,j,k) * Axz(i,j,k)) ) + &
             gupyy_s * ( &
             gupxx_s * Axy(i,j,k) * Axy(i,j,k) + gupyy_s * Ayy(i,j,k) * Ayy(i,j,k) + gupzz_s * Ayz(i,j,k) * Ayz(i,j,k) + &
             TWO * (gupxy_s * Axy(i,j,k) * Ayy(i,j,k) + gupxz_s * Axy(i,j,k) * Ayz(i,j,k) + gupyz_s * Ayy(i,j,k) * Ayz(i,j,k)) ) + &
             gupzz_s * ( &
             gupxx_s * Axz(i,j,k) * Axz(i,j,k) + gupyy_s * Ayz(i,j,k) * Ayz(i,j,k) + gupzz_s * Azz(i,j,k) * Azz(i,j,k) + &
             TWO * (gupxy_s * Axz(i,j,k) * Ayz(i,j,k) + gupxz_s * Axz(i,j,k) * Azz(i,j,k) + gupyz_s * Ayz(i,j,k) * Azz(i,j,k)) ) + &
             TWO * ( &
             gupxy_s * ( &
             gupxx_s * Axx(i,j,k) * Axy(i,j,k) + gupyy_s * Axy(i,j,k) * Ayy(i,j,k) + gupzz_s * Axz(i,j,k) * Ayz(i,j,k) + &
             gupxy_s * (Axx(i,j,k) * Ayy(i,j,k) + Axy(i,j,k) * Axy(i,j,k)) + &
             gupxz_s * (Axx(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Axy(i,j,k)) + &
             gupyz_s * (Axy(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Ayy(i,j,k)) ) + &
             gupxz_s * ( &
             gupxx_s * Axx(i,j,k) * Axz(i,j,k) + gupyy_s * Axy(i,j,k) * Ayz(i,j,k) + gupzz_s * Axz(i,j,k) * Azz(i,j,k) + &
             gupxy_s * (Axx(i,j,k) * Ayz(i,j,k) + Axy(i,j,k) * Axz(i,j,k)) + &
             gupxz_s * (Axx(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Axz(i,j,k)) + &
             gupyz_s * (Axy(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Ayz(i,j,k)) ) + &
             gupyz_s * ( &
             gupxx_s * Axy(i,j,k) * Axz(i,j,k) + gupyy_s * Ayy(i,j,k) * Ayz(i,j,k) + gupzz_s * Ayz(i,j,k) * Azz(i,j,k) + &
             gupxy_s * (Axy(i,j,k) * Ayz(i,j,k) + Ayy(i,j,k) * Axz(i,j,k)) + &
             gupxz_s * (Axy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Axz(i,j,k)) + &
             gupyz_s * (Ayy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Ayz(i,j,k)) ) )) -1.6d1*PI*rho(i,j,k) + EIGHT * PI * S_s
        f_s = - F1o3 *(  gupxx_s * fxx_s + gupyy_s * fyy_s + gupzz_s * fzz_s + &
              TWO* ( gupxy_s * fxy_s + gupxz_s * fxz_s + gupyz_s * fyz_s ) + alpn1_s/chin1_s*f_s)

        fxx_s = alpn1_s * (Rxx_s - EIGHT * PI * Sxx(i,j,k)) - fxx_s
        fxy_s = alpn1_s * (Rxy_s - EIGHT * PI * Sxy(i,j,k)) - fxy_s
        fxz_s = alpn1_s * (Rxz_s - EIGHT * PI * Sxz(i,j,k)) - fxz_s
        fyy_s = alpn1_s * (Ryy_s - EIGHT * PI * Syy(i,j,k)) - fyy_s
        fyz_s = alpn1_s * (Ryz_s - EIGHT * PI * Syz(i,j,k)) - fyz_s
        fzz_s = alpn1_s * (Rzz_s - EIGHT * PI * Szz(i,j,k)) - fzz_s
#else
! Add lapse and S_ij parts to Ricci tensor:

        fxx_s = alpn1_s * (Rxx_s - EIGHT * PI * Sxx(i,j,k)) - fxx_s
        fxy_s = alpn1_s * (Rxy_s - EIGHT * PI * Sxy(i,j,k)) - fxy_s
        fxz_s = alpn1_s * (Rxz_s - EIGHT * PI * Sxz(i,j,k)) - fxz_s
        fyy_s = alpn1_s * (Ryy_s - EIGHT * PI * Syy(i,j,k)) - fyy_s
        fyz_s = alpn1_s * (Ryz_s - EIGHT * PI * Syz(i,j,k)) - fyz_s
        fzz_s = alpn1_s * (Rzz_s - EIGHT * PI * Szz(i,j,k)) - fzz_s

! Compute trace-free part (note: chi^-1 and chi cancel!):

        f_s = F1o3 *(  gupxx_s * fxx_s + gupyy_s * fyy_s + gupzz_s * fzz_s + &
              TWO* ( gupxy_s * fxy_s + gupxz_s * fxz_s + gupyz_s * fyz_s ) )
#endif

        Axx_rhs_s = fxx_s - gxx_s * f_s
        Ayy_rhs_s = fyy_s - gyy_s * f_s
        Azz_rhs_s = fzz_s - gzz_s * f_s
        Axy_rhs_s = fxy_s - gxy(i,j,k) * f_s
        Axz_rhs_s = fxz_s - gxz(i,j,k) * f_s
        Ayz_rhs_s = fyz_s - gyz(i,j,k) * f_s

! Now: store A_il A^l_j into fij:

        fxx_s =       gupxx_s * Axx(i,j,k) * Axx(i,j,k) + gupyy_s * Axy(i,j,k) * Axy(i,j,k) + gupzz_s * Axz(i,j,k) * Axz(i,j,k) + &
             TWO * (gupxy_s * Axx(i,j,k) * Axy(i,j,k) + gupxz_s * Axx(i,j,k) * Axz(i,j,k) + gupyz_s * Axy(i,j,k) * Axz(i,j,k))
        fyy_s =       gupxx_s * Axy(i,j,k) * Axy(i,j,k) + gupyy_s * Ayy(i,j,k) * Ayy(i,j,k) + gupzz_s * Ayz(i,j,k) * Ayz(i,j,k) + &
             TWO * (gupxy_s * Axy(i,j,k) * Ayy(i,j,k) + gupxz_s * Axy(i,j,k) * Ayz(i,j,k) + gupyz_s * Ayy(i,j,k) * Ayz(i,j,k))
        fzz_s =       gupxx_s * Axz(i,j,k) * Axz(i,j,k) + gupyy_s * Ayz(i,j,k) * Ayz(i,j,k) + gupzz_s * Azz(i,j,k) * Azz(i,j,k) + &
             TWO * (gupxy_s * Axz(i,j,k) * Ayz(i,j,k) + gupxz_s * Axz(i,j,k) * Azz(i,j,k) + gupyz_s * Ayz(i,j,k) * Azz(i,j,k))
        fxy_s =       gupxx_s * Axx(i,j,k) * Axy(i,j,k) + gupyy_s * Axy(i,j,k) * Ayy(i,j,k) + gupzz_s * Axz(i,j,k) * Ayz(i,j,k) + &
                    gupxy_s *(Axx(i,j,k) * Ayy(i,j,k) + Axy(i,j,k) * Axy(i,j,k))                            + &
                    gupxz_s *(Axx(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Axy(i,j,k))                            + &
                    gupyz_s *(Axy(i,j,k) * Ayz(i,j,k) + Axz(i,j,k) * Ayy(i,j,k))
        fxz_s =       gupxx_s * Axx(i,j,k) * Axz(i,j,k) + gupyy_s * Axy(i,j,k) * Ayz(i,j,k) + gupzz_s * Axz(i,j,k) * Azz(i,j,k) + &
                    gupxy_s *(Axx(i,j,k) * Ayz(i,j,k) + Axy(i,j,k) * Axz(i,j,k))                            + &
                    gupxz_s *(Axx(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Axz(i,j,k))                            + &
                    gupyz_s *(Axy(i,j,k) * Azz(i,j,k) + Axz(i,j,k) * Ayz(i,j,k))
        fyz_s =       gupxx_s * Axy(i,j,k) * Axz(i,j,k) + gupyy_s * Ayy(i,j,k) * Ayz(i,j,k) + gupzz_s * Ayz(i,j,k) * Azz(i,j,k) + &
                    gupxy_s *(Axy(i,j,k) * Ayz(i,j,k) + Ayy(i,j,k) * Axz(i,j,k))                            + &
                    gupxz_s *(Axy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Axz(i,j,k))                            + &
                    gupyz_s *(Ayy(i,j,k) * Azz(i,j,k) + Ayz(i,j,k) * Ayz(i,j,k))

        f_s = chin1_s
! store D^i D_i Lap in trK_rhs
        trK_rhs_s = f_s*trK_rhs_s

        Axx_rhs_s =           f_s * Axx_rhs_s+ alpn1_s * (trK(i,j,k) * Axx(i,j,k) - TWO * fxx_s)  + &
             TWO * (  Axx(i,j,k) * betaxx(i,j,k) +   Axy(i,j,k) * betayx(i,j,k) +   Axz(i,j,k) * betazx(i,j,k) )- &
               F2o3 * Axx(i,j,k) * div_beta_s

        Ayy_rhs_s =           f_s * Ayy_rhs_s+ alpn1_s * (trK(i,j,k) * Ayy(i,j,k) - TWO * fyy_s)  + &
             TWO * (  Axy(i,j,k) * betaxy(i,j,k) +   Ayy(i,j,k) * betayy(i,j,k) +   Ayz(i,j,k) * betazy(i,j,k) )- &
               F2o3 * Ayy(i,j,k) * div_beta_s

        Azz_rhs_s =           f_s * Azz_rhs_s+ alpn1_s * (trK(i,j,k) * Azz(i,j,k) - TWO * fzz_s)  + &
             TWO * (  Axz(i,j,k) * betaxz(i,j,k) +   Ayz(i,j,k) * betayz(i,j,k) +   Azz(i,j,k) * betazz(i,j,k) )- &
               F2o3 * Azz(i,j,k) * div_beta_s

        Axy_rhs_s =           f_s * Axy_rhs_s+ alpn1_s *( trK(i,j,k) * Axy(i,j,k)  - TWO * fxy_s )+ &
                      Axx(i,j,k) * betaxy(i,j,k)                  +   Axz(i,j,k) * betazy(i,j,k)  + &
                                       Ayy(i,j,k) * betayx(i,j,k) +   Ayz(i,j,k) * betazx(i,j,k)  + &
               F1o3 * Axy(i,j,k) * div_beta_s                -   Axy(i,j,k) * betazz(i,j,k)

        Ayz_rhs_s =           f_s * Ayz_rhs_s+ alpn1_s *( trK(i,j,k) * Ayz(i,j,k)  - TWO * fyz_s )+ &
                      Axy(i,j,k) * betaxz(i,j,k) +   Ayy(i,j,k) * betayz(i,j,k)                   + &
                      Axz(i,j,k) * betaxy(i,j,k)                  +   Azz(i,j,k) * betazy(i,j,k)  + &
               F1o3 * Ayz(i,j,k) * div_beta_s                -   Ayz(i,j,k) * betaxx(i,j,k)

        Axz_rhs_s =           f_s * Axz_rhs_s+ alpn1_s *( trK(i,j,k) * Axz(i,j,k)  - TWO * fxz_s )+ &
                      Axx(i,j,k) * betaxz(i,j,k) +   Axy(i,j,k) * betayz(i,j,k)                   + &
                                       Ayz(i,j,k) * betayx(i,j,k) +   Azz(i,j,k) * betazx(i,j,k)  + &
               F1o3 * Axz(i,j,k) * div_beta_s                -   Axz(i,j,k) * betayy(i,j,k)      !rhs for Aij

! Compute trace of S_ij

        S_s =  f_s * ( gupxx_s * Sxx(i,j,k) + gupyy_s * Syy(i,j,k) + gupzz_s * Szz(i,j,k) + &
           TWO * ( gupxy_s * Sxy(i,j,k) + gupxz_s * Sxz(i,j,k) + gupyz_s * Syz(i,j,k) ) )

        trK_rhs_s = - trK_rhs_s + alpn1_s *( F1o3 * trK(i,j,k) * trK(i,j,k)         + &
                      gupxx_s * fxx_s + gupyy_s * fyy_s + gupzz_s * fzz_s   + &
              TWO * ( gupxy_s * fxy_s + gupxz_s * fxz_s + gupyz_s * fyz_s ) + &
             FOUR * PI * ( rho(i,j,k) + S_s ))                                !rhs for trK

!!!! gauge variable part

#ifdef AMSS_RHS_REORDER
        ! Phase 1 already stored Lap_rhs / beta_rhs; Phase 4 already stored
        ! Gamx_rhs / Gamy_rhs / Gamz_rhs / dtSf_rhs.  Only the unconditional
        ! store-back of these to the arrays is skipped under REORDER -- the
        ! scalar values were written to the arrays earlier, so skip here.
#else
        Lap_rhs(i,j,k) = -TWO*alpn1_s*trK(i,j,k)
        betax_rhs(i,j,k) = FF*dtSfx(i,j,k)
        betay_rhs(i,j,k) = FF*dtSfy(i,j,k)
        betaz_rhs(i,j,k) = FF*dtSfz(i,j,k)

        dtSfx_rhs(i,j,k) = Gamx_rhs_s - eta*dtSfx(i,j,k)
        dtSfy_rhs(i,j,k) = Gamy_rhs_s - eta*dtSfy(i,j,k)
        dtSfz_rhs(i,j,k) = Gamz_rhs_s - eta*dtSfz(i,j,k)
#endif

!~~~~~~> store loop-final values to arrays needed after the loop

        trK_rhs(i,j,k) = trK_rhs_s
#ifdef AMSS_RHS_REORDER
        ! Gamx_rhs/Gamy_rhs/Gamz_rhs already stored in Phase 4 above.
#else
        Gamx_rhs(i,j,k) = Gamx_rhs_s
        Gamy_rhs(i,j,k) = Gamy_rhs_s
        Gamz_rhs(i,j,k) = Gamz_rhs_s
#endif
        Axx_rhs(i,j,k) = Axx_rhs_s
        Axy_rhs(i,j,k) = Axy_rhs_s
        Axz_rhs(i,j,k) = Axz_rhs_s
        Ayy_rhs(i,j,k) = Ayy_rhs_s
        Ayz_rhs(i,j,k) = Ayz_rhs_s
        Azz_rhs(i,j,k) = Azz_rhs_s
#ifdef AMSS_RHS_SKIP_CONSTRAINT_STORES
        ! These arrays are consumed only by the co=0 constraint block below.
        ! The scalar values remain live for the evolution RHS in co=1.
        if (co == 0) then
#endif
        gupxx(i,j,k) = gupxx_s
        gupxy(i,j,k) = gupxy_s
        gupxz(i,j,k) = gupxz_s
        gupyy(i,j,k) = gupyy_s
        gupyz(i,j,k) = gupyz_s
        gupzz(i,j,k) = gupzz_s
        Rxx(i,j,k) = Rxx_s
        Rxy(i,j,k) = Rxy_s
        Rxz(i,j,k) = Rxz_s
        Ryy(i,j,k) = Ryy_s
        Ryz(i,j,k) = Ryz_s
        Rzz(i,j,k) = Rzz_s
        Gamxxx(i,j,k) = Gamxxx_s
        Gamxxy(i,j,k) = Gamxxy_s
        Gamxxz(i,j,k) = Gamxxz_s
        Gamxyy(i,j,k) = Gamxyy_s
        Gamxyz(i,j,k) = Gamxyz_s
        Gamxzz(i,j,k) = Gamxzz_s
        Gamyxx(i,j,k) = Gamyxx_s
        Gamyxy(i,j,k) = Gamyxy_s
        Gamyxz(i,j,k) = Gamyxz_s
        Gamyyy(i,j,k) = Gamyyy_s
        Gamyyz(i,j,k) = Gamyyz_s
        Gamyzz(i,j,k) = Gamyzz_s
        Gamzxx(i,j,k) = Gamzxx_s
        Gamzxy(i,j,k) = Gamzxy_s
        Gamzxz(i,j,k) = Gamzxz_s
        Gamzyy(i,j,k) = Gamzyy_s
        Gamzyz(i,j,k) = Gamzyz_s
        Gamzzz(i,j,k) = Gamzzz_s
#ifdef AMSS_RHS_SKIP_CONSTRAINT_STORES
        end if
#endif

      end do
    end do
#ifdef AMSS_RHS_TILED
      end do
    end do
#endif
  end do
#ifdef AMSS_RHS_OMP_ASSEMBLY
!$omp end parallel do
#endif

#ifdef AMSS_RHS_POINTWISE
#undef betaxx
#undef betaxy
#undef betaxz
#undef betayx
#undef betayy
#undef betayz
#undef betazx
#undef betazy
#undef betazz
#undef chix
#undef chiy
#undef chiz
#undef gxxx
#undef gxxy
#undef gxxz
#undef gxyx
#undef gxyy
#undef gxyz
#undef gxzx
#undef gxzy
#undef gxzz
#undef gyyx
#undef gyyy
#undef gyyz
#undef gyzx
#undef gyzy
#undef gyzz
#undef gzzx
#undef gzzy
#undef gzzz
#undef Lapx
#undef Lapy
#undef Lapz
#undef Kx
#undef Ky
#undef Kz
#undef Gamxx
#undef Gamxy
#undef Gamxz
#undef Gamyx
#undef Gamyy
#undef Gamyz
#undef Gamzx
#undef Gamzy
#undef Gamzz
#endif

!!!!!!!!!advection term part

  SSS(1)=SYM
  SSS(2)=SYM
  SSS(3)=SYM

  AAS(1)=ANTI
  AAS(2)=ANTI
  AAS(3)=SYM

  ASA(1)=ANTI
  ASA(2)=SYM
  ASA(3)=ANTI

  SAA(1)=SYM
  SAA(2)=ANTI
  SAA(3)=ANTI

  ASS(1)=ANTI
  ASS(2)=SYM
  ASS(3)=SYM

  SAS(1)=SYM
  SAS(2)=ANTI
  SAS(3)=SYM

  SSA(1)=SYM
  SSA(2)=SYM
  SSA(3)=ANTI

#ifdef AMSS_BATCH_STENCIL
  if (bs_enabled() .and. Symmetry <= 1) then
    block
      type(field_ptr) :: lfp(24), lfrp(24)
#ifdef AMSS_FUSED_RHS_TAIL
      type(field_ptr) :: kfp(24)
#endif
      real*8 :: ls3(24)
      integer :: lnv
      lnv = 0
#ifdef AMSS_RHS_RAW_DIAG_LOPSIDED
      if (raw_diag_lopsided) then
        ! D(gxx)=D(dxx+1)=D(dxx); the lopsided coefficients sum to zero.
        ! This removes three full-block auxiliary fields from the hot tail.
        lnv = lnv+1; lfp(lnv)%p => dxx;  lfrp(lnv)%p => gxx_rhs;  ls3(lnv) = SYM
      else
        lnv = lnv+1; lfp(lnv)%p => gxx;  lfrp(lnv)%p => gxx_rhs;  ls3(lnv) = SYM
      end if
#else
      lnv = lnv+1; lfp(lnv)%p => gxx;  lfrp(lnv)%p => gxx_rhs;  ls3(lnv) = SYM
#endif
      lnv = lnv+1; lfp(lnv)%p => gxy;  lfrp(lnv)%p => gxy_rhs;  ls3(lnv) = SYM
      lnv = lnv+1; lfp(lnv)%p => gxz;  lfrp(lnv)%p => gxz_rhs;  ls3(lnv) = ANTI
#ifdef AMSS_RHS_RAW_DIAG_LOPSIDED
      if (raw_diag_lopsided) then
        lnv = lnv+1; lfp(lnv)%p => dyy;  lfrp(lnv)%p => gyy_rhs;  ls3(lnv) = SYM
      else
        lnv = lnv+1; lfp(lnv)%p => gyy;  lfrp(lnv)%p => gyy_rhs;  ls3(lnv) = SYM
      end if
#else
      lnv = lnv+1; lfp(lnv)%p => gyy;  lfrp(lnv)%p => gyy_rhs;  ls3(lnv) = SYM
#endif
      lnv = lnv+1; lfp(lnv)%p => gyz;  lfrp(lnv)%p => gyz_rhs;  ls3(lnv) = ANTI
#ifdef AMSS_RHS_RAW_DIAG_LOPSIDED
      if (raw_diag_lopsided) then
        lnv = lnv+1; lfp(lnv)%p => dzz;  lfrp(lnv)%p => gzz_rhs;  ls3(lnv) = SYM
      else
        lnv = lnv+1; lfp(lnv)%p => gzz;  lfrp(lnv)%p => gzz_rhs;  ls3(lnv) = SYM
      end if
#else
      lnv = lnv+1; lfp(lnv)%p => gzz;  lfrp(lnv)%p => gzz_rhs;  ls3(lnv) = SYM
#endif
      lnv = lnv+1; lfp(lnv)%p => Axx;  lfrp(lnv)%p => Axx_rhs;  ls3(lnv) = SYM
      lnv = lnv+1; lfp(lnv)%p => Axy;  lfrp(lnv)%p => Axy_rhs;  ls3(lnv) = SYM
      lnv = lnv+1; lfp(lnv)%p => Axz;  lfrp(lnv)%p => Axz_rhs;  ls3(lnv) = ANTI
      lnv = lnv+1; lfp(lnv)%p => Ayy;  lfrp(lnv)%p => Ayy_rhs;  ls3(lnv) = SYM
      lnv = lnv+1; lfp(lnv)%p => Ayz;  lfrp(lnv)%p => Ayz_rhs;  ls3(lnv) = ANTI
      lnv = lnv+1; lfp(lnv)%p => Azz;  lfrp(lnv)%p => Azz_rhs;  ls3(lnv) = SYM
      lnv = lnv+1; lfp(lnv)%p => chi;  lfrp(lnv)%p => chi_rhs;  ls3(lnv) = SYM
      lnv = lnv+1; lfp(lnv)%p => trK;  lfrp(lnv)%p => trK_rhs;  ls3(lnv) = SYM
      lnv = lnv+1; lfp(lnv)%p => Gamx; lfrp(lnv)%p => Gamx_rhs; ls3(lnv) = SYM
      lnv = lnv+1; lfp(lnv)%p => Gamy; lfrp(lnv)%p => Gamy_rhs; ls3(lnv) = SYM
      lnv = lnv+1; lfp(lnv)%p => Gamz; lfrp(lnv)%p => Gamz_rhs; ls3(lnv) = ANTI
      lnv = lnv+1; lfp(lnv)%p => Lap;   lfrp(lnv)%p => Lap_rhs;   ls3(lnv) = SYM
      lnv = lnv+1; lfp(lnv)%p => betax; lfrp(lnv)%p => betax_rhs; ls3(lnv) = SYM
      lnv = lnv+1; lfp(lnv)%p => betay; lfrp(lnv)%p => betay_rhs; ls3(lnv) = SYM
      lnv = lnv+1; lfp(lnv)%p => betaz; lfrp(lnv)%p => betaz_rhs; ls3(lnv) = ANTI
      lnv = lnv+1; lfp(lnv)%p => dtSfx; lfrp(lnv)%p => dtSfx_rhs; ls3(lnv) = SYM
      lnv = lnv+1; lfp(lnv)%p => dtSfy; lfrp(lnv)%p => dtSfy_rhs; ls3(lnv) = SYM
      lnv = lnv+1; lfp(lnv)%p => dtSfz; lfrp(lnv)%p => dtSfz_rhs; ls3(lnv) = ANTI
#ifdef AMSS_FUSED_RHS_TAIL
      kfp(1)%p => dxx;   kfp(2)%p => gxy;   kfp(3)%p => gxz
      kfp(4)%p => dyy;   kfp(5)%p => gyz;   kfp(6)%p => dzz
      kfp(7)%p => Axx;   kfp(8)%p => Axy;   kfp(9)%p => Axz
      kfp(10)%p => Ayy;  kfp(11)%p => Ayz;  kfp(12)%p => Azz
      kfp(13)%p => chi;  kfp(14)%p => trK;  kfp(15)%p => Gamx
      kfp(16)%p => Gamy; kfp(17)%p => Gamz; kfp(18)%p => Lap
      kfp(19)%p => betax;kfp(20)%p => betay;kfp(21)%p => betaz
      kfp(22)%p => dtSfx;kfp(23)%p => dtSfy;kfp(24)%p => dtSfz
      if (bs_size() == 2) then
        call lopsided_batch(ex,X,Y,Z,lnv,lfp,lfrp,ls3,betax,betay,betaz,Symmetry, &
                            kfp,eps,.true.)
      else
        call lopsided_batch(ex,X,Y,Z,lnv,lfp,lfrp,ls3,betax,betay,betaz,Symmetry)
      end if
#else
      call lopsided_batch(ex,X,Y,Z,lnv,lfp,lfrp,ls3,betax,betay,betaz,Symmetry)
#endif
    end block
  else
  call lopsided(ex,X,Y,Z,gxx,gxx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gxy,gxy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,gxz,gxz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,gyy,gyy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gyz,gyz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,gzz,gzz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Axx,Axx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Axy,Axy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,Axz,Axz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,Ayy,Ayy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Ayz,Ayz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,Azz,Azz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,chi,chi_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,trK,trK_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Gamx,Gamx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,Gamy,Gamy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,Gamz,Gamz_rhs,betax,betay,betaz,Symmetry,SSA)
!!
  call lopsided(ex,X,Y,Z,Lap,Lap_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,betax,betax_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,betay,betay_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,betaz,betaz_rhs,betax,betay,betaz,Symmetry,SSA)

  call lopsided(ex,X,Y,Z,dtSfx,dtSfx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,dtSfy,dtSfy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,dtSfz,dtSfz_rhs,betax,betay,betaz,Symmetry,SSA)
  end if
#else
  call lopsided(ex,X,Y,Z,gxx,gxx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gxy,gxy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,gxz,gxz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,gyy,gyy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,gyz,gyz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,gzz,gzz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Axx,Axx_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Axy,Axy_rhs,betax,betay,betaz,Symmetry,AAS)
  call lopsided(ex,X,Y,Z,Axz,Axz_rhs,betax,betay,betaz,Symmetry,ASA)
  call lopsided(ex,X,Y,Z,Ayy,Ayy_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,Ayz,Ayz_rhs,betax,betay,betaz,Symmetry,SAA)
  call lopsided(ex,X,Y,Z,Azz,Azz_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,chi,chi_rhs,betax,betay,betaz,Symmetry,SSS)
  call lopsided(ex,X,Y,Z,trK,trK_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,Gamx,Gamx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,Gamy,Gamy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,Gamz,Gamz_rhs,betax,betay,betaz,Symmetry,SSA)
!!
  call lopsided(ex,X,Y,Z,Lap,Lap_rhs,betax,betay,betaz,Symmetry,SSS)

  call lopsided(ex,X,Y,Z,betax,betax_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,betay,betay_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,betaz,betaz_rhs,betax,betay,betaz,Symmetry,SSA)

  call lopsided(ex,X,Y,Z,dtSfx,dtSfx_rhs,betax,betay,betaz,Symmetry,ASS)
  call lopsided(ex,X,Y,Z,dtSfy,dtSfy_rhs,betax,betay,betaz,Symmetry,SAS)
  call lopsided(ex,X,Y,Z,dtSfz,dtSfz_rhs,betax,betay,betaz,Symmetry,SSA)
#endif

  if(eps>0)then
! usual Kreiss-Oliger dissipation
#ifdef AMSS_BATCH_STENCIL
#ifdef AMSS_FUSED_RHS_TAIL
  if (bs_enabled() .and. Symmetry <= 1 .and. bs_size() == 2) then
    ! B=2 lopsided_batch already applied KO in the same RHS scan.
    continue
  else if (bs_enabled() .and. Symmetry <= 1) then
#else
  if (bs_enabled() .and. Symmetry <= 1) then
#endif
    block
      type(field_ptr) :: kfp(24), kfrp(24)
      real*8 :: ks3(24)
      integer :: knv
      knv = 0
      knv = knv+1; kfp(knv)%p => chi;  kfrp(knv)%p => chi_rhs;  ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => trK;  kfrp(knv)%p => trK_rhs;  ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => dxx;  kfrp(knv)%p => gxx_rhs; ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => gxy;  kfrp(knv)%p => gxy_rhs; ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => gxz;  kfrp(knv)%p => gxz_rhs; ks3(knv) = ANTI
      knv = knv+1; kfp(knv)%p => dyy;  kfrp(knv)%p => gyy_rhs; ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => gyz;  kfrp(knv)%p => gyz_rhs; ks3(knv) = ANTI
      knv = knv+1; kfp(knv)%p => dzz;  kfrp(knv)%p => gzz_rhs; ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => Axx;  kfrp(knv)%p => Axx_rhs; ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => Axy;  kfrp(knv)%p => Axy_rhs; ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => Axz;  kfrp(knv)%p => Axz_rhs; ks3(knv) = ANTI
      knv = knv+1; kfp(knv)%p => Ayy;  kfrp(knv)%p => Ayy_rhs; ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => Ayz;  kfrp(knv)%p => Ayz_rhs; ks3(knv) = ANTI
      knv = knv+1; kfp(knv)%p => Azz;  kfrp(knv)%p => Azz_rhs; ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => Gamx; kfrp(knv)%p => Gamx_rhs; ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => Gamy; kfrp(knv)%p => Gamy_rhs; ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => Gamz; kfrp(knv)%p => Gamz_rhs; ks3(knv) = ANTI
#if 1
!! bam does not apply dissipation on gauge variables
      knv = knv+1; kfp(knv)%p => Lap;   kfrp(knv)%p => Lap_rhs;   ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => betax; kfrp(knv)%p => betax_rhs; ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => betay; kfrp(knv)%p => betay_rhs; ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => betaz; kfrp(knv)%p => betaz_rhs; ks3(knv) = ANTI
      knv = knv+1; kfp(knv)%p => dtSfx; kfrp(knv)%p => dtSfx_rhs; ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => dtSfy; kfrp(knv)%p => dtSfy_rhs; ks3(knv) = SYM
      knv = knv+1; kfp(knv)%p => dtSfz; kfrp(knv)%p => dtSfz_rhs; ks3(knv) = ANTI
#endif
      call kodis_batch(ex,X,Y,Z,knv,kfp,kfrp,ks3,Symmetry,eps)
    end block
  else
  call kodis(ex,X,Y,Z,chi,chi_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,trK,trK_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dxx,gxx_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxy,gxy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxz,gxz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dyy,gyy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gyz,gyz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dzz,gzz_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axx,Axx_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axy,Axy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axz,Axz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayy,Ayy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayz,Ayz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Azz,Azz_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamx,Gamx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamy,Gamy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamz,Gamz_rhs,SSA,Symmetry,eps)

#if 1
!! bam does not apply dissipation on gauge variables
  call kodis(ex,X,Y,Z,Lap,Lap_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betax,betax_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betay,betay_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betaz,betaz_rhs,SSA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfx,dtSfx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfy,dtSfy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfz,dtSfz_rhs,SSA,Symmetry,eps)
#endif
  end if
#else
  call kodis(ex,X,Y,Z,chi,chi_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,trK,trK_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dxx,gxx_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxy,gxy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gxz,gxz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dyy,gyy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,gyz,gyz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dzz,gzz_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axx,Axx_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axy,Axy_rhs,AAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Axz,Axz_rhs,ASA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayy,Ayy_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Ayz,Ayz_rhs,SAA,Symmetry,eps)
  call kodis(ex,X,Y,Z,Azz,Azz_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamx,Gamx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamy,Gamy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,Gamz,Gamz_rhs,SSA,Symmetry,eps)

#if 1
!! bam does not apply dissipation on gauge variables
  call kodis(ex,X,Y,Z,Lap,Lap_rhs,SSS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betax,betax_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betay,betay_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,betaz,betaz_rhs,SSA,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfx,dtSfx_rhs,ASS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfy,dtSfy_rhs,SAS,Symmetry,eps)
  call kodis(ex,X,Y,Z,dtSfz,dtSfz_rhs,SSA,Symmetry,eps)
#endif
#endif

  endif

  if(co == 0)then
! ham_Res = trR + 2/3 * K^2 - A_ij * A^ij - 16 * PI * rho
! here trR is respect to physical metric
  ham_Res =   gupxx * Rxx + gupyy * Ryy + gupzz * Rzz + &
        TWO* ( gupxy * Rxy + gupxz * Rxz + gupyz * Ryz )

  ham_Res = chin1*ham_Res + F2o3 * trK * trK -(&
       gupxx * ( &
       gupxx * Axx * Axx + gupyy * Axy * Axy + gupzz * Axz * Axz + &
       TWO * (gupxy * Axx * Axy + gupxz * Axx * Axz + gupyz * Axy * Axz) ) + &
       gupyy * ( &
       gupxx * Axy * Axy + gupyy * Ayy * Ayy + gupzz * Ayz * Ayz + &
       TWO * (gupxy * Axy * Ayy + gupxz * Axy * Ayz + gupyz * Ayy * Ayz) ) + &
       gupzz * ( &
       gupxx * Axz * Axz + gupyy * Ayz * Ayz + gupzz * Azz * Azz + &
       TWO * (gupxy * Axz * Ayz + gupxz * Axz * Azz + gupyz * Ayz * Azz) ) + &
       TWO * ( &
       gupxy * ( &
       gupxx * Axx * Axy + gupyy * Axy * Ayy + gupzz * Axz * Ayz + &
       gupxy * (Axx * Ayy + Axy * Axy) + &
       gupxz * (Axx * Ayz + Axz * Axy) + &
       gupyz * (Axy * Ayz + Axz * Ayy) ) + &
       gupxz * ( &
       gupxx * Axx * Axz + gupyy * Axy * Ayz + gupzz * Axz * Azz + &
       gupxy * (Axx * Ayz + Axy * Axz) + &
       gupxz * (Axx * Azz + Axz * Axz) + &
       gupyz * (Axy * Azz + Axz * Ayz) ) + &
       gupyz * ( &
       gupxx * Axy * Axz + gupyy * Ayy * Ayz + gupzz * Ayz * Azz + &
       gupxy * (Axy * Ayz + Ayy * Axz) + &
       gupxz * (Axy * Azz + Ayz * Axz) + &
       gupyz * (Ayy * Azz + Ayz * Ayz) ) ))- F16 * PI * rho

! mov_Res_j = gupkj*(-1/chi d_k chi*A_ij + D_k A_ij) - 2/3 d_j trK - 8 PI s_j where D respect to physical metric
! store D_i A_jk - 1/chi d_i chi*A_jk in gjk_i
  call fderivs(ex,Axx,gxxx,gxxy,gxxz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Axy,gxyx,gxyy,gxyz,X,Y,Z,ANTI,ANTI,SYM ,Symmetry,0)
  call fderivs(ex,Axz,gxzx,gxzy,gxzz,X,Y,Z,ANTI,SYM ,ANTI,Symmetry,0)
  call fderivs(ex,Ayy,gyyx,gyyy,gyyz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)
  call fderivs(ex,Ayz,gyzx,gyzy,gyzz,X,Y,Z,SYM ,ANTI,ANTI,Symmetry,0)
  call fderivs(ex,Azz,gzzx,gzzy,gzzz,X,Y,Z,SYM ,SYM ,SYM ,Symmetry,0)

  gxxx = gxxx - (  Gamxxx * Axx + Gamyxx * Axy + Gamzxx * Axz &
                 + Gamxxx * Axx + Gamyxx * Axy + Gamzxx * Axz) - chix*Axx/chin1
  gxyx = gxyx - (  Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz &
                 + Gamxxx * Axy + Gamyxx * Ayy + Gamzxx * Ayz) - chix*Axy/chin1
  gxzx = gxzx - (  Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz &
                 + Gamxxx * Axz + Gamyxx * Ayz + Gamzxx * Azz) - chix*Axz/chin1
  gyyx = gyyx - (  Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz &
                 + Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz) - chix*Ayy/chin1
  gyzx = gyzx - (  Gamxxz * Axy + Gamyxz * Ayy + Gamzxz * Ayz &
                 + Gamxxy * Axz + Gamyxy * Ayz + Gamzxy * Azz) - chix*Ayz/chin1
  gzzx = gzzx - (  Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz &
                 + Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz) - chix*Azz/chin1
  gxxy = gxxy - (  Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz &
                 + Gamxxy * Axx + Gamyxy * Axy + Gamzxy * Axz) - chiy*Axx/chin1
  gxyy = gxyy - (  Gamxyy * Axx + Gamyyy * Axy + Gamzyy * Axz &
                 + Gamxxy * Axy + Gamyxy * Ayy + Gamzxy * Ayz) - chiy*Axy/chin1
  gxzy = gxzy - (  Gamxyz * Axx + Gamyyz * Axy + Gamzyz * Axz &
                 + Gamxxy * Axz + Gamyxy * Ayz + Gamzxy * Azz) - chiy*Axz/chin1
  gyyy = gyyy - (  Gamxyy * Axy + Gamyyy * Ayy + Gamzyy * Ayz &
                 + Gamxyy * Axy + Gamyyy * Ayy + Gamzyy * Ayz) - chiy*Ayy/chin1
  gyzy = gyzy - (  Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz &
                 + Gamxyy * Axz + Gamyyy * Ayz + Gamzyy * Azz) - chiy*Ayz/chin1
  gzzy = gzzy - (  Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz &
                 + Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz) - chiy*Azz/chin1
  gxxz = gxxz - (  Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz &
                 + Gamxxz * Axx + Gamyxz * Axy + Gamzxz * Axz) - chiz*Axx/chin1
  gxyz = gxyz - (  Gamxyz * Axx + Gamyyz * Axy + Gamzyz * Axz &
                 + Gamxxz * Axy + Gamyxz * Ayy + Gamzxz * Ayz) - chiz*Axy/chin1
  gxzz = gxzz - (  Gamxzz * Axx + Gamyzz * Axy + Gamzzz * Axz &
                 + Gamxxz * Axz + Gamyxz * Ayz + Gamzxz * Azz) - chiz*Axz/chin1
  gyyz = gyyz - (  Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz &
                 + Gamxyz * Axy + Gamyyz * Ayy + Gamzyz * Ayz) - chiz*Ayy/chin1
  gyzz = gyzz - (  Gamxzz * Axy + Gamyzz * Ayy + Gamzzz * Ayz &
                 + Gamxyz * Axz + Gamyyz * Ayz + Gamzyz * Azz) - chiz*Ayz/chin1
  gzzz = gzzz - (  Gamxzz * Axz + Gamyzz * Ayz + Gamzzz * Azz &
                 + Gamxzz * Axz + Gamyzz * Ayz + Gamzzz * Azz) - chiz*Azz/chin1
movx_Res = gupxx*gxxx + gupyy*gxyy + gupzz*gxzz &
          +gupxy*gxyx + gupxz*gxzx + gupyz*gxzy &
          +gupxy*gxxy + gupxz*gxxz + gupyz*gxyz
movy_Res = gupxx*gxyx + gupyy*gyyy + gupzz*gyzz &
          +gupxy*gyyx + gupxz*gyzx + gupyz*gyzy &
          +gupxy*gxyy + gupxz*gxyz + gupyz*gyyz
movz_Res = gupxx*gxzx + gupyy*gyzy + gupzz*gzzz &
          +gupxy*gyzx + gupxz*gzzx + gupyz*gzzy &
          +gupxy*gxzy + gupxz*gxzz + gupyz*gyzz

movx_Res = movx_Res - F2o3*Kx - F8*PI*sx
movy_Res = movy_Res - F2o3*Ky - F8*PI*sy
movz_Res = movz_Res - F2o3*Kz - F8*PI*sz
  endif


  gont = 0

  return

  end function compute_rhs_bssn_fused
