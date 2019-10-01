#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#ifdef TARGET_NEC_SX
C     make sure that the factorized EOS is used on NEC vector computers
# define USE_FACTORIZED_EOS
#endif
C     this was the default, so we keep it that way
#define USE_FACTORIZED_POLY

      subroutine find_rho_2d_ad( 
     I      imin, imax, jmin, jmax, kref, 
     I      tfld, 
     O      tfld_ad, 
     I      sfld, 
     O      sfld_ad, 
     U      rholoc, 
     U      rholoc_ad, 
     I      k, bi, bj, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | o SUBROUTINE FIND_RHO_2D_AD
C     |   Passes sensitivity from rho_ad to t_ad and s_ad
C     |   This is hand coded so that rho_ad does not get zero'd out
C     |     I hope it works!!!
C     *==========================================================*
C     | kRef - determines pressure reference level
C     |        (not used in 'LINEAR' mode)
C     | Note:  k is not used in the forward code, so I'm eliminating it
C     |         here
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global variables ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "EOS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     k    :: Level of Theta/Salt slice
C     kRef :: Pressure reference level
      INTEGER iMin,iMax,jMin,jMax
      INTEGER k,kRef
      _RL tFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sFld  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL tFld_ad  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL sFld_ad  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoLoc (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoLoc_ad(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      INTEGER bi, bj
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
      INTEGER i,j,ip1,ip2
      _RL refTemp,refSalt,sigRef,tP,sP,deltaSig,dRho
      _RL oneMinusKap, facPres, theta_v
      _RL locPres(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoP0  (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL bulkMod(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoNum (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL rhoDen (1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      CHARACTER*(MAX_LEN_MBUF) msgBuf


C     == local adjoint variables ==
      _RL deltasig_ad, facpres_ad, sp_ad, tp_ad, theta_v_ad
      _RL bulkmod_ad(1-olx:snx+olx,1-oly:sny+oly)
      _RL locpres_ad(1-olx:snx+olx,1-oly:sny+oly)
      _RL rhoden_ad(1-olx:snx+olx,1-oly:sny+oly)
      _RL rhonum_ad(1-olx:snx+olx,1-oly:sny+oly)
      _RL rhop0_ad(1-olx:snx+olx,1-oly:sny+oly)
CEOP

C----------------------------------------------
C reset some local adjoint variables...
C----------------------------------------------
      deltasig_ad = 0.d0
      facpres_ad = 0.d0
      sp_ad = 0.d0
      theta_v_ad = 0.d0
      tp_ad = 0.d0
      do ip2 = 1-oly, sny+oly
        do ip1 = 1-olx, snx+olx
          bulkmod_ad(ip1,ip2) = 0.d0
          locpres_ad(ip1,ip2) = 0.d0
          rhoden_ad(ip1,ip2) = 0.d0
          rhonum_ad(ip1,ip2) = 0.d0
          rhop0_ad(ip1,ip2) = 0.d0
        end do
      end do

      if (equationofstate .eq. 'LINEAR') then
        do j = jmin, jmax
          do i = imin, imax
            sfld_ad(i,j) = sfld_ad(i,j)+rholoc_ad(i,j)*rhonil*sbeta
            tfld_ad(i,j) = tfld_ad(i,j)-rholoc_ad(i,j)*rhonil*talpha
          end do
        end do
      else if (equationofstate .eq. 'POLY3') then
        reftemp = eosreft(kref)
        refsalt = eosrefs(kref)
        do j = jmin, jmax
          deltasig_ad = 0.d0
          sp_ad = 0.d0
          tp_ad = 0.d0
          do i = imin, imax
            deltasig_ad = 0.d0
            sp_ad = 0.d0
            tp_ad = 0.d0
            tp = tfld(i,j)-reftemp
            sp = sfld(i,j)-refsalt
            deltasig_ad = deltasig_ad+rholoc_ad(i,j)
            sp_ad = sp_ad+deltasig_ad*((eosc(9,kref)*sp+eosc(5,kref))*
     &sp+eosc(2,kref)+(eosc(9,kref)*sp+eosc(5,kref)+eosc(9,kref)*sp)*sp+
     &(eosc(7,kref)*tp+eosc(8,kref)*sp+eosc(4,kref)+eosc(8,kref)*sp)*tp)
            tp_ad = tp_ad+deltasig_ad*((eosc(6,kref)*tp+eosc(7,kref)*sp+
     &eosc(3,kref))*tp+(eosc(8,kref)*sp+eosc(4,kref))*sp+eosc(1,kref)+(
     &eosc(6,kref)*tp+eosc(7,kref)*sp+eosc(3,kref)+eosc(6,kref)*tp)*tp)
            deltasig_ad = 0.d0
            sfld_ad(i,j) = sfld_ad(i,j)+sp_ad
            sp_ad = 0.d0
            tfld_ad(i,j) = tfld_ad(i,j)+tp_ad
            tp_ad = 0.d0
          end do
        end do
      else if (equationofstate(1:5) .eq. 'JMD95' .or. equationofstate 
     &.eq. 'UNESCO') then
        call pressure_for_eos( bi,bj,imin,imax,jmin,jmax,kref,locpres,
     &mythid )
        call find_rhop0( imin,imax,jmin,jmax,tfld,sfld,rhop0,mythid )
        call find_bulkmod( imin,imax,jmin,jmax,locpres,tfld,sfld,
     &bulkmod,mythid )
        do j = jmin, jmax
          do i = imin, imax
            bulkmod_ad(i,j) = bulkmod_ad(i,j)-rholoc_ad(i,j)*(rhop0(i,j)
     &*(locpres(i,j)*sitobar/bulkmod(i,j)/bulkmod(i,j))/(1.d0-locpres(i,
     &j)*sitobar/bulkmod(i,j))/(1.d0-locpres(i,j)*sitobar/bulkmod(i,j)))
            locpres_ad(i,j) = locpres_ad(i,j)+rholoc_ad(i,j)*(rhop0(i,j)
     &*(sitobar/bulkmod(i,j))/(1.d0-locpres(i,j)*sitobar/bulkmod(i,j))/(
     &1.d0-locpres(i,j)*sitobar/bulkmod(i,j)))
            rhop0_ad(i,j) = rhop0_ad(i,j)+rholoc_ad(i,j)/(1.d0-locpres(
     &i,j)*sitobar/bulkmod(i,j))
          end do
        end do
        call find_bulkmod_ad( imin,imax,jmin,jmax,locpres,locpres_ad,
     &tfld,tfld_ad,sfld,sfld_ad,bulkmod_ad )
        call find_rhop0_ad( imin,imax,jmin,jmax,tfld,tfld_ad,sfld,
     &sfld_ad,rhop0_ad )
        call pressure_for_eos_ad( bi,bj,kref,locpres_ad )
      else if (equationofstate .eq. 'MDJWF') then
        call pressure_for_eos( bi,bj,imin,imax,jmin,jmax,kref,locpres,
     &mythid )
        call find_rhonum( imin,imax,jmin,jmax,locpres,tfld,sfld,rhonum,
     &mythid )
        call find_rhoden( imin,imax,jmin,jmax,locpres,tfld,sfld,rhoden,
     &mythid )
        do j = jmin, jmax
          do i = imin, imax
            rhoden_ad(i,j) = rhoden_ad(i,j)+rholoc_ad(i,j)*rhonum(i,j)
            rhonum_ad(i,j) = rhonum_ad(i,j)+rholoc_ad(i,j)*rhoden(i,j)
          end do
        end do
        call find_rhoden_ad( imin,imax,jmin,jmax,locpres,locpres_ad,
     &tfld,tfld_ad,sfld,sfld_ad,rhoden_ad )
        call find_rhonum_ad( imin,imax,jmin,jmax,locpres,locpres_ad,
     &tfld,tfld_ad,sfld,sfld_ad,rhonum_ad )
        call pressure_for_eos_ad( bi,bj,kref,locpres_ad )
      else if (equationofstate .eq. 'TEOS10') then
        call pressure_for_eos( bi,bj,imin,imax,jmin,jmax,kref,locpres,
     &mythid )
        call find_rhoteos( imin,imax,jmin,jmax,locpres,tfld,sfld,rhonum,
     &rhoden,mythid )
        do j = jmin, jmax
          do i = imin, imax
            rhoden_ad(i,j) = rhoden_ad(i,j)+rholoc_ad(i,j)*rhonum(i,j)
            rhonum_ad(i,j) = rhonum_ad(i,j)+rholoc_ad(i,j)*rhoden(i,j)
          end do
        end do
        call find_rhoteos_ad( imin,imax,jmin,jmax,locpres,locpres_ad,
     &tfld,tfld_ad,sfld,sfld_ad,rhonum_ad,rhoden_ad )
        call pressure_for_eos_ad( bi,bj,kref,locpres_ad )
      else if (equationofstate .eq. 'IDEALG') then
        call pressure_for_eos( bi,bj,imin,imax,jmin,jmax,kref,locpres,
     &mythid )
        oneminuskap = onerl-atm_kappa
        do j = jmin, jmax
          facpres_ad = 0.d0
          theta_v_ad = 0.d0
          do i = imin, imax
            facpres_ad = 0.d0
            theta_v_ad = 0.d0
            if (locpres(i,j) .gt. zerorl .and. tfld(i,j) .gt. zerorl) 
     &then
              facpres = (locpres(i,j)/atm_po)**oneminuskap
              theta_v = tfld(i,j)*(sfld(i,j)*atm_rq+onerl)
              facpres_ad = facpres_ad+rholoc_ad(i,j)*(atm_po/(atm_rd*
     &theta_v))
              theta_v_ad = theta_v_ad-rholoc_ad(i,j)*(atm_po*facpres*
     &atm_rd/(atm_rd*theta_v)/(atm_rd*theta_v))
              sfld_ad(i,j) = sfld_ad(i,j)+theta_v_ad*tfld(i,j)*atm_rq
              tfld_ad(i,j) = tfld_ad(i,j)+theta_v_ad*(sfld(i,j)*atm_rq+
     &onerl)
              theta_v_ad = 0.d0
              locpres_ad(i,j) = locpres_ad(i,j)+facpres_ad/atm_po*
     &oneminuskap*(locpres(i,j)/atm_po)**(oneminuskap-1)
              facpres_ad = 0.d0
            else
            endif
          end do
        end do
        call pressure_for_eos_ad( bi,bj,kref,locpres_ad )
      endif

      end
