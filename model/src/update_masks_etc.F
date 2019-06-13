#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif

C--  File update_masks_etc.F:
C--   Contents
C--   o S/R UPDATE_MASKS_ETC
C--   o FCT SMOOTHMIN_RS( a, b )
C--   o FCT SMOOTHMIN_RL( a, b )
C--   o FCT SMOOTHABS_RS( x )
C--   o FCT SMOOTHABS_RL( x )
Cml   o S/R LIMIT_HFACC_TO_ONE
Cml   o S/R ADLIMIT_HFACC_TO_ONE

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: UPDATE_MASKS_ETC
C     !INTERFACE:
      SUBROUTINE UPDATE_MASKS_ETC( myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE UPDATE_MASKS_ETC
C     | o Re-initialise masks and topography factors after a new
C     |   hFacC has been calculated by the minimizer
C     *==========================================================*
C     | These arrays are used throughout the code and describe
C     | the topography of the domain through masks (0s and 1s)
C     | and fractional height factors (0<hFac<1). The latter
C     | distinguish between the lopped-cell and full-step
C     | topographic representations.
C     |
C     | Updated by Tim Smith 
C     |   Note: This was rewritten closely following ini_masks_etc 
C     |   to work when shelfice package is used. 
C     |   The following is assumed:
C     |   1. selectSigmaCoord.EQ.0: otherwise, you're on your own :) 
C     |   2. useMin4hFacEdges = .false. : this only produces different results
C     |      when useShelfice=.true., which is more general than when
C     |      useMin4hFacEdges=.true.. When useShelfice=.false., this
C     |      method will produce equivalent results to useMin4hFacEdges=.true.
C     |   3. call to ADD_WALLS2MASKS is not included
C     |
C     *==========================================================*
C     | code taken from ini_masks_etc.F
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SURFACE.h"
Cml we need optimcycle for storing the new hFaC(C/W/S) and depth
Cts need tamc header for store directives 
#ifdef ALLOW_AUTODIFF
# include "tamc.h"
# include "tamc_keys.h"
# include "optim.h"
#endif

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myThid -  Number of this instance of INI_MASKS_ETC
      INTEGER myThid

#ifdef ALLOW_DEPTH_CONTROL
C     !FUNCTIONS:
      _RS SMOOTHMIN_RS
      EXTERNAL SMOOTHMIN_RS

C     !LOCAL VARIABLES:
C     == Local variables ==
C     bi,bj   :: Loop counters
C     I,J,K
C     tmpfld  :: Temporary array used to compute & write Total Depth
      INTEGER bi, bj
      INTEGER I, J, K
      _RS tmpfld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tmpVar(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL tmpVar1(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL tmpVar2(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      CHARACTER*(MAX_LEN_MBUF) suff

      INTEGER Im1, Jm1
      _RL hFacCtmp
      _RL hFac1tmp, hFac2tmp
      _RL hFacMnSz

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|


      

Cts ---
C   Calculate hFacC from rLow
C   Here we see the first impact of r_low as a control vector.
Cts ---
CADJ STORE r_low = tapelev_init, key = 1
CADJ STORE ro_surf = tapelev_init, key = 1
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO K=1, Nr
         hFacMnSz=max( hFacMin, min(hFacMinDr*recip_drF(k),1. _d 0) )
         DO J=1-OLy,sNy+OLy
          DO I=1-OLx,sNx+OLx
C      o Non-dimensional distance between grid bound. and domain lower_R bound.
#if defined(ALLOW_DEPTH_CONTROL) && !(defined ALLOW_GENARR2D_CONTROL)

           ! With old control setup, xx_r_low = r_low + xx_depth 
           hFacCtmp = (rF(K)-xx_r_low(I,J,bi,bj))*recip_drF(K)
#else
           ! In genarr2d setup, r_low = r_low_const + xx_depth
           ! i.e. TAF is differentiating w.r.t. anomaly not full field
           hFacCtmp = (rF(K)-R_low(I,J,bi,bj))*recip_drF(K)
#endif /* ALLOW_DEPTH_CONTROL && !ALLOW_GENARR2D_CONTROL */

C      o Select between closed, open, or partial (0, 1, 0-1)
           hFacCtmp = min(max(hFacCtmp,zeroRL),oneRL)

C      o Impose min. fraction and/or size (dimensional
           if ( hFacCtmp.lt.hFacMnSz ) then
            if ( hFacCtmp.lt.hFacMnSz*halfRL ) then
             hFacC(i,j,k,bi,bj) = 0. 
            else
             hFacC(i,j,k,bi,bj) = hFacMnSz
            endif
           else
            hFacC(i,j,k,bi,bj) = hFacCtmp
           endif
          enddo
         enddo
        enddo

Cts ---
C   Still inside bi,bj loop
C   Re-calculate lower-R Boundary position, taking into account hFacC
Cts ---
        DO J=1-OLy,sNy+OLy
         DO I=1-OLx,sNx+OLx
          R_low(i,j,bi,bj) = rF(1)
         ENDDO
        ENDDO
        DO K=Nr,1,-1
         DO J=1-OLy,sNy+OLy
          DO I=1-OLx,sNx+OLx
           R_low(I,J,bi,bj) = R_low(I,J,bi,bj)
     &                      - drF(K)*hFacC(I,J,K,bi,bj)
          ENDDO
         ENDDO
        ENDDO

Cts --- 
C   Still inside bi,bj loop
C   Compute hFacC from Ro_surf
C   This is important when useShelfice=.true. since the call to
C   shelfice_init_depths resets Ro_surf based on ice shelf draft. 
C   This may not be necessary when useShelfice=.false., but will still
C   produce correct results.
Cts ---
        DO k=1, Nr
         hFacMnSz = MAX( hFacMin, MIN(hFacMinDr*recip_drF(k),oneRL) )
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
C      o Non-dimensional distance between grid boundary and model surface
           hFacCtmp = (rF(k)-Ro_surf(i,j,bi,bj))*recip_drF(k)
C      o Reduce the previous fraction : substract the outside part.
           hFacCtmp = hFacC(i,j,k,bi,bj) - MAX( hFacCtmp, zeroRL )
C      o set to zero if empty Column :
           hFacCtmp = MAX( hFacCtmp, zeroRL )
C      o Impose minimum fraction and/or size (dimensional)
           IF ( hFacCtmp.LT.hFacMnSz ) THEN
            IF ( hFacCtmp.LT.hFacMnSz*halfRL ) THEN
             hFacC(i,j,k,bi,bj) = 0.
            ELSE
             hFacC(i,j,k,bi,bj) = hFacMnSz
            ENDIF
           ELSE
             hFacC(i,j,k,bi,bj) = hFacCtmp
           ENDIF
          ENDDO
         ENDDO
        ENDDO

Cts ---
C   Still inside bi,bj loop
C   Recompute Ro_surf, kSurfC, kLowC from hFacC updates
Cts ---
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          tmpVar(i,j) = 0.
          tmpFld(i,j,bi,bj) = 0.
          kSurfC(i,j,bi,bj) = Nr+1
          kLowC (i,j,bi,bj) = 0
         ENDDO
        ENDDO
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           tmpVar(i,j) = tmpVar(i,j) + drF(k)*hFacC(i,j,k,bi,bj)
           tmpFld(i,j,bi,bj) = tmpFld(i,j,bi,bj) + 1.
           IF ( hFacC(i,j,k,bi,bj).NE.zeroRS ) kLowC(i,j,bi,bj) = k
          ENDDO
         ENDDO
        ENDDO
        DO k=Nr,1,-1
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           IF ( hFacC(i,j,k,bi,bj).NE.zeroRS ) kSurfC(i,j,bi,bj) = k
          ENDDO
         ENDDO
        ENDDO
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          Ro_surf(i,j,bi,bj) = R_low(i,j,bi,bj) + tmpVar(i,j)
          maskInC(i,j,bi,bj) = 0.
          IF ( kSurfC(i,j,bi,bj).LE.Nr ) maskInC(i,j,bi,bj) = 1.
         ENDDO
        ENDDO

C-    end bi,bj loops.
       ENDDO
      ENDDO


      IF ( plotLevel.GE.debLevC ) THEN
        _BARRIER
        CALL PLOT_FIELD_XYRS( R_low,
     &         'Model R_low (update_masks_etc)', 1, myThid )
Cts ---
C   Need to print out b/c ice shelf modifies Ro_surf
C   (To be consistent with ini_masks_etc)
Cts ---
        CALL PLOT_FIELD_XYRS( Ro_surf,
     &         'Model Ro_surf (update_masks_etc)', 1, myThid )
      ENDIF

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C     Calculate quantities derived from XY depth map
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
C         Total fluid column thickness (r_unit) :
          tmpfld(i,j,bi,bj) = Ro_surf(i,j,bi,bj) - R_low(i,j,bi,bj)
C         Inverse of fluid column thickness (1/r_unit)
          IF ( tmpfld(i,j,bi,bj) .LE. 0. ) THEN
           recip_Rcol(i,j,bi,bj) = 0.
          ELSE
           recip_Rcol(i,j,bi,bj) = 1. _d 0 / tmpfld(i,j,bi,bj)
          ENDIF
         ENDDO
        ENDDO

Cts ---
C   Still in bi,bj loop
C   Recompute hFacW/S, rSurfW/S, rLow/S with same method as with hFacC
C   above. This is from method 2 in ini_masks_etc      
Cts ---

#ifdef ALLOW_AUTODIFF_TAMC
        act1 = bi - myBxLo(myThid)
        max1 = myBxHi(myThid) - myBxLo(myThid) + 1
        act2 = bj - myByLo(myThid)
        max2 = myByHi(myThid) - myByLo(myThid) + 1
        act3 = myThid - 1
        max3 = nTx*nTy
        act4 = ikey_dynamics - 1 
        iumekey = (act1 + 1) + act2*max1
     &            + act3*max1*max2
     &            + act4*max1*max2*max3

CADJ STORE rLowS(:,:,bi,bj) = 
CADJ &     comlev1_bibj, key=iumekey, byte=isbyte
CADJ STORE rLowW(:,:,bi,bj) = 
CADJ &     comlev1_bibj, key=iumekey, byte=isbyte
CADJ STORE rSurfS(:,:,bi,bj) = 
CADJ &     comlev1_bibj, key=iumekey, byte=isbyte
CADJ STORE rSurfW(:,:,bi,bj) = 
CADJ &     comlev1_bibj, key=iumekey, byte=isbyte
#endif 

         DO k=1, Nr
          hFacMnSz = MAX( hFacMin, MIN(hFacMinDr*recip_drF(k),oneRL) )
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
C      o Non-dimensional distance between grid bound. and domain lower_R bound.
            hFac1tmp = ( rF(k) - rLowW(i,j,bi,bj) )*recip_drF(k)
            hFacCtmp = MIN( hFac1tmp, oneRL )
C      o Impose minimum fraction and/or size (dimensional)
            IF ( hFacCtmp.LT.hFacMnSz*halfRL ) THEN
              hFac1tmp = 0.
            ELSE
              hFac1tmp = MAX( hFacCtmp, hFacMnSz )
            ENDIF
C      o Reduce the previous fraction : substract the outside fraction
C        (i.e., beyond reference (=at rest) surface position rSurfW)
            hFac2tmp = ( rF(k) -rSurfW(i,j,bi,bj) )*recip_drF(k)
            hFacCtmp = hFac1tmp - MAX( hFac2tmp, zeroRL )
C      o Impose minimum fraction and/or size (dimensional)
            IF ( hFacCtmp.LT.hFacMnSz*halfRL ) THEN
              hFacW(i,j,k,bi,bj) = 0.
            ELSE
              hFacW(i,j,k,bi,bj) = MAX( hFacCtmp, hFacMnSz )
            ENDIF
           ENDDO
          ENDDO
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
C      o Non-dimensional distance between grid bound. and domain lower_R bound.
            hFac1tmp = ( rF(k) - rLowS(i,j,bi,bj) )*recip_drF(k)
            hFacCtmp = MIN( hFac1tmp, oneRL )
C      o Impose minimum fraction and/or size (dimensional)
            IF ( hFacCtmp.LT.hFacMnSz*halfRL ) THEN
              hFac1tmp = 0.
            ELSE
              hFac1tmp = MAX( hFacCtmp, hFacMnSz )
            ENDIF
C      o Reduce the previous fraction : substract the outside fraction
C        (i.e., beyond reference (=at rest) surface position rSurfS)
            hFac2tmp = ( rF(k) -rSurfS(i,j,bi,bj) )*recip_drF(k)
            hFacCtmp = hFac1tmp - MAX( hFac2tmp, zeroRL )
C      o Impose minimum fraction and/or size (dimensional)
            IF ( hFacCtmp.LT.hFacMnSz*halfRL ) THEN
              hFacS(i,j,k,bi,bj) = 0.
            ELSE
              hFacS(i,j,k,bi,bj) = MAX( hFacCtmp, hFacMnSz )
            ENDIF
           ENDDO
          ENDDO
         ENDDO

Cts ---
C   Still in bi,bj loop
C   Update rLow, rSurf @ West/South edges accounting for adjusted R_low
C   & Ro_surf due to hFacMin constraint on hFacC.
Cts ---
        DO j=1-OLy,sNy+OLy
         DO i=2-OLx,sNx+OLx
           rLowW(i,j,bi,bj)  =
     &           MAX(   R_low(i-1,j,bi,bj),   R_low(i,j,bi,bj) )
           rSurfW(i,j,bi,bj) =
     &           MIN( Ro_surf(i-1,j,bi,bj), Ro_surf(i,j,bi,bj) )
           rSurfW(i,j,bi,bj) =
     &           MAX( rSurfW(i,j,bi,bj), rLowW(i,j,bi,bj) )
         ENDDO
        ENDDO
        DO j=2-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
           rLowS(i,j,bi,bj)  =
     &           MAX(   R_low(i,j-1,bi,bj),   R_low(i,j,bi,bj) )
           rSurfS(i,j,bi,bj) =
     &           MIN( Ro_surf(i,j-1,bi,bj), Ro_surf(i,j,bi,bj) )
           rSurfS(i,j,bi,bj) =
     &           MAX( rSurfS(i,j,bi,bj), rLowS(i,j,bi,bj) )
         ENDDO
        ENDDO

Cts ---
C   Still in bi,bj loop
C   The following seems only important for useShelfice=.true., but is
C   kept in general in ini_masks_etc, so I'll keep it here too.
Cts ---

C   Adjust rLow & reference rSurf at Western & Southern edges (U & V pts)
C   to get consistent column thickness from Sum_k(hFac*drF) and rSurf-rLow

C   1. Compute total column thickness at Western & Southern edges
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
            tmpVar1(i,j) = 0. _d 0
            tmpVar2(i,j) = 0. _d 0
          ENDDO
         ENDDO
         DO k=1,Nr
          DO j=1-OLy,sNy+OLy
           DO i=1-OLx,sNx+OLx
            tmpVar1(i,j) = tmpVar1(i,j) + drF(k)*hFacW(i,j,k,bi,bj)
            tmpVar2(i,j) = tmpVar2(i,j) + drF(k)*hFacS(i,j,k,bi,bj)
           ENDDO
          ENDDO
         ENDDO

C   2. Adjust both rLow and rSurf at W & S edges (split correction by half)
C   adjust rSurfW and rLowW:
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           tmpVar1(i,j) = rLowW(i,j,bi,bj) + tmpVar1(i,j)
           tmpVar1(i,j) = ( tmpVar1(i,j) -rSurfW(i,j,bi,bj) )*halfRL
          ENDDO
         ENDDO
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           rSurfW(i,j,bi,bj) = rSurfW(i,j,bi,bj) + tmpVar1(i,j)
           rLowW (i,j,bi,bj) = rLowW (i,j,bi,bj) - tmpVar1(i,j)
          ENDDO
         ENDDO

C   3. Adjust rSurfS and rLowS:
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           tmpVar2(i,j) = rLowS(i,j,bi,bj) + tmpVar2(i,j)
           tmpVar2(i,j) = ( tmpVar2(i,j) -rSurfS(i,j,bi,bj) )*halfRL
          ENDDO
         ENDDO
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           rSurfS(i,j,bi,bj) = rSurfS(i,j,bi,bj) + tmpVar2(i,j)
           rLowS (i,j,bi,bj) = rLowS (i,j,bi,bj) - tmpVar2(i,j)
          ENDDO
         ENDDO

C-    end bi,bj loops.
       ENDDO
      ENDDO

Cts ---
C   Note this is done after kSurfW/S & maskInW/S computation in
C   ini_masks_etc. Not sure if important, following ML's version here.
Cts ---
      CALL EXCH_UV_XYZ_RS(hFacW,hFacS,.FALSE.,myThid)
      CALL EXCH_UV_XY_RS(rSurfW,rSurfS,.FALSE.,myThid)
      CALL EXCH_UV_XY_RS(rLowW,rLowS,.FALSE.,myThid)

      IF ( plotLevel.GE.debLevC ) THEN
        _BARRIER
C--   Write to monitor file (standard output)
        CALL PLOT_FIELD_XYZRS( hFacC,'hFacC (update_masks_etc)',
     &                                          Nr, 1, myThid )
        CALL PLOT_FIELD_XYZRS( hFacW,'hFacW (update_masks_etc)',
     &                                          Nr, 1, myThid )
        CALL PLOT_FIELD_XYZRS( hFacS,'hFacS (update_masks_etc)',
     &                                          Nr, 1, myThid )
      ENDIF

C--   Calculate surface k index for interface W & S (U & V points)
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        DO j=1-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx
          kSurfW(i,j,bi,bj) = Nr+1
          kSurfS(i,j,bi,bj) = Nr+1
          DO k=Nr,1,-1
           IF (hFacW(i,j,k,bi,bj).NE.zeroRS) kSurfW(i,j,bi,bj) = k
           IF (hFacS(i,j,k,bi,bj).NE.zeroRS) kSurfS(i,j,bi,bj) = k
          ENDDO
          maskInW(i,j,bi,bj)= 0.
          IF ( kSurfW(i,j,bi,bj).LE.Nr ) maskInW(i,j,bi,bj)= 1.
          maskInS(i,j,bi,bj)= 0.
          IF ( kSurfS(i,j,bi,bj).LE.Nr ) maskInS(i,j,bi,bj)= 1.
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C     Masks and reciprocals of hFac[CWS]
      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO K=1,Nr
         DO J=1-OLy,sNy+OLy
          DO I=1-OLx,sNx+OLx
           IF (hFacC(I,J,K,bi,bj) .NE. 0. ) THEN
            recip_hFacC(I,J,K,bi,bj) = 1. _d 0 / hFacC(I,J,K,bi,bj)
            maskC(I,J,K,bi,bj) = 1.
           ELSE
            recip_hFacC(I,J,K,bi,bj) = 0.
            maskC(I,J,K,bi,bj) = 0.
           ENDIF
           IF (hFacW(I,J,K,bi,bj) .NE. 0. ) THEN
            recip_hFacW(I,J,K,bi,bj) = 1. _d 0 / hFacw(I,J,K,bi,bj)
            maskW(I,J,K,bi,bj) = 1.
           ELSE
            recip_hFacW(I,J,K,bi,bj) = 0.
            maskW(I,J,K,bi,bj) = 0.
           ENDIF
           IF (hFacS(I,J,K,bi,bj) .NE. 0. ) THEN
            recip_hFacS(I,J,K,bi,bj) = 1. _d 0 / hFacS(I,J,K,bi,bj)
            maskS(I,J,K,bi,bj) = 1.
           ELSE
            recip_hFacS(I,J,K,bi,bj) = 0.
            maskS(I,J,K,bi,bj) = 0.
           ENDIF
          ENDDO
         ENDDO
        ENDDO

#ifdef NONLIN_FRSURF
C--   Save initial geometrical hFac factor into h0Fac (fixed in time):
C     Note: In case 1 pkg modifies hFac (from packages_init_fixed, called
C     later in sequence of calls) this pkg would need also to update h0Fac.
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           h0FacC(i,j,k,bi,bj) = _hFacC(i,j,k,bi,bj)
           h0FacW(i,j,k,bi,bj) = _hFacW(i,j,k,bi,bj)
           h0FacS(i,j,k,bi,bj) = _hFacS(i,j,k,bi,bj)
          ENDDO
         ENDDO
        ENDDO
#endif /* NONLIN_FRSURF */
C - end bi,bj loops.
       ENDDO
      ENDDO

#endif /* ALLOW_DEPTH_CONTROL */
      RETURN
      END

#ifdef USE_SMOOTH_MIN
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      _RS FUNCTION SMOOTHMIN_RS( a, b )

      IMPLICIT NONE

      _RS a, b

      _RS SMOOTHABS_RS
      EXTERNAL SMOOTHABS_RS

Cml      smoothMin_R4 = .5*(a+b)
      SMOOTHMIN_RS = .5*( a+b - SMOOTHABS_RS(a-b) )
CML      smoothMin_R4 = MIN(a,b)

      RETURN
      END

      _RL FUNCTION SMOOTHMIN_RL( a, b )

      IMPLICIT NONE

      _RL a, b

      _RL SMOOTHABS_RL
      EXTERNAL SMOOTHABS_RL

Cml      smoothMin_R8 = .5*(a+b)
      SMOOTHMIN_RL = .5*( a+b - SMOOTHABS_RL(a-b) )
Cml      smoothMin_R8 = MIN(a,b)

      RETURN
      END

      _RS FUNCTION SMOOTHABS_RS( x )

      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
C     input parameter
      _RS x
c     local variable
      _RS sf, rsf

      IF ( smoothAbsFuncRange .LT. 0.0 ) THEN
c     limit of smoothMin(a,b) = .5*(a+b)
         SMOOTHABS_RS = 0.
      ELSE
         IF ( smoothAbsFuncRange .NE. 0.0 ) THEN
            sf  = 10.0/smoothAbsFuncRange
            rsf = 1./sf
         ELSE
c     limit of smoothMin(a,b) = min(a,b)
            sf  = 0.
            rsf = 0.
         ENDIF
c
         IF ( x .GT. smoothAbsFuncRange ) THEN
            SMOOTHABS_RS = x
         ELSEIF ( x .LT. -smoothAbsFuncRange ) THEN
            SMOOTHABS_RS = -x
         ELSE
            SMOOTHABS_RS = log(.5*(exp(x*sf)+exp(-x*sf)))*rsf
         ENDIF
      ENDIF

      RETURN
      END

      _RL FUNCTION SMOOTHABS_RL( x )

      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
C     input parameter
      _RL x
c     local variable
      _RL sf, rsf

      IF ( smoothAbsFuncRange .LT. 0.0 ) THEN
c     limit of smoothMin(a,b) = .5*(a+b)
         SMOOTHABS_RL = 0.
      ELSE
         IF ( smoothAbsFuncRange .NE. 0.0 ) THEN
            sf  = 10.0D0/smoothAbsFuncRange
            rsf = 1.D0/sf
         ELSE
c     limit of smoothMin(a,b) = min(a,b)
            sf  = 0.D0
            rsf = 0.D0
         ENDIF
c
         IF ( x .GE. smoothAbsFuncRange ) THEN
            SMOOTHABS_RL = x
         ELSEIF ( x .LE. -smoothAbsFuncRange ) THEN
            SMOOTHABS_RL = -x
         ELSE
            SMOOTHABS_RL = log(.5*(exp(x*sf)+exp(-x*sf)))*rsf
         ENDIF
      ENDIF

      RETURN
      END
#endif /* USE_SMOOTH_MIN */
