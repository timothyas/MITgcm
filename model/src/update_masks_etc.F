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
C--   o FCT REG_hFac_RL( hFac, hFacMin )
C       -> the regular, if statement version of computing hFac
C--   o FCT smooth_hFac_RL( hFac, hFacMin )
C       -> the smooth version of computing hFac
C--   o FCT SMOOTHMIN_RS( a, b )
C--   o FCT SMOOTHMIN_RL( a, b )
C--   o FCT SMOOTHMAX_RS( a, b )
C--   o FCT SMOOTHMAX_RL( a, b )
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
      _RL smooth_hFac_RL
      EXTERNAL smooth_hFac_RL
      _RL reg_hFac_RL
      EXTERNAL reg_hFac_RL
      _RL SMOOTHMIN_RL
      EXTERNAL SMOOTHMIN_RL
      _RL SMOOTHMAX_RL
      EXTERNAL SMOOTHMAX_RL

C     !LOCAL VARIABLES:
C     == Local variables ==
C     bi,bj   :: Loop counters
C     I,J,K
C     tmpfld  :: Temporary array used to compute & write Total Depth
      INTEGER bi, bj
      INTEGER I, J, K
      _RL tmpfld(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tmpVar(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL tmpVarW(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL tmpVarS(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL tmpVar1(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL tmpVar2(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      CHARACTER*(MAX_LEN_MBUF) suff

      INTEGER Im1, Jm1
      _RL hFacCtmp, hFacWtmp, hFacStmp
      _RL hFac1tmp, hFac2tmp
      _RL hFacMnSz

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

Cts ---
C   For depth optimization:
C   Can implement bounds on ctrls so that they don't go above or below
C   max depth, but in case of ice shelf, need to ensure additional topo
C   doesn't go past ice shelf geometry
#ifdef ALLOW_SHELFICE
      do bj=myByLo(myThid), myByHi(myThid) 
       do bi=myBxLo(myThid), myBxHi(myThid) 
        do j=1-OLy,sNy+OLy
         do i=1-OLx,sNx+OLx
          r_low(i,j,bi,bj) = min(r_low(i,j,bi,bj),ro_surf(i,j,bi,bj))
         enddo
        enddo
       enddo
      enddo
#endif


C   Calculate hFac from rLow
CADJ STORE r_low = tapelev_init, key = 1
CADJ STORE ro_surf = tapelev_init, key = 1

C   Re-Initialise rLow & reference rSurf at Western & Southern edges (U & V pts)
C   Note: not final value since these estimates ignore hFacMin constrain
      DO bj=myByLo(myThid), myByHi(myThid)
       DO bi=myBxLo(myThid), myBxHi(myThid)
        i = 1-OLx
        DO j=1-OLy,sNy+OLy
           rLowW (i,j,bi,bj) = rF(1)
           rSurfW(i,j,bi,bj) = rF(1)
        ENDDO
        j = 1-OLy
        DO i=1-OLx,sNx+OLx
           rLowS (i,j,bi,bj) = rF(1)
           rSurfS(i,j,bi,bj) = rF(1)
        ENDDO

        DO j=1-OLy,sNy+OLy
         DO i=2-OLx,sNx+OLx
           rLowW(i,j,bi,bj)  =
#ifdef USE_SMOOTH_MIN
     &  smoothmax_rl(   R_low(i-1,j,bi,bj),   R_low(i,j,bi,bj) )
#else
     &           MAX(   R_low(i-1,j,bi,bj),   R_low(i,j,bi,bj) )
#endif

           rSurfW(i,j,bi,bj) =
#ifdef USE_SMOOTH_MIN
     &  smoothmin_rl( Ro_surf(i-1,j,bi,bj), Ro_surf(i,j,bi,bj) )
#else
     &           MIN( Ro_surf(i-1,j,bi,bj), Ro_surf(i,j,bi,bj) )
#endif
           rSurfW(i,j,bi,bj) =
#ifdef USE_SMOOTH_MIN
     &  smoothmax_rl( rSurfW(i,j,bi,bj), rLowW(i,j,bi,bj) )
#else
     &           MAX( rSurfW(i,j,bi,bj), rLowW(i,j,bi,bj) )
#endif
         ENDDO
        ENDDO
        DO j=2-OLy,sNy+OLy
         DO i=1-OLx,sNx+OLx

           rLowS(i,j,bi,bj)  =
#ifdef USE_SMOOTH_MIN
     &  smoothmax_rl(   R_low(i,j-1,bi,bj),   R_low(i,j,bi,bj) )
#else
     &           MAX(   R_low(i,j-1,bi,bj),   R_low(i,j,bi,bj) )
#endif
           rSurfS(i,j,bi,bj) =
#ifdef USE_SMOOTH_MIN
     &  smoothmin_rl( Ro_surf(i,j-1,bi,bj), Ro_surf(i,j,bi,bj) )
#else
     &           MIN( Ro_surf(i,j-1,bi,bj), Ro_surf(i,j,bi,bj) )
#endif
           rSurfS(i,j,bi,bj) =
#ifdef USE_SMOOTH_MIN
     &  smoothmax_rl( rSurfS(i,j,bi,bj), rLowS(i,j,bi,bj) )
#else
     &           MAX( rSurfS(i,j,bi,bj), rLowS(i,j,bi,bj) )
#endif
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C   Calculate hFac from rLow
CADJ STORE rloww = tapelev_init, key = 1
CADJ STORE rlows = tapelev_init, key = 1
CADJ STORE rsurfw = tapelev_init, key = 1
CADJ STORE rsurfs = tapelev_init, key = 1
CADJ STORE hfacc = tapelev_init, key = 1
CADJ STORE hfacw = tapelev_init, key = 1
CADJ STORE hfacs = tapelev_init, key = 1
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
           hFacWtmp = (rF(K)-rLowW(I,J,bi,bj))*recip_drF(K)
           hFacStmp = (rF(K)-rLowS(I,J,bi,bj))*recip_drF(K)
#endif /* ALLOW_DEPTH_CONTROL && !ALLOW_GENARR2D_CONTROL */

C     o Impose min. fraction and/or size (dimensional)
#ifdef USE_SMOOTH_MIN
           hFacC(i,j,k,bi,bj) = 
     &          smooth_hFac_rl( hFacCtmp, hFacMnSz, smoothHFacEps )
           hFacW(i,j,k,bi,bj) = 
     &          smooth_hFac_rl( hFacWtmp, hFacMnSz, smoothHFacEps )
           hFacS(i,j,k,bi,bj) = 
     &          smooth_hFac_rl( hFacStmp, hFacMnSz, smoothHFacEps )
#else /* ndef smooth_min */
           hFacCtmp = min( max( hFacCtmp, zeroRL ) , oneRL )
           hFacWtmp = min( hFacWtmp, oneRL )
           hFacStmp = min( hFacStmp, oneRL )
           hFacC(i,j,k,bi,bj) = reg_hFac_rl( hFacCtmp, hFacMnSz )
           hFacW(i,j,k,bi,bj) = reg_hFac_rl( hFacWtmp, hFacMnSz )
           hFacS(i,j,k,bi,bj) = reg_hFac_rl( hFacStmp, hFacMnSz )
#endif

          enddo
         enddo
        enddo

Cts ---
C   Still inside bi,bj loop
C   Re-calculate lower-R Boundary position, taking into account hFac
Cts ---
        DO J=1-OLy,sNy+OLy
         DO I=1-OLx,sNx+OLx
          R_low(i,j,bi,bj) = rF(1)
          rLowW(i,j,bi,bj) = rF(1)
          rLowS(i,j,bi,bj) = rF(1)
         ENDDO
        ENDDO
        DO K=Nr,1,-1
         DO J=1-OLy,sNy+OLy
          DO I=1-OLx,sNx+OLx
           R_low(I,J,bi,bj) = R_low(I,J,bi,bj)
     &                      - drF(K)*hFacC(I,J,K,bi,bj)
           rLowW(I,J,bi,bj) = rLowW(I,J,bi,bj)
     &                      - drF(K)*hFacW(I,J,K,bi,bj)
           rLowS(I,J,bi,bj) = rLowS(I,J,bi,bj)
     &                      - drF(K)*hFacS(I,J,K,bi,bj)
          ENDDO
         ENDDO
        ENDDO

Cts --- 
C   Still inside bi,bj loop
C   Compute hFac from rSurf
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
           hFacWtmp = (rF(k)-rSurfW(i,j,bi,bj))*recip_drF(k)
           hFacStmp = (rF(k)-rSurfS(i,j,bi,bj))*recip_drF(k)

C      o Reduce the previous fraction : subtract the outside part.
C      o set to zero if empty Column :
#ifdef USE_SMOOTH_MIN
           hFacCtmp = hFacC(i,j,k,bi,bj) - 
     &                smoothmax_rl( hFacCtmp, zeroRL )
           hFacCtmp = smoothmax_rl( hFacCtmp, zeroRL )

           hFacWtmp = hFacW(i,j,k,bi,bj) - 
     &                smoothmax_rl( hFacWtmp, zeroRL )

           hFacStmp = hFacS(i,j,k,bi,bj) - 
     &                smoothmax_rl( hFacStmp, zeroRL )
#else /* ndef smooth_min */
           hFacCtmp = hFacC(i,j,k,bi,bj) - MAX( hFacCtmp, zeroRL )
           hFacCtmp = MAX( hFacCtmp, zeroRL )

           hFacWtmp = hFacW(i,j,k,bi,bj) - MAX( hFacWtmp, zeroRL )

           hFacStmp = hFacS(i,j,k,bi,bj) - MAX( hFacStmp, zeroRL )
#endif


#ifdef USE_SMOOTH_MIN
           hFacC(i,j,k,bi,bj) = 
     &          smooth_hFac_rl( hFacCtmp, hFacMnSz, smoothHFacEps )

           hFacW(i,j,k,bi,bj) = 
     &          smooth_hFac_rl( hFacWtmp, hFacMnSz, smoothHFacEps )

           hFacS(i,j,k,bi,bj) = 
     &          smooth_hFac_rl( hFacStmp, hFacMnSz, smoothHFacEps )
#else /* ndef smooth_min */
           hFacC(i,j,k,bi,bj) = reg_hFac_rl( hFacCtmp, hFacMnSz )
           hFacW(i,j,k,bi,bj) = reg_hFac_rl( hFacWtmp, hFacMnSz )
           hFacS(i,j,k,bi,bj) = reg_hFac_rl( hFacStmp, hFacMnSz )
#endif
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
          tmpVarW(i,j) = 0.
          tmpVarS(i,j) = 0.
          tmpFld(i,j,bi,bj) = 0.
          kSurfC(i,j,bi,bj) = Nr+1
          kLowC (i,j,bi,bj) = 0
         ENDDO
        ENDDO
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           tmpVar(i,j) = tmpVar(i,j) + drF(k)*hFacC(i,j,k,bi,bj)
           tmpVarW(i,j) = tmpVarW(i,j) + drF(k)*hFacW(i,j,k,bi,bj)
           tmpVarS(i,j) = tmpVarS(i,j) + drF(k)*hFacS(i,j,k,bi,bj)
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
          rSurfW(i,j,bi,bj) = rLowW(i,j,bi,bj) + tmpVarW(i,j)
          rSurfS(i,j,bi,bj) = rLowS(i,j,bi,bj) + tmpVarS(i,j)
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

#if (defined ALLOW_AUTODIFF) && (defined ALLOW_AUTODIFF_MONITOR)
C     Dummy routine for printing ADJhFac
      CALL DUMMY_FOR_HFAC( 'C',0, 0, myThid )
      CALL DUMMY_FOR_HFAC( 'W',0, 0, myThid )
      CALL DUMMY_FOR_HFAC( 'S',0, 0, myThid )
#endif

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

#if (defined ALLOW_AUTODIFF) && (defined ALLOW_AUTODIFF_MONITOR)
C     Dummy routine for printing ADJ recip_hFac
      CALL DUMMY_FOR_HFAC( 'R',0, 0, myThid )
#endif

#endif /* ALLOW_DEPTH_CONTROL */
      RETURN
      END

#ifdef USE_SMOOTH_MIN
C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

      _RL FUNCTION smooth_hFac_rl( hFac, hFacMin, hFacEps )
C     o Impose 1,0, or hFacMin restriction smoothly
C     o Here use spline representation ...

      IMPLICIT NONE
#include "EEPARAMS.h"

      _RL hFac, hFacMin, hFacEps
      _RL low_threshold

      _RL hFac_spline_rl
      EXTERNAL hFac_spline_rl

      low_threshold = (5. _d -1) * hFacMin

      IF ( hFac .LT. low_threshold ) THEN
          smooth_hFac_rl = zeroRL
      ELSEIF ( hFac .LE. hFacMin ) THEN
          smooth_hFac_rl = hFacMin
      ELSEIF ( hFac .LE. hFacMin + hFacEps ) THEN
          smooth_hFac_rl = 
     &          hFac_spline_rl( hFac, hFacMin, hFacEps )
      ELSEIF ( hFac .LE. oneRL ) THEN
          smooth_hFac_rl = hFac
      ELSE
          smooth_hFac_rl = oneRL
      ENDIF

      RETURN
      END

C --------------------------

      _RL FUNCTION hFac_spline_rl( x, x0, eps )

      IMPLICIT NONE
#include "EEPARAMS.h"

      _RL x, x0, eps
      _RL a, b, d, del_x

      a = - oneRL / (eps**2)
      b = (2. _d 0) / eps
      d = x0

      del_x = x - x0

      hFac_spline_rl = a * (del_x**3) + b * (del_x**2) + d

      RETURN
      END


C --------------------------

      _RS FUNCTION SMOOTHMIN_RS( a, b )

      IMPLICIT NONE

      _RS a, b

      _RS SMOOTHABS_RS
      EXTERNAL SMOOTHABS_RS

      SMOOTHMIN_RS = .5*( a+b - SMOOTHABS_RS(a-b) )

      RETURN
      END

C --------------------------

      _RL FUNCTION SMOOTHMIN_RL( a, b )

      IMPLICIT NONE

      _RL a, b

      _RL SMOOTHABS_RL
      EXTERNAL SMOOTHABS_RL

      SMOOTHMIN_RL = .5*( a+b - SMOOTHABS_RL(a-b) )

      RETURN
      END

C --------------------------

      _RS FUNCTION SMOOTHMAX_RS( a, b )

      IMPLICIT NONE

      _RS a, b

      _RS SMOOTHABS_RS
      EXTERNAL SMOOTHABS_RS

      SMOOTHMAX_RS = .5*( a+b + SMOOTHABS_RS(a-b) )

      RETURN
      END

C --------------------------

      _RL FUNCTION SMOOTHMAX_RL( a, b )

      IMPLICIT NONE

      _RL a, b

      _RL SMOOTHABS_RL
      EXTERNAL SMOOTHABS_RL

      SMOOTHMAX_RL = .5*( a+b + SMOOTHABS_RL(a-b) )

      RETURN
      END

C --------------------------

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

C --------------------------

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

C --------------------------
#else

      _RL FUNCTION reg_hFac_rl( hFac, hFacMin )
C     o Impose 1,0, or hFacMin restriction 
C     o as in ini_masks_etc 

      IMPLICIT NONE
#include "EEPARAMS.h"

      _RL hFac, hFacMin

C     o Now apply min size constraint
      if ( hFac.lt.hFacMin ) then
       if ( hFac.lt.hFacMin*halfRL ) then
         reg_hFac_rl = zeroRL
       else
         reg_hFac_rl = hFacMin
       endif
      else
       reg_hFac_rl = hFac
      endif
      END

#endif /* USE_SMOOTH_MIN */
