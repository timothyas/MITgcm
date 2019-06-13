#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: INI_SW_DEPTHS
C     !INTERFACE:
      SUBROUTINE INI_SW_DEPTHS( myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | SUBROUTINE INI_SW_DEPTHS
C     | o Initialise South & West cell boundary values for 
C     |   rLow and rSurf
C     | o Location is key with depth control
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"

C     !INPUT/OUTPUT PARAMETERS:
C     myThid    ::  my Thread Id number
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
      INTEGER bi, bj
      INTEGER i, j, k

cts ---
C   Need to re-do initialization of rLowW/S & rSurfW/S 
C   for correct computation later in
C   the code, which is necessary because these values depend on r_low,
C   which is updated by the previous call to ctrl_depth_init
C   Note: This is particularly important for when useShelfice=.true. but
C   is kept here in general since the computation will be equivalent even
C   when this is false.
Cts ---

C   Initialise rLow & reference rSurf at Western & Southern edges (U & V pts)
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
       ENDDO
      ENDDO

      RETURN
      END

