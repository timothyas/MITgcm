#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#endif

      subroutine smooth_exch_yz_rl (fld_yz, smoothOpNb, mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_exch_yz_rl
C     | o There are no 2D slice exchange routines, so this is a hack
C     |   to do that...
C     |
C     | o Fill a dummy 3D field, run the exchanges, and refill the
C     |   2D slice
C     |
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
#include "SMOOTH.h"
#ifdef ALLOW_OBCS
# include "OBCS_GRID.h"
#endif

C --- Inputs
      _RL fld_yz(1-OLy:sNy+OLy,Nr,nSx,nSy)
      integer smoothOpNb
      integer myThid

#ifdef ALLOW_OBCS
# if (defined ALLOW_OBCS_EAST) || (defined ALLOW_OBCS_WEST)

C --- Local parameters
      character*( 80) fnamegeneric
      character*(MAX_LEN_MBUF) msgBuf
      _RL dummy3D(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      integer i,j,k,bi,bj
      integer itlo,ithi
      integer jtlo,jthi

C --- Loop Counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

C --- 1. Fill the dummy
      DO bj = jtlo,jthi
       DO bi = itlo,ithi

C ---   Initialize
        DO k = 1,Nr
         DO j = 1-OLy,sNy+OLy
          DO i = 1-OLx,sNx+OLx
           dummy3D(i,j,k,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
        ENDDO

C ---   Eastside
        IF ((smooth2DDims(smoothOpNb).EQ.'yze').AND.
     &      (tileHasOBE(bi,bj))) then
         DO j = 1,sNy
          i = OB_Ie(j,bi,bj)
          DO k = 1,Nr
           dummy3D(i,j,k,bi,bj) = fld_yz(j,k,bi,bj)
          ENDDO
         ENDDO

C ---   Westside
        ELSEIF ((smooth2DDims(smoothOpNb).EQ.'yzw').AND.
     &          (tileHasOBW(bi,bj))) then
         DO j = 1,sNy
          i = OB_Iw(j,bi,bj)
          DO k = 1,Nr
           dummy3D(i,j,k,bi,bj) = fld_yz(j,k,bi,bj)
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO

C --- 2. Exchange the dummy
      CALL EXCH_XYZ_RL( dummy3D, myThid )

C --- 3. Refill the input, including overlaps
      DO bj = jtlo,jthi
       DO bi = itlo,ithi

C ---   Eastside
        IF ((smooth2DDims(smoothOpNb).EQ.'yze').AND.
     &      (tileHasOBE(bi,bj))) then
         DO j = 1-OLy,sNy+OLy
          i = OB_Ie(j,bi,bj)
          DO k = 1,Nr
           fld_yz(j,k,bi,bj) = dummy3D(i,j,k,bi,bj)
          ENDDO
         ENDDO

C ---   Westside
        ELSEIF ((smooth2DDims(smoothOpNb).EQ.'yzw').AND.
     &          (tileHasOBW(bi,bj))) then
         DO j = 1-OLy,sNy+OLy
          i = OB_Iw(j,bi,bj)
          DO k = 1,Nr
           fld_yz(j,k,bi,bj) = dummy3D(i,j,k,bi,bj)
          ENDDO
         ENDDO
        ELSE
         DO k = 1,Nr
          DO j = 1-OLy,sNy+OLy
           fld_yz(j,k,bi,bj) = 0. _d 0
          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO

# endif
#endif

      RETURN
      END
