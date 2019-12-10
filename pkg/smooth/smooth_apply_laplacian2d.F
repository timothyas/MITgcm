#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

      subroutine smooth_apply_laplacian2d (
     U     fld_in,
     I     smooth2Dmask,
     I     mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_apply_laplacian2d
C     | o Apply the laplacian operator ONCE to the input field
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
#include "SMOOTH.h"

C --- Inputs
      _RL fld_in(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL smooth2Dmask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      integer myThid

C --- Local parameters
      _RL fld_out(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      integer i,j,k, bi, bj
      integer itlo,ithi
      integer jtlo,jthi

C --- Loop counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)
    
C --- Initialize output
      DO bj = jtlo,jthi
       DO bi = itlo,ithi
        DO j = 1-OLy,sNy+OLy
         DO i = 1-OLx,sNx+OLx
          fld_out(i,j,bi,bj)=0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C --- Apply discretized Laplacian
C     Note: sign is reversed from smooth_diff2d!
      DO bj = jtlo,jthi
       DO bi = itlo,ithi
        DO j = 1,sNy
         DO i = 1,sNx
          if (smooth2Dmask(i,j,bi,bj).NE.0.) then
      fld_out(i,j,bi,bj)=fld_out(i,j,bi,bj)-
     & dyG(i,j,bi,bj)*
     & smooth2Dmask(i,j,bi,bj)*smooth2Dmask(i-1,j,bi,bj)*
     & (fld_in(i,j,bi,bj)-fld_in(i-1,j,bi,bj))*recip_dxC(i,j,bi,bj)

      fld_out(i,j,bi,bj)=fld_out(i,j,bi,bj)-
     & dyG(i+1,j,bi,bj)*
     & smooth2Dmask(i,j,bi,bj)*smooth2Dmask(i+1,j,bi,bj)*
     & (fld_in(i,j,bi,bj)-fld_in(i+1,j,bi,bj))*recip_dxC(i+1,j,bi,bj)

      fld_out(i,j,bi,bj)=fld_out(i,j,bi,bj)-
     & dxG(i,j,bi,bj)*
     & smooth2Dmask(i,j,bi,bj)*smooth2Dmask(i,j-1,bi,bj)*
     & (fld_in(i,j,bi,bj)-fld_in(i,j-1,bi,bj))*recip_dyC(i,j,bi,bj)

      fld_out(i,j,bi,bj)=fld_out(i,j,bi,bj)-
     & dxG(i,j+1,bi,bj)*
     & smooth2Dmask(i,j,bi,bj)*smooth2Dmask(i,j+1,bi,bj)*
     & (fld_in(i,j,bi,bj)-fld_in(i,j+1,bi,bj))*recip_dyC(i,j+1,bi,bj)

          endif
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C --- Replace fld_in with fld_out to return the result
      DO bj = jtlo,jthi
       DO bi = itlo,ithi
        DO j = 1,sNy
         DO i = 1,sNx
          fld_in(i,j,bi,bj)=fld_out(i,j,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CALL EXCH_XY_RL ( fld_in , myThid )

      end
