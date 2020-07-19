#include "SMOOTH_OPTIONS.h"

      subroutine smooth_correl2D (
     U     fld_in,mask_in,smoothOpNb,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_correl2D
C     | o Routine that applies spatial correlation 
C     |   operator to a 2D control field
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
c#include "tamc.h"
#include "SMOOTH.h"

      _RL mask_in(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL fld_in(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      integer smoothOpNb
      integer nbt_in
      character*( 80) fnamegeneric
      integer i,j,bi,bj
      integer itlo,ithi
      integer jtlo,jthi
      integer myThid


      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)


c read smoothing [i.e diffusion] operator:

c division by ~sqrt(area):
      DO bj = jtlo,jthi
       DO bi = itlo,ithi
        DO j = 1,sNy
         DO i = 1,sNx 
      fld_in(i,j,bi,bj)=fld_in(i,j,bi,bj)
     & *sqrt(recip_rA(i,j,bi,bj))
         ENDDO
        ENDDO
       ENDDO 
      ENDDO
      CALL EXCH_XY_RL ( fld_in , myThid )

c do the smoothing:
      nbt_in=smooth2Dnbt(smoothOpNb)/2
      call smooth_diff2D(fld_in,mask_in,nbt_in,mythid)

c division by ~sqrt(var(filter)):
       do bj = jtlo,jthi
        do bi = itlo,ithi
         DO j = 1,sNy
          DO i = 1,sNx
       fld_in(i,j,bi,bj)=fld_in(i,j,bi,bj)
     & *smooth2Dnorm(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      CALL EXCH_XY_RL ( fld_in , myThid )
       
      end
