#include "SMOOTH_OPTIONS.h"

      subroutine smooth_correl3D (
     U     fld_in,smoothOpNb,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_correl3D
C     | o Routine that applies spatial correlation 
C     |   operator to a 3D control field
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
c#include "tamc.h"
#include "SMOOTH.h"


      _RL fld_in(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy)
      integer smoothOpNb
      integer nbt_in
      character*( 80) fnamegeneric
      integer i,j,k,bi,bj
      integer itlo,ithi
      integer jtlo,jthi
      integer myThid


      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

c division by ~sqrt(volume):
      DO bj=jtlo,jthi
       DO bi=itlo,ithi
        DO k=1,Nr
         DO j=1,sNy
          DO i=1,sNx 
           fld_in(i,j,k,bi,bj)=fld_in(i,j,k,bi,bj)
     & *sqrt(recip_rA(i,j,bi,bj)*recip_drF(k))
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
      CALL EXCH_XYZ_RL ( fld_in , myThid )

c do the smoothing:
      nbt_in=smooth3Dnbt(smoothOpNb)/2
      call smooth_diff3D(fld_in,nbt_in,mythid)

c division by ~sqrt(var(filter)):
       do bj = jtlo,jthi
        do bi = itlo,ithi
         DO j = 1,sNy
          DO i = 1,sNx
           DO k = 1,nR
       fld_in(i,j,k,bi,bj)=fld_in(i,j,k,bi,bj)
     & *smooth3Dnorm(i,j,k,bi,bj)
           ENDDO
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      CALL EXCH_XYZ_RL ( fld_in , myThid )

      end
