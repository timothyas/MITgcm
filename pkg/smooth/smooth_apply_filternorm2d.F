#include "SMOOTH_OPTIONS.h"

      subroutine smooth_apply_filternorm2d (
     U     fld_in,smoothOpNb,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_apply_filternorm2d
C     | o Simply read in and apply the normalization field 
C     |   to an input field
C     |
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
#include "SMOOTH.h"

C --- Inputs
      _RL fld_in(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      integer smoothOpNb
      integer myThid

C --- Local parameters
      character*( 80) fnamegeneric
      integer i,j,bi,bj
      integer itlo,ithi
      integer jtlo,jthi

      _RL filtervar_norm(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

C --- Loop Counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

c --- read normalization field [i.e. 1/sqrt(var(filter))]:
      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth2Dnorm',smoothOpNb
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
     &           1, filtervar_norm,1,1,mythid)
      CALL EXCH_XY_RL ( filtervar_norm, myThid )

c --- Do the normalization
      do bj = jtlo,jthi
       do bi = itlo,ithi
        do j = 1,sNy
         do i = 1,sNx
          fld_in(i,j,bi,bj)=
     &      fld_in(i,j,bi,bj)*filtervar_norm(i,j,bi,bj)
         enddo
        enddo
       enddo
      enddo
      CALL EXCH_XY_RL ( fld_in , myThid )
       
      end
