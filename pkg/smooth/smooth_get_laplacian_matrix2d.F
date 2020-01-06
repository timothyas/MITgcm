#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

      subroutine smooth_get_laplacian_matrix2d (
     U     center, north, south, east, west,
     I     smoothOpNb, smooth2DMask, mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_get_laplacian2d
C     | o get the matrix to the discretized 2D (horizontal) 
C     |   laplacian operator
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
#include "SMOOTH.h"

C --- Inputs
C   center, north, south, east, west :
C       fields that define the laplacian stencil at each point
C       for a point (i,j)
C           center = (i,j)
C           north = (i,j+1); south = (i,j-1)
C           east  = (i+1,j); west  = (i-1,j)
C
      _RL center(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL north(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL south(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL east(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL west(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      integer smoothOpNb
      _RL smooth2Dmask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      integer myThid

C --- Local parameters
      character*( 80) fnamegeneric
      integer i,j,bi,bj
      integer itlo,ithi
      integer jtlo,jthi

C --- Loop counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

C --- Diffusion operator (tensor): (Kux, 0; 0, Kvy)
C      write(fnamegeneric(1:80),'(1a,i3.3)')
C     &    'smooth2Doperator',smoothOpNb
C      CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
C     &           1, smooth2D_Kux,1,1,mythid)
C      CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
C     &           1, smooth2D_Kvy,2,1,mythid)
C      CALL EXCH_XY_RL ( smooth2D_Kux, myThid )
C      CALL EXCH_XY_RL ( smooth2D_Kvy, myThid )
    
C --- Initialize output
      DO bj = jtlo,jthi
       DO bi = itlo,ithi
        DO j = 1-OLy,sNy+OLy
         DO i = 1-OLx,sNx+OLx
          center(i,j,bi,bj)=0. _d 0
          north(i,j,bi,bj)=0. _d 0
          south(i,j,bi,bj)=0. _d 0
          east(i,j,bi,bj)=0. _d 0
          west(i,j,bi,bj)=0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C --- Get discretized Laplacian
C     Can operate on overlap region directly because
C     all input fields are "pre-exchanged"
C     Note: sign is reversed from smooth_diff2d!
      DO bj = jtlo,jthi
       DO bi = itlo,ithi
        DO j = 1,sNy
         DO i = 1,sNx
          center(i,j,bi,bj)=center(i,j,bi,bj) -
     &     ( smooth2D_Kux(i,j,bi,bj)*
     &       dyG(i,j,bi,bj)*recip_dxC(i,j,bi,bj)*
     &       smooth2Dmask(i,j,bi,bj)*smooth2Dmask(i-1,j,bi,bj) ) -
     &     ( smooth2D_Kux(i+1,j,bi,bj)*
     &       dyG(i+1,j,bi,bj)*recip_dxC(i+1,j,bi,bj)*
     &       smooth2Dmask(i,j,bi,bj)*smooth2Dmask(i+1,j,bi,bj) ) -
     &     ( smooth2D_Kvy(i,j,bi,bj)*
     &       dxG(i,j,bi,bj)*recip_dyC(i,j,bi,bj)*
     &       smooth2Dmask(i,j,bi,bj)*smooth2Dmask(i,j-1,bi,bj) ) -
     &     ( smooth2D_Kvy(i,j+1,bi,bj)*
     &       dxG(i,j+1,bi,bj)*recip_dyC(i,j+1,bi,bj)*
     &       smooth2Dmask(i,j,bi,bj)*smooth2Dmask(i,j+1,bi,bj) )

          west(i,j,bi,bj)=west(i,j,bi,bj) +
     &     smooth2D_Kux(i,j,bi,bj)*
     &     dyG(i,j,bi,bj)*recip_dxC(i,j,bi,bj)*
     &     smooth2Dmask(i,j,bi,bj)*smooth2Dmask(i-1,j,bi,bj)

          east(i,j,bi,bj)=east(i,j,bi,bj) +
     &     smooth2D_Kux(i+1,j,bi,bj)*
     &     dyG(i+1,j,bi,bj)*recip_dxC(i+1,j,bi,bj)*
     &     smooth2Dmask(i,j,bi,bj)*smooth2Dmask(i+1,j,bi,bj)

          south(i,j,bi,bj)=south(i,j,bi,bj) +
     &     smooth2D_Kvy(i,j,bi,bj)*
     &     dxG(i,j,bi,bj)*recip_dyC(i,j,bi,bj)*
     &     smooth2Dmask(i,j,bi,bj)*smooth2Dmask(i,j-1,bi,bj)

          north(i,j,bi,bj)=north(i,j,bi,bj) +
     &     smooth2D_Kvy(i,j+1,bi,bj)*
     &     dxG(i,j+1,bi,bj)*recip_dyC(i,j+1,bi,bj)*
     &     smooth2Dmask(i,j,bi,bj)*smooth2Dmask(i,j+1,bi,bj)
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CALL EXCH_XY_RL ( center, myThid )
      CALL EXCH_XY_RL ( north, myThid )
      CALL EXCH_XY_RL ( south, myThid )
      CALL EXCH_XY_RL ( east, myThid )
      CALL EXCH_XY_RL ( west, myThid )

      end
