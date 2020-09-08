#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

      subroutine smooth_get_laplacian_2d (
     U     center, east, west, north, south,
     I     smoothOpNb, smooth2DMask, mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_get_laplacian_2d
C     | o get the matrix to the discretized 2D 
C     |   laplacian operator: div( kappa * grad( ) )
C     | o only the diagonal elements of the diffusion tensor:
C     |
C     |             Kux Kuy    Kux  0 
C     |    kappa =  Kvx Kvy ->  0  Kvy
C     |
C     |   are implemented currently
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
#include "SMOOTH.h"

C --- Inputs
C   center, north, south, east, west, up, down :
C       fields that define the laplacian stencil at each point
C       for a point (i,j)
C           center = (i,j)
C           east  = (i+1, j  ); west  = (i-1, j  )
C           north = (i  , j+1); south = (i  , j-1)
      _RL center(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL east  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL west  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL north (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL south (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)

      integer smoothOpNb
      _RL smooth2Dmask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      integer myThid

C --- Local parameters
      _RL xA(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL yA(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      character*( 80) fnamegeneric
      integer i,j,bi,bj
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
          center(i,j,bi,bj)=0. _d 0
          east  (i,j,bi,bj)=0. _d 0
          west  (i,j,bi,bj)=0. _d 0
          north (i,j,bi,bj)=0. _d 0
          south (i,j,bi,bj)=0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C --- Get discretized Laplacian
C     Note: these are the (swapped sign of) the diagonal
C           elements in applied to input field in smooth_rhs
      DO bj = jtlo,jthi
       DO bi = itlo,ithi
        DO j = 1-OLy,sNy+OLy
         DO i = 1-OLx,sNx+OLx

           xA(i,j,bi,bj) = _dyG(i,j,bi,bj) *
     &         _recip_dxC(i,j,bi,bj) * smooth2D_Kux(i,j,bi,bj)

           yA(i,j,bi,bj) = _dxG(i,j,bi,bj) *
     &         _recip_dyC(i,j,bi,bj) * smooth2D_Kvy(i,j,bi,bj)

         ENDDO
        ENDDO

        DO j = 1-OLy+1,sNy+OLy-1
         DO i = 1-OLx+1,sNx+OLx-1

          east(i,j,bi,bj)=east(i,j,bi,bj) + xA(i+1,j,bi,bj)*
     &     smooth2Dmask(i,j,bi,bj)*smooth2Dmask(i+1,j,bi,bj)*
     &     recip_rA(i,j,bi,bj)

          west(i,j,bi,bj)=west(i,j,bi,bj) + xA(i,j,bi,bj)*
     &     smooth2Dmask(i-1,j,bi,bj)*smooth2Dmask(i,j,bi,bj)*
     &     recip_rA(i,j,bi,bj)

          north(i,j,bi,bj)=north(i,j,bi,bj) + yA(i,j+1,bi,bj)*
     &     smooth2Dmask(i,j,bi,bj)*smooth2Dmask(i,j+1,bi,bj)*
     &     recip_rA(i,j,bi,bj)

          south(i,j,bi,bj)=south(i,j,bi,bj) + yA(i,j,bi,bj)*
     &     smooth2Dmask(i,j-1,bi,bj)*smooth2Dmask(i,j,bi,bj)*
     &     recip_rA(i,j,bi,bj)

          center(i,j,bi,bj) = center(i,j,bi,bj) -
     &      east(i,j,bi,bj) - west(i,j,bi,bj) -
     &      north(i,j,bi,bj) - south(i,j,bi,bj)

         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CALL EXCH_XY_RL( center, myThid )
      CALL EXCH_XY_RL( east, myThid )
      CALL EXCH_XY_RL( west, myThid )
      CALL EXCH_XY_RL( north, myThid )
      CALL EXCH_XY_RL( south, myThid )

      end
