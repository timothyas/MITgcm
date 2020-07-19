#include "SMOOTH_OPTIONS.h"

      subroutine smooth_apply_matern_2d (
     U     fld_in,smooth2DMask,smoothOpNb,doNormalize,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_apply_matern_2d
C     | o Routine that applies the following differential operator
C     |   to any field, u
C     |         
C     |         K^2 u - laplacian(u)
C     |         
C     | o This is used for the regularization term in an inverse 
C     |   problem
C     |
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
c#include "tamc.h"
#include "SMOOTH.h"

C --- Inputs
      _RL fld_in(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL smooth2DMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      integer smoothOpNb
      logical doNormalize
      integer myThid

C --- Local parameters
      _RL fld_in_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL center(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL east(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL west(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL north(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL south(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      character*( 80) fnamegeneric
      integer i,j,bi,bj
      integer itlo,ithi
      integer jtlo,jthi


C --- Loop counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

C --- Initialize
      do bj = jtlo,jthi
       do bi = itlo,ithi
        do j = 1,sNy
         do i = 1,sNx
          fld_out(i,j,bi,bj) = 0. _d 0
         enddo
        enddo
       enddo
      enddo

C --- Read in operator
      CALL smooth_read_operator_2d(smoothOpNb,myThid)

C --- Make sure input mask halos are good to go
      CALL EXCH_XY_RL( smooth2DMask, myThid )

C --- Get laplacian matrix elements
      CALL smooth_get_laplacian_2d(
     &  center, east, west, north, south,
     &  smoothOpNb, smooth2DMask, myThid)


c --- 1. Apply inverse of filter normalization (i.e. sqrt(filtervar))
      if (doNormalize) then
       do bj = jtlo,jthi
        do bi = itlo,ithi
         do j = 1,sNy
          do i = 1,sNx
           if ((ABS(smooth2DNorm(i,j,bi,bj)).gt.0. _d 0).and.
     &         (smooth2DMask(i,j,bi,bj).ne.0)) then
            fld_in(i,j,bi,bj)=
     &       fld_in(i,j,bi,bj)/smooth2DNorm(i,j,bi,bj)
           else
            fld_in(i,j,bi,bj)=0. _d 0
           endif
          enddo
         enddo
        enddo
       enddo
       CALL EXCH_XY_RL ( fld_in , myThid )
      endif

C --- 2. Apply laplacian-like operator
      do bj = jtlo,jthi
       do bi = itlo,ithi
        do j = 1,sNy
         do i = 1,sNx

          fld_out(i,j,bi,bj) = fld_in(i,j,bi,bj)*
     &         (smooth2DDelta(i,j,bi,bj) -
     &         center(i,j,bi,bj)) -
     &         east (i,j,bi,bj)*fld_in(i+1,j,bi,bj) -
     &         west (i,j,bi,bj)*fld_in(i-1,j,bi,bj) -
     &         north(i,j,bi,bj)*fld_in(i,j+1,bi,bj) -
     &         south(i,j,bi,bj)*fld_in(i,j-1,bi,bj)

          ! Is this necessary? I can't keep track anymore...
          fld_out(i,j,bi,bj) = fld_out(i,j,bi,bj)*
     &         smooth2DMask(i,j,bi,bj)
         enddo
        enddo
       enddo
      enddo

      CALL EXCH_XYZ_RL ( fld_out, myThid )

c --- 3. Rewrite input <- output
       do bj = jtlo,jthi
        do bi = itlo,ithi
         do j = 1-OLy,sNy+OLy
          do i = 1-OLx,sNx+OLx
           fld_in(i,j,bi,bj) = fld_out(i,j,bi,bj)
          enddo
         enddo
        enddo
       enddo

      end
