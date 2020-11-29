#include "SMOOTH_OPTIONS.h"

      subroutine smooth_inverse_matern_3d (
     U     fld_in,
     I     smooth3DMask,smoothOpNb,doNormalize,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_inverse_matern_3d
C     | o Given the laplace-like differential operator
C     |
C     |         A=( delta - div( kappa * grad() )
C     |
C     |   this employs an iterative method to find the solution
C     |
C     |         u = A^{-1}b
C     |
C     | o This is used as a preconditioner to get an estimate
C     |   of the prior preconditioned misfit Hessian (Gauss Newton
C     |   approximation)
C     |
C     | o if doNormalize==True, result is multiplied by 1/sqrt(var(A^{-1}))
C     |   which is written to file as smooth2Dnorm00<smoothOpNb>
C     |
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
#include "SMOOTH.h"

C --- Inputs
      _RL fld_in      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL smooth3DMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      integer smoothOpNb
      logical doNormalize
      integer myThid

C --- Local parameters
      character*( 80) fnamegeneric
      integer i,j,k,bi,bj
      integer itlo,ithi
      integer jtlo,jthi

      _RL center(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL east  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL west  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL up    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL down  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)


C --- Loop Counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

C --- Update halo retions for mask
C     this is not done right now by ctrl_map_ini_genarr
C     so we can't assume it to be the case
      CALL EXCH_XYZ_RL( smooth3DMask, myThid )
      CALL EXCH_XYZ_RL( fld_in, myThid )

C ------------------------------------------------------------------
C the actual algorithm
C ------------------------------------------------------------------

c --- 0. Apply inv sqrt of cell volume and determinant of mapping:
      if (doNormalize) then
       do bj = jtlo,jthi
        do bi = itlo,ithi
         do k = 1,Nr
          do j = 1-OLy,sNy+OLy
           do i = 1-OLx,sNx+OLx
            fld_in(i,j,k,bi,bj) = fld_in(i,j,k,bi,bj)*
     &          smooth3DRandNorm(i,j,k,bi,bj)
           enddo
          enddo
         enddo
        enddo
       enddo
      endif

c --- 1. Get laplacian matrix elements
      call smooth_get_laplacian_3d(
     &      center, east, west, north, south, up, down,
     &      smoothOpNb, smooth3DMask, myThid)

c --- 2. Create A matrix associated with operator A = (delta - nabla^2)
c        operate directly on overlap region, inputs "pre-exchanged"
      do bj = jtlo,jthi
       do bi = itlo,ithi
        do k = 1,Nr
         do j = 1-OLy,sNy+OLy
          do i = 1-OLx,sNx+OLx
           center(i,j,k,bi,bj) = -center(i,j,k,bi,bj) +
     &          smooth3DDelta(i,j,k,bi,bj)
           east (i,j,k,bi,bj)  = -east (i,j,k,bi,bj)
           west (i,j,k,bi,bj)  = -west (i,j,k,bi,bj)
           north(i,j,k,bi,bj)  = -north(i,j,k,bi,bj)
           south(i,j,k,bi,bj)  = -south(i,j,k,bi,bj)
           up   (i,j,k,bi,bj)  = -up   (i,j,k,bi,bj)
           down (i,j,k,bi,bj)  = -down (i,j,k,bi,bj)
          enddo
         enddo
        enddo
       enddo
      enddo

c --- 3. Get jacobi to solve our problem
      call smooth_jacobi_3d( fld_in, smoothOpNb,
     &      center, east, west, north, south, up, down,
     &      smooth3DMask, myThid )

c --- 4. Apply normalization field [i.e. 1/sqrt(var(filter))]:
      if (doNormalize) then
       do bj = jtlo,jthi
        do bi = itlo,ithi
         do k = 1,Nr
          do j = 1-OLy,sNy+OLy
           do i = 1-OLx,sNx+OLx
            fld_in(i,j,k,bi,bj) = fld_in(i,j,k,bi,bj)*
     &          smooth3DNorm(i,j,k,bi,bj)
           enddo
          enddo
         enddo
        enddo
       enddo
      endif

      end
