#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#endif

      subroutine smooth_inverse_matern_xz (
     U     fld_in,
     I     smooth2DMask,smoothOpNb,doNormalize,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_inverse_matern_xz
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
C     |   which is written to file as smoothXZnorm00<smoothOpNb>
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
      _RL fld_in      (1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL smooth2DMask(1-OLx:sNx+OLx,Nr,nSx,nSy)
      integer smoothOpNb
      logical doNormalize
      integer myThid

#ifdef ALLOW_OBCS
# if (defined ALLOW_OBCS_NORTH) || (defined ALLOW_OBCS_SOUTH)

C --- Local parameters
      character*( 80) fnamegeneric
      integer i,j,k,bi,bj
      integer itlo,ithi
      integer jtlo,jthi

      _RL center(1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL east  (1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL west  (1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL up    (1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL down  (1-OLx:sNx+OLx,Nr,nSx,nSy)

C --- Loop Counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

C --- Update halo retions for mask
C     this is not done right now by ctrl_map_ini_genarr
C     so we can't assume it to be the case
      CALL SMOOTH_EXCH_XZ_RL( smooth2DMask, smoothOpNb, myThid )
      CALL SMOOTH_EXCH_XZ_RL( fld_in, smoothOpNb, myThid )

C ------------------------------------------------------------------
C the actual algorithm
C ------------------------------------------------------------------

c --- 0. Apply inv sqrt of cell volume and determinant of mapping:
      if (doNormalize) then
       do bj = jtlo,jthi
        do bi = itlo,ithi
         IF (((smooth2DDims(smoothOpNb).EQ.'xzn').AND.
     &        (tileHasOBN(bi,bj))).OR.
     &       ((smooth2DDims(smoothOpNb).EQ.'xzs').AND.
     &        (tileHasOBS(bi,bj)))) then
          do k = 1,Nr
           do i = 1-OLx,sNx+OLx
           fld_in(i,k,bi,bj) = fld_in(i,k,bi,bj)*
     &         smoothXZRandNorm(i,k,bi,bj)
           enddo
          enddo
         ENDIF
        enddo
       enddo
      endif

c --- 1. Get laplacian matrix elements
      call smooth_get_laplacian_xz(
     &      center, east, west, up, down,
     &      smoothOpNb, smooth2DMask, myThid)

c --- 2. Create A matrix associated with operator A = (delta - nabla^2)
c        operate directly on overlap region, inputs "pre-exchanged"
      do bj = jtlo,jthi
       do bi = itlo,ithi
        do k = 1,Nr
         do i = 1-OLx,sNx+OLx
          center(i,k,bi,bj) = -center(i,k,bi,bj) +
     &      smoothXZDelta(i,k,bi,bj)
          east (i,k,bi,bj)  = -east(i,k,bi,bj)
          west (i,k,bi,bj)  = -west(i,k,bi,bj)
          up   (i,k,bi,bj)  = -up  (i,k,bi,bj)
          down (i,k,bi,bj)  = -down(i,k,bi,bj)
         enddo
        enddo
       enddo
      enddo

c --- 3. Get jacobi to solve our problem
      call smooth_jacobi_xz( fld_in, smoothOpNb,
     &      center, east, west, up, down,
     &      smooth2DMask, myThid )

c --- 4. Apply normalization field [i.e. 1/sqrt(var(filter))]:
      if (doNormalize) then
       do bj = jtlo,jthi
        do bi = itlo,ithi
         IF (((smooth2DDims(smoothOpNb).EQ.'xzn').AND.
     &        (tileHasOBN(bi,bj))).OR.
     &       ((smooth2DDims(smoothOpNb).EQ.'xzs').AND.
     &        (tileHasOBS(bi,bj)))) then
          do k = 1,Nr
           do i = 1-OLx,sNx+OLx
           fld_in(i,k,bi,bj) = fld_in(i,k,bi,bj)*
     &         smoothXZNorm(i,k,bi,bj)
           enddo
          enddo
         ENDIF
        enddo
       enddo
      endif

# endif
#endif
      end
