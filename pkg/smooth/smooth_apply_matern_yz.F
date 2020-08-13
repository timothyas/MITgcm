#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#endif

      subroutine smooth_apply_matern_yz (
     U     fld_in,smooth2DMask,smoothOpNb,doNormalize,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_apply_matern_yz
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
#ifdef ALLOW_OBCS
# include "OBCS_GRID.h"
#endif

C --- Inputs
      _RL fld_in(1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL smooth2DMask(1-OLy:sNy+OLy,Nr,nSx,nSy)
      integer smoothOpNb
      logical doNormalize
      integer myThid

#ifdef ALLOW_OBCS
# if (defined ALLOW_OBCS_EAST) || (defined ALLOW_OBCS_WEST)

C --- Local parameters
      _RL fld_out(1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL center (1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north  (1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south  (1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL up     (1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL down   (1-OLy:sNy+OLy,Nr,nSx,nSy)
      character*( 80) fnamegeneric
      integer i,j,k,bi,bj
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
        do k = 1,Nr
         do j = 1-OLy,sNy+OLy
          fld_out(j,k,bi,bj) = 0. _d 0
         enddo
        enddo
       enddo
      enddo

C --- Read in operator
      CALL smooth_read_operator_yz(smoothOpNb,myThid)

C --- Make sure input mask halos are good to go
      CALL SMOOTH_EXCH_YZ_RL( smooth2DMask, smoothOpNb, myThid )

C --- Get laplacian matrix elements
      CALL smooth_get_laplacian_yz(
     &  center, north, south, up, down,
     &  smoothOpNb, smooth2DMask, myThid)


c --- 1. Apply inverse of filter normalization (i.e. sqrt(filtervar))
      if (doNormalize) then
       do bj = jtlo,jthi
        do bi = itlo,ithi
         IF (((smooth2DDims(smoothOpNb).EQ.'yze').AND.
     &        (tileHasOBE(bi,bj))).OR.
     &       ((smooth2DDims(smoothOpNb).EQ.'yzw').AND.
     &        (tileHasOBW(bi,bj)))) THEN
          do k = 1,Nr
           do j = 1-OLy,sNy+OLy

            IF ((ABS(smooth2DNorm(j,k,bi,bj)).GT.0. _d 0).AND.
     &          (smooth2DMask(j,k,bi,bj).NE.0)) then

              fld_in(j,k,bi,bj) = fld_in(j,k,bi,bj)/
     &            smoothYZNorm(j,k,bi,bj)
            ELSE
             fld_in(j,k,bi,bj) = 0. _d 0
            ENDIF
           enddo
          enddo
         ENDIF
        enddo
       enddo
      endif

C --- 2. Apply laplacian-like operator
      do bj = jtlo,jthi
       do bi = itlo,ithi
         IF (((smooth2DDims(smoothOpNb).EQ.'yze').AND.
     &        (tileHasOBE(bi,bj))).OR.
     &       ((smooth2DDims(smoothOpNb).EQ.'yzw').AND.
     &        (tileHasOBW(bi,bj)))) THEN
          do k = 1,Nr
           do j = 1,sNy

            fld_out(j,k,bi,bj) = fld_in(j,k,bi,bj)*
     &         (smooth2DDelta(j,k,bi,bj) - center(j,k,bi,bj)) -
     &         north (j,k,bi,bj)*fld_in(j+1,k,bi,bj) -
     &         south (j,k,bi,bj)*fld_in(j-1,k,bi,bj)

             if (k.lt.Nr)
     &         fld_out(j,k,bi,bj) = fld_out(j,k,bi,bj) -
     &          up(j,k,bi,bj)*fld_in(j,k+1,bi,bj)

             if (k.gt.1)
     &         fld_out(j,k,bi,bj) = fld_out(j,k,bi,bj) -
     &          down(j,k,bi,bj)*fld_in(j,k-1,bi,bj)

          enddo
         enddo
        ENDIF
       enddo
      enddo

      CALL SMOOTH_EXCH_YZ_RL ( fld_out, smoothOpNb, myThid )

c --- 3. Rewrite input <- output
       do bj = jtlo,jthi
        do bi = itlo,ithi
         do k = 1,Nr
          do j = 1-OLy,sNy+OLy
           fld_in(j,k,bi,bj) = fld_out(j,k,bi,bj)
          enddo
         enddo
        enddo
       enddo

# endif
#endif
      end
