#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#endif

      subroutine smooth_apply_matern_xz (
     U      fld_in,
     I      smooth2DMask,smoothOpNb,doNormalize,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_apply_matern_xz
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
#include "SMOOTH.h"
#ifdef ALLOW_OBCS
# include "OBCS_GRID.h"
#endif

C --- Inputs
      _RL fld_in(1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL smooth2DMask(1-OLx:sNx+OLx,Nr,nSx,nSy)
      integer smoothOpNb
      logical doNormalize
      integer myThid

#ifdef ALLOW_OBCS
# if (defined ALLOW_OBCS_NORTH) || (defined ALLOW_OBCS_SOUTH)

C --- Local parameters
      _RL fld_out(1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL center (1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL east   (1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL west   (1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL up     (1-OLx:sNx+OLx,Nr,nSx,nSy)
      _RL down   (1-OLx:sNx+OLx,Nr,nSx,nSy)
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
         do i = 1-OLx,sNx+OLx
          fld_out(i,k,bi,bj) = 0. _d 0
         enddo
        enddo
       enddo
      enddo

C --- Read in operator
      CALL smooth_read_operator_xz( smoothOpNb, myThid )

C --- Make sure input mask halos are good to go
      CALL SMOOTH_EXCH_XZ_RL( smooth2DMask, smoothOpNb, myThid )
      CALL SMOOTH_EXCH_XZ_RL( fld_in, smoothOpNb, myThid )

C --- Get laplacian matrix elements
      CALL smooth_get_laplacian_xz(
     &  center, east, west, up, down,
     &  smoothOpNb, smooth2DMask, myThid)


c --- 1. Apply inverse of filter normalization (i.e. sqrt(filtervar))
      if (doNormalize) then
       do bj = jtlo,jthi
        do bi = itlo,ithi
         IF (((smooth2DDims(smoothOpNb).EQ.'xzn').AND.
     &        (tileHasOBN(bi,bj))).OR.
     &       ((smooth2DDims(smoothOpNb).EQ.'xzs').AND.
     &        (tileHasOBS(bi,bj)))) THEN
          do k = 1,Nr
           do i = 1-OLx,sNx+OLx

            IF ((ABS(smooth2DNorm(i,k,bi,bj)).GT.0. _d 0).AND.
     &          (smooth2DMask(i,k,bi,bj).NE.0)) then

              fld_in(i,k,bi,bj) = fld_in(i,k,bi,bj)/
     &            smoothXZNorm(i,k,bi,bj)
            ELSE
             fld_in(i,k,bi,bj) = 0. _d 0
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
         IF (((smooth2DDims(smoothOpNb).EQ.'xzn').AND.
     &        (tileHasOBN(bi,bj))).OR.
     &       ((smooth2DDims(smoothOpNb).EQ.'xzs').AND.
     &        (tileHasOBS(bi,bj)))) THEN
          do k = 1,Nr
           do i = 1,sNx

            fld_out(i,k,bi,bj) = fld_in(i,k,bi,bj)*
     &         (smooth2DDelta(i,k,bi,bj) - center(i,k,bi,bj)) -
     &         east (i,k,bi,bj)*fld_in(i+1,k,bi,bj) -
     &         west (i,k,bi,bj)*fld_in(i-1,k,bi,bj)

             if (k.lt.Nr)
     &         fld_out(i,k,bi,bj) = fld_out(i,k,bi,bj) -
     &          up(i,k,bi,bj)*fld_in(i,k+1,bi,bj)

             if (k.gt.1)
     &         fld_out(i,k,bi,bj) = fld_out(i,k,bi,bj) -
     &          down(i,k,bi,bj)*fld_in(i,k-1,bi,bj)

          enddo
         enddo
        ENDIF
       enddo
      enddo

      CALL SMOOTH_EXCH_XZ_RL ( fld_out, smoothOpNb, myThid )

c --- 3. Rewrite input <- output
C        and divide by randnorm
       do bj = jtlo,jthi
        do bi = itlo,ithi
         do k = 1,Nr
          do i = 1-OLx,sNx+OLx
           fld_in(i,k,bi,bj) = fld_out(i,k,bi,bj)
           if (doNormalize) then
            if (ABS(smoothXZRandNorm(i,k,bi,bj)).gt.0.) then
             fld_in(i,k,bi,bj)=
     &       fld_in(i,k,bi,bj)/smoothXZRandNorm(i,k,bi,bj)
            else
             fld_in(i,k,bi,bj)= 0. _d 0
            endif
           endif
          enddo
         enddo
        enddo
       enddo

# endif
#endif
      end
