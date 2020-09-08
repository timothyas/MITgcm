#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#endif

      subroutine smooth_jacobi_yz (
     U     fld_in,
     I     smoothOpNb,
     I     center, north, south, up, down,
     I     smooth2Dmask,
     I     mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_jacobi_yz
C     | o Invert matrix A via jacobi iterations to solve for u,
C     |   given right hand side b:
C     |          Au = b
C     |
C     | o matrix A is "Laplacian like" and defined by 
C     |   center, north, south, up, down fields
C     |
C     | o b := fld_in -> fld_rhs
C     |   u:= fld_out (fld_in is overwritten with this result at end)
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
C   center, north, south, east, west, up, down :
C       fields that define the laplacian stencil at each point
C       for a point (j,k)
C           center = (j,k)
C           north = (j+1,  k  ); south = (j-1, k  )
C           up    = (j  ,  k+1); down  = (j  , k-1)
C
      _RL fld_in(1-OLy:sNy+OLy,Nr,nSx,nSy)
      integer smoothOpNb
      _RL center(1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north (1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south (1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL up    (1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL down  (1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL smooth2Dmask(1-OLy:sNy+OLy,Nr,nSx,nSy)
      integer myThid

#ifdef ALLOW_OBCS
# if (defined ALLOW_OBCS_EAST) || (defined ALLOW_OBCS_WEST)

C --- Local parameters
      integer i,j,k, bi, bj
      integer itlo,ithi
      integer jtlo,jthi

      _RL TOL, TOL_SQ, fld_in_norm, relErr, absErr
      _RL sor
      integer jacobi_iter

      _RL fld_rhs(1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL fld_out(1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL recip_center(1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL tile_in_norm(nSx,nSy)
      _RL tile_absErr(nSx,nSy)
      character(max_len_mbuf) msgBuf

C --- Loop counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

C --- Initialize solver stuff
      jacobi_iter = 0
      TOL = 1.0 _d -15
      TOL_SQ = TOL*TOL
      relErr = 100. _d 0
      absErr = 0. _d 0 
      sor = smooth2DSOROmega(smoothOpNb)
      
      do bj = jtlo,jthi
       do bi = itlo,ithi
        IF (((smooth2DDims(smoothOpNb).EQ.'yze').AND.
     &       (tileHasOBE(bi,bj))).OR.
     &      ((smooth2DDims(smoothOpNb).EQ.'yzw').AND.
     &       (tileHasOBW(bi,bj)))) then
         DO k = 1,Nr
          DO j = 1-OLy,sNy+OLy
           fld_rhs(j,k,bi,bj) = fld_in(j,k,bi,bj)
           fld_in (j,k,bi,bj) = 1. _d 0
           fld_out(j,k,bi,bj) = 1. _d 0
           if (ABS(center(j,k,bi,bj)).gt.0. _d 0) then
            recip_center(j,k,bi,bj) = 1. _d 0 / center(j,k,bi,bj)
           else
            recip_center(j,k,bi,bj) = 0. _d 0
           endif
          enddo
         enddo
        else
         do k=1,Nr
          do j = 1-OLy,sNy+OLy
           fld_rhs(j,k,bi,bj) = 0. _d 0
           fld_in (j,k,bi,bj) = 0. _d 0
           fld_out(j,k,bi,bj) = 0. _d 0
           recip_center(j,k,bi,bj) = 0. _d 0
          enddo
         enddo
        endif
       enddo
      enddo
    
C --- The jacobi iterations
      do while ((jacobi_iter.lt.smooth2DJacobiMaxIters(smoothOpNb))
     &     .and.(relErr.ge.TOL_SQ))

       ! reset error and update counter
       fld_in_norm = 0. _d 0
       relErr = 0. _d 0
       jacobi_iter = jacobi_iter + 1

      ! --- Odd indices
       do bj = jtlo,jthi
        do bi = itlo,ithi
         tile_in_norm(bi,bj) = 0. _d 0
         tile_absErr(bi,bj) = 0. _d 0
         IF (((smooth2DDims(smoothOpNb).EQ.'yze').AND.
     &        (tileHasOBE(bi,bj))).OR.
     &       ((smooth2DDims(smoothOpNb).EQ.'yzw').AND.
     &        (tileHasOBW(bi,bj)))) then
          do j = 1,sNy
           do k = 1,Nr
            fld_out(j,k,bi,bj) = fld_rhs(j,k,bi,bj) -
     &          north(j,k,bi,bj)*fld_in(j+1,k,bi,bj) -
     &          south(j,k,bi,bj)*fld_out(j-1,k,bi,bj)

            if (k.lt.Nr) fld_out(j,k,bi,bj) = fld_out(j,k,bi,bj) -
     &          down(j,k,bi,bj)*fld_in(j,k+1,bi,bj)

            if (k.gt.1) fld_out(j,k,bi,bj) = fld_out(j,k,bi,bj) -
     &          up(j,k,bi,bj)*fld_out(j,k-1,bi,bj)

            fld_out(j,k,bi,bj) = fld_out(j,k,bi,bj) *
     &          sor*recip_center(j,k,bi,bj) +
     &          (1-sor)*fld_in(j,k,bi,bj)

           enddo
          enddo
         endif
        enddo
       enddo

       ! Update halo region
       CALL SMOOTH_EXCH_YZ_RL ( fld_out, smoothOpNb, myThid )

       do bj = jtlo,jthi
        do bi = itlo,ithi
         IF (((smooth2DDims(smoothOpNb).EQ.'yze').AND.
     &        (tileHasOBE(bi,bj))).OR.
     &       ((smooth2DDims(smoothOpNb).EQ.'yzw').AND.
     &        (tileHasOBW(bi,bj)))) then
          do k = 1,Nr
           do j = 1,sNy

            ! Compute error
            tile_in_norm(bi,bj) = tile_in_norm(bi,bj) +
     &                                 fld_in(j,k,bi,bj)**2
            tile_absErr(bi,bj) = tile_absErr(bi,bj) +
     &          (fld_out(j,k,bi,bj)-fld_in(j,k,bi,bj))**2

            ! Reset initial estimate
            fld_in(j,k,bi,bj)=fld_out(j,k,bi,bj)
           enddo
          enddo
         endif
        enddo
       enddo

       CALL GLOBAL_SUM_TILE_RL(tile_in_norm, fld_in_norm, myThid)
       CALL GLOBAL_SUM_TILE_RL(tile_absErr, absErr, myThid)

       relErr = absErr / fld_in_norm

       ! Update halos
       CALL SMOOTH_EXCH_YZ_RL ( fld_in, smoothOpNb, myThid )

      enddo ! end dowhile <TOL and <n_iter_max

      if ( debugLevel.GE.debLevC ) then

        write(msgBuf,'(2A,I10,A,E22.14)') 'SMOOTH_JACOBI_YZ: ',
     &      'final iter: ', jacobi_iter,
     &      '  final error: ', SQRT(relErr)
        call print_message(msgBuf, standardMessageUnit,
     &      SQUEEZE_RIGHT, myThid)

        write(msgBuf,'(2A,I10,A,E22.14)') 'SMOOTH_JACOBI_YZ: ',
     &      'n_iter_max: ', smooth2DJacobiMaxIters(smoothOpNb),
     &      '  TOL: ', TOL
        call print_message(msgBuf, standardMessageUnit,
     &      SQUEEZE_RIGHT, myThid)
      endif

# endif
#endif

      end
