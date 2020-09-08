#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

      subroutine smooth_jacobi_3d (
     U     fld_in,
     I     smoothOpNb,
     I     center, east, west, north, south, up, down,
     I     smooth3Dmask,
     I     mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_jacobi_3d
C     | o Invert matrix A via jacobi iterations to solve for u,
C     |   given right hand side b:
C     |          Au = b
C     |
C     | o matrix A is "Laplacian like" and defined by
C     |   center, north, south, east, west fields
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

C --- Inputs
C   center, north, south, east, west, up, down :
C       fields that define the laplacian stencil at each point
C       for a point (i,j,k)
C           center = (i,j,k)
C           east  = (i+1, j  ,  k  ); west  = (i-1, j  , k  )
C           north = (i  , j+1,  k  ); south = (i  , j-1, k  )
C           up    = (i  , j  ,  k+1); down  = (i  , j  , k-1)
C
      _RL fld_in(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      integer smoothOpNb
      _RL center(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL east  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL west  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL up    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL down  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL smooth3Dmask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      integer myThid

C --- Local parameters
      integer i,j,k, bi, bj
      integer itlo,ithi
      integer jtlo,jthi

      _RL TOL, TOL_SQ, fld_in_norm, relErr, absErr
      _RL sor
      integer jacobi_iter

      _RL fld_rhs(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL fld_out(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL recip_center(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
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
      sor = smooth3DSOROmega(smoothOpNb)

      do bj = jtlo,jthi
       do bi = itlo,ithi
        do k = 1,Nr
         do j = 1-OLy,sNy+OLy
          do i = 1-OLx,sNx+OLx
           fld_rhs(i,j,k,bi,bj) = fld_in(i,j,k,bi,bj)
           fld_in (i,j,k,bi,bj) = 1. _d 0
           fld_out(i,j,k,bi,bj) = 1. _d 0
           if (ABS(center(i,j,k,bi,bj)).gt.0. _d 0) then
            recip_center(i,j,k,bi,bj) = 1. _d 0 / center(i,j,k,bi,bj)
           else
            recip_center(i,j,k,bi,bj) = 0. _d 0
           endif
          enddo
         enddo
        enddo
       enddo
      enddo

C --- The jacobi iterations
      do while ((jacobi_iter.lt.smooth3DJacobiMaxIters(smoothOpNb))
     &     .and.(relErr.ge.TOL_SQ))

       ! reset error and update counter
       fld_in_norm = 0. _d 0
       relErr = 0. _d 0
       jacobi_iter = jacobi_iter + 1

       do bj = jtlo,jthi
        do bi = itlo,ithi
         tile_in_norm(bi,bj) = 0. _d 0
         tile_absErr(bi,bj) = 0. _d 0
         do k = 1,Nr
          do j = 1,sNy
           do i = 1,sNx
            ! Update "final" estimate
            fld_out(i,j,k,bi,bj) = fld_rhs(i,j,k,bi,bj) -
     &          east(i,j,k,bi,bj) *fld_in(i+1,j,k,bi,bj) -
     &          west(i,j,k,bi,bj) *fld_out(i-1,j,k,bi,bj) -
     &          north(i,j,k,bi,bj)*fld_in(i,j+1,k,bi,bj) -
     &          south(i,j,k,bi,bj)*fld_out(i,j-1,k,bi,bj)


            if (k.lt.Nr) fld_out(i,j,k,bi,bj) = fld_out(i,j,k,bi,bj) -
     &          down(i,j,k,bi,bj)*fld_in(i,j,k+1,bi,bj)

            if (k.gt.1) fld_out(i,j,k,bi,bj) = fld_out(i,j,k,bi,bj) -
     &          up(i,j,k,bi,bj)*fld_out(i,j,k-1,bi,bj)

            fld_out(i,j,k,bi,bj) = fld_out(i,j,k,bi,bj) *
     &          sor*recip_center(i,j,k,bi,bj) +
     &          (1-sor)*fld_in(i,j,k,bi,bj)
           enddo
          enddo
         enddo
        enddo
       enddo

       ! Update halo region
       CALL EXCH_XYZ_RL ( fld_out , myThid )

       do bj = jtlo,jthi
        do bi = itlo,ithi
         do k = 1,Nr
          do j = 1,sNy
           do i = 1,sNx

            ! Compute error
            tile_in_norm(bi,bj) = tile_in_norm(bi,bj) +
     &                                 fld_in(i,j,k,bi,bj)**2
            tile_absErr(bi,bj) = tile_absErr(bi,bj) +
     &          (fld_out(i,j,k,bi,bj)-fld_in(i,j,k,bi,bj))**2

            ! Reset initial estimate
            fld_in(i,j,k,bi,bj)=fld_out(i,j,k,bi,bj)
           enddo
          enddo
         enddo
        enddo
       enddo

       CALL GLOBAL_SUM_TILE_RL(tile_in_norm, fld_in_norm, myThid)
       CALL GLOBAL_SUM_TILE_RL(tile_absErr, absErr, myThid)

       relErr = absErr / fld_in_norm

       ! Update halos
       CALL EXCH_XYZ_RL ( fld_in , myThid )

      enddo ! end dowhile <TOL and <n_iter_max

      if ( debugLevel.GE.debLevC ) then

        write(msgBuf,'(2A,I10,A,E22.14)') 'SMOOTH_JACOBI_3D: ',
     &      'final iter: ', jacobi_iter,
     &      '  final error: ', SQRT(relErr)
        call print_message(msgBuf, standardMessageUnit,
     &      SQUEEZE_RIGHT, myThid)

        write(msgBuf,'(2A,I10,A,E22.14)') 'SMOOTH_JACOBI_3D: ',
     &      'n_iter_max: ', smooth3DJacobiMaxIters(smoothOpNb),
     &      '  TOL: ', TOL
        call print_message(msgBuf, standardMessageUnit,
     &      SQUEEZE_RIGHT, myThid)
      endif

      end
