#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

      subroutine smooth_jacobi_2d (
     I     fld_rhs,
     O     fld_out,
     I     center, north, south, east, west,
     I     smooth2Dmask,
     I     mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_jacobi_2d
C     | o Invert a matrix A via jacobi iterations to solve for u,
C     |   given right hand side b:
C     |          Au = b
C     |
C     | o matrix A is "Laplacian like" and defined by 
C     |   center, north, south, east, west fields
C     |
C     | o b := fld_rhs 
C     |   u:= fld_out
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
#include "SMOOTH.h"

C --- Inputs
C   center, north, south, east, west :
C       fields that define the sparse matrix A at each point
C       for a point (i,j)
C           center = (i,j)
C           north = (i,j+1); south = (i,j-1)
C           east  = (i+1,j); west  = (i-1,j)
      _RL fld_rhs(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL fld_out(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL center(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL north(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL south(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL east(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL west(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL smooth2Dmask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      integer myThid

C --- Local parameters
      integer i,j,k, bi, bj
      integer itlo,ithi
      integer jtlo,jthi

      _RL TOL, TOL_SQ, initial_norm, relErr, absErr
      integer jacobi_iter, n_iter_max

      _RL fld_initial(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL tile_initial_norm(nSx,nSy)
      _RL tile_absErr(nSx,nSy)
      character(max_len_mbuf) msgBuf

C --- Loop counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

C --- Initialize solver stuff
      jacobi_iter = 0
      n_iter_max = 500
      TOL = 1.0 _d -8
      TOL_SQ = TOL*TOL
      relErr = 100. _d 0
      absErr = 0. _d 0 
      
      do bj = jtlo,jthi
       do bi = itlo,ithi
        do j = 1-OLy,sNy+OLy
         do i = 1-OLx,sNx+OLx
          fld_initial(i,j,bi,bj) = 1. _d 0
          fld_out(i,j,bi,bj) = 1. _d 0
         enddo
        enddo
       enddo
      enddo
    
C --- The jacobi iterations
      do while ((jacobi_iter.lt.n_iter_max).and.(relErr.ge.TOL_SQ))

       ! reset error and update counter
       initial_norm = 0. _d 0
       relErr = 0. _d 0
       jacobi_iter = jacobi_iter + 1

       do bj = jtlo,jthi
        do bi = itlo,ithi
         do j = 1,sNy
          do i = 1,sNx
           if ((smooth2Dmask(i,j,bi,bj).NE.0. _d 0).and.
     &         (center(i,j,bi,bj).ne.0. _d 0) ) then

                ! Update "final" estimate
                fld_out(i,j,bi,bj) = fld_rhs(i,j,bi,bj) -
     &              north(i,j,bi,bj)*fld_initial(i,j+1,bi,bj) -
     &              south(i,j,bi,bj)*fld_initial(i,j-1,bi,bj) -
     &              east(i,j,bi,bj) *fld_initial(i+1,j,bi,bj) -
     &              west(i,j,bi,bj) *fld_initial(i-1,j,bi,bj)


                fld_out(i,j,bi,bj) = 
     &              fld_out(i,j,bi,bj)/center(i,j,bi,bj)
           else
            fld_out(i,j,bi,bj)=0. _d 0
           endif
          enddo
         enddo
        enddo
       enddo

       ! Update halo region
       CALL EXCH_XY_RL ( fld_out , myThid )

       do bj = jtlo,jthi
        do bi = itlo,ithi
         do j = 1,sNy
          do i = 1,sNx

           ! Compute error
           tile_initial_norm(bi,bj) = tile_initial_norm(bi,bj) +
     &                                fld_initial(i,j,bi,bj)**2
           tile_absErr(bi,bj) = tile_absErr(bi,bj) +
     &         (fld_out(i,j,bi,bj)-fld_initial(i,j,bi,bj))**2

          enddo
         enddo
        enddo
       enddo

       CALL GLOBAL_SUM_TILE_RL(tile_initial_norm, initial_norm, myThid)
       CALL GLOBAL_SUM_TILE_RL(tile_absErr, absErr, myThid)

       relErr = absErr / initial_norm

       write(msgbuf,'(A,I,A,F)') 'SMOOTH_JACOBI_2D: iter ',
     &  jacobi_iter, ' absErr: ',absErr
       call print_message(msgbuf,standardmessageunit,
     &  SQUEEZE_RIGHT,mythid)

       write(msgbuf,'(A,I,A,F)') 'SMOOTH_JACOBI_2D: iter ',
     &  jacobi_iter, ' initial_norm: ',initial_norm
       call print_message(msgbuf,standardmessageunit,
     &  SQUEEZE_RIGHT,mythid)

       write(msgbuf,'(A,I,A,F)') 'SMOOTH_JACOBI_2D: iter ',
     &  jacobi_iter, ' relErr: ',relErr
       call print_message(msgbuf,standardmessageunit,
     &  SQUEEZE_RIGHT,mythid)


           ! Reset initial estimate
       do bj = jtlo,jthi
        do bi = itlo,ithi
         do j = 1,sNy
          do i = 1,sNx
           fld_initial(i,j,bi,bj)=fld_out(i,j,bi,bj)
          enddo
         enddo
        enddo
       enddo
       CALL EXCH_XY_RL ( fld_initial , myThid )

      enddo ! end dowhile <TOL and <n_iter_max

      end
