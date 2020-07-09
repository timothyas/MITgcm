#include "SMOOTH_OPTIONS.h"

      subroutine smooth_apply_matern_3d (
     U     fld_in,smooth3DMask,smoothOpNb,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_apply_matern_3d
C     | o Routine that applies the following differential operator
C     |   to any field, u
C     |         
C     |         delta u - div( kappa * grad(u) )
C     |
C     |   where delta is a spatially varying field and kappa
C     |   is the tensor described in smooth_get_laplacian_3d
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

C --- Inputs
      _RL fld_in      (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL smooth3DMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      integer myThid

C --- Local parameters
      _RL delta (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL center(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL east  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL west  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL up    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL down  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL fld_out(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      integer smoothOpNb
      character*( 80) fnamegeneric
      integer i,j,k,bi,bj
      integer itlo,ithi
      integer jtlo,jthi


C --- Loop counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

c --- Initialize some containers
c     Not sure what best practices are here...
c     I know these get initialized to standard values
c     when sent to other subroutines, but it
c     seems "unsafe" to not initialize them here
c     since they're created here
      do bj = jtlo,jthi
       do bi = itlo,ithi
        do k = 1,Nr
         do j = 1-OLy,sNy+OLy
          do i = 1-OLx,sNx+OLx
           delta    (i,j,k,bi,bj) = 0. _d 0
           center   (i,j,k,bi,bj) = 0. _d 0
           east     (i,j,k,bi,bj) = 0. _d 0
           west     (i,j,k,bi,bj) = 0. _d 0
           north    (i,j,k,bi,bj) = 0. _d 0
           south    (i,j,k,bi,bj) = 0. _d 0
           up       (i,j,k,bi,bj) = 0. _d 0
           down     (i,j,k,bi,bj) = 0. _d 0
           fld_out  (i,j,k,bi,bj) = 0. _d 0
          enddo
         enddo
        enddo
       enddo
      enddo

C --- Make sure input mask halos are good to go
      CALL EXCH_XYZ_RL( smooth3DMask, myThid )

C --- Read in delta:
C     Provide this file at run time
      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth2Dkappa',smoothOpNb
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
     &           Nr, delta,1,1,mythid)

c --- 1. Apply inverse of filternorm (i.e. filter standard deviation)
      call smooth_apply_filterstd_3d(fld_in, smoothOpNb, myThid)

c --- 2. Get laplacian matrix elements and apply matern operator
      call smooth_get_laplacian_3d(
     &      center, east, west, north, south, up, down,
     &      smoothOpNb, smooth3DMask, myThid)

      do bj = jtlo,jthi
       do bi = itlo,ithi
        do k = 1,Nr
         do j = 1,sNy
          do i = 1,sNx
           fld_out(i,j,k,bi,bj) = fld_in(i,j,k,bi,bj)*
     &          (delta(i,j,k,bi,bj) - center(i,j,k,bi,bj) -
     &          east (i,j,k,bi,bj)*fld_in(i+1,j,k,bi,bj) -
     &          west (i,j,k,bi,bj)*fld_in(i  ,j,k,bi,bj) -
     &          north(i,j,k,bi,bj)*fld_in(i,j+1,k,bi,bj) -
     &          south(i,j,k,bi,bj)*fld_in(i,j  ,k,bi,bj)
           if (k.lt.Nr) fld_out(i,j,k,bi,bj) =
     &          down(i,j,k,bi,bj)*fld_in(i,j,k+1,bi,bj)
           if (k.gt.1) fld_out(i,j,k,bi,bj) =
     &          up  (i,j,k,bi,bj)*fld_in(i,j,k-1,bi,bj)


           ! Is this necessary? I can't keep track anymore...
           fld_out(i,j,k,bi,bj) = fld_out(i,j,k,bi,bj)*
     &          smooth3DMask(i,j,k,bi,bj)
          enddo
         enddo
        enddo
       enddo
      enddo

     CALL EXCH_XYZ_RL ( fld_out, myThid )

c --- 3. Rewrite input <- output
       do bj = jtlo,jthi
        do bi = itlo,ithi
         do k = 1,Nr
          do j = 1-OLy,sNy+OLy
           do i = 1-OLx,sNx+OLx
            fld_in(i,j,k,bi,bj) = fld_out(i,j,k,bi,bj)
           enddo
          enddo
         enddo
        enddo
       enddo

      end
