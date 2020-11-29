#include "SMOOTH_OPTIONS.h"

      subroutine smooth_apply_matern_3d (
     I      fld_in,
     U      fld_reg,
     I      smooth3DMask,smoothOpNb,doNormalize,mythid)

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
C     | o Because run time inputs are processed in ctrl_map_ini_genarr
C     |   this takes fld_in, copies to fld_reg, applies regularization
C     |   and returns that. ctrl_map_ini_genarr then writes that,
C     |   and ctrl_cost_gen2d reads it and adds it to the cost function.
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
      _RL fld_reg     (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL smooth3DMask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      integer smoothOpNb
      logical doNormalize
      integer myThid

C --- Local parameters
      _RL center(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL east  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL west  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL up    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL down  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL fld_out(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      character*( 80) fnamegeneric
      character*(MAX_LEN_MBUF) msgBuf
      integer i,j,k,bi,bj
      integer itlo,ithi
      integer jtlo,jthi
      logical fexist


C --- Loop counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

c --- Initialize
      CALL EXCH_XYZ_RL( fld_in, myThid )
      do bj = jtlo,jthi
       do bi = itlo,ithi
        do k = 1,Nr
         do j = 1-OLy,sNy+OLy
          do i = 1-OLx,sNx+OLx
           fld_out(i,j,k,bi,bj) = 0. _d 0
           fld_reg(i,j,k,bi,bj) = fld_in(i,j,k,bi,bj)
          enddo
         enddo
        enddo
       enddo
      enddo

C --- Make sure input mask halos are good to go
      CALL EXCH_XYZ_RL( smooth3DMask, myThid )

c --- Get laplacian matrix elements
      call smooth_get_laplacian_3d(
     &      center, east, west, north, south, up, down,
     &      smoothOpNb, smooth3DMask, myThid)

c --- 1. Apply inverse of filter normalization ( i.e. sqrt(filtervar))
      if (doNormalize) then
       do bj = jtlo,jthi
        do bi = itlo,ithi
         do k = 1,Nr
          do j = 1-OLy,sNy+OLy
           do i = 1-OLx,sNx+OLx
            if ((ABS(smooth3DNorm(i,j,k,bi,bj)).gt.0. _d 0).and.
     &          (smooth3DMask(i,j,k,bi,bj).ne.0)) then
             fld_reg(i,j,k,bi,bj)=
     &        fld_reg(i,j,k,bi,bj)/smooth3DNorm(i,j,k,bi,bj)
            else
             fld_reg(i,j,k,bi,bj)=0. _d 0
            endif
           enddo
          enddo
         enddo
        enddo
       enddo
      endif

C --- 2. Apply laplacian-like operator
      do bj = jtlo,jthi
       do bi = itlo,ithi
        do k = 1,Nr
         do j = 1,sNy
          do i = 1,sNx

           fld_out(i,j,k,bi,bj) = fld_reg(i,j,k,bi,bj)*
     &          (smooth3DDelta(i,j,k,bi,bj) -
     &          center(i,j,k,bi,bj)) -
     &          east (i,j,k,bi,bj)*fld_reg(i+1,j,k,bi,bj) -
     &          west (i,j,k,bi,bj)*fld_reg(i-1,j,k,bi,bj) -
     &          north(i,j,k,bi,bj)*fld_reg(i,j+1,k,bi,bj) -
     &          south(i,j,k,bi,bj)*fld_reg(i,j-1,k,bi,bj)

           if (k.lt.Nr) fld_out(i,j,k,bi,bj) = fld_out(i,j,k,bi,bj) -
     &          down(i,j,k,bi,bj)*fld_reg(i,j,k+1,bi,bj)
           if (k.gt.1) fld_out(i,j,k,bi,bj) = fld_out(i,j,k,bi,bj) -
     &          up  (i,j,k,bi,bj)*fld_reg(i,j,k-1,bi,bj)


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
            fld_reg(i,j,k,bi,bj) = fld_out(i,j,k,bi,bj)
            if (doNormalize) then
             if (ABS(smooth3DRandNorm(i,j,k,bi,bj)).gt.0.) then
              fld_reg(i,j,k,bi,bj) =
     &        fld_reg(i,j,k,bi,bj)/smooth3DRandNorm(i,j,k,bi,bj)
             else
              fld_reg(i,j,k,bi,bj) = 0. _d 0
             endif
            endif
           enddo
          enddo
         enddo
        enddo
       enddo

      end
