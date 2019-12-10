#include "SMOOTH_OPTIONS.h"

      subroutine smooth_apply_matern_op2d (
     U     fld_in,smooth2DMask,smoothOpNb,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_apply_matern_op2d
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
      integer myThid

C --- Local parameters
      _RL fld_in_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL kappa_squared(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL filtervar_norm(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      integer smoothOpNb
      character*( 80) fnamegeneric
      integer i,j,bi,bj
      integer itlo,ithi
      integer jtlo,jthi


C --- Loop counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

C --- Make sure input mask halos are good to go
      CALL EXCH_XY_RL( smooth2DMask, myThid )

C --- Read in K^2:
C     Provide this file at run time
      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth2Dkappa',smoothOpNb
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
     &           1, kappa_squared,1,1,mythid)

c --- 1. Apply k^2 to fld_in and store in separate field
       do bj = jtlo,jthi
        do bi = itlo,ithi
         do j = 1,sNy
          do i = 1,sNx
            fld_in_k(i,j,bi,bj) = fld_in(i,j,bi,bj)*
     &          kappa_squared(i,j,bi,bj)
          enddo
         enddo
        enddo
       enddo
    
c --- 2. Apply laplacian
      call smooth_apply_laplacian2d(fld_in,smooth2DMask,mythid)

c --- 3. Stick 'em together: K^2(fld) - laplacian(fld)
       do bj = jtlo,jthi
        do bi = itlo,ithi
         do j = 1,sNy
          do i = 1,sNx
           fld_in(i,j,bi,bj) = (fld_in_k(i,j,bi,bj) -
     &         fld_in(i,j,bi,bj))*smooth2DMask(i,j,bi,bj)
          enddo
         enddo
        enddo
       enddo
       CALL EXCH_XY_RL ( fld_in, myThid )
       
      end
