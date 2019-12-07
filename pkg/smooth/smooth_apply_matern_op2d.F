#include "SMOOTH_OPTIONS.h"

      subroutine smooth_apply_matern_op2d (
     U     fld_in,mask_in,smoothOpNb,mythid)

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
      _RL mask_in(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nR,nSx,nSy)
      integer myThid

C --- Local parameters
      _RL fld_in_k(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
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


C --- Read in K^2 = Kux:
C     This is a hack.
C     Provide this file at run time, and make both records
C     the same, because for now will just use Kux
C     To be improved upon when I have more intuition for this...
      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth2Doperator',smoothOpNb
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
     &           1, smooth2D_Kux,1,1,mythid)
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
     &           1, smooth2D_Kvy,2,1,mythid)
      CALL EXCH_XY_RL ( smooth2D_Kux, myThid )
      CALL EXCH_XY_RL ( smooth2D_Kvy, myThid )

c --- read normalization field [i.e. 1/sqrt(var(filter))]:
      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth2Dnorm',smoothOpNb
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
     &           1, smooth2Dnorm,1,1,mythid)
      CALL EXCH_XY_RL ( smooth2Dnorm, myThid )

c --- 1. Apply k^2 to fld_in and store in separate field
       do bj = jtlo,jthi
        do bi = itlo,ithi
         do j = 1,sNy
          do i = 1,sNx
            fld_in_k(i,j,bi,bj) = fld_in(i,j,bi,bj)*
     &          smooth2D_Kux(i,j,bi,bj)
          enddo
         enddo
        enddo
       enddo
    
      CALL EXCH_XY_RL ( fld_in_k, myThid )

c --- 2. Apply laplacian
      call smooth_apply_laplacian2d(fld_in,mask_in,mythid)

c --- 3. Stick 'em together: K^2(fld) - laplacian(fld)
       do bj = jtlo,jthi
        do bi = itlo,ithi
         do j = 1,sNy
          do i = 1,sNx
            fld_in(i,j,bi,bj) = fld_in_k(i,j,bi,bj) -
     &          fld_in(i,j,bi,bj)
          enddo
         enddo
        enddo
       enddo

c --- 4. Normalize: divide by ~sqrt(var(filter)):
       do bj = jtlo,jthi
        do bi = itlo,ithi
         DO j = 1,sNy
          DO i = 1,sNx
       fld_in(i,j,bi,bj)=fld_in(i,j,bi,bj)
     & *smooth2Dnorm(i,j,bi,bj)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      CALL EXCH_XY_RL ( fld_in , myThid )
       
      end
