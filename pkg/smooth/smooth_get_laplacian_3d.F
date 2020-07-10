#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

      subroutine smooth_get_laplacian_3d (
     U     center, east, west, north, south, up, down,
     I     smoothOpNb, smooth3DMask, mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_get_laplacian_3d
C     | o get the matrix to the discretized 3D 
C     |   laplacian operator: div( kappa * grad( ) )
C     | o only the diagonal elements of the diffusion tensor:
C     |
C     |             Kux Kuy Kuz    Kux  0   0
C     |    kappa =  Kvx Kvy Kvz ->  0  Kvy  0
C     |             Kwx Kwy Kwz     0   0  Kwz
C     |
C     |   are implemented currently
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
      _RL center(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL east  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL west  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL up    (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL down  (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

      integer smoothOpNb
      _RL smooth3Dmask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      integer myThid

C --- Local parameters
      _RL xA(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL yA(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL zA(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL recip_vol(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      character*( 80) fnamegeneric
      integer i,j,k,bi,bj
      integer itlo,ithi
      integer jtlo,jthi

C --- Loop counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

C --- Diffusion operator (tensor) as shown above
      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth3Doperator',smoothOpNb
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec, 
     &           Nr,smooth3D_Kwx,1, 1, myThid )
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec, 
     &           Nr,smooth3D_Kwy,2, 1, myThid )
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec, 
     &           Nr,smooth3D_Kwz,3, 1, myThid )
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec, 
     &           Nr,smooth3D_Kux,4, 1, myThid )
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec, 
     &           Nr,smooth3D_Kvy,5, 1, myThid )
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec, 
     &           Nr,smooth3D_Kuz,6, 1, myThid )
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec, 
     &           Nr,smooth3D_Kvz,7, 1, myThid )
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec, 
     &           Nr,smooth3D_Kuy,8, 1, myThid )
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec, 
     &           Nr,smooth3D_Kvx,9, 1, myThid )
      CALL READ_REC_3D_RL(fnamegeneric,smoothprec, 
     &           Nr,smooth3D_kappaR,10, 1, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kwx, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kwy, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kwz, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kux, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kvy, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kuz, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kvz, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kuy, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_Kvx, myThid )
      CALL EXCH_XYZ_RL ( smooth3D_kappaR, myThid )
    
C --- Initialize output
      DO bj = jtlo,jthi
       DO bi = itlo,ithi
        DO k = 1,Nr
         DO j = 1-OLy,sNy+OLy
          DO i = 1-OLx,sNx+OLx
           center(i,j,k,bi,bj)=0. _d 0
           east  (i,j,k,bi,bj)=0. _d 0
           west  (i,j,k,bi,bj)=0. _d 0
           north (i,j,k,bi,bj)=0. _d 0
           south (i,j,k,bi,bj)=0. _d 0
           up    (i,j,k,bi,bj)=0. _d 0
           down  (i,j,k,bi,bj)=0. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C --- Get discretized Laplacian
C     Note: these are the (swapped sign of) the diagonal
C           elements in applied to input field in smooth_rhs
      DO bj = jtlo,jthi
       DO bi = itlo,ithi
        DO k = 1,Nr
         DO j = 1-OLy,sNy+OLy
          DO i = 1-OLx,sNx+OLx

            xA(i,j,k,bi,bj) = _dyG(i,j,bi,bj) * drF(k) *
     &          _hFacW(i,j,k,bi,bj)  * _recip_dxC(i,j,bi,bj) *
     &          smooth3D_Kux(i,j,k,bi,bj)

            yA(i,j,k,bi,bj) = _dxG(i,j,bi,bj) * drF(k) *
     &          _hFacS(i,j,k,bi,bj)  * _recip_dyC(i,j,bi,bj) *
     &          smooth3D_Kvy(i,j,k,bi,bj)

            zA(i,j,k,bi,bj) = _rA(i,j,bi,bj) * 
     &          _hFacC(i,j,k,bi,bj)  * recip_drC(k) *
     &          smooth3D_Kwz(i,j,k,bi,bj)

            recip_vol(i,j,k,bi,bj) = _recip_hFacC(i,j,k,bi,bj)*
     &          recip_rA(i,j,bi,bj) * recip_drF(k)

          ENDDO
         ENDDO
        ENDDO
           

        DO k = 1,Nr
         DO j = 1-OLy+1,sNy+OLy-1
          DO i = 1-OLx+1,sNx+OLx-1

           east(i,j,k,bi,bj)=east(i,j,k,bi,bj) + xA(i+1,j,k,bi,bj)*
     &      smooth3Dmask(i,j,k,bi,bj)*smooth3Dmask(i+1,j,k,bi,bj)*
     &      recip_vol(i,j,k,bi,bj)

           west(i,j,k,bi,bj)=west(i,j,k,bi,bj) + xA(i,j,k,bi,bj)*
     &      smooth3Dmask(i-1,j,k,bi,bj)*smooth3Dmask(i,j,k,bi,bj)*
     &      recip_vol(i,j,k,bi,bj)

           north(i,j,k,bi,bj)=north(i,j,k,bi,bj) + yA(i,j+1,k,bi,bj)*
     &      smooth3Dmask(i,j,k,bi,bj)*smooth3Dmask(i,j+1,k,bi,bj)*
     &      recip_vol(i,j,k,bi,bj)

           south(i,j,k,bi,bj)=south(i,j,k,bi,bj) + yA(i,j,k,bi,bj)*
     &      smooth3Dmask(i,j-1,k,bi,bj)*smooth3Dmask(i,j,k,bi,bj)*
     &      recip_vol(i,j,k,bi,bj)

           if (k.eq.1) then
            up(i,j,k,bi,bj) = 0. _d 0
           else
            up(i,j,k,bi,bj)=up(i,j,k,bi,bj) + zA(i,j,k,bi,bj)*
     &       smooth3Dmask(i,j,k-1,bi,bj)*smooth3Dmask(i,j,k,bi,bj)*
     &       recip_vol(i,j,k,bi,bj)
           endif

           if (k.eq.Nr) then
            down(i,j,k,bi,bj) = 0. _d 0
           else
            down(i,j,k,bi,bj)=down(i,j,k,bi,bj) + zA(i,j,k+1,bi,bj)*
     &       smooth3Dmask(i,j,k,bi,bj)*smooth3Dmask(i,j,k+1,bi,bj)*
     &       recip_vol(i,j,k,bi,bj)
           endif

           center(i,j,k,bi,bj) = center(i,j,k,bi,bj) -
     &       east(i,j,k,bi,bj) - west(i,j,k,bi,bj) -
     &       north(i,j,k,bi,bj) - south(i,j,k,bi,bj) -
     &       up(i,j,k,bi,bj) - down(i,j,k,bi,bj)

          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CALL EXCH_XYZ_RL( center, myThid )
      CALL EXCH_XYZ_RL( east, myThid )
      CALL EXCH_XYZ_RL( west, myThid )
      CALL EXCH_XYZ_RL( north, myThid )
      CALL EXCH_XYZ_RL( south, myThid )
      CALL EXCH_XYZ_RL( up, myThid )
      CALL EXCH_XYZ_RL( down, myThid )

      end
