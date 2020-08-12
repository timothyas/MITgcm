#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_AUTODIFF
# include "AUTODIFF_OPTIONS.h"
#endif

      subroutine smooth_get_laplacian_yz (
     U     center, north, south, up, down,
     I     smoothOpNb, smooth2DMask, mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_get_laplacian_yz
C     | o get the matrix to the discretized 2D 
C     |   laplacian operator: div( kappa * grad( ) )
C     | o only the diagonal elements of the diffusion tensor:
C     |
C     |             Kvy Kvz    Kvy  0 
C     |    kappa =  Kwy Kwz ->  0  Kwz
C     |
C     |   are implemented currently
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
C           north  = (j+1, k  ); south  = (j-1, k  )
C           up     = (j  , k+1); down   = (j  , k-1)
      _RL center(1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL north (1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL south (1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL up    (1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL down  (1-OLy:sNy+OLy,Nr,nSx,nSy)

      integer smoothOpNb
      _RL smooth2Dmask(1-OLy:sNy+OLy,Nr,nSx,nSy)
      integer myThid

#ifdef ALLOW_OBCS
# if (defined ALLOW_OBCS_EAST) || (defined ALLOW_OBCS_WEST)

C --- Local parameters
      _RL yA(1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL ZA(1-OLy:sNy+OLy,Nr,nSx,nSy)
      character*( 80) fnamegeneric
      integer i,j,k,bi,bj
      integer itlo,ithi
      integer jtlo,jthi

C --- Loop counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)

C --- Initialize output
      DO bj = jtlo,jthi
       DO bi = itlo,ithi
        DO k = 1Nr
         DO j = 1-OLy,sNy+OLy
          center(j,k,bi,bj)=0. _d 0
          north (j,k,bi,bj)=0. _d 0
          south (j,k,bi,bj)=0. _d 0
          up    (j,k,bi,bj)=0. _d 0
          down  (j,k,bi,bj)=0. _d 0
          yA    (j,k,bi,bj)=0. _d 0
          zA    (j,k,bi,bj)=0. _d 0
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C --- Get discretized Laplacian
C     Note: these are the (swapped sign of) the diagonal
C           elements in applied to input field in smooth_rhs
      DO bj = jtlo,jthi
       DO bi = itlo,ithi

C ---   Eastside
        IF ((smooth2DDims(smoothOpNb).EQ.'yze').AND.
     &      (tileHasOBE(bi,bj)) then
         DO j = 1-OLy,sNy+OLy
          i = OB_Ie(j,bi,bj)
          DO k = 1,Nr

            yA(j,k,bi,bj) = _dxG(i,j,bi,bj) * drF(k) *
     &          _hFacS(i,j,k,bi,bj)  * _recip_dyC(i,j,bi,bj) *
     &          smoothYZ_Kvy(j,k,bi,bj)

            zA(j,k,bi,bj) = _rA(i,j,bi,bj) * 
     &          _hFacC(i,j,k,bi,bj)  * recip_drC(k) *
     &          smoothYZ_Kwz(j,k,bi,bj)

            recip_vol(j,k,bi,bj) = _recip_hFacC(i,j,k,bi,bj)*
     &          recip_rA(i,j,bi,bj) * recip_drF(k)

          ENDDO
         ENDDO

C ---   Westside
        ELSEIF ((smooth2DDims(smoothOpNb).EQ.'yzw').AND.
     &          (tileHasOBW(bi,bj)) then
         DO j = 1-OLy,sNy+OLy
          i = OB_Iw(j,bi,bj)
          DO k = 1,Nr

            yA(j,k,bi,bj) = _dxG(i,j,bi,bj) * drF(k) *
     &          _hFacS(i,j,k,bi,bj)  * _recip_dyC(i,j,bi,bj) *
     &          smoothYZ_Kvy(j,k,bi,bj)

            zA(j,k,bi,bj) = _rA(i,j,bi,bj) * 
     &          _hFacC(i,j,k,bi,bj)  * recip_drC(k) *
     &          smoothYZ_Kwz(j,k,bi,bj)

            recip_vol(j,k,bi,bj) = _recip_hFacC(i,j,k,bi,bj)*
     &          recip_rA(i,j,bi,bj) * recip_drF(k)

          ENDDO
         ENDDO
        ENDIF
         

        IF (((smooth2DDims(smoothOpNb).EQ.'yze').AND.
     &       (tileHasOBE(bi,bj))).OR.
     &      ((smooth2DDims(smoothOpNb).EQ.'yzw').AND.
     &       (tileHasOBW(bi,bj)))) then
         DO k = 1,Nr
          DO j = 1-OLy+1,sNy+OLy-1

           north(j,k,bi,bj)=north(j,k,bi,bj) + yA(j+1,k,bi,bj)*
     &      smooth2Dmask(j,k,bi,bj)*smooth2Dmask(j+1,k,bi,bj)*
     &      recip_vol(j,k,bi,bj)

           south(j,k,bi,bj)=south(j,k,bi,bj) + yA(j,k,bi,bj)*
     &      smooth2Dmask(j-1,k,bi,bj)*smooth2Dmask(j,k,bi,bj)*
     &      recip_vol(j,k,bi,bj)

           if (k.eq.1) then
            up(j,k,bi,bj) = 0. _d 0
           else
            up(j,k,bi,bj)=up(j,k,bi,bj) + zA(j,k,bi,bj)*
     &       smooth2Dmask(j,k-1,bi,bj)*smooth2Dmask(j,k,bi,bj)*
     &       recip_vol(j,k,bi,bj)
           endif

           if (k.eq.Nr) then
            down(j,k,bi,bj) = 0. _d 0
           else
            down(j,k,bi,bj)=down(j,k,bi,bj) + zA(j,k+1,bi,bj)*
     &       smooth2Dmask(j,k,bi,bj)*smooth2Dmask(j,k+1,bi,bj)*
     &       recip_vol(j,k,bi,bj)
           endif

           center(j,k,bi,bj) = center(j,k,bi,bj) -
     &       north(j,k,bi,bj) - south(j,k,bi,bj) -
     &       up   (j,k,bi,bj) - south(j,k,bi,bj)

          ENDDO
         ENDDO
        ENDIF
       ENDDO
      ENDDO

      CALL EXCH_YZ_RL( center, myThid )
      CALL EXCH_YZ_RL( north, myThid )
      CALL EXCH_YZ_RL( south, myThid )
      CALL EXCH_YZ_RL( up, myThid )
      CALL EXCH_YZ_RL( down, myThid )

# endif
#endif

      end
