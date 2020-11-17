#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_CTRL
# include "CTRL_OPTIONS.h"
#endif
#ifdef ALLOW_SHELFICE
# include "SHELFICE_OPTIONS.h"
#endif

      subroutine smooth_filtervar3D (smoothOpNb,myThid)

C     *==========================================================*
C     | SUBROUTINE smooth_filtervar3D
C     | o Routine that computes the filter variance
C     |   field associated with a diffusion operator, as part
C     |   a 3D spatial correlation operator (smooth_correld3D.F)
C     |   See Weaver and Courtier 01 for details.
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "EESUPPORT.h"
#include "PARAMS.h"
#include "GRID.h"
#include "SMOOTH.h"
#if (defined ALLOW_SHI2D_CONTROL) || (defined ALLOW_SHIFWFLX_CONTROL)
# include "SHELFICE.h"
#endif

C --- Inputs
      integer smoothOpNb, myThid

C --- Local Parameters
      Real*8   port_rand, port_rand_norm
      EXTERNAL PORT_RAND, PORT_RAND_NORM

      Real*8 rand_seed
      integer timeVals(8)

      integer i,j,k, bi, bj, ii, jj, kk
      integer itlo,ithi, jtlo,jthi
      integer diLoc, djLoc,  dkLoc
      integer nbt_in
      character*( 80) fnamegeneric
      character*(MAX_LEN_MBUF) msgBuf
      logical doNormalize,fexist

      _RL smoothTmpFld (1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL smoothTmpMean(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RL smoothTmpVar(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
      _RS smooth3Dmask(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)

C --- Initialize loop counters
      jtlo = mybylo(myThid)
      jthi = mybyhi(myThid)
      itlo = mybxlo(myThid)
      ithi = mybxhi(myThid)

      doNormalize=.FALSE.

c if smooth3Dfilter(smoothOpNb)=0: the filter variance field
c has been computed earlier and is already in the run directory
c so this routine does not do anything

      IF (smooth3Dfilter(smoothOpNb).NE.0) then

      nbt_in=smooth3Dnbt(smoothOpNb)/2

c initialize filter variance field,
c initialize determinant inverse (for matern), and get mask:
      DO bj=jtlo,jthi
       DO bi=itlo,ithi
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           if (smooth3DmaskName(smoothOpNb)(1:5).eq.'maskC') then
            smooth3Dmask(i,j,k,bi,bj) = maskC(i,j,k,bi,bj)
           elseif (smooth3DmaskName(smoothOpNb)(1:5).eq.'maskW') then
            smooth3Dmask(i,j,k,bi,bj) = maskW(i,j,k,bi,bj)
           elseif (smooth3DmaskName(smoothOpNb)(1:5).eq.'maskS') then
            smooth3Dmask(i,j,k,bi,bj) = maskS(i,j,k,bi,bj)
           endif
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C --- Read smoothing [i.e diffusion] operator:
      CALL smooth_read_operator_3d(smoothOpNb,myThid)

      IF (smooth3Dfilter(smoothOpNb).EQ.2) then
c compute the normalization matrix using the approximate method
c
c This method can be quite expensive -- so that the approximate
c method (see below) is usually the prefered one.
c The exact method can be used to check the accuracy
c of the approximate method results (that can be predicted).
c
c note: the exact method requires the adjoint of smooth_diff2D.F (see below)

C --- Check which algorithm is used
       if (smooth3DAlgorithm(smoothOpNb)(1:6).eq.'matern') then
        WRITE(errorMessageUnit,'(A,/,A)' )
     &   "Exact method for matern prior covariance not "
     &   "implemented yet... but would be pretty expensive."
        CALL ALL_PROC_DIE( myThid )
        STOP 'ABNORMAL END: S/R smooth_filtervar3D'

       else


        diLoc=15 !int(5*smooth_L/smooth_dx)
        djLoc=20 !int(5*smooth_L/smooth_dx)
        dkLoc=8

        DO kk=1,dkLoc
        DO ii=1,diLoc,2
        DO jj=1,djLoc,2

        DO bj=jtlo,jthi
         DO bi=itlo,ithi
          DO k=1,Nr
           DO j=1-OLy,sNy+OLy
            DO i=1-OLx,sNx+OLx
             smoothTmpFld(i,j,k,bi,bj)=0.
            ENDDO
           ENDDO
          ENDDO

          DO k=kk,Nr,dkLoc
           DO j=jj,sNy,djLoc
            DO i=ii,sNx,diLoc
             smoothTmpFld(i,j,k,bi,bj)=1.
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO

c note: as we go to adjoint part, we need to have 0 in overlaps
c         so we must NOT have done an exchange for smoothTmpFld

c adjoint:
        WRITE(errorMessageUnit,'(A,/,A)' )
     &   "you need to have adsmooth_diff3D compiled and then:",
     &   "uncomment the line below and comment the stop"
        CALL ALL_PROC_DIE( myThid )
        STOP 'ABNORMAL END: S/R smooth_filtervar3D'
c       call adsmooth_diff3D(smoothTmpFld,nbt_in,myThid)

c division by sqrt(volume)*sqrt(volume) [1 to end adj, 1 to begin fwd]
        DO bj=jtlo,jthi
         DO bi=itlo,ithi
          DO k=1,Nr
           DO j=1,sNy
            DO i=1,sNx
c division by ~sqrt(volume):
          smoothTmpFld(i,j,k,bi,bj)=smoothTmpFld(i,j,k,bi,bj)
     &   *(recip_rA(i,j,bi,bj)*recip_drF(k))
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO

c coming out of adjoint part: overlaps are 0
c going in fwd part: we need to fill them up
        CALL EXCH_XYZ_RL ( smoothTmpFld , myThid )

c fwd:
        call smooth_diff3D(smoothTmpFld,nbt_in,myThid)

c convert variance to normalization factor:
        DO bj=jtlo,jthi
         DO bi=itlo,ithi
          DO k=1,Nr,dkLoc
           DO j=jj,sNy,djLoc
            DO i=ii,sNx,diLoc
            if (smoothTmpFld(i,j,k,bi,bj).NE.0.) then
             smooth3Dnorm(i,j,k,bi,bj)=
     &   1/sqrt(smoothTmpFld(i,j,k,bi,bj))
            endif
            ENDDO
           ENDDO
          ENDDO
         ENDDO
        ENDDO

        ENDDO      !DO ii=1,diLoc
        ENDDO      !DO jj=1,djLoc
        ENDDO      !DO kk=1,dkLoc

       ENDIF ! check smooth algorithm

      ELSEIF (smooth3Dfilter(smoothOpNb).EQ.1) then
c compute the normalization matrix using the approximate method

C --- 1. Initialize mean and variance
      DO bj=jtlo,jthi
       DO bi=itlo,ithi
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           smoothTmpMean(i,j,k,bi,bj)   = 0. _d 0
           smoothTmpVar(i,j,k,bi,bj)   = 0. _d 0
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

c --- 2. Initialize random number generator
c     seed is 1023 * pid * milliseconds from date/time
      rand_seed = 1023. _d 0
      call date_and_time(VALUES=timeVals)
      rand_seed = rand_seed * REAL((1+myPid)*timeVals(8),8)
      smoothTmpFld(1,1,1,1,1)=port_rand(rand_seed)

      DO ii=1,smooth3DNbRand(smoothOpNb)


C --- 3. Initialize field within random sample loop
C     set  points where mask!=0 to random number
C     WC01: divide by sqrt(area)
C     matern: multiply by 1/sqrt(det(F)), provided by user
      DO bj=jtlo,jthi
       DO bi=itlo,ithi
        DO k=1,Nr
         DO j=1,sNy
          DO i=1,sNx

           smoothTmpFld(i,j,k,bi,bj)   = 0. _d 0
           if (smooth3dmask(i,j,k,bi,bj).NE.0) then

             smoothTmpFld(i,j,k,bi,bj)=port_rand_norm(-rand_seed)

             if (smooth3DAlgorithm(smoothOpNb)(1:4).eq.'WC01') then
              smoothTmpFld(i,j,k,bi,bj)=smoothTmpFld(i,j,k,bi,bj)
     &        *sqrt(recip_rA(i,j,bi,bj)*recip_drF(k))

             elseif (smooth3DAlgorithm(smoothOpNb)(1:6).eq.'matern') then
              smoothTmpFld(i,j,k,bi,bj)=smoothTmpFld(i,j,k,bi,bj)*
     &         smooth3DRandNorm(i,j,k,bi,bj)
             endif
           endif
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      CALL EXCH_XYZ_RL ( smoothTmpFld, myThid )

c --- 4. Smooth random field with appropriate algorithm
      if (smooth3DAlgorithm(smoothOpNb)(1:4).eq.'WC01')
     &  call smooth_diff3D(smoothTmpFld,nbt_in,myThid)
      if (smooth3DAlgorithm(smoothOpNb)(1:6).eq.'matern')
     &  call smooth_inverse_matern_3d(
     &      smoothTmpFld,smooth3DMask,smoothOpNb,doNormalize,myThid)

c --- 5. accumulate statistics (to compute the variance later)
      DO bj=jtlo,jthi
       DO bi=itlo,ithi
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
      smoothTmpVar(i,j,k,bi,bj)=smoothTmpVar(i,j,k,bi,bj)
     & +smoothTmpFld(i,j,k,bi,bj)*smoothTmpFld(i,j,k,bi,bj)/
     &  smooth3DNbRand(smoothOpNb)
      smoothTmpMean(i,j,k,bi,bj)=smoothTmpMean(i,j,k,bi,bj)
     & +smoothTmpFld(i,j,k,bi,bj)/smooth3DNbRand(smoothOpNb)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

C --- temporary, write out all samples
      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth3Dfld',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, smoothTmpFld, ii, 1, myThid )

       WRITE(standardMessageUnit,'(A,I4,A,I4)')
     & 'smooth_filtervar3D: ',ii,' members done out of ',
     &  smooth3DNbRand(smoothOpNb)

      ENDDO

c --- 6. compute variance and convert it to normalization factor:
      DO bj=jtlo,jthi
       DO bi=itlo,ithi
        DO k=1,Nr
         DO j=1-OLy,sNy+OLy
          DO i=1-OLx,sNx+OLx
           if (smooth3dmask(i,j,k,bi,bj).NE.0) then
           smooth3Dnorm(i,j,k,bi,bj)=
     & 1/sqrt ( smooth3DNbRand(smoothOpNb)/
     & (smooth3DNbRand(smoothOpNb)-1)* ( smoothTmpVar(i,j,k,bi,bj) -
     & smoothTmpMean(i,j,k,bi,bj)*smoothTmpMean(i,j,k,bi,bj)
     &  )  )
           endif
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      ENDIF

c write smooth3Dnorm_3D to file:
      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth3Dnorm',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, smooth3Dnorm, 1, 1, myThid )
      CALL EXCH_XYZ_RL ( smooth3Dnorm,  myThid )

c write out mean and the last sample
      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth3Dmean',smoothOpNb
      CALL WRITE_REC_3D_RL( fnamegeneric, smoothprec,
     &                      Nr, smoothTmpMean, 1, 1, myThid )

      ENDIF

      RETURN
      END
