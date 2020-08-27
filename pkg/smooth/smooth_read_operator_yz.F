#include "SMOOTH_OPTIONS.h"
#ifdef ALLOW_OBCS
# include "OBCS_OPTIONS.h"
#endif

      subroutine smooth_read_operator_yz (smoothOpNb,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_read_operator_yz
C     | o Read in the tensor Kappa:
C     |
C     |             Kvy Kvz
C     |     Kappa = Kvz Kwz
C     |
C     |   and also the field and delta (matern)
C     |
C     | o for matern algortihm, read three files for three
C     |   entries to the Kappa tensor
C     |     smoothYZKvy<smoothOpNb>
C     |     smoothYZKwz<smoothOpNb>
C     |
C     |   as well as the file smoothYZDelta<smoothOpNb>
C     |
C     |   which form the operator: A = (delta - div(Kappa grad())
C     |
C     | o for matern, all entries except Kvy, Kwz are zero
C     |   and if those two files are not found, they are set
C     |   to 1's. If delta is not found, it is set to zero.
C     |
C     *==========================================================*

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "GRID.h"
#include "PARAMS.h"
#include "SMOOTH.h"

C --- Inputs
      integer smoothOpNb
      integer myThid

#ifdef ALLOW_OBCS
# if (defined ALLOW_OBCS_EAST) || (defined ALLOW_OBCS_WEST)

C --- Local parameters
      character*( 80) fnamegeneric
      character*(MAX_LEN_MBUF) msgBuf
      integer i,j,k,bi,bj
      integer itlo,ithi
      integer jtlo,jthi
      logical fexist

C --- Loop Counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)


      if (smooth2DAlgorithm(smoothOpNb)(1:6).eq.'matern') then

c --- Initialize tensor with default values
        do bj = jtlo,jthi
         do bi = itlo,ithi
          do k = 1,Nr
           do j = 1-OLy,sNy+OLy
            smoothYZ_Kvy(j,k,bi,bj)  = 1. _d 0
            smoothYZ_Kwz(j,k,bi,bj)  = 1. _d 0
            smoothYZDelta(j,k,bi,bj) = 0. _d 0
            smoothYZRandNorm(j,k,bi,bj) = 1. _d 0
            smoothYZNorm(j,k,bi,bj) = 0. _d 0
           enddo
          enddo
         enddo
        enddo

C ---   Kvy
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &        'smooth2DKvy',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        if (fexist) then
          CALL READ_REC_YZ_RL(fnamegeneric, smoothprec, Nr,
     &          smoothYZ_Kvy, 1, 1, myThid )
          CALL SMOOTH_EXCH_YZ_RL( smoothYZ_Kvy, smoothOpNb, myThid )
        else
          WRITE(msgBuf,'(3A)') 'SMOOTH_READ_OPERATOR_YZ: ',
     &      'cannot find file ', fnamegeneric,
     &      ' setting Kvy to ones'
          CALL PRINT_ERROR( msgBuf, myThid )
        endif

C ---   Kwz
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &        'smooth2DKwz',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        if (fexist) then
          CALL READ_REC_YZ_RL(fnamegeneric, smoothprec, Nr,
     &          smoothYZ_Kwz, 1, 1, myThid )
          CALL SMOOTH_EXCH_YZ_RL( smoothYZ_Kwz, smoothOpNb, myThid )
        else
          WRITE(msgBuf,'(3A)') 'SMOOTH_READ_OPERATOR_YZ: ',
     &      'cannot find file ', fnamegeneric,
     &      ' setting Kwz to ones'
          CALL PRINT_ERROR( msgBuf, myThid )
        endif

C ---   delta:
        fexist = .FALSE.
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &      'smooth2DDelta',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        IF (fexist) THEN
          CALL READ_REC_YZ_RL(fnamegeneric,smoothprec,
     &               Nr, smoothYZDelta,1, 1, myThid )
          CALL SMOOTH_EXCH_YZ_RL( smoothYZDelta, smoothOpNb, myThid )
        ELSE
          WRITE(msgBuf,'(4a)')
     &      'SMOOTH_READ_OPERATOR_YZ: ',
     &      'cannot find ', fnamegeneric,
     &      ' setting to zeros'
          CALL PRINT_ERROR(msgBuf, myThid)
        ENDIF

C ---   Look for normalization of white noise process
C       i.e. 1/sqrt(det(F))
        fexist = .FALSE.
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &      'smooth2DRandNorm',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        if (fexist) then
          CALL READ_REC_YZ_RL(fnamegeneric,smoothprec,
     &               Nr, smoothYZRandNorm,1, 1, myThid )
          CALL SMOOTH_EXCH_YZ_RL( smoothYZRandNorm, smoothOpNb, myThid )
        else
          write(msgBuf,'(4a)')
     &      'SMOOTH_READ_OPERATOR_YZ: ',
     &      'cannot find ', fnamegeneric,
     &      ' setting to ones'
          call print_message(msgBuf, standardMessageUnit,
     &        SQUEEZE_RIGHT, myThid)
        endif

      ENDIF ! Check algorithm 'matern'

c ---  read normalization field [i.e. 1/sqrt(var(filter))]:
c     ctrl_check ensures this file exists or is computed, if needed
      fexist = .FALSE.
      write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &    'smooth2Dnorm',smoothOpNb,'.data'
      inquire(file=fnamegeneric,exist=fexist)
      IF (fexist) THEN
        CALL READ_REC_YZ_RL(fnamegeneric,smoothprec,
     &             Nr, smoothYZnorm,1,1,mythid)
        CALL SMOOTH_EXCH_YZ_RL( smoothYZNorm, smoothOpNb, myThid )
      ENDIF


# endif
#endif

      RETURN
      END
