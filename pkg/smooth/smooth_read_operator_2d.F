#include "SMOOTH_OPTIONS.h"

      subroutine smooth_read_operator_2d (smoothOpNb,mythid)

C     *==========================================================*
C     | SUBROUTINE smooth_read_operator_2d
C     | o Read in the tensor Kappa:
C     |
C     |             Kux Kuy
C     |     Kappa = Kvx Kvy
C     |
C     |   and also the field and delta (matern)
C     |
C     | o for WC01, read single file: smooth2Doperator<smoothOpNb>
C     |   with 2 fields in order:
C     |
C     |     Kux, Kvy
C     |
C     | o for matern algortihm, read three files for three
C     |   entries to the Kappa tensor
C     |     smooth2DKux<smoothOpNb>
C     |     smooth2DKvy<smoothOpNb>
C     |
C     |   as well as the file smooth2DDelta<smoothOpNb>
C     |
C     |   which form the operator: A = (delta - div(Kappa grad())
C     |
C     | o for matern, all entries except Kux, Kvy are zero
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

C --- Local parameters
      character*( 80) fnamegeneric
      character*(MAX_LEN_MBUF) msgBuf
      integer i,j,bi,bj
      integer itlo,ithi
      integer jtlo,jthi
      logical fexist

C --- Loop Counters
      jtlo = mybylo(mythid)
      jthi = mybyhi(mythid)
      itlo = mybxlo(mythid)
      ithi = mybxhi(mythid)


C --- WC01: look for single file, smooth2Doperator
      if (smooth2DAlgorithm(smoothOpNb).eq.'WC01') then

c       read smoothing [i.e diffusion] operator:
        write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth2Doperator',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
     &             1, smooth2D_Kux,1,1,mythid)
        CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
     &             1, smooth2D_Kvy,2,1,mythid)


        else
          WRITE(msgBuf,'(3A)') 'SMOOTH_READ_OPERATOR_2D: ',
     &      'cannot find file ', fnamegeneric
          CALL PRINT_ERROR( msgBuf, myThid )
          STOP 'ABNORMAL END: S/R SMOOTH_READ_OPERATOR_2D' 
        endif ! look for smooth2Doperator file

      elseif (smooth2DAlgorithm(smoothOpNb).eq.'matern') then

c --- Initialize tensor with default values
        do bj = jtlo,jthi
         do bi = itlo,ithi
          do j = 1-OLy,sNy+OLy
           do i = 1-OLx,sNx+OLx
            smooth2D_Kux = 1. _d 0
            smooth2D_Kvy = 1. _d 0
            smooth2DDelta = 0. _d 0
           enddo
          enddo
         enddo
        enddo

C ---   Kux
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &        'smooth2DKux',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        if (fexist) then
          CALL READ_REC_3D_RL(fnamegeneric, smoothprec, 1,
     &          smooth2D_Kux, 1, 1, myThid )
        else
          WRITE(msgBuf,'(3A)') 'SMOOTH_READ_OPERATOR_2D: ',
     &      'cannot find file ', fnamegeneric,
     &      ' setting Kux to ones'
          CALL PRINT_ERROR( msgBuf, myThid )
        endif

C ---   Kvy
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &        'smooth2DKvy',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        if (fexist) then
          CALL READ_REC_3D_RL(fnamegeneric, smoothprec, 1,
     &          smooth2D_Kvy, 1, 1, myThid )
        else
          WRITE(msgBuf,'(3A)') 'SMOOTH_READ_OPERATOR_2D: ',
     &      'cannot find file ', fnamegeneric,
     &      ' setting Kvy to ones'
          CALL PRINT_ERROR( msgBuf, myThid )
        endif

C ---   delta:
        fexist = .FALSE.
        write(fnamegeneric(1:80),'(1a,i3.3,a)')
     &      'smooth2DDelta',smoothOpNb,'.data'
        inquire(file=fnamegeneric,exist=fexist)
        IF (fexist) THEN
          CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
     &               1, smooth2DDelta,1, 1, myThid )
          CALL EXCH_XY_RL ( smooth2DDelta, myThid )
        ELSE
          WRITE(msgBuf,'(4a)')
     &      'SMOOTH_READ_OPERATOR_2D: ',
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
          CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
     &               1, smooth2DRandNorm,1, 1, myThid )
          CALL EXCH_XY_RL ( smooth2DRandNorm, myThid )
        else
          write(msgBuf,'(4a)')
     &      'SMOOTH_READ_OPERATOR_2D: ',
     &      'cannot find ', fnamegeneric,
     &      ' setting to ones'
          call print_message(msgBuf, standardMessageUnit,
     &        SQUEEZE_RIGHT, myThid)
        endif

      ENDIF ! Check algorithm 'WC01' or 'matern'

c ---  read normalization field [i.e. 1/sqrt(var(filter))]:
c     ctrl_check ensures this file exists or is computed, if needed
      fexist = .FALSE.
      write(fnamegeneric(1:80),'(1a,i3.3)')
     &    'smooth2Dnorm',smoothOpNb
      inquire(file=fnamegeneric,exist=fexist)
      IF (fexist) THEN
        CALL READ_REC_3D_RL(fnamegeneric,smoothprec,
     &             1, smooth2Dnorm,1,1,mythid)
        CALL EXCH_XY_RL ( smooth2Dnorm, myThid )
      ENDIF

      CALL EXCH_XY_RL ( smooth2D_Kux, myThid )
      CALL EXCH_XY_RL ( smooth2D_Kvy, myThid )


      RETURN
      END
