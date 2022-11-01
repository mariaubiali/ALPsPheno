      logical FUNCTION dummy_cuts(P)
C**************************************************************************
C     INPUT:
C            P(0:3,1)           MOMENTUM OF INCOMING PARTON
C            P(0:3,2)           MOMENTUM OF INCOMING PARTON
C            P(0:3,3)           MOMENTUM OF ...
C            ALL MOMENTA ARE IN THE REST FRAME!!
C            COMMON/JETCUTS/   CUTS ON JETS
C     OUTPUT:
C            TRUE IF EVENTS PASSES ALL CUTS LISTED
C**************************************************************************
      IMPLICIT NONE
c
c     Constants
c
      include 'genps.inc'
      include 'nexternal.inc'
C
C     ARGUMENTS
C
      REAL*8 P(0:3,nexternal)
C
C     PARAMETERS
C
      real*8 PI
      parameter( PI = 3.14159265358979323846d0 )
c
c     particle identification
c
      LOGICAL  IS_A_J(NEXTERNAL),IS_A_L(NEXTERNAL)
      LOGICAL  IS_A_B(NEXTERNAL),IS_A_A(NEXTERNAL),IS_A_ONIUM(NEXTERNAL)
      LOGICAL  IS_A_NU(NEXTERNAL),IS_HEAVY(NEXTERNAL)
      logical  do_cuts(nexternal)
      COMMON /TO_SPECISA/IS_A_J,IS_A_A,IS_A_L,IS_A_B,IS_A_NU,IS_HEAVY,
     . IS_A_ONIUM, do_cuts

c     transverse mass squared, one for each W
      double precision delPhimin, delPhi
      double precision delPhiNum, delPhiDenom
      double precision MTtot1, MTtot2, MTtotW1, MTtotW2, MTmin1, MTmin2
      double precision ppl1(0:3),ppv1(0:3),ppl2(0:3),ppv2(0:3)
      double precision ppmiss(0:3),ppboost(0:3)
      integer i,j,iproc
      include 'maxamps.inc'
      integer idup(nexternal,maxproc,maxsproc)

      dummy_cuts=.true.
        
      DO i=0,3
        ppl1(i)=p(i,3)
        ppv1(i)=p(i,4)
        ppl2(i)=p(i,6)
        ppv2(i)=p(i,7)
        ppmiss(i)=ppv1(i)+ppv2(i)
        ppboost(i)=ppmiss(i)+ppl1(i)+ppl2(i)
      ENDDO

c     delta Phi boost cut < 1.5
      delPhimin = 1.5d0
      delPhiNum = ppmiss(1)*ppboost(1)+ppmiss(2)*ppboost(2)
      delPhiDenom = dsqrt((ppmiss(1)**2 + ppmiss(2)**2)*(ppboost(1)**2 +
     $ ppboost(2)**2))
      delPhi = acos(delPhiNum/delPhiDenom)

c     MT cut, each > 80, to increase tail events
      MTmin1 = 95.0d0
      MTmin2 = 95.0d0
      MTtot1 = 2*dsqrt((ppl1(1)**2 +
     $     ppl1(2)**2)*(ppv1(1)**2 +
     $     ppv1(2)**2)) - 2*(ppl1(1)*ppv1(1) +
     $     ppl1(2)*ppv1(2))

      MTtot2 = 2*dsqrt((ppl2(1)**2 +
     $     ppl2(2)**2)*(ppv2(1)**2 +
     $     ppv2(2)**2)) - 2*(ppl2(1)*ppv2(1) +
     $     ppl2(2)*ppv2(2))


      if (MTtot1.ge.0d0) then
         MTtotW1=dsqrt(MTtot1)
      else
         MTtotW1=0d0
      endif

      if (MTtot2.ge.0d0) then
         MTtotW2=dsqrt(MTtot2)
      else
         MTtotW2=0d0
      endif

      if(delPhi.gt.delPhimin) then
         dummy_cuts=.false.
         return
      endif

      if((MTtotW1.gt.MTmin1).or.(MTtotW2.gt.MTmin2)) then
         dummy_cuts=.false.
         return
      endif

      
      return
      end

      subroutine get_dummy_x1(sjac, X1, R, pbeam1, pbeam2, stot, shat)
      implicit none
      include 'maxparticles.inc'
      include 'run.inc'
c      include 'genps.inc'
      double precision sjac ! jacobian. should be updated not reinit
      double precision X1   ! bjorken X. output
      double precision R    ! random value after grid transfrormation. between 0 and 1
      double precision pbeam1(0:3) ! momentum of the first beam (input and/or output)
      double precision pbeam2(0:3) ! momentum of the second beam (input and/or output)
      double precision stot        ! total energy  (input and /or output)
      double precision shat        ! output

c     global variable to set (or not)
      double precision cm_rap
      logical set_cm_rap
      common/to_cm_rap/set_cm_rap,cm_rap
      
      set_cm_rap=.false. ! then cm_rap will be set as .5d0*dlog(xbk(1)*ebeam(1)/(xbk(2)*ebeam(2)))
                         ! ebeam(1) and ebeam(2) are defined here thanks to 'run.inc'
      shat = x1*ebeam(1)*ebeam(2)
      return 
      end

      subroutine get_dummy_x1_x2(sjac, X, R, pbeam1, pbeam2, stot,shat)
      implicit none
      include 'maxparticles.inc'
      include 'run.inc'
c      include 'genps.inc'
      double precision sjac ! jacobian. should be updated not reinit
      double precision X(2)   ! bjorken X. output
      double precision R(2)    ! random value after grid transfrormation. between 0 and 1
      double precision pbeam1(0:3) ! momentum of the first beam
      double precision pbeam2(0:3) ! momentum of the second beam
      double precision stot        ! total energy
      double precision shat        ! output

c     global variable to set (or not)
      double precision cm_rap
      logical set_cm_rap
      common/to_cm_rap/set_cm_rap,cm_rap
      
      set_cm_rap=.false. ! then cm_rap will be set as .5d0*dlog(xbk(1)*ebeam(1)/(xbk(2)*ebeam(2)))
                         ! ebeam(1) and ebeam(2) are defined here thanks to 'run.inc'
      shat = x(1)*x(2)*ebeam(1)*ebeam(2)
      return 
      end


      logical  function dummy_boostframe()
      implicit none
c
c      
      dummy_boostframe = .false.
      return
      end
      




