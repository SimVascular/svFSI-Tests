c     MUSC 2 Immediate postop
c     Ethan Kung  keo@ucsd.edu

c     Created by Mahdi Esmaily Moghadam 12-01-2010
c     Please report any problem to mesmaily@ucsd.edu, memt63@yahoo.com

c This subroutine initializes the parameters, you need to read the
c comments and specify them carefuly
c--------------------------------------------------------------------
c This is an example for RCR boundary condition with parameters
c Rd, R, and C which are distal, and proximal resistance and capacitor.

      SUBROUTINE INITIALIZE(nTimeStep)
      USE COM
      IMPLICIT NONE
      INTENT(OUT) nTimeStep

      LOGICAl ierr
      INTEGER i, nTimeStep
      REAL(KIND=8), ALLOCATABLE :: tZeroX(:)
c
c********************************************************************
c For instance if pressure in 3D solver is in cgs and here mmHg
c pConv=1334=1.334D3, also if flowrate in 3D solver is mm^3/s and
c here is mL/s qConv=1000=1D3. In the case both solver are using same
c unites you can set these two conversion coefficients to 1D0
      pConv = 1.334D2
      qConv = 1D3

c Only when all the surfaces of you model are coupled with NeumannSrfs
c you may set this to .TRUE.
      pCorr = .FALSE.
      qCorr = .FALSE.

c********************************************************************
c Block for your inputs

c These two value should match to that of defined in solver.inp
      nDirichletSrfs = 0
      nNeumannSrfs   = 21
c Number of unknowns that you need inside your lumped parameter network
      nUnknowns      = 39
c Number of time step between N and N+alpha, for higher values this code
c would be more accurate and more costly
      nTimeStep = 100

c Number of parameters to be printed in AllData file (the first
c nUnknowns columns are the "X" values)
      nXprint = 31

c--------------------------------------------------------------------
c You don't need to change this part                                 !
      ALLOCATE (tZeroX(nUnknowns), srfToXdPtr(nDirichletSrfs))       !
      ALLOCATE (srfToXPtr(nNeumannSrfs))                             !
      tZeroX = 0D0
c--------------------------------------------------------------------

      INCLUDE "initial_values_final.f"

c--------------------------------------------------------------------
c You don't need to change this part                                 !
      INQUIRE (FILE='InitialData', EXIST=ierr)                       !
      IF (.NOT.ierr) THEN                                            !
c         PRINT *, 'Initializing unknowns in LPM'                     !
         OPEN (1, FILE='InitialData',STATUS='NEW',FORM='UNFORMATTED')!
         WRITE (1) 0D0                                               !
         DO i=1, nUnknowns                                           !
            WRITE (1) tZeroX(i)                                      !
         END DO                                                      !
         CLOSE(1)                                                    !
      END IF                                                         !
c--------------------------------------------------------------------

c Surface to X pointer: this defines which Unknown corresponds to which
c suface in "List of Neumann Surfaces" inside solver.inp
c For example, if you have "List of Neumann Surfaces= 2 8 4 ...."
c and you set "srfToXPtr = (/5,3,9,.../)"
C this means X(5) corresponds to surface 2, X(3) <-> surface 8,
c and X(9) <-> surface 4
c Also, Q(1) corresponds to the first element in srfToXPtr, in this
c example, Q(1)<=>X(5), Q(2)<=>X(3)
      srfToXPtr = (/1,2,3,4,5,6,26,27,28,29,30,31,32,
     &             33,34,35,36,37,38,39,7/)

      END SUBROUTINE INITIALIZE

c####################################################################
c Here you should find the f_i=dx_i/dt, based on the following parameters:
c  current x_i:                   x(i)
c  Current time:                  t
c  Flowrates from 3D code:        Q(i)
c  Pressure from Dirichlet faces: P(i)

      SUBROUTINE FINDF(t, x, f, Q, P)
      USE COM
      IMPLICIT NONE
      INTENT(IN) t, Q
      INTENT(OUT) f

      REAL(KIND=8) t, x(nUnknowns), f(nUnknowns), Q(nNeumannSrfs),
     2   P(nDirichletSrfs)

!     These are the dumy variables
      REAL(KIND=8) Qao, Psa, Psv, Tc, Tm, Tsvs, Tsas, Tsad, T1, Tr,
     2   Tmr, Tsr, Aa, Av, Ar, dAr, Pith, Piab, dPith, dPiab,
     3   Rabivc, Rthivc

      INCLUDE "parameters_final.f"

!********************************************************************
!     Time periods
      Tc   = 60.0/fc
      Tm   = MOD(t,Tc)
      Tsvs = Tsvs_temp
      Tsas = const(1)   !*Tc
      Tsad = Tc - Tsas
      T1   = const(2)   !*Tc

!     Atrium
      IF (Tm .LE. T1) THEN
         Aa = 5D-1*(1D0 - COS(2D0*pi*(Tm - T1 + Tsas)/Tsas))
      ELSEIF (T1+Tsad.LE.Tm .AND. Tm.LE.Tc) THEN
         Aa = 5D-1*(1D0 - COS(2D0*pi*(Tm - T1 - Tsad)/Tsas))
      ELSE
         Aa = 0D0
      END IF
      Psa = Aa*(x(20) - Vsa0)/CCsa + csa*(EXP(dsa*(x(20) - Vsa0)) - 1D0)

!     Ventricle
      IF (Tm .LE. Tsvs) THEN
         Av = 5D-1*(1D0 - COS(2D0*pi*Tm/Tsvs))
      ELSE
         Av = 0D0
      END IF
      Psv = Av*(a*(x(22) - Vsv0)**2D0 + b*(x(22) - Vsv0)) +
     2 Csv*(EXP(dsv*(x(22) - Vsv0)) - 1D0)

!     This is the only case that x is directly manipulated, ugly, but the
!     only way!!!
      IF (x(21) .LT. 0D0) x(21) = 0D0

!     Aortic flow
      IF (Psv .GT. x(10)) THEN
         Qao = (SQRT(Rmyo*Rmyo + 4D0*Kao*(Psv - x(10))) - Rmyo)/2D0/Kao
      ELSE
         Qao = 0D0
      END IF

!     Respiration
      Tr = 6D1/fr
      Tmr = MOD(t,Tr)
      Tsr = const(3)*Tr
      IF (Tmr .LT. Tsr) THEN
         Ar  = 5D-1*(1D0 - COS(2D0*pi*Tmr/Tsr))
         dAr = (pi/Tsr)*SIN(2D0*pi*Tmr/Tsr)
      ELSE
         Ar  = 0D0
         dAr = 0D0
      END IF
      Pith  = APith*Ar + P0ith
      dPith = APith*dAr
      Piab  = APiab*Ar + P0iab
      dPiab = APiab*dAr

!     Adjusting parameters based on collapsibility
      IF (x(7) .LT. Pith) THEN
         Csvc = coll*Csvc
      END IF
      IF (x(18) .LT. Piab) THEN
         Cabivc = coll*Cabivc
         Rabivc = R0abivc*(1D0 + Cabivc*(x(18)-Piab)/V0abivc)**(-3D0)
      ELSE
         Rabivc = R0abivc*(1D0 + Cabivc*(x(18)-Piab)/V0abivc)**(-2D0)
      END IF
      IF (x(19) .LT. Pith) THEN
         Cthivc = coll*Cthivc
         Rthivc = R0thivc*(1D0 + Cthivc*(x(19)-Pith)/V0thivc)**(-3D0)
      ELSE
         Rthivc = R0thivc*(1D0 + Cthivc*(x(19)-Pith)/V0thivc)**(-2D0)
      END IF

!     Added pressure to x(1:6)
      offset(1:6) = Q(1:6)*Rp(1:6)
      offset(26:39) = Q(7:20)*Rp(7:20)

!     The main body of equations
      f(1:6)= (Q(1:6) - (x(1:6) - Psa)/Rd(1:6))/C(1:6)
      f(7)  = (Q(21) + (x(8) - x(7))/Rubv)/Csvc + dPith
      f(8)  = (x(9) - (x(8) - x(7))/Rubv)/Cub
      f(9)  = (x(10) - x(9)*Ruba - x(8))/Luba
      f(10) = (Qao - x(9) - x(11))/Cao + dPith
      f(11) = (x(10) - x(11)*Rthao - x(12))/Lthao
      f(12) = (x(11) - x(13) - (x(12) - x(23))/Rla
     2      - (x(12) - x(24))/Rka)/Cthao + dPith
      f(13) = (x(12) - x(13)*Rabao - x(14))/Labao
      f(14) = (x(13) - (x(14) - x(25))/Ria - x(15))/Cabao + dPiab
      f(15) = (x(14) - x(15)*Rlega - x(16))/Llega
      f(16) = (x(15) - (x(16) - x(17))/Rlegc)/Clega
      f(17) = ((x(16) - x(17))/Rlegc - (abs(x(17) - x(18))
     2      + x(17) - x(18))/Rlegv/2D0)/Clegv
      f(18) = ((abs(x(17) - x(18)) + x(17) - x(18))/Rlegv/2D0
     2      - (x(18) - x(19))/Rabivc)/Cabivc + dPiab
      f(19) = ((x(18) - x(19))/Rabivc + (x(24) - x(19))/Rkv
     2      + (x(23) - x(19))/Rlv - (x(19) - Psa)/Rthivc)/Cthivc + dPith
      f(20) = (x(19) - Psa)/Rthivc + SUM((x(1:6) - Psa)/Rd(1:6))
     &      + SUM((x(26:39) - Psa)/Rd(7:20)) - x(21)
      f(21) = (Psa - Psv - Kav*x(21)*x(21))/Lav
      f(22) = x(21) - Qao
      f(23) = ((x(25) - x(23))/Riv + (x(12) - x(23))/Rla
     2      - (x(23) - x(19))/Rlv)/Cl + dPiab
      f(24) = ((x(12) - x(24))/Rka - (x(24) - x(19))/Rkv)/Ck + dPiab
      f(25) = ((x(14) - x(25))/Ria - (x(25) - x(23))/Riv)/Ci + dPiab
      f(26:39)= (Q(7:20) - (x(26:39) - Psa)/Rd(7:20))/C(7:20)

      IF (x(21).EQ.0D0 .AND. Psa.LT.(Psv-Rmyo*Qao)) THEN
         f(21) = 0D0
      END IF

c     Assign the additional parameters to be printed
      Xprint(1)=t !40
      Xprint(2)=Aa !41
      Xprint(3)=Av !42
      Xprint(4)=Psa !43
      Xprint(5)=Psv !44
      Xprint(6)=Aa*(x(20) - Vsa0)/CCsa !45
      Xprint(7)=csa*(EXP(dsa*(x(20) - Vsa0)) - 1D0) !46
      Xprint(8)=Av*(a*(x(22) - Vsv0)**2D0 + b*(x(22) - Vsv0)) !47
      Xprint(9)=Csv*(EXP(dsv*(x(22) - Vsv0)) - 1D0) !48
      Xprint(10)=Qao !49
      Xprint(11)=(x(19)-Psa)/Rthivc !50
      Xprint(12)=Q(1) !51
      Xprint(13)=Q(2) !52
      Xprint(14)=Q(3) !53
      Xprint(15)=Q(4) !54
      Xprint(16)=Q(5) !55
      Xprint(17)=Q(6) !56
      Xprint(18)=Q(7) !57
      Xprint(19)=Q(8) !58
      Xprint(20)=Q(9) !59
      Xprint(21)=Q(10) !60
      Xprint(22)=Q(11) !61
      Xprint(23)=Q(12) !62
      Xprint(24)=Q(13) !63
      Xprint(25)=Q(14) !64
      Xprint(26)=Q(15) !65
      Xprint(27)=Q(16) !66
      Xprint(28)=Q(17) !67
      Xprint(29)=Q(18) !68
      Xprint(30)=Q(19) !69
      Xprint(31)=Q(20) !70
c      Xprint(32)=-Q(21) !71

      RETURN
      END SUBROUTINE FINDF


