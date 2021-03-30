!
! Copyright (c) Stanford University, The Regents of the University of
!               California, and others.
!
! All Rights Reserved.
!
! See Copyright-SimVascular.txt for additional details.
!
! Permission is hereby granted, free of charge, to any person obtaining
! a copy of this software and associated documentation files (the
! "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish,
! distribute, sublicense, and/or sell copies of the Software, and to
! permit persons to whom the Software is furnished to do so, subject
! to the following conditions:
!
! The above copyright notice and this permission notice shall be included
! in all copies or substantial portions of the Software.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!--------------------------------------------------------------------
!
!     This subroutine initializes the parameters, you need to read the
!     comments and specify them carefuly
!
!--------------------------------------------------------------------
!     This is an example of Windkessel model (RCR) coupled to a Neumann
!     BC

      SUBROUTINE cplBC_INI(nFaces, nTimeSteps, qConv, pConv, face)
      IMPLICIT NONE
      INCLUDE "cplBC.h"

      INTEGER, INTENT(IN) :: nFaces
      INTEGER, INTENT(OUT) :: nTimeSteps
      REAL(KIND=8), INTENT(OUT) :: pConv, qConv
      TYPE(cplFaceType), INTENT(OUT) :: face(nFaces)

!--------------------------------------------------------------------
!     Block for your inputs
!     For instance if pressure in 3D solver is in cgs and here mmHg
!     pConv=1334=1.334D3, also if flowrate in 3D solver is mm^3/s and
!     here is mL/s qConv=1000=1D3. In the case both solver are using
!     same unites you can set these two conversion coefficients to 1D0
      pConv = 1.334D2
      qConv = 1D3

!     Number of time step between N and N+alpha
      nTimeSteps = 100
!--------------------------------------------------------------------

      INCLUDE "faces.f"

      END SUBROUTINE cplBC_INI

!####################################################################
!     Here you should find the f_i=dx_i/dt, based on the following
!     parameters:
!     current x_i:         x(i)
!     Flowrates of faces:  Q(i)
!     Pressures of faces:  P(i)
      SUBROUTINE cplBC_FINDF(t, nFaces, nX, X, f, Q, P, offst, nW, Xw)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nFaces, nX, nW
      REAL(KIND=8), INTENT(IN) :: t, Q(nFaces), P(nFaces)
      REAL(KIND=8), INTENT(OUT) :: f(nX), offst(nFaces)
      REAL(KIND=8), INTENT(INOUT) :: x(nX), Xw(nW)

      REAL(KIND=8) Qao, Psa, Psv, Tc, Tm, Tsvs, Tsas, Tsad, T1, Tmr, Tr,
     2 Tsr, Aa, Av, Ar, dAr, Pith, Piab, dPith, dPiab, Rabivc, Rthivc

      INCLUDE "parameters.f"

!     Other model parameters
      Tc   = 60.0/fc
      Tm   = MOD(t,Tc)
      Tsvs = Tsvs_temp
      Tsas = const(1)
      Tsad = Tc - Tsas
      T1   = const(2)

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

      IF (x(21) .LT. 0D0) x(21) = 0D0

!     Aortic flow
      IF (Psv .GT. x(10)) THEN
         Qao = (SQRT(Rmyo*Rmyo + 4D0*Kao*(Psv - x(10))) - Rmyo)/2D0/Kao
      ELSE
         Qao = 0D0
      END IF

!     Respiration
      Tr = 6.0D1/fr
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

!     Parameters based on collapsibility
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

!     Offset due to proximal pulmonary resistances
      offst(1:20) = Q(1:20)*Rp(1:20)

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
     2      + SUM((x(26:39) - Psa)/Rd(7:20)) - x(21)
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

!     Assign the additional parameters to be printed
      Xw(1)=t !40
      Xw(2)=Aa !41
      Xw(3)=Av !42
      Xw(4)=Psa !43
      Xw(5)=Psv !44
      Xw(6)=Aa*(x(20) - Vsa0)/CCsa !45
      Xw(7)=csa*(EXP(dsa*(x(20) - Vsa0)) - 1D0) !46
      Xw(8)=Av*(a*(x(22) - Vsv0)**2D0 + b*(x(22) - Vsv0)) !47
      Xw(9)=Csv*(EXP(dsv*(x(22) - Vsv0)) - 1D0) !48
      Xw(10)=Qao !49
      Xw(11)=(x(19)-Psa)/Rthivc !50
      Xw(12)=Q(1) !51
      Xw(13)=Q(2) !52
      Xw(14)=Q(3) !53
      Xw(15)=Q(4) !54
      Xw(16)=Q(5) !55
      Xw(17)=Q(6) !56
      Xw(18)=Q(7) !57
      Xw(19)=Q(8) !58
      Xw(20)=Q(9) !59
      Xw(21)=Q(10) !60
      Xw(22)=Q(11) !61
      Xw(23)=Q(12) !62
      Xw(24)=Q(13) !63
      Xw(25)=Q(14) !64
      Xw(26)=Q(15) !65
      Xw(27)=Q(16) !66
      Xw(28)=Q(17) !67
      Xw(29)=Q(18) !68
      Xw(30)=Q(19) !69
      Xw(31)=Q(20) !70

      RETURN
      END SUBROUTINE cplBC_FINDF
!####################################################################
