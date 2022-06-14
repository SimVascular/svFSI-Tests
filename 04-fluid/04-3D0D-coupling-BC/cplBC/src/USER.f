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
      pConv = 1D0
      qConv = 1D0

!     Number of time step between N and N+alpha
      nTimeSteps = 100
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!     List of all coupled faces, BC groups and Xptr is specified here
!     face%bGrp : either cplBC_Dir or cplBC_Neu
!     face%name: the name of the face in svFSI input file
!     face%Xptr: the corresponding index of X for this face
!--------------------------------------------------------------------
      face(1)%bGrp  = cplBC_Dir
      face(2)%bGrp  = cplBC_Neu
      face(1)%name  = "lumen_inlet"
      face(2)%name  = "lumen_outlet"
      face(1)%Xptr  = 1
      face(2)%Xptr  = 2

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

!     RCR parameters
      REAL(KIND=8) :: Rp, C, Rd
      REAL(KIND=8), PARAMETER :: pi = ATAN(1D0)*4D0

!     BC
      Rp = 121D0
      C  = 1.5D-4
      Rd = 1212D0

      f(1) = -40D0 * pi * pi * SIN(2D0 * pi * t) 
      f(2)= (1D0/C) * (Q(2) - x(2)/Rd)
      offst(2) = Q(2)*Rp

!     Assign the additional parameters to be printed
!     It is strongly recommended that Xw(1) = t!
      Xw(1)=t
      Xw(2)=offst(2)

      RETURN
      END SUBROUTINE cplBC_FINDF
!####################################################################
