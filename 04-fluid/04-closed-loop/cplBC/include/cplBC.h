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
!     This is a header file for the data structure used for
!     communication with svFSI
!
!--------------------------------------------------------------------

!     cplBC type of coupling to between 3D and OD-LPN models:
!     Dirichlet type coupling, Neumann type coupling
      INTEGER, PARAMETER :: cplBC_Dir = 66112, cplBC_Neu = 66113

      TYPE cplFaceType
         SEQUENCE
!        GenBC_Dir/GenBC_Neu
         INTEGER :: bGrp
!        Pointer to x
         INTEGER :: Xptr
!        Internal cplBC use
         INTEGER :: eqv = 0
!        Flow rates at t
         REAL(KIND=8) Qo
!        Flow rates at t+dt
         REAL(KIND=8) Qn
!        Pressures at t
         REAL(KIND=8) Po
!        Pressures at t+dt
         REAL(KIND=8) Pn
!        Imposed flow/pressure
         REAL(KIND=8) y
!        Name of the face
         CHARACTER(LEN=128) name
      END TYPE cplFaceType
