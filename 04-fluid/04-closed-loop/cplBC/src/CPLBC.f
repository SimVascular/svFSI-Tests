!
! Copyright (c) Stanford University, The Regents of the University of
!               California, and others.
!
! All Rights Reserved.
!
! See Copyright-SimVascular.tXt for additional details.
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
!     Here data are received from 3D domain and it setup the data
!     required for integration of ODE's inside cplBC_FINDF
!
!--------------------------------------------------------------------
!
      PROGRAM cplBC_Integ_X
      IMPLICIT NONE

      INCLUDE "cplBC.h"

      INTEGER i, j, k, n, nFaces, nX, nXprnt, nTimeSteps, fid, istat
      REAL(KIND=8) pConv, qConv, dt, dt_3D, t, trk, r
      CHARACTER(LEN=128) cplBC_commu_name
      CHARACTER flag

      REAL(KIND=8), ALLOCATABLE :: offst(:), X(:), Xrk(:), frk(:,:),
     2   Qrk(:,:), Prk(:,:), Xprnt(:)
      TYPE(cplFaceType), ALLOCATABLE :: inFace(:), face(:)

      fid = 1
      istat = 0

      i = IARGC()
      IF (i .EQ. 0) THEN
         STOP "0D-3D communication file name is required as an argument"
      ELSEIF (i .GT. 1) THEN
         STOP "Too many arguments"
      END IF
      CALL GETARG(1,cplBC_commu_name)

!--------------------------------------------------------------------
!     Reading data required by 0D domain. This includes comparison data.
      OPEN(fid, FILE=cplBC_commu_name, STATUS='OLD', FORM='UNFORMATTED')
      READ(fid) nFaces
      READ(fid) nX
      READ(fid) nXprnt
      READ(fid) dt_3D
      READ(fid) t
      ALLOCATE (inFace(nFaces), face(nFaces), offst(nFaces), X(nX),
     2   Xrk(nX), frk(nX,4), Prk(nFaces,4), Qrk(nFaces,4),
     3   Xprnt(nXprnt))
      READ(fid) X
      DO i=1, nFaces
         READ(fid) inFace(i)%bGrp
         READ(fid) inFace(i)%Qo
         READ(fid) inFace(i)%Qn
         READ(fid) inFace(i)%Po
         READ(fid) inFace(i)%Pn
         READ(fid) inFace(i)%name
         inFace(i)%y  = 0D0
      END DO
      CLOSE(fid)

!     Initialize cplBC variables based on USER inputs
      CALL cplBC_INI(nFaces, nTimeSteps, qConv, pConv, face)

!     Set up equivalence between 3D domain and cplBC user inputs
      DO i=1, nFaces
         DO j=1, nFaces
            IF (inFace(i)%name .EQ. face(j)%name) THEN
               inFace(i)%eqv = j
               inFace(i)%Xptr = face(j)%Xptr
               IF (inFace(i)%bGrp .NE. face(j)%bGrp) THEN
                  PRINT *, "CPLBC-ERROR: face ", i,
     2               " bGrp is different in 0D and 3D domains"
                  istat = -1
                  CALL WRITECOMM()
                  RETURN
               END IF
               EXIT
            END IF
         END DO
         IF (j .GT. nFaces) THEN
            PRINT *, "CPLBC-ERROR: face "//TRIM(inFace(i)%name)//
     2         " not found in cplBC"
            istat = -1
            CALL WRITECOMM()
            RETURN
         END IF
      END DO

!     Scaling flow rates and pressures
      DO i=1, nFaces
         j = inFace(i)%eqv
         face(j)%eqv = i
         face(j)%Qo = inFace(i)%Qo/qConv
         face(j)%Qn = inFace(i)%Qn/qConv
         face(j)%Po = inFace(i)%Po/pConv
         face(j)%Pn = inFace(i)%Pn/pConv
      END DO

!--------------------------------------------------------------------
!     Setting up the system of equations
      offst = 0D0
      dt    = dt_3D/REAL(nTimeSteps,KIND=8)
      DO n=1, nTimeSteps
         DO i=1, 4
            r = REAL(i-1,KIND=8)/3D0
            r = (REAL(n-1,KIND=8) + r)/REAL(nTimeSteps,KIND=8)
            Qrk(:,i) = face(:)%Qo + (face(:)%Qn - face(:)%Qo)*r
            Prk(:,i) = face(:)%Po + (face(:)%Pn - face(:)%Po)*r
         END DO

!        RK-4 1st pass
         trk = t
         Xrk = X
         CALL cplBC_FINDF(trk, nFaces, nX, Xrk, frk(:,1), Qrk(:,1),
     2      Prk(:,1), offst, nXprnt, Xprnt)

!        RK-4 2nd pass
         trk = t + dt/3.0D0
         Xrk = X + dt*frk(:,1)/3.0D0
         CALL cplBC_FINDF(trk, nFaces, nX, Xrk, frk(:,2), Qrk(:,2),
     2      Prk(:,2), offst, nXprnt, Xprnt)

!        RK-4 3rd pass
         trk = t + 2.0D0*dt/3.0D0
         Xrk = X - dt*frk(:,1)/3.0D0 + dt*frk(:,2)
         CALL cplBC_FINDF(trk, nFaces, nX, Xrk, frk(:,3), Qrk(:,3),
     2      Prk(:,3), offst, nXprnt, Xprnt)

!        RK-4 4th pass
         trk = t + dt
         Xrk = X + dt*frk(:,1) - dt*frk(:,2) + dt*frk(:,3)
         CALL cplBC_FINDF(trk, nFaces, nX, Xrk, frk(:,4), Qrk(:,4),
     2      Prk(:,4), offst, nXprnt, Xprnt)

         r = dt/8.0D0
         X = X + r*(frk(:,1) + 3.0D0*(frk(:,2) + frk(:,3)) + frk(:,4))

         DO j=1, nX
            IF (ISNAN(X(j))) THEN
               PRINT *, "CPLBC-ERROR: NaN detected in cplBC"
               istat = -1
               CALL WRITECOMM()
               RETURN
            END IF
         END DO
         t = t + dt
      END DO

      DO i=1, nFaces
         j = inFace(i)%Xptr
         k = inFace(i)%eqv
         IF (inFace(i)%bGrp .EQ. cplBC_Dir) THEN
            inFace(i)%y = qConv*X(j)
         ELSE
            inFace(i)%y = pConv*(X(j) + offst(k))
         END IF
      END DO

      CALL WRITECOMM()

      IF (flag .EQ. 'L') THEN
         OPEN(fid,FILE='cplBC_AllData',STATUS='UNKNOWN',ACCESS='APPEND')
         DO i=1, nX
            WRITE (fid,"(ES14.6E2)",ADVANCE='NO') X(i)
         END DO
         DO i=1, nXprnt
            WRITE (fid,"(ES14.6E2)",ADVANCE='NO') Xprnt(i)
         END DO
         WRITE (fid,*)
         CLOSE(fid)
      END IF

      DEALLOCATE(face, offst, X, Xrk, frk, Prk, Qrk, Xprnt)

      CONTAINS
!-----------------------------------------------------------------------
      SUBROUTINE WRITECOMM()
      IMPLICIT NONE

!     Writing the data required by the 3D domain
      OPEN(fid,FILE=cplBC_commu_name,STATUS='OLD',FORM='UNFORMATTED')
      WRITE(fid) istat
      WRITE(fid) X
      WRITE(fid) Xprnt
      DO i=1, nFaces
         WRITE(fid) inFace(i)%y
      END DO
      CLOSE(fid)

      RETURN
      END SUBROUTINE WRITECOMM
!-----------------------------------------------------------------------
      END PROGRAM cplBC_Integ_X
!######################################################################
