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
!     Additional parameters used in LPN network are defined here
!
!--------------------------------------------------------------------

!     Specifing the constants
      REAL*8, PARAMETER :: ! Order of surfaces: 3~8
     &  Rp(20) = (/
     &  1.431310,
     &  2.007010,
     &  2.555090,
     &  1.245410,
     &  1.016040,
     &  3.002900,
     &  3.482080,
     &  0.322510,
     &  0.847280,
     &  0.617790,
     &  1.035160,
     &  1.437190,
     &  0.297140,
     &  1.723380,
     &  0.135830,
     &  0.201830,
     &  0.554560,
     &  0.834200,
     &  1.193090,
     &  1.308140
     &  /),
     &  C(20) = (/ 
     &  0.009780,
     &  0.006550,
     &  0.005980,
     &  0.016460,
     &  0.017580,
     &  0.003210,
     &  0.002640,
     &  0.117470,
     &  0.030750,
     &  0.052250,
     &  0.015540,
     &  0.013560,
     &  0.211300,
     &  0.009670,
     &  0.579620,
     &  0.298450,
     &  0.080490,
     &  0.031290,
     &  0.023560,
     &  0.013530
     &  /),
     &  Rd(20) = (/ 
     &  12.470410,
     &  16.206790,
     &  18.744990,
     &  7.8476100,
     &  6.7754100,
     &  35.766800,
     &  42.362630,
     &  1.3878500,
     &  4.2446900,
     &  3.3050900,
     &  8.7510400,
     &  9.9326500,
     &  0.8572400,
     &  11.998300,
     &  0.3652900,
     &  0.6721400,
     &  2.0215100,
     &  4.1241700,
     &  5.8012100,
     &  9.6971100
     &  /),
     &  CCsa= 4.00000000D-01,
     &  csa= 2.00000000D-01,
     &  dsa= 3.80000000D-01,
     &  Vsa0= 1.46,
     &  const(3)=(/ 2.20000000D-01,
     &   5.00000000D-02,
     &   4.30000000D-01/),
     &  a= -1.28700000D-01,
     &  b= 8.50000000D+00,
     &  Csv= 1.42200000D+00,
     &  dsv= 7.00000000D-02,
     &  Vsv0= 3.40000000D+00,
     &  Tsvs_temp= 3.10000000D-01,
     &  Rmyo= 4.51000000D-02,
     &  Lav= 2.74000000D-05,
     &  Kav= 8.16000000D-04,
     &  Kao= 8.59000000D-05,
     &  Ruba= 2.01000000D+00,
     &  Luba= 6.56000000D-04,
     &  Cub= 5.70000000D-01,
     &  Rubv= 5.87000000D+00,
     &  Cao= 1.50000000D-01,
     &  Rthao= 2.88000000D-01,
     &  Lthao= 2.28000000D-03,
     &  Cthao= 2.10000000D-02,
     &  Rabao= 2.19000000D+00,
     &  Labao= 2.28000000D-03,
     &  Cabao= 4.45000000D-02,
     &  Rlega= 6.73000000D+00,
     &  Llega= 2.28000000D-03,
     &  Clega= 1.85000000D-02,
     &  Rlegc= 1.58900000D+01,
     &  Clegv= 3.76000000D-01,
     &  Rlegv= 3.71000000D+00,
     &  R0abivc= 5.39000000D-02,
     &  V0abivc= 3.72000000D+00,
     &  R0thivc= 8.38000000D-02,
     &  V0thivc= 4.88000000D+00,
     &  Rla= 3.03000000D+01,
     &  Cl= 6.52000000D-01,
     &  Rlv= 1.90000000D-01,
     &  Rka= 2.13400000D+01,
     &  Ck= 2.72000000D-01,
     &  Rkv= 1.99900000D+00,
     &  Ria= 4.67900000D+01,
     &  Ci= 1.60000000D-01,
     &  Riv= 8.51000000D-01,
     &  P0ith= 0.00000000D+00,
     &  P0iab= 0.00000000D+00,
     &  Apith= 0.00000000D+00,
     &  Apiab= 0.00000000D+00,
     &  fc= 1.20000000D+02,
     &  fr= 3.00000000D+01,
     &  coll= 2.00000000D+01,
     &  pi= 3.14159265D+00
!     The parameters which might change through the code 
      REAL*8 Cthivc, Cabivc, Csvc
      Csvc= 2.56000000D-01
      Cabivc= 2.29000000D-01
      Cthivc= 6.87000000D-02
