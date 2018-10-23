!***********************************************************************

      module varmod
      use stdParams
      use genUtils

      integer, parameter :: nsd = 3

      integer, parameter :: eType_NA = 100, eType_LIN = 101, &
         eType_TRI = 102, eType_TET = 103, eType_BIL = 104, &
         eType_QUD = 105, eType_BIQ = 106, eType_BRK = 107, &
         eType_NRB = 108, eType_WDG = 109

      type faceType
         integer :: nNo
         integer :: nEl
         integer :: eNoN
         integer :: nG
         integer :: eType = eType_NA
         integer :: vtkType
         integer, allocatable :: IEN(:,:)
         integer, allocatable :: gN(:)
         integer, allocatable :: gE(:)
         real(kind=8), allocatable :: w(:)
         real(kind=8), allocatable :: xi(:,:)
         real(kind=8), allocatable :: N(:,:)
         real(kind=8), allocatable :: x(:,:)
         real(kind=8), allocatable :: Nx(:,:,:)
         character(len=strL) :: name
      end type faceType

      type mshType
         integer :: nNo
         integer :: nEl
         integer :: eNoN
         integer :: nG
         integer :: eType = eType_NA
         integer :: vtkType
         integer :: nFa
         integer, allocatable :: IEN(:,:)
         real(kind=8), allocatable :: w(:)
         real(kind=8), allocatable :: xi(:,:)
         real(kind=8), allocatable :: N(:,:)
         real(kind=8), allocatable :: x(:,:)
         real(kind=8), allocatable :: Nx(:,:,:)
         character(len=strL) :: name
         type(faceType), allocatable :: fa(:)
      end type mshType

      interface destroy
         module procedure destroyMsh, destroyFa
      end interface destroy

      integer :: nx, ny, nz
      real(kind=8) :: dx, dy, dz
      real(kind=8), allocatable :: x(:), y(:), z(:)

      type(mshType) :: msh

      contains
!-------------------------------------------------
         subroutine destroyFa(lFa)
         implicit none
         type(faceType), intent(inout) :: lFa

         if (allocated(lFa%IEN)) deallocate(lFa%IEN)
         if (allocated(lFa%gN))  deallocate(lFa%gN)
         if (allocated(lFa%gE))  deallocate(lFa%gE)
         if (allocated(lFa%w))   deallocate(lFa%w)
         if (allocated(lFa%xi))  deallocate(lFa%xi)
         if (allocated(lFa%N))   deallocate(lFa%N)
         if (allocated(lFa%x))   deallocate(lFa%x)
         if (allocated(lFa%Nx))  deallocate(lFa%Nx)

         return
         end subroutine destroyFa
!-------------------------------------------------
         subroutine destroyMsh(lM)
         implicit none
         type(mshType), intent(inout) :: lM
         integer iFa

         if (allocated(lM%IEN))  deallocate(lM%IEN)
         if (allocated(lM%w))    deallocate(lM%w)
         if (allocated(lM%xi))   deallocate(lM%xi)
         if (allocated(lM%N))    deallocate(lM%N)
         if (allocated(lM%x))    deallocate(lM%x)
         if (allocated(lM%Nx))   deallocate(lM%Nx)

         if (allocated(lM%fa)) then
            do iFa=1, lM%nFa
               call destroyFa(lM%fa(iFa))
            end do
            deallocate(lM%fa)
         end if

         return
         end subroutine destroyMsh
!-------------------------------------------------
         pure function cross(U) result(V)
         implicit none
         real(kind=8), intent(in) :: U(:,:)
         real(kind=8) V(size(U,1))

         if (size(U,1) .eq. 2) then
            V(1) =  U(2,1)
            V(2) = -U(1,1)
         else
            V(1) = U(2,1)*U(3,2) - U(3,1)*U(2,2)
            V(2) = U(3,1)*U(1,2) - U(1,1)*U(3,2)
            V(3) = U(1,1)*U(2,2) - U(2,1)*U(1,2)
         end if

         return
         end function cross
!-------------------------------------------------
      end module varmod

!***********************************************************************

      program tet_mesher
      use varmod
      implicit none

      character(len=strL) :: c
      integer :: a, b, e, g, i, j, k, l, Ac, iFa, nNo, nEl, nid(8), &
     &   iord(4,6)
      real(kind=8) :: vol, Jac, Ks(nsd,nsd)

      real(kind=8), allocatable :: xl(:,:), N(:), Nxi(:,:)

      msh%name = "cube"

      i = IARGC()
      if (i .eq. 0) then
         write(stdout,ftab4) "ERROR: num edge elements not specified"
         STOP
      else if (i .gt. 1) then
         write(stdout,ftab4) "ERROR: Too many arguments"
         STOP
      end if

      call getarg(1, c)
      read(c,*,iostat=i) nx
      nx = nx + 1
      dx = 1.0D0/real(nx-1,kind=8)

!     Define cube grid parameters
      ny = nx
      nz = nx
      dy = dx
      dz = dx
      allocate(x(nx), y(ny), z(nz))
      x = 0.0d0
      y = 0.0d0
      z = 0.0d0

!     Cube mesh coordinates
      do i=1, nx
         x(i) = real(i-1,kind=8)*dx
      end do

      do j=1, ny
         y(j) = real(j-1,kind=8)*dy
      end do

      do k=1, nz
         z(k) = real(k-1,kind=8)*dz
      end do

!     Tranform cube coordinates into mesh structure
      msh%nNo = nx*ny*nz
      allocate(msh%x(nsd,msh%nNo))
      do k=1, nz
         do j=1, ny
            do i=1, nx
               Ac = (k-1)*ny*nx + (j-1)*nx + i
               msh%x(1,Ac) = x(i)
               msh%x(2,Ac) = y(j)
               msh%x(3,Ac) = z(k)
            end do
         end do
      end do

!     Define element structure and shape functions
      msh%eNoN = 4
      msh%nEl  = (nx-1)*(ny-1)*(nz-1)*6
      allocate(msh%IEN(msh%eNoN,msh%nEl))
      call selectele(msh)

!     Tetrahedron connectivity
      iord(:,1) = (/6, 4, 1, 2/)
      iord(:,2) = (/6, 4, 2, 3/)
      iord(:,3) = (/6, 4, 3, 7/)
      iord(:,4) = (/6, 4, 7, 8/)
      iord(:,5) = (/6, 4, 8, 5/)
      iord(:,6) = (/6, 4, 5, 1/)

!     Prepare msh%IEN array
      e = 0
      do k=1, nz-1
         do j=1, ny-1
            do i=1, nx-1
!              Get the discrete cube element's nodal indices
               nid(1) = (k-1)*nx*ny + (j-1)*nx + i
               nid(2) = nid(1) + 1
               nid(3) = nid(2) + nx
               nid(4) = nid(3) - 1

               nid(5) = nid(1) + (nx*ny)
               nid(6) = nid(5) + 1
               nid(7) = nid(6) + nx
               nid(8) = nid(7) - 1

!              Decompose cube element into 6 tetrahedrons
               do l=1, 6
                  e = e + 1
                  do a=1, 4
                     b = iord(a,l)
                     msh%IEN(a,e) = nid(b)
                  end do
               end do
            end do
         end do
      end do

!     Compute of volume of the cube
      allocate(xl(nsd,msh%eNoN), N(msh%eNoN), Nxi(nsd,msh%eNoN))
      vol = 0.0D0
      do e=1, msh%nEl
         do a=1, msh%eNoN
            Ac = msh%IEN(a,e)
            xl(:,a) = msh%x(:,Ac)
         end do

         do g=1, msh%nG
            if (g .EQ. 1) THEN
               call GNN(msh%eNoN, msh%Nx(:,:,g), xl, Nxi, Jac, Ks)
            end if
            N = msh%N(:,g)
            vol = vol + msh%w(g)*Jac
         end do
      end do
      write(stdout,ftab1) "Mesh <"//TRIM(msh%name)//"> vol: "// &
     &   TRIM(STR(vol))

      msh%nFa = 6
      allocate(msh%fa(msh%nFa))
      do iFa=1, msh%nFa
         msh%fa(iFa)%eNoN = 3
         call selecteleb(msh%fa(iFa))
      end do

!     Face Z = 0
      iFa = 1
      k   = 1
      nNo = nx*ny
      nEl = (nx-1)*(ny-1)*2
      msh%fa(iFa)%nNo = nNo
      msh%fa(iFa)%nEl = nEl
      write(msh%fa(iFa)%name,'(A)') "Z0"
      call faceAlloc(msh%fa(iFa))
      a = 0
      do j=1, ny
         do i=1, nx
            a  = a + 1
            Ac = (k-1)*nx*ny + (j-1)*nx + i
            msh%fa(iFa)%gN(a)  = Ac
            msh%fa(iFa)%x(:,a) = msh%x(:,Ac)
         end do
      end do

      e = 0
      do j=1, ny-1
         l = (nx-1)*(j-1)
         do i=1, nx-1
            e = e + 1
            msh%fa(iFa)%IEN(1,e) = (j-1)*nx + i + 1
            msh%fa(iFa)%IEN(2,e) = msh%fa(iFa)%IEN(1,e) + nx - 1
            msh%fa(iFa)%IEN(3,e) = msh%fa(iFa)%IEN(1,e) - 1
            msh%fa(iFa)%gE(e)    = 6*l + 1

            e = e + 1
            msh%fa(iFa)%IEN(1,e) = msh%fa(iFa)%IEN(2,e-1)
            msh%fa(iFa)%IEN(2,e) = msh%fa(iFa)%IEN(1,e-1)
            msh%fa(iFa)%IEN(3,e) = msh%fa(iFa)%IEN(2,e) + nx
            msh%fa(iFa)%gE(e)    = 6*l + 2
            l = l + 1
         end do
      end do

!     Face Z = 1
      iFa = 2
      k   = nz
      msh%fa(iFa)%nNo = nNo
      msh%fa(iFa)%nEl = nEl
      write(msh%fa(iFa)%name,'(A)') "Z1"
      call faceAlloc(msh%fa(iFa))
      a = 0
      do j=1, ny
         do i=1, nx
            a  = a + 1
            Ac = (k-1)*nx*ny + (j-1)*nx + i
            msh%fa(iFa)%gN(a)  = Ac
            msh%fa(iFa)%x(:,a) = msh%x(:,Ac)
         end do
      end do

      e = 0
      do j=1, ny-1
         l = (nx-1)*(ny-1)*(k-2) + (nx-1)*(j-1)
         do i=1, nx-1
            e = e + 1
            msh%fa(iFa)%IEN(1,e) = (j-1)*nx + i + 1
            msh%fa(iFa)%IEN(2,e) = msh%fa(iFa)%IEN(1,e) + nx - 1
            msh%fa(iFa)%IEN(3,e) = msh%fa(iFa)%IEN(2,e) + 1
            msh%fa(iFa)%gE(e)    = 6*l + 4

            e = e + 1
            msh%fa(iFa)%IEN(1,e) = msh%fa(iFa)%IEN(2,e-1)
            msh%fa(iFa)%IEN(2,e) = msh%fa(iFa)%IEN(1,e-1)
            msh%fa(iFa)%IEN(3,e) = msh%fa(iFa)%IEN(2,e) - 1
            msh%fa(iFa)%gE(e)    = 6*l + 5
            l = l + 1
         end do
      end do

!     Face X = 0
      iFa = 3
      i   = 1
      nNo = ny*nz
      nEl = (ny-1)*(nz-1)*2
      msh%fa(iFa)%nNo = nNo
      msh%fa(iFa)%nEl = nEl
      write(msh%fa(iFa)%name,'(A)') "X0"
      call faceAlloc(msh%fa(iFa))
      a = 0
      do k=1, nz
         do j=1, ny
            a  = a + 1
            Ac = (k-1)*nx*ny + (j-1)*nx + i
            msh%fa(iFa)%gN(a)  = Ac
            msh%fa(iFa)%x(:,a) = msh%x(:,Ac)
         end do
      end do

      e = 0
      do k=1, nz-1
         l = (nx-1)*(ny-1)*(k-1)
         do j=1, ny-1
            e = e + 1
            msh%fa(iFa)%IEN(1,e) = (k-1)*ny + j + ny
            msh%fa(iFa)%IEN(2,e) = (k-1)*ny + j + 1
            msh%fa(iFa)%IEN(3,e) = msh%fa(iFa)%IEN(1,e) + 1
            msh%fa(iFa)%gE(e)    = 6*l + 5

            e = e + 1
            msh%fa(iFa)%IEN(1,e) = msh%fa(iFa)%IEN(2,e-1)
            msh%fa(iFa)%IEN(2,e) = msh%fa(iFa)%IEN(1,e-1)
            msh%fa(iFa)%IEN(3,e) = msh%fa(iFa)%IEN(1,e) - 1
            msh%fa(iFa)%gE(e)    = 6*l + 6
            l = l + nx - 1
         end do
      end do

!     Face X = 1
      iFa = 4
      i   = nx
      nNo = ny*nz
      nEl = (ny-1)*(nz-1)*2
      msh%fa(iFa)%nNo = nNo
      msh%fa(iFa)%nEl = nEl
      write(msh%fa(iFa)%name,'(A)') "X1"
      call faceAlloc(msh%fa(iFa))
      a = 0
      do k=1, nz
         do j=1, ny
            a  = a + 1
            Ac = (k-1)*nx*ny + (j-1)*nx + i
            msh%fa(iFa)%gN(a)  = Ac
            msh%fa(iFa)%x(:,a) = msh%x(:,Ac)
         end do
      end do

      e = 0
      do k=1, nz-1
         l = (nx-1)*(ny-1)*(k-1) + (nx-2)
         do j=1, ny-1
            e = e + 1
            msh%fa(iFa)%IEN(1,e) = (k-1)*ny + j + 1
            msh%fa(iFa)%IEN(2,e) = msh%fa(iFa)%IEN(1,e) + ny - 1
            msh%fa(iFa)%IEN(3,e) = msh%fa(iFa)%IEN(2,e) + 1
            msh%fa(iFa)%gE(e)    = 6*l + 3

            e = e + 1
            msh%fa(iFa)%IEN(1,e) = msh%fa(iFa)%IEN(2,e-1)
            msh%fa(iFa)%IEN(2,e) = msh%fa(iFa)%IEN(1,e-1)
            msh%fa(iFa)%IEN(3,e) = msh%fa(iFa)%IEN(2,e) - 1
            msh%fa(iFa)%gE(e)    = 6*l + 2
            l = l + nx - 1
         end do
      end do

!     Face Y = 0
      iFa = 5
      j   = 1
      nNo = nx*nz
      nEl = (nx-1)*(nz-1)*2
      msh%fa(iFa)%nNo = nNo
      msh%fa(iFa)%nEl = nEl
      write(msh%fa(iFa)%name,'(A)') "Y0"
      call faceAlloc(msh%fa(iFa))
      a = 0
      do k=1, nz
         do i=1, nx
            a  = a + 1
            Ac = (k-1)*nx*ny + (j-1)*nx + i
            msh%fa(iFa)%gN(a)  = Ac
            msh%fa(iFa)%x(:,a) = msh%x(:,Ac)
         end do
      end do

      e = 0
      do k=1, nz-1
         l = (nx-1)*(ny-1)*(k-1)
         do i=1, nx-1
            e = e + 1
            msh%fa(iFa)%IEN(1,e) = (k-1)*nx + i
            msh%fa(iFa)%IEN(2,e) = msh%fa(iFa)%IEN(1,e) + nx + 1
            msh%fa(iFa)%IEN(3,e) = msh%fa(iFa)%IEN(1,e) + 1
            msh%fa(iFa)%gE(e)    = 6*l + 1

            e = e + 1
            msh%fa(iFa)%IEN(1,e) = msh%fa(iFa)%IEN(2,e-1)
            msh%fa(iFa)%IEN(2,e) = msh%fa(iFa)%IEN(1,e-1)
            msh%fa(iFa)%IEN(3,e) = msh%fa(iFa)%IEN(2,e) + nx
            msh%fa(iFa)%gE(e)    = 6*l + 6
            l = l + 1
         end do
      end do

!     Face Y = 1
      iFa = 6
      j   = ny
      nNo = nx*nz
      nEl = (nx-1)*(nz-1)*2
      msh%fa(iFa)%nNo = nNo
      msh%fa(iFa)%nEl = nEl
      write(msh%fa(iFa)%name,'(A)') "Y1"
      call faceAlloc(msh%fa(iFa))
      a = 0
      do k=1, nz
         do i=1, nx
            a  = a + 1
            Ac = (k-1)*nx*ny + (j-1)*nx + i
            msh%fa(iFa)%gN(a)  = Ac
            msh%fa(iFa)%x(:,a) = msh%x(:,Ac)
         end do
      end do

      e = 0
      do k=1, nz-1
         l = (nx-1)*(ny-1)*(k-1) + (nx-1)*(ny-2)
         do i=1, nx-1
            e = e + 1
            msh%fa(iFa)%IEN(1,e) = (k-1)*nx + i + nx + 1
            msh%fa(iFa)%IEN(2,e) = (k-1)*nx + i
            msh%fa(iFa)%IEN(3,e) = msh%fa(iFa)%IEN(2,e) + 1
            msh%fa(iFa)%gE(e)    = 6*l + 3

            e = e + 1
            msh%fa(iFa)%IEN(1,e) = msh%fa(iFa)%IEN(2,e-1)
            msh%fa(iFa)%IEN(2,e) = msh%fa(iFa)%IEN(1,e-1)
            msh%fa(iFa)%IEN(3,e) = msh%fa(iFa)%IEN(2,e) - 1
            msh%fa(iFa)%gE(e)    = 6*l + 4
            l = l + 1
         end do
      end do

      call VTK(msh)

      call destroy(msh)

      return
      end program tet_mesher

!**************************************************

      subroutine faceAlloc(lFa)
      use varmod
      implicit none
      type(faceType), intent(inout) :: lFa

      allocate(lFa%x(nsd,lFa%nNo))
      allocate(lFa%IEN(lFa%eNoN,lFa%nEl))
      allocate(lFa%gN(lFa%nNo))
      allocate(lFa%gE(lFa%nEl))

      return
      end subroutine faceAlloc

!**************************************************

      subroutine selectele(lM)
      use varmod
      implicit none
      type(mshType), intent(inout) :: lM

      integer :: g

      if (nsd .eq. 3) then
         select case (lM%eNoN)
         case (8)
            lM%eType   = eType_BRK
            lM%nG      = 8
            lM%vtkType = 12
         case (6)
            lM%eType   = eType_WDG
            lM%nG      = 6
            lM%vtkType = 13
         case (4)
            lM%eType   = eType_TET
            lM%nG      = 4
            lM%vtkType = 10
         case default
            write(stdout,ftab4) &
               "Error: unknown combination of nsd & eNoN"
            STOP
         end select
      else
         select case (lM%eNoN)
         case (3)
            lM%eType   = eType_TRI
            lM%nG      = 3
            lM%vtkType = 5
         case (4)
            lM%eType   = eType_BIL
            lM%nG      = 4
            lM%vtkType = 9
         case (9)
            lM%eType   = eType_BIQ
            lM%nG      = 9
            lM%vtkType = 28
         case default
            write(stdout,ftab4) &
               "Error: unknown combination of nsd & eNoN"
            STOP
         end select
      end if

      allocate(lM%w(lM%nG), lM%xi(nsd,lM%nG), lM%N(lM%eNoN,lM%nG), &
         lM%Nx(nsd,lM%eNoN,lM%nG))

      call getGP(nsd, lM%eType, lM%nG, lM%w, lM%xi)

      do g=1, lM%nG
         call getShpF(nsd, lM%eType, lM%eNoN, lM%xi(:,g), lM%N(:,g), &
            lM%Nx(:,:,g))
      end do

      return
      end subroutine selectele

!**************************************************

      subroutine selecteleb(lFa)
      use varmod
      implicit none
      type(faceType), intent(inout) :: lFa

      integer :: g, insd

      insd = nsd - 1
      if (insd .eq. 2) then
         select case (lFa%eNoN)
         case (4)
            lFa%eType   = eType_BIL
            lFa%nG      = 4
            lFa%vtkType = 9
         case (3)
            lFa%eType   = eType_TRI
            lFa%nG      = 3
            lFa%vtkType = 5
         case default
            write(stdout,ftab4) &
               "Error: unknown combination of nsd & eNoNb"
            STOP
         end select
      else
         select case (lFa%eNoN)
         case (2)
            lFa%eType   = eType_LIN
            lFa%nG      = 2
            lFa%vtkType = 3
         case (1)
            lFa%eType   = eType_QUD
            lFa%nG      = 3
            lFa%vtkType = 21
         case default
            write(stdout,ftab4) &
               "Error: unknown combination of nsd & eNoNb"
            STOP
         end select
      end if

      allocate(lFa%w(lFa%nG), lFa%xi(nsd,lFa%nG), &
     &   lFa%N(lFa%eNoN,lFa%nG), lFa%Nx(insd,lFa%eNoN,lFa%nG))

      call getGP(insd, lFa%eType, lFa%nG, lFa%w, lFa%xi)

      do g=1, lFa%nG
         call getShpF(insd, lFa%eType, lFa%eNoN, lFa%xi(:,g), &
     &      lFa%N(:,g), lFa%Nx(:,:,g))
      end do

      return
      end subroutine selecteleb

!**************************************************

      subroutine getGP(insd, eType, nG, w, xi)
      use varmod
      implicit none
      integer, intent(in) :: insd, eType, nG
      real(kind=8), intent(out) :: w(nG), xi(insd,nG)

      real(kind=8) s, t, lz, uz

!     3D elements
      select case (eType)
      case (eType_BRK)
         w = 1D0
         s =  1D0/SQRT(3D0)
         t = -1D0/SQRT(3D0)
         xi(1,1) = s; xi(2,1) = s; xi(3,1) = s
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = s
         xi(1,3) = t; xi(2,3) = s; xi(3,3) = t
         xi(1,4) = s; xi(2,4) = s; xi(3,4) = t
         xi(1,5) = s; xi(2,5) = t; xi(3,5) = s
         xi(1,6) = t; xi(2,6) = t; xi(3,6) = s
         xi(1,7) = t; xi(2,7) = t; xi(3,7) = t
         xi(1,8) = s; xi(2,8) = t; xi(3,8) = t
      case (eType_TET)
         w = 1D0/24D0
         s = (5D0 + 3D0*SQRT(5D0))/2D1
         t = (5D0 -     SQRT(5D0))/2D1
         xi(1,1) = s; xi(2,1) = t; xi(3,1) = t
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = t
         xi(1,3) = t; xi(2,3) = t; xi(3,3) = s
         xi(1,4) = t; xi(2,4) = t; xi(3,4) = t
      case (eType_WDG)
         w  =  1D0/6D0
         s  =  2D0/3D0
         t  =  1D0/6D0
         uz =  1D0/SQRT(3D0)
         lz = -1D0/SQRT(3D0)
         xi(1,1) = s; xi(2,1) = t; xi(3,1) = lz
         xi(1,2) = t; xi(2,2) = s; xi(3,2) = lz
         xi(1,3) = t; xi(2,3) = t; xi(3,3) = lz
         xi(1,4) = s; xi(2,4) = t; xi(3,4) = uz
         xi(1,5) = t; xi(2,5) = s; xi(3,5) = uz
         xi(1,6) = t; xi(2,6) = t; xi(3,6) = uz

!     2D elements
      case (eType_TRI)
         w = 1D0/6D0
         s = 2D0/3D0
         t = 1D0/6D0
         xi(1,1) = s; xi(2,1) = t
         xi(1,2) = t; xi(2,2) = s
         xi(1,3) = t; xi(2,3) = t
      case (eType_BIL)
         w = 1D0
         s =  1D0/SQRT(3D0)
         t = -1D0/SQRT(3D0)
         xi(1,1) = s; xi(2,1) = s
         xi(1,2) = t; xi(2,2) = s
         xi(1,3) = t; xi(2,3) = t
         xi(1,4) = s; xi(2,4) = t
      case (eType_BIQ)
         w(1) = 25D0/81D0; w(2) = 25D0/81D0; w(3) = 25D0/81D0
         w(4) = 25D0/81D0; w(5) = 40D0/81D0; w(6) = 40D0/81D0
         w(7) = 40D0/81D0; w(8) = 40D0/81D0; w(9) = 64D0/81D0
         s    = SQRT(6D-1)
         xi(1,1) =  -s; xi(2,1) =  -s
         xi(1,2) =   s; xi(2,2) =  -s
         xi(1,3) =   s; xi(2,3) =   s
         xi(1,4) =  -s; xi(2,4) =   s
         xi(1,5) = 0D0; xi(2,5) =  -s
         xi(1,6) =   s; xi(2,6) = 0D0
         xi(1,7) = 0D0; xi(2,7) =   s
         xi(1,8) =  -s; xi(2,8) = 0D0
         xi(1,9) = 0D0; xi(2,9) = 0D0

!     1D elements
      case (eType_LIN)
         w = 1D0
         s = 1D0/SQRT(3D0)
         xi(1,1) = -s
         xi(1,2) =  s
      case (eType_QUD)
         w(1) = 5D0/9D0; w(2) = 5D0/9D0; w(3) = 8D0/9D0
         s = SQRT(6D-1)
         xi(1,1) = -s
         xi(1,2) =  s
         xi(1,3) = 0D0
      end select

      return
      end subroutine getGP

!**************************************************

      subroutine getShpF(insd, eType, eNoN, xi, N, Nxi)
      use varmod
      implicit none
      integer, intent(in) :: insd, eType, eNoN
      real(kind=8), intent(out) :: xi(insd), N(eNoN), Nxi(insd,eNoN)

      real(kind=8) :: s, t, mx, my, ux, uy, uz, lx, ly, lz

!     3D elements
      select case (eType)
      case (eType_BRK)
         ux = 1D0 + xi(1); lx = 1D0 - xi(1)
         uy = 1D0 + xi(2); ly = 1D0 - xi(2)
         uz = 1D0 + xi(3); lz = 1D0 - xi(3)
         N(1) = ux*uy*uz/8D0
         N(2) = lx*uy*uz/8D0
         N(3) = lx*uy*lz/8D0
         N(4) = ux*uy*lz/8D0
         N(5) = ux*ly*uz/8D0
         N(6) = lx*ly*uz/8D0
         N(7) = lx*ly*lz/8D0
         N(8) = ux*ly*lz/8D0

         Nxi(1,1) =  uy*uz/8D0
         Nxi(2,1) =  ux*uz/8D0
         Nxi(3,1) =  ux*uy/8D0
         Nxi(1,2) = -uy*uz/8D0
         Nxi(2,2) =  lx*uz/8D0
         Nxi(3,2) =  lx*uy/8D0
         Nxi(1,3) = -uy*lz/8D0
         Nxi(2,3) =  lx*lz/8D0
         Nxi(3,3) = -lx*uy/8D0
         Nxi(1,4) =  uy*lz/8D0
         Nxi(2,4) =  ux*lz/8D0
         Nxi(3,4) = -ux*uy/8D0
         Nxi(1,5) =  ly*uz/8D0
         Nxi(2,5) = -ux*uz/8D0
         Nxi(3,5) =  ux*ly/8D0
         Nxi(1,6) = -ly*uz/8D0
         Nxi(2,6) = -lx*uz/8D0
         Nxi(3,6) =  lx*ly/8D0
         Nxi(1,7) = -ly*lz/8D0
         Nxi(2,7) = -lx*lz/8D0
         Nxi(3,7) = -lx*ly/8D0
         Nxi(1,8) =  ly*lz/8D0
         Nxi(2,8) = -ux*lz/8D0
         Nxi(3,8) = -ux*ly/8D0

      case (eType_TET)
         N(1) = xi(1)
         N(2) = xi(2)
         N(3) = xi(3)
         N(4) = 1D0 - xi(1) - xi(2) - xi(3)

         Nxi(1,1) =  1D0
         Nxi(2,1) =  0D0
         Nxi(3,1) =  0D0
         Nxi(1,2) =  0D0
         Nxi(2,2) =  1D0
         Nxi(3,2) =  0D0
         Nxi(1,3) =  0D0
         Nxi(2,3) =  0D0
         Nxi(3,3) =  1D0
         Nxi(1,4) = -1D0
         Nxi(2,4) = -1D0
         Nxi(3,4) = -1D0

      case (eType_WDG)
         ux = xi(1) ; uy = xi(2) ; uz = 1D0 - ux - uy
         s = (1D0 + xi(3))/2D0; t = (1D0 - xi(3))/2D0
         N(1) = ux*t
         N(2) = uy*t
         N(3) = uz*t
         N(4) = ux*s
         N(5) = uy*s
         N(6) = uz*s

         Nxi(1,1) =  t
         Nxi(2,1) =  0D0
         Nxi(3,1) = -ux/2D0
         Nxi(1,2) =  0D0
         Nxi(2,2) =  t
         Nxi(3,2) = -uy/2D0
         Nxi(1,3) = -t
         Nxi(2,3) = -t
         Nxi(3,3) = -uz/2D0
         Nxi(1,4) =  s
         Nxi(2,4) =  0D0
         Nxi(3,4) =  ux/2D0
         Nxi(1,5) =  0D0
         Nxi(2,5) =  s
         Nxi(3,5) =  uy/2D0
         Nxi(1,6) = -s
         Nxi(2,6) = -s
         Nxi(3,6) =  uz/2D0

!     2D elements
      case (eType_TRI)
         N(1) = xi(1)
         N(2) = xi(2)
         N(3) = 1D0 - xi(1) - xi(2)

         Nxi(1,1) =  1D0
         Nxi(2,1) =  0D0
         Nxi(1,2) =  0D0
         Nxi(2,2) =  1D0
         Nxi(1,3) = -1D0
         Nxi(2,3) = -1D0

      case (eType_BIL)
         ux = 1D0 + xi(1); lx = 1D0 - xi(1)
         uy = 1D0 + xi(2); ly = 1D0 - xi(2)
         N(1) = ux*uy/4D0
         N(2) = lx*uy/4D0
         N(3) = lx*ly/4D0
         N(4) = ux*ly/4D0

         Nxi(1,1) =  uy/4D0
         Nxi(2,1) =  ux/4D0
         Nxi(1,2) = -uy/4D0
         Nxi(2,2) =  lx/4D0
         Nxi(1,3) = -ly/4D0
         Nxi(2,3) = -lx/4D0
         Nxi(1,4) =  ly/4D0
         Nxi(2,4) = -ux/4D0

      case (eType_BIQ)
         ux = 1D0 + xi(1); mx = xi(1); lx = 1D0 - xi(1)
         uy = 1D0 + xi(2); my = xi(2); ly = 1D0 - xi(2)
         N(1) =  mx*lx*my*ly/4D0
         N(2) = -mx*ux*my*ly/4D0
         N(3) =  mx*ux*my*uy/4D0
         N(4) = -mx*lx*my*uy/4D0
         N(5) = -lx*ux*my*ly/2D0
         N(6) =  mx*ux*ly*uy/2D0
         N(7) =  lx*ux*my*uy/2D0
         N(8) = -mx*lx*ly*uy/2D0
         N(9) =  lx*ux*ly*uy

         Nxi(1,1) =  (lx - mx)*my*ly/4D0
         Nxi(2,1) =  (ly - my)*mx*lx/4D0
         Nxi(1,2) = -(ux + mx)*my*ly/4D0
         Nxi(2,2) = -(ly - my)*mx*ux/4D0
         Nxi(1,3) =  (ux + mx)*my*uy/4D0
         Nxi(2,3) =  (uy + my)*mx*ux/4D0
         Nxi(1,4) = -(lx - mx)*my*uy/4D0
         Nxi(2,4) = -(uy + my)*mx*lx/4D0
         Nxi(1,5) = -(lx - ux)*my*ly/2D0
         Nxi(2,5) = -(ly - my)*lx*ux/2D0
         Nxi(1,6) =  (ux + mx)*ly*uy/2D0
         Nxi(2,6) =  (ly - uy)*mx*ux/2D0
         Nxi(1,7) =  (lx - ux)*my*uy/2D0
         Nxi(2,7) =  (uy + my)*lx*ux/2D0
         Nxi(1,8) = -(lx - mx)*ly*uy/2D0
         Nxi(2,8) = -(ly - uy)*mx*lx/2D0
         Nxi(1,9) =  (lx - ux)*ly*uy
         Nxi(2,9) =  (ly - uy)*lx*ux

!     1D elements
      case (eType_LIN)
         N(1) = (1D0 - xi(1))/2D0
         N(2) = (1D0 + xi(1))/2D0

         Nxi(1,1) = -5D-1
         Nxi(1,2) =  5D-1
      case (eType_QUD)
         N(1) = -xi(1)*(1D0 - xi(1))/2D0
         N(2) =  xi(1)*(1D0 + xi(1))/2D0
         N(3) = (1D0 - xi(1))*(1D0 + xi(1))

         Nxi(1,1) = -5D-1 + xi(1)
         Nxi(1,2) =  5D-1 + xi(1)
         Nxi(1,3) = -2D0*xi(1)
      end select

      return
      end subroutine getShpF

!**************************************************

      subroutine GNN(eNoN, Nxi, x, Nx, Jac, ks)
      use varmod, only: nsd
      implicit none

      integer, intent(in) :: eNoN
      real(kind=8), intent(in) :: Nxi(nsd,eNoN), x(nsd,eNoN)
      real(kind=8), intent(out) :: Nx(nsd,eNoN), Jac, ks(nsd,nsd)

      integer :: a
      real(kind=8) :: xXi(nsd,nsd), xiX(nsd,nsd)

      Nx  = 0D0
      xXi = 0D0
      if (nsd .eq. 2) then
         do a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
         end do

         Jac = xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1)

         xiX(1,1) =  xXi(2,2)/Jac
         xiX(1,2) = -xXi(1,2)/Jac
         xiX(2,1) = -xXi(2,1)/Jac
         xiX(2,2) =  xXi(1,1)/Jac

         ks(1,1) = xiX(1,1)*xiX(1,1) + xiX(2,1)*xiX(2,1)
         ks(1,2) = xiX(1,1)*xiX(1,2) + xiX(2,1)*xiX(2,2)
         ks(2,2) = xiX(1,2)*xiX(1,2) + xiX(2,2)*xiX(2,2)
         ks(2,1) = ks(1,2)

         do a=1, eNoN
            Nx(1,a) = Nx(1,a)+ Nxi(1,a)*xiX(1,1) + Nxi(2,a)*xiX(2,1)
            Nx(2,a) = Nx(2,a)+ Nxi(1,a)*xiX(1,2) + Nxi(2,a)*xiX(2,2)
         end do
      else
         do a=1, eNoN
            xXi(:,1) = xXi(:,1) + x(:,a)*Nxi(1,a)
            xXi(:,2) = xXi(:,2) + x(:,a)*Nxi(2,a)
            xXi(:,3) = xXi(:,3) + x(:,a)*Nxi(3,a)
         end do

         Jac = xXi(1,1)*xXi(2,2)*xXi(3,3) + &
               xXi(1,2)*xXi(2,3)*xXi(3,1) + &
               xXi(1,3)*xXi(2,1)*xXi(3,2) - &
               xXi(1,1)*xXi(2,3)*xXi(3,2) - &
               xXi(1,2)*xXi(2,1)*xXi(3,3) - &
               xXi(1,3)*xXi(2,2)*xXi(3,1)

         xiX(1,1) = (xXi(2,2)*xXi(3,3) - xXi(2,3)*xXi(3,2))/Jac
         xiX(1,2) = (xXi(3,2)*xXi(1,3) - xXi(3,3)*xXi(1,2))/Jac
         xiX(1,3) = (xXi(1,2)*xXi(2,3) - xXi(1,3)*xXi(2,2))/Jac
         xiX(2,1) = (xXi(2,3)*xXi(3,1) - xXi(2,1)*xXi(3,3))/Jac
         xiX(2,2) = (xXi(3,3)*xXi(1,1) - xXi(3,1)*xXi(1,3))/Jac
         xiX(2,3) = (xXi(1,3)*xXi(2,1) - xXi(1,1)*xXi(2,3))/Jac
         xiX(3,1) = (xXi(2,1)*xXi(3,2) - xXi(2,2)*xXi(3,1))/Jac
         xiX(3,2) = (xXi(3,1)*xXi(1,2) - xXi(3,2)*xXi(1,1))/Jac
         xiX(3,3) = (xXi(1,1)*xXi(2,2) - xXi(1,2)*xXi(2,1))/Jac

         ks(1,1) = xiX(1,1)*xiX(1,1)+xiX(2,1)*xiX(2,1)+xiX(3,1)*xiX(3,1)
         ks(1,2) = xiX(1,2)*xiX(1,1)+xiX(2,2)*xiX(2,1)+xiX(3,2)*xiX(3,1)
         ks(1,3) = xiX(1,3)*xiX(1,1)+xiX(2,3)*xiX(2,1)+xiX(3,3)*xiX(3,1)
         ks(2,2) = xiX(1,2)*xiX(1,2)+xiX(2,2)*xiX(2,2)+xiX(3,2)*xiX(3,2)
         ks(2,3) = xiX(1,2)*xiX(1,3)+xiX(2,2)*xiX(2,3)+xiX(3,2)*xiX(3,3)
         ks(3,3) = xiX(1,3)*xiX(1,3)+xiX(2,3)*xiX(2,3)+xiX(3,3)*xiX(3,3)
         ks(2,1) = ks(1,2)
         ks(3,1) = ks(1,3)
         ks(3,2) = ks(2,3)

         do a=1, eNoN
            Nx(1,a) = Nx(1,a) + Nxi(1,a)*xiX(1,1) + &
                                Nxi(2,a)*xiX(2,1) + &
                                Nxi(3,a)*xiX(3,1)

            Nx(2,a) = Nx(2,a) + Nxi(1,a)*xiX(1,2) + &
                                Nxi(2,a)*xiX(2,2) + &
                                Nxi(3,a)*xiX(3,2)

            Nx(3,a) = Nx(3,a) + Nxi(1,a)*xiX(1,3) + &
                                Nxi(2,a)*xiX(2,3) + &
                                Nxi(3,a)*xiX(3,3)
         end do
      end if

      return
      end subroutine GNN

!***********************************************************************

      subroutine VTK(lM)
      use varmod
      use vtkXMLMod
      implicit none

      type(mshType), intent(inout) :: lM
      integer :: iFa
      character(len=strL) :: fName, fdir
      logical :: flag

!     Write mesh vtu file
      if (nx .le. 10) then
         write(fdir,'(A)') "./N00"//TRIM(STR(nx-1))
      else if (nx.gt.10 .and. nx.le.100) then
         write(fdir,'(A)') "./N0"//TRIM(STR(nx-1))
      else
         write(fdir,'(A)') "./N"//TRIM(STR(nx-1))
      end if
      call system('mkdir -p '//TRIM(fdir))
      write(fName,'(A)') TRIM(fdir)//"/mesh-complete.mesh.vtu"
      lM%IEN(:,:) = lM%IEN(:,:) - 1
      call writeVTU(lM, fName)
      lM%IEN(:,:) = lM%IEN(:,:) + 1

!     Write face vtp files
      if (allocated(msh%fa)) then
         inquire(file=TRIM(fdir)//"/mesh-surfaces",exist=flag)
         if (.not.flag) call system("mkdir  "//TRIM(fdir)// &
     &      "/mesh-surfaces")
         do iFa=1, lM%nFa
            write(fName,'(A)') TRIM(fdir)//"/mesh-surfaces/"// &
     &         TRIM(lM%fa(iFa)%name)//".vtp"
            lM%fa(iFa)%IEN = lM%fa(iFa)%IEN - 1
            call writeVTP(lM%fa(iFa), fName)
            lM%fa(iFa)%IEN = lM%fa(iFa)%IEN + 1
         end do
      end if

      return
      end subroutine VTK

!**********************************************************************

      subroutine writeVTU(lM, fName)
      use varmod
      use vtkXMLMod

      implicit none

      type(mshType), intent(in) :: lM
      character(len=strL), intent(in) :: fName

      type(vtkXMLType) :: vtu
      integer :: iStat

      call vtkInitWriter(vtu, TRIM(fName), iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTU file write error (init)"
         stop
      end if

      call putVTK_pointCoords(vtu, lM%x, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTU file write error (coords)"
         stop
      end if

      call putVTK_elemIEN(vtu, lM%IEN, lM%vtkType, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTU file write error (ien)"
         stop
      end if

      call vtkWriteToFile(vtu, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTU file write error"
         stop
      end if

      call flushVTK(vtu)

      return
      end subroutine writeVTU

!**********************************************************************

      subroutine writeVTP(lFa, fName)
      use varmod
      use vtkXMLMod
      implicit none

      type(faceType), intent(in) :: lFa
      character(len=strL), intent(in) :: fName

      type(vtkXMLType) :: vtp
      integer :: iStat

      call vtkInitWriter(vtp, TRIM(fName), iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error (init)"
         stop
      end if

      call putVTK_pointCoords(vtp, lFa%x, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error (coords)"
         stop
      end if

      call putVTK_elemIEN(vtp, lFa%IEN, lFa%vtkType, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error (ien)"
         stop
      end if

      if (allocated(lFa%gN)) then
         call putVTK_pointData(vtp, "GlobalNodeID", lFa%gN, iStat)
         if (iStat .lt. 0) then
            write(stdout,ftab4) &
     &         "ERROR: VTP file write error (point data)"
            stop
         end if
      end if

      if (allocated(lFa%gE)) then
         call putVTK_elemData(vtp, "GlobalElementID", lFa%gE, iStat)
         if (iStat .lt. 0) then
            write(stdout,ftab4) &
     &         "ERROR: VTP file write error (cell data)"
            stop
         end if
      end if

      call vtkWriteToFile(vtp, iStat)
      if (iStat .lt. 0) then
         write(stdout,ftab4) "ERROR: VTP file write error"
         stop
      end if

      call flushVTK(vtp)

      return
      end subroutine writeVTP

!***********************************************************************

