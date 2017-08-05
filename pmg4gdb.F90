!............................................................
!. Parallel Multigrid for GDB
!. This is a sample code showing the use of some of the
!. multigrid solvers in the PMGsolver.F90 file
!.
!. Manaure Francisquez
!. 2016
!.
!............................................................
module pars
!.Include the mpi header file
  include 'mpif.h'
!.Include user-defined preprocessor variables and directives
#include "gdbFLAGS.h"

!.Fundamental constants (these never never change):
  real*8, parameter :: pi=4.0d0*atan(1.0_8)
  real*8, parameter :: twopi=2.0d0*pi

!.Program name. Length=exact # of characters used; e.g. character(3) for f3d
  character(3), parameter :: progname = 'gdb'

!.Box size (does not include ghost cells).
! In periodic directions, y(ny0+j)=y(j)+Ly for example
  real*8 :: Lx, Ly, Lz

!.Number of iteration (better time counting)
  integer, parameter :: noit=1

!.Number of MPI blocks along z, y and x (use powers of 2).
!.Total number of processes is MPIblocks(1)*MPIblocks(2)*MPIblocks(3)
!.IMPORTANT:
! 1.Write in rwo-major order (MPI follows C indexing)
!.2.Make nx0G/xprocs,ny0G/yprocs>=4
  integer :: MPIgrid(3)

!.Number of OpenMP threads:
  integer :: OMPthreads    !.Currently OpenMP is not included

!.Global physical grid size of the entire simulation (use powers of 2):
! nx0G=xprocs*nx0, nx0 is the local number of x-points. Likewise for y,z
  integer :: nx0G,ny0G,nz0G

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!..No user-set parameters below this line:
!
!

!.Local number of grid points along x, y and z
  integer :: nx0,ny0,nz0

!.for MPI and ffts:
  integer :: xprocs,yprocs,zprocs,nprocs
  integer :: xprocsm1
  integer :: nx0G2,ny0G2,nz0G2

!.Grid spacing:
!.Assuming x(1)=0.-Lx/2+dx/2 and x(nx0)=Lx/2-dx/2 gives
!.nx0 points between x=0,Lx and thus nx0-1 intervals:
  real*8 :: dx,dy,dz

!.spatial grid:
  real*8, dimension(:), allocatable :: x,y,z,xblock,yblock,zblock

  integer :: nt,ntt    !.time step and frame counters

!.MPI variables
  integer :: ierr,ierror,irc
  integer :: istat(MPI_STATUS_SIZE)
  integer :: myrank,numprocs,nprovided
!.These are used by MPI_CART routines
  integer :: WCOMM,CARTrank
  logical :: CART_BCs(3),reorder=.TRUE.
  integer :: COMM_CART,blockID(3)
!.These are used for creating and using subcommunicators
  integer :: COMM_X,COMM_Y,COMM_XY,remain(3), &
             xID,yID,zID,xyID,xcoord,ycoord,zcoord,xycoord(2)
!.Rank of process to the left and right (along x).
  integer :: left, right
!.Rank of process up and down (along y).
  integer :: up, down

!.Torque scheduler variables
  integer :: args=0
  character(len=8) :: cjobid, crestartid

!.Diffusion coefficients (along x and y) for each quantity
  real*8, allocatable :: Diffs(:,:)

!.Input parameters:
!.Declare here all variables supplied by the input file
!.(add to namelists and BCAST below):
  integer :: nts,nframes,nrst,nfdump
  real*8 :: dt, &
            diffwx,diffpsix,diffnx,diffTex,diffTix,diffvx, &
            diffwy,diffpsiy,diffny,diffTey,diffTiy,diffvy
  real*8 :: de2
#IF (UTIMERS == 1)
!.Timing counters definded by user
  real*8 :: mgloop,mgc1,mgc2
#ENDIF

contains

!*********************************************************
  subroutine compute_pars
  implicit none
!.Number of processes along each direction and total number of processes
  MPIgrid=(/zprocs,yprocs,xprocs/)
!.Local number of grid points along x and z
  nx0=nx0G/xprocs
  ny0=ny0G/yprocs
  nz0=nz0G/zprocs

!.for MPI and ffts:
  nprocs=xprocs*yprocs*zprocs
  xprocsm1=xprocs-1
  nx0G2=nx0G/2+1
  ny0G2=ny0G/2+1
  nz0G2=nz0G/2+1

!.Grid spacing
  dx=Lx/dble(nx0G)
  dy=Ly/dble(ny0G)
  dz=Lz/dble(nz0G)
  end subroutine compute_pars
!*********************************************************
  subroutine alloc_pars
  implicit none
  integer :: i,j,k,ierr
  allocate(x(nx0G),y(ny0G),z(nz0G),xblock(nx0),yblock(ny0),zblock(nz0))
  forall(i=1:nx0G) x(i)=(dble(i-nx0G2)+0.5)*dx
  forall(j=1:ny0G) y(j)=dble(j-ny0G2)*dy
  forall(k=1:nz0G) z(k)=dble(k-nz0G2)*dz

  end subroutine alloc_pars
!*********************************************************
  subroutine read_inputfile
  implicit none
  real*8 :: fLx,fLy,fLz
!.Must call this to get the input params

!.Add here the namelists for parameters (first declared above) from input file:
  namelist/dom/nx0G,ny0G,nz0G,xprocs,yprocs,zprocs,fLx,fLy,fLz
  namelist/comset/OMPthreads
  namelist/list/nts,dt,nframes,nrst,nfdump
  namelist/diff/de2, &
                diffwx,diffpsix,diffnx,difftex,difftix,diffvx, &
                diffwy,diffpsiy,diffny,difftey,difftiy,diffvy

  if (myrank .eq. 0) then
    open(unit=10,file="gdb.in",status="old")
    read (10,dom)
    read (10,comset)
    read (10,list)
    read (10,diff)
    close(unit=10)

    Lx = fLx*twopi
    Ly = fLy*twopi
    Lz = fLz*twopi
  endif

  call MPI_BCAST(nx0G,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ny0G,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nz0G,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(xprocs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(yprocs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(zprocs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Lx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Ly,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(Lz,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(OMPthreads,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nts,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nframes,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nrst,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(nfdump,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(de2,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(diffwx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(diffpsix,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(diffnx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(difftex,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(difftix,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(diffvx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(diffwy,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(diffpsiy,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(diffny,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(difftey,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(difftiy,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(diffvy,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  end subroutine read_inputfile
!*********************************************************
  subroutine dealloc_pars
  implicit none
  deallocate(x,y,z,xblock,yblock,zblock)

  end subroutine dealloc_pars

end module pars
!*****************************************************************
!*****************************************************************
module fields
  use pars
  real*8, dimension(:,:,:), allocatable :: curi,psii
  real*8, dimension(:,:,:), allocatable :: lden,ldeni,deni
  real*8, dimension(:,:,:), allocatable :: lte,ltei
  real*8, dimension(:,:,:), allocatable :: lti,ltii
  real*8, dimension(:,:,:), allocatable :: w,wi,phi,phi0,vort
  real*8, dimension(:,:,:), allocatable :: rpsi,rpsii
  real*8, dimension(:,:,:), allocatable :: psi,psi0
  real*8, dimension(:,:,:), allocatable :: vp,vpi
  real*8, dimension(:), allocatable :: denav
  real*8, dimension(:,:), allocatable :: denyav
contains
  subroutine alloc_fields
  implicit none
  allocatE(curi(nx0,ny0,nz0),psii(nx0,ny0,nz0))
  allocate(lden(nx0,ny0,nz0),ldeni(nx0,ny0,nz0),deni(nx0,ny0,nz0))
  allocate(lte(nx0,ny0,nz0),ltei(nx0,ny0,nz0))
  allocate(lti(nx0,ny0,nz0),ltii(nx0,ny0,nz0))
  allocate(w(nx0,ny0,nz0),wi(nx0,ny0,nz0),phi(nx0,ny0,nz0), &
           phi0(nx0,ny0,nz0),vort(nx0,ny0,nz0))
  allocate(rpsi(nx0,ny0,nz0),rpsii(nx0,ny0,nz0))
  allocate(psi(nx0,ny0,nz0),psi0(nx0,ny0,nz0))
  allocate(vp(nx0,ny0,nz0),vpi(nx0,ny0,nz0))
  allocate(denav(nx0))
  allocate(denyav(nx0,nz0))

  end subroutine alloc_fields
!*********************************************************
  subroutine dealloc_fields
  implicit none
  deallocate(curi,psii)
  deallocate(lden,ldeni,deni)
  deallocate(lte,ltei)
  deallocate(lti,ltii)
  deallocate(w,wi,phi,phi0,vort)
  deallocate(rpsi,rpsii)
  deallocate(psi,psi0)
  deallocate(vp,vpi)
  deallocate(denav)
  deallocate(denyav)
  end subroutine dealloc_fields
end module fields
!*****************************************************************
!*****************************************************************
program pmgdiff
  use PMGsolver
  use pars
  use fields
  implicit none
  real*8, allocatable :: phir(:,:)
  real(DP) :: kx0,ky0,x0,kp,lambda0,sigs
  real(DP), allocatable, dimension(:,:,:) :: lam
!.Diffusion coefficients
  real*8 :: kkxmaxsq,kymaxsq
  real*8 :: kxmax,kymax,knorm,diffexp
  real*8, allocatable :: hdcos(:,:),wrhs(:,:,:)
  real*8 :: mgt1,mgt2
  integer :: i,j,k
!.For MPI stuff
  integer :: IDs(3),ID2s(2)
  real*8, allocatable :: res(:,:,:)
  integer :: pdiff=2
!.for looping over solver
  real*8, allocatable :: wl(:,:,:,:),wrhsl(:,:,:,:)


!.Initialize MPI
  call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,nprovided,ierr)
  call MPI_COMM_DUP(MPI_COMM_WORLD,WCOMM,ierr)
!.Who am I? --- get my rank
  call MPI_COMM_RANK(WCOMM,myrank,ierr)
!.How many processes in the global group?
  call MPI_COMM_SIZE(WCOMM,numprocs,ierr)

!.read the input parameters:
  call read_inputfile
  call compute_pars
  call alloc_pars

!.Check that the number of processes launched is correct
  if (numprocs .ne. nprocs) then
    if (myrank.eq.0) then
      print*,'nprocs parameter in code differs from that in mpirun -np xxx'
      print*,'nprocs in code is',nprocs
    endif
    call MPI_BARRIER(WCOMM,ierr)
    go to 99
  endif

  call alloc_fields

!.MPI_CART boundary conditions along x and z. True=periodic
  CART_BCs(1) = .TRUE.   ! periodic along z
  CART_BCs(2) = .TRUE.   ! periodic along y
  CART_BCs(3) = .FALSE.  ! not periodic along x

!.Create a Cartesian topology
  call MPI_CART_CREATE(MPI_COMM_WORLD,3,MPIgrid,CART_BCs, &
                       reorder,COMM_CART,ierror)
!.Find my coordinate parameters in the Cartesial topology
!.and use them to compute the arrays x,y
  call MPI_CART_COORDS(COMM_CART,myrank,3,blockID,ierror)


!................X, Y and Z SUBCOMMUNICATORS...................!
!.Create x subcommunicators
  remain(1)=.FALSE.     ! Don't keep 1st (z) direction
  remain(2)=.FALSE.     ! Don't keep 2nd (y) direction
  remain(3)=.TRUE.      ! Keep 3rd (x direction)
  call MPI_Cart_sub(COMM_CART, remain, COMM_X, ierror)
  call MPI_Comm_rank(COMM_X, xID, ierr)
  call MPI_Cart_coords(COMM_X, xID, 1, xcoord, ierror)
!.For communication along x, save source and destination rank
  call MPI_CART_SHIFT(COMM_X,0,1,left,right,ierror)
!.MF: These two lines are needed on NERSC. Not sure why.
  if (xID==0) left=MPI_PROC_NULL
  if (xID==xprocsm1) right=MPI_PROC_NULL

!.Create y subcommunicators
  remain(1)=.FALSE.     ! Don't keep 1st (z) direction
  remain(2)=.TRUE.      ! Keep 2nd (y direction)
  remain(3)=.FALSE.     ! Don't keep 3rd (x) direction
  call MPI_Cart_sub(COMM_CART, remain, COMM_Y, ierror)
  call MPI_Comm_rank(COMM_Y, yID, ierr)
  call MPI_Cart_coords(COMM_Y, yID, 1, ycoord, ierror)
!.For communication along y, save source and destination rank
  call MPI_CART_SHIFT(COMM_Y,0,1,down,up,ierror)

!.The following is handy incase myrank.ne.CARTrank
  IDs=(/zID,yID,xID/)
  call MPI_CART_RANK(COMM_CART,IDs,CARTrank,ierror)
!............................................................!
!............XY SUBCOMMUNICATOR FOR 2D xy FFT................!
  remain(1)=.FALSE.      ! Don't keep first (z) direction
  remain(2)=.TRUE.       ! Keep 2nd (y) direction
  remain(3)=.TRUE.       ! Keep 3rd (x) direction
  call MPI_Cart_sub(COMM_CART, remain, COMM_XY, ierror)
  call MPI_Comm_rank(COMM_XY, xyID, ierr)
  call MPI_Cart_coords(COMM_XY, xyID, 2, xycoord, ierror)
!............................................................!

!-Local block's x and z coordinates
  forall(i=1:nx0) xblock(i)=x(i+xID*nx0)
  forall(j=1:ny0) yblock(j)=y(j+yID*ny0)
  forall(k=1:nz0) zblock(k)=z(k+zID*nz0)

!.Test function for Poisson equation with spatially dependent coefficient
  kx0=0.75_dp*TWOPI/Lx
  ky0=2.0_dp*TWOPI/Ly
  lambda0=1.0_dp
  kp=2.0_dp
!  x0=0.2_dp  ! try 0.2 for ny0=1024, Lx=2pi/16
  x0=PI  ! try 0.2 for ny0=1024, Lx=2pi/16
!  sigs=(0.05_dp)**2  ! try 0.5 for ny0=1024, Lx=2pi/16
  sigs=(0.8_dp)**2  ! try 0.5 for ny0=1024, Lx=2pi/16
  allocate(lam(nx0,ny0,nz0),phi0(nx0,ny0,nz0))
  do i=1,nx0
    do j=1,ny0
      lam(i,j,:)=lambda0*(1+((sin(TWOPI*kp*y(j)/Ly))**2)* &
                 exp(-((x(i)-x0)**2)/(2.0_dp*sigs)))
!      lam(i,j,:)=lambda0*(1+exp(-((x(i)-x0)**2)/(2.0_dp*sigs)))
!      lam(i,j,:)=lambda0*(1+exp(-((x(i)-x0)**2)/(2.0_dp*sigs))*exp(-((y(j)-x0)**2)/(2.0_dp*sigs)))
      phi(i,j,:)=cos(kx0*x(i))*sin(ky0*y(j))
      w(i,j,:)=((lam(i,j,:)-lambda0)*(kx0*(x(i)-x0)*tan(kx0*x(i))/sigs &
        +ky0*2.0_dp*TWOPI*kp*cotan(TWOPI*kp*y(j)/Ly)*cotan(ky0*y(j))/Ly) &
        -lam(i,j,:)*(kx0**2+ky0**2))*phi(i,j,:)
!      w(i,j,:)=((lam(i,j,:)-lambda0)*(kx0*(x(i)-x0)*tan(kx0*x(i))/sigs &
!        +ky0*(y(j)-x0)*cotan(ky0*y(j))/sigs) &
!        -lam(i,j,:)*(kx0**2+ky0**2))*phi(i,j,:)
    enddo
  enddo
!  w=-(kx0**2+ky0**2)*phi
!  lam=1.0_dp
  phi0=0.0_dp
!.Right boundary Dirichlet value
  phir=0.0_dp
!  phir=alphad*lmb*0.2_dp

!.Prepare test functions for hyperdiffusion multigrid
  allocate(Diffs(6,2))
!.Rescale and store diffusion coefficients
  kxmax=(2.0d0/3.0d0)*dble(nx0G/2)*twopi/Lx ! using 2/3 de-aliasing rule
  kymax=(2.0d0/3.0d0)*dble(ny0G/2)*twopi/Ly
#IF (HDop==8 .OR. HDop==12)
  if (myrank.eq.0) print*,'x Diffs normalized by kxmax=',kxmax
  if (myrank.eq.0) print*,'y Diffs normalized by kymax=',kymax
  Diffs(1,:)=dt*(/ diffnx/(kxmax**2)   , diffny/(kymax**2)   /)**(HDop/4)
  Diffs(2,:)=dt*(/ diffTex/(kxmax**2)  , diffTey/(kymax**2)  /)**(HDop/4)
  Diffs(3,:)=dt*(/ diffTix/(kxmax**2)  , diffTiy/(kymax**2)  /)**(HDop/4)
  Diffs(4,:)=dt*(/ diffvx/(kxmax**2)   , diffvy/(kymax**2)   /)**(HDop/4)
  Diffs(5,:)=dt*(/ diffwx/(kxmax**2)   , diffwy/(kymax**2)   /)**(HDop/4)
  Diffs(6,:)=dt*(/ diffpsix/(kxmax**2) , diffpsiy/(kymax**2) /)**(HDop/4)
!.Another option is like
!  Diffs(5,:)=dt*(/ diffwx/(kxmax**(HDop/2))   , diffwy/(kymax**(HDop/2))   /) 
#ELSE
  diffexp=HDop !.must be either 2 or 3
  if (HDop==4 .OR. HDop==6) diffexp=HDop/2
!.Two options here:
!.a) Normalize by the smaller k
!  if (kxmax .le. kymax) then
!    knorm = kxmax*kxmax
!  else
!    knorm = kymax*kymax
!  endif
!.b) Normalize by kmaxsq = kxmax**2 + kymax**2
  knorm = kxmax*kxmax + kymax*kymax

  Diffs(1,:)=dt*(/ diffnx/knorm    , diffny/knorm    /)**diffexp
  Diffs(2,:)=dt*(/ diffTex/knorm   , diffTey/knorm   /)**diffexp
  Diffs(3,:)=dt*(/ diffTix/knorm   , diffTiy/knorm   /)**diffexp
  Diffs(4,:)=dt*(/ diffvx/knorm    , diffvy/knorm    /)**diffexp
  Diffs(5,:)=dt*(/ diffwx/knorm    , diffwy/knorm    /)**diffexp
  Diffs(6,:)=dt*(/ diffpsix/knorm  , diffpsiy/knorm  /)**diffexp

!.Another option is like
!  if (kxmax .le. kymax) then
!    knorm=kxmax**2
!  else
!    knorm=kymax**2
!  endif
!  Diffs(5,:)=dt*(/ diffwx/(knorm**diffexp)    , diffwy/(knorm**diffexp)    /)
#ENDIF

!  allocate(phir(ny0,nz0))
!  phir=dble(0.0)
!
  allocate(res(nx0,ny0,nz0))
!.RHS for w hyperdiff
  allocate(wrhs(nx0,ny0,nz0))
  wrhs=w

  allocate(wl(nx0,ny0,nz0,noit),wrhsl(nx0,ny0,nz0,noit))
  do i=1,noit
    wl(:,:,:,i) = w
    wrhsl(:,:,:,i) = wrhs
  enddo

!.Initialize parallel multigrid solvers
  call initialize_PMG(WCOMM,COMM_XY,COMM_X,COMM_Y,Diffs)

  call CPU_TIME(mgt1)
#IF (UTIMERS == 1)
  call CPU_TIME(mgc1)
#ENDIF

  do i=1,noit

  call mgphi(phi0,lam,w,(/2,-1/),phir)
!  call fmgphi(phi0,lam,w,(/2,-1/),phir)

!  call mgpsi(psi0,rpsii,de2,(/-1,-1/))

!!  call mghd(w,wrhs,hdcos(1,:),(/2,2/),1,res)
!  call MGhd(w,wrhs,Diffs(5,:),(/-1,-1/),5)
!  call MGhd(wl(:,:,:,i),wrhsl(:,:,:,i),Diffs(5,:),(/-1,-1/),5,res)

!  deni=dble(1)
!  call fmgphi(psi0,deni,curi,(/-1,-1/),phir)
!  call mgphi(psi0,deni,curi,(/-1,-1/),phir)

  enddo

#IF (UTIMERS == 1)
  call CPU_TIME(mgc2)
  mgloop=mgloop+mgc2-mgc1
#ENDIF
  call CPU_TIME(mgt2)
  if (myrank.eq.0) print*,'time=',(mgt2-mgt1)/dble(noit)

  call MPI_BARRIER(MPI_COMM_WORLD,ierr)

  deallocate(Diffs)
  deallocate(wl,wrhsl)
  call terminate_PMG
  deallocate(lam)
  call dealloc_fields

99 continue
  call dealloc_pars
!.Finilize MPI
  call MPI_FINALIZE(irc)

end program pmgdiff
