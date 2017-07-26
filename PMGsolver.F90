!********************************************************************
!                      PMGsolver
! A parallel multigrid solver for the 2D linear problems
! arising in the plasma fluid model GDB (see B. Zhu, et. al.,
! GDB: a global 3D two-fluid code for plasma turbulence and
! transport in the tokamak edge region, J. of Comp. Phys 2017).
!
! These problems are:
! 1. A Poisson solver to obtain the electrostatic potential phi
!    from the vorticity w
!                 div(lambda*grad(u))=rho
!    where lambda=plasma density, u=phi and rho=w.
!    In the routines below we use the suffix 'ph' to indicate
!    routines to obtain phi.
! 2. An inhomogeneous modified Helmholtz equation to obtain the
!    magnetic flux function psi from the variable rpsi which
!    contains an electron inertia contribution
!                   u-ksq*L(u)=rho
!    where u=psi, rho=rpsi and ksq=(electron skin depth)^2
!    In the routines below we use the suffix 'ps' to indicate
!    routines to obtain psi
! 3. Higher order variants of the (modified) Helmholtz equation
!    arising from implied hyperdiffusion of order 2m
!               u+((-1)^(m+1))*Diff*(L^2m)(u)=rho
!    where Diff is the diffusion coefficient and L is the 2D
!    Laplacian operator, although a non-Laplacian-based diffusion
!    operator is also available respectively. For more details
!    see M. Francisquez, et. al., Multigrid treatment of implicit
!    continuum diffusion, J. of Comp. Phys 2017)
!
! Though these are 2D solvers, we allow for a 3rd dimension so that
! PMGsolver effetively supports solving many xy planes of
! independent 2D problems.
!
! Manaure Francisquez
! mana@dartmouth.edu
! 2017
!
!********************************************************************
! IMPORTANT NOTES
!   i) Not all solvers support the same boundary conditions
!  ii) Modified Helmholtz relaxation has damped relaxation commented out
! ASSUMPTIONS
!   a) Periodic along y and not periodic along x 
!   b) The length of the simulation domain is defined as 2*pi times
!      a factor specified in the input file (fLx along x)
!   c) The problem, the grids and the MPI decomposition are such that
!      only two subdomains along x and/or y can be unified in a single
!      coarsening.
!
! Poisson
!   1. Each process has at least nxmin and nymin points
!      along x and y, respectively.
!   2. nxmin,nymin >= 4
!   3. smallest grid has >=4 points along x
!   4. Boundary condition on the right (along x) is Dirichlet(? check)
!   5. The boundary conditions of lambda are assumed even
! 
! Hyperdiffusion
!   1. Hyperdiffusion may be applied to more than one quantity
!      (e.g. density, temperature, etc.), and a different multigrid
!      parameter set is given for each. nqhd = number of quantities
!      hyperdiffusion is applied to. This variable is inferred from
!      the size of iDiff, and the code looks for nqhd parameter sets
!      in the input file.
!   2. HDop: type of diffusion, some options are (L=Laplacian)
!        =2  u + Diff*(L^2)u=rho
!        =3  u - Diff*(L^3)u=rho
!        =4  u + Diff*(L^2)u=rho split into 2 eqns
!        =6  u - Diff*(L^3)u=rho split into 3 eqns
!        =8  u + Diffx*(d^4 u/dx^4) + Diffy*(d^4 u/dy^4)=rho
!        =12 u - Diffx*(d^6 u/dx^6) - Diffy*(d^6 u/dy^6)=rho
!********************************************************************
module PMGsolver
  implicit none
#include "gdbFLAGS.h"
!.Symbolic names for kind types of 4-byte integers and double precision reals:
  integer, parameter :: I4B = SELECTED_INT_KIND(9)
  integer, parameter :: DP = KIND(1.0D0)

!.Maximum number of grids to visit (=0 if no max)
  integer(I4B) :: LGMAXph,LGMAXps
  integer(I4B), allocatable :: LGMAXhd(:)
!.Semi-coarsening strategy: =2 doubling =4 quadrupling
  integer(I4B) :: NISph,NISps
  integer(I4B), allocatable :: NIShd(:)
!.Type of cycle: GAMMA=1 V-cycle, GAMMA=2 W-cycle
  integer(I4B) :: GAMph,GAMps
  integer(I4B), allocatable :: GAMhd(:)
!.Number of (GAMMA-)cycles (at every grid if using FMG)
  integer(I4B) :: BETph,BETps
  integer(I4B), allocatable :: BEThd(:)
!.Number of relaxation sweeps before and after the course grid correction
  integer(I4B) :: NU1ph,NU2ph,NU1ps,NU2ps
  integer(I4B), allocatable :: NU1hd(:),NU2hd(:)
!.Damping coefficients in relaxation
  real(DP) :: omeph,Oomeph,omeps,Oomeps
  real(DP), allocatable :: omehd(:),Oomehd(:)
  integer(I4B), parameter :: HDopd2=HDop/2,HDopd4=HDop/4
!.Number of elements needed from diagonal neighbors
!.(only needed when HDop=2 or 3)
  integer(I4B), parameter :: ndiag=1+(HDop-2)*2
!.Number of quantities hyperdiffusion is applied to
  integer(I4B) :: nqhd

!.Define a type for an array 3D real arrays
  type allo3D
    real(DP), allocatable :: a(:,:,:)
  end type allo3D
!.Define a type for an array 2D real arrays
  type allo2D
    real(DP), allocatable :: a(:,:)
  end type allo2D
!.Define a type for an array 1D integer arrays
  type Iallo1D
    integer(I4B), allocatable :: a(:)
  end type Iallo1D
!.Define a type for an array 2D integer arrays
  type Iallo2D
    integer(I4B), allocatable :: a(:,:)
  end type Iallo2D
!.Define a type for an array 3D integer arrays
  type Iallo3D
    integer(I4B), allocatable :: a(:,:,:)
  end type Iallo3D
!.Define a type for an array 2D boolean arrays
  type Lallo2D
    logical, allocatable :: a(:,:)
  end type Lallo2D
!.Define a type for an array of arrays of 1D integer arrays
  type IIallo2D
    type(Iallo2D), allocatable :: a(:)
  end type IIallo2D
!.For Poisson equation, one needs lambda coarsened to each
!.grid level. That is a 3D array at each grid level
  type(allo3D), allocatable :: Diff(:)

!.Local (to this rank) fine and coarse dimensions of grids visited
  integer(I4B), allocatable, dimension(:) :: nxFph,nyFph,nxFps,nyFps
  type(Iallo1D), allocatable, dimension(:) :: nxFhd,nyFhd
  integer(I4B), allocatable, dimension(:) :: nxCph,nyCph,nxCps,nyCps
  type(Iallo1D), allocatable, dimension(:) :: nxChd,nyChd
!.Indicate which dimension was coarsened, =0 for both, =1 for x, =2 for y
  integer(I4B), allocatable, dimension(:) :: wdcph,wdcps
  type(Iallo1D), allocatable, dimension(:) :: wdchd
!.Number of grids to be visited
  integer(I4B) :: lgph,lgps
  integer(I4B), allocatable :: lghd(:)
!.Spatial grid details
  real(DP) :: bLx,bLy                   !.Length of the domain box along x and y
  integer(I4B) :: nxG,nyG               !.Number of points in global grid
  integer(I4B) :: nxL,nyL,nzL           !.Number of points in local grid
!.Constant hyper-diffusion coefficient
  real(DP), allocatable :: hDiff(:)
!.k^2 factor in Helmholtz equation = de^2 (electron skin depth squared) in GDB
  real(DP) :: ksq
!.Resolution of current grids, fine and coarse
  integer(I4B) :: nxc,nyc,nxf,nyf
!.Local grid dimensions after unification of subdomains takes place
  integer(I4B) :: nxu,nyu

!........................................................!
!................  MPI related variables ................!
  integer(I4B) :: wrnk0             !.initial rank number in MPI_WORLD_COMM
  integer(I4B) :: iprocs0,jprocs0   !.Number of (finest grid) processes in x and y
  integer(I4B) :: iID0,jID0         !.This process' initial ID in x and y COMMs
  integer(I4B) :: myrank0           !.This process' initial ID in xy COMM
  integer(I4B) :: Wcomm0,XYcomm0    !.Initial World and XY COMMs
  logical :: Lhalf0,Bhalf0          !.Is this rank in Left 1/2 (along x) and/or Bottom 1/2 (along y)
  integer(I4B) :: nxmin,nymin       !.min number of cells per process

!.In the following arrays, entries correspond to differend grid levels
!.Note: x (East-West) and y (North-South) rank

!.rank ID in X and Y COMMs
  integer(I4B), allocatable, dimension(:) :: iIDph,jIDph,iIDps,jIDps
  type(Iallo1D), allocatable, dimension(:) :: iIDhd,jIDhd
!.Number of processors along each direction 
  integer(I4B), allocatable, dimension(:) :: iprocsph,jprocsph,iprocsps, &
                                             jprocsps
  type(Iallo1D), allocatable, dimension(:) :: iprocshd,jprocshd
!.Number of processors along each direction (minus 1)
  integer(I4B) :: iprocsm1,jprocsm1
!.Rank of neighbors, starting with the one above (along y/North) and 
!.proceeding clockwise: 1-North, 2-NE, 3-E, 4-SE, ..., 8-NW
  integer(I4B), allocatable, dimension(:,:) :: nborph,nborps
  type(Iallo2D), allocatable, dimension(:) :: nborhd
!.Array of communicators of procs still active at grid level ng.
!.COMM(ng,1) is the y (North-South) communicator, COMM(ng,2) is the x
!.(East-West) communicator and COMM(ng,3) the 2D communicator for
!.diagonal communication 
  integer(I4B), allocatable, dimension(:,:) :: COMMph,COMMps
  type(Iallo2D), allocatable, dimension(:) :: COMMhd
!.This array contains a boolean for indicating whether this rank is active
!.and whether unification of subdomains takes place or not
!.acu(ng,1) TRUE if this process is still active
!.acu(ng,2) TRUE if it is time to unify subdomains
  logical, allocatable, dimension(:,:) :: acuph,acups
  type(Lallo2D), allocatable, dimension(:) :: acuhd
!.Ranks data is sent to/received from when subdomains are unified
  type(Iallo2D), allocatable :: urnkph(:),urnkps(:)
  type(IIallo2D), allocatable :: urnkhd(:)
!.Number of ranks to whose subdomains will be unified along x and y
  integer(I4B), allocatable :: nurnkph(:,:),nurnkps(:,:)
  type(Iallo2D), allocatable :: nurnkhd(:)

#IF (MGTIMER==1)
!.These are user-defined timers for timing different parts of the code
!.used in profiling
  real(DP) :: initPMG,init2,init1,mg1,mg2,rstdiff1,rstdiff2, &
              prolt1,prolt2,count1,count2,rest1,rest2
  real(DP) :: mgtph,resDtph,postXtph,preXtph,proltph,rNrtph,restph
  real(DP) :: mgtps,resDtps,postXtps,preXtps,proltps,rNrtps,restps
  real(DP) :: mgthd,resDthd,postXthd,preXthd,prolthd,rNrthd,resthd
#ENDIF

contains
  subroutine initialize_PMG(iWcomm,iXYcomm,iXcomm,iYcomm,iDiff)
!.Initialize the PMG solver arrays and communicators needed
!.INPUTS
!.  integer iWcomm: WORLD communicator
!.  integer iXYcomm: 2D XY communicator 
!.  integer iXcomm: 1D X communicator 
!.  integer iYcomm: 1D Y communicator 
!.  real(dp) dimension(:,2) iDiff: Diffusion coefficients along x
!.                                 and y of nqhd quantities
  implicit none
  include 'mpif.h'
  integer(I4B), intent(in) :: iWcomm,iXYcomm,iXcomm,iYcomm
  real(DP), intent(in) :: iDiff(:,:)
  integer(I4B) :: ng,i,j,k,ierr,jnk1(1),jnk2(1),IDs1D(1)
  logical :: jnkl(1)

#IF (MGTIMER==1)
  initPMG=0.0_dp
  call CPU_TIME(init1)
#ENDIF

!.Use a duplicate of the input world COMM for safety/simplicity
  call MPI_COMM_DUP(iWcomm,Wcomm0,ierr) ! World communicator
  call MPI_COMM_RANK(Wcomm0,wrnk0,ierr) ! rank in world comm

  nqhd = size(iDiff,1) !.Number of quantities hyperdiffused

  call read_mginputs !.Read input multigrid parameter set(s)

!.Use a duplicate of the input 2D COMM for safety/simplicity
!.The 2D COMM is used in unifying subdomains, prolongation
!.and diagonal communication
  call MPI_COMM_DUP(iXYcomm,XYcomm0,ierr)

!.Initial grid rank ID's in each COMM
  call MPI_Comm_rank(iXcomm,iID0,ierr)
  call MPI_Comm_rank(iYcomm,jID0,ierr)
  call MPI_CART_RANK(iXYcomm,(/jID0,iID0/),myrank0,ierr)
!.Is this rank in the Left Half of the domain (along x)?
  Lhalf0=.TRUE.
  if (iID0>=iprocs0/2) Lhalf0=.FALSE. 
!.Is this rank in the Bottom Half of the domain (along y)?
  Bhalf0=.TRUE.
  if (jID0>=jprocs0/2) Bhalf0=.FALSE. 

!.Local number of grid points in each direction
  nxL=nxG/iprocs0
  nyL=nyG/jprocs0

!.Determine the grids to be visited by the Poisson and
!.Helmholtz solvers, when to unify subdomains, and initialize necessary arrays
  call pgrids(1,LGMAXph,NISph,lgph,nxFph,nyFph,nxCph,nyCph,wdcph, &
              acuph,iprocsph,jprocsph,urnkph,nurnkph)
  call pgrids(2,LGMAXps,NISps,lgps,nxFps,nyFps,nxCps,nyCps,wdcps, &
              acups,iprocsps,jprocsps,urnkps,nurnkps)
!.Need to allocate space for coarsened lambda at each grid (Diff)
  allocate(Diff(lgph))
  do j=lgph,1,-1
    allocate(Diff(j)%a(nxFph(j),nyFph(j),nzL))
  enddo

!.Communicators at each grid level. The 2D communicator (COMM(ng,3))
!.is used to deal with corner grid points in stencils.
  allocate(COMMph(lgph,3),COMMps(lgps,3))
!.Copy initial (finest) 2D communicator into array of communicators
  call MPI_COMM_DUP(XYcomm0,COMMph(lgph,3),ierr)
  call MPI_COMM_DUP(XYcomm0,COMMps(lgps,3),ierr)
!.Create 2D xy communicators for each coarse grid
  call COMMv2D(COMMph(:,3),XYcomm0,acuph(:,1),iprocsph,jprocsph,lgph)
  call COMMv2D(COMMps(:,3),XYcomm0,acups(:,1),iprocsps,jprocsps,lgps)
!.Save the rank ID of the diagonal neighbors (corners) NE,SE,SW,NW
  allocate(nborph(lgph,8),nborps(lgps,8))
  call diagNeigh(nborph,COMMph(:,3),iprocsph,jprocsph,acuph(:,1),lgph)
  call diagNeigh(nborps,COMMps(:,3),iprocsps,jprocsps,acups(:,1),lgps)

!.Copy initial x and y communicators into array of communicators
  call MPI_COMM_DUP(iYcomm,COMMph(lgph,1),ierr)
  call MPI_COMM_DUP(iXcomm,COMMph(lgph,2),ierr)
  call MPI_COMM_DUP(iYcomm,COMMps(lgps,1),ierr)
  call MPI_COMM_DUP(iXcomm,COMMps(lgps,2),ierr)
!.Create x and y communicators for coarse grids
  call COMMv1D(COMMph(:,1),iYcomm,acuph(:,1),jprocsph,lgph,(/.TRUE./),jID0)
  call COMMv1D(COMMph(:,2),iXcomm,acuph(:,1),iprocsph,lgph,(/.FALSE./),iID0)
  call COMMv1D(COMMps(:,1),iYcomm,acups(:,1),jprocsps,lgps,(/.TRUE./),jID0)
  call COMMv1D(COMMps(:,2),iXcomm,acups(:,1),iprocsps,lgps,(/.FALSE./),iID0)

  allocate(iIDph(lgph),iIDps(lgps))
  allocate(jIDph(lgph),jIDps(lgps))
  do ng=lgph,1,-1         !.Loop over grids visited
    if (acuph(ng,1)) then !.Check if this rank is active (not idle)
!.....This rank's ID in x (East-West) direction
      call MPI_Comm_rank(COMMph(ng,2),iIDph(ng),ierr)
!.....For communication along x, save East and West neighbor rank
      call MPI_CART_SHIFT(COMMph(ng,2),0,1,nborph(ng,7),nborph(ng,3),ierr)
!.....MF: These are needed on NERSC. Not sure why.
      if (iIDph(ng)==0) nborph(ng,7)=MPI_PROC_NULL
      if (iIDph(ng)==iprocsph(lgph)) nborph(ng,3)=MPI_PROC_NULL
!.....This rank's ID in y (North-South) direction
      call MPI_Comm_rank(COMMph(ng,1),jIDph(ng), ierr)
!.....For communication along y, save North and South neighbor rank
      call MPI_CART_SHIFT(COMMph(ng,1),0,1,nborph(ng,5),nborph(ng,1),ierr)
    endif
  enddo
  do ng=lgps,1,-1         !.Loop over grids visited
    if (acups(ng,1)) then !.Check if this rank is active (not idle)
!.....This rank's ID in x (East-West) direction
      call MPI_Comm_rank(COMMps(ng,2),iIDps(ng), ierr)
!.....For communication along x, save East and West neighbor rank
      call MPI_CART_SHIFT(COMMps(ng,2),0,1,nborps(ng,7),nborps(ng,3),ierr)
!.....MF: These are needed on NERSC. Not sure why.
     if (iIDps(ng)==0) nborps(ng,7)=MPI_PROC_NULL
     if (iIDps(ng)==iprocsps(lgps)) nborps(ng,3)=MPI_PROC_NULL
!.....This rank's ID in y (North-South) direction
      call MPI_Comm_rank(COMMps(ng,1),jIDps(ng), ierr)
!.....For communication along y, save North and South neighbor rank
      call MPI_CART_SHIFT(COMMps(ng,1),0,1,nborps(ng,5),nborps(ng,1),ierr)
    endif
  enddo

!.Now do the same for hyperdiffusion quantities
  allocate(hDiff(2)) ! x-y diffusion coefficients used when solver is called
!.Check nxmin and nymin are greater than 8
  if (nxmin<8) then
    nxmin=8
    if (wrnk0==0) print*,'Changing to nxmin=8 for PMG hyperdiff'
  endif
  if (nymin<8) then
    nymin=8
    if (wrnk0==0) print*,'Changing to nymin=8 for PMG hyperdiff'
  endif
!.Allocate arrays for hyperdiffusion parameters of various quantities
  allocate(lghd(nqhd),nxFhd(nqhd),nyFhd(nqhd),nxChd(nqhd),nyChd(nqhd), &
           wdchd(nqhd),acuhd(nqhd),iprocshd(nqhd),jprocshd(nqhd), &
           urnkhd(nqhd),nurnkhd(nqhd))
  allocate(COMMhd(nqhd),nborhd(nqhd),iIDhd(nqhd),jIDhd(nqhd))
  do j=1,nqhd   !.Loop over quantities hyperdiffused
!...Determine the grids to be visited, when to unify subdomains,
!...and initialize necessary arrays
    call pgrids(3,LGMAXhd(j),NIShd(j),lghd(j),nxFhd(j)%a,nyFhd(j)%a, &
      nxChd(j)%a,nyChd(j)%a,wdchd(j)%a,acuhd(j)%a,iprocshd(j)%a, &
      jprocshd(j)%a,urnkhd(j)%a,nurnkhd(j)%a,iDiff(j,:))
!...Communicators at each grid level. The 2D communicator (COMM(ng,3))
!...is used to deal with corner grid points in stencils.
    allocate(COMMhd(j)%a(lghd(j),3))
!...Copy 2D communicator into array of communicators
    call MPI_COMM_DUP(XYcomm0,COMMhd(j)%a(lghd(j),3),ierr)
!...Create 2D xy communicators for each coarse grid
    call COMMv2D(COMMhd(j)%a(:,3),XYcomm0,acuhd(j)%a(:,1), &
                 iprocshd(j)%a,jprocshd(j)%a,lghd(j))
!.Save the rank ID of the diagonal neighbors (corners) NE,SE,SW,NW
    allocate(nborhd(j)%a(lghd(j),8))
    call diagNeigh(nborhd(j)%a,COMMhd(j)%a(:,3),iprocshd(j)%a, &
                   jprocshd(j)%a,acuhd(j)%a(:,1),lghd(j))
!...Copy initial x and y communicators into array of communicators
    call MPI_COMM_DUP(iYcomm,COMMhd(j)%a(lghd(j),1),ierr)
    call MPI_COMM_DUP(iXcomm,COMMhd(j)%a(lghd(j),2),ierr)
!...Create x and y communicators for coarse grids
    call COMMv1D(COMMhd(j)%a(:,1),iYcomm,acuhd(j)%a(:,1), &
                 jprocshd(j)%a,lghd(j),(/.TRUE./),jID0)
    call COMMv1D(COMMhd(j)%a(:,2),iXcomm,acuhd(j)%a(:,1), &
                 iprocshd(j)%a,lghd(j),(/.FALSE./),iID0)
    allocate(iIDhd(j)%a(lghd(j)),jIDhd(j)%a(lghd(j)))
    do ng=lghd(j),1,-1           !.Loop over grids visited
      if (acuhd(j)%a(ng,1)) then !.Check if this rank is active (not idle)
!.......This rank's ID in x (East-West) direction
        call MPI_Comm_rank(COMMhd(j)%a(ng,2),iIDhd(j)%a(ng),ierr)
!.......For communication along x, save East and West neighbor rank
        call MPI_CART_SHIFT(COMMhd(j)%a(ng,2),0,1, &
                            nborhd(j)%a(ng,7),nborhd(j)%a(ng,3),ierr)
!.......MF: These are needed on NERSC. Not sure why.
        if (iIDhd(j)%a(ng)==0) nborhd(j)%a(ng,7)=MPI_PROC_NULL
        if (iIDhd(j)%a(ng)==iprocshd(j)%a(lghd(j))) nborhd(j)%a(ng,3)=MPI_PROC_NULL
!.......This rank's ID in y (North-South) direction
        call MPI_Comm_rank(COMMhd(j)%a(ng,1),jIDhd(j)%a(ng), ierr)
!.......For communication along y, save North and South neighbor rank
        call MPI_CART_SHIFT(COMMhd(j)%a(ng,1),0,1, &
                            nborhd(j)%a(ng,5),nborhd(j)%a(ng,1),ierr)
      endif
    enddo
  enddo

#IF (MGTIMER==1)
  call CPU_TIME(init2)
  initPMG=init2-init1
#ENDIF

  end subroutine initialize_PMG
!*********************************************************
  subroutine read_mginputs
!.Read the input resolution, MPI decomposition, factors defining
!.the domain lengh along x and y, and the multigrid parameter set.
!.The multigrid parameters control which grids are visited, the
!.type and number of iterations to be performed, and when to
!.unify subdomains.
  implicit none
  include 'mpif.h'
  integer(I4B) :: nx0G,ny0G,nz0G,xprocs,yprocs,zprocs
  real(DP) :: fLx,fLy,fLz
  integer(I4B) :: ierr

  allocate(LGMAXhd(nqhd),NIShd(nqhd),GAMhd(nqhd),BEThd(nqhd), &
           NU1hd(nqhd),NU2hd(nqhd),omehd(nqhd),Oomehd(nqhd))
!.Add here the namelists for parameters (first declared above) from input file:
  namelist/dom/nx0G,ny0G,nz0G,xprocs,yprocs,zprocs,fLx,fLy,fLz
  namelist/mgin/LGMAXph,NISph,GAMph,BETph,NU1ph,NU2ph,omeph, &
                LGMAXps,NISps,GAMps,BETps,NU1ps,NU2ps,omeps, &
                LGMAXhd,NIShd,GAMhd,BEThd,NU1hd,NU2hd,omehd, &
                nxmin,nymin

  if (wrnk0 .eq. 0) then
    open(unit=10,file="gdb.in",status="old")
    read (10,dom)
    read (10,mgin)
    close(unit=10)

    nxG=nx0G
    nyG=ny0G
    iprocs0=xprocs
    jprocs0=yprocs
    nzL=nz0G/zprocs

    bLx=fLx*8.0_dp*atan(1.0_8) !.fLx*2*pi
    bLy=fLy*8.0_dp*atan(1.0_8) !.fLy*2*pi
  endif

!.Since only wrnk0 read the inputs, broadcast them to World COMM
  call MPI_BCAST(nxG,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(nyG,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(nzL,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(iprocs0,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(jprocs0,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(bLx,1,MPI_DOUBLE_PRECISION,0,Wcomm0,ierr)
  call MPI_BCAST(bLy,1,MPI_DOUBLE_PRECISION,0,Wcomm0,ierr)
  call MPI_BCAST(LGMAXph,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(NISph,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(GAMph,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(BETph,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(NU1ph,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(NU2ph,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(omeph,1,MPI_DOUBLE_PRECISION,0,Wcomm0,ierr)
  call MPI_BCAST(LGMAXps,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(NISps,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(GAMps,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(BETps,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(NU1ps,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(NU2ps,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(omeps,1,MPI_DOUBLE_PRECISION,0,Wcomm0,ierr)
  call MPI_BCAST(LGMAXhd,nqhd,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(NIShd,nqhd,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(GAMhd,nqhd,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(BEThd,nqhd,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(NU1hd,nqhd,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(NU2hd,nqhd,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(omehd,nqhd,MPI_DOUBLE_PRECISION,0,Wcomm0,ierr)
  call MPI_BCAST(nxmin,1,MPI_INTEGER,0,Wcomm0,ierr)
  call MPI_BCAST(nymin,1,MPI_INTEGER,0,Wcomm0,ierr)

!.One minus the damping parameter in relaxations. Useful in relax routines
  Oomeph=1.0_dp-omeph
  Oomeps=1.0_dp-omeps
  Oomehd=1.0_dp-omehd

  end subroutine read_mginputs
!*********************************************************
  subroutine COMMv1D(COMMv,iCOMM,act,procs,lg,BCs,ID0)
!.Create communicators (along x or y) for each (coarse) level
!.visited and where this process is active. These communicators
!.are used to exchange information in the relaxation, residue,
!.prolongation and restriction routines.
!.INPUTS
!.  integer iCOMM: initial (finest) 1D communicator
!.  logical dimension(:) act: active status of this rank at each grid
!.  integer dimension(:) procs: number of active processes at each grid
!.  integer lg: number of grids to visit
!.  integer BCs: boundary conditions, TRUE for periodic
!.  integer ID0: rank ID in iCOMM
!.OUTPUT
!.  integer dimension(:) COMMv: array of COMMs of active processes in each grid
  implicit none
  include 'mpif.h'
  integer(I4B), intent(in) :: iCOMM,lg,ID0
  integer(I4B), intent(in), allocatable :: procs(:)
  integer(I4B), intent(inout) :: COMMv(:)
  logical, intent(in) :: BCs(:)
  logical, intent(in) :: act(:)
  integer(I4B) :: ng,i,j,k,ierr,IDs1D(1),gap,procs0, &
                  iCOMMgroup,tmpgroup,tmpCOMM
  logical :: reorder=.TRUE.
  integer(I4B), allocatable :: ranks(:)

  procs0=procs(lg) !.initial number of active processes

!.Obtain MPI group of original 1D communicator
  call MPI_Comm_group(iCOMM,iCOMMgroup,ierr)
  do ng=lg-1,1,-1              !.Loop over coarse grids
    allocate(ranks(procs(ng))) !.Ranks in this grid
    k=0                        !.Number of ranks in this grid
    if (procs(ng)>1) then
      gap=procs0/procs(ng)     !.gap (+1) between active processes
!.....This does the left or bottom half of the active processes
      do i=gap,procs0/2,gap
        k=k+1
!        IDs1D=i-1
!        call MPI_CART_RANK(iCOMM,IDs1D,ranks(k),ierr)
        ranks(k)=i-1           !.Rank IDs are zero-indexed
      enddo
!.....This does the right or top half of the active processes
      do i=procs0/2+1,procs0,gap
        k=k+1
!        IDs1D=i-1
!        call MPI_CART_RANK(iCOMM,IDs1D,ranks(k),ierr)
        ranks(k)=i-1           !.Rank IDs are zero-indexed
      enddo
    else
!      IDs1D=iID0
!      call MPI_CART_RANK(iXcomm,IDs1D,xranks(1),ierr)
      ranks=ID0                !.Single active rank
    endif
!...Create new group with only active processes
    call MPI_GROUP_incl(iCOMMgroup,procs(ng),ranks,tmpgroup,ierr)
!...Create new communicator for active processes
!...MPI_COMM_create group available in newer compilers only and is
!...collective only for processes in the group
    call MPI_COMM_create_group(iCOMM,tmpgroup,0,tmpCOMM,ierr)
!    call MPI_COMM_create(iCOMM,tmpgroup,tmpCOMM,ierr)
!...Make a cartesian communicator and save in array of communicators
    if (act(ng)) call MPI_CART_create(tmpCOMM,1,(/procs(ng)/),BCs,reorder,COMMv(ng),ierr)
!...Destroy temporary group and communicator to be used in next iteration
    call MPI_GROUP_free(tmpgroup,ierr)
    if (act(ng)) call MPI_COMM_free(tmpCOMM,ierr)
    deallocate(ranks)
  enddo
  end subroutine COMMv1D
!*********************************************************
  subroutine COMMv2D(COMMv,iCOMM,act,iprocs,jprocs,lg)
!.Create 2D communicators for each (coarse) level visited
!.INPUTS
!.  integer iCOMM: initial (finest) 2D communicator
!.  logical dimension(:) act: active status of this rank at each grid
!.  integer dimension(:) iprocs: number of active processes along x at each grid
!.  integer dimension(:) jprocs: number of active processes along y at each grid
!.  integer lg: number of grids to visit
!.OUTPUT
!.  integer dimension(:) COMMv: array of 2D COMMs of active processes in each grid
  implicit none
  include 'mpif.h'
  integer(I4B), intent(in) :: iCOMM,lg
  integer(I4B), intent(in), allocatable :: iprocs(:),jprocs(:)
  integer(I4B), intent(inout) :: COMMv(:)
  logical :: BCs(2)
  logical, intent(in) :: act(:)
  integer(I4B) :: ng,i,j,k,ierr,IDs2D(2),igap,jgap,iprocs0, &
                  jprocs0,iCOMMgroup,tmpgroup,tmpCOMM
  logical :: reorder=.TRUE.
  integer(I4B), allocatable :: iranks(:),jranks(:),xyranks(:)

  iprocs0=iprocs(lg)      !.initial number of active processes
  jprocs0=jprocs(lg)      !.initial number of active processes
  BCs=(/.TRUE.,.FALSE./)  !.Boundary conditions along y and x

!.Obtain MPI group of original 1D communicator
  call MPI_Comm_group(iCOMM,iCOMMgroup,ierr)
  do ng=lg-1,1,-1                !.Loop over coarse grids
    allocate(iranks(iprocs(ng))) !.Ranks along x in this grid
    k=0                          !.Number of ranks in this grid
    if (iprocs(ng)>1) then
      igap=iprocs0/iprocs(ng)    !.gap (+1) between active processors
!.....This does the left half (along x) of the active processes
      do i=igap,iprocs0/2,igap
        k=k+1
        iranks(k)=i-1            !.Rank IDs are zero-indexed
      enddo
!.....This does the right half (along x) of the active processes
      do i=iprocs0/2+1,iprocs0,igap
        k=k+1
        iranks(k)=i-1            !.Rank IDs are zero-indexed
      enddo
    else
      iranks=iID0                !.Single active rank
    endif
    allocate(jranks(jprocs(ng))) !.Ranks along y in this grid
    k=0
    if (jprocs(ng)>1) then
      jgap=jprocs0/jprocs(ng)    !.gap (+1) between active processors
!.....This does the bottom half (along y) of the active processes
      do j=jgap,jprocs0/2,jgap
        k=k+1
        jranks(k)=j-1            !.Rank IDs are zero-indexed
      enddo
!.....This does the top half (along y) of the active processes
      do j=jprocs0/2+1,jprocs0,jgap
        k=k+1
        jranks(k)=j-1            !.Rank IDs are zero-indexed
      enddo
    else
      jranks=jID0                !.Single active rank
    endif
!...Determine rank IDs of all active processes in 2D communicator
    allocate(xyranks(iprocs(ng)*jprocs(ng)))
    k=0
    do j=1,jprocs(ng)
      do i=1,iprocs(ng)
        k=k+1
        IDs2D=(/jranks(j),iranks(i)/)
        call MPI_CART_RANK(iCOMM,IDs2D,xyranks(k),ierr)
      enddo
    enddo
!...Create new group with only active processes
    call MPI_GROUP_incl(iCOMMgroup,iprocs(ng)*jprocs(ng),xyranks,tmpgroup,ierr)
!...Create new communicator for active processes
!....MPI_COMM_create group available in newer compilers only and is
!....collective only for processes in the group
    call MPI_COMM_create_group(iCOMM,tmpgroup,0,tmpCOMM,ierr)
!    call MPI_COMM_create(iCOMM,tmpgroup,tmpCOMM,ierr)
!...Make a cartesian communicator and save in array of communicators
    if (act(ng)) call MPI_CART_create(tmpCOMM,2,(/jprocs(ng),iprocs(ng)/), &
                   BCs,reorder,COMMv(ng),ierr)
!...Destroy temporary group and communicator to be used in next iteration
    call MPI_GROUP_free(tmpgroup,ierr)
    if (act(ng)) call MPI_COMM_free(tmpCOMM,ierr)
    deallocate(iranks,jranks,xyranks)
  enddo
  end subroutine COMMv2D
!*********************************************************
  subroutine diagNeigh(nbors,XYcomm,iprocs,jprocs,act,lg)
!.Determine the diagonal (NE,SE,SW,NW) neighbors at every grid level
!.INPUTS
!.  integer dimension(:) XYcomm: 2D COMM at each grid
!.  integer dimension(:) iprocs: number of active processes along x at each grid
!.  integer dimension(:) jprocs: number of active processes along y at each grid
!.  logical dimension(:) act: active status of this rank at each grid
!.  integer lg: number of grids to visit
!.OUTPUT
!.  integer dimension(:,:) nbors: array of neibor rank IDs at each grid
!.                                only elements 2,4,6,8 are initialized here
  implicit none
  include 'mpif.h'
  integer(I4B), intent(in) :: XYcomm(:),iprocs(:),jprocs(:)
  logical, intent(in) :: act(:)
  integer(I4B), intent(in) :: lg
  integer(I4B), intent(inout) :: nbors(:,:)
  integer(I4B) :: ng,i,ierr,coord(2),NE(2),SE(2),SW(2),NW(2)
  logical :: jnkl(2)

  do ng=lg,1,-1                    !.Loop over grids visited
!...Initiate them all to NULL, so non-periodic boundaries are NULL
    forall(i=2:8:2) nbors(ng,i)=MPI_PROC_NULL 
    if (act(ng)) then
      call MPI_CART_GET(XYcomm(ng),2,(/jprocs(ng),iprocs(ng)/), &
                        jnkl,coord,ierr) !.This ranks coordinates in 2D COMM
      NE=(/coord(1)+1,coord(2)+1/) !.NorthEast neighbor coordinates 
      SE=(/coord(1)-1,coord(2)+1/) !.SouthEast neighbor coordinates 
      SW=(/coord(1)-1,coord(2)-1/) !.SouthWest neighbor coordinates 
      NW=(/coord(1)+1,coord(2)-1/) !.NorthWest neighbor coordinates 
      if (coord(1)==0) then !.Enforce periodicity in bottom boundary
        SE(1)=jprocs(ng)-1
        SW(1)=jprocs(ng)-1
      endif
      if (coord(1)==jprocs(ng)-1) then !.Enforce periodicity in top boundary
        NE(1)=0
        NW(1)=0
      endif
!.....For non-null diagonal neigbors conver coordinates to a rank ID
      if (coord(2)/=0) then
        call MPI_CART_rank(XYcomm(ng),SW,nbors(ng,6),ierr) !.SW
        call MPI_CART_rank(XYcomm(ng),NW,nbors(ng,8),ierr) !.NW
      endif
      if (coord(2)/=iprocs(ng)-1) then
        call MPI_CART_rank(XYcomm(ng),NE,nbors(ng,2),ierr) !.NE
        call MPI_CART_rank(XYcomm(ng),SE,nbors(ng,4),ierr) !.SE
      endif
    endif
  enddo
  end subroutine diagNeigh
!*********************************************************
  subroutine pgrids(eq,lgMAX,NIS,lg,nxv,nyv,Cnx,Cny,wdc,acu,iprocs, &
                    jprocs,urnk,nurnk,iDiff)
!.Determine the grids visited based on coupling coefficients, when processes
!.are active or not, when to unify subdomains and how many processes remain
!.active at each grid level
!.INPUTS
!.  integer eq: equation, =1 Poisson, =2 modified Helmholtz, =3 hyperdiffusion
!.  integer lgMAX: max number of grids
!.  integer NIS: semi-coarsening factor, =2 double spacing, =4 quadruple spacing
!.  real(dp) dimension(:) iDiff: diffusion coefficients along x and y
!.OUTPUTS
!.  integer dimension(:) nxv: Number of points along x of grids visited
!.  integer dimension(:) nyv: Number of points along y of grids visited
!.  integer dimension(:) Cnx: Number of points along x of next coarser grid at each level
!.  integer dimension(:) Cny: Number of points along y of next coarser grid at each level
!.  integer dimension(:) wdc: "which dimension coarsened" last at each level
!.  logical dimension(:,:) acu: is process active and is subdomain unification
!.                              happening at this level?
!.  integer dimension(:) iprocs: number of processes along x at each level
!.  integer dimension(:) jprocs: number of processes along y at each level
!.  integer dimension(:) urnk: ranks data is sent to/received from when unifying
!.                             subdomains
!.  integer dimension(:,:) nurnk: number of ranks data is sent to/received from
!.                                when unifying subdomains

  implicit none
  include 'mpif.h'
  integer(I4B), intent(in) :: eq, lgMAX, NIS
  integer(I4B), intent(out) :: lg
  integer(I4B), intent(out), allocatable :: iprocs(:),jprocs(:)
  integer(I4B), allocatable, dimension(:), intent(out) :: nxv,nyv,Cnx, &
                                                          Cny,wdc
  integer(I4B), allocatable, intent(out) :: nurnk(:,:)
  type(Iallo2D), allocatable, intent(out) :: urnk(:)
  logical, allocatable, intent(out) :: acu(:,:)
  logical, allocatable :: active(:),unify(:)
  real(DP), intent(in), optional :: iDiff(:)
  real(DP), allocatable :: hx(:),hy(:) !.cell lengths
  real(DP) :: hxs,hys   !.cell lengths squared
  real(DP), allocatable :: cf(:,:)      !.coupling coefficients
  integer(I4B) :: ng                    !.current grid
  integer(I4B) :: igap,jgap,igap0,jgap0 !.gap between processes along x and y
  logical :: oneiproc,onejproc  !.true when # of active procs=1
  logical :: iactive,jactive    !.is process active along x/y?
  logical :: iunify,junify      !.does subdomain unification happen along x/y?
  logical :: iunisend,junisend  !.does this process send its data during unification along x/y?
  integer(I4B) :: minCnx,minCny !.minimum number of points in subdomain along x/y
  integer(I4B) :: inurnk,jnurnk !.number of ranks to unify subdomains with along x/y
  integer(I4B) :: i,j,k,ierror
!.temporary arrays
  real(DP), allocatable :: tmp(:,:),tmphx(:),tmphy(:)
  real(DP) :: tmpcfx,tmpcfy
  integer(I4B), allocatable, dimension(:) :: tmpx,tmpy,tmpw,tmpproc,iugap,jugap
  logical, allocatable :: tmpact(:),tmpuni(:)
  type(Iallo2D), allocatable :: tmpurnk(:)
  integer(I4B), allocatable :: tmpnurnk(:,:)

  nxf=nxL  !.Fine local resolution along x
  nyf=nyL  !.Fine local resolution along y
  nxc=nxL  !.Coarse local resolution along x
  nyc=nyL  !.Coarse local resolution along y

  lg=1                  !.Initial number of grids=1
  iactive=.TRUE.        !.All processes active along x
  jactive=.TRUE.        !.All processes active along y
  allocate(hx(lg),hy(lg)) !.Current fine grid spacing
  allocate(active(lg),unify(lg)) 
  active(lg)=.TRUE.     !.All processes active initially
!.No unification happens at first grid. Requires initial MPI decomposition
!.to be sufficiently coarse (see nxmin and nymin in input file)
  unify(lg)=.FALSE.
!.Gap between active processors. igap=1 every processor along x
!.is active, igap=2 every other processor along x is active, etc.
  igap=1
  jgap=1
!.Number of active processors along each direction
  allocate(iprocs(lg),jprocs(lg))
  iprocs(lg)=iprocs0
  jprocs(lg)=jprocs0

  oneiproc=.FALSE.      !.true when # of active procs=1
  onejproc=.FALSE.
!.Where data is sent to/received from when subdomains are unified
  allocate(urnk(lg),nurnk(lg,2))
  nurnk(lg,:)=1         !.number of processes participating in unification
  allocate(urnk(lg)%a(nurnk(lg,1),nurnk(lg,2)))
  urnk(lg)%a=myrank0 !MPI_PROC_NULL        !.lg entry is not used

!.These define the coarsest grid, which depends on the equation being solved
  if (eq==1 .OR. eq==2) then
    minCnx=2
    minCny=2
  elseif (eq==3) then
#IF (HDop==2 .OR. HDop==8)
    minCnx=3
    minCny=3
#ELIF (HDop==3 .OR. HDop==12)
    minCnx=4
    minCny=4
#ELIF (HDop==4 .OR. HDop==6)
    minCnx=2
    minCny=2
#ENDIF
  endif

  do
!...Begin assuming no unification happens then check later
    iunify=.FALSE.
    junify=.FALSE.
    iunisend=.FALSE.
    junisend=.FALSE.
!...Current cell size along each direction
    hx(1)=bLx/dble(iprocs(1)*nxf)
    hy(1)=bLy/dble(jprocs(1)*nyf)
    if ((nxc*iprocs(1) == minCnx) .OR. (nyc*jprocs(1) == minCny) &
        .OR. (lg==lgMAX)) exit !.reached last grid
    lg=lg+1             !.one more coarser grid
!...Current fine grid spacing (squared)
    hxs=hx(1)**2
    hys=hy(1)**2
!...This allocating and deallocating is just so that quantities
!...are stored in the right order in memory, which is finest
!...last and coarsest first.
    if (lg>2) then
!.....Temp coupling coefficients from previous grid
      allocate(tmp(lg-1,2))
      tmp=cf
      deallocate(cf)
!.....Temp previous grid dimensions
      allocate(tmpx(lg-1),tmpy(lg-1))
      tmpx=nxv
      tmpy=nyv
      deallocate(nxv,nyv)
!.....Temp previous indicator of which dimension was coarsened
      allocate(tmpw(lg-1))
      tmpw=wdc
      deallocate(wdc)
    endif
!...Temp previous cell sizes
    allocate(tmphx(lg-1),tmphy(lg-1))
    tmphx=hx
    tmphy=hy
    deallocate(hx,hy)
!...Allocate arrays with new number of grids
    allocate(cf(lg,2))
    allocate(nxv(lg),nyv(lg))
    allocate(wdc(lg))
    allocate(hx(lg),hy(lg))
    wdc=0 !.only affects the final wdc(1), which isn't used
!...Populate arrays with previous values stored in tmp arrays
    if (lg>2) then
      cf(3:lg,:)=tmp(2:lg-1,:)
      deallocate(tmp)
      nxv(3:lg)=tmpx(2:lg-1)
      nyv(3:lg)=tmpy(2:lg-1)
      deallocate(tmpx,tmpy)
      wdc(2:lg)=tmpw
      wdc(2)=wdc(3) !.temporarily
      deallocate(tmpw)
    endif
    hx(2:lg)=tmphx
    hy(2:lg)=tmphy
    deallocate(tmphx,tmphy)
!...Store current fine grid dimensions
    nxv(2)=nxf
    nyv(2)=nyf
    if (lg>2 .AND. wdc(2)==0) then ! .OR. (wdc() .ne. 0)) then
!.....If already reached a grid with balanced coupling coefficients
!.....coarsen using standard coarsening from now on
      cf(2,1)=1.0_dp
      cf(2,2)=1.0_dp
    else
!.....Otherwise compute the coupling factors + determine coarsening
!.....This identifies along which dimension is the coupling in the
!.....elliptic equation strongest. It is a combination of the difference
!.....in dx vs dy and the coefficients in the equation.

!.....Two options here
!!____________________ Option A __________________________________________!
!!_____ For the Poisson eqn this seems like a better informed approach ____!
!        allocate(tmp(nxf,nyf))
!........tmpD here is the coarsened Diff at every level (currently computed
!........after this subroutine is called)
!        do i=2,nxf-1
!!.........Original, based on exact entries in matrix
!!          tmp(i,:)=(tmpD(ng-lg+2)%a(i+1,:,1)+tmpD(ng-lg+2)%a(i-1,:,1) &
!!                   +2.0_dp*tmpD(ng-lg+2)%a(i,:,1))/hxs
!!.........This one requires less communication
!          tmp(i,:)=(tmpD(ng-lg+2)%a(i+1,:,1)+tmpD(ng-lg+2)%a(i-1,:,1) &
!                   +2.0_dp*tmpD(ng-lg+2)%a(i,:,1))/hxs
!        enddo
!        cf(2,1)=maxval(tmp(2:nxf-1,:))
!!        cf(2,1)=sum(tmp(2:nxf-1,:))/dble((nxf-2)*nyf)
!        do j=2,nyf-1
!          tmp(:,j)=(tmpD(ng-lg+2)%a(:,j+1,1)+tmpD(ng-lg+2)%a(:,j-1,1) &
!                   +2.0_dp*tmpD(ng-lg+2)%a(:,j,1))/hys
!        enddo
!        cf(2,2)=maxval(tmp(:,2:nyf-1))
!!        cf(2,2)=sum(tmp(:,2:nyf-1))/dble(nxf*(nyf-2))
!        deallocate(tmp)
!!____________________ Option B __________________________________________!
!.....It appears that using the same coupling coefficient for both the
!.....Poisson, modified Helmholtz and Laplacian-based hyperdiffusion
!.....works fine (see below). Then each processor can establish the
!.....communication patterns independently because the matrix coefficients
!.....are scalar constants.
#IF (HDop==8 .OR. HDop==12)
      if (eq==3) then
        cf(2,1)=iDiff(1)/(hxs**(HDop/4))
        cf(2,2)=iDiff(2)/(hys**(HDop/4))
      else
        cf(2,1)=1.0_dp/hxs
        cf(2,2)=1.0_dp/hys
      endif
#ELSE
      cf(2,1)=1.0_dp/hxs
      cf(2,2)=1.0_dp/hys
#ENDIF
!!__________________________________________________________________________!
    endif

    if (NIS==2) then
!.....Doubling grid spacing at every level
      if (cf(2,1)/cf(2,2)>1.30_dp) then !. 1.3 from Zubair 2007, requires testing
!.......Semi-coarsen along x
        nxc=nxf/2
        nyc=nyf
        wdc(2)=1
      elseif (cf(2,2)/cf(2,1)>1.30_dp) then !. 1.3 from Zubair 2007, requires testing
!.......Semi-coarsen along y
        nxc=nxf 
        nyc=nyf/2
        wdc(2)=2
      else
!.......Standard coarsening
        nxc=nxf/2
        nyc=nyf/2
        wdc(2)=0
      endif
    else !if (NIS==4) then
!.....Quadrupling grid spacing until problem is isotropic, then doubling
      if (cf(2,1)/cf(2,2)>1.01_dp) then !. 1.3 from Zubair 2007, requires testing
!.......Semi-coarsen along x
        nxc=nxf/4
        nyc=nyf
        wdc(2)=1
!.......If quadrupling the grid spacing overcoarsens the grid along that
!.......dimension so that the next step would require quadrupling along
!.......the other, then just double the grid spacing
#IF (HDop==8 .OR. HDop==12)
        if (eq==3) then
          tmpcfx=iDiff(1)/((bLx/(dble(nxc)*iprocs(1)))**HDopd4)
          tmpcfy=iDiff(2)/((bLy/(dble(nyc)*jprocs(1)))**HDopd4)
        else
          tmpcfx=1.0_dp/(bLx/(dble(nxc)*iprocs(1)))**2
          tmpcfy=1.0_dp/(bLy/(dble(nyc)*jprocs(1)))**2
        endif
#ELSE
        tmpcfx=1.0_dp/(bLx/(dble(nxc)*iprocs(1)))**2
        tmpcfy=1.0_dp/(bLy/(dble(nyc)*jprocs(1)))**2
#ENDIF
        if (tmpcfy/tmpcfx>1.01_dp) nxc=nxc*2
      elseif (cf(2,2)/cf(2,1)>1.01_dp) then !. 1.3 from Zubair 2007, requires testing
!.......Semi-coarsen along y
        nxc=nxf 
        nyc=nyf/4
        wdc(2)=2
!.......If quadrupling the grid spacing overcoarsens the grid along that
!.......dimension so that the next step would require quadrupling along
!.......the other, then just double the grid spacing
#IF (HDop==8 .OR. HDop==12)
        if (eq==3) then
          tmpcfx=iDiff(1)/((bLx/(dble(nxc)*iprocs(1)))**HDopd4)
          tmpcfy=iDiff(2)/((bLy/(dble(nyc)*jprocs(1)))**HDopd4)
        else
          tmpcfx=1.0_dp/(bLx/(dble(nxc)*iprocs(1)))**2
          tmpcfy=1.0_dp/(bLy/(dble(nyc)*jprocs(1)))**2
        endif
#ELSE
        tmpcfx=1.0_dp/(bLx/(dble(nxc)*iprocs(1)))**2
        tmpcfy=1.0_dp/(bLy/(dble(nyc)*jprocs(1)))**2
#ENDIF
        if (tmpcfx/tmpcfy>1.01_dp) nyc=nyc*2
      else
!.......Standard coarsening
        nxc=nxf/2
        nyc=nyf/2
        wdc(2)=0
      endif
    endif

    if (lg>2) then
!.....Temporarily repurpose tmphx/y to save coarse grid number of points
      allocate(tmphx(lg-1),tmphy(lg-1))
      tmphx=Cnx
      tmphy=Cny
      deallocate(Cnx,Cny)
    endif
!...Save local coarse grid dimensions
    allocate(Cnx(lg),Cny(lg))
    Cnx(2)=nxc
    Cny(2)=nyc
    if (lg>2) then
      Cnx(3:lg)=tmphx(2:(lg-1))
      Cny(3:lg)=tmphy(2:(lg-1))
      deallocate(tmphx,tmphy)
    endif

!...Enlarge array with number of processes along x to fit new coarser grid
    allocate(tmpproc(lg-1))
    tmpproc=iprocs
    deallocate(iprocs)
    allocate(iprocs(lg))
    iprocs(2:lg)=tmpproc
    iprocs(1)=iprocs(2)
    deallocatE(tmpproc)

!...Enlarge array with number of processes along y to fit new coarser grid
    allocate(tmpproc(lg-1))
    tmpproc=jprocs
    deallocate(jprocs)
    allocate(jprocs(lg))
    jprocs(2:lg)=tmpproc
    jprocs(1)=jprocs(2)
    deallocatE(tmpproc)

    igap0=igap
    jgap0=jgap
!...Check whether subdomain unification along x has to happen
    if ((.NOT. oneiproc) .AND. (nxc<nxmin)) then
!.....Coarse grid < min # of x cells per process allowed by user.
!.....After relaxation, computing and restricting residue, data is
!.....grouped in smaller number of processes.
      iunify=.TRUE.
      if (iprocs(2)>2) then
!.......The pattern is to agglomerate data towards the center of
!.......the domain, favouring iprocs/2-1 over iprocs/2
        if (nxmin/nxc==4) then
          igap=igap*4   !.Needed for quadruple coarsening (I think)
        else !if (nxmin/nxc==2) then
          igap=igap*2
        endif
!.......Determine if this process remains active in coarser grids
!.......On first half, odd (ID) ranks stay active and even ones don't.
        if (Lhalf0 .AND. (mod(iID0,igap)/=igap-1)) iactive=.FALSE.
!.......On second half, even ranks stay active, and odd ones don't.
        if ((.NOT. Lhalf0) .AND. (mod(iID0,igap)/=0)) iactive=.FALSE.
      else !if (iprocs(1)==2) then
        igap=igap*2
!.......only process with iID=iprocs0/2-1 remains active
        if (iID0/=iprocs0/2-1) iactive=.FALSE. 
      endif
    endif

!...Check whether subdomain unification along y has to happen
    if ((.NOT. onejproc) .AND.(nyc < nymin)) then
!.....Coarse grid has min # of y cells per process allowed by user.
!.....Unify subdomains. The pattern is to agglomerate data towards
!.....the center of the domain, favouring jprocs/2 over jprocs/2+1
      junify=.TRUE.
      if (jprocs(2)>2) then
        if (nymin/nyc==4) then
          jgap=jgap*4   !.Needed for quadruple coarsening (I think)
        else !if (nymin/nyc==2) then
          jgap=jgap*2
        endif

!.......Determine if this process remains active in coarser grids
!.......On first half, odd (ID) ranks stay active and even ones don't.
        if (Bhalf0 .AND. (mod(jID0,jgap)/=jgap-1)) jactive=.FALSE.
!.......On second half, even ranks stay active, and odd ones don't.
        if ((.NOT. Bhalf0) .AND. (mod(jID0,jgap)/=0)) jactive=.FALSE.
      else !if (jprocs(2)==2) then
        jgap=jgap*2
!.......only process with jID=jprocs0/2-1 remains active
        if (jID0/=jprocs0/2-1) jactive=.FALSE.
      endif
    endif

!...Number of processes remaining
    iprocs(1)=iprocs0/igap
    jprocs(1)=jprocs0/jgap
!...Identify when there's only one process left
    if (iprocs(1)==1) oneiproc=.TRUE.
    if (jprocs(1)==1) onejproc=.TRUE.

!...Enlarge active array to fit active processes in new coarser grid
    allocate(tmpact(lg-1))
    tmpact=active
    deallocate(active)
    allocate(active(lg))
    active(1)=.TRUE.
!...If not active along x or y, then this process is inactive
    if ((.NOT. iactive) .OR. (.NOT. jactive)) active(1)=.FALSE.
    active(2:lg)=tmpact
    deallocate(tmpact)

!...Enlarge uni array to fit new coarser grid
    allocate(tmpuni(lg-1))
    tmpuni=unify
    deallocate(unify)
    allocate(unify(lg))
    unify(1)=.FALSE.
!...If unification along x or y happened, set uni=TRUE
!deprecated: if (active(1) .AND. (iunify .OR. junify)) unify(1)=.TRUE.
    if (active(2) .AND. (iunify .OR. junify)) unify(1)=.TRUE.
    unify(2:lg)=tmpuni
    deallocate(tmpuni)

!...Temporarily save urnk and nurnk arrays initialized up to now
    allocate(tmpnurnk(lg-1,2),tmpurnk(lg-1))
    do j=1,lg-1
      allocate(tmpurnk(j)%a(nurnk(j,1),nurnk(j,2)))
    enddo
    tmpurnk=urnk 
    tmpnurnk=nurnk 
    deallocate(urnk,nurnk)
!...Enlarge urnk and nurnk to fit new coarser grid
    allocate(urnk(lg),nurnk(lg,2))
    nurnk(2:lg,:)=tmpnurnk
    do j=2,lg
      allocate(urnk(j)%a(nurnk(j,1),nurnk(j,2)))
      urnk(j)%a=tmpurnk(j-1)%a
      deallocate(tmpurnk(j-1)%a)
    enddo
    deallocate(tmpnurnk,tmpurnk)
!...Determine where data is sent to/received from when
!...subdomains are unified
!...MF: Sorry, the code gets a little messy in the next couple of hundred lines
    nurnk(1,:)=1
    inurnk=1
    jnurnk=1
    if (unify(1)) then  !.Subdomain unification happened

      if (iunify) then  !.Unification happened along x
        if (Lhalf0) then !.Left half
          if (.NOT. oneiproc) then
            if (mod(iID0,igap)/=igap-1) then
!.............These procs send data rightwards to iID0+igap0
              iunisend=.TRUE.
              allocate(iugap(1))
              iugap=igap0       !.Gap between this proc and receiving proc
            else
!.............These procs receive data from iID0-igap0 (on the left)
              allocate(iugap(2))
              iugap=(/-igap0,0/) !.Gap between this proc and sending procs (include self)
              inurnk=2    !.Unification of 2 subdomains along x
            endif
          else !.only one process remains and this is the receiving one
            allocate(iugap(2))
            iugap=(/0,1/)
            inurnk=2    !.Unification of 2 subdomains along x
!***********************************************************************
!deprecated: This stuff was used when I was trying to check for unification
!            along x and/or y at the same time (I think)
!
!            if (.NOT. junify) then
!              allocate(iugap(2))
!              iugap=(/0,1/)
!              inurnk=2    !.Unification of 2 subdomains along x
!            else
!              if (Bhalf0) then ! Lower half
!                if (mod(jID0,jgap)/=jgap-1) then
!!.................These procs send data upwards
!                  iunisend=.TRUE.
!                  allocate(iugap(1))
!                  iugap=0
!                else
!!.................These procs receive data from below
!                  allocate(iugap(2))
!                  iugap=(/0,1/)
!                  nurnk(1,1)=2    ! Unification of 2 subdomains along x
!                endif
!              else    ! Upper half
!                if (mod(jID0,jgap)/=0) then
!!.................These procs send data downwards
!                  iunisend=.TRUE.
!                  allocate(iugap(1))
!                  iugap=0
!                else
!!.................These procs receive data from above
!                  allocate(iugap(2))
!                  iugap=(/0,1/)
!                  nurnk(1,1)=2    ! Unification of 2 subdomains along x
!                endif
!              endif
!            endif
!***********************************************************************
          endif
        else    !.Right half
          if (.NOT. oneiproc) then
            if (mod(iID0,igap)/=0) then
!.............These procs send data leftwards to iID0-igap0
              iunisend=.TRUE.
              allocate(iugap(1)) !.Gap between this proc and receiving proc
              iugap=-igap0
            else
!.............These procs receive data from iID0+igap0 (on the right)
              allocate(iugap(2))
              iugap=(/0,igap0/) !.Gap between this proc and sending procs (include self)
              inurnk=2    !.Unification of 2 subdomains along x
            endif
          else !.only one process remains and this is the one sending data
            iunisend=.TRUE.
            allocate(iugap(1))
            iugap=-1
          endif
        endif
      else !.no unification along x takes place
        allocate(iugap(1))
        iugap=0
      endif !.end if (iunify)

      if (junify) then  !.Unification happened along y
        if (Bhalf0) then !.Lower half
          if (.NOT. onejproc) then
            if (mod(jID0,jgap)/=jgap-1) then
!.............These procs send data upwards to jID0+jgap0 
              junisend=.TRUE.
              allocate(jugap(1)) !.Gap between this proc and receiving proc
              jugap=jgap0
            else
!.............These procs receive data from jID0-jgap0 (below)
              allocate(jugap(2))
              jugap=(/-jgap0,0/) !.Gap between this proc and sending procs (include self)
              jnurnk=2    !.At least one unification (of 2 subdomains) along x
            endif
          else !.only one process remains and this is the receiving one
            allocate(jugap(2))
            jugap=(/0,1/)
            jnurnk=2    !.Unification of 2 subdomains along y
!***********************************************************************
!deprecated: This stuff was used when I was trying to check for unification
!            along x and/or y at the same time (I think)
!
!            if (.NOT. iunify) then
!              allocate(jugap(2))
!              jugap=(/0,1/)
!              jnurnk=2    ! Unification of 2 subdomains along y
!            else
!              if (Lhalf0) then ! Left half
!                if (mod(iID0,igap)/=igap-1) then
!!.................These procs send data rightwards to iID0+igap0
!                  junisend=.TRUE.
!                  allocate(jugap(1))
!                  jugap=0
!                else
!!.................These procs receive data from iID0-igap0 (on the left)
!                  allocate(jugap(2))
!                  jugap=(/0,1/)
!                  nurnk(1,2)=2    ! Unification of 2 subdomains along y
!                endif
!              else    ! Right half
!                if (mod(iID0,igap)/=0) then
!!.................These procs send data leftwards to iID0-igap0
!                  junisend=.TRUE.
!                  allocate(jugap(1))
!                  jugap=0
!                else
!!.................These procs receive data from iID0+igap0 (on the right)
!                  allocate(jugap(2))
!                  jugap=(/0,1/)
!                  nurnk(1,2)=2    ! Unification of 2 subdomains along y
!                endif
!              endif
!            endif
!***********************************************************************
          endif
        else    !.Upper half
          if (.NOT. onejproc) then
            if (mod(jID0,jgap)/=0) then
!.............These procs send data downwards to jID0-jgap0
              junisend=.TRUE.
              allocate(jugap(1)) !.Gap between this proc and receiving proc
              jugap=-jgap0
            else
!.............These procs receive data from jID0+jgap0 (above)
              allocate(jugap(2))
              jugap=(/0,jgap0/) !.Gap between this proc and sending procs (include self)
              jnurnk=2    !.Unification of 2 subdomains along y
            endif
          else !.only one process remains and this is the one sending data
            allocate(jugap(1))
            jugap=-1
          endif
        endif
      else !.no unification along y takes place
        allocate(jugap(1))
        jugap=0
      endif !.end if (junify)

!.................................................................................!
!.This section is a horrible but important HACK
!.Luckily this is only done at initialization
!.
      if ( (iunify .AND. (.NOT. junify)) .OR. &
           ((.NOT. iunify) .AND. junify) ) then
!.......Subdomain unification takes place along a single dimension
        nurnk(1,1)=inurnk
        nurnk(1,2)=jnurnk
      else
!.......Subdomain unification takes place along x and y at the same time
        if (Lhalf0 .AND. Bhalf0) then
!.........Bottom left quadrant
          if (( (.NOT. oneiproc) .AND. (.NOT. onejproc) .AND. &
                (((mod(iID0,igap)==igap-1) .AND. (mod(jID0,jgap)==jgap-1)) .OR. &
                 ((mod(iID0,igap)/=igap-1) .AND. (mod(jID0,jgap)/=jgap-1))) ) .OR. &
              (oneiproc .AND. (mod(jID0,jgap)==jgap-1)) .OR. &
              (onejproc .AND. (mod(iID0,igap)==igap-1)) ) then
            nurnk(1,1)=inurnk
            nurnk(1,2)=jnurnk
          elseif ((mod(iID0,igap)/=igap-1) .AND. (mod(jID0,jgap)==jgap-1)) then
            nurnk(1,1)=inurnk
            nurnk(1,2)=1
            deallocate(jugap)
            allocate(jugap(1))
            jugap=0
          elseif ((mod(iID0,igap)==igap-1) .AND. (mod(jID0,jgap)/=jgap-1)) then
            nurnk(1,1)=1
            nurnk(1,2)=jnurnk
            deallocate(iugap)
            allocate(iugap(1))
            iugap=0
          endif
        elseif ((.NOT. Lhalf0) .AND. Bhalf0) then
!.........Bottom right quadrant
          if (( (.NOT. onejproc) .AND. (((mod(iID0,igap)==0) .AND. (mod(jID0,jgap)==jgap-1)) &
               .OR. ((mod(iID0,igap)/=0) .AND. (mod(jID0,jgap)/=jgap-1)))) .OR. &
              (onejproc .AND. (mod(iID0,igap)==0))) then
            nurnk(1,1)=inurnk
            nurnk(1,2)=jnurnk
          elseif ((mod(iID0,igap)/=0) .AND. (mod(jID0,jgap)==jgap-1)) then
            nurnk(1,1)=inurnk
            nurnk(1,2)=1
            deallocate(jugap)
            allocate(jugap(1))
            jugap=0
          elseif ((mod(iID0,igap)==0) .AND. (mod(jID0,jgap)/=jgap-1)) then
            nurnk(1,1)=1
            nurnk(1,2)=jnurnk
            deallocate(iugap)
            allocate(iugap(1))
            iugap=0
          endif
        elseif (Lhalf0 .AND. (.NOT. Bhalf0)) then
!.........Top left quadrant
          if (((.NOT. oneiproc) .AND. (((mod(iID0,igap)==igap-1) .AND. (mod(jID0,jgap)==0)) &
               .OR. ((mod(iID0,igap)/=igap-1) .AND. (mod(jID0,jgap)/=0)))) .OR. &
              (oneiproc .AND. (mod(jID0,jgap)==0))) then
            nurnk(1,1)=inurnk
            nurnk(1,2)=jnurnk
          elseif ((mod(iID0,igap)/=igap-1) .AND. (mod(jID0,jgap)==0)) then
            nurnk(1,1)=inurnk
            nurnk(1,2)=1
            deallocate(jugap)
            allocate(jugap(1))
            jugap=0
          elseif ((mod(iID0,igap)==igap-1) .AND. (mod(jID0,jgap)/=0)) then
            nurnk(1,1)=1
            nurnk(1,2)=jnurnk
            deallocate(iugap)
            allocate(iugap(1))
            iugap=0
          endif
        else !if ((.NOT. Lhalf0) .AND. (.NOT. Bhalf0)) then
!.........Top right quadrant
          if (((mod(iID0,igap)==0) .AND. (mod(jID0,jgap)==0)) .OR. &
             ((mod(iID0,igap)/=0) .AND. (mod(jID0,jgap)/=0))) then
            nurnk(1,1)=inurnk
            nurnk(1,2)=jnurnk
          elseif ((mod(iID0,igap)/=0) .AND. (mod(jID0,jgap)==0)) then
            nurnk(1,1)=inurnk
            nurnk(1,2)=1
            deallocate(jugap)
            allocate(jugap(1))
            jugap=0
          elseif ((mod(iID0,igap)==0) .AND. (mod(jID0,jgap)/=0)) then
            nurnk(1,1)=1
            nurnk(1,2)=jnurnk
            deallocate(iugap)
            allocate(iugap(1))
            iugap=0
          endif
        endif
      endif
!....................... END OF HACK ..............................................!

!.....Convert coordinates of sending/receiving procs given by gaps
!.....into rank IDs in 2D communicator
      allocate(urnk(1)%a(nurnk(1,1),nurnk(1,2))) 
      do j=1,nurnk(1,2)
        do i=1,nurnk(1,1)
          call MPI_CART_RANK(XYcomm0,(/jID0+jugap(j),iID0+iugap(i)/),urnk(1)%a(i,j),ierror)
        enddo
      enddo
      deallocate(iugap,jugap)

    else
!.....In this case no subdomain unification takes place
      allocate(urnk(1)%a(1,1)) 
      urnk(1)%a=myrank0 !MPI_PROC_NULL
    endif !.end if (unify(1))
    
!...Next fine grid dimension is current coarse grid dimension
    nxf=nxc
    nyf=nyc
!...If unification takes place along one direction, active processes
!...have twice as many points as nxc,nyc in the next grid
    if (iunify) nxf=2*nxf
    if (junify) nyf=2*nyf

  enddo

!.Only the finest grid is used (relaxation only)
  if (lgMAX==1) then
    allocate(nxv(1),nyv(1),Cnx(1),Cny(1))
    allocate(wdc(1),cf(1,2))
    wdc=0
  endif
  nxv(1)=nxf
  nyv(1)=nyf
  Cnx(1)=nxf
  Cny(1)=nyf

!.Join active and unify into a single array
  allocate(acu(lg,2))
  acu(:,1)=active
  acu(:,2)=unify

!.This is not necessary because the last coupling coefficient is not used
!  hxs=(bLx/dble(nxf))**2
!  hys=(bLy/dble(nyf))**2
  hxs=hx(1)**2
  hys=hy(1)**2
  if (wdc(1)==0) then
    cf(1,1)=1.0_dp
    cf(1,2)=1.0_dp
  else
!.....Two options here
!!____________________ Option A __________________________________________!
!!_____ For the Poisson eqn this seems like a better informed approach ____!
!      allocate(tmp(nxf,nyf))
!      do i=2,nxf-1
!        tmp(i,:)=(Diff(1)%a(i+1,:,1)+Diff(1)%a(i-1,:,1)+2.0_dp*Diff(1)%a(i,:,1))/hxs
!      enddo
!      cf(1,1)=maxval(tmp(2:nxf-1,:))
!!      cf(1,1)=sum(tmp(2:nxf-1,:))/dble((nxf-2)*nyf)
!      do j=2,nyf-1
!        tmp(:,j)=(Diff(1)%a(:,j+1,1)+Diff(1)%a(:,j-1,1)+2.0_dp*Diff(1)%a(:,j,1))/hys
!      enddo
!      cf(1,2)=maxval(tmp(:,2:nyf-1))
!!      cf(1,2)=sum(tmp(:,2:nyf-1))/dble(nxf*(nyf-2))
!      deallocate(tmp)
!!____________________ Option B __________________________________________!
!...This seems to work fine
#IF (HDop==8 .OR. HDop==12)
    if (eq==3) then
      cf(1,1)=iDiff(1)/(hxs**(HDop/4))
      cf(1,2)=iDiff(2)/(hys**(HDop/4))
    else
      cf(1,1)=1.0_dp/hxs
      cf(1,2)=1.0_dp/hys
    endif
#ELSE
    cf(1,1)=1.0_dp/hxs
    cf(1,2)=1.0_dp/hys
#ENDIF
!!__________________________________________________________________________!
  endif

#IF (DispG==1)
!.Display parameters at each grid for processes selected by i,j below
if (eq==3 .AND. LGMAX==3) then
do j=1,1
do i=0,3
  if (iID0==i .AND. jID0==j) then
    call MPI_CART_RANK(XYcomm0,(/jID0,iID0/),k,ierror)
    write(*,'("myrank=",I3)') k
    do ng=lg,1,-1
    if (nurnk(ng,1)*nurnk(ng,2)==1) then
    write(*,'("ng=",I2," | nx=",I4," ny=",I4," | nxC=",I4," nyC=",I4," | uni=", &
      L2," act=",L2," | iprocs=",i3," jprocs=",i3," | nurnk=",2I3," urnk=",I4)') &
      ng,nxv(ng),nyv(ng),Cnx(ng),Cny(ng),unify(ng),active(ng), &
      iprocs(ng),jprocs(ng),nurnk(ng,:),urnk(ng)%a(:,:)
    elseif (nurnk(ng,1)*nurnk(ng,2)==2) then
    write(*,'("ng=",I2," | nx=",I4," ny=",I4," | nxC=",I4," nyC=",I4," | uni=", &
      L2," act=",L2," | iprocs=",i3," jprocs=",i3," | nurnk=",2I3," urnk=",2I4)') &
      ng,nxv(ng),nyv(ng),Cnx(ng),Cny(ng),unify(ng),active(ng), &
      iprocs(ng),jprocs(ng),nurnk(ng,:),urnk(ng)%a(:,:)
    elseif (nurnk(ng,1)*nurnk(ng,2)==4) then
    write(*,'("ng=",I2," | nx=",I4," ny=",I4," | nxC=",I4," nyC=",I4," | uni=", &
      L2," act=",L2," | iprocs=",i3," jprocs=",i3," | nurnk=",2I3," urnk=",4I4)') &
      ng,nxv(ng),nyv(ng),Cnx(ng),Cny(ng),unify(ng),active(ng), &
      iprocs(ng),jprocs(ng),nurnk(ng,:),urnk(ng)%a(:,:)
    endif
    enddo
    print*,'-----------------------------------------------------'
  endif
  call MPI_BARRIER(XYcomm0,ierror)
enddo
enddo
endif
#ENDIF

  end subroutine pgrids
!*********************************************************
  subroutine terminate_PMG
!.Deallocate all the arrays used by PMG solver
  implicit none
  integer(I4B) :: ng,i,j,ierr
  
  do ng=lgph,1,-1
    deallocate(Diff(ng)%a)
  enddo
  deallocate(Diff)
  deallocate(wdcph,nxFph,nyFph,nxCph,nyCph)
  deallocate(wdcps,nxFps,nyFps,nxCps,nyCps)
  deallocate(iIDph,jIDph,iIDps,jIDps)
  deallocate(iprocsph,jprocsph,iprocsps,jprocsps)
  deallocate(nborph,nborps)
  do ng=lgph,1,-1
    if (acuph(ng,1)) then
      do i=1,3
        call MPI_COMM_FREE(COMMph(ng,i),ierr)
      enddo
    endif
    deallocate(urnkph(ng)%a)
  enddo
  do ng=lgps,1,-1
    if (acups(ng,1)) then
      do i=1,3
        call MPI_COMM_FREE(COMMps(ng,i),ierr)
      enddo
    endif
    deallocate(urnkps(ng)%a)
  enddo
  deallocate(COMMph,COMMps,acuph,acups)
  deallocate(urnkph,urnkps,nurnkph,nurnkps)
  do j=1,nqhd
    deallocate(wdchd(j)%a,nxFhd(j)%a,nyFhd(j)%a,nxChd(j)%a,nyChd(j)%a)
    deallocate(iIDhd(j)%a,jIDhd(j)%a,iprocshd(j)%a,jprocshd(j)%a)
    deallocate(nborhd(j)%a)
    do ng=lghd(j),1,-1
      if (acuhd(j)%a(ng,1)) then
        do i=1,3
          call MPI_COMM_FREE(COMMhd(j)%a(ng,i),ierr)
        enddo
      endif
      deallocate(urnkhd(j)%a(ng)%a)
    enddo
    deallocate(COMMhd(j)%a,acuhd(j)%a,urnkhd(j)%a,nurnkhd(j)%a)
  enddo
  deallocate(lghd,wdchd,nxFhd,nyFhd,nxChd,nyChd)
  deallocate(iIDhd,jIDhd,iprocshd,jprocshd,nborhd)
  deallocate(COMMhd,acuhd,urnkhd,nurnkhd)
  deallocate(LGMAXhd,NIShd,GAMhd,BEThd,NU1hd,NU2hd)
  deallocate(hDiff)
  end subroutine terminate_PMG
!*********************************************************
  recursive subroutine Gcycphi(ng,u,lambda,rhs,xBC,ur)
!.Gamma cycles to solve Poisson equation and obtain phi
!                 div(lambda*grad(u))=rho
!.where lambda=plasma density, u=phi and rho=w.
!.INPUTS
!. integer ng: current grid
!. real(dp) dimension(:,:,:) u: current solution at grid ng
!. real(dp) dimension(:,:,:) lambda: spatially dependent coefficient
!. real(dp) dimension(:,:,:) rhs: the right-hand side of Poisson equation
!. integer(2) xBC: boundary conditions along x
!. real(dp) dimension(:,:) ur: Dirichlet value for Right boundary
!.OUTPUTS
!. real(dp) dimension(:,:,:) u: updated solution at grid ng
  integer(I4B), intent(in) :: ng
  real(DP), intent(inout) :: u(:,:,:)
  type(allo3d), allocatable :: lambda(:)
  real(DP), allocatable, intent(in) :: rhs(:,:,:)
  real(DP), intent(in), optional :: ur(:,:)
  integer(I4B), intent(in) :: xBC(:)
  integer(I4B) :: jpost, jpre, i,j
  real(DP), allocatable :: res(:,:,:),v(:,:,:)
  real(DP), allocatable :: ur0(:,:)

  iprocsm1=iprocsph(ng)-1 !.ID of last rank in x COMM
  jprocsm1=jprocsph(ng)-1 !.ID of last rank in y COMM

  if (ng==1) then        !.Bottom of gamma cycle: Solve on coarsest grid
#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do i=1,NU2ph
      call relaxphi(ng,u,lambda(ng)%a,rhs,xBC,ur)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    postXtph=postXtph+count2-count1
#ENDIF
  else                  !.On downward stroke towards coarser grids
#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpre=1,NU1ph      !.Pre-smoothing
      call relaxphi(ng,u,lambda(ng)%a,rhs,xBC,ur)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    preXtph=preXtph+count2-count1
#ENDIF

!...Number of grid points in coarse grid 
    nxc=nxCph(ng)
    nyc=nyCph(ng)
!...Number of grid points in coarse grid after subdomain unification
    nxc=nxCph(ng)
    nxu=nxFph(ng-1)
    nyu=nyFph(ng-1)
    allocate(res(nxu,nyu,nzL),v(nxu,nyu,nzL))
    allocate(ur0(nyu,nzL))
#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
!...Restricted residual is the next RHS 
    res=rstrct(resphi(ng,u,lambda(ng)%a,rhs,xBC,ur),wdcph(ng), &
      xBC,acuph(ng-1,:),urnkph(ng-1)%a,nurnkph(ng-1,:),ur)
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    rNrtph=rNrtph+count2-count1
#ENDIF

!...Next fine grid is the coarse, or unified subdomain, grid
    nxf=nxu
    nyf=nyu

    v=0.0_dp               !.Zero for initial guess in next relaxation
    ur0=0.0_dp             !.Zero right side boundary value for Dirichlect BC

    if (acuph(ng-1,1)) then
!.....Only active processes proceed to next coarser grid
      do i=1,GAMph
!.......Recursive call for the coarse grid correction
        call Gcycphi(ng-1,v,lambda,res,(/xBC(1),-1/),ur0)
      enddo
    endif

!...Redefine quantities that were modified in coarser grids
    iprocsm1=iprocsph(ng)-1
    jprocsm1=jprocsph(ng)-1
    nxc=nxCph(ng)
    nyc=nyCph(ng)
    nxf=nxFph(ng)
    nyf=nyFph(ng)

#IF (MGTIMER==1)
    call CPU_TIME(prolt1)
#ENDIF
!...Upward stroke -> finer grid: prolong coarse grid error + correct solution
    u=u+prolng(v,wdcph(ng),xBC,COMMph(ng,:),iIDph(ng),jprocsph(ng), &
      nborph(ng,:),acuph(ng-1,:),urnkph(ng-1)%a,nurnkph(ng-1,:),ur0)
#IF (MGTIMER==1)
    call CPU_TIME(prolt2)
    proltph=proltph+prolt2-prolt1
#ENDIF

#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpost=1,NU2ph      !.Post-smoothing
      call relaxphi(ng,u,lambda(ng)%a,rhs,xBC,ur)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    postXtph=postXtph+count2-count1
#ENDIF
  endif

  end subroutine Gcycphi
!*********************************************************
  subroutine relaxphi(ng,u,lam,rhs,xBC,ur)
!.Red-black Gauss-Seidel relaxation of
!       div(lam*grad(u))=rhs
!.using a 2nd order accurate finite difference stencil.
!.This is discretized such that the operator is self-adjoint.
!.INPUTS
!. integer ng: current grid
!. real(dp) dimension(:,:,:) u: current solution at grid ng
!. real(dp) dimension(:,:,:) lam: spatially dependent coefficient
!. real(dp) dimension(:,:,:) rhs: the right-hand side of Poisson equation
!. integer(2) xBC: boundary conditions along x
!. real(dp) dimension(:,:) ur: Dirichlet value for Right boundary
!.OUTPUTS
!. real(dp) dimension(:,:,:) u: updated solution at grid ng
  implicit none
  include 'mpif.h'
  real(DP), intent(inout) :: u(:,:,:)
  real(DP), intent(in), optional :: ur(:,:)
  real(DP), allocatable, intent(in) :: lam(:,:,:),rhs(:,:,:)
  integer(I4B), intent(in) :: ng,xBC(:)
  real(DP), allocatable :: Lsu(:),loc1D(:) !,Rsu(:,:)
  integer(I4B) :: i,j,k,il,nxfd2,nyfd2,Nrank,Srank,Erank,Wrank
  real(DP) :: fx,fy,fc
  real(DP) :: Lqu,Rqu,Lru,Rru
  real(DP), allocatable, dimension(:,:) :: sbufW,rbufE,sbufE,rbufW, &
                                           sbufN,rbufS,sbufS,rbufN
  integer(I4B) :: req(8),stat(MPI_STATUS_SIZE,8),ierr

  nxfd2=nxf/2
  nyfd2=nyf/2

!.Neigboring ranks to communicate with
  Nrank=nborph(ng,1)
  Erank=nborph(ng,3)
  Srank=nborph(ng,5)
  Wrank=nborph(ng,7)

!.Boundary conditions are implemented with some amount of flexibility.
!.The right boundary for example, if u(nxf+1) is needed, it will be
!.written as
!.  u(nxf+1) = Rq*u(nxf) + Rr*u(nxf-1) + Rs
!.so that symmetry, extrapolation or Dirichlet conditions can be
!.easily imposed.
!......... BCs of u ................................!
!.LHS BC of u
  allocate(Lsu(nzL),loc1D(nzL))
  if (xBC(1)==0) then
!...linear extrapolation w/o boundary value
    Lqu =  2.0_dp
    Lru = -1.0_dp
    Lsu =  0.0_dp
  elseif (xBC(1)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Lqu = -1.0_dp
    Lru =  0.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u(1,:,k))
    call MPI_ALLreduce(loc1D,Lsu,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMph(ng,1),ierr)
    Lsu =  2.0_dp*Lsu/dble(nyf*jprocsph(ng))
  else ! if (xBC(1)==1 .OR. xBC(1)==-1) then !.symmetry
    Lqu = dble(xBC(1))
    Lru = 0.0_dp
    Lsu = 0.0_dp
  endif
!.RHS BC of u
!  allocate(Rsu(nyf,nzL))
  if (xBC(2)==0) then
!...linear extrapolation w/o boundary value
    Rqu =  2.0_dp
    Rru = -1.0_dp
!    Rsu =  0.0_dp
!!...linear extrapolation w/ boundary value
!    Rqu = -1.0_dp
!    Rru =  0.0_dp
!    Rsu =  2.0_dp*ur
!  elseif (xBC(2)==2) then
!!...ky=0 component has even symmetry, ky.neq.0 component is odd
!    Rqu = -1.0_dp
!    Rru =  0.0_dp
!    forall(k=1:nzL) loc1D(k)=sum(u(nxf,:,k))
!    call MPI_ALLreduce(loc1D,Rsu,nzL,MPI_DOUBLE_PRECISION, &
!                       MPI_SUM,COMMph(ng,1),ierr)
!    Rsu =  2.0_dp*Rsu/dble(nyf*jprocsph(ng))
  else ! if (xBC(2)==1 .OR. xBC(2)==-1) then !.symmetry
    Rqu = dble(xBC(2))
    Rru = 0.0_dp
!    Rsu = 0.0_dp
  endif
!......... BCs of lambda ................................!
!.LHS BC of lambda assumed even
!.RHS BC of lambda assumed even
!........................................................!
  fx = 0.5_dp*((dble(nxf*iprocsph(ng))/bLx)**2)
  fy = 0.5_dp*((dble(nyf*jprocsph(ng))/bLy)**2)
  fc = 2.0_dp*(fx+fy)
  if (nxf>2) then ! Used in one sided stencil at RHS Dirichlet boundary
    il=nxf-2
  else
    il=nxf-1
  endif

!.Send buffer names indicate direction data is sent to
!.and receive buffers direction data is received from.
! CAREFUL!!!!
!.When sending/receiving lam data, the names don't apply. Try to use
!.the opposite of what the name says: e.g. sbufE sends lam data West,
!.and rbufS receives lam data from North.
  allocate(sbufW(nyfd2,nzL),rbufE(nyfd2,nzL),sbufE(nyfd2,nzL),rbufW(nyfd2,nzL))
  allocate(sbufN(nxfd2,nzL),rbufS(nxfd2,nzL),sbufS(nxfd2,nzL),rbufN(nxfd2,nzL))
!.Comunicate black grid points needed
!.Send bottom row of u to South rank
  sbufS=u(2:nxf:2,1,:)
  call MPI_ISEND(sbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(1),ierr)
!.Receive u row from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(2),ierr)
!.Send bottom row of lam to South rank
  sbufN=lam(2:nxf:2,1,:)
  call MPI_ISEND(sbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(3),ierr)
!.Receive lam row from North rank, for lam(i,j+1)
  call MPI_IRECV(rbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(4),ierr)
!.Send first column of u to West rank
  sbufW=u(1,2:nyf:2,:)
  call MPI_ISEND(sbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(5),ierr)
!.Receive u column from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(6),ierr)
!.Send first column of lam to West rank
  sbufE=lam(1,2:nyf:2,:)
  call MPI_ISEND(sbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(7),ierr)
!.Receive lam column from East rank, for lam(i+1,j)
  call MPI_IRECV(rbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(8),ierr)

!.First do the even-even and odd-odd squares of the grid,
!.i.e., the red squares of the checker-board
  do j=2,nyf-2,2
    do i=2,nxf-2,2
      u(i,j,:)=Oomeph*u(i,j,:)+omeph*(fx*((lam(i+1,j,:)+lam(i,j,:))*u(i+1,j,:)+(lam(i,j,:)+lam(i-1,j,:))*u(i-1,j,:)) &
        +fy*((lam(i,j+1,:)+lam(i,j,:))*u(i,j+1,:)+(lam(i,j,:)+lam(i,j-1,:))*u(i,j-1,:)) &
        -rhs(i,j,:))/(fx*(lam(i+1,j,:)+lam(i-1,j,:))+fy*(lam(i,j+1,:)+lam(i,j-1,:))+fc*lam(i,j,:))
    enddo
  enddo
!.Top, even i
  call MPI_WAITALL(2,(/req(2),req(4)/),(/stat(:,2),stat(:,4)/),ierr)
  do i=2,nxf-2,2
    u(i,nyf,:)=Oomeph*u(i,nyf,:)+omeph*(fx*((lam(i+1,nyf,:)+lam(i,nyf,:))*u(i+1,nyf,:) &
      +(lam(i,nyf,:)+lam(i-1,nyf,:))*u(i-1,nyf,:))+fy*((rbufS(i/2,:)+lam(i,nyf,:))*rbufN(i/2,:) &
      +(lam(i,nyf,:)+lam(i,nyf-1,:))*u(i,nyf-1,:))-rhs(i,nyf,:))/(fx*(lam(i+1,nyf,:)+lam(i-1,nyf,:)) &
      +fy*(rbufS(i/2,:)+lam(i,nyf-1,:))+fc*lam(i,nyf,:))
  enddo
!.Right boundary, even j
    call MPI_WAITALL(2,(/req(6),req(8)/),(/stat(:,6),stat(:,8)/),ierr)
  if (iIDph(ng) .eq. iprocsm1) then !.Right most process (along x)
    do j=2,nyf-2,2
      u(nxf,j,:)=Oomeph*u(nxf,j,:)+omeph*(fx*0.25_dp*(2.0_dp*lam(nxf,j,:)*(6.0_dp*ur(j,:)+u(il,j,:)) &
        +(-2.0_dp*lam(nxf,j,:)+4.0_dp*lam(nxf-1,j,:))*u(nxf-1,j,:)) &
        +fy*((lam(nxf,j+1,:)+lam(nxf,j,:))*u(nxf,j+1,:)+(lam(nxf,j,:)+lam(nxf,j-1,:))*u(nxf,j-1,:))-rhs(nxf,j,:)) &
        /(fx*(lam(nxf,j,:)+lam(nxf-1,j,:))+fy*(lam(nxf,j+1,:)+lam(nxf,j-1,:))+fc*lam(nxf,j,:))
    enddo
    u(nxf,nyf,:)=Oomeph*u(nxf,nyf,:)+omeph*(fx*0.25_dp*(2.0_dp*lam(nxf,nyf,:)*(6.0_dp*ur(nyf,:)+u(il,nyf,:)) &
      +(-2.0_dp*lam(nxf,nyf,:)+4.0_dp*lam(nxf-1,nyf,:))*u(nxf-1,nyf,:)) &
      +fy*((rbufS(nxfd2,:)+lam(nxf,nyf,:))*rbufN(nxfd2,:)+(lam(nxf,nyf,:)+lam(nxf,nyf-1,:))*u(nxf,nyf-1,:))-rhs(nxf,nyf,:)) &
      /(fx*(lam(nxf,nyf,:)+lam(nxf-1,nyf,:))+fy*(rbufS(nxfd2,:)+lam(nxf,nyf-1,:))+fc*lam(nxf,nyf,:)) !.global NE corner
  else
    do j=2,nyf-2,2
      u(nxf,j,:)=Oomeph*u(nxf,j,:)+omeph*(fx*((rbufW(j/2,:)+lam(nxf,j,:))*rbufE(j/2,:)+(lam(nxf,j,:)+lam(nxf-1,j,:))*u(nxf-1,j,:)) &
        +fy*((lam(nxf,j+1,:)+lam(nxf,j,:))*u(nxf,j+1,:)+(lam(nxf,j,:)+lam(nxf,j-1,:))*u(nxf,j-1,:)) &
        -rhs(nxf,j,:))/(fx*(rbufW(j/2,:)+lam(nxf-1,j,:))+fy*(lam(nxf,j+1,:)+lam(nxf,j-1,:))+fc*lam(nxf,j,:))
    enddo
    u(nxf,nyf,:)=Oomeph*u(nxf,nyf,:)+omeph*(fx*((rbufW(nyfd2,:)+lam(nxf,nyf,:))*rbufE(nyfd2,:) &
      +(lam(nxf,nyf,:)+lam(nxf-1,nyf,:))*u(nxf-1,nyf,:))+fy*((rbufS(nxfd2,:)+lam(nxf,nyf,:))*rbufN(nxfd2,:) &
      +(lam(nxf,nyf,:)+lam(nxf,nyf-1,:))*u(nxf,nyf-1,:))-rhs(nxf,nyf,:))/(fx*(rbufW(nyfd2,:)+lam(nxf-1,nyf,:)) &
      +fy*(rbufS(nxfd2,:)+lam(nxf,nyf-1,:))+fc*lam(nxf,nyf,:)) !.local NE corner
  endif

  call MPI_WAITALL(2,(/req(1),req(3)/),(/stat(:,1),stat(:,3)/),ierr)
  call MPI_WAITALL(2,(/req(5),req(7)/),(/stat(:,5),stat(:,7)/),ierr)
!.Send top row of array to North rank
  sbufN=u(1:nxf-1:2,nyf,:)
  call MPI_ISEND(sbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(1),ierr)
!.Receive u row from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(2),ierr)
!.Send top row of lam to North rank
  sbufS=lam(1:nxf-1:2,nyf,:)
  call MPI_ISEND(sbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(3),ierr)
!.Receive lam row from South rank, for lam(i,j-1)
  call MPI_IRECV(rbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(4),ierr)
!.Send last column of array to East rank
  sbufE=u(nxf,1:nyf-1:2,:)
  call MPI_ISEND(sbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(5),ierr)
!.Receive u column from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(6),ierr)
!.Send last column of lam to East rank
  sbufW=lam(nxf,1:nyf-1:2,:)
  call MPI_ISEND(sbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(7),ierr)
!.Receive lam column from West rank, for lam(i-1,j)
  call MPI_IRECV(rbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(8),ierr)
  do j=3,nyf-1,2
    do i=3,nxf-1,2
      u(i,j,:)=Oomeph*u(i,j,:)+omeph*(fx*((lam(i+1,j,:)+lam(i,j,:))*u(i+1,j,:)+(lam(i,j,:)+lam(i-1,j,:))*u(i-1,j,:)) &
        +fy*((lam(i,j+1,:)+lam(i,j,:))*u(i,j+1,:)+(lam(i,j,:)+lam(i,j-1,:))*u(i,j-1,:)) &
        -rhs(i,j,:))/(fx*(lam(i+1,j,:)+lam(i-1,j,:))+fy*(lam(i,j+1,:)+lam(i,j-1,:))+fc*lam(i,j,:))
    enddo
  enddo
!.Bottom, odd i
  call MPI_WAITALL(2,(/req(2),req(4)/),(/stat(:,2),stat(:,4)/),ierr)
  do i=3,nxf-1,2
    u(i,1,:)=Oomeph*u(i,1,:)+omeph*(fx*((lam(i+1,1,:)+lam(i,1,:))*u(i+1,1,:)+(lam(i,1,:)+lam(i-1,1,:))*u(i-1,1,:)) &
      +fy*((lam(i,2,:)+lam(i,1,:))*u(i,2,:)+(lam(i,1,:)+rbufN((i+1)/2,:))*rbufS((i+1)/2,:)) &
      -rhs(i,1,:))/(fx*(lam(i+1,1,:)+lam(i-1,1,:))+fy*(lam(i,2,:)+rbufN((i+1)/2,:))+fc*lam(i,1,:))
  enddo
!.Left boundary, odd j
  call MPI_WAITALL(2,(/req(6),req(8)/),(/stat(:,6),stat(:,8)/),ierr)
  if (iIDph(ng)==0) then !.Left most process (along x)
    u(1,1,:)=Oomeph*u(1,1,:)+omeph*(fx*((lam(2,1,:)+(1.0_dp+2.0_dp*Lru)*lam(1,1,:))*u(2,1,:)+2.0_dp*lam(1,1,:)*Lsu) &
      +fy*((lam(1,2,:)+lam(1,1,:))*u(1,2,:)+(lam(1,1,:)+rbufN(1,:))*rbufS(1,:)) &
      -rhs(1,1,:))/(fx*(lam(2,1,:)+(1.0_dp-2.0_dp*Lqu)*lam(1,1,:))+fy*(lam(1,2,:)+rbufN(1,:))+fc*lam(1,1,:)) !.global SW corner
    do j=3,nyf-1,2
      u(1,j,:)=Oomeph*u(1,j,:)+omeph*(fx*((lam(2,j,:)+(1.0_dp+2.0_dp*Lru)*lam(1,j,:))*u(2,j,:)+2.0_dp*lam(1,j,:)*Lsu) &
        +fy*((lam(1,j+1,:)+lam(1,j,:))*u(1,j+1,:)+(lam(1,j,:)+lam(1,j-1,:))*u(1,j-1,:)) &
        -rhs(1,j,:))/(fx*(lam(2,j,:)+(1.0_dp-2.0_dp*Lqu)*lam(1,j,:))+fy*(lam(1,j+1,:)+lam(1,j-1,:))+fc*lam(1,j,:))
    enddo
  else
    u(1,1,:)=Oomeph*u(1,1,:)+omeph*(fx*((lam(2,1,:)+lam(1,1,:))*u(2,1,:)+(lam(1,1,:)+rbufE(1,:))*rbufW(1,:)) &
      +fy*((lam(1,2,:)+lam(1,1,:))*u(1,2,:)+(lam(1,1,:)+rbufN(1,:))*rbufS(1,:)) &
      -rhs(1,1,:))/(fx*(lam(2,1,:)+rbufE(1,:))+fy*(lam(1,2,:)+rbufN(1,:))+fc*lam(1,1,:)) !.local SW corner
    do j=3,nyf-1,2
      u(1,j,:)=Oomeph*u(1,j,:)+omeph*(fx*((lam(2,j,:)+lam(1,j,:))*u(2,j,:)+(lam(1,j,:)+rbufE((j+1)/2,:))*rbufW((j+1)/2,:)) &
        +fy*((lam(1,j+1,:)+lam(1,j,:))*u(1,j+1,:)+(lam(1,j,:)+lam(1,j-1,:))*u(1,j-1,:)) &
        -rhs(1,j,:))/(fx*(lam(2,j,:)+rbufE((j+1)/2,:))+fy*(lam(1,j+1,:)+lam(1,j-1,:))+fc*lam(1,j,:))
    enddo
  endif

    call MPI_WAITALL(2,(/req(1),req(3)/),(/stat(:,1),stat(:,3)/),ierr)
    call MPI_WAITALL(2,(/req(5),req(7)/),(/stat(:,5),stat(:,7)/),ierr)
!.Now do even-odd and odd-even squares of the grid, i.e, the black
!.squares of the checker-board
!.Comunicate red grid points needed
!.Send bottom row of u to South rank
  sbufS=u(1:nxf-1:2,1,:)
  call MPI_ISEND(sbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(1),ierr)
!.Receive u row from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(2),ierr)
!.Send bottom row of lam to South rank
  sbufN=lam(1:nxf-1:2,1,:)
  call MPI_ISEND(sbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(3),ierr)
!.Receive lam row from North rank, for lam(i,j+1)
  call MPI_IRECV(rbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(4),ierr)
!.Send last column of array to East rank
  sbufE=u(nxf,2:nyf:2,:)
  call MPI_ISEND(sbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(5),ierr)
!.Receive u column from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(6),ierr)
!.Send last column of lam to East rank
  sbufW=lam(nxf,2:nyf:2,:)
  call MPI_ISEND(sbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(7),ierr)
!.Receive lam column from West rank, for lam(i-1,j)
  call MPI_IRECV(rbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(8),ierr)
  do j=2,nyf-2,2
    do i=3,nxf-1,2
      u(i,j,:)=Oomeph*u(i,j,:)+omeph*(fx*((lam(i+1,j,:)+lam(i,j,:))*u(i+1,j,:)+(lam(i,j,:)+lam(i-1,j,:))*u(i-1,j,:)) &
        +fy*((lam(i,j+1,:)+lam(i,j,:))*u(i,j+1,:)+(lam(i,j,:)+lam(i,j-1,:))*u(i,j-1,:)) &
        -rhs(i,j,:))/(fx*(lam(i+1,j,:)+lam(i-1,j,:))+fy*(lam(i,j+1,:)+lam(i,j-1,:))+fc*lam(i,j,:))
    enddo
  enddo
!.Top, odd i
  call MPI_WAITALL(2,(/req(2),req(4)/),(/stat(:,2),stat(:,4)/),ierr)
  do i=3,nxf-1,2
    u(i,nyf,:)=Oomeph*u(i,nyf,:)+omeph*(fx*((lam(i+1,nyf,:)+lam(i,nyf,:))*u(i+1,nyf,:) &
      +(lam(i,nyf,:)+lam(i-1,nyf,:))*u(i-1,nyf,:))+fy*((rbufS((i+1)/2,:)+lam(i,nyf,:))*rbufN((i+1)/2,:) &
      +(lam(i,nyf,:)+lam(i,nyf-1,:))*u(i,nyf-1,:))-rhs(i,nyf,:))/(fx*(lam(i+1,nyf,:) &
      +lam(i-1,nyf,:))+fy*(rbufS((i+1)/2,:)+lam(i,nyf-1,:))+fc*lam(i,nyf,:))
  enddo
!.Left boundary, even j
  call MPI_WAITALL(2,(/req(6),req(8)/),(/stat(:,6),stat(:,8)/),ierr)
  if (iIDph(ng)==0) then !.Left most process (along x)
    do j=2,nyf-2,2
      u(1,j,:)=Oomeph*u(1,j,:)+omeph*(fx*((lam(2,j,:)+(1.0_dp+2.0_dp*Lru)*lam(1,j,:))*u(2,j,:)+2.0_dp*lam(1,j,:)*Lsu) &
        +fy*((lam(1,j+1,:)+lam(1,j,:))*u(1,j+1,:)+(lam(1,j,:)+lam(1,j-1,:))*u(1,j-1,:)) &
        -rhs(1,j,:))/(fx*(lam(2,j,:)+(1.0_dp-2.0_dp*Lqu)*lam(1,j,:))+fy*(lam(1,j+1,:)+lam(1,j-1,:))+fc*lam(1,j,:))
    enddo
    u(1,nyf,:)=Oomeph*u(1,nyf,:)+omeph*(fx*((lam(2,nyf,:)+(1.0_dp+2.0_dp*Lru)*lam(1,nyf,:))*u(2,nyf,:)+2.0_dp*lam(1,nyf,:)*Lsu) &
      +fy*((rbufS(1,:)+lam(1,nyf,:))*rbufN(1,:)+(lam(1,nyf,:)+lam(1,nyf-1,:))*u(1,nyf-1,:)) &
      -rhs(1,nyf,:))/(fx*(lam(2,nyf,:)+(1.0_dp-2.0_dp*Lqu)*lam(1,nyf,:))+fy*(rbufS(1,:)+lam(1,nyf-1,:))+fc*lam(1,nyf,:)) !.global NW corner
  else
    do j=2,nyf-2,2
      u(1,j,:)=Oomeph*u(1,j,:)+omeph*(fx*((lam(2,j,:)+lam(1,j,:))*u(2,j,:)+(lam(1,j,:)+rbufE(j/2,:))*rbufW(j/2,:)) &
        +fy*((lam(1,j+1,:)+lam(1,j,:))*u(1,j+1,:)+(lam(1,j,:)+lam(1,j-1,:))*u(1,j-1,:)) &
        -rhs(1,j,:))/(fx*(lam(2,j,:)+rbufE(j/2,:))+fy*(lam(1,j+1,:)+lam(1,j-1,:))+fc*lam(1,j,:))
    enddo
    u(1,nyf,:)=Oomeph*u(1,nyf,:)+omeph*(fx*((lam(2,nyf,:)+lam(1,nyf,:))*u(2,nyf,:)+(lam(1,nyf,:)+rbufE(nyfd2,:))*rbufW(nyfd2,:)) &
      +fy*((rbufS(1,:)+lam(1,nyf,:))*rbufN(1,:)+(lam(1,nyf,:)+lam(1,nyf-1,:))*u(1,nyf-1,:)) &
      -rhs(1,nyf,:))/(fx*(lam(2,nyf,:)+rbufE(nyfd2,:))+fy*(rbufS(1,:)+lam(1,nyf-1,:))+fc*lam(1,nyf,:)) !.local NW corner
  endif

  call MPI_WAITALL(2,(/req(1),req(3)/),(/stat(:,1),stat(:,3)/),ierr)
  call MPI_WAITALL(2,(/req(5),req(7)/),(/stat(:,5),stat(:,7)/),ierr)
!.Send top row of array to North rank
  sbufN=u(2:nxf:2,nyf,:)
  call MPI_ISEND(sbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(1),ierr)
!.Receive u row from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(2),ierr)
!.Send top row of lam to North rank
  sbufS=lam(2:nxf:2,nyf,:)
  call MPI_ISEND(sbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(3),ierr)
!.Receive lam row from South rank, for lam(i,j-1)
  call MPI_IRECV(rbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(4),ierr)
!.Send first column of u to West rank
  sbufW=u(1,1:nyf-1:2,:)
  call MPI_ISEND(sbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(5),ierr)
!.Receive u column from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(6),ierr)
!.Send first column of lam to West rank
  sbufE=lam(1,1:nyf-1:2,:)
  call MPI_ISEND(sbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(7),ierr)
!.Receive lam column from East rank, for lam(i+1,j)
  call MPI_IRECV(rbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(8),ierr)
  do j=3,nyf-1,2
    do i=2,nxf-2,2
      u(i,j,:)=Oomeph*u(i,j,:)+omeph*(fx*((lam(i+1,j,:)+lam(i,j,:))*u(i+1,j,:)+(lam(i,j,:)+lam(i-1,j,:))*u(i-1,j,:)) &
        +fy*((lam(i,j+1,:)+lam(i,j,:))*u(i,j+1,:)+(lam(i,j,:)+lam(i,j-1,:))*u(i,j-1,:)) &
        -rhs(i,j,:))/(fx*(lam(i+1,j,:)+lam(i-1,j,:))+fy*(lam(i,j+1,:)+lam(i,j-1,:))+fc*lam(i,j,:))
    enddo
  enddo
!.Bottom, even i
  call MPI_WAITALL(2,(/req(2),req(4)/),(/stat(:,2),stat(:,4)/),ierr)
  do i=2,nxf-2,2
    u(i,1,:)=Oomeph*u(i,1,:)+omeph*(fx*((lam(i+1,1,:)+lam(i,1,:))*u(i+1,1,:)+(lam(i,1,:)+lam(i-1,1,:))*u(i-1,1,:)) &
      +fy*((lam(i,2,:)+lam(i,1,:))*u(i,2,:)+(lam(i,1,:)+rbufN(i/2,:))*rbufS(i/2,:)) &
      -rhs(i,1,:))/(fx*(lam(i+1,1,:)+lam(i-1,1,:))+fy*(lam(i,2,:)+rbufN(i/2,:))+fc*lam(i,1,:))
  enddo
!.Right boundary, odd j
    call MPI_WAITALL(2,(/req(6),req(8)/),(/stat(:,6),stat(:,8)/),ierr)
  if (iIDph(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf,1,:)=Oomeph*u(nxf,1,:)+omeph*(fx*0.25_dp*(2.0_dp*lam(nxf,1,:)*(6.0_dp*ur(1,:)+u(il,1,:)) &
      +(-2.0_dp*lam(nxf,1,:)+4.0_dp*lam(nxf-1,1,:))*u(nxf-1,1,:)) &
      +fy*((lam(nxf,2,:)+lam(nxf,1,:))*u(nxf,2,:)+(lam(nxf,1,:)+rbufN(nxfd2,:))*rbufS(nxfd2,:))-rhs(nxf,1,:)) &
      /(fx*(lam(nxf,1,:)+lam(nxf-1,1,:))+fy*(lam(nxf,2,:)+rbufN(nxfd2,:))+fc*lam(nxf,1,:)) !.global SE corner
    do j=3,nyf-1,2
      u(nxf,j,:)=Oomeph*u(nxf,j,:)+omeph*(fx*0.25_dp*(2.0_dp*lam(nxf,j,:)*(6.0_dp*ur(j,:)+u(il,j,:)) &
        +(-2.0_dp*lam(nxf,j,:)+4.0_dp*lam(nxf-1,j,:))*u(nxf-1,j,:)) &
        +fy*((lam(nxf,j+1,:)+lam(nxf,j,:))*u(nxf,j+1,:)+(lam(nxf,j,:)+lam(nxf,j-1,:))*u(nxf,j-1,:))-rhs(nxf,j,:)) &
        /(fx*(lam(nxf,j,:)+lam(nxf-1,j,:))+fy*(lam(nxf,j+1,:)+lam(nxf,j-1,:))+fc*lam(nxf,j,:))
    enddo
  else
    u(nxf,1,:)=Oomeph*u(nxf,1,:)+omeph*(fx*((rbufW(1,:)+lam(nxf,1,:))*rbufE(1,:)+(lam(nxf,1,:)+lam(nxf-1,1,:))*u(nxf-1,1,:)) &
      +fy*((lam(nxf,2,:)+lam(nxf,1,:))*u(nxf,2,:)+(lam(nxf,1,:)+rbufN(nxfd2,:))*rbufS(nxfd2,:)) &
      -rhs(nxf,1,:))/(fx*(rbufW(1,:)+lam(nxf-1,1,:))+fy*(lam(nxf,2,:)+rbufN(nxfd2,:))+fc*lam(nxf,1,:)) !.local SE corner
    do j=3,nyf-1,2
      u(nxf,j,:)=Oomeph*u(nxf,j,:)+omeph*(fx*((rbufW((j+1)/2,:)+lam(nxf,j,:))*rbufE((j+1)/2,:) &
        +(lam(nxf,j,:)+lam(nxf-1,j,:))*u(nxf-1,j,:))+fy*((lam(nxf,j+1,:)+lam(nxf,j,:))*u(nxf,j+1,:) &
        +(lam(nxf,j,:)+lam(nxf,j-1,:))*u(nxf,j-1,:))-rhs(nxf,j,:)) &
        /(fx*(rbufW((j+1)/2,:)+lam(nxf-1,j,:))+fy*(lam(nxf,j+1,:)+lam(nxf,j-1,:))+fc*lam(nxf,j,:))
    enddo
  endif

  call MPI_WAITALL(2,(/req(1),req(3)/),(/stat(:,1),stat(:,3)/),ierr)
  call MPI_WAITALL(2,(/req(5),req(7)/),(/stat(:,5),stat(:,7)/),ierr)
  end subroutine relaxphi
!*********************************************************
  function resphi(ng,u,lam,rhs,xBC,ur)
!.Returns minus the residual of the Poisson equation.
!       res = rhs - div(lam*grad(u))
!.using a 2nd order accurate finite difference stencil.
!.This is discretized such that the operator is self-adjoint.
!.INPUTS
!. integer ng: current grid
!. real(dp) dimension(:,:,:) u: current solution at grid ng
!. real(dp) dimension(:,:,:) lam: spatially dependent coefficient
!. real(dp) dimension(:,:,:) rhs: the right-hand side of Poisson equation
!. integer(2) xBC: boundary conditions along x
!. real(dp) dimension(:,:) ur: Dirichlet value for Right boundary
!.OUTPUTS
!. real(dp) dimension(:,:,:) resphi: negated residual
  implicit none
  include 'mpif.h'
  real(DP), dimension(:,:,:), intent(in) :: u, lam, rhs
  real(DP), intent(in), optional :: ur(:,:)
  integer(I4B), intent(in) :: ng,xBC(:)
  real(DP), allocatable :: resphi(:,:,:)
  real(DP), allocatable :: Lsu(:),loc1D(:) !Lsu(:,:),Rsu(:,:)
  integer(I4B) :: i,j,k,il,Nrank,Srank,Erank,Wrank
  real(DP) :: fx,fy,fc,Lqu,Lru,Rqu,Rru,Rsu
  real(DP), allocatable, dimension(:,:) :: sbufW,rbufE,sbufE,rbufW, &
                                           sbufN,rbufS,sbufS,rbufN
  real(DP), allocatable, dimension(:) :: SEu,SElam,NWu,NWlam
  integer(I4B) :: req(8),stat(MPI_STATUS_SIZE,8)
  integer(I4B) :: ierr

#IF (MGTIMER==1)
    call CPU_TIME(rest1)
#ENDIF

!.Neigboring ranks to communicate with
  Nrank=nborph(ng,1)
  Erank=nborph(ng,3)
  Srank=nborph(ng,5)
  Wrank=nborph(ng,7)

  allocate(resphi(nxf,nyf,nzL))
!.Boundary conditions are implemented with some amount of flexibility.
!.The right boundary for example, if u(nxf+1) is needed, it will be
!.written as
!.  u(nxf+1) = Rq*u(nxf) + Rr*u(nxf-1) + Rs
!.so that symmetry, extrapolation or Dirichlet conditions can be
!.easily imposed.
!.LHS BC of u
!  allocate(Lsu(nyf,nzL))
  allocate(Lsu(nzL),loc1D(nzL))
  if (xBC(1)==0) then
!...linear extrapolation w/o boundary value
    Lqu =  2.0_dp
    Lru = -1.0_dp
    Lsu =  0.0_dp
  elseif (xBC(1)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Lqu = -1.0_dp
    Lru =  0.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u(1,:,k))
    call MPI_ALLreduce(loc1D,Lsu,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMph(ng,1),ierr)
    Lsu =  2.0_dp*Lsu/dble(nyf*jprocsph(ng))
  else ! if (xBC(1)==1 .OR. xBC(1)==-1) then !.symmetry
    Lqu = dble(xBC(1))
    Lru = 0.0_dp
    Lsu = 0.0_dp
  endif
!.RHS BC of u
!  allocate(Rsu(nyf,nzL))
  if (xBC(2)==0) then
!...linear extrapolation w/o boundary value
    Rqu =  2.0_dp
    Rru = -1.0_dp
    Rsu =  0.0_dp
!...linear extrapolation w/ boundary value
!    Rqu = -1.0_dp
!    Rru =  0.0_dp
!    Rsu =  2.0_dp*ur
  else ! if (xBC(2)==1 .OR. xBC(2)==-1) then !.symmetry
    Rqu = dble(xBC(2))
    Rru = 0.0_dp
    Rsu = 0.0_dp
  endif
  fx = 0.5_dp*((dble(nxf*iprocsph(ng))/bLx)**2)         ! 1/(hx^2)
  fy = 0.5_dp*((dble(nyf*jprocsph(ng))/bLy)**2)         ! 1/(hy^2)
  fc = 2.0_dp*(fx+fy)
  if (nxf>2) then ! Used in one-sided stencil at RHS Dirichlet boundary
    il = nxf-2
  else
    il = nxf-1
  endif

!.When sending/receiving lam data, the names don't apply. Try to use
!.the opposite of what the name says: e.g. sbufE sends lam data West,
!.and rbufS receives lam data from North.
  allocate(sbufW(nyf,nzL),rbufE(nyf,nzL),sbufE(nyf,nzL),rbufW(nyf,nzL))
  allocate(sbufN(nxf,nzL),rbufS(nxf,nzL),sbufS(nxf,nzL),rbufN(nxf,nzL))
  allocate(SEu(nzL),SElam(nzL),NWu(nzL),NWlam(nzL))
!.Send top row of u to North rank
  sbufN=u(:,nyf,:)
  call MPI_ISEND(sbufN,nxf*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(1),ierr)
!.Send top row of lam to North rank
  sbufS=lam(:,nyf,:)
  call MPI_ISEND(sbufS,nxf*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(2),ierr)
!.Send right column of u to East rank
  sbufE=u(nxf,:,:)
  call MPI_ISEND(sbufE,nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(3),ierr)
!.Send right column of lam to East rank
  sbufW=lam(nxf,:,:)
  call MPI_ISEND(sbufW,nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(4),ierr)
!.Receive u row from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxf*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(5),ierr)
!.Receive lam row from South rank, for lam(i,j-1)
  call MPI_IRECV(rbufN,nxf*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(6),ierr)
!.Receive u column from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(7),ierr)
!.Receive lam column from West rank, for lam(i-1,j)
  call MPI_IRECV(rbufE,nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(8),ierr)
!.Interior points
  do j=2,nyf-1
    do i=2,nxf-1
      resphi(i,j,:)=-(fx*((lam(i+1,j,:)+lam(i,j,:))*u(i+1,j,:) &
        +(lam(i,j,:)+lam(i-1,j,:))*u(i-1,j,:))+fy*((lam(i,j+1,:)+lam(i,j,:))*u(i,j+1,:) &
        +(lam(i,j,:)+lam(i,j-1,:))*u(i,j-1,:))-(fx*(lam(i+1,j,:)+lam(i-1,j,:)) &
        +fy*(lam(i,j+1,:)+lam(i,j-1,:))+fc*lam(i,j,:))*u(i,j,:))+rhs(i,j,:)
    enddo
  enddo

!.Boundary points
  call MPI_WAITALL(2,(/req(5),req(6)/),(/stat(:,5),stat(:,6)/),ierr)
  do i=2,nxf-1 !.Bottom boundary
    resphi(i,1,:)=-(fx*((lam(i+1,1,:)+lam(i,1,:))*u(i+1,1,:) &
      +(lam(i,1,:)+lam(i-1,1,:))*u(i-1,1,:))+fy*((lam(i,2,:)+lam(i,1,:))*u(i,2,:) &
      +(lam(i,1,:)+rbufN(i,:))*rbufS(i,:))-(fx*(lam(i+1,1,:)+lam(i-1,1,:)) &
      +fy*(lam(i,2,:)+rbufN(i,:))+fc*lam(i,1,:))*u(i,1,:))+rhs(i,1,:)
  enddo
  SEu=rbufS(nxf,:)
  SElam=rbufN(nxf,:)
  call MPI_WAITALL(2,(/req(7),req(8)/),(/stat(:,7),stat(:,8)/),ierr)
  if (iIDph(ng) == 0) then !.Left most process (along x)
    resphi(1,1,:)=-(fx*((lam(2,1,:)+lam(1,1,:))*u(2,1,:) &  !.global SW corner
      +2.0_dp*lam(1,1,:)*(Lqu*u(1,1,:)+Lru*u(2,1,:)+Lsu))+fy*((lam(1,2,:)+lam(1,1,:))*u(1,2,:) &
      +(lam(1,1,:)+rbufN(1,:))*rbufS(1,:))-(fx*(lam(2,1,:)+lam(1,1,:)) &
      +fy*(lam(1,2,:)+rbufN(1,:))+fc*lam(1,1,:))*u(1,1,:))+rhs(1,1,:)
    do j=2,nyf-1 !.global left boundary
      resphi(1,j,:)=-(fx*((lam(2,j,:)+lam(1,j,:))*u(2,j,:) &
        +2.0_dp*lam(1,j,:)*(Lqu*u(1,j,:)+Lru*u(2,j,:)+Lsu))+fy*((lam(1,j+1,:)+lam(1,j,:))*u(1,j+1,:) &
        +(lam(1,j,:)+lam(1,j-1,:))*u(1,j-1,:))-(fx*(lam(2,j,:)+lam(1,j,:)) &
        +fy*(lam(1,j+1,:)+lam(1,j-1,:))+fc*lam(1,j,:))*u(1,j,:))+rhs(1,j,:)
    enddo
  else
    resphi(1,1,:)=-(fx*((lam(2,1,:)+lam(1,1,:))*u(2,1,:) & !.local SW corner
      +(lam(1,1,:)+rbufE(1,:))*rbufW(1,:))+fy*((lam(1,2,:)+lam(1,1,:))*u(1,2,:) &
      +(lam(1,1,:)+rbufN(1,:))*rbufS(1,:))-(fx*(lam(2,1,:)+rbufE(1,:)) &
      +fy*(lam(1,2,:)+rbufN(1,:))+fc*lam(1,1,:))*u(1,1,:))+rhs(1,1,:)
    do j=2,nyf-1 !.local left boundary
      resphi(1,j,:)=-(fx*((lam(2,j,:)+lam(1,j,:))*u(2,j,:) &
        +(lam(1,j,:)+rbufE(j,:))*rbufW(j,:))+fy*((lam(1,j+1,:)+lam(1,j,:))*u(1,j+1,:) &
        +(lam(1,j,:)+lam(1,j-1,:))*u(1,j-1,:))-(fx*(lam(2,j,:)+rbufE(j,:)) &
        +fy*(lam(1,j+1,:)+lam(1,j-1,:))+fc*lam(1,j,:))*u(1,j,:))+rhs(1,j,:)
    enddo
    NWu=rbufW(nyf,:)
    NWlam=rbufE(nyf,:)
  endif

  call MPI_WAITALL(2,(/req(1),req(2)/),(/stat(:,1),stat(:,2)/),ierr)
  call MPI_WAITALL(2,(/req(3),req(4)/),(/stat(:,3),stat(:,4)/),ierr)
!.Send bottom row of u to South rank
  sbufS=u(:,1,:)
  call MPI_ISEND(sbufS,nxf*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(1),ierr)
!.Send bottom row of lam to South rank
  sbufN=lam(:,1,:)
  call MPI_ISEND(sbufN,nxf*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(2),ierr)
!.Send left column of u to West rank
  sbufW=u(1,:,:)
  call MPI_ISEND(sbufW,nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(3),ierr)
!.Send left column of lam to West rank
  sbufE=lam(1,:,:)
  call MPI_ISEND(sbufE,nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(4),ierr)
!.Receive u row from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxf*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(5),ierr)
!.Receive lam row from North rank, for lam(i,j+1)
  call MPI_IRECV(rbufS,nxf*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(6),ierr)
!.Receive u column from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(7),ierr)
!.Receive lam column from East rank, for lam(i+1,j)
  call MPI_IRECV(rbufW,nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(8),ierr)
  call MPI_WAITALL(2,(/req(5),req(6)/),(/stat(:,5),stat(:,6)/),ierr)
  do i=2,nxf-1 !.Top boundary
    resphi(i,nyf,:)=-(fx*((lam(i+1,nyf,:)+lam(i,nyf,:))*u(i+1,nyf,:) &
      +(lam(i,nyf,:)+lam(i-1,nyf,:))*u(i-1,nyf,:))+fy*((rbufS(i,:)+lam(i,nyf,:))*rbufN(i,:) &
      +(lam(i,nyf,:)+lam(i,nyf-1,:))*u(i,nyf-1,:))-(fx*(lam(i+1,nyf,:)+lam(i-1,nyf,:)) &
      +fy*(rbufS(i,:)+lam(i,nyf-1,:))+fc*lam(i,nyf,:))*u(i,nyf,:))+rhs(i,nyf,:)
  enddo
  if (iIDph(ng)==0) then !.Left most process (along x)
!...global NW corner
    resphi(1,nyf,:)=-(fx*((lam(2,nyf,:)+lam(1,nyf,:))*u(2,nyf,:) &
      +2.0_dp*lam(1,nyf,:)*(Lqu*u(1,nyf,:)+Lru*u(2,nyf,:)+Lsu))+fy*((rbufS(1,:)+lam(1,nyf,:))*rbufN(1,:) &
      +(lam(1,nyf,:)+lam(1,nyf-1,:))*u(1,nyf-1,:))-(fx*(lam(2,nyf,:)+lam(1,nyf,:)) &
      +fy*(rbufS(1,:)+lam(1,nyf-1,:))+fc*lam(1,nyf,:))*u(1,nyf,:))+rhs(1,nyf,:)
  else !.local NW corner
    resphi(1,nyf,:)=-(fx*((lam(2,nyf,:)+lam(1,nyf,:))*u(2,nyf,:) &
      +(lam(1,nyf,:)+NWlam)*NWu)+fy*((rbufS(1,:)+lam(1,nyf,:))*rbufN(1,:) &
      +(lam(1,nyf,:)+lam(1,nyf-1,:))*u(1,nyf-1,:))-(fx*(lam(2,nyf,:)+NWlam) &
      +fy*(rbufS(1,:)+lam(1,nyf-1,:))+fc*lam(1,nyf,:))*u(1,nyf,:))+rhs(1,nyf,:)
  endif

  call MPI_WAITALL(2,(/req(7),req(8)/),(/stat(:,7),stat(:,8)/),ierr)
  if (iIDph(ng) == iprocsm1) then !.Right most process (along x)
    resphi(nxf,1,:)=-(fx*0.25_dp*(2.0_dp*lam(nxf,1,:)*(6.0_dp*ur(1,:)+u(il,1,:)) &
      +(-2.0_dp*lam(nxf,1,:)+4.0_dp*lam(nxf-1,1,:))*u(nxf-1,1,:)) &
      +fy*((lam(nxf,2,:)+lam(nxf,1,:))*u(nxf,2,:)+(lam(nxf,1,:)+SElam)*SEu) &
      -(fx*(lam(nxf,1,:)+lam(nxf-1,1,:))+fy*(lam(nxf,2,:)+SElam) &
      +fc*lam(nxf,1,:))*u(nxf,1,:))+rhs(nxf,1,:) !.global SE corner
    do j=2,nyf-1 !.global right boundary
      resphi(nxf,j,:)=-(fx*0.25_dp*(2.0_dp*lam(nxf,j,:)*(6.0_dp*ur(j,:)+u(il,j,:)) &
        +(-2.0_dp*lam(nxf,j,:)+4.0_dp*lam(nxf-1,j,:))*u(nxf-1,j,:)) &
        +fy*((lam(nxf,j+1,:)+lam(nxf,j,:))*u(nxf,j+1,:)+(lam(nxf,j,:)+lam(nxf,j-1,:))*u(nxf,j-1,:)) &
        -(fx*(lam(nxf,j,:)+lam(nxf-1,j,:))+fy*(lam(nxf,j+1,:)+lam(nxf,j-1,:)) &
        +fc*lam(nxf,j,:))*u(nxf,j,:))+rhs(nxf,j,:)
    enddo
    resphi(nxf,nyf,:)=-(fx*0.25_dp*(2.0_dp*lam(nxf,nyf,:)*(6.0_dp*ur(nyf,:)+u(il,nyf,:)) &
      +(-2.0_dp*lam(nxf,nyf,:)+4.0_dp*lam(nxf-1,nyf,:))*u(nxf-1,nyf,:)) &
      +fy*((rbufS(nxf,:)+lam(nxf,nyf,:))*rbufN(nxf,:)+(lam(nxf,nyf,:)+lam(nxf,nyf-1,:))*u(nxf,nyf-1,:)) &
      -(fx*(lam(nxf,nyf,:)+lam(nxf-1,nyf,:))+fy*(rbufS(nxf,:)+lam(nxf,nyf-1,:)) &
      +fc*lam(nxf,nyf,:))*u(nxf,nyf,:))+rhs(nxf,nyf,:) !.global NE corner
  else
    resphi(nxf,1,:)=-(fx*((rbufW(1,:)+lam(nxf,1,:))*rbufE(1,:) & !.local SE corner
      +(lam(nxf,1,:)+lam(nxf-1,1,:))*u(nxf-1,1,:))+fy*((lam(nxf,2,:)+lam(nxf,1,:))*u(nxf,2,:) &
      +(lam(nxf,1,:)+SElam)*SEu)-(fx*(rbufW(1,:)+lam(nxf-1,1,:)) &
      +fy*(lam(nxf,2,:)+SElam)+fc*lam(nxf,1,:))*u(nxf,1,:))+rhs(nxf,1,:)
    do j=2,nyf-1 !.local right boundary
      resphi(nxf,j,:)=-(fx*((rbufW(j,:)+lam(nxf,j,:))*rbufE(j,:) &
        +(lam(nxf,j,:)+lam(nxf-1,j,:))*u(nxf-1,j,:))+fy*((lam(nxf,j+1,:)+lam(nxf,j,:))*u(nxf,j+1,:) &
        +(lam(nxf,j,:)+lam(nxf,j-1,:))*u(nxf,j-1,:))-(fx*(rbufW(j,:)+lam(nxf-1,j,:)) &
        +fy*(lam(nxf,j+1,:)+lam(nxf,j-1,:))+fc*lam(nxf,j,:))*u(nxf,j,:))+rhs(nxf,j,:)
    enddo
    resphi(nxf,nyf,:)=-(fx*((rbufW(nyf,:)+lam(nxf,nyf,:))*rbufE(nyf,:) &
      +(lam(nxf,nyf,:)+lam(nxf-1,nyf,:))*u(nxf-1,nyf,:))+fy*((rbufS(nxf,:)+lam(nxf,nyf,:))*rbufN(nxf,:) &
      +(lam(nxf,nyf,:)+lam(nxf,nyf-1,:))*u(nxf,nyf-1,:))-(fx*(rbufW(nyf,:)+lam(nxf-1,nyf,:)) &
      +fy*(rbufS(nxf,:)+lam(nxf,nyf-1,:))+fc*lam(nxf,nyf,:))*u(nxf,nyf,:))+rhs(nxf,nyf,:) !.global NE corner
  endif
  call MPI_WAITALL(2,(/req(1),req(2)/),(/stat(:,1),stat(:,2)/),ierr)
  call MPI_WAITALL(2,(/req(3),req(4)/),(/stat(:,3),stat(:,4)/),ierr)

#IF (MGTIMER==1)
    call CPU_TIME(rest2)
    restph=restph+rest2-rest1
#ENDIF
  end function resphi
!*********************************************************
  recursive subroutine Gcycpsi(ng,u,rhs,xBC)
!.Gamma cycles to solve modified Helmholtz equation and obtain psi
!                   u-ksq*L(u)=rho
!.where u=psi, rho=rpsi and ksq=(electron skin depth)^2
!.INPUTS
!. integer ng: current grid
!. real(dp) dimension(:,:,:) u: current solution at grid ng
!. real(dp) dimension(:,:,:) rhs: the right-hand side of mod Helmholtz equation
!. integer(2) xBC: boundary conditions along x
!.OUTPUTS
!. real(dp) dimension(:,:,:) u: updated solution at grid ng
  implicit none
  integer(I4B), intent(in) :: ng
  real(DP), intent(inout) :: u(:,:,:)
  real(DP), intent(in) :: rhs(:,:,:)
  integer(I4B), intent(in) :: xBC(:)
  integer(I4B) :: jpost, jpre, i
  real(DP), allocatable :: res(:,:,:),v(:,:,:)

  iprocsm1=iprocsps(ng)-1 !.ID of last rank in x COMM
  jprocsm1=jprocsps(ng)-1 !.ID of last rank in y COMM

  if (ng==1) then        !.Bottom of gamma cycle: Solve on coarsest grid
#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpost=1,NU2ps
      call relaxpsi(ng,u,rhs,xBC)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    postXtps=postXtps+count2-count1
#ENDIF
  else                  !.On downward stroke towards coarser grids
#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpre=1,NU1ps      !.Pre-smoothing
      call relaxpsi(ng,u,rhs,xBC)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    preXtps=preXtps+count2-count1
#ENDIF

!...Number of grid points in coarse grid
    nxc=nxCps(ng)
    nyc=nyCps(ng)
!...Number of grid points in coarse grid after subdomain unification
    nxu=nxFps(ng-1)
    nyu=nyFps(ng-1)
    allocate(res(nxu,nyu,nzL),v(nxu,nyu,nzL))

#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
!...Restricted residual is the next RHS 
    res=rstrct(respsi(ng,u,rhs,xBC),wdcps(ng),xBC, &
      acups(ng-1,:),urnkps(ng-1)%a,nurnkps(ng-1,:))
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    rNrtps=rNrtps+count2-count1
#ENDIF

!...Next fine grid is the coarse, or unified subdomain, grid
    nxf=nxu
    nyf=nyu

    v=0.0_dp !.Zero for initial guess in next relaxation

    if (acups(ng-1,1)) then
!.....Only active processes proceed to next coarser grid
      do i=1,GAMps
        call Gcycpsi(ng-1,v,res,xBC)  !.Recursive call for coarse grid correction
      enddo
    endif

!...Redefine quantities that were modified in coarser grids
    iprocsm1=iprocsps(ng)-1
    jprocsm1=jprocsps(ng)-1
    nxc=nxCps(ng)
    nyc=nyCps(ng)
    nxf=nxFps(ng)
    nyf=nyFps(ng)

#IF (MGTIMER==1)
    call CPU_TIME(prolt1)
#ENDIF
!...Upward stroke -> finer grid: prolong coarse grid error + correct solution
    u=u+prolng(v,wdcps(ng),xBC,COMMps(ng,:),iIDps(ng),jprocsps(ng), &
      nborps(ng,:),acups(ng-1,:),urnkps(ng-1)%a,nurnkps(ng-1,:))
#IF (MGTIMER==1)
    call CPU_TIME(prolt2)
    proltps=proltps+prolt2-prolt1
#ENDIF

#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpost=1,NU2ps      !.Post-smoothing
      call relaxpsi(ng,u,rhs,xBC)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    postXtps=postXtps+count2-count1
#ENDIF
  endif

  end subroutine Gcycpsi
!*********************************************************
  subroutine relaxpsi(ng,u,rhs,xBC)
!.Red-black Gauss-Seidel relaxation of
!                   u-ksq*L(u)=rho
!.using a 2nd order accurate finite difference stencil.
!.INPUTS
!. integer ng: current grid
!. real(dp) dimension(:,:,:) u: current solution at grid ng
!. real(dp) dimension(:,:,:) rhs: the right-hand side of mod Helmholtz equation
!. integer(2) xBC: boundary conditions along x
!.OUTPUTS
!. real(dp) dimension(:,:,:) u: updated solution at grid ng
  implicit none
  include 'mpif.h'
  real(DP), intent(inout) :: u(:,:,:)
  real(DP), intent(in) :: rhs(:,:,:)
  integer(I4B), intent(in) :: ng,xBC(:)
  integer(I4B) :: i,j,k,nxfd2,nyfd2,Nrank,Srank,Erank,Wrank
  real(DP) :: fx,fy,fc,Lq,Rq,delta,deltaL,deltaR
  real(DP), allocatable :: Ls(:),Rs(:),loc1D(:)
  real(DP), allocatable, dimension(:,:) :: sbufW,rbufE,sbufE,rbufW, &
                                           sbufN,rbufS,sbufS,rbufN
  integer(I4B) :: req(4),stat(MPI_STATUS_SIZE,4),ierr
  integer(I4B) :: req2(4),stat2(MPI_STATUS_SIZE,4)

  nxfd2=nxf/2
  nyfd2=nyf/2
!.Neigboring ranks to communicate with
  Nrank=nborps(ng,1)
  Erank=nborps(ng,3)
  Srank=nborps(ng,5)
  Wrank=nborps(ng,7)

!.Boundary conditions are implemented with some amount of flexibility.
!.The right boundary for example, if u(nxf+1) is needed, it will be
!.written as
!.  u(nxf+1) = Rq*u(nxf) + Rr*u(nxf-1) + Rs
!.so that symmetry, extrapolation or Dirichlet conditions can be
!.easily imposed.
!......... BCs of u ................................!
  allocate(Ls(nzL),Rs(nzL),loc1D(nzL))
!.LHS BC of u
  if (xBC(1)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Lq = -1.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u(1,:,k))
    call MPI_ALLreduce(loc1D,ls,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMps(ng,1),ierr)
    ls =  2.0_dp*ls/dble(nyf*jprocsps(ng))
  else ! if (xBC(1)==1 .OR. xBC(1)==-1) then !.symmetry
    Lq = dble(xBC(1))
    ls = 0.0_dp
  endif
!.RHS BC of u
  if (xBC(2)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Rq = -1.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u(nxf,:,k))
    call MPI_ALLreduce(loc1D,rs,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMps(ng,1),ierr)
    rs =  2.0_dp*rs/dble(nyf*jprocsps(ng))
  else ! if (xBC(2)==1 .OR. xBC(2)==-1) then !.symmetry
    Rq = dble(xBC(2))
    rs = 0.0_dp
  endif
!........................................................!
  fx = -ksq*(dble(nxf)*iprocsps(ng)/bLx)**2
  fy = -ksq*(dble(nyf)*jprocsps(ng)/bLy)**2
  fc =  2.0_dp*(fx+fy)
!  delta  = omeps/(1.0_dp-fc)
!  deltaL = omeps/(1.0_dp-(fc-fx*Lq))
!  deltaR = omeps/(1.0_dp-(fc-fx*Rq))
  delta  = 1.0_dp/(1.0_dp-fc)
  deltaL = 1.0_dp/(1.0_dp-(fc-fx*Lq))
  deltaR = 1.0_dp/(1.0_dp-(fc-fx*Rq))

!.Send buffer names indicate direction data is sent to
!.and receive buffers direction data is received from.
!.CAREFUL!!!!
  allocate(sbufS(nxfd2,nzL),rbufN(nxfd2,nzL),sbufW(nyfd2,nzL),rbufE(nyfd2,nzL))
!.First do the even-even and odd-odd squares of the grid,
!.i.e., the red squares of the checker-board
!.Comunicate black grid points needed
!.Send bottom row of u to South rank
  sbufS=u(2:nxf:2,1,:)
  call MPI_ISEND(sbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMps(ng,1),req(1),ierr)
!.Receive u row from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMps(ng,1),req(2),ierr)
!.Send first column of u to West rank
  sbufW=u(1,2:nyf:2,:)
  call MPI_ISEND(sbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMps(ng,2),req(3),ierr)
!.Receive u column from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMps(ng,2),req(4),ierr)
  do j=2,nyf-2,2
    do i=2,nxf-2,2
!      u(i,j,:)=Oomeps*u(i,j,:)+(rhs(i,j,:)-fx*(u(i+1,j,:)+u(i-1,j,:)) &
      u(i,j,:)=(rhs(i,j,:)-fx*(u(i+1,j,:)+u(i-1,j,:)) &
        -fy*(u(i,j+1,:)+u(i,j-1,:)))*delta
    enddo
  enddo
!.Top, even i
  call MPI_WAIT(req(2),stat(:,2),ierr)
  do i=2,nxf-2,2
!    u(i,nyf,:)=Oomeps*u(i,nyf,:)+(rhs(i,nyf,:)-fx*(u(i+1,nyf,:)+u(i-1,nyf,:)) &
    u(i,nyf,:)=(rhs(i,nyf,:)-fx*(u(i+1,nyf,:)+u(i-1,nyf,:)) &
      -fy*(rbufN(i/2,:)+u(i,nyf-1,:)))*delta
  enddo
!.Right boundary, even j
  call MPI_WAIT(req(4),stat(:,4),ierr)
  if (iIDps(ng) .eq. iprocsm1) then !.Right most process (along x)
    do j=2,nyf-2,2
!      u(nxf,j,:)=Oomeps*u(nxf,j,:)+(rhs(nxf,j,:)-fx*(rs+u(nxf-1,j,:)) &
      u(nxf,j,:)=(rhs(nxf,j,:)-fx*(rs+u(nxf-1,j,:)) &
        -fy*(u(nxf,j+1,:)+u(nxf,j-1,:)))*deltaR
    enddo
!    u(nxf,nyf,:)=Oomeps*u(nxf,nyf,:)+(rhs(nxf,nyf,:)-fx*(rs+u(nxf-1,nyf,:)) &
    u(nxf,nyf,:)=(rhs(nxf,nyf,:)-fx*(rs+u(nxf-1,nyf,:)) &
      -fy*(rbufN(nxfd2,:)+u(nxf,nyf-1,:)))*deltaR !.global NE corner
  else
    do j=2,nyf-2,2
!      u(nxf,j,:)=Oomeps*u(nxf,j,:)+(rhs(nxf,j,:)-fx*(rbufE(j/2,:)+u(nxf-1,j,:)) &
      u(nxf,j,:)=(rhs(nxf,j,:)-fx*(rbufE(j/2,:)+u(nxf-1,j,:)) &
        -fy*(u(nxf,j+1,:)+u(nxf,j-1,:)))*delta
    enddo
    u(nxf,nyf,:)=(rhs(nxf,nyf,:)-fx*(rbufE(nyfd2,:)+u(nxf-1,nyf,:)) &
      -fy*(rbufN(nxfd2,:)+u(nxf,nyf-1,:)))*delta !.local NE corner
  endif

  call MPI_WAITALL(2,(/req(1),req(3)/),(/stat(:,1),stat(:,3)/),ierr)
  allocate(sbufN(nxfd2,nzL),rbufS(nxfd2,nzL),sbufE(nyfd2,nzL),rbufW(nyfd2,nzL))
!.Send top row of array to North rank
  sbufN=u(1:nxf-1:2,nyf,:)
  call MPI_ISEND(sbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMps(ng,1),req2(1),ierr)
!.Receive u row from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMps(ng,1),req2(2),ierr)
!.Send last column of array to East rank
  sbufE=u(nxf,1:nyf-1:2,:)
  call MPI_ISEND(sbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMps(ng,2),req2(3),ierr)
!.Receive u column from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMps(ng,2),req2(4),ierr)
  do j=3,nyf-1,2
    do i=3,nxf-1,2
!      u(i,j,:)=Oomeps*u(i,j,:)+(rhs(i,j,:)-fx*(u(i+1,j,:)+u(i-1,j,:)) &
      u(i,j,:)=(rhs(i,j,:)-fx*(u(i+1,j,:)+u(i-1,j,:)) &
        -fy*(u(i,j+1,:)+u(i,j-1,:)))*delta
    enddo
  enddo
!.Bottom, odd i
  call MPI_WAIT(req2(2),stat2(:,2),ierr)
  do i=3,nxf-1,2
!    u(i,1,:)=Oomeps*u(i,1,:)+(rhs(i,1,:)-fx*(u(i+1,1,:)+u(i-1,1,:)) &
    u(i,1,:)=(rhs(i,1,:)-fx*(u(i+1,1,:)+u(i-1,1,:)) &
      -fy*(u(i,2,:)+rbufS((i+1)/2,:)))*delta
  enddo
!.Left boundary, odd j
  call MPI_WAIT(req2(4),stat2(:,4),ierr)
  if (iIDps(ng)==0) then !.Left most process (along x)
!    u(1,1,:)=Oomeps*u(1,1,:)+(rhs(1,1,:)-fx*(u(2,1,:)+Ls) &
    u(1,1,:)=(rhs(1,1,:)-fx*(u(2,1,:)+Ls) &
      -fy*(u(1,2,:)+rbufS(1,:)))*deltaL !.global SW corner
    do j=3,nyf-1,2
!      u(1,j,:)=Oomeps*u(1,j,:)+(rhs(1,j,:)-fx*(u(2,j,:)+Ls) &
      u(1,j,:)=(rhs(1,j,:)-fx*(u(2,j,:)+Ls) &
        -fy*(u(1,j+1,:)+u(1,j-1,:)))*deltaL
    enddo
  else
!    u(1,1,:)=Oomeps*u(1,1,:)+(rhs(1,1,:)-fx*(u(2,1,:)+rbufW(1,:)) &
    u(1,1,:)=(rhs(1,1,:)-fx*(u(2,1,:)+rbufW(1,:)) &
      -fy*(u(1,2,:)+rbufS(1,:)))*delta
    do j=3,nyf-1,2
!      u(1,j,:)=Oomeps*u(1,j,:)+(rhs(1,j,:)-fx*(u(2,j,:)+rbufW((j+1)/2,:)) &
      u(1,j,:)=(rhs(1,j,:)-fx*(u(2,j,:)+rbufW((j+1)/2,:)) &
        -fy*(u(1,j+1,:)+u(1,j-1,:)))*delta !.local SW corner
    enddo
  endif
!  call MPI_WAITALL(4,(/req(1),req(3),req2(1),req2(3)/), &
!    (/stat(:,1),stat(:,3),stat2(:,1),stat2(:,3)/),ierr)
  call MPI_WAITALL(2,(/req2(1),req2(3)/), &
    (/stat2(:,1),stat2(:,3)/),ierr)

!.Now do even-odd and odd-even squares of the grid, i.e, the black
!.squares of the checker-board
!.Comunicate red grid points needed
!.Send bottom row of u to South rank
  sbufS=u(1:nxf-1:2,1,:)
  call MPI_ISEND(sbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMps(ng,1),req(1),ierr)
!.Receive u row from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMps(ng,1),req(2),ierr)
!.Send last column of array to East rank
  sbufE=u(nxf,2:nyf:2,:)
  call MPI_ISEND(sbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMps(ng,2),req(3),ierr)
!.Receive u column from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMps(ng,2),req(4),ierr)
  do j=2,nyf-2,2
    do i=3,nxf-1,2
!      u(i,j,:)=Oomeps*u(i,j,:)+(rhs(i,j,:)-fx*(u(i+1,j,:)+u(i-1,j,:)) &
      u(i,j,:)=(rhs(i,j,:)-fx*(u(i+1,j,:)+u(i-1,j,:)) &
        -fy*(u(i,j+1,:)+u(i,j-1,:)))*delta
    enddo
  enddo
!.Top, odd i
  call MPI_WAIT(req(2),stat(:,2),ierr)
  do i=3,nxf-1,2
!    u(i,nyf,:)=Oomeps*u(i,nyf,:)+(rhs(i,nyf,:)-fx*(u(i+1,nyf,:)+u(i-1,nyf,:)) &
    u(i,nyf,:)=(rhs(i,nyf,:)-fx*(u(i+1,nyf,:)+u(i-1,nyf,:)) &
      -fy*(rbufN((i+1)/2,:)+u(i,nyf-1,:)))*delta
  enddo
!.Left boundary, even j
  call MPI_WAIT(req(4),stat(:,4),ierr)
  if (iIDps(ng)==0) then !.Left most process (along x)
    do j=2,nyf-2,2
!      u(1,j,:)=Oomeps*u(1,j,:)+(rhs(1,j,:)-fx*(u(2,j,:)+Ls) &
      u(1,j,:)=(rhs(1,j,:)-fx*(u(2,j,:)+Ls) &
        -fy*(u(1,j+1,:)+u(1,j-1,:)))*deltaL
    enddo
!    u(1,nyf,:)=Oomeps*u(1,nyf,:)+(rhs(1,nyf,:)-fx*(u(2,nyf,:)+Ls) &
    u(1,nyf,:)=(rhs(1,nyf,:)-fx*(u(2,nyf,:)+Ls) &
      -fy*(rbufN(1,:)+u(1,nyf-1,:)))*deltaL !.global NW corner
  else
    do j=2,nyf-2,2
!      u(1,j,:)=Oomeps*u(1,j,:)+(rhs(1,j,:)-fx*(u(2,j,:)+rbufW(j/2,:)) &
      u(1,j,:)=(rhs(1,j,:)-fx*(u(2,j,:)+rbufW(j/2,:)) &
        -fy*(u(1,j+1,:)+u(1,j-1,:)))*delta
    enddo
!    u(1,nyf,:)=Oomeps*u(1,nyf,:)+(rhs(1,nyf,:)-fx*(u(2,nyf,:)+rbufW(nyfd2,:)) &
    u(1,nyf,:)=(rhs(1,nyf,:)-fx*(u(2,nyf,:)+rbufW(nyfd2,:)) &
      -fy*(rbufN(1,:)+u(1,nyf-1,:)))*delta !.local NW corner
  endif
!  call MPI_WAITALL(2,(/req(1),req(3)/),(/stat(:,1),stat(:,3)/),ierr)

!.Send top row of array to North rank
  sbufN=u(2:nxf:2,nyf,:)
  call MPI_ISEND(sbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMps(ng,1),req2(1),ierr)
!.Receive u row from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMps(ng,1),req2(2),ierr)
!.Send first column of u to West rank
  sbufW=u(1,1:nyf-1:2,:)
  call MPI_ISEND(sbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMps(ng,2),req2(3),ierr)
!.Receive u column from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMps(ng,2),req2(4),ierr)
  do j=3,nyf-1,2
    do i=2,nxf-2,2
!      u(i,j,:)=Oomeps*u(i,j,:)+(rhs(i,j,:)-fx*(u(i+1,j,:)+u(i-1,j,:)) &
      u(i,j,:)=(rhs(i,j,:)-fx*(u(i+1,j,:)+u(i-1,j,:)) &
        -fy*(u(i,j+1,:)+u(i,j-1,:)))*delta
    enddo
  enddo
!.Bottom, even i
  call MPI_WAIT(req2(2),stat2(:,2),ierr)
  do i=2,nxf-2,2
!    u(i,1,:)=Oomeps*u(i,1,:)+(rhs(i,1,:)-fx*(u(i+1,1,:)+u(i-1,1,:)) &
    u(i,1,:)=(rhs(i,1,:)-fx*(u(i+1,1,:)+u(i-1,1,:)) &
      -fy*(u(i,2,:)+rbufS(i/2,:)))*delta
  enddo
!.Right boundary, odd j
  call MPI_WAIT(req2(4),stat2(:,4),ierr)
  if (iIDps(ng) .eq. iprocsm1) then !.Right most process (along x)
!    u(nxf,1,:)=Oomeps*u(nxf,1,:)+(rhs(nxf,1,:)-fx*(rs+u(nxf-1,1,:)) &
    u(nxf,1,:)=(rhs(nxf,1,:)-fx*(rs+u(nxf-1,1,:)) &
      -fy*(u(nxf,2,:)+rbufS(nxfd2,:)))*deltaR !.global SE corner
    do j=3,nyf-1,2
!      u(nxf,j,:)=Oomeps*u(nxf,j,:)+(rhs(nxf,j,:)-fx*(rs+u(nxf-1,j,:)) &
      u(nxf,j,:)=(rhs(nxf,j,:)-fx*(rs+u(nxf-1,j,:)) &
        -fy*(u(nxf,j+1,:)+u(nxf,j-1,:)))*deltaR
    enddo
  else
!    u(nxf,1,:)=Oomeps*u(nxf,1,:)+(rhs(nxf,1,:)-fx*(rbufE(1,:)+u(nxf-1,1,:)) &
    u(nxf,1,:)=(rhs(nxf,1,:)-fx*(rbufE(1,:)+u(nxf-1,1,:)) &
      -fy*(u(nxf,2,:)+rbufS(nxfd2,:)))*delta
    do j=3,nyf-1,2
!      u(nxf,j,:)=Oomeps*u(nxf,j,:)+(rhs(nxf,j,:)-fx*(rbufE((j+1)/2,:)+u(nxf-1,j,:)) &
      u(nxf,j,:)=(rhs(nxf,j,:)-fx*(rbufE((j+1)/2,:)+u(nxf-1,j,:)) &
        -fy*(u(nxf,j+1,:)+u(nxf,j-1,:)))*delta !.local SE corner
    enddo
  endif
  call MPI_WAITALL(4,(/req(1),req(3),req2(1),req2(3)/), &
    (/stat(:,1),stat(:,3),stat2(:,1),stat2(:,3)/),ierr)

  end subroutine relaxpsi
!*********************************************************
  function respsi(ng,u,rhs,xBC)
!.Returns minus the residual of the modified Helmholtz equation.
!                   res = rho - (u-ksq*L(u))
!.INPUTS
!. integer ng: current grid
!. real(dp) dimension(:,:,:) u: current solution at grid ng
!. real(dp) dimension(:,:,:) rhs: the right-hand side of mod Helmholtz equation
!. integer(2) xBC: boundary conditions along x
!.OUTPUTS
!. real(dp) dimension(:,:,:) u: updated solution at grid ng
  implicit none
  include 'mpif.h'
  real(DP), allocatable :: respsi(:,:,:)
  real(DP), intent(in) :: u(:,:,:),rhs(:,:,:)
  integer(I4B), intent(in) :: ng,xBC(:)
  integer(I4B) :: i,j,k,Nrank,Srank,Erank,Wrank
  real(DP) :: fx,fy,fc,Lq,Rq
  real(DP), allocatable :: Ls(:),Rs(:),loc1D(:)
  real(DP), allocatable, dimension(:,:) :: sbufW,rbufE,sbufE,rbufW, &
                                           sbufN,rbufS,sbufS,rbufN
  integer(I4B) :: req(8),stat(MPI_STATUS_SIZE,8)
  integer(I4B) :: ierr

#IF (MGTIMER==1)
    call CPU_TIME(rest1)
#ENDIF

!.Neigboring ranks to communicate with
  Nrank=nborps(ng,1)
  Erank=nborps(ng,3)
  Srank=nborps(ng,5)
  Wrank=nborps(ng,7)

  allocate(respsi(nxf,nyf,nzL),Ls(nzL),Rs(nzL),loc1D(nzL))
!.Boundary conditions are implemented with some amount of flexibility.
!.The right boundary for example, if u(nxf+1) is needed, it will be
!.written as
!.  u(nxf+1) = Rq*u(nxf) + Rr*u(nxf-1) + Rs
!.so that symmetry, extrapolation or Dirichlet conditions can be
!.easily imposed.
!.LHS BC of u
  if (xBC(1)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Lq = -1.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u(1,:,k))
    call MPI_ALLreduce(loc1D,ls,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMps(ng,1),ierr)
    ls =  2.0_dp*ls/dble(nyf*jprocsps(ng))
  else ! if (xBC(1)==1 .OR. xBC(1)==-1) then !.symmetry
    Lq = dble(xBC(1))
    ls = 0.0_dp
  endif
!.RHS BC of u
  if (xBC(2)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Rq = -1.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u(nxf,:,k))
    call MPI_ALLreduce(loc1D,rs,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMps(ng,1),ierr)
    rs =  2.0_dp*rs/dble(nyf*jprocsps(ng))
  else ! if (xBC(2)==1 .OR. xBC(2)==-1) then !.symmetry
    Rq = dble(xBC(2))
    rs = 0.0_dp
  endif
  fx = -ksq*(dble(nxf)*iprocsps(ng)/bLx)**2         !.1/(hx^2)
  fy = -ksq*(dble(nyf)*jprocsps(ng)/bLy)**2         !.1/(hy^2)
  fc =  2.0_dp*(fx+fy)

  allocate(sbufN(nxf,nzL),rbufS(nxf,nzL),sbufE(nyf,nzL),rbufW(nyf,nzL))
  allocate(sbufS(nxf,nzL),rbufN(nxf,nzL),sbufW(nyf,nzL),rbufE(nyf,nzL))
!.Send top row of u to North rank
  sbufN=u(:,nyf,:)
  call MPI_ISEND(sbufN,nxf*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMps(ng,1),req(1),ierr)
!.Send bottom row of u to South rank
  sbufS=u(:,1,:)
  call MPI_ISEND(sbufS,nxf*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMps(ng,1),req(2),ierr)
!.Send right column of u to East rank
  sbufE=u(nxf,:,:)
  call MPI_ISEND(sbufE,nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMps(ng,2),req(3),ierr)
!.Send left column of u to West rank
  sbufW=u(1,:,:)
  call MPI_ISEND(sbufW,nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMps(ng,2),req(4),ierr)
!.Receive u row from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxf*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMps(ng,1),req(5),ierr)
!.Receive u row from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxf*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMps(ng,1),req(6),ierr)
!.Receive u column from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMps(ng,2),req(7),ierr)
!.Receive u column from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMps(ng,2),req(8),ierr)
!.Interior points
  do j=2,nyf-1
    do i=2,nxf-1
      respsi(i,j,:)=rhs(i,j,:)-(fx*(u(i+1,j,:)+u(i-1,j,:)) &
        +fy*(u(i,j+1,:)+u(i,j-1,:))+(1.0_dp-fc)*u(i,j,:))
    enddo
  enddo

!.Boundary points
  call MPI_WAIT(req(5),stat(:,5),ierr)
  do i=2,nxf-1 !.bottom boundary
    respsi(i,1,:)=rhs(i,1,:)-(fx*(u(i+1,1,:)+u(i-1,1,:)) &
      +fy*(u(i,2,:)+rbufS(i,:))+(1.0_dp-fc)*u(i,1,:))
  enddo
  call MPI_WAIT(req(6),stat(:,6),ierr)
  do i=2,nxf-1 !.top boundary
    respsi(i,nyf,:)=rhs(i,nyf,:)-(fx*(u(i+1,nyf,:)+u(i-1,nyf,:)) &
      +fy*(rbufN(i,:)+u(i,nyf-1,:))+(1.0_dp-fc)*u(i,nyf,:))
  enddo
!.Left boundary
  call MPI_WAIT(req(7),stat(:,7),ierr)
  if (iIDps(ng) == 0) then !.Left most process (along x)
    respsi(1,1,:)=rhs(1,1,:)-(fx*(u(2,1,:)+Lq*u(1,1,:)+Ls) &
      +fy*(u(1,2,:)+rbufS(1,:))+(1.0_dp-fc)*u(1,1,:)) !.global SW corner
    do j=2,nyf-1
      respsi(1,j,:)=rhs(1,j,:)-(fx*(u(2,j,:)+Lq*u(1,j,:)+Ls) &
        +fy*(u(1,j+1,:)+u(1,j-1,:))+(1.0_dp-fc)*u(1,j,:))
    enddo
    respsi(1,nyf,:)=rhs(1,nyf,:)-(fx*(u(2,nyf,:)+Lq*u(1,nyf,:)+Ls) &
      +fy*(rbufN(1,:)+u(1,nyf-1,:))+(1.0_dp-fc)*u(1,nyf,:)) !.global NW corner
  else
    respsi(1,1,:)=rhs(1,1,:)-(fx*(u(2,1,:)+rbufW(1,:)) &
      +fy*(u(1,2,:)+rbufS(1,:))+(1.0_dp-fc)*u(1,1,:)) !.local SW corner
    do j=2,nyf-1
      respsi(1,j,:)=rhs(1,j,:)-(fx*(u(2,j,:)+rbufW(j,:)) &
        +fy*(u(1,j+1,:)+u(1,j-1,:))+(1.0_dp-fc)*u(1,j,:))
    enddo
    respsi(1,nyf,:)=rhs(1,nyf,:)-(fx*(u(2,nyf,:)+rbufW(nyf,:)) &
      +fy*(rbufN(1,:)+u(1,nyf-1,:))+(1.0_dp-fc)*u(1,nyf,:)) !.local NW corner
  endif
!.Right boundary
  call MPI_WAIT(req(8),stat(:,8),ierr)
  if (iIDps(ng) == iprocsm1) then !.Right most process (along x)
    respsi(nxf,1,:)=rhs(nxf,1,:)-(fx*(Rq*u(nxf,1,:)+Rs+u(nxf-1,1,:)) &
      +fy*(u(nxf,2,:)+rbufS(nxf,:))+(1.0_dp-fc)*u(nxf,1,:)) !.global SE corner
    do j=2,nyf-1
      respsi(nxf,j,:)=rhs(nxf,j,:)-(fx*(Rq*u(nxf,j,:)+Rs+u(nxf-1,j,:)) &
        +fy*(u(nxf,j+1,:)+u(nxf,j-1,:))+(1.0_dp-fc)*u(nxf,j,:))
    enddo
    respsi(nxf,nyf,:)=rhs(nxf,nyf,:)-(fx*(Rq*u(nxf,nyf,:)+Rs+u(nxf-1,nyf,:)) &
      +fy*(rbufN(nxf,:)+u(nxf,nyf-1,:))+(1.0_dp-fc)*u(nxf,nyf,:)) !.global NE corner
  else
    respsi(nxf,1,:)=rhs(nxf,1,:)-(fx*(rbufE(1,:)+u(nxf-1,1,:)) &
      +fy*(u(nxf,2,:)+rbufS(nxf,:))+(1.0_dp-fc)*u(nxf,1,:)) !.local SE corner
    do j=2,nyf-1
      respsi(nxf,j,:)=rhs(nxf,j,:)-(fx*(rbufE(j,:)+u(nxf-1,j,:)) &
        +fy*(u(nxf,j+1,:)+u(nxf,j-1,:))+(1.0_dp-fc)*u(nxf,j,:))
    enddo
    respsi(nxf,nyf,:)=rhs(nxf,nyf,:)-(fx*(rbufE(nyf,:)+u(nxf-1,nyf,:)) &
      +fy*(rbufN(nxf,:)+u(nxf,nyf-1,:))+(1.0_dp-fc)*u(nxf,nyf,:)) !.local NE corner
  endif
  call MPI_WAITALL(4,(/req(1),req(2),req(3),req(4)/), &
    (/stat(:,1),stat(:,2),stat(:,3),stat(:,4)/),ierr)

#IF (MGTIMER==1)
    call CPU_TIME(rest2)
    restps=restps+rest2-rest1
#ENDIF
  end function respsi
!*********************************************************
#IF (HDop==4)
  recursive subroutine Gcychd(ng,u1,u2,rhs1,rhs2,xBC,qn)
!.Gamma cycles to solve equation from implicit 4th order
!.Laplacian-based hyperdiffusion equation
!.                   u + Diff*(L^2)u=rho
!.(L=Laplacian) splitting the problem into a system of
!.two 2nd order eqns.
!.INPUTS
!. integer ng: current grid
!. real(dp) dimension(:,:,:) u1, u2: current solution at grid ng
!. real(dp) dimension(:,:,:) rhs1, rhs2: the right-hand side rho vector
!. integer(2) xBC: boundary conditions along x
!. integer qn: which quantity is being hyperdiffused (1 to nqhd)
!.OUTPUTS
!. real(dp) dimension(:,:,:) u1, u2: updated solution at grid ng
  implicit none
  integer(I4B), intent(in) :: ng, xBC(:), qn
  real(DP), intent(inout) :: u1(:,:,:),u2(:,:,:)
  real(DP), intent(in) :: rhs1(:,:,:),rhs2(:,:,:)
  integer(I4B) :: jpost, jpre, i
  real(DP), allocatable :: res1f(:,:,:),res1c(:,:,:),v1(:,:,:),v2(:,:,:), &
                           res2f(:,:,:),res2c(:,:,:)

  iprocsm1=iprocshd(qn)%a(ng)-1 !.ID of last rank in x COMM
  jprocsm1=jprocshd(qn)%a(ng)-1 !.ID of last rank in y COMM

  if (ng==1) then        !.Bottom of gamma cycle: Solve on coarsest grid
#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpost=1,NU2hd(qn)
      call relaxhd(ng,u1,u2,rhs1,rhs2,xBC,qn)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    postXthd=postXthd+count2-count1
#ENDIF
  else                  !.On downward stroke towards coarser grids
#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpre=1,NU1hd(qn)      !.Pre-smoothing
      call relaxhd(ng,u1,u2,rhs1,rhs2,xBC,qn)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    preXthd=preXthd+count2-count1
#ENDIF

!...Number of grid points in coarse grid
    nxc=nxChd(qn)%a(ng)
    nyc=nyChd(qn)%a(ng)
!...Number of grid points in coarse grid after subdomain unification
    nxu=nxFhd(qn)%a(ng-1)
    nyu=nyFhd(qn)%a(ng-1)
    allocate(res1f(nxf,nyf,nzL),res1c(nxu,nyu,nzL),v1(nxu,nyu,nzL))
    allocate(res2f(nxf,nyf,nzL),res2c(nxu,nyu,nzL),v2(nxu,nyu,nzL))

#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    call reshd(ng,res1f,res2f,u1,u2,rhs1,rhs2,xBC,qn)   !.Fine grid residual of each equation

!...Restriction of the residual is the next RHS
    res1c=rstrct(res1f,wdchd(qn)%a(ng),xBC,acuhd(qn)%a(ng-1,:), &
                 urnkhd(qn)%a(ng-1)%a,nurnkhd(qn)%a(ng-1,:))
    res2c=rstrct(res2f,wdchd(qn)%a(ng),xBC,acuhd(qn)%a(ng-1,:), &
                 urnkhd(qn)%a(ng-1)%a,nurnkhd(qn)%a(ng-1,:))
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    rNrthd=rNrthd+count2-count1
#ENDIF

!...Next fine grid is the coarse, or unified subdomain, grid
    nxf=nxu
    nyf=nyu

!...Zero for initial guess in next relaxation
    v1=0.0_dp
    v2=0.0_dp
    if (acuhd(qn)%a(ng-1,1)) then
!.....Only active processes proceed to next coarser grid
      do i=1,GAMhd(qn)
        call Gcychd(ng-1,v1,v2,res1c,res2c,xBC,qn)  !.Recursive call for the coarse grid correction
      enddo
    endif

!...Redefine quantities that were modified in coarser grids
    iprocsm1=iprocshd(qn)%a(ng)-1
    jprocsm1=jprocshd(qn)%a(ng)-1
    nxc=nxChd(qn)%a(ng)
    nyc=nyChd(qn)%a(ng)
    nxf=nxFhd(qn)%a(ng)
    nyf=nyFhd(qn)%a(ng)

#IF (MGTIMER==1)
    call CPU_TIME(prolt1)
#ENDIF
!...Upward stroke -> finer grid: prolong coarse grid error and
!...correct each solution
    u1=u1+prolng(v1,wdchd(qn)%a(ng),xBC,COMMhd(qn)%a(ng,:), &
       iIDhd(qn)%a(ng),jprocshd(qn)%a(ng),nborhd(qn)%a(ng,:), &
       acuhd(qn)%a(ng-1,:),urnkhd(qn)%a(ng-1)%a,nurnkhd(qn)%a(ng-1,:))
    u2=u2+prolng(v2,wdchd(qn)%a(ng),xBC,COMMhd(qn)%a(ng,:), &
       iIDhd(qn)%a(ng),jprocshd(qn)%a(ng),nborhd(qn)%a(ng,:), &
       acuhd(qn)%a(ng-1,:),urnkhd(qn)%a(ng-1)%a,nurnkhd(qn)%a(ng-1,:))
#IF (MGTIMER==1)
    call CPU_TIME(prolt2)
    prolthd=prolthd+prolt2-prolt1
#ENDIF

#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpost=1,NU2hd(qn)      !.Post-smoothing
      call relaxhd(ng,u1,u2,rhs1,rhs2,xBC,qn)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    postXthd=postXthd+count2-count1
#ENDIF
  endif

  end subroutine Gcychd
!*********************************************************
  subroutine relaxhd(ng,u1,u2,rhs1,rhs2,xBC,qn)
!.Red-black Gauss-Seidel relaxation system of 2 equations
!.                  u2 + Diff*L u1 = rhs1 (= rho)
!.                       L u2 - u1 = rhs2 (= 0)
!.INPUTS
!. integer ng: current grid
!. real(dp) dimension(:,:,:) u1, u2: current solution at grid ng
!. real(dp) dimension(:,:,:) rhs1, rhs2: the right-hand side of mod Helmholtz equation
!. integer(2) xBC: boundary conditions along x
!. integer qn: quantity being hyperdiffused (1 to nqhd)
!.OUTPUTS
!. real(dp) dimension(:,:,:) u1, u2: updated solution at grid ng
  implicit none
  include 'mpif.h'
  real(DP), intent(inout) :: u1(:,:,:),u2(:,:,:)
  real(DP), intent(in) :: rhs1(:,:,:),rhs2(:,:,:)
  integer(I4B), intent(in) :: ng,xBC(:),qn
  integer(I4B) :: i,j,k,Nrank,Srank,Erank,Wrank,nxfd2,nyfd2
  integer(I4B) :: ierr,req(8),stat(MPI_STATUS_SIZE,8)
  real(DP), allocatable :: Ls1(:),Rs1(:),Ls2(:),Rs2(:),loc1D(:),b1(:),b2(:)
  real(DP), allocatable, dimension(:,:,:) :: sbufN,sbufE,sbufS,sbufW, &
                                             rbufN,rbufE,rbufS,rbufW
  real(DP) :: fx,fy,fc,fxD,fyD,fcD,delta,Lq,Rq,deltaL,deltaR

  nxfd2=nxf/2
  nyfd2=nyf/2

!.Neigboring ranks to communicate with
  Nrank = nborhd(qn)%a(ng,1)
  Erank = nborhd(qn)%a(ng,3)
  Srank = nborhd(qn)%a(ng,5)
  Wrank = nborhd(qn)%a(ng,7)

!.Boundary conditions are implemented with some amount of flexibility.
!.The right boundary for example, if u(nxf+1) is needed, it will be
!.written as
!.  u(nxf+1) = Rq*u(nxf) + Rr*u(nxf-1) + Rs
!.so that symmetry, extrapolation or Dirichlet conditions can be
!.easily imposed.
!......... BCs of u ................................!
  allocate(Ls1(nzL),Rs1(nzL),loc1D(nzL))
  allocate(Ls2(nzL),Rs2(nzL))
!.LHS BC of u
  if (xBC(1)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Lq  = -1.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u1(1,:,k))
    call MPI_ALLreduce(loc1D,Ls1,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Ls1 =  2.0_dp*Ls1/dble(nyf*jprocshd(qn)%a(ng))
    forall(k=1:nzL) loc1D(k)=sum(u2(1,:,k))
    call MPI_ALLreduce(loc1D,Ls2,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Ls2 =  2.0_dp*Ls2/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC==1 .OR. xBC(1)==-1) then !.symmetry
    Lq  = dble(xBC(1))
    Ls1 = 0.0_dp
    Ls2 = 0.0_dp
  endif
!.RHS BC of u
  if (xBC(2)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Rq  = -1.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u1(nxf,:,k))
    call MPI_ALLreduce(loc1D,Rs1,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Rs1 =  2.0_dp*Rs1/dble(nyf*jprocshd(qn)%a(ng))
    forall(k=1:nzL) loc1D(k)=sum(u2(nxf,:,k))
    call MPI_ALLreduce(loc1D,Rs2,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Rs2 =  2.0_dp*Rs2/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC(2)==1 .OR. xBC(2)==-1) then !.symmetry
    Rq  = dble(xBC(2))
    Rs1 = 0.0_dp
    Rs2 = 0.0_dp
  endif
!........................................................!

  allocate(b1(nzL),b2(nzL))
  allocate(sbufW(nyfd2,nzL,HDopd2),sbufS(nxfd2,nzL,HDopd2))
  allocate(rbufN(nxfd2,nzL,HDopd2),rbufE(nyfd2,nzL,HDopd2))
!.Send bottom row to South rank, receive row from North rank
  forall(i=2:nxf:2,k=1:nzL) sbufS(i/2,k,:)=(/u1(i,1,k),u2(i,1,k)/)
  call MPI_ISEND(sbufS,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(1),ierr)
!.Receive row from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(2),ierr)
!.Send left column to West rank, receive column from East rank
  forall(j=2:nyf:2,k=1:nzL) sbufW(j/2,k,:)=(/u1(1,j,k),u2(1,j,k)/)
  call MPI_ISEND(sbufW,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(3),ierr)
!.Receive column from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(4),ierr)
  
  fx    = (dble(nxf*iprocshd(qn)%a(ng))/bLx)**2
  fy    = (dble(nyf*jprocshd(qn)%a(ng))/bLy)**2
  fc    = 2.0_dp*(fx+fy)
  fxD   = fx*hDiff(1)
  fyD   = fy*hDiff(2)
!  fcD   = fc*Diff
  fcD   = 2.0_dp*(hDiff(1)*fx+hDiff(2)*fy)
  delta = omehd(qn)/(fcD*fc+1.0_dp)
  deltaL= omehd(qn)/((fx*Lq-fc)*(fxD*Lq-fcD)+1.0_dp)
  deltaR= omehd(qn)/((fx*Rq-fc)*(fxD*Rq-fcD)+1.0_dp)

!.Relax inner points with even i and even j
  do j=2,nyf-2,2
    do i=2,nxf-2,2
      b1=rhs1(i,j,:)-fxD*(u1(i+1,j,:)+u1(i-1,j,:))-fyD*(u1(i,j+1,:)+u1(i,j-1,:))
      b2=rhs2(i,j,:)-fx*(u2(i+1,j,:)+u2(i-1,j,:))-fy*(u2(i,j+1,:)+u2(i,j-1,:))
      u1(i,j,:)=Oomehd(qn)*u1(i,j,:)+(-fc*b1-b2)*delta
      u2(i,j,:)=Oomehd(qn)*u2(i,j,:)+(-fcD*b2+b1)*delta
    enddo
  enddo
!.Top boundary, even i
  call MPI_WAIT(req(2),stat(:,2),ierr)
  do i=2,nxf-2,2
    b1=rhs1(i,nyf,:)-fxD*(u1(i+1,nyf,:)+u1(i-1,nyf,:))-fyD*(rbufN(i/2,:,1)+u1(i,nyf-1,:))
    b2=rhs2(i,nyf,:)-fx*(u2(i+1,nyf,:)+u2(i-1,nyf,:))-fy*(rbufN(i/2,:,2)+u2(i,nyf-1,:))
    u1(i,nyf,:)=Oomehd(qn)*u1(i,nyf,:)+(-fc*b1-b2)*delta
    u2(i,nyf,:)=Oomehd(qn)*u2(i,nyf,:)+(-fcD*b2+b1)*delta
  enddo
!.Right boundary, even j
  call MPI_WAIT(req(4),stat(:,4),ierr)
  if (iIDhd(qn)%a(ng) .eq. iprocsm1) then !.Right most process (along x)
    do j=2,nyf-2,2
      b1=rhs1(nxf,j,:)-fxD*(Rs1+u1(nxf-1,j,:))-fyD*(u1(nxf,j+1,:)+u1(nxf,j-1,:))
      b2=rhs2(nxf,j,:)-fx*(Rs2+u2(nxf-1,j,:))-fy*(u2(nxf,j+1,:)+u2(nxf,j-1,:))
      u1(nxf,j,:)=Oomehd(qn)*u1(nxf,j,:)+((fx*Rq-fc)*b1-b2)*deltaR
      u2(nxf,j,:)=Oomehd(qn)*u2(nxf,j,:)+((fxD*Rq-fcD)*b2+b1)*deltaR
    enddo
!...global NE corner
    b1=rhs1(nxf,nyf,:)-fxD*(Rs1+u1(nxf-1,nyf,:))-fyD*(rbufN(nxfd2,:,1)+u1(nxf,nyf-1,:))
    b2=rhs2(nxf,nyf,:)-fx*(Rs2+u2(nxf-1,nyf,:))-fy*(rbufN(nxfd2,:,2)+u2(nxf,nyf-1,:))
    u1(nxf,nyf,:)=Oomehd(qn)*u1(nxf,nyf,:)+((fx*Rq-fc)*b1-b2)*deltaR
    u2(nxf,nyf,:)=Oomehd(qn)*u2(nxf,nyf,:)+((fxD*Rq-fcD)*b2+b1)*deltaR
  else
    do j=2,nyf-2,2
      b1=rhs1(nxf,j,:)-fxD*(rbufE(j/2,:,1)+u1(nxf-1,j,:))-fyD*(u1(nxf,j+1,:)+u1(nxf,j-1,:))
      b2=rhs2(nxf,j,:)-fx*(rbufE(j/2,:,2)+u2(nxf-1,j,:))-fy*(u2(nxf,j+1,:)+u2(nxf,j-1,:))
      u1(nxf,j,:)=Oomehd(qn)*u1(nxf,j,:)+(-fc*b1-b2)*delta
      u2(nxf,j,:)=Oomehd(qn)*u2(nxf,j,:)+(-fcD*b2+b1)*delta
    enddo
!...local NE corner
    b1=rhs1(nxf,nyf,:)-fxD*(rbufE(nyfd2,:,1)+u1(nxf-1,nyf,:))-fyD*(rbufN(nxfd2,:,1)+u1(nxf,nyf-1,:))
    b2=rhs2(nxf,nyf,:)-fx*(rbufE(nyfd2,:,2)+u2(nxf-1,nyf,:))-fy*(rbufN(nxfd2,:,2)+u2(nxf,nyf-1,:))
    u1(nxf,nyf,:)=Oomehd(qn)*u1(nxf,nyf,:)+(-fc*b1-b2)*delta
    u2(nxf,nyf,:)=Oomehd(qn)*u2(nxf,nyf,:)+(-fcD*b2+b1)*delta
  endif

  allocate(sbufN(nxfd2,nzL,HDopd2),sbufE(nyfd2,nzL,HDopd2))
  allocate(rbufW(nyfd2,nzL,HDopd2),rbufS(nxfd2,nzL,HDopd2))
!.Send top row to North rank, receive row from South rank
  forall(i=1:nxf-1:2,k=1:nzL) sbufN((i+1)/2,k,:)=(/u1(i,nyf,k),u2(i,nyf,k)/)
  call MPI_ISEND(sbufN,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(5),ierr)
!.Receive rows from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(6),ierr)
!.Send right column to East rank, receive column from West rank
  forall(j=1:nyf-1:2,k=1:nzL) sbufE((j+1)/2,k,:)=(/u1(nxf,j,k),u2(nxf,j,k)/)
  call MPI_ISEND(sbufE,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(7),ierr)
!.Receive columns from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(8),ierr)

!.Relax inner points with odd i and odd j
  do j=3,nyf-1,2
    do i=3,nxf-1,2
      b1=rhs1(i,j,:)-fxD*(u1(i+1,j,:)+u1(i-1,j,:))-fyD*(u1(i,j+1,:)+u1(i,j-1,:))
      b2=rhs2(i,j,:)-fx*(u2(i+1,j,:)+u2(i-1,j,:))-fy*(u2(i,j+1,:)+u2(i,j-1,:))
      u1(i,j,:)=Oomehd(qn)*u1(i,j,:)+(-fc*b1-b2)*delta
      u2(i,j,:)=Oomehd(qn)*u2(i,j,:)+(-fcD*b2+b1)*delta
    enddo
  enddo
!.Bottom boundary, odd i
  call MPI_WAIT(req(6),stat(:,6),ierr)
  do i=3,nxf-1,2
    b1=rhs1(i,1,:)-fxD*(u1(i+1,1,:)+u1(i-1,1,:))-fyD*(u1(i,2,:)+rbufS((i+1)/2,:,1))
    b2=rhs2(i,1,:)-fx*(u2(i+1,1,:)+u2(i-1,1,:))-fy*(u2(i,2,:)+rbufS((i+1)/2,:,2))
    u1(i,1,:)=Oomehd(qn)*u1(i,1,:)+(-fc*b1-b2)*delta
    u2(i,1,:)=Oomehd(qn)*u2(i,1,:)+(-fcD*b2+b1)*delta
  enddo
!.Left boundary, odd j
  call MPI_WAIT(req(8),stat(:,8),ierr)
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
!...global SW corner
    b1=rhs1(1,1,:)-fxD*(u1(2,1,:)+Ls1)-fyD*(u1(1,2,:)+rbufS(1,:,1))
    b2=rhs2(1,1,:)-fx*(u2(2,1,:)+Ls2)-fy*(u2(1,2,:)+rbufS(1,:,2))
    u1(1,1,:)=Oomehd(qn)*u1(1,1,:)+((fx*Lq-fc)*b1-b2)*deltaL
    u2(1,1,:)=Oomehd(qn)*u2(1,1,:)+((fxD*Lq-fcD)*b2+b1)*deltaL
    do j=3,nyf-1,2
      b1=rhs1(1,j,:)-fxD*(u1(2,j,:)+Ls1)-fyD*(u1(1,j+1,:)+u1(1,j-1,:))
      b2=rhs2(1,j,:)-fx*(u2(2,j,:)+Ls2)-fy*(u2(1,j+1,:)+u2(1,j-1,:))
      u1(1,j,:)=Oomehd(qn)*u1(1,j,:)+((fx*Lq-fc)*b1-b2)*deltaL
      u2(1,j,:)=Oomehd(qn)*u2(1,j,:)+((fxD*Lq-fcD)*b2+b1)*deltaL
    enddo
  else
!...local SW corner
    b1=rhs1(1,1,:)-fxD*(u1(2,1,:)+rbufW(1,:,1))-fyD*(u1(1,2,:)+rbufS(1,:,1))
    b2=rhs2(1,1,:)-fx*(u2(2,1,:)+rbufW(1,:,2))-fy*(u2(1,2,:)+rbufS(1,:,2))
    u1(1,1,:)=Oomehd(qn)*u1(1,1,:)+(-fc*b1-b2)*delta
    u2(1,1,:)=Oomehd(qn)*u2(1,1,:)+(-fcD*b2+b1)*delta
    do j=3,nyf-1,2
      b1=rhs1(1,j,:)-fxD*(u1(2,j,:)+rbufW((j+1)/2,:,1))-fyD*(u1(1,j+1,:)+u1(1,j-1,:))
      b2=rhs2(1,j,:)-fx*(u2(2,j,:)+rbufW((j+1)/2,:,2))-fy*(u2(1,j+1,:)+u2(1,j-1,:))
      u1(1,j,:)=Oomehd(qn)*u1(1,j,:)+(-fc*b1-b2)*delta
      u2(1,j,:)=Oomehd(qn)*u2(1,j,:)+(-fcD*b2+b1)*delta
    enddo
  endif

  call MPI_WAITALL(4,(/req(1),req(3),req(5),req(7)/), &
                     (/stat(:,1),stat(:,3),stat(:,5),stat(:,7)/),ierr)

!.Now do even-odd and odd-even squares of the grid, i.e, the black
!.squares of the checker-board
!.Send bottom row to South rank, receive row from North rank
  forall(i=1:nxf-1:2,k=1:nzL) sbufS((i+1)/2,k,:)=(/u1(i,1,k),u2(i,1,k)/)
  call MPI_ISEND(sbufS,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(1),ierr)
!.Receive row from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(2),ierr)
!.Send right column to East rank, receive column from West rank
  forall(j=2:nyf:2,k=1:nzL) sbufE(j/2,k,:)=(/u1(nxf,j,k),u2(nxf,j,k)/)
  call MPI_ISEND(sbufE,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(3),ierr)
!.Receive columns from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(4),ierr)
!.Relax inner points with odd i and even j
  do j=2,nyf-2,2
    do i=3,nxf-1,2
      b1=rhs1(i,j,:)-fxD*(u1(i+1,j,:)+u1(i-1,j,:))-fyD*(u1(i,j+1,:)+u1(i,j-1,:))
      b2=rhs2(i,j,:)-fx*(u2(i+1,j,:)+u2(i-1,j,:))-fy*(u2(i,j+1,:)+u2(i,j-1,:))
      u1(i,j,:)=Oomehd(qn)*u1(i,j,:)+(-fc*b1-b2)*delta
      u2(i,j,:)=Oomehd(qn)*u2(i,j,:)+(-fcD*b2+b1)*delta
    enddo
  enddo
!.Top boundary, odd i
  call MPI_WAIT(req(2),stat(:,2),ierr)
  do i=3,nxf-1,2
    b1=rhs1(i,nyf,:)-fxD*(u1(i+1,nyf,:)+u1(i-1,nyf,:))-fyD*(rbufN((i+1)/2,:,1)+u1(i,nyf-1,:))
    b2=rhs2(i,nyf,:)-fx*(u2(i+1,nyf,:)+u2(i-1,nyf,:))-fy*(rbufN((i+1)/2,:,2)+u2(i,nyf-1,:))
    u1(i,nyf,:)=Oomehd(qn)*u1(i,nyf,:)+(-fc*b1-b2)*delta
    u2(i,nyf,:)=Oomehd(qn)*u2(i,nyf,:)+(-fcD*b2+b1)*delta
  enddo
!.Left boundary, even j
  call MPI_WAIT(req(4),stat(:,4),ierr)
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    do j=2,nyf-2,2
      b1=rhs1(1,j,:)-fxD*(u1(2,j,:)+Ls1)-fyD*(u1(1,j+1,:)+u1(1,j-1,:))
      b2=rhs2(1,j,:)-fx*(u2(2,j,:)+Ls2)-fy*(u2(1,j+1,:)+u2(1,j-1,:))
      u1(1,j,:)=Oomehd(qn)*u1(1,j,:)+((fx*Lq-fc)*b1-b2)*deltaL
      u2(1,j,:)=Oomehd(qn)*u2(1,j,:)+((fxD*Lq-fcD)*b2+b1)*deltaL
    enddo
!...global NW corner
    b1=rhs1(1,nyf,:)-fxD*(u1(2,nyf,:)+Ls1)-fyD*(rbufN(1,:,1)+u1(1,nyf-1,:))
    b2=rhs2(1,nyf,:)-fx*(u2(2,nyf,:)+Ls2)-fy*(rbufN(1,:,2)+u2(1,nyf-1,:))
    u1(1,nyf,:)=Oomehd(qn)*u1(1,nyf,:)+((fx*Lq-fc)*b1-b2)*deltaL
    u2(1,nyf,:)=Oomehd(qn)*u2(1,nyf,:)+((fxD*Lq-fcD)*b2+b1)*deltaL
  else
    do j=2,nyf-2,2
      b1=rhs1(1,j,:)-fxD*(u1(2,j,:)+rbufW(j/2,:,1))-fyD*(u1(1,j+1,:)+u1(1,j-1,:))
      b2=rhs2(1,j,:)-fx*(u2(2,j,:)+rbufW(j/2,:,2))-fy*(u2(1,j+1,:)+u2(1,j-1,:))
      u1(1,j,:)=Oomehd(qn)*u1(1,j,:)+(-fc*b1-b2)*delta
      u2(1,j,:)=Oomehd(qn)*u2(1,j,:)+(-fcD*b2+b1)*delta
    enddo
!...local NW corner
    b1=rhs1(1,nyf,:)-fxD*(u1(2,nyf,:)+rbufW(nyfd2,:,1))-fyD*(rbufN(1,:,1)+u1(1,nyf-1,:))
    b2=rhs2(1,nyf,:)-fx*(u2(2,nyf,:)+rbufW(nyfd2,:,2))-fy*(rbufN(1,:,2)+u2(1,nyf-1,:))
    u1(1,nyf,:)=Oomehd(qn)*u1(1,nyf,:)+(-fc*b1-b2)*delta
    u2(1,nyf,:)=Oomehd(qn)*u2(1,nyf,:)+(-fcD*b2+b1)*delta
  endif

!.Send top row to North rank, receive row from South rank
  forall(i=2:nxf:2,k=1:nzL) sbufN(i/2,k,:)=(/u1(i,nyf,k),u2(i,nyf,k)/)
  call MPI_ISEND(sbufN,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(5),ierr)
!.Receive rows from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(6),ierr)
!.Send left column to West rank, receive column from East rank
  forall(j=1:nyf-1:2,k=1:nzL) sbufW((j+1)/2,k,:)=(/u1(1,j,k),u2(1,j,k)/)
  call MPI_ISEND(sbufW,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(7),ierr)
!.Receive column from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(8),ierr)
!.Relax inner points with even i and odd j
  do j=3,nyf-1,2
    do i=2,nxf-2,2
      b1=rhs1(i,j,:)-fxD*(u1(i+1,j,:)+u1(i-1,j,:))-fyD*(u1(i,j+1,:)+u1(i,j-1,:))
      b2=rhs2(i,j,:)-fx*(u2(i+1,j,:)+u2(i-1,j,:))-fy*(u2(i,j+1,:)+u2(i,j-1,:))
      u1(i,j,:)=Oomehd(qn)*u1(i,j,:)+(-fc*b1-b2)*delta
      u2(i,j,:)=Oomehd(qn)*u2(i,j,:)+(-fcD*b2+b1)*delta
    enddo
  enddo
!.Bottom boundary, even i
  call MPI_WAIT(req(6),stat(:,6),ierr)
  do i=2,nxf-2,2
    b1=rhs1(i,1,:)-fxD*(u1(i+1,1,:)+u1(i-1,1,:))-fyD*(u1(i,2,:)+rbufS(i/2,:,1))
    b2=rhs2(i,1,:)-fx*(u2(i+1,1,:)+u2(i-1,1,:))-fy*(u2(i,2,:)+rbufS(i/2,:,2))
    u1(i,1,:)=Oomehd(qn)*u1(i,1,:)+(-fc*b1-b2)*delta
    u2(i,1,:)=Oomehd(qn)*u2(i,1,:)+(-fcD*b2+b1)*delta
  enddo
!.Right boundary, odd j
  call MPI_WAIT(req(8),stat(:,8),ierr)
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
!...global SE corner
    b1=rhs1(nxf,1,:)-fxD*(Rs1+u1(nxf-1,1,:))-fyD*(u1(nxf,2,:)+rbufS(nxfd2,:,1))
    b2=rhs2(nxf,1,:)-fx*(Rs2+u2(nxf-1,1,:))-fy*(u2(nxf,2,:)+rbufS(nxfd2,:,2))
    u1(nxf,1,:)=Oomehd(qn)*u1(nxf,1,:)+((fx*Rq-fc)*b1-b2)*deltaR
    u2(nxf,1,:)=Oomehd(qn)*u2(nxf,1,:)+((fxD*Rq-fcD)*b2+b1)*deltaR
    do j=3,nyf-1,2
      b1=rhs1(nxf,j,:)-fxD*(Rs1+u1(nxf-1,j,:))-fyD*(u1(nxf,j+1,:)+u1(nxf,j-1,:))
      b2=rhs2(nxf,j,:)-fx*(Rs2+u2(nxf-1,j,:))-fy*(u2(nxf,j+1,:)+u2(nxf,j-1,:))
      u1(nxf,j,:)=Oomehd(qn)*u1(nxf,j,:)+((fx*Rq-fc)*b1-b2)*deltaR
      u2(nxf,j,:)=Oomehd(qn)*u2(nxf,j,:)+((fxD*Rq-fcD)*b2+b1)*deltaR
    enddo
  else
!...local SE corner
    b1=rhs1(nxf,1,:)-fxD*(rbufE(1,:,1)+u1(nxf-1,1,:))-fyD*(u1(nxf,2,:)+rbufS(nxfd2,:,1))
    b2=rhs2(nxf,1,:)-fx*(rbufE(1,:,2)+u2(nxf-1,1,:))-fy*(u2(nxf,2,:)+rbufS(nxfd2,:,2))
    u1(nxf,1,:)=Oomehd(qn)*u1(nxf,1,:)+(-fc*b1-b2)*delta
    u2(nxf,1,:)=Oomehd(qn)*u2(nxf,1,:)+(-fcD*b2+b1)*delta
    do j=3,nyf-1,2
      b1=rhs1(nxf,j,:)-fxD*(rbufE((j+1)/2,:,1)+u1(nxf-1,j,:))-fyD*(u1(nxf,j+1,:)+u1(nxf,j-1,:))
      b2=rhs2(nxf,j,:)-fx*(rbufE((j+1)/2,:,2)+u2(nxf-1,j,:))-fy*(u2(nxf,j+1,:)+u2(nxf,j-1,:))
      u1(nxf,j,:)=Oomehd(qn)*u1(nxf,j,:)+(-fc*b1-b2)*delta
      u2(nxf,j,:)=Oomehd(qn)*u2(nxf,j,:)+(-fcD*b2+b1)*delta
    enddo
  endif

  call MPI_WAITALL(4,(/req(1),req(3),req(5),req(7)/), &
                     (/stat(:,1),stat(:,3),stat(:,5),stat(:,7)/),ierr)

  end subroutine relaxhd
!*********************************************************
  subroutine reshd(ng,res1,res2,u1,u2,rhs1,rhs2,xBC,qn)
!.Returns minus the residual of each equation in the split approach
!.to solving the 4th order hyperdiffusion problem
!.      res1 = rhs1 - (u2 + Diff*L u1)
!.      res2 = rhs2 - (L u2 - u1) 
!.(L=Laplacian) using a 2nd order accurate finite difference stencil.
!.INPUTS
!. integer ng: current grid
!. real(dp) dimension(:,:,:) u1, u2: current solution at grid ng
!. real(dp) dimension(:,:,:) rhs1, rhs2: the right-hand side of Poisson equation
!. integer(2) xBC: boundary conditions along x
!. integer qn: quantity being hyperdiffused (1 to nqhd)
!.OUTPUTS
!. real(dp) dimension(:,:,:) res1, res2: negated residual
  implicit none
  include 'mpif.h'
  real(DP), intent(out) :: res1(:,:,:),res2(:,:,:)
  real(DP), intent(in) :: u1(:,:,:),u2(:,:,:),rhs1(:,:,:),rhs2(:,:,:)
  integer(I4B), intent(in) :: ng,xBC(:),qn
  integer(I4B) :: i,j,k,Nrank,Srank,Erank,Wrank
  real(DP) :: fx,fy,fc,fxD,fyD,fcD,Lq,Rq
  real(DP), allocatable :: Ls1(:),Ls2(:),Rs1(:),Rs2(:),loc1D(:)
  real(DP), allocatable, dimension(:,:,:) :: sbufW,rbufE,sbufE,rbufW, &
                                             sbufN,rbufS,sbufS,rbufN
  integer(I4B) :: req(8),stat(MPI_STATUS_SIZE,8),ierr

#IF (MGTIMER==1)
    call CPU_TIME(rest1)
#ENDIF

!.Neigboring ranks to communicate with
  Nrank =nborhd(qn)%a(ng,1)
  Erank =nborhd(qn)%a(ng,3)
  Srank =nborhd(qn)%a(ng,5)
  Wrank =nborhd(qn)%a(ng,7)

!.Boundary conditions are implemented with some amount of flexibility.
!.The right boundary for example, if u(nxf+1) is needed, it will be
!.written as
!.  u(nxf+1) = Rq*u(nxf) + Rr*u(nxf-1) + Rs
!.so that symmetry, extrapolation or Dirichlet conditions can be
!.easily imposed.
!......... BCs of u ................................!
  allocate(Ls1(nzL),Rs1(nzL),loc1D(nzL))
  allocate(Ls2(nzL),Rs2(nzL))
!.LHS BC of u
  if (xBC(1)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Lq  = -1.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u1(1,:,k))
    call MPI_ALLreduce(loc1D,Ls1,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Ls1 =  2.0_dp*Ls1/dble(nyf*jprocshd(qn)%a(ng))
    forall(k=1:nzL) loc1D(k)=sum(u2(1,:,k))
    call MPI_ALLreduce(loc1D,Ls2,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Ls2 = 2.0_dp*Ls2/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC==1 .OR. xBC(1)==-1) then !.symmetry
    Lq  = dble(xBC(1))
    Ls1 = 0.0_dp
    Ls2 = 0.0_dp
  endif
!.RHS BC of u
  if (xBC(2)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Rq  = -1.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u1(nxf,:,k))
    call MPI_ALLreduce(loc1D,Rs1,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Rs1 =  2.0_dp*Rs1/dble(nyf*jprocshd(qn)%a(ng))
    forall(k=1:nzL) loc1D(k)=sum(u2(nxf,:,k))
    call MPI_ALLreduce(loc1D,Rs2,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Rs2 =  2.0_dp*Rs2/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC(2)==1 .OR. xBC(2)==-1) then !.symmetry
    Rq  = dble(xBC(2))
    Rs1 = 0.0_dp
    Rs2 = 0.0_dp
  endif
!........................................................!

  allocate(sbufW(nyf,nzL,HDopd2),sbufS(nxf,nzL,HDopd2), &
           sbufN(nxf,nzL,HDopd2),sbufE(nyf,nzL,HDopd2))
  allocate(rbufW(nyf,nzL,HDopd2),rbufS(nxf,nzL,HDopd2), &
           rbufN(nxf,nzL,HDopd2),rbufE(nyf,nzL,HDopd2))
!.Send top row to North rank, receive row from South rank
  forall(i=1:nxf,k=1:nzL) sbufN(i,k,:)=(/u1(i,nyf,k),u2(i,nyf,k)/)
  call MPI_ISEND(sbufN,nxf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(1),ierr)
!.Send right column to East rank, receive column from West rank
  forall(j=1:nyf,k=1:nzL) sbufE(j,k,:)=(/u1(nxf,j,k),u2(nxf,j,k)/)
  call MPI_ISEND(sbufE,nyf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(2),ierr)
!.Send bottom row to South rank, receive row from North rank
  forall(i=1:nxf,k=1:nzL) sbufS(i,k,:)=(/u1(i,1,k),u2(i,1,k)/)
  call MPI_ISEND(sbufS,nxf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(3),ierr)
!.Send left column to West rank, receive column from East rank
  forall(j=1:nyf,k=1:nzL) sbufW(j,k,:)=(/u1(1,j,k),u2(1,j,k)/)
  call MPI_ISEND(sbufW,nyf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(4),ierr)
!.Receive rows from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(5),ierr)
!.Receive columns from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,nyf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(6),ierr)
!.Receive row from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(7),ierr)
!.Receive column from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,nyf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(8),ierr)

  fx  = (dble(nxf*iprocshd(qn)%a(ng))/bLx)**2         !.1/(hx^2)
  fy  = (dble(nyf*jprocshd(qn)%a(ng))/bLy)**2         !.1/(hy^2)
  fc  = 2.0_dp*(fx+fy)
  fxD = hDiff(1)*fx
  fyD = hDiff(2)*fy
!  fcD = Diff*fc
  fcD = 2.0_dp*(hDiff(1)*fx+hDiff(2)*fy)

!.Interior points
  do j=2,nyf-1
    do i=2,nxf-1
      res1(i,j,:)=rhs1(i,j,:)-(u2(i,j,:)+fxD*(u1(i+1,j,:)+u1(i-1,j,:)) &
        +fyD*(u1(i,j+1,:)+u1(i,j-1,:))-fcD*u1(i,j,:))
      res2(i,j,:)=rhs2(i,j,:)-(fx*(u2(i+1,j,:)+u2(i-1,j,:)) &
        +fy*(u2(i,j+1,:)+u2(i,j-1,:))-fc*u2(i,j,:)-u1(i,j,:))
    enddo
  enddo

!.Bottom boundary
  call MPI_WAIT(req(5),stat(:,5),ierr)
  do i=2,nxf-1
    res1(i,1,:)=rhs1(i,1,:)-(u2(i,1,:)+fxD*(u1(i+1,1,:)+u1(i-1,1,:)) &
      +fyD*(u1(i,2,:)+rbufS(i,:,1))-fcD*u1(i,1,:))
    res2(i,1,:)=rhs2(i,1,:)-(fx*(u2(i+1,1,:)+u2(i-1,1,:)) &
      +fy*(u2(i,2,:)+rbufS(i,:,2))-fc*u2(i,1,:)-u1(i,1,:))
  enddo
!.Left boundary
  call MPI_WAIT(req(6),stat(:,6),ierr)
  if (iIDhd(qn)%a(ng) == 0) then !.Left most process (along x)
!...global SW corner
    res1(1,1,:)=rhs1(1,1,:)-(u2(1,1,:)+fxD*(u1(2,1,:)+Lq*u1(1,1,:)+Ls1) &
      +fyD*(u1(1,2,:)+rbufS(1,:,1))-fcD*u1(1,1,:))
    res2(1,1,:)=rhs2(1,1,:)-(fx*(u2(2,1,:)+Lq*u2(1,1,:)+Ls2) &
      +fy*(u2(1,2,:)+rbufS(1,:,2))-fc*u2(1,1,:)-u1(1,1,:))
    do j=2,nyf-1
      res1(1,j,:)=rhs1(1,j,:)-(u2(1,j,:)+fxD*(u1(2,j,:)+Lq*u1(1,j,:)+Ls1) &
        +fyD*(u1(1,j+1,:)+u1(1,j-1,:))-fcD*u1(1,j,:))
      res2(1,j,:)=rhs2(1,j,:)-(fx*(u2(2,j,:)+Lq*u2(1,j,:)+Ls2) &
        +fy*(u2(1,j+1,:)+u2(1,j-1,:))-fc*u2(1,j,:)-u1(1,j,:))
    enddo
  else
!...local SW corner
    res1(1,1,:)=rhs1(1,1,:)-(u2(1,1,:)+fxD*(u1(2,1,:)+rbufW(1,:,1)) &
      +fyD*(u1(1,2,:)+rbufS(1,:,1))-fcD*u1(1,1,:))
    res2(1,1,:)=rhs2(1,1,:)-(fx*(u2(2,1,:)+rbufW(1,:,2)) &
      +fy*(u2(1,2,:)+rbufS(1,:,2))-fc*u2(1,1,:)-u1(1,1,:))
    do j=2,nyf-1
      res1(1,j,:)=rhs1(1,j,:)-(u2(1,j,:)+fxD*(u1(2,j,:)+rbufW(j,:,1)) &
        +fyD*(u1(1,j+1,:)+u1(1,j-1,:))-fcD*u1(1,j,:))
      res2(1,j,:)=rhs2(1,j,:)-(fx*(u2(2,j,:)+rbufW(j,:,2)) &
        +fy*(u2(1,j+1,:)+u2(1,j-1,:))-fc*u2(1,j,:)-u1(1,j,:))
    enddo
  endif

  call MPI_WAIT(req(7),stat(:,7),ierr)
  if (iIDhd(qn)%a(ng) == 0) then !.Left most process (along x)
!...global NW corner
    res1(1,nyf,:)=rhs1(1,nyf,:)-(u2(1,nyf,:)+fxD*(u1(2,nyf,:)+Lq*u1(1,nyf,:)+Ls1) &
      +fyD*(rbufN(1,:,1)+u1(1,nyf-1,:))-fcD*u1(1,nyf,:))
    res2(1,nyf,:)=rhs2(1,nyf,:)-(fx*(u2(2,nyf,:)+Lq*u2(1,nyf,:)+Ls2) &
      +fy*(rbufN(1,:,2)+u2(1,nyf-1,:))-fc*u2(1,nyf,:)-u1(1,nyf,:))
  else
!...local NW corner
    res1(1,nyf,:)=rhs1(1,nyf,:)-(u2(1,nyf,:)+fxD*(u1(2,nyf,:)+rbufW(nyf,:,1)) &
      +fyD*(rbufN(1,:,1)+u1(1,nyf-1,:))-fcD*u1(1,nyf,:))
    res2(1,nyf,:)=rhs2(1,nyf,:)-(fx*(u2(2,nyf,:)+rbufW(nyf,:,2)) &
      +fy*(rbufN(1,:,2)+u2(1,nyf-1,:))-fc*u2(1,nyf,:)-u1(1,nyf,:))
  endif
!.Top boundary
  do i=2,nxf-1
    res1(i,nyf,:)=rhs1(i,nyf,:)-(u2(i,nyf,:)+fxD*(u1(i+1,nyf,:)+u1(i-1,nyf,:)) &
      +fyD*(rbufN(i,:,1)+u1(i,nyf-1,:))-fcD*u1(i,nyf,:))
    res2(i,nyf,:)=rhs2(i,nyf,:)-(fx*(u2(i+1,nyf,:)+u2(i-1,nyf,:)) &
      +fy*(rbufN(i,:,2)+u2(i,nyf-1,:))-fc*u2(i,nyf,:)-u1(i,nyf,:))
  enddo

!.Right boundary
  call MPI_WAIT(req(8),stat(:,8),ierr)
  if (iIDhd(qn)%a(ng) == iprocsm1) then !.Right most process (along x)
!...global NE corner
    res1(nxf,nyf,:)=rhs1(nxf,nyf,:)-(u2(nxf,nyf,:)+fxD*(Rq*u1(nxf,nyf,:)+Rs1 &
      +u1(nxf-1,nyf,:))+fyD*(rbufN(nxf,:,1)+u1(nxf,nyf-1,:))-fcD*u1(nxf,nyf,:))
    res2(nxf,nyf,:)=rhs2(nxf,nyf,:)-(fx*(Rq*u2(nxf,nyf,:)+Rs2+u2(nxf-1,nyf,:)) &
      +fy*(rbufN(nxf,:,2)+u2(nxf,nyf-1,:))-fc*u2(nxf,nyf,:)-u1(nxf,nyf,:))
!...global SE corner
    res1(nxf,1,:)=rhs1(nxf,1,:)-(u2(nxf,1,:)+fxD*(Rq*u1(nxf,1,:)+Rs1+u1(nxf-1,1,:)) &
      +fyD*(u1(nxf,2,:)+rbufS(nxf,:,1))-fcD*u1(nxf,1,:))
    res2(nxf,1,:)=rhs2(nxf,1,:)-(fx*(Rq*u2(nxf,1,:)+Rs2+u2(nxf-1,1,:)) &
      +fy*(u2(nxf,2,:)+rbufS(nxf,:,2))-fc*u2(nxf,1,:)-u1(nxf,1,:))
    do j=2,nyf-1
      res1(nxf,j,:)=rhs1(nxf,j,:)-(u2(nxf,j,:)+fxD*(Rq*u1(nxf,j,:)+Rs1+u1(nxf-1,j,:)) &
        +fyD*(u1(nxf,j+1,:)+u1(nxf,j-1,:))-fcD*u1(nxf,j,:))
      res2(nxf,j,:)=rhs2(nxf,j,:)-(fx*(Rq*u2(nxf,j,:)+Rs2+u2(nxf-1,j,:)) &
        +fy*(u2(nxf,j+1,:)+u2(nxf,j-1,:))-fc*u2(nxf,j,:)-u1(nxf,j,:))
    enddo
  else
!...local NE corner
    res1(nxf,nyf,:)=rhs1(nxf,nyf,:)-(u2(nxf,nyf,:)+fxD*(rbufE(nyf,:,1)+u1(nxf-1,nyf,:)) &
      +fyD*(rbufN(nxf,:,1)+u1(nxf,nyf-1,:))-fcD*u1(nxf,nyf,:))
    res2(nxf,nyf,:)=rhs2(nxf,nyf,:)-(fx*(rbufE(nyf,:,2)+u2(nxf-1,nyf,:)) &
      +fy*(rbufN(nxf,:,2)+u2(nxf,nyf-1,:))-fc*u2(nxf,nyf,:)-u1(nxf,nyf,:))
!...local SE corner
    res1(nxf,1,:)=rhs1(nxf,1,:)-(u2(nxf,1,:)+fxD*(rbufE(1,:,1)+u1(nxf-1,1,:)) &
      +fyD*(u1(nxf,2,:)+rbufS(nxf,:,1))-fcD*u1(nxf,1,:))
    res2(nxf,1,:)=rhs2(nxf,1,:)-(fx*(rbufE(1,:,2)+u2(nxf-1,1,:)) &
      +fy*(u2(nxf,2,:)+rbufS(nxf,:,2))-fc*u2(nxf,1,:)-u1(nxf,1,:))
    do j=2,nyf-1
      res1(nxf,j,:)=rhs1(nxf,j,:)-(u2(nxf,j,:)+fxD*(rbufE(j,:,1)+u1(nxf-1,j,:)) &
        +fyD*(u1(nxf,j+1,:)+u1(nxf,j-1,:))-fcD*u1(nxf,j,:))
      res2(nxf,j,:)=rhs2(nxf,j,:)-(fx*(rbufE(j,:,2)+u2(nxf-1,j,:)) &
        +fy*(u2(nxf,j+1,:)+u2(nxf,j-1,:))-fc*u2(nxf,j,:)-u1(nxf,j,:))
    enddo
  endif
  call MPI_WAITALL(4,(/req(1),req(2),req(3),req(4)/), &
                     (/stat(:,1),stat(:,2),stat(:,3),stat(:,4)/),ierr)

#IF (MGTIMER==1)
    call CPU_TIME(rest2)
    resthd=resthd+rest2-rest1
#ENDIF
  end subroutine reshd
!*********************************************************
#ELIF (HDop==6)
  recursive subroutine Gcychd(ng,u1,u2,u3,rhs1,rhs2,rhs3,xBC,qn)
!.Gamma cycles to solve equation from implicit 6th order
!.Laplacian-based hyperdiffusion equation
!.                   u + Diff*(L^3)u=rho
!.(L=Laplacian) splitting the problem into a system of
!.three 2nd order eqns.
!.INPUTS
!. integer ng: current grid
!. real(dp) dimension(:,:,:) u1, u2, u3: current solution at grid ng
!. real(dp) dimension(:,:,:) rhs1, rhs2, rhs3: the right-hand side rho vector
!. integer(2) xBC: boundary conditions along x
!. integer qn: which quantity is being hyperdiffused (1 to nqhd)
!.OUTPUTS
!. real(dp) dimension(:,:,:) u1, u2, u3: updated solution at grid ng
  implicit none
  integer(I4B), intent(in) :: ng, xBC(:), qn
  real(DP), intent(inout) :: u1(:,:,:),u2(:,:,:),u3(:,:,:)
  real(DP), intent(in) :: rhs1(:,:,:),rhs2(:,:,:),rhs3(:,:,:)
  integer(I4B) :: jpost, jpre, i
  real(DP), allocatable :: res1f(:,:,:),res1c(:,:,:), &
                           res2f(:,:,:),res2c(:,:,:), &
                           res3f(:,:,:),res3c(:,:,:), &
                           v1(:,:,:),v2(:,:,:),v3(:,:,:)

  iprocsm1=iprocshd(qn)%a(ng)-1 !.ID of last rank in x COMM
  jprocsm1=jprocshd(qn)%a(ng)-1 !.ID of last rank in y COMM

  if (ng==1) then        !.Bottom of gamma cycle: Solve on coarsest grid
#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpost=1,NU2hd(qn)
      call relaxhd(ng,u1,u2,u3,rhs1,rhs2,rhs3,xBC,qn)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    postXthd=postXthd+count2-count1
#ENDIF
  else                  !.On downward stroke towards coarser grids
#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpre=1,NU1hd(qn)      !.Pre-smoothing
      call relaxhd(ng,u1,u2,u3,rhs1,rhs2,rhs3,xBC,qn)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    preXthd=preXthd+count2-count1
#ENDIF

!...Number of grid points in coarse grid
    nxc=nxChd(qn)%a(ng)
    nyc=nyChd(qn)%a(ng)
!...Number of grid points in coarse grid after subdomain unification
    nxu=nxFhd(qn)%a(ng-1)
    nyu=nyFhd(qn)%a(ng-1)
    allocate(res1f(nxf,nyf,nzL),res1c(nxu,nyu,nzL),v1(nxu,nyu,nzL))
    allocate(res2f(nxf,nyf,nzL),res2c(nxu,nyu,nzL),v2(nxu,nyu,nzL))
    allocate(res3f(nxf,nyf,nzL),res3c(nxu,nyu,nzL),v3(nxu,nyu,nzL))

#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    call reshd(ng,res1f,res2f,res3f,u1,u2,u3,rhs1,rhs2,rhs3,xBC,qn)   !.Fine grid residual of each equation

!...Restriction of the residual is the next RHS
    res1c=rstrct(res1f,wdchd(qn)%a(ng),xBC,acuhd(qn)%a(ng-1,:), &
                 urnkhd(qn)%a(ng-1)%a,nurnkhd(qn)%a(ng-1,:))
    res2c=rstrct(res2f,wdchd(qn)%a(ng),xBC,acuhd(qn)%a(ng-1,:), &
                 urnkhd(qn)%a(ng-1)%a,nurnkhd(qn)%a(ng-1,:))
    res3c=rstrct(res3f,wdchd(qn)%a(ng),xBC,acuhd(qn)%a(ng-1,:), &
                 urnkhd(qn)%a(ng-1)%a,nurnkhd(qn)%a(ng-1,:))
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    rNrthd=rNrthd+count2-count1
#ENDIF

!...Next fine grid is the coarse, or unified subdomain, grid
    nxf=nxu
    nyf=nyu

!...Zero for initial guess in next relaxation
    v1=0.0_dp
    v2=0.0_dp
    v3=0.0_dp
    if (acuhd(qn)%a(ng-1,1)) then
!.....Only active processes proceed to next coarser grid
      do i=1,GAMhd(qn)
        call Gcychd(ng-1,v1,v2,v3,res1c,res2c,res3c,xBC,qn)  !.Recursive call for the coarse grid correction
      enddo
    endif

!...Redefine quantities that were modified in coarser grids
    iprocsm1=iprocshd(qn)%a(ng)-1
    jprocsm1=jprocshd(qn)%a(ng)-1
    nxc=nxChd(qn)%a(ng)
    nyc=nyChd(qn)%a(ng)
    nxf=nxFhd(qn)%a(ng)
    nyf=nyFhd(qn)%a(ng)

#IF (MGTIMER==1)
    call CPU_TIME(prolt1)
#ENDIF
!...Upward stroke -> finer grid: prolong coarse grid error
!...and correct each solution
    u1=u1+prolng(v1,wdchd(qn)%a(ng),xBC,COMMhd(qn)%a(ng,:), &
       iIDhd(qn)%a(ng),jprocshd(qn)%a(ng),nborhd(qn)%a(ng,:), &
       acuhd(qn)%a(ng-1,:),urnkhd(qn)%a(ng-1)%a,nurnkhd(qn)%a(ng-1,:))
    u2=u2+prolng(v2,wdchd(qn)%a(ng),xBC,COMMhd(qn)%a(ng,:), &
       iIDhd(qn)%a(ng),jprocshd(qn)%a(ng),nborhd(qn)%a(ng,:), &
       acuhd(qn)%a(ng-1,:),urnkhd(qn)%a(ng-1)%a,nurnkhd(qn)%a(ng-1,:))
    u3=u3+prolng(v3,wdchd(qn)%a(ng),xBC,COMMhd(qn)%a(ng,:), &
       iIDhd(qn)%a(ng),jprocshd(qn)%a(ng),nborhd(qn)%a(ng,:), &
       acuhd(qn)%a(ng-1,:),urnkhd(qn)%a(ng-1)%a,nurnkhd(qn)%a(ng-1,:))
#IF (MGTIMER==1)
    call CPU_TIME(prolt2)
    prolthd=prolthd+prolt2-prolt1
#ENDIF

#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpost=1,NU2hd(qn)      !.Post-smoothing
      call relaxhd(ng,u1,u2,u3,rhs1,rhs2,rhs3,xBC,qn)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    postXthd=postXthd+count2-count1
#ENDIF
  endif

  end subroutine Gcychd
!*********************************************************
  subroutine relaxhd(ng,u1,u2,u3,rhs1,rhs2,rhs3,xBC,qn)
!.Red-black Gauss-Seidel relaxation system of 3 equations
!.                  u3 - Diff*L u1 = rhs1 (= rho)
!.                       L u2 - u1 = rhs2 (= 0)
!.                       L u3 - u2 = rhs3 (= 0)
!.INPUTS
!. integer ng: current grid
!. real(dp) dimension(:,:,:) u1, u2, u3: current solution at grid ng
!. real(dp) dimension(:,:,:) rhs1, rhs2, rhs3: the right-hand side vector
!. integer(2) xBC: boundary conditions along x
!. integer qn: quantity being hyperdiffused (1 to nqhd)
!.OUTPUTS
!. real(dp) dimension(:,:,:) u1, u2: updated solution at grid ng
  implicit none
  include 'mpif.h'
  real(DP), intent(inout) :: u1(:,:,:),u2(:,:,:),u3(:,:,:)
  real(DP), intent(in) :: rhs1(:,:,:),rhs2(:,:,:),rhs3(:,:,:)
  integer(I4B), intent(in) :: ng,xBC(:),qn
  integer(I4B) :: i,j,k,Nrank,Srank,Erank,Wrank,nxfd2,nyfd2
  integer(I4B) :: ierr,req(8),stat(MPI_STATUS_SIZE,8)
  real(DP), allocatable :: Ls1(:),Rs1(:),Ls2(:),Rs2(:), &
                           Ls3(:),Rs3(:),loc1D(:),b1(:),b2(:),b3(:)
  real(DP), allocatable, dimension(:,:,:) :: sbufN,sbufE,sbufS,sbufW, &
                                             rbufN,rbufE,rbufS,rbufW
  real(DP) :: fx,fy,fc,fcsq,fr,frsq,fl,flsq,fxD,fyD,fcD, &
              ffcD,frD,flD,ffrD,fflD,Lq,Rq,delta,deltaL,deltaR

  nxfd2=nxf/2
  nyfd2=nyf/2

!.Neigboring ranks to communicate with
  Nrank = nborhd(qn)%a(ng,1)
  Erank = nborhd(qn)%a(ng,3)
  Srank = nborhd(qn)%a(ng,5)
  Wrank = nborhd(qn)%a(ng,7)

!.Boundary conditions are implemented with some amount of flexibility.
!.The right boundary for example, if u(nxf+1) is needed, it will be
!.written as
!.  u(nxf+1) = Rq*u(nxf) + Rr*u(nxf-1) + Rs
!.so that symmetry, extrapolation or Dirichlet conditions can be
!.easily imposed.
!......... BCs of u ................................!
  allocate(Ls1(nzL),Rs1(nzL),loc1D(nzL))
  allocate(Ls2(nzL),Rs2(nzL),Ls3(nzL),Rs3(nzL))
!.LHS BC of u
  if (xBC(1)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Lq  = -1.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u1(1,:,k))
    call MPI_ALLreduce(loc1D,Ls1,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Ls1 =  2.0_dp*Ls1/dble(nyf*jprocshd(qn)%a(ng))
    forall(k=1:nzL) loc1D(k)=sum(u2(1,:,k))
    call MPI_ALLreduce(loc1D,Ls2,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Ls2 =  2.0_dp*Ls2/dble(nyf*jprocshd(qn)%a(ng))
    forall(k=1:nzL) loc1D(k)=sum(u3(1,:,k))
    call MPI_ALLreduce(loc1D,Ls3,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Ls3 =  2.0_dp*Ls3/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC==1 .OR. xBC(1)==-1) then !.symmetry
    Lq  = dble(xBC(1))
    Ls1 = 0.0_dp
    Ls2 = 0.0_dp
    Ls3 = 0.0_dp
  endif
!.RHS BC of u
  if (xBC(2)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Rq  = -1.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u1(nxf,:,k))
    call MPI_ALLreduce(loc1D,Rs1,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Rs1 =  2.0_dp*Rs1/dble(nyf*jprocshd(qn)%a(ng))
    forall(k=1:nzL) loc1D(k)=sum(u2(nxf,:,k))
    call MPI_ALLreduce(loc1D,Rs2,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Rs2 =  2.0_dp*Rs2/dble(nyf*jprocshd(qn)%a(ng))
    forall(k=1:nzL) loc1D(k)=sum(u3(nxf,:,k))
    call MPI_ALLreduce(loc1D,Rs3,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Rs3 =  2.0_dp*Rs3/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC(2)==1 .OR. xBC(2)==-1) then !.symmetry
    Rq  = dble(xBC(2))
    Rs1 = 0.0_dp
    Rs2 = 0.0_dp
    Rs3 = 0.0_dp
  endif
!........................................................!

  allocate(sbufW(nyfd2,nzL,HDopd2),sbufS(nxfd2,nzL,HDopd2))
  allocate(rbufN(nxfd2,nzL,HDopd2),rbufE(nyfd2,nzL,HDopd2))
!.Send bottom row to South rank, receive row from North rank
  forall(i=2:nxf:2,k=1:nzL) sbufS(i/2,k,:)=(/u1(i,1,k),u2(i,1,k),u3(i,1,k)/)
  call MPI_ISEND(sbufS,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(1),ierr)
!.Receive row from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(2),ierr)
!.Send left column to West rank, receive column from East rank
  forall(j=2:nyf:2,k=1:nzL) sbufW(j/2,k,:)=(/u1(1,j,k),u2(1,j,k),u3(1,j,k)/)
  call MPI_ISEND(sbufW,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(3),ierr)
!.Receive column from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(4),ierr)
  
  fx     = (dble(nxf*iprocshd(qn)%a(ng))/bLx)**2
  fy     = (dble(nyf*jprocshd(qn)%a(ng))/bLy)**2
  fc     = 2.0_dp*(fx+fy)
  fc     = 2.0_dp*(fx+fy)
  fcsq   = fc*fc
  fr     = fx*Rq-fc
  frsq   = fr*fr
  fl     = fx*Lq-fc
  flsq   = fl*fl
  fxD    = fx*hDiff(1)
  fyD    = fy*hDiff(2)
!  fcD    = fc*Diff
  fcD    = 2.0_dp*(hDiff(1)*fx+hDiff(2)*fy)
  ffcD   = fc*fcD
  frD    = fxD*Rq-fcD
  flD    = fxD*Lq-fcD
  ffrD   = fr*frD
  fflD   = fl*flD
  delta  = omehd(qn)/(-fcD*fcsq-1.0_dp)
  deltaL = omehd(qn)/(flsq*flD-1.0_dp)
  deltaR = omehd(qn)/(frsq*frD-1.0_dp)


  allocate(b1(nzL),b2(nzL),b3(nzL))
!.Inner points, even i and even j
  do j=2,nyf-2,2
    do i=2,nxf-2,2
      b1=rhs1(i,j,:)+fxD*(u1(i+1,j,:)+u1(i-1,j,:))+fyD*(u1(i,j+1,:)+u1(i,j-1,:))
      b2=rhs2(i,j,:)-fx*(u2(i+1,j,:)+u2(i-1,j,:))-fy*(u2(i,j+1,:)+u2(i,j-1,:))
      b3=rhs3(i,j,:)-fx*(u3(i+1,j,:)+u3(i-1,j,:))-fy*(u3(i,j+1,:)+u3(i,j-1,:))
      u1(i,j,:)=Oomehd(qn)*u1(i,j,:)+(-fcsq*b1+b2-fc*b3)*delta
      u2(i,j,:)=Oomehd(qn)*u2(i,j,:)+(fc*b1+ffcD*b2+b3)*delta
      u3(i,j,:)=Oomehd(qn)*u3(i,j,:)+(-b1-fcD*b2+ffcD*b3)*delta
    enddo
  enddo
!.Top boundary, even i
  call MPI_WAIT(req(2),stat(:,2),ierr)
  do i=2,nxf-2,2
    b1=rhs1(i,nyf,:)+fxD*(u1(i+1,nyf,:)+u1(i-1,nyf,:))+fyD*(rbufN(i/2,:,1)+u1(i,nyf-1,:))
    b2=rhs2(i,nyf,:)-fx*(u2(i+1,nyf,:)+u2(i-1,nyf,:))-fy*(rbufN(i/2,:,2)+u2(i,nyf-1,:))
    b3=rhs3(i,nyf,:)-fx*(u3(i+1,nyf,:)+u3(i-1,nyf,:))-fy*(rbufN(i/2,:,3)+u3(i,nyf-1,:))
    u1(i,nyf,:)=Oomehd(qn)*u1(i,nyf,:)+(-fcsq*b1+b2-fc*b3)*delta
    u2(i,nyf,:)=Oomehd(qn)*u2(i,nyf,:)+(fc*b1+ffcD*b2+b3)*delta
    u3(i,nyf,:)=Oomehd(qn)*u3(i,nyf,:)+(-b1-fcD*b2+ffcD*b3)*delta
  enddo
!.Right boundary, even j
  call MPI_WAIT(req(4),stat(:,4),ierr)
  if (iIDhd(qn)%a(ng) .eq. iprocsm1) then !.Right most process (along x)
    do j=2,nyf-2,2
      b1=rhs1(nxf,j,:)+fxD*(Rs1+u1(nxf-1,j,:))+fyD*(u1(nxf,j+1,:)+u1(nxf,j-1,:))
      b2=rhs2(nxf,j,:)-fx*(Rs2+u2(nxf-1,j,:))-fy*(u2(nxf,j+1,:)+u2(nxf,j-1,:))
      b3=rhs3(nxf,j,:)-fx*(Rs3+u3(nxf-1,j,:))-fy*(u3(nxf,j+1,:)+u3(nxf,j-1,:))
      u1(nxf,j,:)=Oomehd(qn)*u1(nxf,j,:)+(-frsq*b1+b2+fr*b3)*deltaR
      u2(nxf,j,:)=Oomehd(qn)*u2(nxf,j,:)+(-fr*b1+ffrD*b2+b3)*deltaR
      u3(nxf,j,:)=Oomehd(qn)*u3(nxf,j,:)+(-b1+frD*b2+ffrD*b3)*deltaR
    enddo
!...global NE corner
    b1=rhs1(nxf,nyf,:)+fxD*(Rs1+u1(nxf-1,nyf,:))+fyD*(rbufN(nxfd2,:,1)+u1(nxf,nyf-1,:))
    b2=rhs2(nxf,nyf,:)-fx*(Rs2+u2(nxf-1,nyf,:))-fy*(rbufN(nxfd2,:,2)+u2(nxf,nyf-1,:))
    b3=rhs3(nxf,nyf,:)-fx*(Rs3+u3(nxf-1,nyf,:))-fy*(rbufN(nxfd2,:,3)+u3(nxf,nyf-1,:))
    u1(nxf,nyf,:)=Oomehd(qn)*u1(nxf,nyf,:)+(-frsq*b1+b2+fr*b3)*deltaR
    u2(nxf,nyf,:)=Oomehd(qn)*u2(nxf,nyf,:)+(-fr*b1+ffrD*b2+b3)*deltaR
    u3(nxf,nyf,:)=Oomehd(qn)*u3(nxf,nyf,:)+(-b1+frD*b2+ffrD*b3)*deltaR
  else
    do j=2,nyf-2,2
      b1=rhs1(nxf,j,:)+fxD*(rbufE(j/2,:,1)+u1(nxf-1,j,:))+fyD*(u1(nxf,j+1,:)+u1(nxf,j-1,:))
      b2=rhs2(nxf,j,:)-fx*(rbufE(j/2,:,2)+u2(nxf-1,j,:))-fy*(u2(nxf,j+1,:)+u2(nxf,j-1,:))
      b3=rhs3(nxf,j,:)-fx*(rbufE(j/2,:,3)+u3(nxf-1,j,:))-fy*(u3(nxf,j+1,:)+u3(nxf,j-1,:))
      u1(nxf,j,:)=Oomehd(qn)*u1(nxf,j,:)+(-fcsq*b1+b2-fc*b3)*delta
      u2(nxf,j,:)=Oomehd(qn)*u2(nxf,j,:)+(fc*b1+ffcD*b2+b3)*delta
      u3(nxf,j,:)=Oomehd(qn)*u3(nxf,j,:)+(-b1-fcD*b2+ffcD*b3)*delta
    enddo
!...local NE corner
    b1=rhs1(nxf,nyf,:)+fxD*(rbufE(nyfd2,:,1)+u1(nxf-1,nyf,:))+fyD*(rbufN(nxfd2,:,1)+u1(nxf,nyf-1,:))
    b2=rhs2(nxf,nyf,:)-fx*(rbufE(nyfd2,:,2)+u2(nxf-1,nyf,:))-fy*(rbufN(nxfd2,:,2)+u2(nxf,nyf-1,:))
    b3=rhs3(nxf,nyf,:)-fx*(rbufE(nyfd2,:,3)+u3(nxf-1,nyf,:))-fy*(rbufN(nxfd2,:,3)+u3(nxf,nyf-1,:))
    u1(nxf,nyf,:)=Oomehd(qn)*u1(nxf,nyf,:)+(-fcsq*b1+b2-fc*b3)*delta
    u2(nxf,nyf,:)=Oomehd(qn)*u2(nxf,nyf,:)+(fc*b1+ffcD*b2+b3)*delta
    u3(nxf,nyf,:)=Oomehd(qn)*u3(nxf,nyf,:)+(-b1-fcD*b2+ffcD*b3)*delta
  endif

  allocate(sbufN(nxfd2,nzL,HDopd2),sbufE(nyfd2,nzL,HDopd2))
  allocate(rbufW(nyfd2,nzL,HDopd2),rbufS(nxfd2,nzL,HDopd2))
!.Send top row to North rank, receive row from South rank
  forall(i=1:nxf-1:2,k=1:nzL) sbufN((i+1)/2,k,:)=(/u1(i,nyf,k),u2(i,nyf,k),u3(i,nyf,k)/)
  call MPI_ISEND(sbufN,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(5),ierr)
!.Receive rows from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(6),ierr)
!.Send right column to East rank, receive column from West rank
  forall(j=1:nyf-1:2,k=1:nzL) sbufE((j+1)/2,k,:)=(/u1(nxf,j,k),u2(nxf,j,k),u3(nxf,j,k)/)
  call MPI_ISEND(sbufE,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(7),ierr)
!.Receive columns from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(8),ierr)

!.Inner points, odd i and odd j
  do j=3,nyf-1,2
    do i=3,nxf-1,2
      b1=rhs1(i,j,:)+fxD*(u1(i+1,j,:)+u1(i-1,j,:))+fyD*(u1(i,j+1,:)+u1(i,j-1,:))
      b2=rhs2(i,j,:)-fx*(u2(i+1,j,:)+u2(i-1,j,:))-fy*(u2(i,j+1,:)+u2(i,j-1,:))
      b3=rhs3(i,j,:)-fx*(u3(i+1,j,:)+u3(i-1,j,:))-fy*(u3(i,j+1,:)+u3(i,j-1,:))
      u1(i,j,:)=Oomehd(qn)*u1(i,j,:)+(-fcsq*b1+b2-fc*b3)*delta
      u2(i,j,:)=Oomehd(qn)*u2(i,j,:)+(fc*b1+ffcD*b2+b3)*delta
      u3(i,j,:)=Oomehd(qn)*u3(i,j,:)+(-b1-fcD*b2+ffcD*b3)*delta
    enddo
  enddo
!.Bottom boundary, odd i
  call MPI_WAIT(req(6),stat(:,6),ierr)
  do i=3,nxf-1,2
    b1=rhs1(i,1,:)+fxD*(u1(i+1,1,:)+u1(i-1,1,:))+fyD*(u1(i,2,:)+rbufS((i+1)/2,:,1))
    b2=rhs2(i,1,:)-fx*(u2(i+1,1,:)+u2(i-1,1,:))-fy*(u2(i,2,:)+rbufS((i+1)/2,:,2))
    b3=rhs3(i,1,:)-fx*(u3(i+1,1,:)+u3(i-1,1,:))-fy*(u3(i,2,:)+rbufS((i+1)/2,:,3))
    u1(i,1,:)=Oomehd(qn)*u1(i,1,:)+(-fcsq*b1+b2-fc*b3)*delta
    u2(i,1,:)=Oomehd(qn)*u2(i,1,:)+(fc*b1+ffcD*b2+b3)*delta
    u3(i,1,:)=Oomehd(qn)*u3(i,1,:)+(-b1-fcD*b2+ffcD*b3)*delta
  enddo
!.Left boundary, odd j
  call MPI_WAIT(req(8),stat(:,8),ierr)
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
!...global SW corner
    b1=rhs1(1,1,:)+fxD*(u1(2,1,:)+Ls1)+fyD*(u1(1,2,:)+rbufS(1,:,1))
    b2=rhs2(1,1,:)-fx*(u2(2,1,:)+Ls2)-fy*(u2(1,2,:)+rbufS(1,:,2))
    b3=rhs3(1,1,:)-fx*(u3(2,1,:)+Ls3)-fy*(u3(1,2,:)+rbufS(1,:,3))
    u1(1,1,:)=Oomehd(qn)*u1(1,1,:)+(-flsq*b1+b2+fl*b3)*deltaL
    u2(1,1,:)=Oomehd(qn)*u2(1,1,:)+(-fl*b1+fflD*b2+b3)*deltaL
    u3(1,1,:)=Oomehd(qn)*u3(1,1,:)+(-b1+flD*b2+fflD*b3)*deltaL
    do j=3,nyf-1,2
      b1=rhs1(1,j,:)+fxD*(u1(2,j,:)+Ls1)+fyD*(u1(1,j+1,:)+u1(1,j-1,:))
      b2=rhs2(1,j,:)-fx*(u2(2,j,:)+Ls2)-fy*(u2(1,j+1,:)+u2(1,j-1,:))
      b3=rhs3(1,j,:)-fx*(u3(2,j,:)+Ls3)-fy*(u3(1,j+1,:)+u3(1,j-1,:))
      u1(1,j,:)=Oomehd(qn)*u1(1,j,:)+(-flsq*b1+b2+fl*b3)*deltaL
      u2(1,j,:)=Oomehd(qn)*u2(1,j,:)+(-fl*b1+fflD*b2+b3)*deltaL
      u3(1,j,:)=Oomehd(qn)*u3(1,j,:)+(-b1+flD*b2+fflD*b3)*deltaL
    enddo
  else
!...local SW corner
    b1=rhs1(1,1,:)+fxD*(u1(2,1,:)+rbufW(1,:,1))+fyD*(u1(1,2,:)+rbufS(1,:,1))
    b2=rhs2(1,1,:)-fx*(u2(2,1,:)+rbufW(1,:,2))-fy*(u2(1,2,:)+rbufS(1,:,2))
    b3=rhs3(1,1,:)-fx*(u3(2,1,:)+rbufW(1,:,3))-fy*(u3(1,2,:)+rbufS(1,:,3))
    u1(1,1,:)=Oomehd(qn)*u1(1,1,:)+(-fcsq*b1+b2-fc*b3)*delta
    u2(1,1,:)=Oomehd(qn)*u2(1,1,:)+(fc*b1+ffcD*b2+b3)*delta
    u3(1,1,:)=Oomehd(qn)*u3(1,1,:)+(-b1-fcD*b2+ffcD*b3)*delta
    do j=3,nyf-1,2
      b1=rhs1(1,j,:)+fxD*(u1(2,j,:)+rbufW((j+1)/2,:,1))+fyD*(u1(1,j+1,:)+u1(1,j-1,:))
      b2=rhs2(1,j,:)-fx*(u2(2,j,:)+rbufW((j+1)/2,:,2))-fy*(u2(1,j+1,:)+u2(1,j-1,:))
      b3=rhs3(1,j,:)-fx*(u3(2,j,:)+rbufW((j+1)/2,:,3))-fy*(u3(1,j+1,:)+u3(1,j-1,:))
      u1(1,j,:)=Oomehd(qn)*u1(1,j,:)+(-fcsq*b1+b2-fc*b3)*delta
      u2(1,j,:)=Oomehd(qn)*u2(1,j,:)+(fc*b1+ffcD*b2+b3)*delta
      u3(1,j,:)=Oomehd(qn)*u3(1,j,:)+(-b1-fcD*b2+ffcD*b3)*delta
    enddo
  endif

  call MPI_WAITALL(4,(/req(1),req(3),req(5),req(7)/), &
                     (/stat(:,1),stat(:,3),stat(:,5),stat(:,7)/),ierr)

!.Now do even-odd and odd-even squares of the grid, i.e, the black
!.squares of the checker-board
!.Send bottom row to South rank, receive row from North rank
  forall(i=1:nxf-1:2,k=1:nzL) sbufS((i+1)/2,k,:)=(/u1(i,1,k),u2(i,1,k),u3(i,1,k)/)
  call MPI_ISEND(sbufS,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(1),ierr)
!.Receive row from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(2),ierr)
!.Send right column to East rank, receive column from West rank
  forall(j=2:nyf:2,k=1:nzL) sbufE(j/2,k,:)=(/u1(nxf,j,k),u2(nxf,j,k),u3(nxf,j,k)/)
  call MPI_ISEND(sbufE,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(3),ierr)
!.Receive columns from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(4),ierr)
!.Inner points, odd i and even j
  do j=2,nyf-2,2
    do i=3,nxf-1,2
      b1=rhs1(i,j,:)+fxD*(u1(i+1,j,:)+u1(i-1,j,:))+fyD*(u1(i,j+1,:)+u1(i,j-1,:))
      b2=rhs2(i,j,:)-fx*(u2(i+1,j,:)+u2(i-1,j,:))-fy*(u2(i,j+1,:)+u2(i,j-1,:))
      b3=rhs3(i,j,:)-fx*(u3(i+1,j,:)+u3(i-1,j,:))-fy*(u3(i,j+1,:)+u3(i,j-1,:))
      u1(i,j,:)=Oomehd(qn)*u1(i,j,:)+(-fcsq*b1+b2-fc*b3)*delta
      u2(i,j,:)=Oomehd(qn)*u2(i,j,:)+(fc*b1+ffcD*b2+b3)*delta
      u3(i,j,:)=Oomehd(qn)*u3(i,j,:)+(-b1-fcD*b2+ffcD*b3)*delta
    enddo
  enddo
!.Top boundary, odd i
  call MPI_WAIT(req(2),stat(:,2),ierr)
  do i=3,nxf-1,2
    b1=rhs1(i,nyf,:)+fxD*(u1(i+1,nyf,:)+u1(i-1,nyf,:))+fyD*(rbufN((i+1)/2,:,1)+u1(i,nyf-1,:))
    b2=rhs2(i,nyf,:)-fx*(u2(i+1,nyf,:)+u2(i-1,nyf,:))-fy*(rbufN((i+1)/2,:,2)+u2(i,nyf-1,:))
    b3=rhs3(i,nyf,:)-fx*(u3(i+1,nyf,:)+u3(i-1,nyf,:))-fy*(rbufN((i+1)/2,:,3)+u3(i,nyf-1,:))
    u1(i,nyf,:)=Oomehd(qn)*u1(i,nyf,:)+(-fcsq*b1+b2-fc*b3)*delta
    u2(i,nyf,:)=Oomehd(qn)*u2(i,nyf,:)+(fc*b1+ffcD*b2+b3)*delta
    u3(i,nyf,:)=Oomehd(qn)*u3(i,nyf,:)+(-b1-fcD*b2+ffcD*b3)*delta
  enddo
!.Left boundary, even j
  call MPI_WAIT(req(4),stat(:,4),ierr)
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    do j=2,nyf-2,2
      b1=rhs1(1,j,:)+fxD*(u1(2,j,:)+Ls1)+fyD*(u1(1,j+1,:)+u1(1,j-1,:))
      b2=rhs2(1,j,:)-fx*(u2(2,j,:)+Ls2)-fy*(u2(1,j+1,:)+u2(1,j-1,:))
      b3=rhs3(1,j,:)-fx*(u3(2,j,:)+Ls3)-fy*(u3(1,j+1,:)+u3(1,j-1,:))
      u1(1,j,:)=Oomehd(qn)*u1(1,j,:)+(-flsq*b1+b2+fl*b3)*deltaL
      u2(1,j,:)=Oomehd(qn)*u2(1,j,:)+(-fl*b1+fflD*b2+b3)*deltaL
      u3(1,j,:)=Oomehd(qn)*u3(1,j,:)+(-b1+flD*b2+fflD*b3)*deltaL
    enddo
!...global NW corner
    b1=rhs1(1,nyf,:)+fxD*(u1(2,nyf,:)+Ls1)+fyD*(rbufN(1,:,1)+u1(1,nyf-1,:))
    b2=rhs2(1,nyf,:)-fx*(u2(2,nyf,:)+Ls2)-fy*(rbufN(1,:,2)+u2(1,nyf-1,:))
    b3=rhs3(1,nyf,:)-fx*(u3(2,nyf,:)+Ls3)-fy*(rbufN(1,:,3)+u3(1,nyf-1,:))
    u1(1,nyf,:)=Oomehd(qn)*u1(1,nyf,:)+(-flsq*b1+b2+fl*b3)*deltaL
    u2(1,nyf,:)=Oomehd(qn)*u2(1,nyf,:)+(-fl*b1+fflD*b2+b3)*deltaL
    u3(1,nyf,:)=Oomehd(qn)*u3(1,nyf,:)+(-b1+flD*b2+fflD*b3)*deltaL
  else
    do j=2,nyf-2,2
      b1=rhs1(1,j,:)+fxD*(u1(2,j,:)+rbufW(j/2,:,1))+fyD*(u1(1,j+1,:)+u1(1,j-1,:))
      b2=rhs2(1,j,:)-fx*(u2(2,j,:)+rbufW(j/2,:,2))-fy*(u2(1,j+1,:)+u2(1,j-1,:))
      b3=rhs3(1,j,:)-fx*(u3(2,j,:)+rbufW(j/2,:,3))-fy*(u3(1,j+1,:)+u3(1,j-1,:))
      u1(1,j,:)=Oomehd(qn)*u1(1,j,:)+(-fcsq*b1+b2-fc*b3)*delta
      u2(1,j,:)=Oomehd(qn)*u2(1,j,:)+(fc*b1+ffcD*b2+b3)*delta
      u3(1,j,:)=Oomehd(qn)*u3(1,j,:)+(-b1-fcD*b2+ffcD*b3)*delta
    enddo
!...local NW corner
    b1=rhs1(1,nyf,:)+fxD*(u1(2,nyf,:)+rbufW(nyf/2,:,1))+fyD*(rbufN(1,:,1)+u1(1,nyf-1,:))
    b2=rhs2(1,nyf,:)-fx*(u2(2,nyf,:)+rbufW(nyf/2,:,2))-fy*(rbufN(1,:,2)+u2(1,nyf-1,:))
    b3=rhs3(1,nyf,:)-fx*(u3(2,nyf,:)+rbufW(nyf/2,:,3))-fy*(rbufN(1,:,3)+u3(1,nyf-1,:))
    u1(1,nyf,:)=Oomehd(qn)*u1(1,nyf,:)+(-fcsq*b1+b2-fc*b3)*delta
    u2(1,nyf,:)=Oomehd(qn)*u2(1,nyf,:)+(fc*b1+ffcD*b2+b3)*delta
    u3(1,nyf,:)=Oomehd(qn)*u3(1,nyf,:)+(-b1-fcD*b2+ffcD*b3)*delta
  endif

!.Send top row to North rank, receive row from South rank
  forall(i=2:nxf:2,k=1:nzL) sbufN(i/2,k,:)=(/u1(i,nyf,k),u2(i,nyf,k),u3(i,nyf,k)/)
  call MPI_ISEND(sbufN,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(5),ierr)
!.Receive rows from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(6),ierr)
!.Send left column to West rank, receive column from East rank
  forall(j=1:nyf-1:2,k=1:nzL) sbufW((j+1)/2,k,:)=(/u1(1,j,k),u2(1,j,k),u3(1,j,k)/)
  call MPI_ISEND(sbufW,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(7),ierr)
!.Receive column from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,nyfd2*nzL*HDopd2,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(8),ierr)
!.Inner points, even i and odd j
  do j=3,nyf-1,2
    do i=2,nxf-2,2
      b1=rhs1(i,j,:)+fxD*(u1(i+1,j,:)+u1(i-1,j,:))+fyD*(u1(i,j+1,:)+u1(i,j-1,:))
      b2=rhs2(i,j,:)-fx*(u2(i+1,j,:)+u2(i-1,j,:))-fy*(u2(i,j+1,:)+u2(i,j-1,:))
      b3=rhs3(i,j,:)-fx*(u3(i+1,j,:)+u3(i-1,j,:))-fy*(u3(i,j+1,:)+u3(i,j-1,:))
      u1(i,j,:)=Oomehd(qn)*u1(i,j,:)+(-fcsq*b1+b2-fc*b3)*delta
      u2(i,j,:)=Oomehd(qn)*u2(i,j,:)+(fc*b1+ffcD*b2+b3)*delta
      u3(i,j,:)=Oomehd(qn)*u3(i,j,:)+(-b1-fcD*b2+ffcD*b3)*delta
    enddo
  enddo
!.Bottom boundary, even i
  call MPI_WAIT(req(6),stat(:,6),ierr)
  do i=2,nxf-2,2
    b1=rhs1(i,1,:)+fxD*(u1(i+1,1,:)+u1(i-1,1,:))+fyD*(u1(i,2,:)+rbufS(i/2,:,1))
    b2=rhs2(i,1,:)-fx*(u2(i+1,1,:)+u2(i-1,1,:))-fy*(u2(i,2,:)+rbufS(i/2,:,2))
    b3=rhs3(i,1,:)-fx*(u3(i+1,1,:)+u3(i-1,1,:))-fy*(u3(i,2,:)+rbufS(i/2,:,3))
    u1(i,1,:)=Oomehd(qn)*u1(i,1,:)+(-fcsq*b1+b2-fc*b3)*delta
    u2(i,1,:)=Oomehd(qn)*u2(i,1,:)+(fc*b1+ffcD*b2+b3)*delta
    u3(i,1,:)=Oomehd(qn)*u3(i,1,:)+(-b1-fcD*b2+ffcD*b3)*delta
  enddo
!.Right boundary, odd j
  call MPI_WAIT(req(8),stat(:,8),ierr)
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
!...global SE corner
    b1=rhs1(nxf,1,:)+fxD*(Rs1+u1(nxf-1,1,:))+fyD*(u1(nxf,2,:)+rbufS(nxfd2,:,1))
    b2=rhs2(nxf,1,:)-fx*(Rs2+u2(nxf-1,1,:))-fy*(u2(nxf,2,:)+rbufS(nxfd2,:,2))
    b3=rhs3(nxf,1,:)-fx*(Rs3+u3(nxf-1,1,:))-fy*(u3(nxf,2,:)+rbufS(nxfd2,:,3))
    u1(nxf,1,:)=Oomehd(qn)*u1(nxf,1,:)+(-frsq*b1+b2+fr*b3)*deltaR
    u2(nxf,1,:)=Oomehd(qn)*u2(nxf,1,:)+(-fr*b1+ffrD*b2+b3)*deltaR
    u3(nxf,1,:)=Oomehd(qn)*u3(nxf,1,:)+(-b1+frD*b2+ffrD*b3)*deltaR
    do j=3,nyf-1,2
      b1=rhs1(nxf,j,:)+fxD*(Rs1+u1(nxf-1,j,:))+fyD*(u1(nxf,j+1,:)+u1(nxf,j-1,:))
      b2=rhs2(nxf,j,:)-fx*(Rs2+u2(nxf-1,j,:))-fy*(u2(nxf,j+1,:)+u2(nxf,j-1,:))
      b3=rhs3(nxf,j,:)-fx*(Rs3+u3(nxf-1,j,:))-fy*(u3(nxf,j+1,:)+u3(nxf,j-1,:))
      u1(nxf,j,:)=Oomehd(qn)*u1(nxf,j,:)+(-frsq*b1+b2+fr*b3)*deltaR
      u2(nxf,j,:)=Oomehd(qn)*u2(nxf,j,:)+(-fr*b1+ffrD*b2+b3)*deltaR
      u3(nxf,j,:)=Oomehd(qn)*u3(nxf,j,:)+(-b1+frD*b2+ffrD*b3)*deltaR
    enddo
  else
!...local SE corner
    b1=rhs1(nxf,1,:)+fxD*(rbufE(1,:,1)+u1(nxf-1,1,:))+fyD*(u1(nxf,2,:)+rbufS(nxfd2,:,1))
    b2=rhs2(nxf,1,:)-fx*(rbufE(1,:,2)+u2(nxf-1,1,:))-fy*(u2(nxf,2,:)+rbufS(nxfd2,:,2))
    b3=rhs3(nxf,1,:)-fx*(rbufE(1,:,3)+u3(nxf-1,1,:))-fy*(u3(nxf,2,:)+rbufS(nxfd2,:,3))
    u1(nxf,1,:)=Oomehd(qn)*u1(nxf,1,:)+(-fcsq*b1+b2-fc*b3)*delta
    u2(nxf,1,:)=Oomehd(qn)*u2(nxf,1,:)+(fc*b1+ffcD*b2+b3)*delta
    u3(nxf,1,:)=Oomehd(qn)*u3(nxf,1,:)+(-b1-fcD*b2+ffcD*b3)*delta
    do j=3,nyf-1,2
      b1=rhs1(nxf,j,:)+fxD*(rbufE((j+1)/2,:,1)+u1(nxf-1,j,:))+fyD*(u1(nxf,j+1,:)+u1(nxf,j-1,:))
      b2=rhs2(nxf,j,:)-fx*(rbufE((j+1)/2,:,2)+u2(nxf-1,j,:))-fy*(u2(nxf,j+1,:)+u2(nxf,j-1,:))
      b3=rhs3(nxf,j,:)-fx*(rbufE((j+1)/2,:,3)+u3(nxf-1,j,:))-fy*(u3(nxf,j+1,:)+u3(nxf,j-1,:))
      u1(nxf,j,:)=Oomehd(qn)*u1(nxf,j,:)+(-fcsq*b1+b2-fc*b3)*delta
      u2(nxf,j,:)=Oomehd(qn)*u2(nxf,j,:)+(fc*b1+ffcD*b2+b3)*delta
      u3(nxf,j,:)=Oomehd(qn)*u3(nxf,j,:)+(-b1-fcD*b2+ffcD*b3)*delta
    enddo
  endif

  call MPI_WAITALL(4,(/req(1),req(3),req(5),req(7)/), &
                     (/stat(:,1),stat(:,3),stat(:,5),stat(:,7)/),ierr)
  end subroutine relaxhd
!*********************************************************
  subroutine reshd(ng,res1,res2,res3,u1,u2,u3,rhs1,rhs2,rhs3,xBC,qn)
!.Returns minus the residual of each equation in the split approach
!.to solving the 6th order hyperdiffusion problem
!.      res1 = rhs1 - (u3 - Diff*L u1)
!.      res2 = rhs2 - (L u2 - u1) 
!.      res3 = rhs3 - (L u3 - u2) 
!.(L=Laplacian) using a 2nd order accurate finite difference stencil.
!.INPUTS
!. integer ng: current grid
!. real(dp) dimension(:,:,:) u1, u2, u3: current solution at grid ng
!. real(dp) dimension(:,:,:) rhs1, rhs2, rhs3: the right-hand side of Poisson equation
!. integer(2) xBC: boundary conditions along x
!. integer qn: quantity being hyperdiffused (1 to nqhd)
!.OUTPUTS
!. real(dp) dimension(:,:,:) res1, res2, res3: negated residual
  implicit none
  include 'mpif.h'
  real(DP), intent(out) :: res1(:,:,:),res2(:,:,:),res3(:,:,:)
  real(DP), intent(in) :: u1(:,:,:),u2(:,:,:),u3(:,:,:), &
                          rhs1(:,:,:),rhs2(:,:,:),rhs3(:,:,:)
  integer(I4B), intent(in) :: ng,xBC(:),qn
  integer(I4B) :: i,j,k,Nrank,Srank,Erank,Wrank
  real(DP) :: fx,fy,fc,fxD,fyD,fcD,Lq,Rq
  real(DP), allocatable :: Ls1(:),Ls2(:),Ls3(:),Rs1(:),Rs2(:),Rs3(:),loc1D(:)
  real(DP), allocatable, dimension(:,:,:) :: sbufW,rbufE,sbufE,rbufW, &
                                             sbufN,rbufS,sbufS,rbufN
  integer(I4B) :: req(8),stat(MPI_STATUS_SIZE,8),ierr

#IF (MGTIMER==1)
    call CPU_TIME(rest1)
#ENDIF

!.Neigboring ranks to communicate with
  Nrank =nborhd(qn)%a(ng,1)
  Erank =nborhd(qn)%a(ng,3)
  Srank =nborhd(qn)%a(ng,5)
  Wrank =nborhd(qn)%a(ng,7)

!.Boundary conditions are implemented with some amount of flexibility.
!.The right boundary for example, if u(nxf+1) is needed, it will be
!.written as
!.  u(nxf+1) = Rq*u(nxf) + Rr*u(nxf-1) + Rs
!.so that symmetry, extrapolation or Dirichlet conditions can be
!.easily imposed.
!......... BCs of u ................................!
  allocate(Ls1(nzL),Rs1(nzL),loc1D(nzL))
  allocate(Ls2(nzL),Rs2(nzL),Ls3(nzL),Rs3(nzL))
!.LHS BC of u
  if (xBC(1)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Lq  = -1.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u1(1,:,k))
    call MPI_ALLreduce(loc1D,Ls1,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Ls1 =  2.0_dp*Ls1/dble(nyf*jprocshd(qn)%a(ng))
    forall(k=1:nzL) loc1D(k)=sum(u2(1,:,k))
    call MPI_ALLreduce(loc1D,Ls2,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Ls2 =  2.0_dp*Ls2/dble(nyf*jprocshd(qn)%a(ng))
    forall(k=1:nzL) loc1D(k)=sum(u3(1,:,k))
    call MPI_ALLreduce(loc1D,Ls3,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Ls3 =  2.0_dp*Ls3/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC==1 .OR. xBC(1)==-1) then ! symmetry
    Lq  = dble(xBC(1))
    Ls1 = 0.0_dp
    Ls2 = 0.0_dp
    Ls3 = 0.0_dp
  endif
!.RHS BC of u
  if (xBC(2)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Rq  = -1.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u1(nxf,:,k))
    call MPI_ALLreduce(loc1D,Rs1,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Rs1 =  2.0_dp*Rs1/dble(nyf*jprocshd(qn)%a(ng))
    forall(k=1:nzL) loc1D(k)=sum(u2(nxf,:,k))
    call MPI_ALLreduce(loc1D,Rs2,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Rs2 =  2.0_dp*Rs2/dble(nyf*jprocshd(qn)%a(ng))
    forall(k=1:nzL) loc1D(k)=sum(u3(nxf,:,k))
    call MPI_ALLreduce(loc1D,Rs3,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Rs3 =  2.0_dp*Rs3/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC(2)==1 .OR. xBC(2)==-1) then ! symmetry
    Rq  = dble(xBC(2))
    Rs1 = 0.0_dp
    Rs2 = 0.0_dp
    Rs3 = 0.0_dp
  endif
!........................................................!

  allocate(sbufW(nyf,nzL,HDopd2),sbufS(nxf,nzL,HDopd2), &
           sbufN(nxf,nzL,HDopd2),sbufE(nyf,nzL,HDopd2))
  allocate(rbufW(nyf,nzL,HDopd2),rbufS(nxf,nzL,HDopd2), &
           rbufN(nxf,nzL,HDopd2),rbufE(nyf,nzL,HDopd2))
!.Send top row to North rank, receive row from South rank
  forall(i=1:nxf,k=1:nzL) sbufN(i,k,:)=(/u1(i,nyf,k),u2(i,nyf,k),u3(i,nyf,k)/)
  call MPI_ISEND(sbufN,nxf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(1),ierr)
!.Send right column to East rank, receive column from West rank
  forall(j=1:nyf,k=1:nzL) sbufE(j,k,:)=(/u1(nxf,j,k),u2(nxf,j,k),u3(nxf,j,k)/)
  call MPI_ISEND(sbufE,nyf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(2),ierr)
!.Send bottom row to South rank, receive row from North rank
  forall(i=1:nxf,k=1:nzL) sbufS(i,k,:)=(/u1(i,1,k),u2(i,1,k),u3(i,1,k)/)
  call MPI_ISEND(sbufS,nxf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(3),ierr)
!.Send left column to West rank, receive column from East rank
  forall(j=1:nyf,k=1:nzL) sbufW(j,k,:)=(/u1(1,j,k),u2(1,j,k),u3(1,j,k)/)
  call MPI_ISEND(sbufW,nyf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(4),ierr)
!.Receive rows from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(5),ierr)
!.Receive columns from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,nyf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(6),ierr)
!.Receive row from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(7),ierr)
!.Receive column from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,nyf*nzL*HDopd2,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(8),ierr)

  fx  = (dble(nxf*iprocshd(qn)%a(ng))/bLx)**2         !.1/(hx^2)
  fy  = (dble(nyf*jprocshd(qn)%a(ng))/bLy)**2         !.1/(hy^2)
  fc  = 2.0_dp*(fx+fy)
  fxD = hDiff(1)*fx
  fyD = hDiff(2)*fy
!  fcD = Diff*fc
  fcD = 2.0_dp*(hDiff(1)*fx+hDiff(2)*fy)

!.Interior points
  do j=2,nyf-1
    do i=2,nxf-1
      res1(i,j,:)=rhs1(i,j,:)-(u3(i,j,:)-fxD*(u1(i+1,j,:)+u1(i-1,j,:)) &
        -fyD*(u1(i,j+1,:)+u1(i,j-1,:))+fcD*u1(i,j,:))
      res2(i,j,:)=rhs2(i,j,:)-(fx*(u2(i+1,j,:)+u2(i-1,j,:)) &
        +fy*(u2(i,j+1,:)+u2(i,j-1,:))-fc*u2(i,j,:)-u1(i,j,:))
      res3(i,j,:)=rhs3(i,j,:)-(fx*(u3(i+1,j,:)+u3(i-1,j,:)) &
        +fy*(u3(i,j+1,:)+u3(i,j-1,:))-fc*u3(i,j,:)-u2(i,j,:))
    enddo
  enddo

!.Boundary points
  call MPI_WAIT(req(5),stat(:,5),ierr)
  do i=2,nxf-1
    res1(i,1,:)=rhs1(i,1,:)-(u3(i,1,:)-fxD*(u1(i+1,1,:)+u1(i-1,1,:)) &
      -fyD*(u1(i,2,:)+rbufS(i,:,1))+fcD*u1(i,1,:))
    res2(i,1,:)=rhs2(i,1,:)-(fx*(u2(i+1,1,:)+u2(i-1,1,:)) &
      +fy*(u2(i,2,:)+rbufS(i,:,2))-fc*u2(i,1,:)-u1(i,1,:))
    res3(i,1,:)=rhs3(i,1,:)-(fx*(u3(i+1,1,:)+u3(i-1,1,:)) &
      +fy*(u3(i,2,:)+rbufS(i,:,3))-fc*u3(i,1,:)-u2(i,1,:))
  enddo
!.Left boundary
  call MPI_WAIT(req(6),stat(:,6),ierr)
  if (iIDhd(qn)%a(ng) == 0) then !.Left most process (along x)
!...global SW corner
    res1(1,1,:)=rhs1(1,1,:)-(u3(1,1,:)-fxD*(u1(2,1,:)+Lq*u1(1,1,:)+Ls1) &
      -fyD*(u1(1,2,:)+rbufS(1,:,1))+fcD*u1(1,1,:))
    res2(1,1,:)=rhs2(1,1,:)-(fx*(u2(2,1,:)+Lq*u2(1,1,:)+Ls2) &
      +fy*(u2(1,2,:)+rbufS(1,:,2))-fc*u2(1,1,:)-u1(1,1,:))
    res3(1,1,:)=rhs3(1,1,:)-(fx*(u3(2,1,:)+Lq*u3(1,1,:)+Ls3) &
      +fy*(u3(1,2,:)+rbufS(1,:,3))-fc*u3(1,1,:)-u2(1,1,:))
    do j=2,nyf-1
      res1(1,j,:)=rhs1(1,j,:)-(u3(1,j,:)-fxD*(u1(2,j,:)+Lq*u1(1,j,:)+Ls1) &
        -fyD*(u1(1,j+1,:)+u1(1,j-1,:))+fcD*u1(1,j,:))
      res2(1,j,:)=rhs2(1,j,:)-(fx*(u2(2,j,:)+Lq*u2(1,j,:)+Ls2) &
        +fy*(u2(1,j+1,:)+u2(1,j-1,:))-fc*u2(1,j,:)-u1(1,j,:))
      res3(1,j,:)=rhs3(1,j,:)-(fx*(u3(2,j,:)+Lq*u3(1,j,:)+Ls3) &
        +fy*(u3(1,j+1,:)+u3(1,j-1,:))-fc*u3(1,j,:)-u2(1,j,:))
    enddo
  else
!...local SW corner
    res1(1,1,:)=rhs1(1,1,:)-(u3(1,1,:)-fxD*(u1(2,1,:)+rbufW(1,:,1)) &
      -fyD*(u1(1,2,:)+rbufS(1,:,1))+fcD*u1(1,1,:))
    res2(1,1,:)=rhs2(1,1,:)-(fx*(u2(2,1,:)+rbufW(1,:,2)) &
      +fy*(u2(1,2,:)+rbufS(1,:,2))-fc*u2(1,1,:)-u1(1,1,:))
    res3(1,1,:)=rhs3(1,1,:)-(fx*(u3(2,1,:)+rbufW(1,:,3)) &
      +fy*(u3(1,2,:)+rbufS(1,:,3))-fc*u3(1,1,:)-u2(1,1,:))
    do j=2,nyf-1
      res1(1,j,:)=rhs1(1,j,:)-(u3(1,j,:)-fxD*(u1(2,j,:)+rbufW(j,:,1)) &
        -fyD*(u1(1,j+1,:)+u1(1,j-1,:))+fcD*u1(1,j,:))
      res2(1,j,:)=rhs2(1,j,:)-(fx*(u2(2,j,:)+rbufW(j,:,2)) &
        +fy*(u2(1,j+1,:)+u2(1,j-1,:))-fc*u2(1,j,:)-u1(1,j,:))
      res3(1,j,:)=rhs3(1,j,:)-(fx*(u3(2,j,:)+rbufW(j,:,3)) &
        +fy*(u3(1,j+1,:)+u3(1,j-1,:))-fc*u3(1,j,:)-u2(1,j,:))
    enddo
  endif

  call MPI_WAIT(req(7),stat(:,7),ierr)
  if (iIDhd(qn)%a(ng) == 0) then !.Left most process (along x)
!...global NW corner
    res1(1,nyf,:)=rhs1(1,nyf,:)-(u3(1,nyf,:)-fxD*(u1(2,nyf,:)+Lq*u1(1,nyf,:)+Ls1) &
      -fyD*(rbufN(1,:,1)+u1(1,nyf-1,:))+fcD*u1(1,nyf,:))
    res2(1,nyf,:)=rhs2(1,nyf,:)-(fx*(u2(2,nyf,:)+Lq*u2(1,nyf,:)+Ls2) &
      +fy*(rbufN(1,:,2)+u2(1,nyf-1,:))-fc*u2(1,nyf,:)-u1(1,nyf,:))
    res3(1,nyf,:)=rhs3(1,nyf,:)-(fx*(u3(2,nyf,:)+Lq*u3(1,nyf,:)+Ls3) &
      +fy*(rbufN(1,:,3)+u3(1,nyf-1,:))-fc*u3(1,nyf,:)-u2(1,nyf,:))
  else
!...local NW corner
    res1(1,nyf,:)=rhs1(1,nyf,:)-(u3(1,nyf,:)-fxD*(u1(2,nyf,:)+rbufW(nyf,:,1)) &
      -fyD*(rbufN(1,:,1)+u1(1,nyf-1,:))+fcD*u1(1,nyf,:))
    res2(1,nyf,:)=rhs2(1,nyf,:)-(fx*(u2(2,nyf,:)+rbufW(nyf,:,2)) &
      +fy*(rbufN(1,:,2)+u2(1,nyf-1,:))-fc*u2(1,nyf,:)-u1(1,nyf,:))
    res3(1,nyf,:)=rhs3(1,nyf,:)-(fx*(u3(2,nyf,:)+rbufW(nyf,:,3)) &
      +fy*(rbufN(1,:,3)+u3(1,nyf-1,:))-fc*u3(1,nyf,:)-u2(1,nyf,:))
  endif
!.Top boundary
  do i=2,nxf-1
    res1(i,nyf,:)=rhs1(i,nyf,:)-(u3(i,nyf,:)-fxD*(u1(i+1,nyf,:)+u1(i-1,nyf,:)) &
      -fyD*(rbufN(i,:,1)+u1(i,nyf-1,:))+fcD*u1(i,nyf,:))
    res2(i,nyf,:)=rhs2(i,nyf,:)-(fx*(u2(i+1,nyf,:)+u2(i-1,nyf,:)) &
      +fy*(rbufN(i,:,2)+u2(i,nyf-1,:))-fc*u2(i,nyf,:)-u1(i,nyf,:))
    res3(i,nyf,:)=rhs3(i,nyf,:)-(fx*(u3(i+1,nyf,:)+u3(i-1,nyf,:)) &
      +fy*(rbufN(i,:,3)+u3(i,nyf-1,:))-fc*u3(i,nyf,:)-u2(i,nyf,:))
  enddo

!.Right boundary
  call MPI_WAIT(req(8),stat(:,8),ierr)
  if (iIDhd(qn)%a(ng) == iprocsm1) then !.Right most process (along x)
!...global NE corner
    res1(nxf,nyf,:)=rhs1(nxf,nyf,:)-(u3(nxf,nyf,:)-fxD*(Rq*u1(nxf,nyf,:)+Rs1 &
      +u1(nxf-1,nyf,:))-fyD*(rbufN(nxf,:,1)+u1(nxf,nyf-1,:))+fcD*u1(nxf,nyf,:))
    res2(nxf,nyf,:)=rhs2(nxf,nyf,:)-(fx*(Rq*u2(nxf,nyf,:)+Rs2+u2(nxf-1,nyf,:)) &
      +fy*(rbufN(nxf,:,2)+u2(nxf,nyf-1,:))-fc*u2(nxf,nyf,:)-u1(nxf,nyf,:))
    res3(nxf,nyf,:)=rhs3(nxf,nyf,:)-(fx*(Rq*u3(nxf,nyf,:)+Rs3+u3(nxf-1,nyf,:)) &
      +fy*(rbufN(nxf,:,3)+u3(nxf,nyf-1,:))-fc*u3(nxf,nyf,:)-u2(nxf,nyf,:))
!...global SE corner
    res1(nxf,1,:)=rhs1(nxf,1,:)-(u3(nxf,1,:)-fxD*(Rq*u1(nxf,1,:)+Rs1+u1(nxf-1,1,:)) &
      -fyD*(u1(nxf,2,:)+rbufS(nxf,:,1))+fcD*u1(nxf,1,:))
    res2(nxf,1,:)=rhs2(nxf,1,:)-(fx*(Rq*u2(nxf,1,:)+Rs2+u2(nxf-1,1,:)) &
      +fy*(u2(nxf,2,:)+rbufS(nxf,:,2))-fc*u2(nxf,1,:)-u1(nxf,1,:))
    res3(nxf,1,:)=rhs3(nxf,1,:)-(fx*(Rq*u3(nxf,1,:)+Rs3+u3(nxf-1,1,:)) &
      +fy*(u3(nxf,2,:)+rbufS(nxf,:,3))-fc*u3(nxf,1,:)-u2(nxf,1,:))
    do j=2,nyf-1
      res1(nxf,j,:)=rhs1(nxf,j,:)-(u3(nxf,j,:)-fxD*(Rq*u1(nxf,j,:)+Rs1+u1(nxf-1,j,:)) &
        -fyD*(u1(nxf,j+1,:)+u1(nxf,j-1,:))+fcD*u1(nxf,j,:))
      res2(nxf,j,:)=rhs2(nxf,j,:)-(fx*(Rq*u2(nxf,j,:)+Rs2+u2(nxf-1,j,:)) &
        +fy*(u2(nxf,j+1,:)+u2(nxf,j-1,:))-fc*u2(nxf,j,:)-u1(nxf,j,:))
      res3(nxf,j,:)=rhs3(nxf,j,:)-(fx*(Rq*u3(nxf,j,:)+Rs3+u3(nxf-1,j,:)) &
        +fy*(u3(nxf,j+1,:)+u3(nxf,j-1,:))-fc*u3(nxf,j,:)-u2(nxf,j,:))
    enddo
  else
!...local NE corner
    res1(nxf,nyf,:)=rhs1(nxf,nyf,:)-(u3(nxf,nyf,:)-fxD*(rbufE(nyf,:,1)+u1(nxf-1,nyf,:)) &
      -fyD*(rbufN(nxf,:,1)+u1(nxf,nyf-1,:))+fcD*u1(nxf,nyf,:))
    res2(nxf,nyf,:)=rhs2(nxf,nyf,:)-(fx*(rbufE(nyf,:,2)+u2(nxf-1,nyf,:)) &
      +fy*(rbufN(nxf,:,2)+u2(nxf,nyf-1,:))-fc*u2(nxf,nyf,:)-u1(nxf,nyf,:))
    res3(nxf,nyf,:)=rhs3(nxf,nyf,:)-(fx*(rbufE(nyf,:,3)+u3(nxf-1,nyf,:)) &
      +fy*(rbufN(nxf,:,3)+u3(nxf,nyf-1,:))-fc*u3(nxf,nyf,:)-u2(nxf,nyf,:))
!...local SE corner
    res1(nxf,1,:)=rhs1(nxf,1,:)-(u3(nxf,1,:)-fxD*(rbufE(1,:,1)+u1(nxf-1,1,:)) &
      -fyD*(u1(nxf,2,:)+rbufS(nxf,:,1))+fcD*u1(nxf,1,:))
    res2(nxf,1,:)=rhs2(nxf,1,:)-(fx*(rbufE(1,:,2)+u2(nxf-1,1,:)) &
      +fy*(u2(nxf,2,:)+rbufS(nxf,:,2))-fc*u2(nxf,1,:)-u1(nxf,1,:))
    res3(nxf,1,:)=rhs3(nxf,1,:)-(fx*(rbufE(1,:,3)+u3(nxf-1,1,:)) &
      +fy*(u3(nxf,2,:)+rbufS(nxf,:,3))-fc*u3(nxf,1,:)-u2(nxf,1,:))
    do j=2,nyf-1
      res1(nxf,j,:)=rhs1(nxf,j,:)-(u3(nxf,j,:)-fxD*(rbufE(j,:,1)+u1(nxf-1,j,:)) &
        -fyD*(u1(nxf,j+1,:)+u1(nxf,j-1,:))+fcD*u1(nxf,j,:))
      res2(nxf,j,:)=rhs2(nxf,j,:)-(fx*(rbufE(j,:,2)+u2(nxf-1,j,:)) &
        +fy*(u2(nxf,j+1,:)+u2(nxf,j-1,:))-fc*u2(nxf,j,:)-u1(nxf,j,:))
      res3(nxf,j,:)=rhs3(nxf,j,:)-(fx*(rbufE(j,:,3)+u3(nxf-1,j,:)) &
        +fy*(u3(nxf,j+1,:)+u3(nxf,j-1,:))-fc*u3(nxf,j,:)-u2(nxf,j,:))
    enddo
  endif
  call MPI_WAITALL(4,(/req(1),req(2),req(3),req(4)/), &
                     (/stat(:,1),stat(:,2),stat(:,3),stat(:,4)/),ierr)

#IF (MGTIMER==1)
    call CPU_TIME(rest2)
    resthd=resthd+rest2-rest1
#ENDIF
  end subroutine reshd
!*********************************************************
#ELSE
!.If (HDop==2 .OR. HDop==3 .OR. HDop==8 .OR. HDop==12)
  recursive subroutine Gcychd(ng,u,rhs,xBC,qn)
!.Gamma cycles to solve equation from implicit hyperdiffusion
!.equation arising from Laplacian based 4th order diffusion
!.(HDop==2)
!.                   u + Diff*(L^2)u=rho,
!.or Laplacian-based 6th order diffusion (HDop==3)
!.                   u - Diff*(L^3)u=rho,
!.or anisotropic without cross terms (ANISOWC) 4th order
!.diffusion (HDop==8)
!.    u + ( Diffx*(d^4/dx^4) + Diffy*(d^4/dy^4) )u=rho
!.or ANISOWC 6th order diffusion (HDop==12)
!.    u - ( Diffx*(d^6/dx^6) + Diffy*(d^6/dy^6) )u=rho
!.using the direct application of their finite difference stencil.
!.INPUTS
!. integer ng: current grid
!. real(dp) dimension(:,:,:) u: current solution at grid ng
!. real(dp) dimension(:,:,:) rhs: the right-hand side rho vector
!. integer(2) xBC: boundary conditions along x
!. integer qn: which quantity is being hyperdiffused (1 to nqhd)
!.OUTPUTS
!. real(dp) dimension(:,:,:) u: updated solution at grid ng
  integer(I4B), intent(in) :: ng, xBC(:), qn
  real(DP), intent(inout) :: u(:,:,:)
  real(DP), intent(in) :: rhs(:,:,:)
  integer(I4B) :: jpost, jpre, i
  real(DP), allocatable :: res(:,:,:),v(:,:,:)

  iprocsm1=iprocshd(qn)%a(ng)-1 !.ID of last rank in x COMM
  jprocsm1=jprocshd(qn)%a(ng)-1 !.ID of last rank in y COMM

  if (ng==1) then        !.Bottom of gamma cycle: Solve on coarsest grid
#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpost=1,NU2hd(qn)
      call relaxhd(ng,u,rhs,xBC,qn)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    postXthd=postXthd+count2-count1
#ENDIF
  else                  !.On downward stroke towards coarser grids
#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpre=1,NU1hd(qn)      !.Pre-smoothing
      call relaxhd(ng,u,rhs,xBC,qn)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    preXthd=preXthd+count2-count1
#ENDIF

!...Number of grid points in coarse grid
    nxc=nxChd(qn)%a(ng)
    nyc=nyChd(qn)%a(ng)
!...Number of grid points in coarse grid after subdomain unification
    nxu=nxFhd(qn)%a(ng-1)
    nyu=nyFhd(qn)%a(ng-1)
    allocate(res(nxu,nyu,nzL),v(nxu,nyu,nzL))

#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
!...Restricted residual is the next RHS
    res=rstrct(reshd(ng,u,rhs,xBC,qn),wdchd(qn)%a(ng),xBC, &
      acuhd(qn)%a(ng-1,:),urnkhd(qn)%a(ng-1)%a,nurnkhd(qn)%a(ng-1,:))
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    rNrthd=rNrthd+count2-count1
#ENDIF

!...Next fine grid is the coarse, or unified subdomain, grid
    nxf=nxu
    nyf=nyu

    v=0.0_dp !.Zero for initial guess in next relaxation
    if (acuhd(qn)%a(ng-1,1)) then
!.....Only active processes proceed to next coarser grid
      do i=1,GAMhd(qn)
        call Gcychd(ng-1,v,res,xBC,qn)  !.Recursive call for the coarse grid correction
      enddo
    endif

!...Redefine quantities that were modified in coarser grids
    iprocsm1=iprocshd(qn)%a(ng)-1
    jprocsm1=jprocshd(qn)%a(ng)-1
    nxc=nxChd(qn)%a(ng)
    nyc=nyChd(qn)%a(ng)
    nxf=nxFhd(qn)%a(ng)
    nyf=nyFhd(qn)%a(ng)

#IF (MGTIMER==1)
    call CPU_TIME(prolt1)
#ENDIF
!...Upward stroke -> finer grid: prolong coarse grid error + correct solution
    u=u+prolng(v,wdchd(qn)%a(ng),xBC,COMMhd(qn)%a(ng,:), &
      iIDhd(qn)%a(ng),jprocshd(qn)%a(ng),nborhd(qn)%a(ng,:), &
      acuhd(qn)%a(ng-1,:),urnkhd(qn)%a(ng-1)%a,nurnkhd(qn)%a(ng-1,:))
#IF (MGTIMER==1)
    call CPU_TIME(prolt2)
    prolthd=prolthd+prolt2-prolt1
#ENDIF

#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpost=1,NU2hd(qn)      !.Post-smoothing
      call relaxhd(ng,u,rhs,xBC,qn)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    postXthd=postXthd+count2-count1
#ENDIF
  endif

  end subroutine Gcychd
!*********************************************************
#IF (HDop==2 .OR. HDop==3)
!*********************************************************
  subroutine relaxhd(ng,u,rhs,xBC,qn)
!.Local damped Gauss-Seidel relaxation of Laplacian-based
!.4th order hyperdiffusion equation (HDop==2)
!.                   u + Diff*(L^2)u=rho,
!.or Laplacian-based 6th order diffusion (HDop==3)
!.                   u - Diff*(L^3)u=rho,
!.directly using its finite difference stencil.
!.INPUTS
!. integer ng: current grid
!. real(dp) dimension(:,:,:) u: current solution at grid ng
!. real(dp) dimension(:,:,:) rhs: the right-hand side vector
!. integer(2) xBC: boundary conditions along x
!. integer qn: quantity being hyperdiffused (1 to nqhd)
!.OUTPUTS
!. real(dp) dimension(:,:,:) u: updated solution at grid ng
  implicit none
  include 'mpif.h'
  real(DP), intent(inout) :: u(:,:,:)
  real(DP), intent(in) :: rhs(:,:,:)
  integer(I4B), intent(in) :: ng,xBC(:),qn
  integer(I4B) :: i,j,k,Nrank,Srank,Erank,Wrank,NErank,SErank,SWrank,NWrank
  integer(I4B) :: ierr,req(16),stat(MPI_STATUS_SIZE,16)
  real(DP), allocatable :: Ls(:,:),Rs(:,:),loc2D(:,:)
  real(DP), allocatable, dimension(:,:,:) :: sbufN,sbufE,sbufS,sbufW, &
                                             rbufN,rbufE,rbufS,rbufW
  real(DP), allocatable, dimension(:,:) :: sbufNE,sbufSE,sbufSW,sbufNW, &
                                           rbufNE,rbufSE,rbufSW,rbufNW
  real(DP) :: rdx,rdy,fx,fxx,fy,fyy,fxy,fc,delta,Lq,Rq
#IF (HDop==3)
  real(DP) :: fx6,fx4,fy6,fy4,fxxy,fxyy
#ENDIF
  real(DP), allocatable :: deltaL(:),deltaR(:)

!.Neigboring ranks to communicate with
  Nrank  = nborhd(qn)%a(ng,1)
  NErank = nborhd(qn)%a(ng,2)
  Erank  = nborhd(qn)%a(ng,3)
  SErank = nborhd(qn)%a(ng,4)
  Srank  = nborhd(qn)%a(ng,5)
  SWrank = nborhd(qn)%a(ng,6)
  Wrank  = nborhd(qn)%a(ng,7)
  NWrank = nborhd(qn)%a(ng,8)

!.Boundary conditions are implemented with some amount of flexibility.
!.The right boundary for example, if u(nxf+1) is needed, it will be
!.written as
!.  u(nxf+1) = Rq*u(nxf) + Rr*u(nxf-1) + Rs
!.so that symmetry, extrapolation or Dirichlet conditions can be
!.easily imposed.
!......... BCs of u ................................!
  allocate(Ls(HDop,nzL),Rs(HDop,nzL),loc2D(HDop,nzL))
!.LHS BC of u
  if (xBC(1)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Lq = -1.0_dp
    forall(i=1:HDop,k=1:nzL) loc2D(i,k)=sum(u(i,:,k))
    call MPI_ALLreduce(loc2D,Ls,HDop*nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Ls = 2.0_dp*Ls/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC==1 .OR. xBC(1)==-1) then ! symmetry
    Lq = dble(xBC(1))
    Ls = 0.0_dp
  endif
!.RHS BC of u
  if (xBC(2)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Rq = -1.0_dp
    forall(i=1:HDop,k=1:nzL) loc2D(i,k)=sum(u(nxf-i+1,:,k))
    call MPI_ALLreduce(loc2D,rs,HDop*nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Rs = 2.0_dp*Rs/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC(2)==1 .OR. xBC(2)==-1) then ! symmetry
    Rq = dble(xBC(2))
    Rs = 0.0_dp
  endif
!........................................................!

  allocate(sbufN(nxf,HDop,nzL),sbufE(HDop,nyf,nzL),sbufW(HDop,nyf,nzL),sbufS(nxf,HDop,nzL))
  allocate(rbufN(nxf,HDop,nzL),rbufE(HDop,nyf,nzL),rbufW(HDop,nyf,nzL),rbufS(nxf,HDop,nzL))
  allocate(sbufNE(ndiag,nzL),sbufSE(ndiag,nzL),sbufSW(ndiag,nzL),sbufNW(ndiag,nzL))
  allocate(rbufNE(ndiag,nzL),rbufSE(ndiag,nzL),rbufSW(ndiag,nzL),rbufNW(ndiag,nzL))
!.Send top 2 rows to North rank, receive rows from South rank
  sbufN=u(:,nyf-HDop+1:nyf,:)
  call MPI_ISEND(sbufN,nxf*HDop*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(1),ierr)
!.Send right 2 columns to East rank, receive columns from West rank
  sbufE=u(nxf-HDop+1:nxf,:,:)
  call MPI_ISEND(sbufE,HDop*nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(2),ierr)
!.Send left 2 columns to West rank, receive columns from East rank
  sbufW=u(1:HDop,:,:)
  call MPI_ISEND(sbufW,HDop*nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(3),ierr)
!.Send bottom 2 rows to South rank, receive rows from North rank
  sbufS=u(:,1:HDop,:)
  call MPI_ISEND(sbufS,nxf*HDop*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(4),ierr)
#IF (HDop==2)
  sbufNE(1,:)=u(nxf,nyf,:)
  sbufSE(1,:)=u(nxf,1,:)
  sbufSW(1,:)=u(1,1,:)
  sbufNW(1,:)=u(1,nyf,:)
#ELSE
!.First element is the directly diagonal also needed for HDop==2
!.The other two follow the clockwise order, starting at 12 o'clock
  forall(k=1:nzL) sbufNE(:,k)=(/u(nxf,nyf,k),u(nxf,nyf-1,k),u(nxf-1,nyf,k)/)
  forall(k=1:nzL) sbufSE(:,k)=(/u(nxf,1,k),u(nxf,2,k),u(nxf-1,1,k)/)
  forall(k=1:nzL) sbufSW(:,k)=(/u(1,1,k),u(1,2,k),u(2,1,k)/)
  forall(k=1:nzL) sbufNW(:,k)=(/u(1,nyf,k),u(2,nyf,k),u(1,nyf-1,k)/)
#ENDIF
!.Send NE corner to NE rank, receive from SW rank
  call MPI_ISEND(sbufNE,ndiag*nzL,MPI_DOUBLE_PRECISION,NErank, &
                 0,COMMhd(qn)%a(ng,3),req(9),ierr)
!.Send SE corner to SE rank, receive from NW rank
  call MPI_ISEND(sbufSE,ndiag*nzL,MPI_DOUBLE_PRECISION,SErank, &
                 0,COMMhd(qn)%a(ng,3),req(10),ierr)
!.Send SW corner to SW rank, receive from NE rank
  call MPI_ISEND(sbufSW,ndiag*nzL,MPI_DOUBLE_PRECISION,SWrank, &
                 0,COMMhd(qn)%a(ng,3),req(11),ierr)
!.Send NW corner to NW rank, receive from SE rank
  call MPI_ISEND(sbufNW,ndiag*nzL,MPI_DOUBLE_PRECISION,NWrank, &
                 0,COMMhd(qn)%a(ng,3),req(12),ierr)
!.Receive rows from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxf*HDop*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(5),ierr)
!.Receive columns from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,HDop*nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(6),ierr)
!.Receive columns from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,HDop*nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(7),ierr)
!.Receive rows from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxf*HDop*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(8),ierr)
!.Receive NE corner from SW rank, for u(i-1,j-1)
  call MPI_IRECV(rbufSW,ndiag*nzL,MPI_DOUBLE_PRECISION,SWrank, &
                 0,COMMhd(qn)%a(ng,3),req(13),ierr)
!.Receive NW corner from SE rank, for u(i+1,j-1)
  call MPI_IRECV(rbufSE,ndiag*nzL,MPI_DOUBLE_PRECISION,SErank, &
                 0,COMMhd(qn)%a(ng,3),req(14),ierr)
!.Receive SE corner from NW rank, for u(i-1,j+1)
  call MPI_IRECV(rbufNW,ndiag*nzL,MPI_DOUBLE_PRECISION,NWrank, &
                 0,COMMhd(qn)%a(ng,3),req(15),ierr)
!.Receive SW corner from NE rank, for u(i+1,j+1)
  call MPI_IRECV(rbufNE,ndiag*nzL,MPI_DOUBLE_PRECISION,NErank, &
                 0,COMMhd(qn)%a(ng,3),req(16),ierr)
  
  call MPI_WAITALL(3,(/req(5),req(6),req(13)/),(/stat(:,5),stat(:,6),stat(:,13)/),ierr)

  rdx =  dble(nxf*iprocshd(qn)%a(ng))/bLx
  rdy =  dble(nyf*jprocshd(qn)%a(ng))/bLy
#IF (HDop==2)
  fxx =  hDiff(1)*(rdx**4)
  fyy =  hDiff(2)*(rdy**4)
  fxy =  (hDiff(1)+hDiff(2))*((rdx*rdy)**2)
  fx  = -4.0_dp*fxx-2.0_dp*fxy
  fy  = -4.0_dp*fyy-2.0_dp*fxy
  fc  =  1.0_dp+2.0_dp*(3.0_dp*(fxx+fyy)+2.0_dp*fxy)

  allocate(deltaL(1),deltaR(1))
  delta  = omehd(qn)/fc
  deltaL = omehd(qn)/(fc+fx*Lq)
  deltaR = omehd(qn)/(fc+fx*Rq)

!........................... FIRST ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
!...global SW corner
    u(1,1,:)=Oomehd(qn)*u(1,1,:)+(rhs(1,1,:) &
!    u(1,1,:)=(rhs(1,1,:) &
      -fxx*(u(3,1,:)+Lq*u(2,1,:)+Ls(2,:))-fx*(u(2,1,:)+Ls(1,:)) &
      -fyy*(u(1,3,:)+rbufS(1,1,:))-fy*(u(1,2,:)+rbufS(1,2,:)) &
      -fxy*(u(2,2,:)+rbufS(2,2,:)+Lq*rbufS(1,2,:)+Ls(1,:)+Lq*u(1,2,:)+Ls(1,:)))*deltaL
    u(2,1,:)=Oomehd(qn)*u(2,1,:)+(rhs(2,1,:) &
!    u(2,1,:)=(rhs(2,1,:) &
      -fxx*(u(4,1,:)+Lq*u(1,1,:)+Ls(1,:))-fx*(u(3,1,:)+u(1,1,:)) &
      -fyy*(u(2,3,:)+rbufS(2,1,:))-fy*(u(2,2,:)+rbufS(2,2,:)) &
      -fxy*(u(3,2,:)+rbufS(3,2,:)+rbufS(1,2,:)+u(1,2,:)))*delta
  else
!...local SW corner
    u(1,1,:)=Oomehd(qn)*u(1,1,:)+(rhs(1,1,:) &
!    u(1,1,:)=(rhs(1,1,:) &
      -fxx*(u(3,1,:)+rbufW(1,1,:))-fx*(u(2,1,:)+rbufW(2,1,:)) &
      -fyy*(u(1,3,:)+rbufS(1,1,:))-fy*(u(1,2,:)+rbufS(1,2,:)) &
      -fxy*(u(2,2,:)+rbufS(2,2,:)+rbufSW(1,:)+rbufW(2,2,:)))*delta
    u(2,1,:)=Oomehd(qn)*u(2,1,:)+(rhs(2,1,:) &
!    u(2,1,:)=(rhs(2,1,:) &
      -fxx*(u(4,1,:)+rbufW(2,1,:))-fx*(u(3,1,:)+u(1,1,:)) &
      -fyy*(u(2,3,:)+rbufS(2,1,:))-fy*(u(2,2,:)+rbufS(2,2,:)) &
      -fxy*(u(3,2,:)+rbufS(3,2,:)+rbufS(1,2,:)+u(1,2,:)))*delta
  endif
  do i=3,nxf-2 !.away from x-boundaries
    u(i,1,:)=Oomehd(qn)*u(i,1,:)+(rhs(i,1,:) &
!    u(i,1,:)=(rhs(i,1,:) &
      -fxx*(u(i+2,1,:)+u(i-2,1,:))-fx*(u(i+1,1,:)+u(i-1,1,:)) &
      -fyy*(u(i,3,:)+rbufS(i,1,:))-fy*(u(i,2,:)+rbufS(i,2,:)) &
      -fxy*(u(i+1,2,:)+rbufS(i+1,2,:)+rbufS(i-1,2,:)+u(i-1,2,:)))*delta
  enddo
  call MPI_WAITALL(2,(/req(7),req(14)/),(/stat(:,7),stat(:,14)/),ierr)
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-1,1,:)=Oomehd(qn)*u(nxf-1,1,:)+(rhs(nxf-1,1,:) &
!    u(nxf-1,1,:)=(rhs(nxf-1,1,:) &
      -fxx*(Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-3,1,:))-fx*(u(nxf,1,:)+u(nxf-2,1,:)) &
      -fyy*(u(nxf-1,3,:)+rbufS(nxf-1,1,:))-fy*(u(nxf-1,2,:)+rbufS(nxf-1,2,:)) &
      -fxy*(u(nxf,2,:)+rbufS(nxf,2,:)+rbufS(nxf-2,2,:)+u(nxf-2,2,:)))*delta
!...global SE corner
    u(nxf,1,:)=Oomehd(qn)*u(nxf,1,:)+(rhs(nxf,1,:) &
!    u(nxf,1,:)=(rhs(nxf,1,:) &
      -fxx*(Rq*u(nxf-1,1,:)+Rs(2,:)+u(nxf-2,1,:))-fx*(rs(1,:)+u(nxf-1,1,:)) &
      -fyy*(u(nxf,3,:)+rbufS(nxf,1,:))-fy*(u(nxf,2,:)+rbufS(nxf,2,:)) &
      -fxy*(Rq*u(nxf,2,:)+Rs(1,:)+Rq*rbufS(nxf,2,:)+Rs(1,:)+rbufS(nxf-1,2,:)+u(nxf-1,2,:)))*deltaR
  else
    u(nxf-1,1,:)=Oomehd(qn)*u(nxf-1,1,:)+(rhs(nxf-1,1,:) &
!    u(nxf-1,1,:)=(rhs(nxf-1,1,:) &
      -fxx*(rbufE(1,1,:)+u(nxf-3,1,:))-fx*(u(nxf,1,:)+u(nxf-2,1,:)) &
      -fyy*(u(nxf-1,3,:)+rbufS(nxf-1,1,:))-fy*(u(nxf-1,2,:)+rbufS(nxf-1,2,:)) &
      -fxy*(u(nxf,2,:)+rbufS(nxf,2,:)+rbufS(nxf-2,2,:)+u(nxf-2,2,:)))*delta
!...local SE corner
    u(nxf,1,:)=Oomehd(qn)*u(nxf,1,:)+(rhs(nxf,1,:) &
!    u(nxf,1,:)=(rhs(nxf,1,:) &
      -fxx*(rbufE(2,1,:)+u(nxf-2,1,:))-fx*(rbufE(1,1,:)+u(nxf-1,1,:)) &
      -fyy*(u(nxf,3,:)+rbufS(nxf,1,:))-fy*(u(nxf,2,:)+rbufS(nxf,2,:)) &
      -fxy*(rbufE(1,2,:)+rbufSE(1,:)+rbufS(nxf-1,2,:)+u(nxf-1,2,:)))*delta
  endif
!........................... SECOND ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,2,:)=Oomehd(qn)*u(1,2,:)+(rhs(1,2,:) &
!    u(1,2,:)=(rhs(1,2,:) &
      -fxx*(u(3,2,:)+Lq*u(2,2,:)+Ls(2,:))-fx*(u(2,2,:)+Ls(1,:)) &
      -fyy*(u(1,4,:)+rbufS(1,2,:))-fy*(u(1,3,:)+u(1,1,:)) &
      -fxy*(u(2,3,:)+u(2,1,:)+Lq*u(1,1,:)+Ls(1,:)+Lq*u(1,3,:)+Ls(1,:)))*deltaL
    u(2,2,:)=Oomehd(qn)*u(2,2,:)+(rhs(2,2,:) &
!    u(2,2,:)=(rhs(2,2,:) &
      -fxx*(u(4,2,:)+Lq*u(1,2,:)+Ls(1,:))-fx*(u(3,2,:)+u(1,2,:)) &
      -fyy*(u(2,4,:)+rbufS(2,2,:))-fy*(u(2,3,:)+u(2,1,:)) &
      -fxy*(u(3,3,:)+u(3,1,:)+u(1,1,:)+u(1,3,:)))*delta
  else
    u(1,2,:)=Oomehd(qn)*u(1,2,:)+(rhs(1,2,:) &
!    u(1,2,:)=(rhs(1,2,:) &
      -fxx*(u(3,2,:)+rbufW(1,2,:))-fx*(u(2,2,:)+rbufW(2,2,:)) &
      -fyy*(u(1,4,:)+rbufS(1,2,:))-fy*(u(1,3,:)+u(1,1,:)) &
      -fxy*(u(2,3,:)+u(2,1,:)+rbufW(2,1,:)+rbufW(2,3,:)))*delta
    u(2,2,:)=Oomehd(qn)*u(2,2,:)+(rhs(2,2,:) &
!    u(2,2,:)=(rhs(2,2,:) &
      -fxx*(u(4,2,:)+rbufW(2,2,:))-fx*(u(3,2,:)+u(1,2,:)) &
      -fyy*(u(2,4,:)+rbufS(2,2,:))-fy*(u(2,3,:)+u(2,1,:)) &
      -fxy*(u(3,3,:)+u(3,1,:)+u(1,1,:)+u(1,3,:)))*delta
  endif
  do i=3,nxf-2 !.away from x-boundaries
    u(i,2,:)=Oomehd(qn)*u(i,2,:)+(rhs(i,2,:) &
!    u(i,2,:)=(rhs(i,2,:) &
      -fxx*(u(i+2,2,:)+u(i-2,2,:))-fx*(u(i+1,2,:)+u(i-1,2,:)) &
      -fyy*(u(i,4,:)+rbufS(i,2,:))-fy*(u(i,3,:)+u(i,1,:)) &
      -fxy*(u(i+1,3,:)+u(i+1,1,:)+u(i-1,1,:)+u(i-1,3,:)))*delta
  enddo
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-1,2,:)=Oomehd(qn)*u(nxf-1,2,:)+(rhs(nxf-1,2,:) &
!    u(nxf-1,2,:)=(rhs(nxf-1,2,:) &
      -fxx*(Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-3,2,:))-fx*(u(nxf,2,:)+u(nxf-2,2,:)) &
      -fyy*(u(nxf-1,4,:)+rbufS(nxf-1,2,:))-fy*(u(nxf-1,3,:)+u(nxf-1,1,:)) &
      -fxy*(u(nxf,3,:)+u(nxf,1,:)+u(nxf-2,1,:)+u(nxf-2,3,:)))*delta
    u(nxf,2,:)=Oomehd(qn)*u(nxf,2,:)+(rhs(nxf,2,:) &
!    u(nxf,2,:)=(rhs(nxf,2,:) &
      -fxx*(Rq*u(nxf-1,2,:)+Rs(2,:)+u(nxf-2,2,:))-fx*(rs(1,:)+u(nxf-1,2,:)) &
      -fyy*(u(nxf,4,:)+rbufS(nxf,2,:))-fy*(u(nxf,3,:)+u(nxf,1,:)) &
      -fxy*(Rq*u(nxf,3,:)+Rs(1,:)+Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-1,1,:)+u(nxf-1,3,:)))*deltaR
  else
    u(nxf-1,2,:)=Oomehd(qn)*u(nxf-1,2,:)+(rhs(nxf-1,2,:) &
!    u(nxf-1,2,:)=(rhs(nxf-1,2,:) &
      -fxx*(rbufE(1,2,:)+u(nxf-3,2,:))-fx*(u(nxf,2,:)+u(nxf-2,2,:)) &
      -fyy*(u(nxf-1,4,:)+rbufS(nxf-1,2,:))-fy*(u(nxf-1,3,:)+u(nxf-1,1,:)) &
      -fxy*(u(nxf,3,:)+u(nxf,1,:)+u(nxf-2,1,:)+u(nxf-2,3,:)))*delta
    u(nxf,2,:)=Oomehd(qn)*u(nxf,2,:)+(rhs(nxf,2,:) &
!    u(nxf,2,:)=(rhs(nxf,2,:) &
      -fxx*(rbufE(2,2,:)+u(nxf-2,2,:))-fx*(rbufE(1,2,:)+u(nxf-1,2,:)) &
      -fyy*(u(nxf,4,:)+rbufS(nxf,2,:))-fy*(u(nxf,3,:)+u(nxf,1,:)) &
      -fxy*(rbufE(1,3,:)+rbufE(1,1,:)+u(nxf-1,1,:)+u(nxf-1,3,:)))*delta
  endif
!................... AWAY FROM BOTTOM AND TOP BOUNDARIES ...................!
  do j=3,nyf-2
!...Left boundary
    if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
      u(1,j,:)=Oomehd(qn)*u(1,j,:)+(rhs(1,j,:) &
!      u(1,j,:)=(rhs(1,j,:) &
        -fxx*(u(3,j,:)+Lq*u(2,j,:)+Ls(2,:))-fx*(u(2,j,:)+Ls(1,:)) &
        -fyy*(u(1,j+2,:)+u(1,j-2,:))-fy*(u(1,j+1,:)+u(1,j-1,:)) &
        -fxy*(u(2,j+1,:)+u(2,j-1,:)+Lq*u(1,j-1,:)+Ls(1,:)+Lq*u(1,j+1,:)+Ls(1,:)))*deltaL
      u(2,j,:)=Oomehd(qn)*u(2,j,:)+(rhs(2,j,:) &
!      u(2,j,:)=(rhs(2,j,:) &
        -fxx*(u(4,j,:)+Lq*u(1,j,:)+Ls(1,:))-fx*(u(3,j,:)+u(1,j,:)) &
        -fyy*(u(2,j+2,:)+u(2,j-2,:))-fy*(u(2,j+1,:)+u(2,j-1,:)) &
        -fxy*(u(3,j+1,:)+u(3,j-1,:)+u(1,j-1,:)+u(1,j+1,:)))*delta
    else
      u(1,j,:)=Oomehd(qn)*u(1,j,:)+(rhs(1,j,:) &
!      u(1,j,:)=(rhs(1,j,:) &
        -fxx*(u(3,j,:)+rbufW(1,j,:))-fx*(u(2,j,:)+rbufW(2,j,:)) &
        -fyy*(u(1,j+2,:)+u(1,j-2,:))-fy*(u(1,j+1,:)+u(1,j-1,:)) &
        -fxy*(u(2,j+1,:)+u(2,j-1,:)+rbufW(2,j-1,:)+rbufW(2,j+1,:)))*delta
      u(2,j,:)=Oomehd(qn)*u(2,j,:)+(rhs(2,j,:) &
!      u(2,j,:)=(rhs(2,j,:) &
        -fxx*(u(4,j,:)+rbufW(2,j,:))-fx*(u(3,j,:)+u(1,j,:)) &
        -fyy*(u(2,j+2,:)+u(2,j-2,:))-fy*(u(2,j+1,:)+u(2,j-1,:)) &
        -fxy*(u(3,j+1,:)+u(3,j-1,:)+u(1,j-1,:)+u(1,j+1,:)))*delta
    endif
!...Inner points
    do i=3,nxf-2
      u(i,j,:)=Oomehd(qn)*u(i,j,:)+(rhs(i,j,:) &
!      u(i,j,:)=(rhs(i,j,:) &
        -fxx*(u(i+2,j,:)+u(i-2,j,:))-fx*(u(i+1,j,:)+u(i-1,j,:)) &
        -fyy*(u(i,j+2,:)+u(i,j-2,:))-fy*(u(i,j+1,:)+u(i,j-1,:)) &
        -fxy*(u(i+1,j+1,:)+u(i+1,j-1,:)+u(i-1,j-1,:)+u(i-1,j+1,:)))*delta
    enddo
!...Right boundary
    if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
      u(nxf-1,j,:)=Oomehd(qn)*u(nxf-1,j,:)+(rhs(nxf-1,j,:) &
!      u(nxf-1,j,:)=(rhs(nxf-1,j,:) &
        -fxx*(Rq*u(nxf,j,:)+Rs(1,:)+u(nxf-3,j,:))-fx*(u(nxf,j,:)+u(nxf-2,j,:)) &
        -fyy*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:))-fy*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:)) &
        -fxy*(u(nxf,j+1,:)+u(nxf,j-1,:)+u(nxf-2,j-1,:)+u(nxf-2,j+1,:)))*delta
      u(nxf,j,:)=Oomehd(qn)*u(nxf,j,:)+(rhs(nxf,j,:) &
!      u(nxf,j,:)=(rhs(nxf,j,:) &
        -fxx*(Rq*u(nxf-1,j,:)+Rs(2,:)+u(nxf-2,j,:))-fx*(rs(1,:)+u(nxf-1,j,:)) &
        -fyy*(u(nxf,j+2,:)+u(nxf,j-2,:))-fy*(u(nxf,j+1,:)+u(nxf,j-1,:)) &
        -fxy*(Rq*u(nxf,j+1,:)+Rs(1,:)+Rq*u(nxf,j-1,:)+Rs(1,:)+u(nxf-1,j-1,:)+u(nxf-1,j+1,:)))*deltaR
    else
      u(nxf-1,j,:)=Oomehd(qn)*u(nxf-1,j,:)+(rhs(nxf-1,j,:) &
!      u(nxf-1,j,:)=(rhs(nxf-1,j,:) &
        -fxx*(rbufE(1,j,:)+u(nxf-3,j,:))-fx*(u(nxf,j,:)+u(nxf-2,j,:)) &
        -fyy*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:))-fy*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:)) &
        -fxy*(u(nxf,j+1,:)+u(nxf,j-1,:)+u(nxf-2,j-1,:)+u(nxf-2,j+1,:)))*delta
      u(nxf,j,:)=Oomehd(qn)*u(nxf,j,:)+(rhs(nxf,j,:) &
!      u(nxf,j,:)=(rhs(nxf,j,:) &
        -fxx*(rbufE(2,j,:)+u(nxf-2,j,:))-fx*(rbufE(1,j,:)+u(nxf-1,j,:)) &
        -fyy*(u(nxf,j+2,:)+u(nxf,j-2,:))-fy*(u(nxf,j+1,:)+u(nxf,j-1,:)) &
        -fxy*(rbufE(1,j+1,:)+rbufE(1,j-1,:)+u(nxf-1,j-1,:)+u(nxf-1,j+1,:)))*delta
    endif
  enddo
  call MPI_WAIT(req(8),stat(:,8),ierr)
!........................... SECOND TO TOP ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,nyf-1,:)=Oomehd(qn)*u(1,nyf-1,:)+(rhs(1,nyf-1,:) &
!    u(1,nyf-1,:)=(rhs(1,nyf-1,:) &
      -fxx*(u(3,nyf-1,:)+Lq*u(2,nyf-1,:)+Ls(2,:))-fx*(u(2,nyf-1,:)+Ls(1,:)) &
      -fyy*(rbufN(1,1,:)+u(1,nyf-3,:))-fy*(u(1,nyf,:)+u(1,nyf-2,:)) &
      -fxy*(u(2,nyf,:)+u(2,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:)+Lq*u(1,nyf,:)+Ls(1,:)))*deltaL
    u(2,nyf-1,:)=Oomehd(qn)*u(2,nyf-1,:)+(rhs(2,nyf-1,:) &
!    u(2,nyf-1,:)=(rhs(2,nyf-1,:) &
      -fxx*(u(4,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:))-fx*(u(3,nyf-1,:)+u(1,nyf-1,:)) &
      -fyy*(rbufN(2,1,:)+u(2,nyf-3,:))-fy*(u(2,nyf,:)+u(2,nyf-2,:)) &
      -fxy*(u(3,nyf,:)+u(3,nyf-2,:)+u(1,nyf-2,:)+u(1,nyf,:)))*delta
  else
    u(1,nyf-1,:)=Oomehd(qn)*u(1,nyf-1,:)+(rhs(1,nyf-1,:) &
!    u(1,nyf-1,:)=(rhs(1,nyf-1,:) &
      -fxx*(u(3,nyf-1,:)+rbufW(1,nyf-1,:))-fx*(u(2,nyf-1,:)+rbufW(2,nyf-1,:)) &
      -fyy*(rbufN(1,1,:)+u(1,nyf-3,:))-fy*(u(1,nyf,:)+u(1,nyf-2,:)) &
      -fxy*(u(2,nyf,:)+u(2,nyf-2,:)+rbufW(2,nyf-2,:)+rbufW(2,nyf,:)))*delta
    u(2,nyf-1,:)=Oomehd(qn)*u(2,nyf-1,:)+(rhs(2,nyf-1,:) &
!    u(2,nyf-1,:)=(rhs(2,nyf-1,:) &
      -fxx*(u(4,nyf-1,:)+rbufW(2,nyf-1,:))-fx*(u(3,nyf-1,:)+u(1,nyf-1,:)) &
      -fyy*(rbufN(2,1,:)+u(2,nyf-3,:))-fy*(u(2,nyf,:)+u(2,nyf-2,:)) &
      -fxy*(u(3,nyf,:)+u(3,nyf-2,:)+u(1,nyf-2,:)+u(1,nyf,:)))*delta
  endif
  do i=3,nxf-2 !.away from x-boundaries
    u(i,nyf-1,:)=Oomehd(qn)*u(i,nyf-1,:)+(rhs(i,nyf-1,:) &
!    u(i,nyf-1,:)=(rhs(i,nyf-1,:) &
      -fxx*(u(i+2,nyf-1,:)+u(i-2,nyf-1,:))-fx*(u(i+1,nyf-1,:)+u(i-1,nyf-1,:)) &
      -fyy*(rbufN(i,1,:)+u(i,nyf-3,:))-fy*(u(i,nyf,:)+u(i,nyf-2,:)) &
      -fxy*(u(i+1,nyf,:)+u(i+1,nyf-2,:)+u(i-1,nyf-2,:)+u(i-1,nyf,:)))*delta
  enddo
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-1,nyf-1,:)=Oomehd(qn)*u(nxf-1,nyf-1,:)+(rhs(nxf-1,nyf-1,:) &
!    u(nxf-1,nyf-1,:)=(rhs(nxf-1,nyf-1,:) &
      -fxx*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-3,nyf-1,:))-fx*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)) &
      -fyy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:))-fy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:)) &
      -fxy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)+u(nxf-2,nyf,:)))*delta
    u(nxf,nyf-1,:)=Oomehd(qn)*u(nxf,nyf-1,:)+(rhs(nxf,nyf-1,:) &
!    u(nxf,nyf-1,:)=(rhs(nxf,nyf-1,:) &
      -fxx*(Rq*u(nxf-1,nyf-1,:)+Rs(2,:)+u(nxf-2,nyf-1,:))-fx*(rs(1,:)+u(nxf-1,nyf-1,:)) &
      -fyy*(rbufN(nxf,1,:)+u(nxf,nyf-3,:))-fy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)) &
      -fxy*(Rq*u(nxf,nyf,:)+Rs(1,:)+Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-1,nyf-2,:)+u(nxf-1,nyf,:)))*deltaR
  else
    u(nxf-1,nyf-1,:)=Oomehd(qn)*u(nxf-1,nyf-1,:)+(rhs(nxf-1,nyf-1,:) &
!    u(nxf-1,nyf-1,:)=(rhs(nxf-1,nyf-1,:) &
      -fxx*(rbufE(1,nyf-1,:)+u(nxf-3,nyf-1,:))-fx*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)) &
      -fyy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:))-fy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:)) &
      -fxy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)+u(nxf-2,nyf,:)))*delta
    u(nxf,nyf-1,:)=Oomehd(qn)*u(nxf,nyf-1,:)+(rhs(nxf,nyf-1,:) &
!    u(nxf,nyf-1,:)=(rhs(nxf,nyf-1,:) &
      -fxx*(rbufE(2,nyf-1,:)+u(nxf-2,nyf-1,:))-fx*(rbufE(1,nyf-1,:)+u(nxf-1,nyf-1,:)) &
      -fyy*(rbufN(nxf,1,:)+u(nxf,nyf-3,:))-fy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)) &
      -fxy*(rbufE(1,nyf,:)+rbufE(1,nyf-2,:)+u(nxf-1,nyf-2,:)+u(nxf-1,nyf,:)))*delta
  endif
!........................... TOP ROW ...........................!
  call MPI_WAIT(req(15),stat(:,15),ierr) !.NW corner grid point, u(i-1,j+1)
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
!...global NW corner
    u(1,nyf,:)=Oomehd(qn)*u(1,nyf,:)+(rhs(1,nyf,:) &
!    u(1,nyf,:)=(rhs(1,nyf,:) &
      -fxx*(u(3,nyf,:)+Lq*u(2,nyf,:)+Ls(2,:))-fx*(u(2,nyf,:)+Ls(1,:)) &
      -fyy*(rbufN(1,2,:)+u(1,nyf-2,:))-fy*(rbufN(1,1,:)+u(1,nyf-1,:)) &
      -fxy*(rbufN(2,1,:)+u(2,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:)+Lq*rbufN(1,1,:)+Ls(1,:)))*deltaL
    u(2,nyf,:)=Oomehd(qn)*u(2,nyf,:)+(rhs(2,nyf,:) &
!    u(2,nyf,:)=(rhs(2,nyf,:) &
      -fxx*(u(4,nyf,:)+Lq*u(1,nyf,:)+Ls(1,:))-fx*(u(3,nyf,:)+u(1,nyf,:)) &
      -fyy*(rbufN(2,2,:)+u(2,nyf-2,:))-fy*(rbufN(2,1,:)+u(2,nyf-1,:)) &
      -fxy*(rbufN(3,1,:)+u(3,nyf-1,:)+u(1,nyf-1,:)+rbufN(1,1,:)))*delta
  else
!...local NW corner
    u(1,nyf,:)=Oomehd(qn)*u(1,nyf,:)+(rhs(1,nyf,:) &
!    u(1,nyf,:)=(rhs(1,nyf,:) &
      -fxx*(u(3,nyf,:)+rbufW(1,nyf,:))-fx*(u(2,nyf,:)+rbufW(2,nyf,:)) &
      -fyy*(rbufN(1,2,:)+u(1,nyf-2,:))-fy*(rbufN(1,1,:)+u(1,nyf-1,:)) &
      -fxy*(rbufN(2,1,:)+u(2,nyf-1,:)+rbufW(2,nyf-1,:)+rbufNW(1,:)))*delta
    u(2,nyf,:)=Oomehd(qn)*u(2,nyf,:)+(rhs(2,nyf,:) &
!    u(2,nyf,:)=(rhs(2,nyf,:) &
      -fxx*(u(4,nyf,:)+rbufW(2,nyf,:))-fx*(u(3,nyf,:)+u(1,nyf,:)) &
      -fyy*(rbufN(2,2,:)+u(2,nyf-2,:))-fy*(rbufN(2,1,:)+u(2,nyf-1,:)) &
      -fxy*(rbufN(3,1,:)+u(3,nyf-1,:)+u(1,nyf-1,:)+rbufN(1,1,:)))*delta
  endif
  do i=3,nxf-2 !.away from x-boundaries
    u(i,nyf,:)=Oomehd(qn)*u(i,nyf,:)+(rhs(i,nyf,:) &
!    u(i,nyf,:)=(rhs(i,nyf,:) &
      -fxx*(u(i+2,nyf,:)+u(i-2,nyf,:))-fx*(u(i+1,nyf,:)+u(i-1,nyf,:)) &
      -fyy*(rbufN(i,2,:)+u(i,nyf-2,:))-fy*(rbufN(i,1,:)+u(i,nyf-1,:)) &
      -fxy*(rbufN(i+1,1,:)+u(i+1,nyf-1,:)+u(i-1,nyf-1,:)+rbufN(i-1,1,:)))*delta
  enddo
  call MPI_WAIT(req(16),stat(:,16),ierr) !.NE corner grid point, u(i+1,j+1)
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-1,nyf,:)=Oomehd(qn)*u(nxf-1,nyf,:)+(rhs(nxf-1,nyf,:) &
      -fxx*(Rq*u(nxf,nyf,:)+Rs(1,:)+u(nxf-3,nyf,:))-fx*(u(nxf,nyf,:)+u(nxf-2,nyf,:)) &
      -fyy*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:))-fy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:)) &
      -fxy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)+rbufN(nxf-2,1,:)))*delta
!...global NE corner
    u(nxf,nyf,:)=Oomehd(qn)*u(nxf,nyf,:)+(rhs(nxf,nyf,:) &
!    u(nxf,nyf,:)=(rhs(nxf,nyf,:) &
      -fxx*(Rq*u(nxf-1,nyf,:)+Rs(2,:)+u(nxf-2,nyf,:))-fx*(rs(1,:)+u(nxf-1,nyf,:)) &
      -fyy*(rbufN(nxf,2,:)+u(nxf,nyf-2,:))-fy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)) &
      -fxy*(Rq*rbufN(nxf,1,:)+Rs(1,:)+Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-1,nyf-1,:)+rbufN(nxf-1,1,:)))*deltaR
  else
    u(nxf-1,nyf,:)=Oomehd(qn)*u(nxf-1,nyf,:)+(rhs(nxf-1,nyf,:) &
!    u(nxf-1,nyf,:)=(rhs(nxf-1,nyf,:) &
      -fxx*(rbufE(1,nyf,:)+u(nxf-3,nyf,:))-fx*(u(nxf,nyf,:)+u(nxf-2,nyf,:)) &
      -fyy*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:))-fy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:)) &
      -fxy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)+rbufN(nxf-2,1,:)))*delta
!...local NE corner
    u(nxf,nyf,:)=Oomehd(qn)*u(nxf,nyf,:)+(rhs(nxf,nyf,:) &
!    u(nxf,nyf,:)=(rhs(nxf,nyf,:) &
      -fxx*(rbufE(2,nyf,:)+u(nxf-2,nyf,:))-fx*(rbufE(1,nyf,:)+u(nxf-1,nyf,:)) &
      -fyy*(rbufN(nxf,2,:)+u(nxf,nyf-2,:))-fy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)) &
      -fxy*(rbufNE(1,:)+rbufE(1,nyf-1,:)+u(nxf-1,nyf-1,:)+rbufN(nxf-1,1,:)))*delta
  endif

#ELSE
!IF (HDop==3)
  fx6  =  hDiff(1)*(rdx**6)
  fy6  =  hDiff(2)*(rdy**6)
  fx4  =  hDiff(1)*(rdx**4)*(rdy**2)
  fy4  =  hDiff(2)*(rdx**2)*(rdy**4)
  fc   =  1.0_dp+4.0_dp*(5.0_dp*(fx6+fy6)+9.0_dp*(fx4+fy4))
  fx   =  15.0_dp*fx6+18.0_dp*fy4+24.0_dp*fx4
  fxx  = -6.0_dp*(fx6+fx4)
  fy   =  15.0_dp*fy6+18.0_dp*fx4+24.0_dp*fy4
  fyy  = -6.0_dp*(fy6+fy4)
  fxy  = -12.0_dp*(fx4+fy4)
  fxxy =  3.0_dp*fx4
  fxyy =  3.0_dp*fy4

  allocate(deltaL(2),deltaR(2))
  delta  = omehd(qn)/fc
  deltaL = (/omehd(qn)/(fc-fx*Lq),omehd(qn)/(fc-fx6*Lq)/)
  deltaR = (/omehd(qn)/(fc-fx*Rq),omehd(qn)/(fc-fx6*Rq)/)

!........................... FIRST ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
!...global SW corner
    u(1,1,:)=Oomehd(qn)*u(1,1,:)+(rhs(1,1,:) &
!    u(1,1,:)=(rhs(1,1,:) &
      +fx6*(u(4,1,:)+Lq*u(3,1,:)+Ls(3,:))+fxx*(u(3,1,:)+Lq*u(2,1,:)+Ls(2,:))+fx*(u(2,1,:)+Ls(1,:)) &
      +fy6*(u(1,4,:)+rbufS(1,1,:))+fyy*(u(1,3,:)+rbufS(1,2,:))+fy*(u(1,2,:)+rbufS(1,3,:)) &
      +fxy*(u(2,2,:)+rbufS(2,3,:)+Lq*rbufS(1,3,:)+Ls(1,:)+Lq*u(1,2,:)+Ls(1,:)) &
      +fxxy*(u(3,2,:)+rbufS(3,3,:)+Lq*rbufS(2,3,:)+Ls(2,:)+Lq*u(2,2,:)+Ls(2,:)) &
      +fxyy*(u(2,3,:)+rbufS(2,2,:)+Lq*rbufS(1,2,:)+Ls(1,:)+Lq*u(1,3,:)+Ls(1,:)))*deltaL(1)
    u(2,1,:)=Oomehd(qn)*u(2,1,:)+(rhs(2,1,:) &
!    u(2,1,:)=(rhs(2,1,:) &
      +fx6*(u(5,1,:)+Ls(2,:))+fxx*(u(4,1,:)+Lq*u(1,1,:)+Ls(1,:))+fx*(u(3,1,:)+u(1,1,:)) &
      +fy6*(u(2,4,:)+rbufS(2,1,:))+fyy*(u(2,3,:)+rbufS(2,2,:))+fy*(u(2,2,:)+rbufS(2,3,:)) &
      +fxy*(u(3,2,:)+rbufS(3,3,:)+rbufS(1,3,:)+u(1,2,:)) &
      +fxxy*(u(4,2,:)+rbufS(4,3,:)+Lq*rbufS(1,3,:)+Ls(1,:)+Lq*u(1,2,:)+Ls(1,:)) &
      +fxyy*(u(3,3,:)+rbufS(3,2,:)+rbufS(1,2,:)+u(1,3,:)))*deltaL(2)
    u(3,1,:)=Oomehd(qn)*u(3,1,:)+(rhs(3,1,:) &
!    u(3,1,:)=(rhs(3,1,:) &
      +fx6*(u(6,1,:)+Lq*u(1,1,:)+Ls(1,:))+fxx*(u(5,1,:)+u(1,1,:))+fx*(u(4,1,:)+u(2,1,:)) &
      +fy6*(u(3,4,:)+rbufS(3,1,:))+fyy*(u(3,3,:)+rbufS(3,2,:))+fy*(u(3,2,:)+rbufS(3,3,:)) &
      +fxy*(u(4,2,:)+rbufS(4,3,:)+rbufS(2,3,:)+u(2,2,:)) &
      +fxxy*(u(5,2,:)+rbufS(5,3,:)+rbufS(1,3,:)+u(1,2,:)) &
      +fxyy*(u(4,3,:)+rbufS(4,2,:)+rbufS(2,2,:)+u(2,3,:)))*delta
  else
!...local SW corner
    u(1,1,:)=Oomehd(qn)*u(1,1,:)+(rhs(1,1,:) &
!    u(1,1,:)=(rhs(1,1,:) &
      +fx6*(u(4,1,:)+rbufW(1,1,:))+fxx*(u(3,1,:)+rbufW(2,1,:))+fx*(u(2,1,:)+rbufW(3,1,:)) &
      +fy6*(u(1,4,:)+rbufS(1,1,:))+fyy*(u(1,3,:)+rbufS(1,2,:))+fy*(u(1,2,:)+rbufS(1,3,:)) &
      +fxy*(u(2,2,:)+rbufS(2,3,:)+rbufSW(1,:)+rbufW(3,2,:)) &
      +fxxy*(u(3,2,:)+rbufS(3,3,:)+rbufSW(3,:)+rbufW(2,2,:)) &
      +fxyy*(u(2,3,:)+rbufS(2,2,:)+rbufSW(2,:)+rbufW(3,3,:)))*delta
    u(2,1,:)=Oomehd(qn)*u(2,1,:)+(rhs(2,1,:) &
!    u(2,1,:)=(rhs(2,1,:) &
      +fx6*(u(5,1,:)+rbufW(2,1,:))+fxx*(u(4,1,:)+rbufW(3,1,:))+fx*(u(3,1,:)+u(1,1,:)) &
      +fy6*(u(2,4,:)+rbufS(2,1,:))+fyy*(u(2,3,:)+rbufS(2,2,:))+fy*(u(2,2,:)+rbufS(2,3,:)) &
      +fxy*(u(3,2,:)+rbufS(3,3,:)+rbufS(1,3,:)+u(1,2,:)) &
      +fxxy*(u(4,2,:)+rbufS(4,3,:)+rbufSW(1,:)+rbufW(3,2,:)) &
      +fxyy*(u(3,3,:)+rbufS(3,2,:)+rbufS(1,2,:)+u(1,3,:)))*delta
    u(3,1,:)=Oomehd(qn)*u(3,1,:)+(rhs(3,1,:) &
!    u(3,1,:)=(rhs(3,1,:) &
      +fx6*(u(6,1,:)+rbufW(3,1,:))+fxx*(u(5,1,:)+u(1,1,:))+fx*(u(4,1,:)+u(2,1,:)) &
      +fy6*(u(3,4,:)+rbufS(3,1,:))+fyy*(u(3,3,:)+rbufS(3,2,:))+fy*(u(3,2,:)+rbufS(3,3,:)) &
      +fxy*(u(4,2,:)+rbufS(4,3,:)+rbufS(2,3,:)+u(2,2,:)) &
      +fxxy*(u(5,2,:)+rbufS(5,3,:)+rbufS(1,3,:)+u(1,2,:)) &
      +fxyy*(u(4,3,:)+rbufS(4,2,:)+rbufS(2,2,:)+u(2,3,:)))*delta
  endif
  do i=4,nxf-3 !.away from x-boundaries
    u(i,1,:)=Oomehd(qn)*u(i,1,:)+(rhs(i,1,:) &
!    u(i,1,:)=(rhs(i,1,:) &
      +fx6*(u(i+3,1,:)+u(i-3,1,:))+fxx*(u(i+2,1,:)+u(i-2,1,:))+fx*(u(i+1,1,:)+u(i-1,1,:)) &
      +fy6*(u(i,4,:)+rbufS(i,1,:))+fyy*(u(i,3,:)+rbufS(i,2,:))+fy*(u(i,2,:)+rbufS(i,3,:)) &
      +fxy*(u(i+1,2,:)+rbufS(i+1,3,:)+rbufS(i-1,3,:)+u(i-1,2,:)) &
      +fxxy*(u(i+2,2,:)+rbufS(i+2,3,:)+rbufS(i-2,3,:)+u(i-2,2,:)) &
      +fxyy*(u(i+1,3,:)+rbufS(i+1,2,:)+rbufS(i-1,2,:)+u(i-1,3,:)))*delta
  enddo
  call MPI_WAITALL(2,(/req(7),req(14)/),(/stat(:,7),stat(:,14)/),ierr)
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-2,1,:)=Oomehd(qn)*u(nxf-2,1,:)+(rhs(nxf-2,1,:) &
!    u(nxf-2,1,:)=(rhs(nxf-2,1,:) &
      +fx6*(Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-5,1,:))+fxx*(u(nxf,1,:)+u(nxf-4,1,:))+fx*(u(nxf-1,1,:)+u(nxf-3,1,:)) &
      +fy6*(u(nxf-2,4,:)+rbufS(nxf-2,1,:))+fyy*(u(nxf-2,3,:)+rbufS(nxf-2,2,:))+fy*(u(nxf-2,2,:)+rbufS(nxf-2,3,:)) &
      +fxy*(u(nxf-1,2,:)+rbufS(nxf-1,3,:)+rbufS(nxf-3,3,:)+u(nxf-3,2,:)) &
      +fxxy*(u(nxf,2,:)+rbufS(nxf,3,:)+rbufS(nxf-4,3,:)+u(nxf-4,2,:)) &
      +fxyy*(u(nxf-1,3,:)+rbufS(nxf-1,2,:)+rbufS(nxf-3,2,:)+u(nxf-3,3,:)))*delta
    u(nxf-1,1,:)=Oomehd(qn)*u(nxf-1,1,:)+(rhs(nxf-1,1,:) &
!    u(nxf-1,1,:)=(rhs(nxf-1,1,:) &
      +fx6*(rs(2,:)+u(nxf-4,1,:))+fxx*(Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-3,1,:))+fx*(u(nxf,1,:)+u(nxf-2,1,:)) &
      +fy6*(u(nxf-1,4,:)+rbufS(nxf-1,1,:))+fyy*(u(nxf-1,3,:)+rbufS(nxf-1,2,:))+fy*(u(nxf-1,2,:)+rbufS(nxf-1,3,:)) &
      +fxy*(u(nxf,2,:)+rbufS(nxf,3,:)+rbufS(nxf-2,3,:)+u(nxf-2,2,:)) &
      +fxxy*(Rq*u(nxf,2,:)+Rs(1,:)+Rq*rbufS(nxf,3,:)+Rs(1,:)+rbufS(nxf-3,3,:)+u(nxf-3,2,:)) &
      +fxyy*(u(nxf,3,:)+rbufS(nxf,2,:)+rbufS(nxf-2,2,:)+u(nxf-2,3,:)))*deltaR(2)
!...global SE corner
    u(nxf,1,:)=Oomehd(qn)*u(nxf,1,:)+(rhs(nxf,1,:) &
!    u(nxf,1,:)=(rhs(nxf,1,:) &
      +fx6*(Rq*u(nxf-2,1,:)+Rs(3,:)+u(nxf-3,1,:))+fxx*(Rq*u(nxf-1,1,:)+Rs(2,:)+u(nxf-2,1,:))+fx*(rs(1,:)+u(nxf-1,1,:)) &
      +fy6*(u(nxf,4,:)+rbufS(nxf,1,:))+fyy*(u(nxf,3,:)+rbufS(nxf,2,:))+fy*(u(nxf,2,:)+rbufS(nxf,3,:)) &
      +fxy*(Rq*u(nxf,2,:)+Rs(1,:)+Rq*rbufS(nxf,3,:)+Rs(1,:)+rbufS(nxf-1,3,:)+u(nxf-1,2,:)) &
      +fxxy*(Rq*u(nxf-1,2,:)+Rs(2,:)+Rq*rbufS(nxf-1,3,:)+Rs(2,:)+rbufS(nxf-2,3,:)+u(nxf-2,2,:)) &
      +fxyy*(Rq*u(nxf,3,:)+Rs(1,:)+Rq*rbufS(nxf,2,:)+Rs(1,:)+rbufS(nxf-1,2,:)+u(nxf-1,3,:)))*deltaR(1)
  else
    u(nxf-2,1,:)=Oomehd(qn)*u(nxf-2,1,:)+(rhs(nxf-2,1,:) &
!    u(nxf-2,1,:)=(rhs(nxf-2,1,:) &
      +fx6*(rbufE(1,1,:)+u(nxf-5,1,:))+fxx*(u(nxf,1,:)+u(nxf-4,1,:))+fx*(u(nxf-1,1,:)+u(nxf-3,1,:)) &
      +fy6*(u(nxf-2,4,:)+rbufS(nxf-2,1,:))+fyy*(u(nxf-2,3,:)+rbufS(nxf-2,2,:))+fy*(u(nxf-2,2,:)+rbufS(nxf-2,3,:)) &
      +fxy*(u(nxf-1,2,:)+rbufS(nxf-1,3,:)+rbufS(nxf-3,3,:)+u(nxf-3,2,:)) &
      +fxxy*(u(nxf,2,:)+rbufS(nxf,3,:)+rbufS(nxf-4,3,:)+u(nxf-4,2,:)) &
      +fxyy*(u(nxf-1,3,:)+rbufS(nxf-1,2,:)+rbufS(nxf-3,2,:)+u(nxf-3,3,:)))*delta
    u(nxf-1,1,:)=Oomehd(qn)*u(nxf-1,1,:)+(rhs(nxf-1,1,:) &
!    u(nxf-1,1,:)=(rhs(nxf-1,1,:) &
      +fx6*(rbufE(2,1,:)+u(nxf-4,1,:))+fxx*(rbufE(1,1,:)+u(nxf-3,1,:))+fx*(u(nxf,1,:)+u(nxf-2,1,:)) &
      +fy6*(u(nxf-1,4,:)+rbufS(nxf-1,1,:))+fyy*(u(nxf-1,3,:)+rbufS(nxf-1,2,:))+fy*(u(nxf-1,2,:)+rbufS(nxf-1,3,:)) &
      +fxy*(u(nxf,2,:)+rbufS(nxf,3,:)+rbufS(nxf-2,3,:)+u(nxf-2,2,:)) &
      +fxxy*(rbufE(1,2,:)+rbufSE(1,:)+rbufS(nxf-3,3,:)+u(nxf-3,2,:)) &
      +fxyy*(u(nxf,3,:)+rbufS(nxf,2,:)+rbufS(nxf-2,2,:)+u(nxf-2,3,:)))*delta
!...local SE corner
    u(nxf,1,:)=Oomehd(qn)*u(nxf,1,:)+(rhs(nxf,1,:) &
!    u(nxf,1,:)=(rhs(nxf,1,:) &
      +fx6*(rbufE(3,1,:)+u(nxf-3,1,:))+fxx*(rbufE(2,1,:)+u(nxf-2,1,:))+fx*(rbufE(1,1,:)+u(nxf-1,1,:)) &
      +fy6*(u(nxf,4,:)+rbufS(nxf,1,:))+fyy*(u(nxf,3,:)+rbufS(nxf,2,:))+fy*(u(nxf,2,:)+rbufS(nxf,3,:)) &
      +fxy*(rbufE(1,2,:)+rbufSE(1,:)+rbufS(nxf-1,3,:)+u(nxf-1,2,:)) &
      +fxxy*(rbufE(2,2,:)+rbufSE(2,:)+rbufS(nxf-2,3,:)+u(nxf-2,2,:)) &
      +fxyy*(rbufE(1,3,:)+rbufSE(3,:)+rbufS(nxf-1,2,:)+u(nxf-1,3,:)))*delta
  endif
!........................... SECOND ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,2,:)=Oomehd(qn)*u(1,2,:)+(rhs(1,2,:) &
!    u(1,2,:)=(rhs(1,2,:) &
      +fx6*(u(4,2,:)+Lq*u(3,2,:)+Ls(3,:))+fxx*(u(3,2,:)+Lq*u(2,2,:)+Ls(2,:))+fx*(u(2,2,:)+Ls(1,:)) &
      +fy6*(u(1,5,:)+rbufS(1,2,:))+fyy*(u(1,4,:)+rbufS(1,3,:))+fy*(u(1,3,:)+u(1,1,:)) &
      +fxy*(u(2,3,:)+u(2,1,:)+Lq*u(1,1,:)+Ls(1,:)+Lq*u(1,3,:)+Ls(1,:)) &
      +fxxy*(u(3,3,:)+u(3,1,:)+Lq*u(2,1,:)+Ls(2,:)+Lq*u(2,3,:)+Ls(2,:)) &
      +fxyy*(u(2,4,:)+rbufS(2,3,:)+Lq*rbufS(1,3,:)+Ls(1,:)+Lq*u(1,4,:)+Ls(1,:)))*deltaL(1)
    u(2,2,:)=Oomehd(qn)*u(2,2,:)+(rhs(2,2,:) &
!    u(2,2,:)=(rhs(2,2,:) &
      +fx6*(u(5,2,:)+Ls(2,:))+fxx*(u(4,2,:)+Lq*u(1,2,:)+Ls(1,:))+fx*(u(3,2,:)+u(1,2,:)) &
      +fy6*(u(2,5,:)+rbufS(2,2,:))+fyy*(u(2,4,:)+rbufS(2,3,:))+fy*(u(2,3,:)+u(2,1,:)) &
      +fxy*(u(3,3,:)+u(3,1,:)+u(1,1,:)+u(1,3,:)) &
      +fxxy*(u(4,3,:)+u(4,1,:)+Lq*u(1,1,:)+Ls(1,:)+Lq*u(1,3,:)+Ls(1,:)) &
      +fxyy*(u(3,4,:)+rbufS(3,3,:)+rbufS(1,3,:)+u(1,4,:)))*deltaL(2)
    u(3,2,:)=Oomehd(qn)*u(3,2,:)+(rhs(3,2,:) &
!    u(3,2,:)=(rhs(3,2,:) &
      +fx6*(u(6,2,:)+Lq*u(1,2,:)+Ls(1,:))+fxx*(u(5,2,:)+u(1,2,:))+fx*(u(4,2,:)+u(2,2,:)) &
      +fy6*(u(3,5,:)+rbufS(3,2,:))+fyy*(u(3,4,:)+rbufS(3,3,:))+fy*(u(3,3,:)+u(3,1,:)) &
      +fxy*(u(4,3,:)+u(4,1,:)+u(2,1,:)+u(2,3,:)) &
      +fxxy*(u(5,3,:)+u(5,1,:)+u(1,1,:)+u(1,3,:)) &
      +fxyy*(u(4,4,:)+rbufS(4,3,:)+rbufS(2,3,:)+u(2,4,:)))*delta
  else
    u(1,2,:)=Oomehd(qn)*u(1,2,:)+(rhs(1,2,:) &
!    u(1,2,:)=(rhs(1,2,:) &
      +fx6*(u(4,2,:)+rbufW(1,2,:))+fxx*(u(3,2,:)+rbufW(2,2,:))+fx*(u(2,2,:)+rbufW(3,2,:)) &
      +fy6*(u(1,5,:)+rbufS(1,2,:))+fyy*(u(1,4,:)+rbufS(1,3,:))+fy*(u(1,3,:)+u(1,1,:)) &
      +fxy*(u(2,3,:)+u(2,1,:)+rbufW(3,1,:)+rbufW(3,3,:)) &
      +fxxy*(u(3,3,:)+u(3,1,:)+rbufW(2,1,:)+rbufW(2,3,:)) &
      +fxyy*(u(2,4,:)+rbufS(2,3,:)+rbufSW(1,:)+rbufW(3,4,:)))*delta
    u(2,2,:)=Oomehd(qn)*u(2,2,:)+(rhs(2,2,:) &
!    u(2,2,:)=(rhs(2,2,:) &
      +fx6*(u(5,2,:)+rbufW(2,2,:))+fxx*(u(4,2,:)+rbufW(3,2,:))+fx*(u(3,2,:)+u(1,2,:)) &
      +fy6*(u(2,5,:)+rbufS(2,2,:))+fyy*(u(2,4,:)+rbufS(2,3,:))+fy*(u(2,3,:)+u(2,1,:)) &
      +fxy*(u(3,3,:)+u(3,1,:)+u(1,1,:)+u(1,3,:)) &
      +fxxy*(u(4,3,:)+u(4,1,:)+rbufW(3,1,:)+rbufW(3,3,:)) &
      +fxyy*(u(3,4,:)+rbufS(3,3,:)+rbufS(1,3,:)+u(1,4,:)))*delta
    u(3,2,:)=Oomehd(qn)*u(3,2,:)+(rhs(3,2,:) &
!    u(3,2,:)=(rhs(3,2,:) &
      +fx6*(u(6,2,:)+rbufW(3,2,:))+fxx*(u(5,2,:)+u(1,2,:))+fx*(u(4,2,:)+u(2,2,:)) &
      +fy6*(u(3,5,:)+rbufS(3,2,:))+fyy*(u(3,4,:)+rbufS(3,3,:))+fy*(u(3,3,:)+u(3,1,:)) &
      +fxy*(u(4,3,:)+u(4,1,:)+u(2,1,:)+u(2,3,:)) &
      +fxxy*(u(5,3,:)+u(5,1,:)+u(1,1,:)+u(1,3,:)) &
      +fxyy*(u(4,4,:)+rbufS(4,3,:)+rbufS(2,3,:)+u(2,4,:)))*delta
  endif
  do i=4,nxf-3 !.away from x-boundaries
    u(i,2,:)=Oomehd(qn)*u(i,2,:)+(rhs(i,2,:) &
!    u(i,2,:)=(rhs(i,2,:) &
      +fx6*(u(i+3,2,:)+u(i-3,2,:))+fxx*(u(i+2,2,:)+u(i-2,2,:))+fx*(u(i+1,2,:)+u(i-1,2,:)) &
      +fy6*(u(i,5,:)+rbufS(i,2,:))+fyy*(u(i,4,:)+rbufS(i,3,:))+fy*(u(i,3,:)+u(i,1,:)) &
      +fxy*(u(i+1,3,:)+u(i+1,1,:)+u(i-1,1,:)+u(i-1,3,:)) &
      +fxxy*(u(i+2,3,:)+u(i+2,1,:)+u(i-2,1,:)+u(i-2,3,:)) &
      +fxyy*(u(i+1,4,:)+rbufS(i+1,3,:)+rbufS(i-1,3,:)+u(i-1,4,:)))*delta
  enddo
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-2,2,:)=Oomehd(qn)*u(nxf-2,2,:)+(rhs(nxf-2,2,:) &
!    u(nxf-2,2,:)=(rhs(nxf-2,2,:) &
      +fx6*(Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-5,2,:))+fxx*(u(nxf,2,:)+u(nxf-4,2,:))+fx*(u(nxf-1,2,:)+u(nxf-3,2,:)) &
      +fy6*(u(nxf-2,5,:)+rbufS(nxf-2,2,:))+fyy*(u(nxf-2,4,:)+rbufS(nxf-2,3,:))+fy*(u(nxf-2,3,:)+u(nxf-2,1,:)) &
      +fxy*(u(nxf-1,3,:)+u(nxf-1,1,:)+u(nxf-3,1,:)+u(nxf-3,3,:)) &
      +fxxy*(u(nxf,3,:)+u(nxf,1,:)+u(nxf-4,1,:)+u(nxf-4,3,:)) &
      +fxyy*(u(nxf-1,4,:)+rbufS(nxf-1,3,:)+rbufS(nxf-3,3,:)+u(nxf-3,4,:)))*delta
    u(nxf-1,2,:)=Oomehd(qn)*u(nxf-1,2,:)+(rhs(nxf-1,2,:) &
!    u(nxf-1,2,:)=(rhs(nxf-1,2,:) &
      +fx6*(rs(2,:)+u(nxf-4,2,:))+fxx*(Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-3,2,:))+fx*(u(nxf,2,:)+u(nxf-2,2,:)) &
      +fy6*(u(nxf-1,5,:)+rbufS(nxf-1,2,:))+fyy*(u(nxf-1,4,:)+rbufS(nxf-1,3,:))+fy*(u(nxf-1,3,:)+u(nxf-1,1,:)) &
      +fxy*(u(nxf,3,:)+u(nxf,1,:)+u(nxf-2,1,:)+u(nxf-2,3,:)) &
      +fxxy*(Rq*u(nxf,3,:)+Rs(1,:)+Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-3,1,:)+u(nxf-3,3,:)) &
      +fxyy*(u(nxf,4,:)+rbufS(nxf,3,:)+rbufS(nxf-2,3,:)+u(nxf-2,4,:)))*deltaR(2)
    u(nxf,2,:)=Oomehd(qn)*u(nxf,2,:)+(rhs(nxf,2,:) &
!    u(nxf,2,:)=(rhs(nxf,2,:) &
      +fx6*(Rq*u(nxf-2,2,:)+Rs(3,:)+u(nxf-3,2,:))+fxx*(Rq*u(nxf-1,2,:)+Rs(2,:)+u(nxf-2,2,:))+fx*(rs(1,:)+u(nxf-1,2,:)) &
      +fy6*(u(nxf,5,:)+rbufS(nxf,2,:))+fyy*(u(nxf,4,:)+rbufS(nxf,3,:))+fy*(u(nxf,3,:)+u(nxf,1,:)) &
      +fxy*(Rq*u(nxf,3,:)+Rs(1,:)+Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-1,1,:)+u(nxf-1,3,:)) &
      +fxxy*(Rq*u(nxf-1,3,:)+Rs(2,:)+Rq*u(nxf-1,1,:)+Rs(2,:)+u(nxf-2,1,:)+u(nxf-2,3,:)) &
      +fxyy*(Rq*u(nxf,4,:)+Rs(1,:)+Rq*rbufS(nxf,3,:)+Rs(1,:)+rbufS(nxf-1,3,:)+u(nxf-1,4,:)))*deltaR(1)
  else
    u(nxf-2,2,:)=Oomehd(qn)*u(nxf-2,2,:)+(rhs(nxf-2,2,:) &
!    u(nxf-2,2,:)=(rhs(nxf-2,2,:) &
      +fx6*(rbufE(1,2,:)+u(nxf-5,2,:))+fxx*(u(nxf,2,:)+u(nxf-4,2,:))+fx*(u(nxf-1,2,:)+u(nxf-3,2,:)) &
      +fy6*(u(nxf-2,5,:)+rbufS(nxf-2,2,:))+fyy*(u(nxf-2,4,:)+rbufS(nxf-2,3,:))+fy*(u(nxf-2,3,:)+u(nxf-2,1,:)) &
      +fxy*(u(nxf-1,3,:)+u(nxf-1,1,:)+u(nxf-3,1,:)+u(nxf-3,3,:)) &
      +fxxy*(u(nxf,3,:)+u(nxf,1,:)+u(nxf-4,1,:)+u(nxf-4,3,:)) &
      +fxyy*(u(nxf-1,4,:)+rbufS(nxf-1,3,:)+rbufS(nxf-3,3,:)+u(nxf-3,4,:)))*delta
    u(nxf-1,2,:)=Oomehd(qn)*u(nxf-1,2,:)+(rhs(nxf-1,2,:) &
!    u(nxf-1,2,:)=(rhs(nxf-1,2,:) &
      +fx6*(rbufE(2,2,:)+u(nxf-4,2,:))+fxx*(rbufE(1,2,:)+u(nxf-3,2,:))+fx*(u(nxf,2,:)+u(nxf-2,2,:)) &
      +fy6*(u(nxf-1,5,:)+rbufS(nxf-1,2,:))+fyy*(u(nxf-1,4,:)+rbufS(nxf-1,3,:))+fy*(u(nxf-1,3,:)+u(nxf-1,1,:)) &
      +fxy*(u(nxf,3,:)+u(nxf,1,:)+u(nxf-2,1,:)+u(nxf-2,3,:)) &
      +fxxy*(rbufE(1,3,:)+rbufE(1,1,:)+u(nxf-3,1,:)+u(nxf-3,3,:)) &
      +fxyy*(u(nxf,4,:)+rbufS(nxf,3,:)+rbufS(nxf-2,3,:)+u(nxf-2,4,:)))*delta
    u(nxf,2,:)=Oomehd(qn)*u(nxf,2,:)+(rhs(nxf,2,:) &
!    u(nxf,2,:)=(rhs(nxf,2,:) &
      +fx6*(rbufE(3,2,:)+u(nxf-3,2,:))+fxx*(rbufE(2,2,:)+u(nxf-2,2,:))+fx*(rbufE(1,2,:)+u(nxf-1,2,:)) &
      +fy6*(u(nxf,5,:)+rbufS(nxf,2,:))+fyy*(u(nxf,4,:)+rbufS(nxf,3,:))+fy*(u(nxf,3,:)+u(nxf,1,:)) &
      +fxy*(rbufE(1,3,:)+rbufE(1,1,:)+u(nxf-1,1,:)+u(nxf-1,3,:)) &
      +fxxy*(rbufE(2,3,:)+rbufE(2,1,:)+u(nxf-2,1,:)+u(nxf-2,3,:)) &
      +fxyy*(rbufE(1,4,:)+rbufSE(1,:)+rbufS(nxf-1,3,:)+u(nxf-1,4,:)))*delta
  endif
!........................... THIRD ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,3,:)=Oomehd(qn)*u(1,3,:)+(rhs(1,3,:) &
!    u(1,3,:)=(rhs(1,3,:) &
      +fx6*(u(4,3,:)+Lq*u(3,3,:)+Ls(3,:))+fxx*(u(3,3,:)+Lq*u(2,3,:)+Ls(2,:))+fx*(u(2,3,:)+Ls(1,:)) &
      +fy6*(u(1,6,:)+rbufS(1,3,:))+fyy*(u(1,5,:)+u(1,1,:))+fy*(u(1,4,:)+u(1,2,:)) &
      +fxy*(u(2,4,:)+u(2,2,:)+Lq*u(1,2,:)+Ls(1,:)+Lq*u(1,4,:)+Ls(1,:)) &
      +fxxy*(u(3,4,:)+u(3,2,:)+Lq*u(2,2,:)+Ls(2,:)+Lq*u(2,4,:)+Ls(2,:)) &
      +fxyy*(u(2,5,:)+u(2,1,:)+Lq*u(1,1,:)+Ls(1,:)+Lq*u(1,5,:)+Ls(1,:)))*deltaL(1)
    u(2,3,:)=Oomehd(qn)*u(2,3,:)+(rhs(2,3,:) &
!    u(2,3,:)=(rhs(2,3,:) &
      +fx6*(u(5,3,:)+Ls(2,:))+fxx*(u(4,3,:)+Lq*u(1,3,:)+Ls(1,:))+fx*(u(3,3,:)+u(1,3,:)) &
      +fy6*(u(2,6,:)+rbufS(2,3,:))+fyy*(u(2,5,:)+u(2,1,:))+fy*(u(2,4,:)+u(2,2,:)) &
      +fxy*(u(3,4,:)+u(3,2,:)+u(1,2,:)+u(1,4,:)) &
      +fxxy*(u(4,4,:)+u(4,2,:)+Lq*u(1,2,:)+Ls(1,:)+Lq*u(1,4,:)+Ls(1,:)) &
      +fxyy*(u(3,5,:)+u(3,1,:)+u(1,1,:)+u(1,5,:)))*deltaL(2)
    u(3,3,:)=Oomehd(qn)*u(3,3,:)+(rhs(3,3,:) &
!    u(3,3,:)=(rhs(3,3,:) &
      +fx6*(u(6,3,:)+Lq*u(1,3,:)+Ls(1,:))+fxx*(u(5,3,:)+u(1,3,:))+fx*(u(4,3,:)+u(2,3,:)) &
      +fy6*(u(3,6,:)+rbufS(3,3,:))+fyy*(u(3,5,:)+u(3,1,:))+fy*(u(3,4,:)+u(3,2,:)) &
      +fxy*(u(4,4,:)+u(4,2,:)+u(2,2,:)+u(2,4,:)) &
      +fxxy*(u(5,4,:)+u(5,2,:)+u(1,2,:)+u(1,4,:)) &
      +fxyy*(u(4,5,:)+u(4,1,:)+u(2,1,:)+u(2,5,:)))*delta
  else
    u(1,3,:)=Oomehd(qn)*u(1,3,:)+(rhs(1,3,:) &
!    u(1,3,:)=(rhs(1,3,:) &
      +fx6*(u(4,3,:)+rbufW(1,3,:))+fxx*(u(3,3,:)+rbufW(2,3,:))+fx*(u(2,3,:)+rbufW(3,3,:)) &
      +fy6*(u(1,6,:)+rbufS(1,3,:))+fyy*(u(1,5,:)+u(1,1,:))+fy*(u(1,4,:)+u(1,2,:)) &
      +fxy*(u(2,4,:)+u(2,2,:)+rbufW(3,2,:)+rbufW(3,4,:)) &
      +fxxy*(u(3,4,:)+u(3,2,:)+rbufW(2,2,:)+rbufW(2,4,:)) &
      +fxyy*(u(2,5,:)+u(2,1,:)+rbufW(3,1,:)+rbufW(3,5,:)))*delta
    u(2,3,:)=Oomehd(qn)*u(2,3,:)+(rhs(2,3,:) &
!    u(2,3,:)=(rhs(2,3,:) &
      +fx6*(u(5,3,:)+rbufW(2,3,:))+fxx*(u(4,3,:)+rbufW(3,3,:))+fx*(u(3,3,:)+u(1,3,:)) &
      +fy6*(u(2,6,:)+rbufS(2,3,:))+fyy*(u(2,5,:)+u(2,1,:))+fy*(u(2,4,:)+u(2,2,:)) &
      +fxy*(u(3,4,:)+u(3,2,:)+u(1,2,:)+u(1,4,:)) &
      +fxxy*(u(4,4,:)+u(4,2,:)+rbufW(3,2,:)+rbufW(3,4,:)) &
      +fxyy*(u(3,5,:)+u(3,1,:)+u(1,1,:)+u(1,5,:)))*delta
    u(3,3,:)=Oomehd(qn)*u(3,3,:)+(rhs(3,3,:) &
!    u(3,3,:)=(rhs(3,3,:) &
      +fx6*(u(6,3,:)+rbufW(3,3,:))+fxx*(u(5,3,:)+u(1,3,:))+fx*(u(4,3,:)+u(2,3,:)) &
      +fy6*(u(3,6,:)+rbufS(3,3,:))+fyy*(u(3,5,:)+u(3,1,:))+fy*(u(3,4,:)+u(3,2,:)) &
      +fxy*(u(4,4,:)+u(4,2,:)+u(2,2,:)+u(2,4,:)) &
      +fxxy*(u(5,4,:)+u(5,2,:)+u(1,2,:)+u(1,4,:)) &
      +fxyy*(u(4,5,:)+u(4,1,:)+u(2,1,:)+u(2,5,:)))*delta
  endif
  do i=4,nxf-3 !.away from x-boundaries
    u(i,3,:)=Oomehd(qn)*u(i,3,:)+(rhs(i,3,:) &
!    u(i,3,:)=(rhs(i,3,:) &
      +fx6*(u(i+3,3,:)+u(i-3,3,:))+fxx*(u(i+2,3,:)+u(i-2,3,:))+fx*(u(i+1,3,:)+u(i-1,3,:)) &
      +fy6*(u(i,6,:)+rbufS(i,3,:))+fyy*(u(i,5,:)+u(i,1,:))+fy*(u(i,4,:)+u(i,2,:)) &
      +fxy*(u(i+1,4,:)+u(i+1,2,:)+u(i-1,2,:)+u(i-1,4,:)) &
      +fxxy*(u(i+2,4,:)+u(i+2,2,:)+u(i-2,2,:)+u(i-2,4,:)) &
      +fxyy*(u(i+1,5,:)+u(i+1,1,:)+u(i-1,1,:)+u(i-1,5,:)))*delta
  enddo
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-2,3,:)=Oomehd(qn)*u(nxf-2,3,:)+(rhs(nxf-2,3,:) &
!    u(nxf-2,3,:)=(rhs(nxf-2,3,:) &
      +fx6*(Rq*u(nxf,3,:)+Rs(1,:)+u(nxf-5,3,:))+fxx*(u(nxf,3,:)+u(nxf-4,3,:))+fx*(u(nxf-1,3,:)+u(nxf-3,3,:)) &
      +fy6*(u(nxf-2,6,:)+rbufS(nxf-2,3,:))+fyy*(u(nxf-2,5,:)+u(nxf-2,1,:))+fy*(u(nxf-2,4,:)+u(nxf-2,2,:)) &
      +fxy*(u(nxf-1,4,:)+u(nxf-1,2,:)+u(nxf-3,2,:)+u(nxf-3,4,:)) &
      +fxxy*(u(nxf,4,:)+u(nxf,2,:)+u(nxf-4,2,:)+u(nxf-4,4,:)) &
      +fxyy*(u(nxf-1,5,:)+u(nxf-1,1,:)+u(nxf-3,1,:)+u(nxf-3,5,:)))*delta
    u(nxf-1,3,:)=Oomehd(qn)*u(nxf-1,3,:)+(rhs(nxf-1,3,:) &
!    u(nxf-1,3,:)=(rhs(nxf-1,3,:) &
      +fx6*(rs(2,:)+u(nxf-4,3,:))+fxx*(Rq*u(nxf,3,:)+Rs(1,:)+u(nxf-3,3,:))+fx*(u(nxf,3,:)+u(nxf-2,3,:)) &
      +fy6*(u(nxf-1,6,:)+rbufS(nxf-1,3,:))+fyy*(u(nxf-1,5,:)+u(nxf-1,1,:))+fy*(u(nxf-1,4,:)+u(nxf-1,2,:)) &
      +fxy*(u(nxf,4,:)+u(nxf,2,:)+u(nxf-2,2,:)+u(nxf-2,4,:)) &
      +fxxy*(Rq*u(nxf,4,:)+Rs(1,:)+Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-3,2,:)+u(nxf-3,4,:)) &
      +fxyy*(u(nxf,5,:)+u(nxf,1,:)+u(nxf-2,1,:)+u(nxf-2,5,:)))*deltaR(2)
    u(nxf,3,:)=Oomehd(qn)*u(nxf,3,:)+(rhs(nxf,3,:) &
!    u(nxf,3,:)=(rhs(nxf,3,:) &
      +fx6*(Rq*u(nxf-2,3,:)+Rs(3,:)+u(nxf-3,3,:))+fxx*(Rq*u(nxf-1,3,:)+Rs(2,:)+u(nxf-2,3,:))+fx*(rs(1,:)+u(nxf-1,3,:)) &
      +fy6*(u(nxf,6,:)+rbufS(nxf,3,:))+fyy*(u(nxf,5,:)+u(nxf,1,:))+fy*(u(nxf,4,:)+u(nxf,2,:)) &
      +fxy*(Rq*u(nxf,4,:)+Rs(1,:)+Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-1,2,:)+u(nxf-1,4,:)) &
      +fxxy*(Rq*u(nxf-1,4,:)+Rs(2,:)+Rq*u(nxf-1,2,:)+Rs(2,:)+u(nxf-2,2,:)+u(nxf-2,4,:)) &
      +fxyy*(Rq*u(nxf,5,:)+Rs(1,:)+Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-1,1,:)+u(nxf-1,5,:)))*deltaR(1)
  else
    u(nxf-2,3,:)=Oomehd(qn)*u(nxf-2,3,:)+(rhs(nxf-2,3,:) &
!    u(nxf-2,3,:)=(rhs(nxf-2,3,:) &
      +fx6*(rbufE(1,3,:)+u(nxf-5,3,:))+fxx*(u(nxf,3,:)+u(nxf-4,3,:))+fx*(u(nxf-1,3,:)+u(nxf-3,3,:)) &
      +fy6*(u(nxf-2,6,:)+rbufS(nxf-2,3,:))+fyy*(u(nxf-2,5,:)+u(nxf-2,1,:))+fy*(u(nxf-2,4,:)+u(nxf-2,2,:)) &
      +fxy*(u(nxf-1,4,:)+u(nxf-1,2,:)+u(nxf-3,2,:)+u(nxf-3,4,:)) &
      +fxxy*(u(nxf,4,:)+u(nxf,2,:)+u(nxf-4,2,:)+u(nxf-4,4,:)) &
      +fxyy*(u(nxf-1,5,:)+u(nxf-1,1,:)+u(nxf-3,1,:)+u(nxf-3,5,:)))*delta
    u(nxf-1,3,:)=Oomehd(qn)*u(nxf-1,3,:)+(rhs(nxf-1,3,:) &
!    u(nxf-1,3,:)=(rhs(nxf-1,3,:) &
      +fx6*(rbufE(2,3,:)+u(nxf-4,3,:))+fxx*(rbufE(1,3,:)+u(nxf-3,3,:))+fx*(u(nxf,3,:)+u(nxf-2,3,:)) &
      +fy6*(u(nxf-1,6,:)+rbufS(nxf-1,3,:))+fyy*(u(nxf-1,5,:)+u(nxf-1,1,:))+fy*(u(nxf-1,4,:)+u(nxf-1,2,:)) &
      +fxy*(u(nxf,4,:)+u(nxf,2,:)+u(nxf-2,2,:)+u(nxf-2,4,:)) &
      +fxxy*(rbufE(1,4,:)+rbufE(1,2,:)+u(nxf-3,2,:)+u(nxf-3,4,:)) &
      +fxyy*(u(nxf,5,:)+u(nxf,1,:)+u(nxf-2,1,:)+u(nxf-2,5,:)))*delta
    u(nxf,3,:)=Oomehd(qn)*u(nxf,3,:)+(rhs(nxf,3,:) &
!    u(nxf,3,:)=(rhs(nxf,3,:) &
      +fx6*(rbufE(3,3,:)+u(nxf-3,3,:))+fxx*(rbufE(2,3,:)+u(nxf-2,3,:))+fx*(rbufE(1,3,:)+u(nxf-1,3,:)) &
      +fy6*(u(nxf,6,:)+rbufS(nxf,3,:))+fyy*(u(nxf,5,:)+u(nxf,1,:))+fy*(u(nxf,4,:)+u(nxf,2,:)) &
      +fxy*(rbufE(1,4,:)+rbufE(1,2,:)+u(nxf-1,2,:)+u(nxf-1,4,:)) &
      +fxxy*(rbufE(2,4,:)+rbufE(2,2,:)+u(nxf-2,2,:)+u(nxf-2,4,:)) &
      +fxyy*(rbufE(1,5,:)+rbufE(1,1,:)+u(nxf-1,1,:)+u(nxf-1,5,:)))*delta
  endif
!................... AWAY FROM BOTTOM AND TOP BOUNDARIES ...................!
  do j=4,nyf-3
!...Left boundary
    if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
      u(1,j,:)=Oomehd(qn)*u(1,j,:)+(rhs(1,j,:) &
!      u(1,j,:)=(rhs(1,j,:) &
        +fx6*(u(4,j,:)+Lq*u(3,j,:)+Ls(3,:))+fxx*(u(3,j,:)+Lq*u(2,j,:)+Ls(2,:))+fx*(u(2,j,:)+Ls(1,:)) &
        +fy6*(u(1,j+3,:)+u(1,j-3,:))+fyy*(u(1,j+2,:)+u(1,j-2,:))+fy*(u(1,j+1,:)+u(1,j-1,:)) &
        +fxy*(u(2,j+1,:)+u(2,j-1,:)+Lq*u(1,j-1,:)+Ls(1,:)+Lq*u(1,j+1,:)+Ls(1,:)) &
        +fxxy*(u(3,j+1,:)+u(3,j-1,:)+Lq*u(2,j-1,:)+Ls(2,:)+Lq*u(2,j+1,:)+Ls(2,:)) &
        +fxyy*(u(2,j+2,:)+u(2,j-2,:)+Lq*u(1,j-2,:)+Ls(1,:)+Lq*u(1,j+2,:)+Ls(1,:)))*deltaL(1)
      u(2,j,:)=Oomehd(qn)*u(2,j,:)+(rhs(2,j,:) &
!      u(2,j,:)=(rhs(2,j,:) &
        +fx6*(u(5,j,:)+Ls(2,:))+fxx*(u(4,j,:)+Lq*u(1,j,:)+Ls(1,:))+fx*(u(3,j,:)+u(1,j,:)) &
        +fy6*(u(2,j+3,:)+u(2,j-3,:))+fyy*(u(2,j+2,:)+u(2,j-2,:))+fy*(u(2,j+1,:)+u(2,j-1,:)) &
        +fxy*(u(3,j+1,:)+u(3,j-1,:)+u(1,j-1,:)+u(1,j+1,:)) &
        +fxxy*(u(4,j+1,:)+u(4,j-1,:)+Lq*u(1,j-1,:)+Ls(1,:)+Lq*u(1,j+1,:)+Ls(1,:)) &
        +fxyy*(u(3,j+2,:)+u(3,j-2,:)+u(1,j-2,:)+u(1,j+2,:)))*deltaL(2)
      u(3,j,:)=Oomehd(qn)*u(3,j,:)+(rhs(3,j,:) &
!      u(3,j,:)=(rhs(3,j,:) &
        +fx6*(u(6,j,:)+Lq*u(1,j,:)+Ls(1,:))+fxx*(u(5,j,:)+u(1,j,:))+fx*(u(4,j,:)+u(2,j,:)) &
        +fy6*(u(3,j+3,:)+u(3,j-3,:))+fyy*(u(3,j+2,:)+u(3,j-2,:))+fy*(u(3,j+1,:)+u(3,j-1,:)) &
        +fxy*(u(4,j+1,:)+u(4,j-1,:)+u(2,j-1,:)+u(2,j+1,:)) &
        +fxxy*(u(5,j+1,:)+u(5,j-1,:)+u(1,j-1,:)+u(1,j+1,:)) &
        +fxyy*(u(4,j+2,:)+u(4,j-2,:)+u(2,j-2,:)+u(2,j+2,:)))*delta
    else
      u(1,j,:)=Oomehd(qn)*u(1,j,:)+(rhs(1,j,:) &
!      u(1,j,:)=(rhs(1,j,:) &
        +fx6*(u(4,j,:)+rbufW(1,j,:))+fxx*(u(3,j,:)+rbufW(2,j,:))+fx*(u(2,j,:)+rbufW(3,j,:)) &
        +fy6*(u(1,j+3,:)+u(1,j-3,:))+fyy*(u(1,j+2,:)+u(1,j-2,:))+fy*(u(1,j+1,:)+u(1,j-1,:)) &
        +fxy*(u(2,j+1,:)+u(2,j-1,:)+rbufW(3,j-1,:)+rbufW(3,j+1,:)) &
        +fxxy*(u(3,j+1,:)+u(3,j-1,:)+rbufW(2,j-1,:)+rbufW(2,j+1,:)) &
        +fxyy*(u(2,j+2,:)+u(2,j-2,:)+rbufW(3,j-2,:)+rbufW(3,j+2,:)))*delta
      u(2,j,:)=Oomehd(qn)*u(2,j,:)+(rhs(2,j,:) &
!      u(2,j,:)=(rhs(2,j,:) &
        +fx6*(u(5,j,:)+rbufW(2,j,:))+fxx*(u(4,j,:)+rbufW(3,j,:))+fx*(u(3,j,:)+u(1,j,:)) &
        +fy6*(u(2,j+3,:)+u(2,j-3,:))+fyy*(u(2,j+2,:)+u(2,j-2,:))+fy*(u(2,j+1,:)+u(2,j-1,:)) &
        +fxy*(u(3,j+1,:)+u(3,j-1,:)+u(1,j-1,:)+u(1,j+1,:)) &
        +fxxy*(u(4,j+1,:)+u(4,j-1,:)+rbufW(3,j-1,:)+rbufW(3,j+1,:)) &
        +fxyy*(u(3,j+2,:)+u(3,j-2,:)+u(1,j-2,:)+u(1,j+2,:)))*delta
      u(3,j,:)=Oomehd(qn)*u(3,j,:)+(rhs(3,j,:) &
!      u(3,j,:)=(rhs(3,j,:) &
        +fx6*(u(6,j,:)+rbufW(3,j,:))+fxx*(u(5,j,:)+u(1,j,:))+fx*(u(4,j,:)+u(2,j,:)) &
        +fy6*(u(3,j+3,:)+u(3,j-3,:))+fyy*(u(3,j+2,:)+u(3,j-2,:))+fy*(u(3,j+1,:)+u(3,j-1,:)) &
        +fxy*(u(4,j+1,:)+u(4,j-1,:)+u(2,j-1,:)+u(2,j+1,:)) &
        +fxxy*(u(5,j+1,:)+u(5,j-1,:)+u(1,j-1,:)+u(1,j+1,:)) &
        +fxyy*(u(4,j+2,:)+u(4,j-2,:)+u(2,j-2,:)+u(2,j+2,:)))*delta
    endif
!...Inner points
    do i=4,nxf-3
      u(i,j,:)=Oomehd(qn)*u(i,j,:)+(rhs(i,j,:) &
!      u(i,j,:)=(rhs(i,j,:) &
        +fx6*(u(i+3,j,:)+u(i-3,j,:))+fxx*(u(i+2,j,:)+u(i-2,j,:))+fx*(u(i+1,j,:)+u(i-1,j,:)) &
        +fy6*(u(i,j+3,:)+u(i,j-3,:))+fyy*(u(i,j+2,:)+u(i,j-2,:))+fy*(u(i,j+1,:)+u(i,j-1,:)) &
        +fxy*(u(i+1,j+1,:)+u(i+1,j-1,:)+u(i-1,j-1,:)+u(i-1,j+1,:)) &
        +fxxy*(u(i+2,j+1,:)+u(i+2,j-1,:)+u(i-2,j-1,:)+u(i-2,j+1,:)) &
        +fxyy*(u(i+1,j+2,:)+u(i+1,j-2,:)+u(i-1,j-2,:)+u(i-1,j+2,:)))*delta
    enddo
!...Right boundary
    if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
      u(nxf-2,j,:)=Oomehd(qn)*u(nxf-2,j,:)+(rhs(nxf-2,j,:) &
!      u(nxf-2,j,:)=(rhs(nxf-2,j,:) &
        +fx6*(Rq*u(nxf,j,:)+Rs(1,:)+u(nxf-5,j,:))+fxx*(u(nxf,j,:)+u(nxf-4,j,:))+fx*(u(nxf-1,j,:)+u(nxf-3,j,:)) &
        +fy6*(u(nxf-2,j+3,:)+u(nxf-2,j-3,:))+fyy*(u(nxf-2,j+2,:)+u(nxf-2,j-2,:))+fy*(u(nxf-2,j+1,:)+u(nxf-2,j-1,:)) &
        +fxy*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:)+u(nxf-3,j-1,:)+u(nxf-3,j+1,:)) &
        +fxxy*(u(nxf,j+1,:)+u(nxf,j-1,:)+u(nxf-4,j-1,:)+u(nxf-4,j+1,:)) &
        +fxyy*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:)+u(nxf-3,j-2,:)+u(nxf-3,j+2,:)))*delta
      u(nxf-1,j,:)=Oomehd(qn)*u(nxf-1,j,:)+(rhs(nxf-1,j,:) &
!      u(nxf-1,j,:)=(rhs(nxf-1,j,:) &
        +fx6*(rs(2,:)+u(nxf-4,j,:))+fxx*(Rq*u(nxf,j,:)+Rs(1,:)+u(nxf-3,j,:))+fx*(u(nxf,j,:)+u(nxf-2,j,:)) &
        +fy6*(u(nxf-1,j+3,:)+u(nxf-1,j-3,:))+fyy*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:))+fy*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:)) &
        +fxy*(u(nxf,j+1,:)+u(nxf,j-1,:)+u(nxf-2,j-1,:)+u(nxf-2,j+1,:)) &
        +fxxy*(Rq*u(nxf,j+1,:)+Rs(1,:)+Rq*u(nxf,j-1,:)+Rs(1,:)+u(nxf-3,j-1,:)+u(nxf-3,j+1,:)) &
        +fxyy*(u(nxf,j+2,:)+u(nxf,j-2,:)+u(nxf-2,j-2,:)+u(nxf-2,j+2,:)))*deltaR(2)
      u(nxf,j,:)=Oomehd(qn)*u(nxf,j,:)+(rhs(nxf,j,:) &
!      u(nxf,j,:)=(rhs(nxf,j,:) &
        +fx6*(Rq*u(nxf-2,j,:)+Rs(3,:)+u(nxf-3,j,:))+fxx*(Rq*u(nxf-1,j,:)+Rs(2,:)+u(nxf-2,j,:))+fx*(rs(1,:)+u(nxf-1,j,:)) &
        +fy6*(u(nxf,j+3,:)+u(nxf,j-3,:))+fyy*(u(nxf,j+2,:)+u(nxf,j-2,:))+fy*(u(nxf,j+1,:)+u(nxf,j-1,:)) &
        +fxy*(Rq*u(nxf,j+1,:)+Rs(1,:)+Rq*u(nxf,j-1,:)+Rs(1,:)+u(nxf-1,j-1,:)+u(nxf-1,j+1,:)) &
        +fxxy*(Rq*u(nxf-1,j+1,:)+Rs(2,:)+Rq*u(nxf-1,j-1,:)+Rs(2,:)+u(nxf-2,j-1,:)+u(nxf-2,j+1,:)) &
        +fxyy*(Rq*u(nxf,j+2,:)+Rs(1,:)+Rq*u(nxf,j-2,:)+Rs(1,:)+u(nxf-1,j-2,:)+u(nxf-1,j+2,:)))*deltaR(1)
    else
      u(nxf-2,j,:)=Oomehd(qn)*u(nxf-2,j,:)+(rhs(nxf-2,j,:) &
!      u(nxf-2,j,:)=(rhs(nxf-2,j,:) &
        +fx6*(rbufE(1,j,:)+u(nxf-5,j,:))+fxx*(u(nxf,j,:)+u(nxf-4,j,:))+fx*(u(nxf-1,j,:)+u(nxf-3,j,:)) &
        +fy6*(u(nxf-2,j+3,:)+u(nxf-2,j-3,:))+fyy*(u(nxf-2,j+2,:)+u(nxf-2,j-2,:))+fy*(u(nxf-2,j+1,:)+u(nxf-2,j-1,:)) &
        +fxy*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:)+u(nxf-3,j-1,:)+u(nxf-3,j+1,:)) &
        +fxxy*(u(nxf,j+1,:)+u(nxf,j-1,:)+u(nxf-4,j-1,:)+u(nxf-4,j+1,:)) &
        +fxyy*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:)+u(nxf-3,j-2,:)+u(nxf-3,j+2,:)))*delta
      u(nxf-1,j,:)=Oomehd(qn)*u(nxf-1,j,:)+(rhs(nxf-1,j,:) &
!      u(nxf-1,j,:)=(rhs(nxf-1,j,:) &
        +fx6*(rbufE(2,j,:)+u(nxf-4,j,:))+fxx*(rbufE(1,j,:)+u(nxf-3,j,:))+fx*(u(nxf,j,:)+u(nxf-2,j,:)) &
        +fy6*(u(nxf-1,j+3,:)+u(nxf-1,j-3,:))+fyy*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:))+fy*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:)) &
        +fxy*(u(nxf,j+1,:)+u(nxf,j-1,:)+u(nxf-2,j-1,:)+u(nxf-2,j+1,:)) &
        +fxxy*(rbufE(1,j+1,:)+rbufE(1,j-1,:)+u(nxf-3,j-1,:)+u(nxf-3,j+1,:)) &
        +fxyy*(u(nxf,j+2,:)+u(nxf,j-2,:)+u(nxf-2,j-2,:)+u(nxf-2,j+2,:)))*delta
      u(nxf,j,:)=Oomehd(qn)*u(nxf,j,:)+(rhs(nxf,j,:) &
!      u(nxf,j,:)=(rhs(nxf,j,:) &
        +fx6*(rbufE(3,j,:)+u(nxf-3,j,:))+fxx*(rbufE(2,j,:)+u(nxf-2,j,:))+fx*(rbufE(1,j,:)+u(nxf-1,j,:)) &
        +fy6*(u(nxf,j+3,:)+u(nxf,j-3,:))+fyy*(u(nxf,j+2,:)+u(nxf,j-2,:))+fy*(u(nxf,j+1,:)+u(nxf,j-1,:)) &
        +fxy*(rbufE(1,j+1,:)+rbufE(1,j-1,:)+u(nxf-1,j-1,:)+u(nxf-1,j+1,:)) &
        +fxxy*(rbufE(2,j+1,:)+rbufE(2,j-1,:)+u(nxf-2,j-1,:)+u(nxf-2,j+1,:)) &
        +fxyy*(rbufE(1,j+2,:)+rbufE(1,j-2,:)+u(nxf-1,j-2,:)+u(nxf-1,j+2,:)))*delta
    endif
  enddo
  call MPI_WAIT(req(8),stat(:,8),ierr)
!........................... THIRD TO LAST ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,nyf-2,:)=Oomehd(qn)*u(1,nyf-2,:)+(rhs(1,nyf-2,:) &
!    u(1,nyf-2,:)=(rhs(1,nyf-2,:) &
      +fx6*(u(4,nyf-2,:)+Lq*u(3,nyf-2,:)+Ls(3,:))+fxx*(u(3,nyf-2,:)+Lq*u(2,nyf-2,:)+Ls(2,:))+fx*(u(2,nyf-2,:)+Ls(1,:)) &
      +fy6*(rbufN(1,1,:)+u(1,nyf-5,:))+fyy*(u(1,nyf,:)+u(1,nyf-4,:))+fy*(u(1,nyf-1,:)+u(1,nyf-3,:)) &
      +fxy*(u(2,nyf-1,:)+u(2,nyf-3,:)+Lq*u(1,nyf-3,:)+Ls(1,:)+Lq*u(1,nyf-1,:)+Ls(1,:)) &
      +fxxy*(u(3,nyf-1,:)+u(3,nyf-3,:)+Lq*u(2,nyf-3,:)+Ls(2,:)+Lq*u(2,nyf-1,:)+Ls(2,:)) &
      +fxyy*(u(2,nyf,:)+u(2,nyf-4,:)+Lq*u(1,nyf-4,:)+Ls(1,:)+Lq*u(1,nyf,:)+Ls(1,:)))*deltaL(1)
    u(2,nyf-2,:)=Oomehd(qn)*u(2,nyf-2,:)+(rhs(2,nyf-2,:) &
!    u(2,nyf-2,:)=(rhs(2,nyf-2,:) &
      +fx6*(u(5,nyf-2,:)+Ls(2,:))+fxx*(u(4,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:))+fx*(u(3,nyf-2,:)+u(1,nyf-2,:)) &
      +fy6*(rbufN(2,1,:)+u(2,nyf-5,:))+fyy*(u(2,nyf,:)+u(2,nyf-4,:))+fy*(u(2,nyf-1,:)+u(2,nyf-3,:)) &
      +fxy*(u(3,nyf-1,:)+u(3,nyf-3,:)+u(1,nyf-3,:)+u(1,nyf-1,:)) &
      +fxxy*(u(4,nyf-1,:)+u(4,nyf-3,:)+Lq*u(1,nyf-3,:)+Ls(1,:)+Lq*u(1,nyf-1,:)+Ls(1,:)) &
      +fxyy*(u(3,nyf,:)+u(3,nyf-4,:)+u(1,nyf-4,:)+u(1,nyf,:)))*deltaL(2)
    u(3,nyf-2,:)=Oomehd(qn)*u(3,nyf-2,:)+(rhs(3,nyf-2,:) &
!    u(3,nyf-2,:)=(rhs(3,nyf-2,:) &
      +fx6*(u(6,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:))+fxx*(u(5,nyf-2,:)+u(1,nyf-2,:))+fx*(u(4,nyf-2,:)+u(2,nyf-2,:)) &
      +fy6*(rbufN(3,1,:)+u(3,nyf-5,:))+fyy*(u(3,nyf,:)+u(3,nyf-4,:))+fy*(u(3,nyf-1,:)+u(3,nyf-3,:)) &
      +fxy*(u(4,nyf-1,:)+u(4,nyf-3,:)+u(2,nyf-3,:)+u(2,nyf-1,:)) &
      +fxxy*(u(5,nyf-1,:)+u(5,nyf-3,:)+u(1,nyf-3,:)+u(1,nyf-1,:)) &
      +fxyy*(u(4,nyf,:)+u(4,nyf-4,:)+u(2,nyf-4,:)+u(2,nyf,:)))*delta
  else
    u(1,nyf-2,:)=Oomehd(qn)*u(1,nyf-2,:)+(rhs(1,nyf-2,:) &
!    u(1,nyf-2,:)=(rhs(1,nyf-2,:) &
      +fx6*(u(4,nyf-2,:)+rbufW(1,nyf-2,:))+fxx*(u(3,nyf-2,:)+rbufW(2,nyf-2,:))+fx*(u(2,nyf-2,:)+rbufW(3,nyf-2,:)) &
      +fy6*(rbufN(1,1,:)+u(1,nyf-5,:))+fyy*(u(1,nyf,:)+u(1,nyf-4,:))+fy*(u(1,nyf-1,:)+u(1,nyf-3,:)) &
      +fxy*(u(2,nyf-1,:)+u(2,nyf-3,:)+rbufW(3,nyf-3,:)+rbufW(3,nyf-1,:)) &
      +fxxy*(u(3,nyf-1,:)+u(3,nyf-3,:)+rbufW(2,nyf-3,:)+rbufW(2,nyf-1,:)) &
      +fxyy*(u(2,nyf,:)+u(2,nyf-4,:)+rbufW(3,nyf-4,:)+rbufW(3,nyf,:)))*delta
    u(2,nyf-2,:)=Oomehd(qn)*u(2,nyf-2,:)+(rhs(2,nyf-2,:) &
!    u(2,nyf-2,:)=(rhs(2,nyf-2,:) &
      +fx6*(u(5,nyf-2,:)+rbufW(2,nyf-2,:))+fxx*(u(4,nyf-2,:)+rbufW(3,nyf-2,:))+fx*(u(3,nyf-2,:)+u(1,nyf-2,:)) &
      +fy6*(rbufN(2,1,:)+u(2,nyf-5,:))+fyy*(u(2,nyf,:)+u(2,nyf-4,:))+fy*(u(2,nyf-1,:)+u(2,nyf-3,:)) &
      +fxy*(u(3,nyf-1,:)+u(3,nyf-3,:)+u(1,nyf-3,:)+u(1,nyf-1,:)) &
      +fxxy*(u(4,nyf-1,:)+u(4,nyf-3,:)+rbufW(3,nyf-3,:)+rbufW(3,nyf-1,:)) &
      +fxyy*(u(3,nyf,:)+u(3,nyf-4,:)+u(1,nyf-4,:)+u(1,nyf,:)))*delta
    u(3,nyf-2,:)=Oomehd(qn)*u(3,nyf-2,:)+(rhs(3,nyf-2,:) &
!    u(3,nyf-2,:)=(rhs(3,nyf-2,:) &
      +fx6*(u(6,nyf-2,:)+rbufW(3,nyf-2,:))+fxx*(u(5,nyf-2,:)+u(1,nyf-2,:))+fx*(u(4,nyf-2,:)+u(2,nyf-2,:)) &
      +fy6*(rbufN(3,1,:)+u(3,nyf-5,:))+fyy*(u(3,nyf,:)+u(3,nyf-4,:))+fy*(u(3,nyf-1,:)+u(3,nyf-3,:)) &
      +fxy*(u(4,nyf-1,:)+u(4,nyf-3,:)+u(2,nyf-3,:)+u(2,nyf-1,:)) &
      +fxxy*(u(5,nyf-1,:)+u(5,nyf-3,:)+u(1,nyf-3,:)+u(1,nyf-1,:)) &
      +fxyy*(u(4,nyf,:)+u(4,nyf-4,:)+u(2,nyf-4,:)+u(2,nyf,:)))*delta
  endif
  do i=4,nxf-3 !.away from x-boundaries
    u(i,nyf-2,:)=Oomehd(qn)*u(i,nyf-2,:)+(rhs(i,nyf-2,:) &
!    u(i,nyf-2,:)=(rhs(i,nyf-2,:) &
      +fx6*(u(i+3,nyf-2,:)+u(i-3,nyf-2,:))+fxx*(u(i+2,nyf-2,:)+u(i-2,nyf-2,:))+fx*(u(i+1,nyf-2,:)+u(i-1,nyf-2,:)) &
      +fy6*(rbufN(i,1,:)+u(i,nyf-5,:))+fyy*(u(i,nyf,:)+u(i,nyf-4,:))+fy*(u(i,nyf-1,:)+u(i,nyf-3,:)) &
      +fxy*(u(i+1,nyf-1,:)+u(i+1,nyf-3,:)+u(i-1,nyf-3,:)+u(i-1,nyf-1,:)) &
      +fxxy*(u(i+2,nyf-1,:)+u(i+2,nyf-3,:)+u(i-2,nyf-3,:)+u(i-2,nyf-1,:)) &
      +fxyy*(u(i+1,nyf,:)+u(i+1,nyf-4,:)+u(i-1,nyf-4,:)+u(i-1,nyf,:)))*delta
  enddo
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-2,nyf-2,:)=Oomehd(qn)*u(nxf-2,nyf-2,:)+(rhs(nxf-2,nyf-2,:) &
!    u(nxf-2,nyf-2,:)=(rhs(nxf-2,nyf-2,:) &
      +fx6*(Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-5,nyf-2,:))+fxx*(u(nxf,nyf-2,:)+u(nxf-4,nyf-2,:))+fx*(u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:)) &
      +fy6*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-5,:))+fyy*(u(nxf-2,nyf,:)+u(nxf-2,nyf-4,:))+fy*(u(nxf-2,nyf-1,:)+u(nxf-2,nyf-3,:)) &
      +fxy*(u(nxf-1,nyf-1,:)+u(nxf-1,nyf-3,:)+u(nxf-3,nyf-3,:)+u(nxf-3,nyf-1,:)) &
      +fxxy*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:)+u(nxf-4,nyf-3,:)+u(nxf-4,nyf-1,:)) &
      +fxyy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-4,:)+u(nxf-3,nyf-4,:)+u(nxf-3,nyf,:)))*delta
    u(nxf-1,nyf-2,:)=Oomehd(qn)*u(nxf-1,nyf-2,:)+(rhs(nxf-1,nyf-2,:) &
!    u(nxf-1,nyf-2,:)=(rhs(nxf-1,nyf-2,:) &
      +fx6*(rs(2,:)+u(nxf-4,nyf-2,:))+fxx*(Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-3,nyf-2,:))+fx*(u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)) &
      +fy6*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-5,:))+fyy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-4,:))+fy*(u(nxf-1,nyf-1,:)+u(nxf-1,nyf-3,:)) &
      +fxy*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:)+u(nxf-2,nyf-3,:)+u(nxf-2,nyf-1,:)) &
      +fxxy*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+Rq*u(nxf,nyf-3,:)+Rs(1,:)+u(nxf-3,nyf-3,:)+u(nxf-3,nyf-1,:)) &
      +fxyy*(u(nxf,nyf,:)+u(nxf,nyf-4,:)+u(nxf-2,nyf-4,:)+u(nxf-2,nyf,:)))*deltaR(2)
    u(nxf,nyf-2,:)=Oomehd(qn)*u(nxf,nyf-2,:)+(rhs(nxf,nyf-2,:) &
!    u(nxf,nyf-2,:)=(rhs(nxf,nyf-2,:) &
      +fx6*(Rq*u(nxf-2,nyf-2,:)+Rs(3,:)+u(nxf-3,nyf-2,:))+fxx*(Rq*u(nxf-1,nyf-2,:)+Rs(2,:)+u(nxf-2,nyf-2,:))+fx*(rs(1,:)+u(nxf-1,nyf-2,:)) &
      +fy6*(rbufN(nxf,1,:)+u(nxf,nyf-5,:))+fyy*(u(nxf,nyf,:)+u(nxf,nyf-4,:))+fy*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:)) &
      +fxy*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+Rq*u(nxf,nyf-3,:)+Rs(1,:)+u(nxf-1,nyf-3,:)+u(nxf-1,nyf-1,:)) &
      +fxxy*(Rq*u(nxf-1,nyf-1,:)+Rs(2,:)+Rq*u(nxf-1,nyf-3,:)+Rs(2,:)+u(nxf-2,nyf-3,:)+u(nxf-2,nyf-1,:)) &
      +fxyy*(Rq*u(nxf,nyf,:)+Rs(1,:)+Rq*u(nxf,nyf-4,:)+Rs(1,:)+u(nxf-1,nyf-4,:)+u(nxf-1,nyf,:)))*deltaR(1)
  else
    u(nxf-2,nyf-2,:)=Oomehd(qn)*u(nxf-2,nyf-2,:)+(rhs(nxf-2,nyf-2,:) &
!    u(nxf-2,nyf-2,:)=(rhs(nxf-2,nyf-2,:) &
      +fx6*(rbufE(1,nyf-2,:)+u(nxf-5,nyf-2,:))+fxx*(u(nxf,nyf-2,:)+u(nxf-4,nyf-2,:))+fx*(u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:)) &
      +fy6*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-5,:))+fyy*(u(nxf-2,nyf,:)+u(nxf-2,nyf-4,:))+fy*(u(nxf-2,nyf-1,:)+u(nxf-2,nyf-3,:)) &
      +fxy*(u(nxf-1,nyf-1,:)+u(nxf-1,nyf-3,:)+u(nxf-3,nyf-3,:)+u(nxf-3,nyf-1,:)) &
      +fxxy*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:)+u(nxf-4,nyf-3,:)+u(nxf-4,nyf-1,:)) &
      +fxyy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-4,:)+u(nxf-3,nyf-4,:)+u(nxf-3,nyf,:)))*delta
    u(nxf-1,nyf-2,:)=Oomehd(qn)*u(nxf-1,nyf-2,:)+(rhs(nxf-1,nyf-2,:) &
!    u(nxf-1,nyf-2,:)=(rhs(nxf-1,nyf-2,:) &
      +fx6*(rbufE(2,nyf-2,:)+u(nxf-4,nyf-2,:))+fxx*(rbufE(1,nyf-2,:)+u(nxf-3,nyf-2,:))+fx*(u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)) &
      +fy6*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-5,:))+fyy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-4,:))+fy*(u(nxf-1,nyf-1,:)+u(nxf-1,nyf-3,:)) &
      +fxy*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:)+u(nxf-2,nyf-3,:)+u(nxf-2,nyf-1,:)) &
      +fxxy*(rbufE(1,nyf-1,:)+rbufE(1,nyf-3,:)+u(nxf-3,nyf-3,:)+u(nxf-3,nyf-1,:)) &
      +fxyy*(u(nxf,nyf,:)+u(nxf,nyf-4,:)+u(nxf-2,nyf-4,:)+u(nxf-2,nyf,:)))*delta
    u(nxf,nyf-2,:)=Oomehd(qn)*u(nxf,nyf-2,:)+(rhs(nxf,nyf-2,:) &
!    u(nxf,nyf-2,:)=(rhs(nxf,nyf-2,:) &
      +fx6*(rbufE(3,nyf-2,:)+u(nxf-3,nyf-2,:))+fxx*(rbufE(2,nyf-2,:)+u(nxf-2,nyf-2,:))+fx*(rbufE(1,nyf-2,:)+u(nxf-1,nyf-2,:)) &
      +fy6*(rbufN(nxf,1,:)+u(nxf,nyf-5,:))+fyy*(u(nxf,nyf,:)+u(nxf,nyf-4,:))+fy*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:)) &
      +fxy*(rbufE(1,nyf-1,:)+rbufE(1,nyf-3,:)+u(nxf-1,nyf-3,:)+u(nxf-1,nyf-1,:)) &
      +fxxy*(rbufE(2,nyf-1,:)+rbufE(2,nyf-3,:)+u(nxf-2,nyf-3,:)+u(nxf-2,nyf-1,:)) &
      +fxyy*(rbufE(1,nyf,:)+rbufE(1,nyf-4,:)+u(nxf-1,nyf-4,:)+u(nxf-1,nyf,:)))*delta
  endif
!........................... SECOND TO LAST ROW ...........................!
  call MPI_WAIT(req(15),stat(:,15),ierr) !.NW corner grid points
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,nyf-1,:)=Oomehd(qn)*u(1,nyf-1,:)+(rhs(1,nyf-1,:) &
!    u(1,nyf-1,:)=(rhs(1,nyf-1,:) &
      +fx6*(u(4,nyf-1,:)+Lq*u(3,nyf-1,:)+Ls(3,:))+fxx*(u(3,nyf-1,:)+Lq*u(2,nyf-1,:)+Ls(2,:))+fx*(u(2,nyf-1,:)+Ls(1,:)) &
      +fy6*(rbufN(1,2,:)+u(1,nyf-4,:))+fyy*(rbufN(1,1,:)+u(1,nyf-3,:))+fy*(u(1,nyf,:)+u(1,nyf-2,:)) &
      +fxy*(u(2,nyf,:)+u(2,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:)+Lq*u(1,nyf,:)+Ls(1,:)) &
      +fxxy*(u(3,nyf,:)+u(3,nyf-2,:)+Lq*u(2,nyf-2,:)+Ls(2,:)+Lq*u(2,nyf,:)+Ls(2,:)) &
      +fxyy*(rbufN(2,1,:)+u(2,nyf-3,:)+Lq*u(1,nyf-3,:)+Ls(1,:)+Lq*rbufN(1,1,:)+Ls(1,:)))*deltaL(1)
    u(2,nyf-1,:)=Oomehd(qn)*u(2,nyf-1,:)+(rhs(2,nyf-1,:) &
!    u(2,nyf-1,:)=(rhs(2,nyf-1,:) &
      +fx6*(u(5,nyf-1,:)+Ls(2,:))+fxx*(u(4,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:))+fx*(u(3,nyf-1,:)+u(1,nyf-1,:)) &
      +fy6*(rbufN(2,2,:)+u(2,nyf-4,:))+fyy*(rbufN(2,1,:)+u(2,nyf-3,:))+fy*(u(2,nyf,:)+u(2,nyf-2,:)) &
      +fxy*(u(3,nyf,:)+u(3,nyf-2,:)+u(1,nyf-2,:)+u(1,nyf,:)) &
      +fxxy*(u(4,nyf,:)+u(4,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:)+Lq*u(1,nyf,:)+Ls(1,:)) &
      +fxyy*(rbufN(3,1,:)+u(3,nyf-3,:)+u(1,nyf-3,:)+rbufN(1,1,:)))*deltaL(2)
    u(3,nyf-1,:)=Oomehd(qn)*u(3,nyf-1,:)+(rhs(3,nyf-1,:) &
!    u(3,nyf-1,:)=(rhs(3,nyf-1,:) &
      +fx6*(u(6,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:))+fxx*(u(5,nyf-1,:)+u(1,nyf-1,:))+fx*(u(4,nyf-1,:)+u(2,nyf-1,:)) &
      +fy6*(rbufN(3,2,:)+u(3,nyf-4,:))+fyy*(rbufN(3,1,:)+u(3,nyf-3,:))+fy*(u(3,nyf,:)+u(3,nyf-2,:)) &
      +fxy*(u(4,nyf,:)+u(4,nyf-2,:)+u(2,nyf-2,:)+u(2,nyf,:)) &
      +fxxy*(u(5,nyf,:)+u(5,nyf-2,:)+u(1,nyf-2,:)+u(1,nyf,:)) &
      +fxyy*(rbufN(4,1,:)+u(4,nyf-3,:)+u(2,nyf-3,:)+rbufN(2,1,:)))*delta
  else
    u(1,nyf-1,:)=Oomehd(qn)*u(1,nyf-1,:)+(rhs(1,nyf-1,:) &
!    u(1,nyf-1,:)=(rhs(1,nyf-1,:) &
      +fx6*(u(4,nyf-1,:)+rbufW(1,nyf-1,:))+fxx*(u(3,nyf-1,:)+rbufW(2,nyf-1,:))+fx*(u(2,nyf-1,:)+rbufW(3,nyf-1,:)) &
      +fy6*(rbufN(1,2,:)+u(1,nyf-4,:))+fyy*(rbufN(1,1,:)+u(1,nyf-3,:))+fy*(u(1,nyf,:)+u(1,nyf-2,:)) &
      +fxy*(u(2,nyf,:)+u(2,nyf-2,:)+rbufW(3,nyf-2,:)+rbufW(3,nyf,:)) &
      +fxxy*(u(3,nyf,:)+u(3,nyf-2,:)+rbufW(2,nyf-2,:)+rbufW(2,nyf,:)) &
      +fxyy*(rbufN(2,1,:)+u(2,nyf-3,:)+rbufW(3,nyf-3,:)+rbufNW(1,:)))*delta
    u(2,nyf-1,:)=Oomehd(qn)*u(2,nyf-1,:)+(rhs(2,nyf-1,:) &
!    u(2,nyf-1,:)=(rhs(2,nyf-1,:) &
      +fx6*(u(5,nyf-1,:)+rbufW(2,nyf-1,:))+fxx*(u(4,nyf-1,:)+rbufW(3,nyf-1,:))+fx*(u(3,nyf-1,:)+u(1,nyf-1,:)) &
      +fy6*(rbufN(2,2,:)+u(2,nyf-4,:))+fyy*(rbufN(2,1,:)+u(2,nyf-3,:))+fy*(u(2,nyf,:)+u(2,nyf-2,:)) &
      +fxy*(u(3,nyf,:)+u(3,nyf-2,:)+u(1,nyf-2,:)+u(1,nyf,:)) &
      +fxxy*(u(4,nyf,:)+u(4,nyf-2,:)+rbufW(3,nyf-2,:)+rbufW(3,nyf,:)) &
      +fxyy*(rbufN(3,1,:)+u(3,nyf-3,:)+u(1,nyf-3,:)+rbufN(1,1,:)))*delta
    u(3,nyf-1,:)=Oomehd(qn)*u(3,nyf-1,:)+(rhs(3,nyf-1,:) &
!    u(3,nyf-1,:)=(rhs(3,nyf-1,:) &
      +fx6*(u(6,nyf-1,:)+rbufW(3,nyf-1,:))+fxx*(u(5,nyf-1,:)+u(1,nyf-1,:))+fx*(u(4,nyf-1,:)+u(2,nyf-1,:)) &
      +fy6*(rbufN(3,2,:)+u(3,nyf-4,:))+fyy*(rbufN(3,1,:)+u(3,nyf-3,:))+fy*(u(3,nyf,:)+u(3,nyf-2,:)) &
      +fxy*(u(4,nyf,:)+u(4,nyf-2,:)+u(2,nyf-2,:)+u(2,nyf,:)) &
      +fxxy*(u(5,nyf,:)+u(5,nyf-2,:)+u(1,nyf-2,:)+u(1,nyf,:)) &
      +fxyy*(rbufN(4,1,:)+u(4,nyf-3,:)+u(2,nyf-3,:)+rbufN(2,1,:)))*delta
  endif
  do i=4,nxf-3 !.away from x-boundaries
    u(i,nyf-1,:)=Oomehd(qn)*u(i,nyf-1,:)+(rhs(i,nyf-1,:) &
!    u(i,nyf-1,:)=(rhs(i,nyf-1,:) &
      +fx6*(u(i+3,nyf-1,:)+u(i-3,nyf-1,:))+fxx*(u(i+2,nyf-1,:)+u(i-2,nyf-1,:))+fx*(u(i+1,nyf-1,:)+u(i-1,nyf-1,:)) &
      +fy6*(rbufN(i,2,:)+u(i,nyf-4,:))+fyy*(rbufN(i,1,:)+u(i,nyf-3,:))+fy*(u(i,nyf,:)+u(i,nyf-2,:)) &
      +fxy*(u(i+1,nyf,:)+u(i+1,nyf-2,:)+u(i-1,nyf-2,:)+u(i-1,nyf,:)) &
      +fxxy*(u(i+2,nyf,:)+u(i+2,nyf-2,:)+u(i-2,nyf-2,:)+u(i-2,nyf,:)) &
      +fxyy*(rbufN(i+1,1,:)+u(i+1,nyf-3,:)+u(i-1,nyf-3,:)+rbufN(i-1,1,:)))*delta
  enddo
  call MPI_WAIT(req(16),stat(:,16),ierr) !.NE corner grid point
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-2,nyf-1,:)=Oomehd(qn)*u(nxf-2,nyf-1,:)+(rhs(nxf-2,nyf-1,:) &
!    u(nxf-2,nyf-1,:)=(rhs(nxf-2,nyf-1,:) &
      +fx6*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-5,nyf-1,:))+fxx*(u(nxf,nyf-1,:)+u(nxf-4,nyf-1,:))+fx*(u(nxf-1,nyf-1,:)+u(nxf-3,nyf-1,:)) &
      +fy6*(rbufN(nxf-2,2,:)+u(nxf-2,nyf-4,:))+fyy*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-3,:))+fy*(u(nxf-2,nyf,:)+u(nxf-2,nyf-2,:)) &
      +fxy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:)+u(nxf-3,nyf,:)) &
      +fxxy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)+u(nxf-4,nyf-2,:)+u(nxf-4,nyf,:)) &
      +fxyy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:)+u(nxf-3,nyf-3,:)+rbufN(nxf-3,1,:)))*delta
    u(nxf-1,nyf-1,:)=Oomehd(qn)*u(nxf-1,nyf-1,:)+(rhs(nxf-1,nyf-1,:) &
!    u(nxf-1,nyf-1,:)=(rhs(nxf-1,nyf-1,:) &
      +fx6*(rs(2,:)+u(nxf-4,nyf-1,:))+fxx*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-3,nyf-1,:))+fx*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)) &
      +fy6*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-4,:))+fyy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:))+fy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:)) &
      +fxy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)+u(nxf-2,nyf,:)) &
      +fxxy*(Rq*u(nxf,nyf,:)+Rs(1,:)+Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-3,nyf-2,:)+u(nxf-3,nyf,:)) &
      +fxyy*(rbufN(nxf,1,:)+u(nxf,nyf-3,:)+u(nxf-2,nyf-3,:)+rbufN(nxf-2,1,:)))*deltaR(2)
    u(nxf,nyf-1,:)=Oomehd(qn)*u(nxf,nyf-1,:)+(rhs(nxf,nyf-1,:) &
!    u(nxf,nyf-1,:)=(rhs(nxf,nyf-1,:) &
      +fx6*(Rq*u(nxf-2,nyf-1,:)+Rs(3,:)+u(nxf-3,nyf-1,:))+fxx*(Rq*u(nxf-1,nyf-1,:)+Rs(2,:)+u(nxf-2,nyf-1,:))+fx*(rs(1,:)+u(nxf-1,nyf-1,:)) &
      +fy6*(rbufN(nxf,2,:)+u(nxf,nyf-4,:))+fyy*(rbufN(nxf,1,:)+u(nxf,nyf-3,:))+fy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)) &
      +fxy*(Rq*u(nxf,nyf,:)+Rs(1,:)+Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-1,nyf-2,:)+u(nxf-1,nyf,:)) &
      +fxxy*(Rq*u(nxf-1,nyf,:)+Rs(2,:)+Rq*u(nxf-1,nyf-2,:)+Rs(2,:)+u(nxf-2,nyf-2,:)+u(nxf-2,nyf,:)) &
      +fxyy*(Rq*rbufN(nxf,1,:)+Rs(1,:)+Rq*u(nxf,nyf-3,:)+Rs(1,:)+u(nxf-1,nyf-3,:)+rbufN(nxf-1,1,:)))*deltaR(1)
  else
    u(nxf-2,nyf-1,:)=Oomehd(qn)*u(nxf-2,nyf-1,:)+(rhs(nxf-2,nyf-1,:) &
!    u(nxf-2,nyf-1,:)=(rhs(nxf-2,nyf-1,:) &
      +fx6*(rbufE(1,nyf-1,:)+u(nxf-5,nyf-1,:))+fxx*(u(nxf,nyf-1,:)+u(nxf-4,nyf-1,:))+fx*(u(nxf-1,nyf-1,:)+u(nxf-3,nyf-1,:)) &
      +fy6*(rbufN(nxf-2,2,:)+u(nxf-2,nyf-4,:))+fyy*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-3,:))+fy*(u(nxf-2,nyf,:)+u(nxf-2,nyf-2,:)) &
      +fxy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:)+u(nxf-3,nyf,:)) &
      +fxxy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)+u(nxf-4,nyf-2,:)+u(nxf-4,nyf,:)) &
      +fxyy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:)+u(nxf-3,nyf-3,:)+rbufN(nxf-3,1,:)))*delta
    u(nxf-1,nyf-1,:)=Oomehd(qn)*u(nxf-1,nyf-1,:)+(rhs(nxf-1,nyf-1,:) &
!    u(nxf-1,nyf-1,:)=(rhs(nxf-1,nyf-1,:) &
      +fx6*(rbufE(2,nyf-1,:)+u(nxf-4,nyf-1,:))+fxx*(rbufE(1,nyf-1,:)+u(nxf-3,nyf-1,:))+fx*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)) &
      +fy6*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-4,:))+fyy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:))+fy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:)) &
      +fxy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)+u(nxf-2,nyf,:)) &
      +fxxy*(rbufE(1,nyf,:)+rbufE(1,nyf-2,:)+u(nxf-3,nyf-2,:)+u(nxf-3,nyf,:)) &
      +fxyy*(rbufN(nxf,1,:)+u(nxf,nyf-3,:)+u(nxf-2,nyf-3,:)+rbufN(nxf-2,1,:)))*delta
    u(nxf,nyf-1,:)=Oomehd(qn)*u(nxf,nyf-1,:)+(rhs(nxf,nyf-1,:) &
!    u(nxf,nyf-1,:)=(rhs(nxf,nyf-1,:) &
      +fx6*(rbufE(3,nyf-1,:)+u(nxf-3,nyf-1,:))+fxx*(rbufE(2,nyf-1,:)+u(nxf-2,nyf-1,:))+fx*(rbufE(1,nyf-1,:)+u(nxf-1,nyf-1,:)) &
      +fy6*(rbufN(nxf,2,:)+u(nxf,nyf-4,:))+fyy*(rbufN(nxf,1,:)+u(nxf,nyf-3,:))+fy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)) &
      +fxy*(rbufE(1,nyf,:)+rbufE(1,nyf-2,:)+u(nxf-1,nyf-2,:)+u(nxf-1,nyf,:)) &
      +fxxy*(rbufE(2,nyf,:)+rbufE(2,nyf-2,:)+u(nxf-2,nyf-2,:)+u(nxf-2,nyf,:)) &
      +fxyy*(rbufNE(1,:)+rbufE(1,nyf-3,:)+u(nxf-1,nyf-3,:)+rbufN(nxf-1,1,:)))*delta
  endif
!........................... TOP TO LAST ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
!...global NW corner
    u(1,nyf,:)=Oomehd(qn)*u(1,nyf,:)+(rhs(1,nyf,:) &
!    u(1,nyf,:)=(rhs(1,nyf,:) &
      +fx6*(u(4,nyf,:)+Lq*u(3,nyf,:)+Ls(3,:))+fxx*(u(3,nyf,:)+Lq*u(2,nyf,:)+Ls(2,:))+fx*(u(2,nyf,:)+Ls(1,:)) &
      +fy6*(rbufN(1,3,:)+u(1,nyf-3,:))+fyy*(rbufN(1,2,:)+u(1,nyf-2,:))+fy*(rbufN(1,1,:)+u(1,nyf-1,:)) &
      +fxy*(rbufN(2,1,:)+u(2,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:)+Lq*rbufN(1,1,:)+Ls(1,:)) &
      +fxxy*(rbufN(3,1,:)+u(3,nyf-1,:)+Lq*u(2,nyf-1,:)+Ls(2,:)+Lq*rbufN(2,1,:)+Ls(2,:)) &
      +fxyy*(rbufN(2,2,:)+u(2,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:)+Lq*rbufN(1,2,:)+Ls(1,:)))*deltaL(1)
    u(2,nyf,:)=Oomehd(qn)*u(2,nyf,:)+(rhs(2,nyf,:) &
!    u(2,nyf,:)=(rhs(2,nyf,:) &
      +fx6*(u(5,nyf,:)+Ls(2,:))+fxx*(u(4,nyf,:)+Lq*u(1,nyf,:)+Ls(1,:))+fx*(u(3,nyf,:)+u(1,nyf,:)) &
      +fy6*(rbufN(2,3,:)+u(2,nyf-3,:))+fyy*(rbufN(2,2,:)+u(2,nyf-2,:))+fy*(rbufN(2,1,:)+u(2,nyf-1,:)) &
      +fxy*(rbufN(3,1,:)+u(3,nyf-1,:)+u(1,nyf-1,:)+rbufN(1,1,:)) &
      +fxxy*(rbufN(4,1,:)+u(4,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:)+Lq*rbufN(1,1,:)+Ls(1,:)) &
      +fxyy*(rbufN(3,2,:)+u(3,nyf-2,:)+u(1,nyf-2,:)+rbufN(1,2,:)))*deltaL(2)
    u(3,nyf,:)=Oomehd(qn)*u(3,nyf,:)+(rhs(3,nyf,:) &
!    u(3,nyf,:)=(rhs(3,nyf,:) &
      +fx6*(u(6,nyf,:)+Lq*u(1,nyf,:)+Ls(1,:))+fxx*(u(5,nyf,:)+u(1,nyf,:))+fx*(u(4,nyf,:)+u(2,nyf,:)) &
      +fy6*(rbufN(3,3,:)+u(3,nyf-3,:))+fyy*(rbufN(3,2,:)+u(3,nyf-2,:))+fy*(rbufN(3,1,:)+u(3,nyf-1,:)) &
      +fxy*(rbufN(4,1,:)+u(4,nyf-1,:)+u(2,nyf-1,:)+rbufN(2,1,:)) &
      +fxxy*(rbufN(5,1,:)+u(5,nyf-1,:)+u(1,nyf-1,:)+rbufN(1,1,:)) &
      +fxyy*(rbufN(4,2,:)+u(4,nyf-2,:)+u(2,nyf-2,:)+rbufN(2,2,:)))*delta
  else
!...local NW corner
    u(1,nyf,:)=Oomehd(qn)*u(1,nyf,:)+(rhs(1,nyf,:) &
!    u(1,nyf,:)=(rhs(1,nyf,:) &
      +fx6*(u(4,nyf,:)+rbufW(1,nyf,:))+fxx*(u(3,nyf,:)+rbufW(2,nyf,:))+fx*(u(2,nyf,:)+rbufW(3,nyf,:)) &
      +fy6*(rbufN(1,3,:)+u(1,nyf-3,:))+fyy*(rbufN(1,2,:)+u(1,nyf-2,:))+fy*(rbufN(1,1,:)+u(1,nyf-1,:)) &
      +fxy*(rbufN(2,1,:)+u(2,nyf-1,:)+rbufW(3,nyf-1,:)+rbufNW(1,:)) &
      +fxxy*(rbufN(3,1,:)+u(3,nyf-1,:)+rbufW(2,nyf-1,:)+rbufNW(3,:)) &
      +fxyy*(rbufN(2,2,:)+u(2,nyf-2,:)+rbufW(3,nyf-2,:)+rbufNW(2,:)))*delta
    u(2,nyf,:)=Oomehd(qn)*u(2,nyf,:)+(rhs(2,nyf,:) &
!    u(2,nyf,:)=(rhs(2,nyf,:) &
      +fx6*(u(5,nyf,:)+rbufW(2,nyf,:))+fxx*(u(4,nyf,:)+rbufW(3,nyf,:))+fx*(u(3,nyf,:)+u(1,nyf,:)) &
      +fy6*(rbufN(2,3,:)+u(2,nyf-3,:))+fyy*(rbufN(2,2,:)+u(2,nyf-2,:))+fy*(rbufN(2,1,:)+u(2,nyf-1,:)) &
      +fxy*(rbufN(3,1,:)+u(3,nyf-1,:)+u(1,nyf-1,:)+rbufN(1,1,:)) &
      +fxxy*(rbufN(4,1,:)+u(4,nyf-1,:)+rbufW(3,nyf-1,:)+rbufNW(1,:)) &
      +fxyy*(rbufN(3,2,:)+u(3,nyf-2,:)+u(1,nyf-2,:)+rbufN(1,2,:)))*delta
    u(3,nyf,:)=Oomehd(qn)*u(3,nyf,:)+(rhs(3,nyf,:) &
!    u(3,nyf,:)=(rhs(3,nyf,:) &
      +fx6*(u(6,nyf,:)+rbufW(3,nyf,:))+fxx*(u(5,nyf,:)+u(1,nyf,:))+fx*(u(4,nyf,:)+u(2,nyf,:)) &
      +fy6*(rbufN(3,3,:)+u(3,nyf-3,:))+fyy*(rbufN(3,2,:)+u(3,nyf-2,:))+fy*(rbufN(3,1,:)+u(3,nyf-1,:)) &
      +fxy*(rbufN(4,1,:)+u(4,nyf-1,:)+u(2,nyf-1,:)+rbufN(2,1,:)) &
      +fxxy*(rbufN(5,1,:)+u(5,nyf-1,:)+u(1,nyf-1,:)+rbufN(1,1,:)) &
      +fxyy*(rbufN(4,2,:)+u(4,nyf-2,:)+u(2,nyf-2,:)+rbufN(2,2,:)))*delta
  endif
  do i=4,nxf-3 !.away from x-boundaries
    u(i,nyf,:)=Oomehd(qn)*u(i,nyf,:)+(rhs(i,nyf,:) &
!    u(i,nyf,:)=(rhs(i,nyf,:) &
      +fx6*(u(i+3,nyf,:)+u(i-3,nyf,:))+fxx*(u(i+2,nyf,:)+u(i-2,nyf,:))+fx*(u(i+1,nyf,:)+u(i-1,nyf,:)) &
      +fy6*(rbufN(i,3,:)+u(i,nyf-3,:))+fyy*(rbufN(i,2,:)+u(i,nyf-2,:))+fy*(rbufN(i,1,:)+u(i,nyf-1,:)) &
      +fxy*(rbufN(i+1,1,:)+u(i+1,nyf-1,:)+u(i-1,nyf-1,:)+rbufN(i-1,1,:)) &
      +fxxy*(rbufN(i+2,1,:)+u(i+2,nyf-1,:)+u(i-2,nyf-1,:)+rbufN(i-2,1,:)) &
      +fxyy*(rbufN(i+1,2,:)+u(i+1,nyf-2,:)+u(i-1,nyf-2,:)+rbufN(i-1,2,:)))*delta
  enddo
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-2,nyf,:)=Oomehd(qn)*u(nxf-2,nyf,:)+(rhs(nxf-2,nyf,:) &
!    u(nxf-2,nyf,:)=(rhs(nxf-2,nyf,:) &
      +fx6*(Rq*u(nxf,nyf,:)+Rs(1,:)+u(nxf-5,nyf,:))+fxx*(u(nxf,nyf,:)+u(nxf-4,nyf,:))+fx*(u(nxf-1,nyf,:)+u(nxf-3,nyf,:)) &
      +fy6*(rbufN(nxf-2,3,:)+u(nxf-2,nyf-3,:))+fyy*(rbufN(nxf-2,2,:)+u(nxf-2,nyf-2,:))+fy*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-1,:)) &
      +fxy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:)+u(nxf-3,nyf-1,:)+rbufN(nxf-3,1,:)) &
      +fxxy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)+u(nxf-4,nyf-1,:)+rbufN(nxf-4,1,:)) &
      +fxyy*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:)+rbufN(nxf-3,2,:)))*delta
    u(nxf-1,nyf,:)=Oomehd(qn)*u(nxf-1,nyf,:)+(rhs(nxf-1,nyf,:) &
!    u(nxf-1,nyf,:)=(rhs(nxf-1,nyf,:) &
      +fx6*(rs(2,:)+u(nxf-4,nyf,:))+fxx*(Rq*u(nxf,nyf,:)+Rs(1,:)+u(nxf-3,nyf,:))+fx*(u(nxf,nyf,:)+u(nxf-2,nyf,:)) &
      +fy6*(rbufN(nxf-1,3,:)+u(nxf-1,nyf-3,:))+fyy*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:))+fy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:)) &
      +fxy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)+rbufN(nxf-2,1,:)) &
      +fxxy*(Rq*rbufN(nxf,1,:)+Rs(1,:)+Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-3,nyf-1,:)+rbufN(nxf-3,1,:)) &
      +fxyy*(rbufN(nxf,2,:)+u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)+rbufN(nxf-2,2,:)))*deltaR(2)
!...global NE corner
    u(nxf,nyf,:)=Oomehd(qn)*u(nxf,nyf,:)+(rhs(nxf,nyf,:) &
!    u(nxf,nyf,:)=(rhs(nxf,nyf,:) &
      +fx6*(Rq*u(nxf-2,nyf,:)+Rs(3,:)+u(nxf-3,nyf,:))+fxx*(Rq*u(nxf-1,nyf,:)+Rs(2,:)+u(nxf-2,nyf,:))+fx*(rs(1,:)+u(nxf-1,nyf,:)) &
      +fy6*(rbufN(nxf,3,:)+u(nxf,nyf-3,:))+fyy*(rbufN(nxf,2,:)+u(nxf,nyf-2,:))+fy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)) &
      +fxy*(Rq*rbufN(nxf,1,:)+Rs(1,:)+Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-1,nyf-1,:)+rbufN(nxf-1,1,:)) &
      +fxxy*(Rq*rbufN(nxf-1,1,:)+Rs(2,:)+Rq*u(nxf-1,nyf-1,:)+Rs(2,:)+u(nxf-2,nyf-1,:)+rbufN(nxf-2,1,:)) &
      +fxyy*(Rq*rbufN(nxf,2,:)+Rs(1,:)+Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-1,nyf-2,:)+rbufN(nxf-1,2,:)))*deltaR(1)
  else
    u(nxf-2,nyf,:)=Oomehd(qn)*u(nxf-2,nyf,:)+(rhs(nxf-2,nyf,:) &
!    u(nxf-2,nyf,:)=(rhs(nxf-2,nyf,:) &
      +fx6*(rbufE(1,nyf,:)+u(nxf-5,nyf,:))+fxx*(u(nxf,nyf,:)+u(nxf-4,nyf,:))+fx*(u(nxf-1,nyf,:)+u(nxf-3,nyf,:)) &
      +fy6*(rbufN(nxf-2,3,:)+u(nxf-2,nyf-3,:))+fyy*(rbufN(nxf-2,2,:)+u(nxf-2,nyf-2,:))+fy*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-1,:)) &
      +fxy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:)+u(nxf-3,nyf-1,:)+rbufN(nxf-3,1,:)) &
      +fxxy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)+u(nxf-4,nyf-1,:)+rbufN(nxf-4,1,:)) &
      +fxyy*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:)+rbufN(nxf-3,2,:)))*delta
    u(nxf-1,nyf,:)=Oomehd(qn)*u(nxf-1,nyf,:)+(rhs(nxf-1,nyf,:) &
!    u(nxf-1,nyf,:)=(rhs(nxf-1,nyf,:) &
      +fx6*(rbufE(2,nyf,:)+u(nxf-4,nyf,:))+fxx*(rbufE(1,nyf,:)+u(nxf-3,nyf,:))+fx*(u(nxf,nyf,:)+u(nxf-2,nyf,:)) &
      +fy6*(rbufN(nxf-1,3,:)+u(nxf-1,nyf-3,:))+fyy*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:))+fy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:)) &
      +fxy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)+rbufN(nxf-2,1,:)) &
      +fxxy*(rbufNE(1,:)+rbufE(1,nyf-1,:)+u(nxf-3,nyf-1,:)+rbufN(nxf-3,1,:)) &
      +fxyy*(rbufN(nxf,2,:)+u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)+rbufN(nxf-2,2,:)))*delta
!...local NE corner
    u(nxf,nyf,:)=Oomehd(qn)*u(nxf,nyf,:)+(rhs(nxf,nyf,:) &
!    u(nxf,nyf,:)=(rhs(nxf,nyf,:) &
      +fx6*(rbufE(3,nyf,:)+u(nxf-3,nyf,:))+fxx*(rbufE(2,nyf,:)+u(nxf-2,nyf,:))+fx*(rbufE(1,nyf,:)+u(nxf-1,nyf,:)) &
      +fy6*(rbufN(nxf,3,:)+u(nxf,nyf-3,:))+fyy*(rbufN(nxf,2,:)+u(nxf,nyf-2,:))+fy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)) &
      +fxy*(rbufNE(1,:)+rbufE(1,nyf-1,:)+u(nxf-1,nyf-1,:)+rbufN(nxf-1,1,:)) &
      +fxxy*(rbufNE(3,:)+rbufE(2,nyf-1,:)+u(nxf-2,nyf-1,:)+rbufN(nxf-2,1,:)) &
      +fxyy*(rbufNE(2,:)+rbufE(1,nyf-2,:)+u(nxf-1,nyf-2,:)+rbufN(nxf-1,2,:)))*delta
  endif
#ENDIF

  call MPI_WAITALL(8,(/req(1),req(2),req(3),req(4),req(9),req(10),req(11),req(12)/), &
    (/stat(:,1),stat(:,2),stat(:,3),stat(:,4),stat(:,9),stat(:,10),stat(:,11),stat(:,12)/),ierr)

  end subroutine relaxhd
!*********************************************************
  function reshd(ng,u,rhs,xBC,qn)
!.Returns minus the residual. Input quantities
!.are u, av, rhs, while the residual is returned in reshd.
!.reshd,u,rhs are (nxf,nyf) arrays
  implicit none
  include 'mpif.h'
  real(DP), allocatable :: reshd(:,:,:)
  real(DP), intent(in) :: u(:,:,:),rhs(:,:,:)
  integer(I4B), intent(in) :: ng,xBC(:),qn
  integer(I4B) :: i,j,k,Nrank,Srank,Erank,Wrank,NErank,SErank,SWrank,NWrank
  real(DP) :: rdx,rdy,fx,fxx,fy,fyy,fxy,fc,Lq,Rq
#IF (HDop==3)
  real(DP) :: fx6,fy6,fx4,fy4,fxxy,fxyy
#ENDIF
  real(DP), allocatable :: Ls(:,:),rs(:,:),loc2D(:,:)
  real(DP), allocatable, dimension(:,:,:) :: sbufW,rbufE,sbufE,rbufW, &
                                             sbufN,rbufS,sbufS,rbufN
  real(DP), allocatable, dimension(:,:) :: sbufNE,sbufSE,sbufSW,sbufNW, &
                                           rbufNE,rbufSE,rbufSW,rbufNW
  integer(I4B) :: req(16),stat(MPI_STATUS_SIZE,16),ierr

#IF (MGTIMER==1)
    call CPU_TIME(rest1)
#ENDIF

!.Neigboring ranks to communicate with
  Nrank  = nborhd(qn)%a(ng,1)
  NErank = nborhd(qn)%a(ng,2)
  Erank  = nborhd(qn)%a(ng,3)
  SErank = nborhd(qn)%a(ng,4)
  Srank  = nborhd(qn)%a(ng,5)
  SWrank = nborhd(qn)%a(ng,6)
  Wrank  = nborhd(qn)%a(ng,7)
  NWrank = nborhd(qn)%a(ng,8)

  allocate(reshd(nxf,nyf,nzL))
!.Boundary conditions are implemented with some amount of flexibility.
!.The right boundary for example, if u(nxf+1) is needed, it will be
!.written as
!.  u(nxf+1) = Rq*u(nxf) + Rr*u(nxf-1) + Rs
!.so that symmetry, extrapolation or Dirichlet conditions can be
!.easily imposed.
!......... BCs of u ................................!
  allocate(Ls(HDop,nzL),rs(HDop,nzL),loc2D(HDop,nzL))
!.LHS BC of u
  if (xBC(1)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Lq = -1.0_dp
    forall(i=1:HDop,k=1:nzL) loc2D(i,k)=sum(u(i,:,k))
    call MPI_ALLreduce(loc2D,Ls,HDop*nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    ls =  2.0_dp*Ls/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC(1)==1 .OR. xBC(1)==-1) then ! symmetry
    Lq = dble(xBC(1))
    Ls = 0.0_dp
  endif
!.RHS BC of u
  if (xBC(2)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Rq = -1.0_dp
    forall(i=1:HDop,k=1:nzL) loc2D(i,k)=sum(u(nxf-i+1,:,k))
    call MPI_ALLreduce(loc2D,Rs,HDop*nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    Rs =  2.0_dp*Rs/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC(2)==1 .OR. xBC(2)==-1) then ! symmetry
    Rq = dble(xBC(2))
    Rs = 0.0_dp
  endif
!........................................................!

  allocate(sbufN(nxf,HDop,nzL),rbufS(nxf,HDop,nzL),sbufE(HDop,nyf,nzL),rbufW(HDop,nyf,nzL))
  allocate(sbufS(nxf,HDop,nzL),rbufN(nxf,HDop,nzL),sbufW(HDop,nyf,nzL),rbufE(HDop,nyf,nzL))
  allocate(sbufNE(ndiag,nzL),sbufSE(ndiag,nzL),sbufSW(ndiag,nzL),sbufNW(ndiag,nzL))
  allocate(rbufNE(ndiag,nzL),rbufSE(ndiag,nzL),rbufSW(ndiag,nzL),rbufNW(ndiag,nzL))
!.Send top 2 rows of u to North rank
  sbufN=u(:,nyf-HDop+1:nyf,:)
  call MPI_ISEND(sbufN,nxf*HDop*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(1),ierr)
!.Send bottom 2 rows of u to South rank
  sbufS=u(:,1:HDop,:)
  call MPI_ISEND(sbufS,nxf*HDop*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(2),ierr)
!.Send right 2 columns of u to East rank
  sbufE=u(nxf-HDop+1:nxf,:,:)
  call MPI_ISEND(sbufE,HDop*nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(3),ierr)
!.Send 2 left columns of u to West rank
  sbufW=u(1:HDop,:,:)
  call MPI_ISEND(sbufW,HDop*nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(4),ierr)
#IF (HDop==2)
  sbufNE(1,:)=u(nxf,nyf,:)
  sbufSE(1,:)=u(nxf,1,:)
  sbufSW(1,:)=u(1,1,:)
  sbufNW(1,:)=u(1,nyf,:)
#ELSE
!.First element is the directly diagonal also needed for HDop==2
!.The other two follow the clockwise order, starting at 12 o'clock
  forall(k=1:nzL) sbufNE(:,k)=(/u(nxf,nyf,k),u(nxf,nyf-1,k),u(nxf-1,nyf,k)/)
  forall(k=1:nzL) sbufSE(:,k)=(/u(nxf,1,k),u(nxf,2,k),u(nxf-1,1,k)/)
  forall(k=1:nzL) sbufSW(:,k)=(/u(1,1,k),u(1,2,k),u(2,1,k)/)
  forall(k=1:nzL) sbufNW(:,k)=(/u(1,nyf,k),u(2,nyf,k),u(1,nyf-1,k)/)
#ENDIF
!.Send NE corner to NE rank, receive from SW rank
  call MPI_ISEND(sbufNE,ndiag*nzL,MPI_DOUBLE_PRECISION,NErank, &
                 0,COMMhd(qn)%a(ng,3),req(9),ierr)
!.Send SE corner to SE rank, receive from NW rank
  call MPI_ISEND(sbufSE,ndiag*nzL,MPI_DOUBLE_PRECISION,SErank, &
                 0,COMMhd(qn)%a(ng,3),req(10),ierr)
!.Send SW corner to SW rank, receive from NE rank
  call MPI_ISEND(sbufSW,ndiag*nzL,MPI_DOUBLE_PRECISION,SWrank, &
                 0,COMMhd(qn)%a(ng,3),req(11),ierr)
!.Send NW corner to NW rank, receive from SE rank
  call MPI_ISEND(sbufNW,ndiag*nzL,MPI_DOUBLE_PRECISION,NWrank, &
                 0,COMMhd(qn)%a(ng,3),req(12),ierr)
!.Receive u rows from South rank, for u(i,j-2) and u(i,j-1)
  call MPI_IRECV(rbufS,nxf*HDop*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(5),ierr)
!.Receive u rows from North rank, for u(i,j+2) and u(i,j+1)
  call MPI_IRECV(rbufN,nxf*HDop*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(6),ierr)
!.Receive u columns from West rank, for u(i-2,j) and u(i-1,j)
  call MPI_IRECV(rbufW,HDop*nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(7),ierr)
!.Receive u columns from East rank, for u(i+2,j) and u(i+1,j)
  call MPI_IRECV(rbufE,HDop*nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(8),ierr)
!.Receive NE corner from SW rank, for u(i-1,j-1)
  call MPI_IRECV(rbufSW,ndiag*nzL,MPI_DOUBLE_PRECISION,SWrank, &
                 0,COMMhd(qn)%a(ng,3),req(13),ierr)
!.Receive NW corner from SE rank, for u(i+1,j-1)
  call MPI_IRECV(rbufSE,ndiag*nzL,MPI_DOUBLE_PRECISION,SErank, &
                 0,COMMhd(qn)%a(ng,3),req(14),ierr)
!.Receive SE corner from NW rank, for u(i-1,j+1)
  call MPI_IRECV(rbufNW,ndiag*nzL,MPI_DOUBLE_PRECISION,NWrank, &
                 0,COMMhd(qn)%a(ng,3),req(15),ierr)
!.Receive SW corner from NE rank, for u(i+1,j+1)
  call MPI_IRECV(rbufNE,ndiag*nzL,MPI_DOUBLE_PRECISION,NErank, &
                 0,COMMhd(qn)%a(ng,3),req(16),ierr)

  rdx = dble(nxf*iprocshd(qn)%a(ng))/bLx
  rdy = dble(nyf*jprocshd(qn)%a(ng))/bLy
#IF (HDop==2)
  fxx = hDiff(1)*(rdx**4)
  fyy = hDiff(2)*(rdy**4)
  fxy = (hDiff(1)+hDiff(2))*((rdx*rdy)**2)
  fx  = -4.0_dp*fxx-2.0_dp*fxy
  fy  = -4.0_dp*fyy-2.0_dp*fxy
  fc  = 1.0_dp+2.0_dp*(3.0_dp*(fxx+fyy)+2.0_dp*fxy)

!.Interior points
  do j=3,nyf-2
    do i=3,nxf-2
      reshd(i,j,:)=rhs(i,j,:)-(fxx*(u(i+2,j,:)+u(i-2,j,:))+fx*(u(i+1,j,:)+u(i-1,j,:)) &
        +fyy*(u(i,j+2,:)+u(i,j-2,:))+fy*(u(i,j+1,:)+u(i,j-1,:)) &
        +fxy*(u(i+1,j+1,:)+u(i+1,j-1,:)+u(i-1,j-1,:)+u(i-1,j+1,:))+fc*u(i,j,:))
    enddo
  enddo

!.Boundary points
  call MPI_WAIT(req(5),stat(:,5),ierr)
  do i=3,nxf-2 !.Bottom boundary
    reshd(i,1,:)=rhs(i,1,:)-(fxx*(u(i+2,1,:)+u(i-2,1,:))+fx*(u(i+1,1,:)+u(i-1,1,:)) &
      +fyy*(u(i,3,:)+rbufS(i,1,:))+fy*(u(i,2,:)+rbufS(i,2,:)) &
      +fxy*(u(i+1,2,:)+rbufS(i+1,2,:)+rbufS(i-1,2,:)+u(i-1,2,:))+fc*u(i,1,:))
    reshd(i,2,:)=rhs(i,2,:)-(fxx*(u(i+2,2,:)+u(i-2,2,:))+fx*(u(i+1,2,:)+u(i-1,2,:)) &
      +fyy*(u(i,4,:)+rbufS(i,2,:))+fy*(u(i,3,:)+u(i,1,:)) &
      +fxy*(u(i+1,3,:)+u(i+1,1,:)+u(i-1,1,:)+u(i-1,3,:))+fc*u(i,2,:))
  enddo
  call MPI_WAIT(req(6),stat(:,6),ierr)
  do i=3,nxf-2 !.Top boundary
    reshd(i,nyf-1,:)=rhs(i,nyf-1,:)-(fxx*(u(i+2,nyf-1,:)+u(i-2,nyf-1,:))+fx*(u(i+1,nyf-1,:)+u(i-1,nyf-1,:)) &
      +fyy*(rbufN(i,1,:)+u(i,nyf-3,:))+fy*(u(i,nyf,:)+u(i,nyf-2,:)) &
      +fxy*(u(i+1,nyf,:)+u(i+1,nyf-2,:)+u(i-1,nyf-2,:)+u(i-1,nyf,:))+fc*u(i,nyf-1,:))
    reshd(i,nyf,:)=rhs(i,nyf,:)-(fxx*(u(i+2,nyf,:)+u(i-2,nyf,:))+fx*(u(i+1,nyf,:)+u(i-1,nyf,:)) &
      +fyy*(rbufN(i,2,:)+u(i,nyf-2,:))+fy*(rbufN(i,1,:)+u(i,nyf-1,:)) &
      +fxy*(rbufN(i+1,1,:)+u(i+1,nyf-1,:)+u(i-1,nyf-1,:)+rbufN(i-1,1,:))+fc*u(i,nyf,:))
  enddo
!.Left boundary
  call MPI_WAITALL(3,(/req(7),req(13),req(15)/),(/stat(:,7),stat(:,13),stat(:,15)/),ierr)
  if (iIDhd(qn)%a(ng) == 0) then !.Left most process (along x)
!...global SW corner
    reshd(1,1,:)=rhs(1,1,:)-(fxx*(u(3,1,:)+Lq*u(2,1,:)+Ls(2,:))+fx*(u(2,1,:)+Lq*u(1,1,:)+Ls(1,:)) &
      +fyy*(u(1,3,:)+rbufS(1,1,:))+fy*(u(1,2,:)+rbufS(1,2,:)) &
      +fxy*(u(2,2,:)+rbufS(2,2,:)+Lq*rbufS(1,2,:)+Ls(1,:)+Lq*u(1,2,:)+Ls(1,:))+fc*u(1,1,:))
    reshd(2,1,:)=rhs(2,1,:)-(fxx*(u(4,1,:)+Lq*u(1,1,:)+Ls(1,:))+fx*(u(3,1,:)+u(1,1,:)) &
      +fyy*(u(2,3,:)+rbufS(2,1,:))+fy*(u(2,2,:)+rbufS(2,2,:)) &
      +fxy*(u(3,2,:)+rbufS(3,2,:)+rbufS(1,2,:)+u(1,2,:))+fc*u(2,1,:))
    reshd(1,2,:)=rhs(1,2,:)-(fxx*(u(3,2,:)+Lq*u(2,2,:)+Ls(2,:))+fx*(u(2,2,:)+Lq*u(1,2,:)+Ls(1,:)) &
      +fyy*(u(1,4,:)+rbufS(1,2,:))+fy*(u(1,3,:)+u(1,1,:)) &
      +fxy*(u(2,3,:)+u(2,1,:)+Lq*u(1,1,:)+Ls(1,:)+Lq*u(1,3,:)+Ls(1,:))+fc*u(1,2,:))
    reshd(2,2,:)=rhs(2,2,:)-(fxx*(u(4,2,:)+Lq*u(1,2,:)+Ls(1,:))+fx*(u(3,2,:)+u(1,2,:)) &
      +fyy*(u(2,4,:)+rbufS(2,2,:))+fy*(u(2,3,:)+u(2,1,:)) &
      +fxy*(u(3,3,:)+u(3,1,:)+u(1,1,:)+u(1,3,:))+fc*u(2,2,:))
    do j=3,nyf-2 !.away from x-boundaries
      reshd(1,j,:)=rhs(1,j,:)-(fxx*(u(3,j,:)+Lq*u(2,j,:)+Ls(2,:))+fx*(u(2,j,:)+Lq*u(1,j,:)+Ls(1,:)) &
        +fyy*(u(1,j+2,:)+u(1,j-2,:))+fy*(u(1,j+1,:)+u(1,j-1,:)) &
        +fxy*(u(2,j+1,:)+u(2,j-1,:)+Lq*u(1,j-1,:)+Ls(1,:)+Lq*u(1,j+1,:)+Ls(1,:))+fc*u(1,j,:))
      reshd(2,j,:)=rhs(2,j,:)-(fxx*(u(4,j,:)+Lq*u(1,j,:)+Ls(1,:))+fx*(u(3,j,:)+u(1,j,:)) &
        +fyy*(u(2,j+2,:)+u(2,j-2,:))+fy*(u(2,j+1,:)+u(2,j-1,:)) &
        +fxy*(u(3,j+1,:)+u(3,j-1,:)+u(1,j-1,:)+u(1,j+1,:))+fc*u(2,j,:))
    enddo
!...global NW corner
    reshd(1,nyf-1,:)=rhs(1,nyf-1,:)-(fxx*(u(3,nyf-1,:)+Lq*u(2,nyf-1,:)+Ls(2,:))+fx*(u(2,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:)) &
      +fyy*(rbufN(1,1,:)+u(1,nyf-3,:))+fy*(u(1,nyf,:)+u(1,nyf-2,:)) &
      +fxy*(u(2,nyf,:)+u(2,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:)+Lq*u(1,nyf,:)+Ls(1,:))+fc*u(1,nyf-1,:))
    reshd(2,nyf-1,:)=rhs(2,nyf-1,:)-(fxx*(u(4,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:))+fx*(u(3,nyf-1,:)+u(1,nyf-1,:)) &
      +fyy*(rbufN(2,1,:)+u(2,nyf-3,:))+fy*(u(2,nyf,:)+u(2,nyf-2,:)) &
      +fxy*(u(3,nyf,:)+u(3,nyf-2,:)+u(1,nyf-2,:)+u(1,nyf,:))+fc*u(2,nyf-1,:))
    reshd(1,nyf,:)=rhs(1,nyf,:)-(fxx*(u(3,nyf,:)+Lq*u(2,nyf,:)+Ls(2,:))+fx*(u(2,nyf,:)+Lq*u(1,nyf,:)+Ls(1,:)) &
      +fyy*(rbufN(1,2,:)+u(1,nyf-2,:))+fy*(rbufN(1,1,:)+u(1,nyf-1,:)) &
      +fxy*(rbufN(2,1,:)+u(2,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:)+Lq*rbufN(1,1,:)+Ls(1,:))+fc*u(1,nyf,:))
    reshd(2,nyf,:)=rhs(2,nyf,:)-(fxx*(u(4,nyf,:)+Lq*u(1,nyf,:)+Ls(1,:))+fx*(u(3,nyf,:)+u(1,nyf,:)) &
      +fyy*(rbufN(2,2,:)+u(2,nyf-2,:))+fy*(rbufN(2,1,:)+u(2,nyf-1,:)) &
      +fxy*(rbufN(3,1,:)+u(3,nyf-1,:)+u(1,nyf-1,:)+rbufN(1,1,:))+fc*u(2,nyf,:))
  else
!...local SW corner
    reshd(1,1,:)=rhs(1,1,:)-(fxx*(u(3,1,:)+rbufW(1,1,:))+fx*(u(2,1,:)+rbufW(2,1,:)) &
      +fyy*(u(1,3,:)+rbufS(1,1,:))+fy*(u(1,2,:)+rbufS(1,2,:)) &
      +fxy*(u(2,2,:)+rbufS(2,2,:)+rbufSW(1,:)+rbufW(2,2,:))+fc*u(1,1,:))
    reshd(2,1,:)=rhs(2,1,:)-(fxx*(u(4,1,:)+rbufW(2,1,:))+fx*(u(3,1,:)+u(1,1,:)) &
      +fyy*(u(2,3,:)+rbufS(2,1,:))+fy*(u(2,2,:)+rbufS(2,2,:)) &
      +fxy*(u(3,2,:)+rbufS(3,2,:)+rbufS(1,2,:)+u(1,2,:))+fc*u(2,1,:))
    reshd(1,2,:)=rhs(1,2,:)-(fxx*(u(3,2,:)+rbufW(1,2,:))+fx*(u(2,2,:)+rbufW(2,2,:)) &
      +fyy*(u(1,4,:)+rbufS(1,2,:))+fy*(u(1,3,:)+u(1,1,:)) &
      +fxy*(u(2,3,:)+u(2,1,:)+rbufW(2,1,:)+rbufW(2,3,:))+fc*u(1,2,:))
    reshd(2,2,:)=rhs(2,2,:)-(fxx*(u(4,2,:)+rbufW(2,2,:))+fx*(u(3,2,:)+u(1,2,:)) &
      +fyy*(u(2,4,:)+rbufS(2,2,:))+fy*(u(2,3,:)+u(2,1,:)) &
      +fxy*(u(3,3,:)+u(3,1,:)+u(1,1,:)+u(1,3,:))+fc*u(2,2,:))
    do j=3,nyf-2
      reshd(1,j,:)=rhs(1,j,:)-(fxx*(u(3,j,:)+rbufW(1,j,:))+fx*(u(2,j,:)+rbufW(2,j,:)) &
        +fyy*(u(1,j+2,:)+u(1,j-2,:))+fy*(u(1,j+1,:)+u(1,j-1,:)) &
        +fxy*(u(2,j+1,:)+u(2,j-1,:)+rbufW(2,j-1,:)+rbufW(2,j+1,:))+fc*u(1,j,:))
      reshd(2,j,:)=rhs(2,j,:)-(fxx*(u(4,j,:)+rbufW(2,j,:))+fx*(u(3,j,:)+u(1,j,:)) &
        +fyy*(u(2,j+2,:)+u(2,j-2,:))+fy*(u(2,j+1,:)+u(2,j-1,:)) &
        +fxy*(u(3,j+1,:)+u(3,j-1,:)+u(1,j-1,:)+u(1,j+1,:))+fc*u(2,j,:))
    enddo
!...local NW corner
    reshd(1,nyf-1,:)=rhs(1,nyf-1,:)-(fxx*(u(3,nyf-1,:)+rbufW(1,nyf-1,:))+fx*(u(2,nyf-1,:)+rbufW(2,nyf-1,:)) &
      +fyy*(rbufN(1,1,:)+u(1,nyf-3,:))+fy*(u(1,nyf,:)+u(1,nyf-2,:)) &
      +fxy*(u(2,nyf,:)+u(2,nyf-2,:)+rbufW(2,nyf-2,:)+rbufW(2,nyf,:))+fc*u(1,nyf-1,:))
    reshd(2,nyf-1,:)=rhs(2,nyf-1,:)-(fxx*(u(4,nyf-1,:)+rbufW(2,nyf-1,:))+fx*(u(3,nyf-1,:)+u(1,nyf-1,:)) &
      +fyy*(rbufN(2,1,:)+u(2,nyf-3,:))+fy*(u(2,nyf,:)+u(2,nyf-2,:)) &
      +fxy*(u(3,nyf,:)+u(3,nyf-2,:)+u(1,nyf-2,:)+u(1,nyf,:))+fc*u(2,nyf-1,:))
    reshd(1,nyf,:)=rhs(1,nyf,:)-(fxx*(u(3,nyf,:)+rbufW(1,nyf,:))+fx*(u(2,nyf,:)+rbufW(2,nyf,:)) &
      +fyy*(rbufN(1,2,:)+u(1,nyf-2,:))+fy*(rbufN(1,1,:)+u(1,nyf-1,:)) &
      +fxy*(rbufN(2,1,:)+u(2,nyf-1,:)+rbufW(2,nyf-1,:)+rbufNW(1,:))+fc*u(1,nyf,:))
    reshd(2,nyf,:)=rhs(2,nyf,:)-(fxx*(u(4,nyf,:)+rbufW(2,nyf,:))+fx*(u(3,nyf,:)+u(1,nyf,:)) &
      +fyy*(rbufN(2,2,:)+u(2,nyf-2,:))+fy*(rbufN(2,1,:)+u(2,nyf-1,:)) &
      +fxy*(rbufN(3,1,:)+u(3,nyf-1,:)+u(1,nyf-1,:)+rbufN(1,1,:))+fc*u(2,nyf,:))
  endif
!.Right boundary
  call MPI_WAITALL(3,(/req(8),req(14),req(16)/),(/stat(:,8),stat(:,14),stat(:,16)/),ierr)
  if (iIDhd(qn)%a(ng) == iprocsm1) then !.Right most process (along x)
!...SE corner
    reshd(nxf-1,1,:)=rhs(nxf-1,1,:)-(fxx*(Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-3,1,:))+fx*(u(nxf,1,:)+u(nxf-2,1,:)) &
      +fyy*(u(nxf-1,3,:)+rbufS(nxf-1,1,:))+fy*(u(nxf-1,2,:)+rbufS(nxf-1,2,:)) &
      +fxy*(u(nxf,2,:)+rbufS(nxf,2,:)+rbufS(nxf-2,2,:)+u(nxf-2,2,:))+fc*u(nxf-1,1,:))
    reshd(nxf,1,:)=rhs(nxf,1,:)-(fxx*(Rq*u(nxf-1,1,:)+Rs(2,:)+u(nxf-2,1,:))+fx*(Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-1,1,:)) &
      +fyy*(u(nxf,3,:)+rbufS(nxf,1,:))+fy*(u(nxf,2,:)+rbufS(nxf,2,:)) &
      +fxy*(Rq*u(nxf,2,:)+Rs(1,:)+Rq*rbufS(nxf,2,:)+Rs(1,:)+rbufS(nxf-1,2,:)+u(nxf-1,2,:))+fc*u(nxf,1,:))
    reshd(nxf-1,2,:)=rhs(nxf-1,2,:)-(fxx*(Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-3,2,:))+fx*(u(nxf,2,:)+u(nxf-2,2,:)) &
      +fyy*(u(nxf-1,4,:)+rbufS(nxf-1,2,:))+fy*(u(nxf-1,3,:)+u(nxf-1,1,:)) &
      +fxy*(u(nxf,3,:)+u(nxf,1,:)+u(nxf-2,1,:)+u(nxf-2,3,:))+fc*u(nxf-1,2,:))
    reshd(nxf,2,:)=rhs(nxf,2,:)-(fxx*(Rq*u(nxf-1,2,:)+Rs(2,:)+u(nxf-2,2,:))+fx*(Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-1,2,:)) &
      +fyy*(u(nxf,4,:)+rbufS(nxf,2,:))+fy*(u(nxf,3,:)+u(nxf,1,:)) &
      +fxy*(Rq*u(nxf,3,:)+Rs(1,:)+Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-1,1,:)+u(nxf-1,3,:))+fc*u(nxf,2,:))
    do j=3,nyf-2
      reshd(nxf-1,j,:)=rhs(nxf-1,j,:)-(fxx*(Rq*u(nxf,j,:)+Rs(1,:)+u(nxf-3,j,:))+fx*(u(nxf,j,:)+u(nxf-2,j,:)) &
        +fyy*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:))+fy*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:)) &
        +fxy*(u(nxf,j+1,:)+u(nxf,j-1,:)+u(nxf-2,j-1,:)+u(nxf-2,j+1,:))+fc*u(nxf-1,j,:))
      reshd(nxf,j,:)=rhs(nxf,j,:)-(fxx*(Rq*u(nxf-1,j,:)+Rs(2,:)+u(nxf-2,j,:))+fx*(Rq*u(nxf,j,:)+Rs(1,:)+u(nxf-1,j,:)) &
        +fyy*(u(nxf,j+2,:)+u(nxf,j-2,:))+fy*(u(nxf,j+1,:)+u(nxf,j-1,:)) &
        +fxy*(Rq*u(nxf,j+1,:)+Rs(1,:)+Rq*u(nxf,j-1,:)+Rs(1,:)+u(nxf-1,j-1,:)+u(nxf-1,j+1,:))+fc*u(nxf,j,:))
    enddo
!...NE corner
    reshd(nxf-1,nyf-1,:)=rhs(nxf-1,nyf-1,:)-(fxx*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-3,nyf-1,:))+fx*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)) &
      +fyy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:))+fy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:)) &
      +fxy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)+u(nxf-2,nyf,:))+fc*u(nxf-1,nyf-1,:))
    reshd(nxf,nyf-1,:)=rhs(nxf,nyf-1,:)-(fxx*(Rq*u(nxf-1,nyf-1,:)+Rs(2,:)+u(nxf-2,nyf-1,:))+fx*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-1,nyf-1,:)) &
      +fyy*(rbufN(nxf,1,:)+u(nxf,nyf-3,:))+fy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)) &
      +fxy*(Rq*u(nxf,nyf,:)+Rs(1,:)+Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-1,nyf-2,:)+u(nxf-1,nyf,:))+fc*u(nxf,nyf-1,:))
    reshd(nxf-1,nyf,:)=rhs(nxf-1,nyf,:)-(fxx*(Rq*u(nxf,nyf,:)+Rs(1,:)+u(nxf-3,nyf,:))+fx*(u(nxf,nyf,:)+u(nxf-2,nyf,:)) &
      +fyy*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:))+fy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:)) &
      +fxy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)+rbufN(nxf-2,1,:))+fc*u(nxf-1,nyf,:))
    reshd(nxf,nyf,:)=rhs(nxf,nyf,:)-(fxx*(Rq*u(nxf-1,nyf,:)+Rs(2,:)+u(nxf-2,nyf,:))+fx*(Rq*u(nxf,nyf,:)+Rs(1,:)+u(nxf-1,nyf,:)) &
      +fyy*(rbufN(nxf,2,:)+u(nxf,nyf-2,:))+fy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)) &
      +fxy*(Rq*rbufN(nxf,1,:)+Rs(1,:)+Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-1,nyf-1,:)+rbufN(nxf-1,1,:))+fc*u(nxf,nyf,:))
  else
    reshd(nxf-1,1,:)=rhs(nxf-1,1,:)-(fxx*(rbufE(1,1,:)+u(nxf-3,1,:))+fx*(u(nxf,1,:)+u(nxf-2,1,:)) &
      +fyy*(u(nxf-1,3,:)+rbufS(nxf-1,1,:))+fy*(u(nxf-1,2,:)+rbufS(nxf-1,2,:)) &
      +fxy*(u(nxf,2,:)+rbufS(nxf,2,:)+rbufS(nxf-2,2,:)+u(nxf-2,2,:))+fc*u(nxf-1,1,:))
    reshd(nxf,1,:)=rhs(nxf,1,:)-(fxx*(rbufE(2,1,:)+u(nxf-2,1,:))+fx*(rbufE(1,1,:)+u(nxf-1,1,:)) &
      +fyy*(u(nxf,3,:)+rbufS(nxf,1,:))+fy*(u(nxf,2,:)+rbufS(nxf,2,:)) &
      +fxy*(rbufE(1,2,:)+rbufSE(1,:)+rbufS(nxf-1,2,:)+u(nxf-1,2,:))+fc*u(nxf,1,:))
    reshd(nxf-1,2,:)=rhs(nxf-1,2,:)-(fxx*(rbufE(1,2,:)+u(nxf-3,2,:))+fx*(u(nxf,2,:)+u(nxf-2,2,:)) &
      +fyy*(u(nxf-1,4,:)+rbufS(nxf-1,2,:))+fy*(u(nxf-1,3,:)+u(nxf-1,1,:)) &
      +fxy*(u(nxf,3,:)+u(nxf,1,:)+u(nxf-2,1,:)+u(nxf-2,3,:))+fc*u(nxf-1,2,:))
    reshd(nxf,2,:)=rhs(nxf,2,:)-(fxx*(rbufE(2,2,:)+u(nxf-2,2,:))+fx*(rbufE(1,2,:)+u(nxf-1,2,:)) &
      +fyy*(u(nxf,4,:)+rbufS(nxf,2,:))+fy*(u(nxf,3,:)+u(nxf,1,:)) &
      +fxy*(rbufE(1,3,:)+rbufE(1,1,:)+u(nxf-1,1,:)+u(nxf-1,3,:))+fc*u(nxf,2,:))
    do j=3,nyf-2
      reshd(nxf-1,j,:)=rhs(nxf-1,j,:)-(fxx*(rbufE(1,j,:)+u(nxf-3,j,:))+fx*(u(nxf,j,:)+u(nxf-2,j,:)) &
        +fyy*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:))+fy*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:)) &
        +fxy*(u(nxf,j+1,:)+u(nxf,j-1,:)+u(nxf-2,j-1,:)+u(nxf-2,j+1,:))+fc*u(nxf-1,j,:))
      reshd(nxf,j,:)=rhs(nxf,j,:)-(fxx*(rbufE(2,j,:)+u(nxf-2,j,:))+fx*(rbufE(1,j,:)+u(nxf-1,j,:)) &
        +fyy*(u(nxf,j+2,:)+u(nxf,j-2,:))+fy*(u(nxf,j+1,:)+u(nxf,j-1,:)) &
        +fxy*(rbufE(1,j+1,:)+rbufE(1,j-1,:)+u(nxf-1,j-1,:)+u(nxf-1,j+1,:))+fc*u(nxf,j,:))
    enddo
    reshd(nxf-1,nyf-1,:)=rhs(nxf-1,nyf-1,:)-(fxx*(rbufE(1,nyf-1,:)+u(nxf-3,nyf-1,:))+fx*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)) &
      +fyy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:))+fy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:)) &
      +fxy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)+u(nxf-2,nyf,:))+fc*u(nxf-1,nyf-1,:))
    reshd(nxf,nyf-1,:)=rhs(nxf,nyf-1,:)-(fxx*(rbufE(2,nyf-1,:)+u(nxf-2,nyf-1,:))+fx*(rbufE(1,nyf-1,:)+u(nxf-1,nyf-1,:)) &
      +fyy*(rbufN(nxf,1,:)+u(nxf,nyf-3,:))+fy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)) &
      +fxy*(rbufE(1,nyf,:)+rbufE(1,nyf-2,:)+u(nxf-1,nyf-2,:)+u(nxf-1,nyf,:))+fc*u(nxf,nyf-1,:))
    reshd(nxf-1,nyf,:)=rhs(nxf-1,nyf,:)-(fxx*(rbufE(1,nyf,:)+u(nxf-3,nyf,:))+fx*(u(nxf,nyf,:)+u(nxf-2,nyf,:)) &
      +fyy*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:))+fy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:)) &
      +fxy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)+rbufN(nxf-2,1,:))+fc*u(nxf-1,nyf,:))
    reshd(nxf,nyf,:)=rhs(nxf,nyf,:)-(fxx*(rbufE(2,nyf,:)+u(nxf-2,nyf,:))+fx*(rbufE(1,nyf,:)+u(nxf-1,nyf,:)) &
      +fyy*(rbufN(nxf,2,:)+u(nxf,nyf-2,:))+fy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)) &
      +fxy*(rbufNE(1,:)+rbufE(1,nyf-1,:)+u(nxf-1,nyf-1,:)+rbufN(nxf-1,1,:))+fc*u(nxf,nyf,:))
  endif

#ELSE
!IF (HDop==3)

  fx6 = hDiff(1)*(rdx**6)
  fy6 = hDiff(2)*(rdy**6)
  fx4 = hDiff(1)*(rdx**4)*(rdy**2)
  fy4 = hDiff(2)*(rdx**2)*(rdy**4)
  fc  = 1.0_dp+4.0_dp*(5.0_dp*(fx6+fy6)+9.0_dp*(fx4+fy4))
  fx  = 15.0_dp*fx6+18.0_dp*fy4+24.0_dp*fx4
  fxx =-6.0_dp*(fx6+fx4)
  fy  = 15.0_dp*fy6+18.0_dp*fx4+24.0_dp*fy4
  fyy =-6.0_dp*(fy6+fy4)
  fxy =-12.0_dp*(fx4+fy4)
  fxxy= 3.0_dp*fx4
  fxyy= 3.0_dp*fy4

!.Interior points
  do j=4,nyf-3
    do i=4,nxf-3
      reshd(i,j,:)=rhs(i,j,:)-(fc*u(i,j,:) &
        -fx6*(u(i+3,j,:)+u(i-3,j,:))-fxx*(u(i+2,j,:)+u(i-2,j,:))-fx*(u(i+1,j,:)+u(i-1,j,:)) &
        -fy6*(u(i,j+3,:)+u(i,j-3,:))-fyy*(u(i,j+2,:)+u(i,j-2,:))-fy*(u(i,j+1,:)+u(i,j-1,:)) &
        -fxy*(u(i+1,j+1,:)+u(i+1,j-1,:)+u(i-1,j-1,:)+u(i-1,j+1,:)) &
        -fxxy*(u(i+2,j+1,:)+u(i+2,j-1,:)+u(i-2,j-1,:)+u(i-2,j+1,:)) &
        -fxyy*(u(i+1,j+2,:)+u(i+1,j-2,:)+u(i-1,j-2,:)+u(i-1,j+2,:)))
    enddo
  enddo

!.Boundary points
  call MPI_WAIT(req(5),stat(:,5),ierr)
  do i=4,nxf-3 ! Bottom boundary
    reshd(i,1,:)=rhs(i,1,:)-(fc*u(i,1,:) &
      -fx6*(u(i+3,1,:)+u(i-3,1,:))-fxx*(u(i+2,1,:)+u(i-2,1,:))-fx*(u(i+1,1,:)+u(i-1,1,:)) &
      -fy6*(u(i,4,:)+rbufS(i,1,:))-fyy*(u(i,3,:)+rbufS(i,2,:))-fy*(u(i,2,:)+rbufS(i,3,:)) &
      -fxy*(u(i+1,2,:)+rbufS(i+1,3,:)+rbufS(i-1,3,:)+u(i-1,2,:)) &
      -fxxy*(u(i+2,2,:)+rbufS(i+2,3,:)+rbufS(i-2,3,:)+u(i-2,2,:)) &
      -fxyy*(u(i+1,3,:)+rbufS(i+1,2,:)+rbufS(i-1,2,:)+u(i-1,3,:)))
    reshd(i,2,:)=rhs(i,2,:)-(fc*u(i,2,:) &
      -fx6*(u(i+3,2,:)+u(i-3,2,:))-fxx*(u(i+2,2,:)+u(i-2,2,:))-fx*(u(i+1,2,:)+u(i-1,2,:)) &
      -fy6*(u(i,5,:)+rbufS(i,2,:))-fyy*(u(i,4,:)+rbufS(i,3,:))-fy*(u(i,3,:)+u(i,1,:)) &
      -fxy*(u(i+1,3,:)+u(i+1,1,:)+u(i-1,1,:)+u(i-1,3,:)) &
      -fxxy*(u(i+2,3,:)+u(i+2,1,:)+u(i-2,1,:)+u(i-2,3,:)) &
      -fxyy*(u(i+1,4,:)+rbufS(i+1,3,:)+rbufS(i-1,3,:)+u(i-1,4,:)))
    reshd(i,3,:)=rhs(i,3,:)-(fc*u(i,3,:) &
      -fx6*(u(i+3,3,:)+u(i-3,3,:))-fxx*(u(i+2,3,:)+u(i-2,3,:))-fx*(u(i+1,3,:)+u(i-1,3,:)) &
      -fy6*(u(i,6,:)+rbufS(i,3,:))-fyy*(u(i,5,:)+u(i,1,:))-fy*(u(i,4,:)+u(i,2,:)) &
      -fxy*(u(i+1,4,:)+u(i+1,2,:)+u(i-1,2,:)+u(i-1,4,:)) &
      -fxxy*(u(i+2,4,:)+u(i+2,2,:)+u(i-2,2,:)+u(i-2,4,:)) &
      -fxyy*(u(i+1,5,:)+u(i+1,1,:)+u(i-1,1,:)+u(i-1,5,:)))
  enddo
  call MPI_WAIT(req(6),stat(:,6),ierr)
  do i=4,nxf-3 ! Top boundary
    reshd(i,nyf-2,:)=rhs(i,nyf-2,:)-(fc*u(i,nyf-2,:) &
      -fx6*(u(i+3,nyf-2,:)+u(i-3,nyf-2,:))-fxx*(u(i+2,nyf-2,:)+u(i-2,nyf-2,:))-fx*(u(i+1,nyf-2,:)+u(i-1,nyf-2,:)) &
      -fy6*(rbufN(i,1,:)+u(i,nyf-5,:))-fyy*(u(i,nyf,:)+u(i,nyf-4,:))-fy*(u(i,nyf-1,:)+u(i,nyf-3,:)) &
      -fxy*(u(i+1,nyf-1,:)+u(i+1,nyf-3,:)+u(i-1,nyf-3,:)+u(i-1,nyf-1,:)) &
      -fxxy*(u(i+2,nyf-1,:)+u(i+2,nyf-3,:)+u(i-2,nyf-3,:)+u(i-2,nyf-1,:)) &
      -fxyy*(u(i+1,nyf,:)+u(i+1,nyf-4,:)+u(i-1,nyf-4,:)+u(i-1,nyf,:)))
    reshd(i,nyf-1,:)=rhs(i,nyf-1,:)-(fc*u(i,nyf-1,:) &
      -fx6*(u(i+3,nyf-1,:)+u(i-3,nyf-1,:))-fxx*(u(i+2,nyf-1,:)+u(i-2,nyf-1,:))-fx*(u(i+1,nyf-1,:)+u(i-1,nyf-1,:)) &
      -fy6*(rbufN(i,2,:)+u(i,nyf-4,:))-fyy*(rbufN(i,1,:)+u(i,nyf-3,:))-fy*(u(i,nyf,:)+u(i,nyf-2,:)) &
      -fxy*(u(i+1,nyf,:)+u(i+1,nyf-2,:)+u(i-1,nyf-2,:)+u(i-1,nyf,:)) &
      -fxxy*(u(i+2,nyf,:)+u(i+2,nyf-2,:)+u(i-2,nyf-2,:)+u(i-2,nyf,:)) &
      -fxyy*(rbufN(i+1,1,:)+u(i+1,nyf-3,:)+u(i-1,nyf-3,:)+rbufN(i-1,1,:)))
    reshd(i,nyf,:)=rhs(i,nyf,:)-(fc*u(i,nyf,:) &
      -fx6*(u(i+3,nyf,:)+u(i-3,nyf,:))-fxx*(u(i+2,nyf,:)+u(i-2,nyf,:))-fx*(u(i+1,nyf,:)+u(i-1,nyf,:)) &
      -fy6*(rbufN(i,3,:)+u(i,nyf-3,:))-fyy*(rbufN(i,2,:)+u(i,nyf-2,:))-fy*(rbufN(i,1,:)+u(i,nyf-1,:)) &
      -fxy*(rbufN(i+1,1,:)+u(i+1,nyf-1,:)+u(i-1,nyf-1,:)+rbufN(i-1,1,:)) &
      -fxxy*(rbufN(i+2,1,:)+u(i+2,nyf-1,:)+u(i-2,nyf-1,:)+rbufN(i-2,1,:)) &
      -fxyy*(rbufN(i+1,2,:)+u(i+1,nyf-2,:)+u(i-1,nyf-2,:)+rbufN(i-1,2,:)))
  enddo
!.Left boundary
  call MPI_WAITALL(3,(/req(7),req(13),req(15)/),(/stat(:,7),stat(:,13),stat(:,15)/),ierr)
  if (iIDhd(qn)%a(ng) == 0) then !.Left most process (along x)
!...SW corner
    reshd(1,1,:)=rhs(1,1,:)-(fc*u(1,1,:) &
      -fx6*(u(4,1,:)+Lq*u(3,1,:)+Ls(3,:))-fxx*(u(3,1,:)+Lq*u(2,1,:)+Ls(2,:))-fx*(u(2,1,:)+Lq*u(1,1,:)+Ls(1,:)) &
      -fy6*(u(1,4,:)+rbufS(1,1,:))-fyy*(u(1,3,:)+rbufS(1,2,:))-fy*(u(1,2,:)+rbufS(1,3,:)) &
      -fxy*(u(2,2,:)+rbufS(2,3,:)+Lq*rbufS(1,3,:)+Ls(1,:)+Lq*u(1,2,:)+Ls(1,:)) &
      -fxxy*(u(3,2,:)+rbufS(3,3,:)+Lq*rbufS(2,3,:)+Ls(2,:)+Lq*u(2,2,:)+Ls(2,:)) &
      -fxyy*(u(2,3,:)+rbufS(2,2,:)+Lq*rbufS(1,2,:)+Ls(1,:)+Lq*u(1,3,:)+Ls(1,:)))
    reshd(2,1,:)=rhs(2,1,:)-(fc*u(2,1,:) &
      -fx6*(u(5,1,:)+Lq*u(2,1,:)+Ls(2,:))-fxx*(u(4,1,:)+Lq*u(1,1,:)+Ls(1,:))-fx*(u(3,1,:)+u(1,1,:)) &
      -fy6*(u(2,4,:)+rbufS(2,1,:))-fyy*(u(2,3,:)+rbufS(2,2,:))-fy*(u(2,2,:)+rbufS(2,3,:)) &
      -fxy*(u(3,2,:)+rbufS(3,3,:)+rbufS(1,3,:)+u(1,2,:)) &
      -fxxy*(u(4,2,:)+rbufS(4,3,:)+Lq*rbufS(1,3,:)+Ls(1,:)+Lq*u(1,2,:)+Ls(1,:)) &
      -fxyy*(u(3,3,:)+rbufS(3,2,:)+rbufS(1,2,:)+u(1,3,:)))
    reshd(3,1,:)=rhs(3,1,:)-(fc*u(3,1,:) &
      -fx6*(u(6,1,:)+Lq*u(1,1,:)+Ls(1,:))-fxx*(u(5,1,:)+u(1,1,:))-fx*(u(4,1,:)+u(2,1,:)) &
      -fy6*(u(3,4,:)+rbufS(3,1,:))-fyy*(u(3,3,:)+rbufS(3,2,:))-fy*(u(3,2,:)+rbufS(3,3,:)) &
      -fxy*(u(4,2,:)+rbufS(4,3,:)+rbufS(2,3,:)+u(2,2,:)) &
      -fxxy*(u(5,2,:)+rbufS(5,3,:)+rbufS(1,3,:)+u(1,2,:)) &
      -fxyy*(u(4,3,:)+rbufS(4,2,:)+rbufS(2,2,:)+u(2,3,:)))
    reshd(1,2,:)=rhs(1,2,:)-(fc*u(1,2,:) &
      -fx6*(u(4,2,:)+Lq*u(3,2,:)+Ls(3,:))-fxx*(u(3,2,:)+Lq*u(2,2,:)+Ls(2,:))-fx*(u(2,2,:)+Lq*u(1,2,:)+Ls(1,:)) &
      -fy6*(u(1,5,:)+rbufS(1,2,:))-fyy*(u(1,4,:)+rbufS(1,3,:))-fy*(u(1,3,:)+u(1,1,:)) &
      -fxy*(u(2,3,:)+u(2,1,:)+Lq*u(1,1,:)+Ls(1,:)+Lq*u(1,3,:)+Ls(1,:)) &
      -fxxy*(u(3,3,:)+u(3,1,:)+Lq*u(2,1,:)+Ls(2,:)+Lq*u(2,3,:)+Ls(2,:)) &
      -fxyy*(u(2,4,:)+rbufS(2,3,:)+Lq*rbufS(1,3,:)+Ls(1,:)+Lq*u(1,4,:)+Ls(1,:)))
    reshd(2,2,:)=rhs(2,2,:)-(fc*u(2,2,:) &
      -fx6*(u(5,2,:)+Lq*u(2,2,:)+Ls(2,:))-fxx*(u(4,2,:)+Lq*u(1,2,:)+Ls(1,:))-fx*(u(3,2,:)+u(1,2,:)) &
      -fy6*(u(2,5,:)+rbufS(2,2,:))-fyy*(u(2,4,:)+rbufS(2,3,:))-fy*(u(2,3,:)+u(2,1,:)) &
      -fxy*(u(3,3,:)+u(3,1,:)+u(1,1,:)+u(1,3,:)) &
      -fxxy*(u(4,3,:)+u(4,1,:)+Lq*u(1,1,:)+Ls(1,:)+Lq*u(1,3,:)+Ls(1,:)) &
      -fxyy*(u(3,4,:)+rbufS(3,3,:)+rbufS(1,3,:)+u(1,4,:)))
    reshd(3,2,:)=rhs(3,2,:)-(fc*u(3,2,:) &
      -fx6*(u(6,2,:)+Lq*u(1,2,:)+Ls(1,:))-fxx*(u(5,2,:)+u(1,2,:))-fx*(u(4,2,:)+u(2,2,:)) &
      -fy6*(u(3,5,:)+rbufS(3,2,:))-fyy*(u(3,4,:)+rbufS(3,3,:))-fy*(u(3,3,:)+u(3,1,:)) &
      -fxy*(u(4,3,:)+u(4,1,:)+u(2,1,:)+u(2,3,:)) &
      -fxxy*(u(5,3,:)+u(5,1,:)+u(1,1,:)+u(1,3,:)) &
      -fxyy*(u(4,4,:)+rbufS(4,3,:)+rbufS(2,3,:)+u(2,4,:)))
    reshd(1,3,:)=rhs(1,3,:)-(fc*u(1,3,:) &
      -fx6*(u(4,3,:)+Lq*u(3,3,:)+Ls(3,:))-fxx*(u(3,3,:)+Lq*u(2,3,:)+Ls(2,:))-fx*(u(2,3,:)+Lq*u(1,3,:)+Ls(1,:)) &
      -fy6*(u(1,6,:)+rbufS(1,3,:))-fyy*(u(1,5,:)+u(1,1,:))-fy*(u(1,4,:)+u(1,2,:)) &
      -fxy*(u(2,4,:)+u(2,2,:)+Lq*u(1,2,:)+Ls(1,:)+Lq*u(1,4,:)+Ls(1,:)) &
      -fxxy*(u(3,4,:)+u(3,2,:)+Lq*u(2,2,:)+Ls(2,:)+Lq*u(2,4,:)+Ls(2,:)) &
      -fxyy*(u(2,5,:)+u(2,1,:)+Lq*u(1,1,:)+Ls(1,:)+Lq*u(1,5,:)+Ls(1,:)))
    reshd(2,3,:)=rhs(2,3,:)-(fc*u(2,3,:) &
      -fx6*(u(5,3,:)+Lq*u(2,3,:)+Ls(2,:))-fxx*(u(4,3,:)+Lq*u(1,3,:)+Ls(1,:))-fx*(u(3,3,:)+u(1,3,:)) &
      -fy6*(u(2,6,:)+rbufS(2,3,:))-fyy*(u(2,5,:)+u(2,1,:))-fy*(u(2,4,:)+u(2,2,:)) &
      -fxy*(u(3,4,:)+u(3,2,:)+u(1,2,:)+u(1,4,:)) &
      -fxxy*(u(4,4,:)+u(4,2,:)+Lq*u(1,2,:)+Ls(1,:)+Lq*u(1,4,:)+Ls(1,:)) &
      -fxyy*(u(3,5,:)+u(3,1,:)+u(1,1,:)+u(1,5,:)))
    reshd(3,3,:)=rhs(3,3,:)-(fc*u(3,3,:) &
      -fx6*(u(6,3,:)+Lq*u(1,3,:)+Ls(1,:))-fxx*(u(5,3,:)+u(1,3,:))-fx*(u(4,3,:)+u(2,3,:)) &
      -fy6*(u(3,6,:)+rbufS(3,3,:))-fyy*(u(3,5,:)+u(3,1,:))-fy*(u(3,4,:)+u(3,2,:)) &
      -fxy*(u(4,4,:)+u(4,2,:)+u(2,2,:)+u(2,4,:)) &
      -fxxy*(u(5,4,:)+u(5,2,:)+u(1,2,:)+u(1,4,:)) &
      -fxyy*(u(4,5,:)+u(4,1,:)+u(2,1,:)+u(2,5,:)))
    do j=4,nyf-3 !.Left boundary
      reshd(1,j,:)=rhs(1,j,:)-(fc*u(1,j,:) &
        -fx6*(u(4,j,:)+Lq*u(3,j,:)+Ls(3,:))-fxx*(u(3,j,:)+Lq*u(2,j,:)+Ls(2,:))-fx*(u(2,j,:)+Lq*u(1,j,:)+Ls(1,:)) &
        -fy6*(u(1,j+3,:)+u(1,j-3,:))-fyy*(u(1,j+2,:)+u(1,j-2,:))-fy*(u(1,j+1,:)+u(1,j-1,:)) &
        -fxy*(u(2,j+1,:)+u(2,j-1,:)+Lq*u(1,j-1,:)+Ls(1,:)+Lq*u(1,j+1,:)+Ls(1,:)) &
        -fxxy*(u(3,j+1,:)+u(3,j-1,:)+Lq*u(2,j-1,:)+Ls(2,:)+Lq*u(2,j+1,:)+Ls(2,:)) &
        -fxyy*(u(2,j+2,:)+u(2,j-2,:)+Lq*u(1,j-2,:)+Ls(1,:)+Lq*u(1,j+2,:)+Ls(1,:)))
      reshd(2,j,:)=rhs(2,j,:)-(fc*u(2,j,:) &
        -fx6*(u(5,j,:)+Lq*u(2,j,:)+Ls(2,:))-fxx*(u(4,j,:)+Lq*u(1,j,:)+Ls(1,:))-fx*(u(3,j,:)+u(1,j,:)) &
        -fy6*(u(2,j+3,:)+u(2,j-3,:))-fyy*(u(2,j+2,:)+u(2,j-2,:))-fy*(u(2,j+1,:)+u(2,j-1,:)) &
        -fxy*(u(3,j+1,:)+u(3,j-1,:)+u(1,j-1,:)+u(1,j+1,:)) &
        -fxxy*(u(4,j+1,:)+u(4,j-1,:)+Lq*u(1,j-1,:)+Ls(1,:)+Lq*u(1,j+1,:)+Ls(1,:)) &
        -fxyy*(u(3,j+2,:)+u(3,j-2,:)+u(1,j-2,:)+u(1,j+2,:)))
      reshd(3,j,:)=rhs(3,j,:)-(fc*u(3,j,:) &
        -fx6*(u(6,j,:)+Lq*u(1,j,:)+Ls(1,:))-fxx*(u(5,j,:)+u(1,j,:))-fx*(u(4,j,:)+u(2,j,:)) &
        -fy6*(u(3,j+3,:)+u(3,j-3,:))-fyy*(u(3,j+2,:)+u(3,j-2,:))-fy*(u(3,j+1,:)+u(3,j-1,:)) &
        -fxy*(u(4,j+1,:)+u(4,j-1,:)+u(2,j-1,:)+u(2,j+1,:)) &
        -fxxy*(u(5,j+1,:)+u(5,j-1,:)+u(1,j-1,:)+u(1,j+1,:)) &
        -fxyy*(u(4,j+2,:)+u(4,j-2,:)+u(2,j-2,:)+u(2,j+2,:)))
    enddo
!...NW corner
    reshd(1,nyf-2,:)=rhs(1,nyf-2,:)-(fc*u(1,nyf-2,:) &
      -fx6*(u(4,nyf-2,:)+Lq*u(3,nyf-2,:)+Ls(3,:))-fxx*(u(3,nyf-2,:)+Lq*u(2,nyf-2,:)+Ls(2,:))-fx*(u(2,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:)) &
      -fy6*(rbufN(1,1,:)+u(1,nyf-5,:))-fyy*(u(1,nyf,:)+u(1,nyf-4,:))-fy*(u(1,nyf-1,:)+u(1,nyf-3,:)) &
      -fxy*(u(2,nyf-1,:)+u(2,nyf-3,:)+Lq*u(1,nyf-3,:)+Ls(1,:)+Lq*u(1,nyf-1,:)+Ls(1,:)) &
      -fxxy*(u(3,nyf-1,:)+u(3,nyf-3,:)+Lq*u(2,nyf-3,:)+Ls(2,:)+Lq*u(2,nyf-1,:)+Ls(2,:)) &
      -fxyy*(u(2,nyf,:)+u(2,nyf-4,:)+Lq*u(1,nyf-4,:)+Ls(1,:)+Lq*u(1,nyf,:)+Ls(1,:)))
    reshd(2,nyf-2,:)=rhs(2,nyf-2,:)-(fc*u(2,nyf-2,:) &
      -fx6*(u(5,nyf-2,:)+Lq*u(2,nyf-2,:)+Ls(2,:))-fxx*(u(4,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:))-fx*(u(3,nyf-2,:)+u(1,nyf-2,:)) &
      -fy6*(rbufN(2,1,:)+u(2,nyf-5,:))-fyy*(u(2,nyf,:)+u(2,nyf-4,:))-fy*(u(2,nyf-1,:)+u(2,nyf-3,:)) &
      -fxy*(u(3,nyf-1,:)+u(3,nyf-3,:)+u(1,nyf-3,:)+u(1,nyf-1,:)) &
      -fxxy*(u(4,nyf-1,:)+u(4,nyf-3,:)+Lq*u(1,nyf-3,:)+Ls(1,:)+Lq*u(1,nyf-1,:)+Ls(1,:)) &
      -fxyy*(u(3,nyf,:)+u(3,nyf-4,:)+u(1,nyf-4,:)+u(1,nyf,:)))
    reshd(3,nyf-2,:)=rhs(3,nyf-2,:)-(fc*u(3,nyf-2,:) &
      -fx6*(u(6,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:))-fxx*(u(5,nyf-2,:)+u(1,nyf-2,:))-fx*(u(4,nyf-2,:)+u(2,nyf-2,:)) &
      -fy6*(rbufN(3,1,:)+u(3,nyf-5,:))-fyy*(u(3,nyf,:)+u(3,nyf-4,:))-fy*(u(3,nyf-1,:)+u(3,nyf-3,:)) &
      -fxy*(u(4,nyf-1,:)+u(4,nyf-3,:)+u(2,nyf-3,:)+u(2,nyf-1,:)) &
      -fxxy*(u(5,nyf-1,:)+u(5,nyf-3,:)+u(1,nyf-3,:)+u(1,nyf-1,:)) &
      -fxyy*(u(4,nyf,:)+u(4,nyf-4,:)+u(2,nyf-4,:)+u(2,nyf,:)))
    reshd(1,nyf-1,:)=rhs(1,nyf-1,:)-(fc*u(1,nyf-1,:) &
      -fx6*(u(4,nyf-1,:)+Lq*u(3,nyf-1,:)+Ls(3,:))-fxx*(u(3,nyf-1,:)+Lq*u(2,nyf-1,:)+Ls(2,:))-fx*(u(2,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:)) &
      -fy6*(rbufN(1,2,:)+u(1,nyf-4,:))-fyy*(rbufN(1,1,:)+u(1,nyf-3,:))-fy*(u(1,nyf,:)+u(1,nyf-2,:)) &
      -fxy*(u(2,nyf,:)+u(2,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:)+Lq*u(1,nyf,:)+Ls(1,:)) &
      -fxxy*(u(3,nyf,:)+u(3,nyf-2,:)+Lq*u(2,nyf-2,:)+Ls(2,:)+Lq*u(2,nyf,:)+Ls(2,:)) &
      -fxyy*(rbufN(2,1,:)+u(2,nyf-3,:)+Lq*u(1,nyf-3,:)+Ls(1,:)+Lq*rbufN(1,1,:)+Ls(1,:)))
    reshd(2,nyf-1,:)=rhs(2,nyf-1,:)-(fc*u(2,nyf-1,:) &
      -fx6*(u(5,nyf-1,:)+Lq*u(2,nyf-1,:)+Ls(2,:))-fxx*(u(4,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:))-fx*(u(3,nyf-1,:)+u(1,nyf-1,:)) &
      -fy6*(rbufN(2,2,:)+u(2,nyf-4,:))-fyy*(rbufN(2,1,:)+u(2,nyf-3,:))-fy*(u(2,nyf,:)+u(2,nyf-2,:)) &
      -fxy*(u(3,nyf,:)+u(3,nyf-2,:)+u(1,nyf-2,:)+u(1,nyf,:)) &
      -fxxy*(u(4,nyf,:)+u(4,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:)+Lq*u(1,nyf,:)+Ls(1,:)) &
      -fxyy*(rbufN(3,1,:)+u(3,nyf-3,:)+u(1,nyf-3,:)+rbufN(1,1,:)))
    reshd(3,nyf-1,:)=rhs(3,nyf-1,:)-(fc*u(3,nyf-1,:) &
      -fx6*(u(6,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:))-fxx*(u(5,nyf-1,:)+u(1,nyf-1,:))-fx*(u(4,nyf-1,:)+u(2,nyf-1,:)) &
      -fy6*(rbufN(3,2,:)+u(3,nyf-4,:))-fyy*(rbufN(3,1,:)+u(3,nyf-3,:))-fy*(u(3,nyf,:)+u(3,nyf-2,:)) &
      -fxy*(u(4,nyf,:)+u(4,nyf-2,:)+u(2,nyf-2,:)+u(2,nyf,:)) &
      -fxxy*(u(5,nyf,:)+u(5,nyf-2,:)+u(1,nyf-2,:)+u(1,nyf,:)) &
      -fxyy*(rbufN(4,1,:)+u(4,nyf-3,:)+u(2,nyf-3,:)+rbufN(2,1,:)))
    reshd(1,nyf,:)=rhs(1,nyf,:)-(fc*u(1,nyf,:) &
      -fx6*(u(4,nyf,:)+Lq*u(3,nyf,:)+Ls(3,:))-fxx*(u(3,nyf,:)+Lq*u(2,nyf,:)+Ls(2,:))-fx*(u(2,nyf,:)+Lq*u(1,nyf,:)+Ls(1,:)) &
      -fy6*(rbufN(1,3,:)+u(1,nyf-3,:))-fyy*(rbufN(1,2,:)+u(1,nyf-2,:))-fy*(rbufN(1,1,:)+u(1,nyf-1,:)) &
      -fxy*(rbufN(2,1,:)+u(2,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:)+Lq*rbufN(1,1,:)+Ls(1,:)) &
      -fxxy*(rbufN(3,1,:)+u(3,nyf-1,:)+Lq*u(2,nyf-1,:)+Ls(2,:)+Lq*rbufN(2,1,:)+Ls(2,:)) &
      -fxyy*(rbufN(2,2,:)+u(2,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:)+Lq*rbufN(1,2,:)+Ls(1,:)))
    reshd(2,nyf,:)=rhs(2,nyf,:)-(fc*u(2,nyf,:) &
      -fx6*(u(5,nyf,:)+Lq*u(2,nyf,:)+Ls(2,:))-fxx*(u(4,nyf,:)+Lq*u(1,nyf,:)+Ls(1,:))-fx*(u(3,nyf,:)+u(1,nyf,:)) &
      -fy6*(rbufN(2,3,:)+u(2,nyf-3,:))-fyy*(rbufN(2,2,:)+u(2,nyf-2,:))-fy*(rbufN(2,1,:)+u(2,nyf-1,:)) &
      -fxy*(rbufN(3,1,:)+u(3,nyf-1,:)+u(1,nyf-1,:)+rbufN(1,1,:)) &
      -fxxy*(rbufN(4,1,:)+u(4,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:)+Lq*rbufN(1,1,:)+Ls(1,:)) &
      -fxyy*(rbufN(3,2,:)+u(3,nyf-2,:)+u(1,nyf-2,:)+rbufN(1,2,:)))
    reshd(3,nyf,:)=rhs(3,nyf,:)-(fc*u(3,nyf,:) &
      -fx6*(u(6,nyf,:)+Lq*u(1,nyf,:)+Ls(1,:))-fxx*(u(5,nyf,:)+u(1,nyf,:))-fx*(u(4,nyf,:)+u(2,nyf,:)) &
      -fy6*(rbufN(3,3,:)+u(3,nyf-3,:))-fyy*(rbufN(3,2,:)+u(3,nyf-2,:))-fy*(rbufN(3,1,:)+u(3,nyf-1,:)) &
      -fxy*(rbufN(4,1,:)+u(4,nyf-1,:)+u(2,nyf-1,:)+rbufN(2,1,:)) &
      -fxxy*(rbufN(5,1,:)+u(5,nyf-1,:)+u(1,nyf-1,:)+rbufN(1,1,:)) &
      -fxyy*(rbufN(4,2,:)+u(4,nyf-2,:)+u(2,nyf-2,:)+rbufN(2,2,:)))
  else
    reshd(1,1,:)=rhs(1,1,:)-(fc*u(1,1,:) &
      -fx6*(u(4,1,:)+rbufW(1,1,:))-fxx*(u(3,1,:)+rbufW(2,1,:))-fx*(u(2,1,:)+rbufW(3,1,:)) &
      -fy6*(u(1,4,:)+rbufS(1,1,:))-fyy*(u(1,3,:)+rbufS(1,2,:))-fy*(u(1,2,:)+rbufS(1,3,:)) &
      -fxy*(u(2,2,:)+rbufS(2,3,:)+rbufSW(1,:)+rbufW(3,2,:)) &
      -fxxy*(u(3,2,:)+rbufS(3,3,:)+rbufSW(3,:)+rbufW(2,2,:)) &
      -fxyy*(u(2,3,:)+rbufS(2,2,:)+rbufSW(2,:)+rbufW(3,3,:)))
    reshd(2,1,:)=rhs(2,1,:)-(fc*u(2,1,:) &
      -fx6*(u(5,1,:)+rbufW(2,1,:))-fxx*(u(4,1,:)+rbufW(3,1,:))-fx*(u(3,1,:)+u(1,1,:)) &
      -fy6*(u(2,4,:)+rbufS(2,1,:))-fyy*(u(2,3,:)+rbufS(2,2,:))-fy*(u(2,2,:)+rbufS(2,3,:)) &
      -fxy*(u(3,2,:)+rbufS(3,3,:)+rbufS(1,3,:)+u(1,2,:)) &
      -fxxy*(u(4,2,:)+rbufS(4,3,:)+rbufSW(1,:)+rbufW(3,2,:)) &
      -fxyy*(u(3,3,:)+rbufS(3,2,:)+rbufS(1,2,:)+u(1,3,:)))
    reshd(3,1,:)=rhs(3,1,:)-(fc*u(3,1,:) &
      -fx6*(u(6,1,:)+rbufW(3,1,:))-fxx*(u(5,1,:)+u(1,1,:))-fx*(u(4,1,:)+u(2,1,:)) &
      -fy6*(u(3,4,:)+rbufS(3,1,:))-fyy*(u(3,3,:)+rbufS(3,2,:))-fy*(u(3,2,:)+rbufS(3,3,:)) &
      -fxy*(u(4,2,:)+rbufS(4,3,:)+rbufS(2,3,:)+u(2,2,:)) &
      -fxxy*(u(5,2,:)+rbufS(5,3,:)+rbufS(1,3,:)+u(1,2,:)) &
      -fxyy*(u(4,3,:)+rbufS(4,2,:)+rbufS(2,2,:)+u(2,3,:)))
    reshd(1,2,:)=rhs(1,2,:)-(fc*u(1,2,:) &
      -fx6*(u(4,2,:)+rbufW(1,2,:))-fxx*(u(3,2,:)+rbufW(2,2,:))-fx*(u(2,2,:)+rbufW(3,2,:)) &
      -fy6*(u(1,5,:)+rbufS(1,2,:))-fyy*(u(1,4,:)+rbufS(1,3,:))-fy*(u(1,3,:)+u(1,1,:)) &
      -fxy*(u(2,3,:)+u(2,1,:)+rbufW(3,1,:)+rbufW(3,3,:)) &
      -fxxy*(u(3,3,:)+u(3,1,:)+rbufW(2,1,:)+rbufW(2,3,:)) &
      -fxyy*(u(2,4,:)+rbufS(2,3,:)+rbufSW(1,:)+rbufW(3,4,:)))
    reshd(2,2,:)=rhs(2,2,:)-(fc*u(2,2,:) &
      -fx6*(u(5,2,:)+rbufW(2,2,:))-fxx*(u(4,2,:)+rbufW(3,2,:))-fx*(u(3,2,:)+u(1,2,:)) &
      -fy6*(u(2,5,:)+rbufS(2,2,:))-fyy*(u(2,4,:)+rbufS(2,3,:))-fy*(u(2,3,:)+u(2,1,:)) &
      -fxy*(u(3,3,:)+u(3,1,:)+u(1,1,:)+u(1,3,:)) &
      -fxxy*(u(4,3,:)+u(4,1,:)+rbufW(3,1,:)+rbufW(3,3,:)) &
      -fxyy*(u(3,4,:)+rbufS(3,3,:)+rbufS(1,3,:)+u(1,4,:)))
    reshd(3,2,:)=rhs(3,2,:)-(fc*u(3,2,:) &
      -fx6*(u(6,2,:)+rbufW(3,2,:))-fxx*(u(5,2,:)+u(1,2,:))-fx*(u(4,2,:)+u(2,2,:)) &
      -fy6*(u(3,5,:)+rbufS(3,2,:))-fyy*(u(3,4,:)+rbufS(3,3,:))-fy*(u(3,3,:)+u(3,1,:)) &
      -fxy*(u(4,3,:)+u(4,1,:)+u(2,1,:)+u(2,3,:)) &
      -fxxy*(u(5,3,:)+u(5,1,:)+u(1,1,:)+u(1,3,:)) &
      -fxyy*(u(4,4,:)+rbufS(4,3,:)+rbufS(2,3,:)+u(2,4,:)))
    reshd(1,3,:)=rhs(1,3,:)-(fc*u(1,3,:) &
      -fx6*(u(4,3,:)+rbufW(1,3,:))-fxx*(u(3,3,:)+rbufW(2,3,:))-fx*(u(2,3,:)+rbufW(3,3,:)) &
      -fy6*(u(1,6,:)+rbufS(1,3,:))-fyy*(u(1,5,:)+u(1,1,:))-fy*(u(1,4,:)+u(1,2,:)) &
      -fxy*(u(2,4,:)+u(2,2,:)+rbufW(3,2,:)+rbufW(3,4,:)) &
      -fxxy*(u(3,4,:)+u(3,2,:)+rbufW(2,2,:)+rbufW(2,4,:)) &
      -fxyy*(u(2,5,:)+u(2,1,:)+rbufW(3,1,:)+rbufW(3,5,:)))
    reshd(2,3,:)=rhs(2,3,:)-(fc*u(2,3,:) &
      -fx6*(u(5,3,:)+rbufW(2,3,:))-fxx*(u(4,3,:)+rbufW(3,3,:))-fx*(u(3,3,:)+u(1,3,:)) &
      -fy6*(u(2,6,:)+rbufS(2,3,:))-fyy*(u(2,5,:)+u(2,1,:))-fy*(u(2,4,:)+u(2,2,:)) &
      -fxy*(u(3,4,:)+u(3,2,:)+u(1,2,:)+u(1,4,:)) &
      -fxxy*(u(4,4,:)+u(4,2,:)+rbufW(3,2,:)+rbufW(3,4,:)) &
      -fxyy*(u(3,5,:)+u(3,1,:)+u(1,1,:)+u(1,5,:)))
    reshd(3,3,:)=rhs(3,3,:)-(fc*u(3,3,:) &
      -fx6*(u(6,3,:)+rbufW(3,3,:))-fxx*(u(5,3,:)+u(1,3,:))-fx*(u(4,3,:)+u(2,3,:)) &
      -fy6*(u(3,6,:)+rbufS(3,3,:))-fyy*(u(3,5,:)+u(3,1,:))-fy*(u(3,4,:)+u(3,2,:)) &
      -fxy*(u(4,4,:)+u(4,2,:)+u(2,2,:)+u(2,4,:)) &
      -fxxy*(u(5,4,:)+u(5,2,:)+u(1,2,:)+u(1,4,:)) &
      -fxyy*(u(4,5,:)+u(4,1,:)+u(2,1,:)+u(2,5,:)))
    do j=4,nyf-3
      reshd(1,j,:)=rhs(1,j,:)-(fc*u(1,j,:) &
        -fx6*(u(4,j,:)+rbufW(1,j,:))-fxx*(u(3,j,:)+rbufW(2,j,:))-fx*(u(2,j,:)+rbufW(3,j,:)) &
        -fy6*(u(1,j+3,:)+u(1,j-3,:))-fyy*(u(1,j+2,:)+u(1,j-2,:))-fy*(u(1,j+1,:)+u(1,j-1,:)) &
        -fxy*(u(2,j+1,:)+u(2,j-1,:)+rbufW(3,j-1,:)+rbufW(3,j+1,:)) &
        -fxxy*(u(3,j+1,:)+u(3,j-1,:)+rbufW(2,j-1,:)+rbufW(2,j+1,:)) &
        -fxyy*(u(2,j+2,:)+u(2,j-2,:)+rbufW(3,j-2,:)+rbufW(3,j+2,:)))
      reshd(2,j,:)=rhs(2,j,:)-(fc*u(2,j,:) &
        -fx6*(u(5,j,:)+rbufW(2,j,:))-fxx*(u(4,j,:)+rbufW(3,j,:))-fx*(u(3,j,:)+u(1,j,:)) &
        -fy6*(u(2,j+3,:)+u(2,j-3,:))-fyy*(u(2,j+2,:)+u(2,j-2,:))-fy*(u(2,j+1,:)+u(2,j-1,:)) &
        -fxy*(u(3,j+1,:)+u(3,j-1,:)+u(1,j-1,:)+u(1,j+1,:)) &
        -fxxy*(u(4,j+1,:)+u(4,j-1,:)+rbufW(3,j-1,:)+rbufW(3,j+1,:)) &
        -fxyy*(u(3,j+2,:)+u(3,j-2,:)+u(1,j-2,:)+u(1,j+2,:)))
      reshd(3,j,:)=rhs(3,j,:)-(fc*u(3,j,:) &
        -fx6*(u(6,j,:)+rbufW(3,j,:))-fxx*(u(5,j,:)+u(1,j,:))-fx*(u(4,j,:)+u(2,j,:)) &
        -fy6*(u(3,j+3,:)+u(3,j-3,:))-fyy*(u(3,j+2,:)+u(3,j-2,:))-fy*(u(3,j+1,:)+u(3,j-1,:)) &
        -fxy*(u(4,j+1,:)+u(4,j-1,:)+u(2,j-1,:)+u(2,j+1,:)) &
        -fxxy*(u(5,j+1,:)+u(5,j-1,:)+u(1,j-1,:)+u(1,j+1,:)) &
        -fxyy*(u(4,j+2,:)+u(4,j-2,:)+u(2,j-2,:)+u(2,j+2,:)))
    enddo
    reshd(1,nyf-2,:)=rhs(1,nyf-2,:)-(fc*u(1,nyf-2,:) &
      -fx6*(u(4,nyf-2,:)+rbufW(1,nyf-2,:))-fxx*(u(3,nyf-2,:)+rbufW(2,nyf-2,:))-fx*(u(2,nyf-2,:)+rbufW(3,nyf-2,:)) &
      -fy6*(rbufN(1,1,:)+u(1,nyf-5,:))-fyy*(u(1,nyf,:)+u(1,nyf-4,:))-fy*(u(1,nyf-1,:)+u(1,nyf-3,:)) &
      -fxy*(u(2,nyf-1,:)+u(2,nyf-3,:)+rbufW(3,nyf-3,:)+rbufW(3,nyf-1,:)) &
      -fxxy*(u(3,nyf-1,:)+u(3,nyf-3,:)+rbufW(2,nyf-3,:)+rbufW(2,nyf-1,:)) &
      -fxyy*(u(2,nyf,:)+u(2,nyf-4,:)+rbufW(3,nyf-4,:)+rbufW(3,nyf,:)))
    reshd(2,nyf-2,:)=rhs(2,nyf-2,:)-(fc*u(2,nyf-2,:) &
      -fx6*(u(5,nyf-2,:)+rbufW(2,nyf-2,:))-fxx*(u(4,nyf-2,:)+rbufW(3,nyf-2,:))-fx*(u(3,nyf-2,:)+u(1,nyf-2,:)) &
      -fy6*(rbufN(2,1,:)+u(2,nyf-5,:))-fyy*(u(2,nyf,:)+u(2,nyf-4,:))-fy*(u(2,nyf-1,:)+u(2,nyf-3,:)) &
      -fxy*(u(3,nyf-1,:)+u(3,nyf-3,:)+u(1,nyf-3,:)+u(1,nyf-1,:)) &
      -fxxy*(u(4,nyf-1,:)+u(4,nyf-3,:)+rbufW(3,nyf-3,:)+rbufW(3,nyf-1,:)) &
      -fxyy*(u(3,nyf,:)+u(3,nyf-4,:)+u(1,nyf-4,:)+u(1,nyf,:)))
    reshd(3,nyf-2,:)=rhs(3,nyf-2,:)-(fc*u(3,nyf-2,:) &
      -fx6*(u(6,nyf-2,:)+rbufW(3,nyf-2,:))-fxx*(u(5,nyf-2,:)+u(1,nyf-2,:))-fx*(u(4,nyf-2,:)+u(2,nyf-2,:)) &
      -fy6*(rbufN(3,1,:)+u(3,nyf-5,:))-fyy*(u(3,nyf,:)+u(3,nyf-4,:))-fy*(u(3,nyf-1,:)+u(3,nyf-3,:)) &
      -fxy*(u(4,nyf-1,:)+u(4,nyf-3,:)+u(2,nyf-3,:)+u(2,nyf-1,:)) &
      -fxxy*(u(5,nyf-1,:)+u(5,nyf-3,:)+u(1,nyf-3,:)+u(1,nyf-1,:)) &
      -fxyy*(u(4,nyf,:)+u(4,nyf-4,:)+u(2,nyf-4,:)+u(2,nyf,:)))
    reshd(1,nyf-1,:)=rhs(1,nyf-1,:)-(fc*u(1,nyf-1,:) &
      -fx6*(u(4,nyf-1,:)+rbufW(1,nyf-1,:))-fxx*(u(3,nyf-1,:)+rbufW(2,nyf-1,:))-fx*(u(2,nyf-1,:)+rbufW(3,nyf-1,:)) &
      -fy6*(rbufN(1,2,:)+u(1,nyf-4,:))-fyy*(rbufN(1,1,:)+u(1,nyf-3,:))-fy*(u(1,nyf,:)+u(1,nyf-2,:)) &
      -fxy*(u(2,nyf,:)+u(2,nyf-2,:)+rbufW(3,nyf-2,:)+rbufW(3,nyf,:)) &
      -fxxy*(u(3,nyf,:)+u(3,nyf-2,:)+rbufW(2,nyf-2,:)+rbufW(2,nyf,:)) &
      -fxyy*(rbufN(2,1,:)+u(2,nyf-3,:)+rbufW(3,nyf-3,:)+rbufNW(1,:)))
    reshd(2,nyf-1,:)=rhs(2,nyf-1,:)-(fc*u(2,nyf-1,:) &
      -fx6*(u(5,nyf-1,:)+rbufW(2,nyf-1,:))-fxx*(u(4,nyf-1,:)+rbufW(3,nyf-1,:))-fx*(u(3,nyf-1,:)+u(1,nyf-1,:)) &
      -fy6*(rbufN(2,2,:)+u(2,nyf-4,:))-fyy*(rbufN(2,1,:)+u(2,nyf-3,:))-fy*(u(2,nyf,:)+u(2,nyf-2,:)) &
      -fxy*(u(3,nyf,:)+u(3,nyf-2,:)+u(1,nyf-2,:)+u(1,nyf,:)) &
      -fxxy*(u(4,nyf,:)+u(4,nyf-2,:)+rbufW(3,nyf-2,:)+rbufW(3,nyf,:)) &
      -fxyy*(rbufN(3,1,:)+u(3,nyf-3,:)+u(1,nyf-3,:)+rbufN(1,1,:)))
    reshd(3,nyf-1,:)=rhs(3,nyf-1,:)-(fc*u(3,nyf-1,:) &
      -fx6*(u(6,nyf-1,:)+rbufW(3,nyf-1,:))-fxx*(u(5,nyf-1,:)+u(1,nyf-1,:))-fx*(u(4,nyf-1,:)+u(2,nyf-1,:)) &
      -fy6*(rbufN(3,2,:)+u(3,nyf-4,:))-fyy*(rbufN(3,1,:)+u(3,nyf-3,:))-fy*(u(3,nyf,:)+u(3,nyf-2,:)) &
      -fxy*(u(4,nyf,:)+u(4,nyf-2,:)+u(2,nyf-2,:)+u(2,nyf,:)) &
      -fxxy*(u(5,nyf,:)+u(5,nyf-2,:)+u(1,nyf-2,:)+u(1,nyf,:)) &
      -fxyy*(rbufN(4,1,:)+u(4,nyf-3,:)+u(2,nyf-3,:)+rbufN(2,1,:)))
    reshd(1,nyf,:)=rhs(1,nyf,:)-(fc*u(1,nyf,:) &
      -fx6*(u(4,nyf,:)+rbufW(1,nyf,:))-fxx*(u(3,nyf,:)+rbufW(2,nyf,:))-fx*(u(2,nyf,:)+rbufW(3,nyf,:)) &
      -fy6*(rbufN(1,3,:)+u(1,nyf-3,:))-fyy*(rbufN(1,2,:)+u(1,nyf-2,:))-fy*(rbufN(1,1,:)+u(1,nyf-1,:)) &
      -fxy*(rbufN(2,1,:)+u(2,nyf-1,:)+rbufW(3,nyf-1,:)+rbufNW(1,:)) &
      -fxxy*(rbufN(3,1,:)+u(3,nyf-1,:)+rbufW(2,nyf-1,:)+rbufNW(3,:)) &
      -fxyy*(rbufN(2,2,:)+u(2,nyf-2,:)+rbufW(3,nyf-2,:)+rbufNW(2,:)))
    reshd(2,nyf,:)=rhs(2,nyf,:)-(fc*u(2,nyf,:) &
      -fx6*(u(5,nyf,:)+rbufW(2,nyf,:))-fxx*(u(4,nyf,:)+rbufW(3,nyf,:))-fx*(u(3,nyf,:)+u(1,nyf,:)) &
      -fy6*(rbufN(2,3,:)+u(2,nyf-3,:))-fyy*(rbufN(2,2,:)+u(2,nyf-2,:))-fy*(rbufN(2,1,:)+u(2,nyf-1,:)) &
      -fxy*(rbufN(3,1,:)+u(3,nyf-1,:)+u(1,nyf-1,:)+rbufN(1,1,:)) &
      -fxxy*(rbufN(4,1,:)+u(4,nyf-1,:)+rbufW(3,nyf-1,:)+rbufNW(1,:)) &
      -fxyy*(rbufN(3,2,:)+u(3,nyf-2,:)+u(1,nyf-2,:)+rbufN(1,2,:)))
    reshd(3,nyf,:)=rhs(3,nyf,:)-(fc*u(3,nyf,:) &
      -fx6*(u(6,nyf,:)+rbufW(3,nyf,:))-fxx*(u(5,nyf,:)+u(1,nyf,:))-fx*(u(4,nyf,:)+u(2,nyf,:)) &
      -fy6*(rbufN(3,3,:)+u(3,nyf-3,:))-fyy*(rbufN(3,2,:)+u(3,nyf-2,:))-fy*(rbufN(3,1,:)+u(3,nyf-1,:)) &
      -fxy*(rbufN(4,1,:)+u(4,nyf-1,:)+u(2,nyf-1,:)+rbufN(2,1,:)) &
      -fxxy*(rbufN(5,1,:)+u(5,nyf-1,:)+u(1,nyf-1,:)+rbufN(1,1,:)) &
      -fxyy*(rbufN(4,2,:)+u(4,nyf-2,:)+u(2,nyf-2,:)+rbufN(2,2,:)))
  endif
!.Right boundary
  call MPI_WAITALL(3,(/req(8),req(14),req(16)/),(/stat(:,8),stat(:,14),stat(:,16)/),ierr)
  if (iIDhd(qn)%a(ng) == iprocsm1) then !.Right most process (along x)
!...SE corner
    reshd(nxf-2,1,:)=rhs(nxf-2,1,:)-(fc*u(nxf-2,1,:) &
      -fx6*(Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-5,1,:))-fxx*(u(nxf,1,:)+u(nxf-4,1,:))-fx*(u(nxf-1,1,:)+u(nxf-3,1,:)) &
      -fy6*(u(nxf-2,4,:)+rbufS(nxf-2,1,:))-fyy*(u(nxf-2,3,:)+rbufS(nxf-2,2,:))-fy*(u(nxf-2,2,:)+rbufS(nxf-2,3,:)) &
      -fxy*(u(nxf-1,2,:)+rbufS(nxf-1,3,:)+rbufS(nxf-3,3,:)+u(nxf-3,2,:)) &
      -fxxy*(u(nxf,2,:)+rbufS(nxf,3,:)+rbufS(nxf-4,3,:)+u(nxf-4,2,:)) &
      -fxyy*(u(nxf-1,3,:)+rbufS(nxf-1,2,:)+rbufS(nxf-3,2,:)+u(nxf-3,3,:)))
    reshd(nxf-1,1,:)=rhs(nxf-1,1,:)-(fc*u(nxf-1,1,:) &
      -fx6*(Rq*u(nxf-1,1,:)+Rs(2,:)+u(nxf-4,1,:))-fxx*(Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-3,1,:))-fx*(u(nxf,1,:)+u(nxf-2,1,:)) &
      -fy6*(u(nxf-1,4,:)+rbufS(nxf-1,1,:))-fyy*(u(nxf-1,3,:)+rbufS(nxf-1,2,:))-fy*(u(nxf-1,2,:)+rbufS(nxf-1,3,:)) &
      -fxy*(u(nxf,2,:)+rbufS(nxf,3,:)+rbufS(nxf-2,3,:)+u(nxf-2,2,:)) &
      -fxxy*(Rq*u(nxf,2,:)+Rs(1,:)+Rq*rbufS(nxf,3,:)+Rs(1,:)+rbufS(nxf-3,3,:)+u(nxf-3,2,:)) &
      -fxyy*(u(nxf,3,:)+rbufS(nxf,2,:)+rbufS(nxf-2,2,:)+u(nxf-2,3,:)))
    reshd(nxf,1,:)=rhs(nxf,1,:)-(fc*u(nxf,1,:) &
      -fx6*(Rq*u(nxf-2,1,:)+Rs(3,:)+u(nxf-3,1,:))-fxx*(Rq*u(nxf-1,1,:)+Rs(2,:)+u(nxf-2,1,:))-fx*(Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-1,1,:)) &
      -fy6*(u(nxf,4,:)+rbufS(nxf,1,:))-fyy*(u(nxf,3,:)+rbufS(nxf,2,:))-fy*(u(nxf,2,:)+rbufS(nxf,3,:)) &
      -fxy*(Rq*u(nxf,2,:)+Rs(1,:)+Rq*rbufS(nxf,3,:)+Rs(1,:)+rbufS(nxf-1,3,:)+u(nxf-1,2,:)) &
      -fxxy*(Rq*u(nxf-1,2,:)+Rs(2,:)+Rq*rbufS(nxf-1,3,:)+Rs(2,:)+rbufS(nxf-2,3,:)+u(nxf-2,2,:)) &
      -fxyy*(Rq*u(nxf,3,:)+Rs(1,:)+Rq*rbufS(nxf,2,:)+Rs(1,:)+rbufS(nxf-1,2,:)+u(nxf-1,3,:)))
    reshd(nxf-2,2,:)=rhs(nxf-2,2,:)-(fc*u(nxf-2,2,:) &
      -fx6*(Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-5,2,:))-fxx*(u(nxf,2,:)+u(nxf-4,2,:))-fx*(u(nxf-1,2,:)+u(nxf-3,2,:)) &
      -fy6*(u(nxf-2,5,:)+rbufS(nxf-2,2,:))-fyy*(u(nxf-2,4,:)+rbufS(nxf-2,3,:))-fy*(u(nxf-2,3,:)+u(nxf-2,1,:)) &
      -fxy*(u(nxf-1,3,:)+u(nxf-1,1,:)+u(nxf-3,1,:)+u(nxf-3,3,:)) &
      -fxxy*(u(nxf,3,:)+u(nxf,1,:)+u(nxf-4,1,:)+u(nxf-4,3,:)) &
      -fxyy*(u(nxf-1,4,:)+rbufS(nxf-1,3,:)+rbufS(nxf-3,3,:)+u(nxf-3,4,:)))
    reshd(nxf-1,2,:)=rhs(nxf-1,2,:)-(fc*u(nxf-1,2,:) &
      -fx6*(Rq*u(nxf-1,2,:)+Rs(2,:)+u(nxf-4,2,:))-fxx*(Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-3,2,:))-fx*(u(nxf,2,:)+u(nxf-2,2,:)) &
      -fy6*(u(nxf-1,5,:)+rbufS(nxf-1,2,:))-fyy*(u(nxf-1,4,:)+rbufS(nxf-1,3,:))-fy*(u(nxf-1,3,:)+u(nxf-1,1,:)) &
      -fxy*(u(nxf,3,:)+u(nxf,1,:)+u(nxf-2,1,:)+u(nxf-2,3,:)) &
      -fxxy*(Rq*u(nxf,3,:)+Rs(1,:)+Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-3,1,:)+u(nxf-3,3,:)) &
      -fxyy*(u(nxf,4,:)+rbufS(nxf,3,:)+rbufS(nxf-2,3,:)+u(nxf-2,4,:)))
    reshd(nxf,2,:)=rhs(nxf,2,:)-(fc*u(nxf,2,:) &
      -fx6*(Rq*u(nxf-2,2,:)+Rs(3,:)+u(nxf-3,2,:))-fxx*(Rq*u(nxf-1,2,:)+Rs(2,:)+u(nxf-2,2,:))-fx*(Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-1,2,:)) &
      -fy6*(u(nxf,5,:)+rbufS(nxf,2,:))-fyy*(u(nxf,4,:)+rbufS(nxf,3,:))-fy*(u(nxf,3,:)+u(nxf,1,:)) &
      -fxy*(Rq*u(nxf,3,:)+Rs(1,:)+Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-1,1,:)+u(nxf-1,3,:)) &
      -fxxy*(Rq*u(nxf-1,3,:)+Rs(2,:)+Rq*u(nxf-1,1,:)+Rs(2,:)+u(nxf-2,1,:)+u(nxf-2,3,:)) &
      -fxyy*(Rq*u(nxf,4,:)+Rs(1,:)+Rq*rbufS(nxf,3,:)+Rs(1,:)+rbufS(nxf-1,3,:)+u(nxf-1,4,:)))
    reshd(nxf-2,3,:)=rhs(nxf-2,3,:)-(fc*u(nxf-2,3,:) &
      -fx6*(Rq*u(nxf,3,:)+Rs(1,:)+u(nxf-5,3,:))-fxx*(u(nxf,3,:)+u(nxf-4,3,:))-fx*(u(nxf-1,3,:)+u(nxf-3,3,:)) &
      -fy6*(u(nxf-2,6,:)+rbufS(nxf-2,3,:))-fyy*(u(nxf-2,5,:)+u(nxf-2,1,:))-fy*(u(nxf-2,4,:)+u(nxf-2,2,:)) &
      -fxy*(u(nxf-1,4,:)+u(nxf-1,2,:)+u(nxf-3,2,:)+u(nxf-3,4,:)) &
      -fxxy*(u(nxf,4,:)+u(nxf,2,:)+u(nxf-4,2,:)+u(nxf-4,4,:)) &
      -fxyy*(u(nxf-1,5,:)+u(nxf-1,1,:)+u(nxf-3,1,:)+u(nxf-3,5,:)))
    reshd(nxf-1,3,:)=rhs(nxf-1,3,:)-(fc*u(nxf-1,3,:) &
      -fx6*(Rq*u(nxf-1,3,:)+Rs(2,:)+u(nxf-4,3,:))-fxx*(Rq*u(nxf,3,:)+Rs(1,:)+u(nxf-3,3,:))-fx*(u(nxf,3,:)+u(nxf-2,3,:)) &
      -fy6*(u(nxf-1,6,:)+rbufS(nxf-1,3,:))-fyy*(u(nxf-1,5,:)+u(nxf-1,1,:))-fy*(u(nxf-1,4,:)+u(nxf-1,2,:)) &
      -fxy*(u(nxf,4,:)+u(nxf,2,:)+u(nxf-2,2,:)+u(nxf-2,4,:)) &
      -fxxy*(Rq*u(nxf,4,:)+Rs(1,:)+Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-3,2,:)+u(nxf-3,4,:)) &
      -fxyy*(u(nxf,5,:)+u(nxf,1,:)+u(nxf-2,1,:)+u(nxf-2,5,:)))
    reshd(nxf,3,:)=rhs(nxf,3,:)-(fc*u(nxf,3,:) &
      -fx6*(Rq*u(nxf-2,3,:)+Rs(3,:)+u(nxf-3,3,:))-fxx*(Rq*u(nxf-1,3,:)+Rs(2,:)+u(nxf-2,3,:))-fx*(Rq*u(nxf,3,:)+Rs(1,:)+u(nxf-1,3,:)) &
      -fy6*(u(nxf,6,:)+rbufS(nxf,3,:))-fyy*(u(nxf,5,:)+u(nxf,1,:))-fy*(u(nxf,4,:)+u(nxf,2,:)) &
      -fxy*(Rq*u(nxf,4,:)+Rs(1,:)+Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-1,2,:)+u(nxf-1,4,:)) &
      -fxxy*(Rq*u(nxf-1,4,:)+Rs(2,:)+Rq*u(nxf-1,2,:)+Rs(2,:)+u(nxf-2,2,:)+u(nxf-2,4,:)) &
      -fxyy*(Rq*u(nxf,5,:)+Rs(1,:)+Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-1,1,:)+u(nxf-1,5,:)))
    do j=4,nyf-3 !.Right boundary
      reshd(nxf-2,j,:)=rhs(nxf-2,j,:)-(fc*u(nxf-2,j,:) &
        -fx6*(Rq*u(nxf,j,:)+Rs(1,:)+u(nxf-5,j,:))-fxx*(u(nxf,j,:)+u(nxf-4,j,:))-fx*(u(nxf-1,j,:)+u(nxf-3,j,:)) &
        -fy6*(u(nxf-2,j+3,:)+u(nxf-2,j-3,:))-fyy*(u(nxf-2,j+2,:)+u(nxf-2,j-2,:))-fy*(u(nxf-2,j+1,:)+u(nxf-2,j-1,:)) &
        -fxy*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:)+u(nxf-3,j-1,:)+u(nxf-3,j+1,:)) &
        -fxxy*(u(nxf,j+1,:)+u(nxf,j-1,:)+u(nxf-4,j-1,:)+u(nxf-4,j+1,:)) &
        -fxyy*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:)+u(nxf-3,j-2,:)+u(nxf-3,j+2,:)))
      reshd(nxf-1,j,:)=rhs(nxf-1,j,:)-(fc*u(nxf-1,j,:) &
        -fx6*(Rq*u(nxf-1,j,:)+Rs(2,:)+u(nxf-4,j,:))-fxx*(Rq*u(nxf,j,:)+Rs(1,:)+u(nxf-3,j,:))-fx*(u(nxf,j,:)+u(nxf-2,j,:)) &
        -fy6*(u(nxf-1,j+3,:)+u(nxf-1,j-3,:))-fyy*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:))-fy*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:)) &
        -fxy*(u(nxf,j+1,:)+u(nxf,j-1,:)+u(nxf-2,j-1,:)+u(nxf-2,j+1,:)) &
        -fxxy*(Rq*u(nxf,j+1,:)+Rs(1,:)+Rq*u(nxf,j-1,:)+Rs(1,:)+u(nxf-3,j-1,:)+u(nxf-3,j+1,:)) &
        -fxyy*(u(nxf,j+2,:)+u(nxf,j-2,:)+u(nxf-2,j-2,:)+u(nxf-2,j+2,:)))
      reshd(nxf,j,:)=rhs(nxf,j,:)-(fc*u(nxf,j,:) &
        -fx6*(Rq*u(nxf-2,j,:)+Rs(3,:)+u(nxf-3,j,:))-fxx*(Rq*u(nxf-1,j,:)+Rs(2,:)+u(nxf-2,j,:))-fx*(Rq*u(nxf,j,:)+Rs(1,:)+u(nxf-1,j,:)) &
        -fy6*(u(nxf,j+3,:)+u(nxf,j-3,:))-fyy*(u(nxf,j+2,:)+u(nxf,j-2,:))-fy*(u(nxf,j+1,:)+u(nxf,j-1,:)) &
        -fxy*(Rq*u(nxf,j+1,:)+Rs(1,:)+Rq*u(nxf,j-1,:)+Rs(1,:)+u(nxf-1,j-1,:)+u(nxf-1,j+1,:)) &
        -fxxy*(Rq*u(nxf-1,j+1,:)+Rs(2,:)+Rq*u(nxf-1,j-1,:)+Rs(2,:)+u(nxf-2,j-1,:)+u(nxf-2,j+1,:)) &
        -fxyy*(Rq*u(nxf,j+2,:)+Rs(1,:)+Rq*u(nxf,j-2,:)+Rs(1,:)+u(nxf-1,j-2,:)+u(nxf-1,j+2,:)))
    enddo
!...NE corner
    reshd(nxf-2,nyf-2,:)=rhs(nxf-2,nyf-2,:)-(fc*u(nxf-2,nyf-2,:) &
      -fx6*(Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-5,nyf-2,:))-fxx*(u(nxf,nyf-2,:)+u(nxf-4,nyf-2,:))-fx*(u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:)) &
      -fy6*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-5,:))-fyy*(u(nxf-2,nyf,:)+u(nxf-2,nyf-4,:))-fy*(u(nxf-2,nyf-1,:)+u(nxf-2,nyf-3,:)) &
      -fxy*(u(nxf-1,nyf-1,:)+u(nxf-1,nyf-3,:)+u(nxf-3,nyf-3,:)+u(nxf-3,nyf-1,:)) &
      -fxxy*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:)+u(nxf-4,nyf-3,:)+u(nxf-4,nyf-1,:)) &
      -fxyy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-4,:)+u(nxf-3,nyf-4,:)+u(nxf-3,nyf,:)))
    reshd(nxf-1,nyf-2,:)=rhs(nxf-1,nyf-2,:)-(fc*u(nxf-1,nyf-2,:) &
      -fx6*(Rq*u(nxf-1,nyf-2,:)+Rs(2,:)+u(nxf-4,nyf-2,:))-fxx*(Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-3,nyf-2,:))-fx*(u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)) &
      -fy6*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-5,:))-fyy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-4,:))-fy*(u(nxf-1,nyf-1,:)+u(nxf-1,nyf-3,:)) &
      -fxy*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:)+u(nxf-2,nyf-3,:)+u(nxf-2,nyf-1,:)) &
      -fxxy*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+Rq*u(nxf,nyf-3,:)+Rs(1,:)+u(nxf-3,nyf-3,:)+u(nxf-3,nyf-1,:)) &
      -fxyy*(u(nxf,nyf,:)+u(nxf,nyf-4,:)+u(nxf-2,nyf-4,:)+u(nxf-2,nyf,:)))
    reshd(nxf,nyf-2,:)=rhs(nxf,nyf-2,:)-(fc*u(nxf,nyf-2,:) &
      -fx6*(Rq*u(nxf-2,nyf-2,:)+Rs(3,:)+u(nxf-3,nyf-2,:))-fxx*(Rq*u(nxf-1,nyf-2,:)+Rs(2,:)+u(nxf-2,nyf-2,:))-fx*(Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-1,nyf-2,:)) &
      -fy6*(rbufN(nxf,1,:)+u(nxf,nyf-5,:))-fyy*(u(nxf,nyf,:)+u(nxf,nyf-4,:))-fy*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:)) &
      -fxy*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+Rq*u(nxf,nyf-3,:)+Rs(1,:)+u(nxf-1,nyf-3,:)+u(nxf-1,nyf-1,:)) &
      -fxxy*(Rq*u(nxf-1,nyf-1,:)+Rs(2,:)+Rq*u(nxf-1,nyf-3,:)+Rs(2,:)+u(nxf-2,nyf-3,:)+u(nxf-2,nyf-1,:)) &
      -fxyy*(Rq*u(nxf,nyf,:)+Rs(1,:)+Rq*u(nxf,nyf-4,:)+Rs(1,:)+u(nxf-1,nyf-4,:)+u(nxf-1,nyf,:)))
    reshd(nxf-2,nyf-1,:)=rhs(nxf-2,nyf-1,:)-(fc*u(nxf-2,nyf-1,:) &
      -fx6*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-5,nyf-1,:))-fxx*(u(nxf,nyf-1,:)+u(nxf-4,nyf-1,:))-fx*(u(nxf-1,nyf-1,:)+u(nxf-3,nyf-1,:)) &
      -fy6*(rbufN(nxf-2,2,:)+u(nxf-2,nyf-4,:))-fyy*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-3,:))-fy*(u(nxf-2,nyf,:)+u(nxf-2,nyf-2,:)) &
      -fxy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:)+u(nxf-3,nyf,:)) &
      -fxxy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)+u(nxf-4,nyf-2,:)+u(nxf-4,nyf,:)) &
      -fxyy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:)+u(nxf-3,nyf-3,:)+rbufN(nxf-3,1,:)))
    reshd(nxf-1,nyf-1,:)=rhs(nxf-1,nyf-1,:)-(fc*u(nxf-1,nyf-1,:) &
      -fx6*(Rq*u(nxf-1,nyf-1,:)+Rs(2,:)+u(nxf-4,nyf-1,:))-fxx*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-3,nyf-1,:))-fx*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)) &
      -fy6*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-4,:))-fyy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:))-fy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:)) &
      -fxy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)+u(nxf-2,nyf,:)) &
      -fxxy*(Rq*u(nxf,nyf,:)+Rs(1,:)+Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-3,nyf-2,:)+u(nxf-3,nyf,:)) &
      -fxyy*(rbufN(nxf,1,:)+u(nxf,nyf-3,:)+u(nxf-2,nyf-3,:)+rbufN(nxf-2,1,:)))
    reshd(nxf,nyf-1,:)=rhs(nxf,nyf-1,:)-(fc*u(nxf,nyf-1,:) &
      -fx6*(Rq*u(nxf-2,nyf-1,:)+Rs(3,:)+u(nxf-3,nyf-1,:))-fxx*(Rq*u(nxf-1,nyf-1,:)+Rs(2,:)+u(nxf-2,nyf-1,:))-fx*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-1,nyf-1,:)) &
      -fy6*(rbufN(nxf,2,:)+u(nxf,nyf-4,:))-fyy*(rbufN(nxf,1,:)+u(nxf,nyf-3,:))-fy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)) &
      -fxy*(Rq*u(nxf,nyf,:)+Rs(1,:)+Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-1,nyf-2,:)+u(nxf-1,nyf,:)) &
      -fxxy*(Rq*u(nxf-1,nyf,:)+Rs(2,:)+Rq*u(nxf-1,nyf-2,:)+Rs(2,:)+u(nxf-2,nyf-2,:)+u(nxf-2,nyf,:)) &
      -fxyy*(Rq*rbufN(nxf,1,:)+Rs(1,:)+Rq*u(nxf,nyf-3,:)+Rs(1,:)+u(nxf-1,nyf-3,:)+rbufN(nxf-1,1,:)))
    reshd(nxf-2,nyf,:)=rhs(nxf-2,nyf,:)-(fc*u(nxf-2,nyf,:) &
      -fx6*(Rq*u(nxf,nyf,:)+Rs(1,:)+u(nxf-5,nyf,:))-fxx*(u(nxf,nyf,:)+u(nxf-4,nyf,:))-fx*(u(nxf-1,nyf,:)+u(nxf-3,nyf,:)) &
      -fy6*(rbufN(nxf-2,3,:)+u(nxf-2,nyf-3,:))-fyy*(rbufN(nxf-2,2,:)+u(nxf-2,nyf-2,:))-fy*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-1,:)) &
      -fxy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:)+u(nxf-3,nyf-1,:)+rbufN(nxf-3,1,:)) &
      -fxxy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)+u(nxf-4,nyf-1,:)+rbufN(nxf-4,1,:)) &
      -fxyy*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:)+rbufN(nxf-3,2,:)))
    reshd(nxf-1,nyf,:)=rhs(nxf-1,nyf,:)-(fc*u(nxf-1,nyf,:) &
      -fx6*(Rq*u(nxf-1,nyf,:)+Rs(2,:)+u(nxf-4,nyf,:))-fxx*(Rq*u(nxf,nyf,:)+Rs(1,:)+u(nxf-3,nyf,:))-fx*(u(nxf,nyf,:)+u(nxf-2,nyf,:)) &
      -fy6*(rbufN(nxf-1,3,:)+u(nxf-1,nyf-3,:))-fyy*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:))-fy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:)) &
      -fxy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)+rbufN(nxf-2,1,:)) &
      -fxxy*(Rq*rbufN(nxf,1,:)+Rs(1,:)+Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-3,nyf-1,:)+rbufN(nxf-3,1,:)) &
      -fxyy*(rbufN(nxf,2,:)+u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)+rbufN(nxf-2,2,:)))
    reshd(nxf,nyf,:)=rhs(nxf,nyf,:)-(fc*u(nxf,nyf,:) &
      -fx6*(Rq*u(nxf-2,nyf,:)+Rs(3,:)+u(nxf-3,nyf,:))-fxx*(Rq*u(nxf-1,nyf,:)+Rs(2,:)+u(nxf-2,nyf,:))-fx*(Rq*u(nxf,nyf,:)+Rs(1,:)+u(nxf-1,nyf,:)) &
      -fy6*(rbufN(nxf,3,:)+u(nxf,nyf-3,:))-fyy*(rbufN(nxf,2,:)+u(nxf,nyf-2,:))-fy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)) &
      -fxy*(Rq*rbufN(nxf,1,:)+Rs(1,:)+Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-1,nyf-1,:)+rbufN(nxf-1,1,:)) &
      -fxxy*(Rq*rbufN(nxf-1,1,:)+Rs(2,:)+Rq*u(nxf-1,nyf-1,:)+Rs(2,:)+u(nxf-2,nyf-1,:)+rbufN(nxf-2,1,:)) &
      -fxyy*(Rq*rbufN(nxf,2,:)+Rs(1,:)+Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-1,nyf-2,:)+rbufN(nxf-1,2,:)))
  else
    reshd(nxf-2,1,:)=rhs(nxf-2,1,:)-(fc*u(nxf-2,1,:) &
      -fx6*(rbufE(1,1,:)+u(nxf-5,1,:))-fxx*(u(nxf,1,:)+u(nxf-4,1,:))-fx*(u(nxf-1,1,:)+u(nxf-3,1,:)) &
      -fy6*(u(nxf-2,4,:)+rbufS(nxf-2,1,:))-fyy*(u(nxf-2,3,:)+rbufS(nxf-2,2,:))-fy*(u(nxf-2,2,:)+rbufS(nxf-2,3,:)) &
      -fxy*(u(nxf-1,2,:)+rbufS(nxf-1,3,:)+rbufS(nxf-3,3,:)+u(nxf-3,2,:)) &
      -fxxy*(u(nxf,2,:)+rbufS(nxf,3,:)+rbufS(nxf-4,3,:)+u(nxf-4,2,:)) &
      -fxyy*(u(nxf-1,3,:)+rbufS(nxf-1,2,:)+rbufS(nxf-3,2,:)+u(nxf-3,3,:)))
    reshd(nxf-1,1,:)=rhs(nxf-1,1,:)-(fc*u(nxf-1,1,:) &
      -fx6*(rbufE(2,1,:)+u(nxf-4,1,:))-fxx*(rbufE(1,1,:)+u(nxf-3,1,:))-fx*(u(nxf,1,:)+u(nxf-2,1,:)) &
      -fy6*(u(nxf-1,4,:)+rbufS(nxf-1,1,:))-fyy*(u(nxf-1,3,:)+rbufS(nxf-1,2,:))-fy*(u(nxf-1,2,:)+rbufS(nxf-1,3,:)) &
      -fxy*(u(nxf,2,:)+rbufS(nxf,3,:)+rbufS(nxf-2,3,:)+u(nxf-2,2,:)) &
      -fxxy*(rbufE(1,2,:)+rbufSE(1,:)+rbufS(nxf-3,3,:)+u(nxf-3,2,:)) &
      -fxyy*(u(nxf,3,:)+rbufS(nxf,2,:)+rbufS(nxf-2,2,:)+u(nxf-2,3,:)))
    reshd(nxf,1,:)=rhs(nxf,1,:)-(fc*u(nxf,1,:) &
      -fx6*(rbufE(3,1,:)+u(nxf-3,1,:))-fxx*(rbufE(2,1,:)+u(nxf-2,1,:))-fx*(rbufE(1,1,:)+u(nxf-1,1,:)) &
      -fy6*(u(nxf,4,:)+rbufS(nxf,1,:))-fyy*(u(nxf,3,:)+rbufS(nxf,2,:))-fy*(u(nxf,2,:)+rbufS(nxf,3,:)) &
      -fxy*(rbufE(1,2,:)+rbufSE(1,:)+rbufS(nxf-1,3,:)+u(nxf-1,2,:)) &
      -fxxy*(rbufE(2,2,:)+rbufSE(2,:)+rbufS(nxf-2,3,:)+u(nxf-2,2,:)) &
      -fxyy*(rbufE(1,3,:)+rbufSE(3,:)+rbufS(nxf-1,2,:)+u(nxf-1,3,:)))
    reshd(nxf-2,2,:)=rhs(nxf-2,2,:)-(fc*u(nxf-2,2,:) &
      -fx6*(rbufE(1,2,:)+u(nxf-5,2,:))-fxx*(u(nxf,2,:)+u(nxf-4,2,:))-fx*(u(nxf-1,2,:)+u(nxf-3,2,:)) &
      -fy6*(u(nxf-2,5,:)+rbufS(nxf-2,2,:))-fyy*(u(nxf-2,4,:)+rbufS(nxf-2,3,:))-fy*(u(nxf-2,3,:)+u(nxf-2,1,:)) &
      -fxy*(u(nxf-1,3,:)+u(nxf-1,1,:)+u(nxf-3,1,:)+u(nxf-3,3,:)) &
      -fxxy*(u(nxf,3,:)+u(nxf,1,:)+u(nxf-4,1,:)+u(nxf-4,3,:)) &
      -fxyy*(u(nxf-1,4,:)+rbufS(nxf-1,3,:)+rbufS(nxf-3,3,:)+u(nxf-3,4,:)))
    reshd(nxf-1,2,:)=rhs(nxf-1,2,:)-(fc*u(nxf-1,2,:) &
      -fx6*(rbufE(2,2,:)+u(nxf-4,2,:))-fxx*(rbufE(1,2,:)+u(nxf-3,2,:))-fx*(u(nxf,2,:)+u(nxf-2,2,:)) &
      -fy6*(u(nxf-1,5,:)+rbufS(nxf-1,2,:))-fyy*(u(nxf-1,4,:)+rbufS(nxf-1,3,:))-fy*(u(nxf-1,3,:)+u(nxf-1,1,:)) &
      -fxy*(u(nxf,3,:)+u(nxf,1,:)+u(nxf-2,1,:)+u(nxf-2,3,:)) &
      -fxxy*(rbufE(1,3,:)+rbufE(1,1,:)+u(nxf-3,1,:)+u(nxf-3,3,:)) &
      -fxyy*(u(nxf,4,:)+rbufS(nxf,3,:)+rbufS(nxf-2,3,:)+u(nxf-2,4,:)))
    reshd(nxf,2,:)=rhs(nxf,2,:)-(fc*u(nxf,2,:) &
      -fx6*(rbufE(3,2,:)+u(nxf-3,2,:))-fxx*(rbufE(2,2,:)+u(nxf-2,2,:))-fx*(rbufE(1,2,:)+u(nxf-1,2,:)) &
      -fy6*(u(nxf,5,:)+rbufS(nxf,2,:))-fyy*(u(nxf,4,:)+rbufS(nxf,3,:))-fy*(u(nxf,3,:)+u(nxf,1,:)) &
      -fxy*(rbufE(1,3,:)+rbufE(1,1,:)+u(nxf-1,1,:)+u(nxf-1,3,:)) &
      -fxxy*(rbufE(2,3,:)+rbufE(2,1,:)+u(nxf-2,1,:)+u(nxf-2,3,:)) &
      -fxyy*(rbufE(1,4,:)+rbufSE(1,:)+rbufS(nxf-1,3,:)+u(nxf-1,4,:)))
    reshd(nxf-2,3,:)=rhs(nxf-2,3,:)-(fc*u(nxf-2,3,:) &
      -fx6*(rbufE(1,3,:)+u(nxf-5,3,:))-fxx*(u(nxf,3,:)+u(nxf-4,3,:))-fx*(u(nxf-1,3,:)+u(nxf-3,3,:)) &
      -fy6*(u(nxf-2,6,:)+rbufS(nxf-2,3,:))-fyy*(u(nxf-2,5,:)+u(nxf-2,1,:))-fy*(u(nxf-2,4,:)+u(nxf-2,2,:)) &
      -fxy*(u(nxf-1,4,:)+u(nxf-1,2,:)+u(nxf-3,2,:)+u(nxf-3,4,:)) &
      -fxxy*(u(nxf,4,:)+u(nxf,2,:)+u(nxf-4,2,:)+u(nxf-4,4,:)) &
      -fxyy*(u(nxf-1,5,:)+u(nxf-1,1,:)+u(nxf-3,1,:)+u(nxf-3,5,:)))
    reshd(nxf-1,3,:)=rhs(nxf-1,3,:)-(fc*u(nxf-1,3,:) &
      -fx6*(rbufE(2,3,:)+u(nxf-4,3,:))-fxx*(rbufE(1,3,:)+u(nxf-3,3,:))-fx*(u(nxf,3,:)+u(nxf-2,3,:)) &
      -fy6*(u(nxf-1,6,:)+rbufS(nxf-1,3,:))-fyy*(u(nxf-1,5,:)+u(nxf-1,1,:))-fy*(u(nxf-1,4,:)+u(nxf-1,2,:)) &
      -fxy*(u(nxf,4,:)+u(nxf,2,:)+u(nxf-2,2,:)+u(nxf-2,4,:)) &
      -fxxy*(rbufE(1,4,:)+rbufE(1,2,:)+u(nxf-3,2,:)+u(nxf-3,4,:)) &
      -fxyy*(u(nxf,5,:)+u(nxf,1,:)+u(nxf-2,1,:)+u(nxf-2,5,:)))
    reshd(nxf,3,:)=rhs(nxf,3,:)-(fc*u(nxf,3,:) &
      -fx6*(rbufE(3,3,:)+u(nxf-3,3,:))-fxx*(rbufE(2,3,:)+u(nxf-2,3,:))-fx*(rbufE(1,3,:)+u(nxf-1,3,:)) &
      -fy6*(u(nxf,6,:)+rbufS(nxf,3,:))-fyy*(u(nxf,5,:)+u(nxf,1,:))-fy*(u(nxf,4,:)+u(nxf,2,:)) &
      -fxy*(rbufE(1,4,:)+rbufE(1,2,:)+u(nxf-1,2,:)+u(nxf-1,4,:)) &
      -fxxy*(rbufE(2,4,:)+rbufE(2,2,:)+u(nxf-2,2,:)+u(nxf-2,4,:)) &
      -fxyy*(rbufE(1,5,:)+rbufE(1,1,:)+u(nxf-1,1,:)+u(nxf-1,5,:)))
    do j=4,nyf-3
      reshd(nxf-2,j,:)=rhs(nxf-2,j,:)-(fc*u(nxf-2,j,:) &
        -fx6*(rbufE(1,j,:)+u(nxf-5,j,:))-fxx*(u(nxf,j,:)+u(nxf-4,j,:))-fx*(u(nxf-1,j,:)+u(nxf-3,j,:)) &
        -fy6*(u(nxf-2,j+3,:)+u(nxf-2,j-3,:))-fyy*(u(nxf-2,j+2,:)+u(nxf-2,j-2,:))-fy*(u(nxf-2,j+1,:)+u(nxf-2,j-1,:)) &
        -fxy*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:)+u(nxf-3,j-1,:)+u(nxf-3,j+1,:)) &
        -fxxy*(u(nxf,j+1,:)+u(nxf,j-1,:)+u(nxf-4,j-1,:)+u(nxf-4,j+1,:)) &
        -fxyy*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:)+u(nxf-3,j-2,:)+u(nxf-3,j+2,:)))
      reshd(nxf-1,j,:)=rhs(nxf-1,j,:)-(fc*u(nxf-1,j,:) &
        -fx6*(rbufE(2,j,:)+u(nxf-4,j,:))-fxx*(rbufE(1,j,:)+u(nxf-3,j,:))-fx*(u(nxf,j,:)+u(nxf-2,j,:)) &
        -fy6*(u(nxf-1,j+3,:)+u(nxf-1,j-3,:))-fyy*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:))-fy*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:)) &
        -fxy*(u(nxf,j+1,:)+u(nxf,j-1,:)+u(nxf-2,j-1,:)+u(nxf-2,j+1,:)) &
        -fxxy*(rbufE(1,j+1,:)+rbufE(1,j-1,:)+u(nxf-3,j-1,:)+u(nxf-3,j+1,:)) &
        -fxyy*(u(nxf,j+2,:)+u(nxf,j-2,:)+u(nxf-2,j-2,:)+u(nxf-2,j+2,:)))
      reshd(nxf,j,:)=rhs(nxf,j,:)-(fc*u(nxf,j,:) &
        -fx6*(rbufE(3,j,:)+u(nxf-3,j,:))-fxx*(rbufE(2,j,:)+u(nxf-2,j,:))-fx*(rbufE(1,j,:)+u(nxf-1,j,:)) &
        -fy6*(u(nxf,j+3,:)+u(nxf,j-3,:))-fyy*(u(nxf,j+2,:)+u(nxf,j-2,:))-fy*(u(nxf,j+1,:)+u(nxf,j-1,:)) &
        -fxy*(rbufE(1,j+1,:)+rbufE(1,j-1,:)+u(nxf-1,j-1,:)+u(nxf-1,j+1,:)) &
        -fxxy*(rbufE(2,j+1,:)+rbufE(2,j-1,:)+u(nxf-2,j-1,:)+u(nxf-2,j+1,:)) &
        -fxyy*(rbufE(1,j+2,:)+rbufE(1,j-2,:)+u(nxf-1,j-2,:)+u(nxf-1,j+2,:)))
    enddo
    reshd(nxf-2,nyf-2,:)=rhs(nxf-2,nyf-2,:)-(fc*u(nxf-2,nyf-2,:) &
      -fx6*(rbufE(1,nyf-2,:)+u(nxf-5,nyf-2,:))-fxx*(u(nxf,nyf-2,:)+u(nxf-4,nyf-2,:))-fx*(u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:)) &
      -fy6*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-5,:))-fyy*(u(nxf-2,nyf,:)+u(nxf-2,nyf-4,:))-fy*(u(nxf-2,nyf-1,:)+u(nxf-2,nyf-3,:)) &
      -fxy*(u(nxf-1,nyf-1,:)+u(nxf-1,nyf-3,:)+u(nxf-3,nyf-3,:)+u(nxf-3,nyf-1,:)) &
      -fxxy*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:)+u(nxf-4,nyf-3,:)+u(nxf-4,nyf-1,:)) &
      -fxyy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-4,:)+u(nxf-3,nyf-4,:)+u(nxf-3,nyf,:)))
    reshd(nxf-1,nyf-2,:)=rhs(nxf-1,nyf-2,:)-(fc*u(nxf-1,nyf-2,:) &
      -fx6*(rbufE(2,nyf-2,:)+u(nxf-4,nyf-2,:))-fxx*(rbufE(1,nyf-2,:)+u(nxf-3,nyf-2,:))-fx*(u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)) &
      -fy6*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-5,:))-fyy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-4,:))-fy*(u(nxf-1,nyf-1,:)+u(nxf-1,nyf-3,:)) &
      -fxy*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:)+u(nxf-2,nyf-3,:)+u(nxf-2,nyf-1,:)) &
      -fxxy*(rbufE(1,nyf-1,:)+rbufE(1,nyf-3,:)+u(nxf-3,nyf-3,:)+u(nxf-3,nyf-1,:)) &
      -fxyy*(u(nxf,nyf,:)+u(nxf,nyf-4,:)+u(nxf-2,nyf-4,:)+u(nxf-2,nyf,:)))
    reshd(nxf,nyf-2,:)=rhs(nxf,nyf-2,:)-(fc*u(nxf,nyf-2,:) &
      -fx6*(rbufE(3,nyf-2,:)+u(nxf-3,nyf-2,:))-fxx*(rbufE(2,nyf-2,:)+u(nxf-2,nyf-2,:))-fx*(rbufE(1,nyf-2,:)+u(nxf-1,nyf-2,:)) &
      -fy6*(rbufN(nxf,1,:)+u(nxf,nyf-5,:))-fyy*(u(nxf,nyf,:)+u(nxf,nyf-4,:))-fy*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:)) &
      -fxy*(rbufE(1,nyf-1,:)+rbufE(1,nyf-3,:)+u(nxf-1,nyf-3,:)+u(nxf-1,nyf-1,:)) &
      -fxxy*(rbufE(2,nyf-1,:)+rbufE(2,nyf-3,:)+u(nxf-2,nyf-3,:)+u(nxf-2,nyf-1,:)) &
      -fxyy*(rbufE(1,nyf,:)+rbufE(1,nyf-4,:)+u(nxf-1,nyf-4,:)+u(nxf-1,nyf,:)))
    reshd(nxf-2,nyf-1,:)=rhs(nxf-2,nyf-1,:)-(fc*u(nxf-2,nyf-1,:) &
      -fx6*(rbufE(1,nyf-1,:)+u(nxf-5,nyf-1,:))-fxx*(u(nxf,nyf-1,:)+u(nxf-4,nyf-1,:))-fx*(u(nxf-1,nyf-1,:)+u(nxf-3,nyf-1,:)) &
      -fy6*(rbufN(nxf-2,2,:)+u(nxf-2,nyf-4,:))-fyy*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-3,:))-fy*(u(nxf-2,nyf,:)+u(nxf-2,nyf-2,:)) &
      -fxy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:)+u(nxf-3,nyf,:)) &
      -fxxy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)+u(nxf-4,nyf-2,:)+u(nxf-4,nyf,:)) &
      -fxyy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:)+u(nxf-3,nyf-3,:)+rbufN(nxf-3,1,:)))
    reshd(nxf-1,nyf-1,:)=rhs(nxf-1,nyf-1,:)-(fc*u(nxf-1,nyf-1,:) &
      -fx6*(rbufE(2,nyf-1,:)+u(nxf-4,nyf-1,:))-fxx*(rbufE(1,nyf-1,:)+u(nxf-3,nyf-1,:))-fx*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)) &
      -fy6*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-4,:))-fyy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:))-fy*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:)) &
      -fxy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)+u(nxf-2,nyf,:)) &
      -fxxy*(rbufE(1,nyf,:)+rbufE(1,nyf-2,:)+u(nxf-3,nyf-2,:)+u(nxf-3,nyf,:)) &
      -fxyy*(rbufN(nxf,1,:)+u(nxf,nyf-3,:)+u(nxf-2,nyf-3,:)+rbufN(nxf-2,1,:)))
    reshd(nxf,nyf-1,:)=rhs(nxf,nyf-1,:)-(fc*u(nxf,nyf-1,:) &
      -fx6*(rbufE(3,nyf-1,:)+u(nxf-3,nyf-1,:))-fxx*(rbufE(2,nyf-1,:)+u(nxf-2,nyf-1,:))-fx*(rbufE(1,nyf-1,:)+u(nxf-1,nyf-1,:)) &
      -fy6*(rbufN(nxf,2,:)+u(nxf,nyf-4,:))-fyy*(rbufN(nxf,1,:)+u(nxf,nyf-3,:))-fy*(u(nxf,nyf,:)+u(nxf,nyf-2,:)) &
      -fxy*(rbufE(1,nyf,:)+rbufE(1,nyf-2,:)+u(nxf-1,nyf-2,:)+u(nxf-1,nyf,:)) &
      -fxxy*(rbufE(2,nyf,:)+rbufE(2,nyf-2,:)+u(nxf-2,nyf-2,:)+u(nxf-2,nyf,:)) &
      -fxyy*(rbufNE(1,:)+rbufE(1,nyf-3,:)+u(nxf-1,nyf-3,:)+rbufN(nxf-1,1,:)))
    reshd(nxf-2,nyf,:)=rhs(nxf-2,nyf,:)-(fc*u(nxf-2,nyf,:) &
      -fx6*(rbufE(1,nyf,:)+u(nxf-5,nyf,:))-fxx*(u(nxf,nyf,:)+u(nxf-4,nyf,:))-fx*(u(nxf-1,nyf,:)+u(nxf-3,nyf,:)) &
      -fy6*(rbufN(nxf-2,3,:)+u(nxf-2,nyf-3,:))-fyy*(rbufN(nxf-2,2,:)+u(nxf-2,nyf-2,:))-fy*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-1,:)) &
      -fxy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:)+u(nxf-3,nyf-1,:)+rbufN(nxf-3,1,:)) &
      -fxxy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)+u(nxf-4,nyf-1,:)+rbufN(nxf-4,1,:)) &
      -fxyy*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:)+rbufN(nxf-3,2,:)))
    reshd(nxf-1,nyf,:)=rhs(nxf-1,nyf,:)-(fc*u(nxf-1,nyf,:) &
      -fx6*(rbufE(2,nyf,:)+u(nxf-4,nyf,:))-fxx*(rbufE(1,nyf,:)+u(nxf-3,nyf,:))-fx*(u(nxf,nyf,:)+u(nxf-2,nyf,:)) &
      -fy6*(rbufN(nxf-1,3,:)+u(nxf-1,nyf-3,:))-fyy*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:))-fy*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:)) &
      -fxy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:)+rbufN(nxf-2,1,:)) &
      -fxxy*(rbufNE(1,:)+rbufE(1,nyf-1,:)+u(nxf-3,nyf-1,:)+rbufN(nxf-3,1,:)) &
      -fxyy*(rbufN(nxf,2,:)+u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:)+rbufN(nxf-2,2,:)))
    reshd(nxf,nyf,:)=rhs(nxf,nyf,:)-(fc*u(nxf,nyf,:) &
      -fx6*(rbufE(3,nyf,:)+u(nxf-3,nyf,:))-fxx*(rbufE(2,nyf,:)+u(nxf-2,nyf,:))-fx*(rbufE(1,nyf,:)+u(nxf-1,nyf,:)) &
      -fy6*(rbufN(nxf,3,:)+u(nxf,nyf-3,:))-fyy*(rbufN(nxf,2,:)+u(nxf,nyf-2,:))-fy*(rbufN(nxf,1,:)+u(nxf,nyf-1,:)) &
      -fxy*(rbufNE(1,:)+rbufE(1,nyf-1,:)+u(nxf-1,nyf-1,:)+rbufN(nxf-1,1,:)) &
      -fxxy*(rbufNE(3,:)+rbufE(2,nyf-1,:)+u(nxf-2,nyf-1,:)+rbufN(nxf-2,1,:)) &
      -fxyy*(rbufNE(2,:)+rbufE(1,nyf-2,:)+u(nxf-1,nyf-2,:)+rbufN(nxf-1,2,:)))
  endif

#ENDIF

  call MPI_WAITALL(8,(/req(1),req(2),req(3),req(4),req(9),req(10),req(11),req(12)/), &
    (/stat(:,1),stat(:,2),stat(:,3),stat(:,4),stat(:,9),stat(:,10),stat(:,11),stat(:,12)/),ierr)

#IF (MGTIMER==1)
    call CPU_TIME(rest2)
    resthd=resthd+rest2-rest1
#ENDIF
  end function reshd
!*********************************************************
#ELIF (HDop==8 .OR. HDop==12)
!*********************************************************
  subroutine relaxhd(ng,u,rhs,xBC,qn)
! Red-black-whie-orange Gauss-Seidel relaxation. The current value of the
! solution u is updated, using the right-hand-side function rhs.
  implicit none
  include 'mpif.h'
  real(DP), intent(inout) :: u(:,:,:)
  real(DP), intent(in) :: rhs(:,:,:)
  integer(I4B), intent(in) :: ng,xBC(:),qn
  integer(I4B) :: i,j,k,Nrank,Srank,Erank,Wrank
  integer(I4B) :: ierr,req(8),stat(MPI_STATUS_SIZE,8)
  real(DP), allocatable :: Ls(:,:),rs(:,:),loc2D(:,:)
  real(DP), allocatable, dimension(:,:,:) :: sbufN,sbufE,sbufS,sbufW, &
                                             rbufN,rbufE,rbufS,rbufW
  real(DP) :: fx,fy,delta,Lq,Rq
  real(DP), allocatable :: deltaL(:),deltaR(:)

  Nrank=nborhd(qn)%a(ng,1)
  Erank=nborhd(qn)%a(ng,3)
  Srank=nborhd(qn)%a(ng,5)
  Wrank=nborhd(qn)%a(ng,7)

!.Boundary conditions are implemented with some amount of flexibility.
!.The right boundary for example, if u(nxf+1) is needed, it will be
!.written as
!.  u(nxf+1) = Rq*u(nxf) + Rr*u(nxf-1) + Rs
!.so that symmetry, extrapolation or Dirichlet conditions can be
!.easily imposed.
!......... BCs of u ................................!
  allocate(Ls(HDopd4,nzL),rs(HDopd4,nzL),loc2D(HDopd4,nzL))
!.LHS BC of u
  if (xBC(1)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Lq=-1.0_dp
    forall(i=1:HDopd4,k=1:nzL) loc2D(i,k)=sum(u(i,:,k))
    call MPI_ALLreduce(loc2D,ls,HDopd4*nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    ls=2.0_dp*ls/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC==1 .OR. xBC(1)==-1) then ! symmetry
    Lq=dble(xBC(1))
    ls=0.0_dp
  endif
!.RHS BC of u
  if (xBC(2)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Rq=-1.0_dp
    forall(i=1:HDopd4,k=1:nzL) loc2D(i,k)=sum(u(nxf-i+1,:,k))
    call MPI_ALLreduce(loc2D,rs,HDopd4*nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    rs=2.0_dp*rs/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC(2)==1 .OR. xBC(2)==-1) then ! symmetry
    Rq=dble(xBC(2))
    rs=0.0_dp
  endif
!........................................................!
  fx=hDiff(1)*((dble(nxf*iprocshd(qn)%a(ng))/bLx)**HDopd2)
  fy=hDiff(2)*((dble(nyf*jprocshd(qn)%a(ng))/bLy)**HDopd2)

  allocate(sbufN(nxf,HDopd4,nzL),sbufE(HDopd4,nyf,nzL),sbufW(HDopd4,nyf,nzL),sbufS(nxf,HDopd4,nzL))
  allocate(rbufN(nxf,HDopd4,nzL),rbufE(HDopd4,nyf,nzL),rbufW(HDopd4,nyf,nzL),rbufS(nxf,HDopd4,nzL))
!.Send top 2 rows to North rank, receive rows from South rank
  sbufN=u(:,nyf-HDopd4+1:nyf,:)
  call MPI_ISEND(sbufN,nxf*HDopd4*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(1),ierr)
!.Send right 2 columns to East rank, receive columns from West rank
  sbufE=u(nxf-HDopd4+1:nxf,:,:)
  call MPI_ISEND(sbufE,HDopd4*nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(2),ierr)
!.Send left 2 columns to West rank, receive columns from East rank
  sbufW=u(1:HDopd4,:,:)
  call MPI_ISEND(sbufW,HDopd4*nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(3),ierr)
!.Send bottom 2 rows to South rank, receive rows from North rank
  sbufS=u(:,1:HDopd4,:)
  call MPI_ISEND(sbufS,nxf*HDopd4*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(4),ierr)
!.Receive rows from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxf*HDopd4*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(5),ierr)
!.Receive columns from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,HDopd4*nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(6),ierr)
!.Receive columns from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,HDopd4*nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(7),ierr)
!.Receive rows from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxf*HDopd4*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(8),ierr)
  
  call MPI_WAITALL(2,(/req(5),req(6)/),(/stat(:,5),stat(:,6)/),ierr)
#IF (HDop==8)
  delta=omehd(qn)/(1.0_dp+6.0_dp*(fx+fy))
  allocate(deltaL(1),deltaR(1))
  deltaL=omehd(qn)/(1.0_dp+6.0_dp*(fx+fy)-4.0_dp*fx*Lq)
  deltaR=omehd(qn)/(1.0_dp+6.0_dp*(fx+fy)-4.0_dp*fx*Rq)

!........................... FIRST ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,1,:)=Oomehd(qn)*u(1,1,:)+(rhs(1,1,:) &
!    u(1,1,:)=(rhs(1,1,:) &
      -fx*(u(3,1,:)-4.0_dp*(u(2,1,:)+Ls(1,:))+Lq*u(2,1,:)+Ls(2,:)) &
      -fy*(u(1,3,:)-4.0_dp*(u(1,2,:)+rbufS(1,2,:))+rbufS(1,1,:)))*deltaL
    u(2,1,:)=Oomehd(qn)*u(2,1,:)+(rhs(2,1,:) &
!    u(2,1,:)=(rhs(2,1,:) &
      -fx*(u(4,1,:)-4.0_dp*(u(3,1,:)+u(1,1,:))+Lq*u(1,1,:)+Ls(1,:)) &
      -fy*(u(2,3,:)-4.0_dp*(u(2,2,:)+rbufS(2,2,:))+rbufS(2,1,:)))*delta
  else
    u(1,1,:)=Oomehd(qn)*u(1,1,:)+(rhs(1,1,:) &
!    u(1,1,:)=(rhs(1,1,:) &
      -fx*(u(3,1,:)-4.0_dp*(u(2,1,:)+rbufW(2,1,:))+rbufW(1,1,:)) &
      -fy*(u(1,3,:)-4.0_dp*(u(1,2,:)+rbufS(1,2,:))+rbufS(1,1,:)))*delta
    u(2,1,:)=Oomehd(qn)*u(2,1,:)+(rhs(2,1,:) &
!    u(2,1,:)=(rhs(2,1,:) &
      -fx*(u(4,1,:)-4.0_dp*(u(3,1,:)+u(1,1,:))+rbufW(2,1,:)) &
      -fy*(u(2,3,:)-4.0_dp*(u(2,2,:)+rbufS(2,2,:))+rbufS(2,1,:)))*delta
  endif
  do i=3,nxf-2
    u(i,1,:)=Oomehd(qn)*u(i,1,:)+(rhs(i,1,:) &
!    u(i,1,:)=(rhs(i,1,:) &
      -fx*(u(i+2,1,:)-4.0_dp*(u(i+1,1,:)+u(i-1,1,:))+u(i-2,1,:)) &
      -fy*(u(i,3,:)-4.0_dp*(u(i,2,:)+rbufS(i,2,:))+rbufS(i,1,:)))*delta
  enddo
  call MPI_WAIT(req(7),stat(:,7),ierr)
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-1,1,:)=Oomehd(qn)*u(nxf-1,1,:)+(rhs(nxf-1,1,:) &
!    u(nxf-1,1,:)=(rhs(nxf-1,1,:) &
      -fx*(Rq*u(nxf,1,:)+Rs(1,:)-4.0_dp*(u(nxf,1,:)+u(nxf-2,1,:))+u(nxf-3,1,:)) &
      -fy*(u(nxf-1,3,:)-4.0_dp*(u(nxf-1,2,:)+rbufS(nxf-1,2,:))+rbufS(nxf-1,1,:)))*delta
    u(nxf,1,:)=Oomehd(qn)*u(nxf,1,:)+(rhs(nxf,1,:) &
!    u(nxf,1,:)=(rhs(nxf,1,:) &
      -fx*(Rq*u(nxf-1,1,:)+Rs(2,:)-4.0_dp*(rs(1,:)+u(nxf-1,1,:))+u(nxf-2,1,:)) &
      -fy*(u(nxf,3,:)-4.0_dp*(u(nxf,2,:)+rbufS(nxf,2,:))+rbufS(nxf,1,:)))*deltaR
  else
    u(nxf-1,1,:)=Oomehd(qn)*u(nxf-1,1,:)+(rhs(nxf-1,1,:) &
!    u(nxf-1,1,:)=(rhs(nxf-1,1,:) &
      -fx*(rbufE(1,1,:)-4.0_dp*(u(nxf,1,:)+u(nxf-2,1,:))+u(nxf-3,1,:)) &
      -fy*(u(nxf-1,3,:)-4.0_dp*(u(nxf-1,2,:)+rbufS(nxf-1,2,:))+rbufS(nxf-1,1,:)))*delta
    u(nxf,1,:)=Oomehd(qn)*u(nxf,1,:)+(rhs(nxf,1,:) &
!    u(nxf,1,:)=(rhs(nxf,1,:) &
      -fx*(rbufE(2,1,:)-4.0_dp*(rbufE(1,1,:)+u(nxf-1,1,:))+u(nxf-2,1,:)) &
      -fy*(u(nxf,3,:)-4.0_dp*(u(nxf,2,:)+rbufS(nxf,2,:))+rbufS(nxf,1,:)))*delta
  endif
!........................... SECOND ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,2,:)=Oomehd(qn)*u(1,2,:)+(rhs(1,2,:) &
!    u(1,2,:)=(rhs(1,2,:) &
      -fx*(u(3,2,:)-4.0_dp*(u(2,2,:)+Ls(1,:))+Lq*u(2,2,:)+Ls(2,:)) &
      -fy*(u(1,4,:)-4.0_dp*(u(1,3,:)+u(1,1,:))+rbufS(1,2,:)))*deltaL
    u(2,2,:)=Oomehd(qn)*u(2,2,:)+(rhs(2,2,:) &
!    u(2,2,:)=(rhs(2,2,:) &
      -fx*(u(4,2,:)-4.0_dp*(u(3,2,:)+u(1,2,:))+Lq*u(1,2,:)+Ls(1,:)) &
      -fy*(u(2,4,:)-4.0_dp*(u(2,3,:)+u(2,1,:))+rbufS(2,2,:)))*delta
  else
    u(1,2,:)=Oomehd(qn)*u(1,2,:)+(rhs(1,2,:) &
!    u(1,2,:)=(rhs(1,2,:) &
      -fx*(u(3,2,:)-4.0_dp*(u(2,2,:)+rbufW(2,2,:))+rbufW(1,2,:)) &
      -fy*(u(1,4,:)-4.0_dp*(u(1,3,:)+u(1,1,:))+rbufS(1,2,:)))*delta
    u(2,2,:)=Oomehd(qn)*u(2,2,:)+(rhs(2,2,:) &
!    u(2,2,:)=(rhs(2,2,:) &
      -fx*(u(4,2,:)-4.0_dp*(u(3,2,:)+u(1,2,:))+rbufW(2,2,:)) &
      -fy*(u(2,4,:)-4.0_dp*(u(2,3,:)+u(2,1,:))+rbufS(2,2,:)))*delta
  endif
  do i=3,nxf-2
    u(i,2,:)=Oomehd(qn)*u(i,2,:)+(rhs(i,2,:) &
!    u(i,2,:)=(rhs(i,2,:) &
    -fx*(u(i+2,2,:)-4.0_dp*(u(i+1,2,:)+u(i-1,2,:))+u(i-2,2,:)) &
    -fy*(u(i,4,:)-4.0_dp*(u(i,3,:)+u(i,1,:))+rbufS(i,2,:)))*delta
  enddo
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-1,2,:)=Oomehd(qn)*u(nxf-1,2,:)+(rhs(nxf-1,2,:) &
!    u(nxf-1,2,:)=(rhs(nxf-1,2,:) &
      -fx*(Rq*u(nxf,2,:)+Rs(1,:)-4.0_dp*(u(nxf,2,:)+u(nxf-2,2,:))+u(nxf-3,2,:)) &
      -fy*(u(nxf-1,4,:)-4.0_dp*(u(nxf-1,3,:)+u(nxf-1,1,:))+rbufS(nxf-1,2,:)))*delta
    u(nxf,2,:)=Oomehd(qn)*u(nxf,2,:)+(rhs(nxf,2,:) &
!    u(nxf,2,:)=(rhs(nxf,2,:) &
      -fx*(Rq*u(nxf-1,2,:)+Rs(2,:)-4.0_dp*(rs(1,:)+u(nxf-1,2,:))+u(nxf-2,2,:)) &
      -fy*(u(nxf,4,:)-4.0_dp*(u(nxf,3,:)+u(nxf,1,:))+rbufS(nxf,2,:)))*deltaR
  else
    u(nxf-1,2,:)=Oomehd(qn)*u(nxf-1,2,:)+(rhs(nxf-1,2,:) &
!    u(nxf-1,2,:)=(rhs(nxf-1,2,:) &
      -fx*(rbufE(1,2,:)-4.0_dp*(u(nxf,2,:)+u(nxf-2,2,:))+u(nxf-3,2,:)) &
      -fy*(u(nxf-1,4,:)-4.0_dp*(u(nxf-1,3,:)+u(nxf-1,1,:))+rbufS(nxf-1,2,:)))*delta
    u(nxf,2,:)=Oomehd(qn)*u(nxf,2,:)+(rhs(nxf,2,:) &
!    u(nxf,2,:)=(rhs(nxf,2,:) &
      -fx*(rbufE(2,2,:)-4.0_dp*(rbufE(1,2,:)+u(nxf-1,2,:))+u(nxf-2,2,:)) &
      -fy*(u(nxf,4,:)-4.0_dp*(u(nxf,3,:)+u(nxf,1,:))+rbufS(nxf,2,:)))*delta
  endif
!................... AWAY FROM BOTTOM AND TOP BOUNDARIES ...................!
  do j=3,nyf-2
!...Left boundary
    if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
     u(1,j,:)=Oomehd(qn)*u(1,j,:)+(rhs(1,j,:) &
!      u(1,j,:)=(rhs(1,j,:) &
        -fx*(u(3,j,:)-4.0_dp*(u(2,j,:)+Ls(1,:))+Lq*u(2,j,:)+Ls(2,:)) &
        -fy*(u(1,j+2,:)-4.0_dp*(u(1,j+1,:)+u(1,j-1,:))+u(1,j-2,:)))*deltaL
      u(2,j,:)=Oomehd(qn)*u(2,j,:)+(rhs(2,j,:) &
!      u(2,j,:)=(rhs(2,j,:) &
        -fx*(u(4,j,:)-4.0_dp*(u(3,j,:)+u(1,j,:))+Lq*u(1,j,:)+Ls(1,:)) &
        -fy*(u(2,j+2,:)-4.0_dp*(u(2,j+1,:)+u(2,j-1,:))+u(2,j-2,:)))*delta
    else
      u(1,j,:)=Oomehd(qn)*u(1,j,:)+(rhs(1,j,:) &
!      u(1,j,:)=(rhs(1,j,:) &
        -fx*(u(3,j,:)-4.0_dp*(u(2,j,:)+rbufW(2,j,:))+rbufW(1,j,:)) &
        -fy*(u(1,j+2,:)-4.0_dp*(u(1,j+1,:)+u(1,j-1,:))+u(1,j-2,:)))*delta
      u(2,j,:)=Oomehd(qn)*u(2,j,:)+(rhs(2,j,:) &
!      u(2,j,:)=(rhs(2,j,:) &
        -fx*(u(4,j,:)-4.0_dp*(u(3,j,:)+u(1,j,:))+rbufW(2,j,:)) &
        -fy*(u(2,j+2,:)-4.0_dp*(u(2,j+1,:)+u(2,j-1,:))+u(2,j-2,:)))*delta
    endif
    do i=3,nxf-2
      u(i,j,:)=Oomehd(qn)*u(i,j,:)+(rhs(i,j,:) &
!      u(i,j,:)=(rhs(i,j,:) &
        -fx*(u(i+2,j,:)-4.0_dp*(u(i+1,j,:)+u(i-1,j,:))+u(i-2,j,:)) &
        -fy*(u(i,j+2,:)-4.0_dp*(u(i,j+1,:)+u(i,j-1,:))+u(i,j-2,:)))*delta
    enddo
!...Right boundary
    if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
      u(nxf-1,j,:)=Oomehd(qn)*u(nxf-1,j,:)+(rhs(nxf-1,j,:) &
!      u(nxf-1,j,:)=(rhs(nxf-1,j,:) &
        -fx*(Rq*u(nxf,j,:)+Rs(1,:)-4.0_dp*(u(nxf,j,:)+u(nxf-2,j,:))+u(nxf-3,j,:)) &
        -fy*(u(nxf-1,j+2,:)-4.0_dp*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:))+u(nxf-1,j-2,:)))*delta
      u(nxf,j,:)=Oomehd(qn)*u(nxf,j,:)+(rhs(nxf,j,:) &
!      u(nxf,j,:)=(rhs(nxf,j,:) &
        -fx*(Rq*u(nxf-1,j,:)+Rs(2,:)-4.0_dp*(rs(1,:)+u(nxf-1,j,:))+u(nxf-2,j,:)) &
        -fy*(u(nxf,j+2,:)-4.0_dp*(u(nxf,j+1,:)+u(nxf,j-1,:))+u(nxf,j-2,:)))*deltaR
    else
      u(nxf-1,j,:)=Oomehd(qn)*u(nxf-1,j,:)+(rhs(nxf-1,j,:) &
!      u(nxf-1,j,:)=(rhs(nxf-1,j,:) &
        -fx*(rbufE(1,j,:)-4.0_dp*(u(nxf,j,:)+u(nxf-2,j,:))+u(nxf-3,j,:)) &
        -fy*(u(nxf-1,j+2,:)-4.0_dp*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:))+u(nxf-1,j-2,:)))*delta
      u(nxf,j,:)=Oomehd(qn)*u(nxf,j,:)+(rhs(nxf,j,:) &
!      u(nxf,j,:)=(rhs(nxf,j,:) &
        -fx*(rbufE(2,j,:)-4.0_dp*(rbufE(1,j,:)+u(nxf-1,j,:))+u(nxf-2,j,:)) &
        -fy*(u(nxf,j+2,:)-4.0_dp*(u(nxf,j+1,:)+u(nxf,j-1,:))+u(nxf,j-2,:)))*delta
    endif
  enddo
  call MPI_WAIT(req(8),stat(:,8),ierr)
!........................... SECOND TO TOP ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,nyf-1,:)=Oomehd(qn)*u(1,nyf-1,:)+(rhs(1,nyf-1,:) &
!    u(1,nyf-1,:)=(rhs(1,nyf-1,:) &
      -fx*(u(3,nyf-1,:)-4.0_dp*(u(2,nyf-1,:)+Ls(1,:))+Lq*u(2,nyf-1,:)+Ls(2,:)) &
      -fy*(rbufN(1,1,:)-4.0_dp*(u(1,nyf,:)+u(1,nyf-2,:))+u(1,nyf-3,:)))*deltaL
    u(2,nyf-1,:)=Oomehd(qn)*u(2,nyf-1,:)+(rhs(2,nyf-1,:) &
!    u(2,nyf-1,:)=(rhs(2,nyf-1,:) &
      -fx*(u(4,nyf-1,:)-4.0_dp*(u(3,nyf-1,:)+u(1,nyf-1,:))+Lq*u(1,nyf-1,:)+Ls(1,:)) &
      -fy*(rbufN(2,1,:)-4.0_dp*(u(2,nyf,:)+u(2,nyf-2,:))+u(2,nyf-3,:)))*delta
  else
    u(1,nyf-1,:)=Oomehd(qn)*u(1,nyf-1,:)+(rhs(1,nyf-1,:) &
!    u(1,nyf-1,:)=(rhs(1,nyf-1,:) &
      -fx*(u(3,nyf-1,:)-4.0_dp*(u(2,nyf-1,:)+rbufW(2,nyf-1,:))+rbufW(1,nyf-1,:)) &
      -fy*(rbufN(1,1,:)-4.0_dp*(u(1,nyf,:)+u(1,nyf-2,:))+u(1,nyf-3,:)))*delta
    u(2,nyf-1,:)=Oomehd(qn)*u(2,nyf-1,:)+(rhs(2,nyf-1,:) &
!    u(2,nyf-1,:)=(rhs(2,nyf-1,:) &
      -fx*(u(4,nyf-1,:)-4.0_dp*(u(3,nyf-1,:)+u(1,nyf-1,:))+rbufW(2,nyf-1,:)) &
      -fy*(rbufN(2,1,:)-4.0_dp*(u(2,nyf,:)+u(2,nyf-2,:))+u(2,nyf-3,:)))*delta
  endif
  do i=3,nxf-2
    u(i,nyf-1,:)=Oomehd(qn)*u(i,nyf-1,:)+(rhs(i,nyf-1,:) &
!    u(i,nyf-1,:)=(rhs(i,nyf-1,:) &
      -fx*(u(i+2,nyf-1,:)-4.0_dp*(u(i+1,nyf-1,:)+u(i-1,nyf-1,:))+u(i-2,nyf-1,:)) &
      -fy*(rbufN(i,1,:)-4.0_dp*(u(i,nyf,:)+u(i,nyf-2,:))+u(i,nyf-3,:)))*delta
  enddo
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-1,nyf-1,:)=Oomehd(qn)*u(nxf-1,nyf-1,:)+(rhs(nxf-1,nyf-1,:) &
!    u(nxf-1,nyf-1,:)=(rhs(nxf-1,nyf-1,:) &
     -fx*(Rq*u(nxf,nyf-1,:)+Rs(1,:)-4.0_dp*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:))+u(nxf-3,nyf-1,:)) &
     -fy*(rbufN(nxf-1,1,:)-4.0_dp*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:))+u(nxf-1,nyf-3,:)))*delta
    u(nxf,nyf-1,:)=Oomehd(qn)*u(nxf,nyf-1,:)+(rhs(nxf,nyf-1,:) &
!    u(nxf,nyf-1,:)=(rhs(nxf,nyf-1,:) &
     -fx*(Rq*u(nxf-1,nyf-1,:)+Rs(2,:)-4.0_dp*(rs(1,:)+u(nxf-1,nyf-1,:))+u(nxf-2,nyf-1,:)) &
     -fy*(rbufN(nxf,1,:)-4.0_dp*(u(nxf,nyf,:)+u(nxf,nyf-2,:))+u(nxf,nyf-3,:)))*deltaR
  else
    u(nxf-1,nyf-1,:)=Oomehd(qn)*u(nxf-1,nyf-1,:)+(rhs(nxf-1,nyf-1,:) &
!    u(nxf-1,nyf-1,:)=(rhs(nxf-1,nyf-1,:) &
      -fx*(rbufE(1,nyf-1,:)-4.0_dp*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:))+u(nxf-3,nyf-1,:)) &
      -fy*(rbufN(nxf-1,1,:)-4.0_dp*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:))+u(nxf-1,nyf-3,:)))*delta
    u(nxf,nyf-1,:)=Oomehd(qn)*u(nxf,nyf-1,:)+(rhs(nxf,nyf-1,:) &
!    u(nxf,nyf-1,:)=(rhs(nxf,nyf-1,:) &
      -fx*(rbufE(2,nyf-1,:)-4.0_dp*(rbufE(1,nyf-1,:)+u(nxf-1,nyf-1,:))+u(nxf-2,nyf-1,:)) &
      -fy*(rbufN(nxf,1,:)-4.0_dp*(u(nxf,nyf,:)+u(nxf,nyf-2,:))+u(nxf,nyf-3,:)))*delta
  endif
!........................... TOP ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,nyf,:)=Oomehd(qn)*u(1,nyf,:)+(rhs(1,nyf,:) &
!    u(1,nyf,:)=(rhs(1,nyf,:) &
      -fx*(u(3,nyf,:)-4.0_dp*(u(2,nyf,:)+Ls(1,:))+Lq*u(2,nyf,:)+Ls(2,:)) &
      -fy*(rbufN(1,2,:)-4.0_dp*(rbufN(1,1,:)+u(1,nyf-1,:))+u(1,nyf-2,:)))*deltaL
    u(2,nyf,:)=Oomehd(qn)*u(2,nyf,:)+(rhs(2,nyf,:) &
!    u(2,nyf,:)=(rhs(2,nyf,:) &
     -fx*(u(4,nyf,:)-4.0_dp*(u(3,nyf,:)+u(1,nyf,:))+Lq*u(1,nyf,:)+Ls(1,:)) &
     -fy*(rbufN(2,2,:)-4.0_dp*(rbufN(2,1,:)+u(2,nyf-1,:))+u(2,nyf-2,:)))*delta
  else
    u(1,nyf,:)=Oomehd(qn)*u(1,nyf,:)+(rhs(1,nyf,:) &
!    u(1,nyf,:)=(rhs(1,nyf,:) &
      -fx*(u(3,nyf,:)-4.0_dp*(u(2,nyf,:)+rbufW(2,nyf,:))+rbufW(1,nyf,:)) &
      -fy*(rbufN(1,2,:)-4.0_dp*(rbufN(1,1,:)+u(1,nyf-1,:))+u(1,nyf-2,:)))*delta
    u(2,nyf,:)=Oomehd(qn)*u(2,nyf,:)+(rhs(2,nyf,:) &
!    u(2,nyf,:)=(rhs(2,nyf,:) &
      -fx*(u(4,nyf,:)-4.0_dp*(u(3,nyf,:)+u(1,nyf,:))+rbufW(2,nyf,:)) &
      -fy*(rbufN(2,2,:)-4.0_dp*(rbufN(2,1,:)+u(2,nyf-1,:))+u(2,nyf-2,:)))*delta
  endif
  do i=3,nxf-2
    u(i,nyf,:)=Oomehd(qn)*u(i,nyf,:)+(rhs(i,nyf,:) &
!    u(i,nyf,:)=(rhs(i,nyf,:) &
      -fx*(u(i+2,nyf,:)-4.0_dp*(u(i+1,nyf,:)+u(i-1,nyf,:))+u(i-2,nyf,:)) &
      -fy*(rbufN(i,2,:)-4.0_dp*(rbufN(i,1,:)+u(i,nyf-1,:))+u(i,nyf-2,:)))*delta
  enddo
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-1,nyf,:)=Oomehd(qn)*u(nxf-1,nyf,:)+(rhs(nxf-1,nyf,:) &
!    u(nxf-1,nyf,:)=(rhs(nxf-1,nyf,:) &
      -fx*(Rq*u(nxf,nyf,:)+Rs(1,:)-4.0_dp*(u(nxf,nyf,:)+u(nxf-2,nyf,:))+u(nxf-3,nyf,:)) &
      -fy*(rbufN(nxf-1,2,:)-4.0_dp*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:))+u(nxf-1,nyf-2,:)))*delta
    u(nxf,nyf,:)=Oomehd(qn)*u(nxf,nyf,:)+(rhs(nxf,nyf,:) &
!    u(nxf,nyf,:)=(rhs(nxf,nyf,:) &
      -fx*(Rq*u(nxf-1,nyf,:)+Rs(2,:)-4.0_dp*(rs(1,:)+u(nxf-1,nyf,:))+u(nxf-2,nyf,:)) &
      -fy*(rbufN(nxf,2,:)-4.0_dp*(rbufN(nxf,1,:)+u(nxf,nyf-1,:))+u(nxf,nyf-2,:)))*deltaR
  else
    u(nxf-1,nyf,:)=Oomehd(qn)*u(nxf-1,nyf,:)+(rhs(nxf-1,nyf,:) &
!    u(nxf-1,nyf,:)=(rhs(nxf-1,nyf,:) &
      -fx*(rbufE(1,nyf,:)-4.0_dp*(u(nxf,nyf,:)+u(nxf-2,nyf,:))+u(nxf-3,nyf,:)) &
      -fy*(rbufN(nxf-1,2,:)-4.0_dp*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:))+u(nxf-1,nyf-2,:)))*delta
    u(nxf,nyf,:)=Oomehd(qn)*u(nxf,nyf,:)+(rhs(nxf,nyf,:) &
!    u(nxf,nyf,:)=(rhs(nxf,nyf,:) &
      -fx*(rbufE(2,nyf,:)-4.0_dp*(rbufE(1,nyf,:)+u(nxf-1,nyf,:))+u(nxf-2,nyf,:)) &
      -fy*(rbufN(nxf,2,:)-4.0_dp*(rbufN(nxf,1,:)+u(nxf,nyf-1,:))+u(nxf,nyf-2,:)))*delta
  endif

#ELIF (HDop==12)
  delta=omehd(qn)/(1.0_dp+20.0_dp*(fx+fy))
  allocate(deltaL(2),deltaR(2))
  deltaL(1)=omehd(qn)/(1.0_dp+20.0_dp*(fx+fy)-15.0_dp*fx*Lq)
  deltaL(2)=omehd(qn)/(1.0_dp+20.0_dp*(fx+fy)-fx*Lq)
  deltaR(1)=omehd(qn)/(1.0_dp+20.0_dp*(fx+fy)-15.0_dp*fx*Rq)
  deltaR(2)=omehd(qn)/(1.0_dp+20.0_dp*(fx+fy)-fx*Rq)

!........................... FIRST ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,1,:)=Oomehd(qn)*u(1,1,:)+(rhs(1,1,:) &
!    u(1,1,:)=(rhs(1,1,:) &
      +fx*(u(4,1,:)-6.0_dp*(u(3,1,:)+Lq*u(2,1,:)+Ls(2,:))+15.0_dp*(u(2,1,:)+Ls(1,:))+Lq*u(3,1,:)+Ls(3,:)) &
      +fy*(u(1,4,:)-6.0_dp*(u(1,3,:)+rbufS(1,2,:))+15.0_dp*(u(1,2,:)+rbufS(1,3,:))+rbufS(1,1,:)))*deltaL(1)
    u(2,1,:)=Oomehd(qn)*u(2,1,:)+(rhs(2,1,:) &
!    u(2,1,:)=(rhs(2,1,:) &
      +fx*(u(5,1,:)-6.0_dp*(u(4,1,:)+Lq*u(1,1,:)+Ls(1,:))+15.0_dp*(u(3,1,:)+u(1,1,:))+Ls(2,:)) &
      +fy*(u(2,4,:)-6.0_dp*(u(2,3,:)+rbufS(2,2,:))+15.0_dp*(u(2,2,:)+rbufS(2,3,:))+rbufS(2,1,:)))*deltaL(2)
    u(3,1,:)=Oomehd(qn)*u(3,1,:)+(rhs(3,1,:) &
!    u(3,1,:)=(rhs(3,1,:) &
      +fx*(u(6,1,:)-6.0_dp*(u(5,1,:)+u(1,1,:))+15.0_dp*(u(4,1,:)+u(2,1,:))+Lq*u(1,1,:)+Ls(1,:)) &
      +fy*(u(3,4,:)-6.0_dp*(u(3,3,:)+rbufS(3,2,:))+15.0_dp*(u(3,2,:)+rbufS(3,3,:))+rbufS(3,1,:)))*delta
  else
    u(1,1,:)=Oomehd(qn)*u(1,1,:)+(rhs(1,1,:) &
!    u(1,1,:)=(rhs(1,1,:) &
      +fx*(u(4,1,:)-6.0_dp*(u(3,1,:)+rbufW(2,1,:))+15.0_dp*(u(2,1,:)+rbufW(3,1,:))+rbufW(1,1,:)) &
      +fy*(u(1,4,:)-6.0_dp*(u(1,3,:)+rbufS(1,2,:))+15.0_dp*(u(1,2,:)+rbufS(1,3,:))+rbufS(1,1,:)))*delta
    u(2,1,:)=Oomehd(qn)*u(2,1,:)+(rhs(2,1,:) &
!    u(2,1,:)=(rhs(2,1,:) &
      +fx*(u(5,1,:)-6.0_dp*(u(4,1,:)+rbufW(3,1,:))+15.0_dp*(u(3,1,:)+u(1,1,:))+rbufW(2,1,:)) &
      +fy*(u(2,4,:)-6.0_dp*(u(2,3,:)+rbufS(2,2,:))+15.0_dp*(u(2,2,:)+rbufS(2,3,:))+rbufS(2,1,:)))*delta
    u(3,1,:)=Oomehd(qn)*u(3,1,:)+(rhs(3,1,:) &
!    u(3,1,:)=(rhs(3,1,:) &
      +fx*(u(6,1,:)-6.0_dp*(u(5,1,:)+u(1,1,:))+15.0_dp*(u(4,1,:)+u(2,1,:))+rbufW(3,1,:)) &
      +fy*(u(3,4,:)-6.0_dp*(u(3,3,:)+rbufS(3,2,:))+15.0_dp*(u(3,2,:)+rbufS(3,3,:))+rbufS(3,1,:)))*delta
  endif
  do i=4,nxf-3
    u(i,1,:)=Oomehd(qn)*u(i,1,:)+(rhs(i,1,:) &
!    u(i,1,:)=(rhs(i,1,:) &
      +fx*(u(i+3,1,:)-6.0_dp*(u(i+2,1,:)+u(i-2,1,:))+15.0_dp*(u(i+1,1,:)+u(i-1,1,:))+u(i-3,1,:)) &
      +fy*(u(i,4,:)-6.0_dp*(u(i,3,:)+rbufS(i,2,:))+15.0_dp*(u(i,2,:)+rbufS(i,3,:))+rbufS(i,1,:)))*delta
  enddo
  call MPI_WAIT(req(7),stat(:,7),ierr)
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-2,1,:)=Oomehd(qn)*u(nxf-2,1,:)+(rhs(nxf-2,1,:) &
!    u(nxf-2,1,:)=(rhs(nxf-2,1,:) &
      +fx*(Rq*u(nxf,1,:)+Rs(1,:)-6.0_dp*(u(nxf,1,:)+u(nxf-4,1,:))+15.0_dp*(u(nxf-1,1,:)+u(nxf-3,1,:))+u(nxf-5,1,:)) &
      +fy*(u(nxf-2,4,:)-6.0_dp*(u(nxf-2,3,:)+rbufS(nxf-2,2,:))+15.0_dp*(u(nxf-2,2,:)+rbufS(nxf-2,3,:))+rbufS(nxf-2,1,:)))*delta
    u(nxf-1,1,:)=Oomehd(qn)*u(nxf-1,1,:)+(rhs(nxf-1,1,:) &
!    u(nxf-1,1,:)=(rhs(nxf-1,1,:) &
      +fx*(rs(2,:)-6.0_dp*(Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-3,1,:))+15.0_dp*(u(nxf,1,:)+u(nxf-2,1,:))+u(nxf-4,1,:)) &
      +fy*(u(nxf-1,4,:)-6.0_dp*(u(nxf-1,3,:)+rbufS(nxf-1,2,:))+15.0_dp*(u(nxf-1,2,:)+rbufS(nxf-1,3,:))+rbufS(nxf-1,1,:)))*deltaR(2)
    u(nxf,1,:)=Oomehd(qn)*u(nxf,1,:)+(rhs(nxf,1,:) &
!    u(nxf,1,:)=(rhs(nxf,1,:) &
      +fx*(Rq*u(nxf-2,1,:)+Rs(3,:)-6.0_dp*(Rq*u(nxf-1,1,:)+Rs(2,:)+u(nxf-2,1,:))+15.0_dp*(rs(1,:)+u(nxf-1,1,:))+u(nxf-3,1,:)) &
      +fy*(u(nxf,4,:)-6.0_dp*(u(nxf,3,:)+rbufS(nxf,2,:))+15.0_dp*(u(nxf,2,:)+rbufS(nxf,3,:))+rbufS(nxf,1,:)))*deltaR(1)
  else
    u(nxf-2,1,:)=Oomehd(qn)*u(nxf-2,1,:)+(rhs(nxf-2,1,:) &
!    u(nxf-2,1,:)=(rhs(nxf-2,1,:) &
      +fx*(rbufE(1,1,:)-6.0_dp*(u(nxf,1,:)+u(nxf-4,1,:))+15.0_dp*(u(nxf-1,1,:)+u(nxf-3,1,:))+u(nxf-5,1,:)) &
      +fy*(u(nxf-2,4,:)-6.0_dp*(u(nxf-2,3,:)+rbufS(nxf-2,2,:))+15.0_dp*(u(nxf-2,2,:)+rbufS(nxf-2,3,:))+rbufS(nxf-2,1,:)))*delta
    u(nxf-1,1,:)=Oomehd(qn)*u(nxf-1,1,:)+(rhs(nxf-1,1,:) &
!    u(nxf-1,1,:)=(rhs(nxf-1,1,:) &
      +fx*(rbufE(2,1,:)-6.0_dp*(rbufE(1,1,:)+u(nxf-3,1,:))+15.0_dp*(u(nxf,1,:)+u(nxf-2,1,:))+u(nxf-4,1,:)) &
      +fy*(u(nxf-1,4,:)-6.0_dp*(u(nxf-1,3,:)+rbufS(nxf-1,2,:))+15.0_dp*(u(nxf-1,2,:)+rbufS(nxf-1,3,:))+rbufS(nxf-1,1,:)))*delta
    u(nxf,1,:)=Oomehd(qn)*u(nxf,1,:)+(rhs(nxf,1,:) &
!    u(nxf,1,:)=(rhs(nxf,1,:) &
      +fx*(rbufE(3,1,:)-6.0_dp*(rbufE(2,1,:)+u(nxf-2,1,:))+15.0_dp*(rbufE(1,1,:)+u(nxf-1,1,:))+u(nxf-3,1,:)) &
      +fy*(u(nxf,4,:)-6.0_dp*(u(nxf,3,:)+rbufS(nxf,2,:))+15.0_dp*(u(nxf,2,:)+rbufS(nxf,3,:))+rbufS(nxf,1,:)))*delta
  endif
!........................... SECOND ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,2,:)=Oomehd(qn)*u(1,2,:)+(rhs(1,2,:) &
!    u(1,2,:)=(rhs(1,2,:) &
      +fx*(u(4,2,:)-6.0_dp*(u(3,2,:)+Lq*u(2,2,:)+Ls(2,:))+15.0_dp*(u(2,2,:)+Ls(1,:))+Lq*u(3,2,:)+Ls(3,:)) &
      +fy*(u(1,5,:)-6.0_dp*(u(1,4,:)+rbufS(1,3,:))+15.0_dp*(u(1,3,:)+u(1,1,:))+rbufS(1,2,:)))*deltaL(1)
    u(2,2,:)=Oomehd(qn)*u(2,2,:)+(rhs(2,2,:) &
!    u(2,2,:)=(rhs(2,2,:) &
      +fx*(u(5,2,:)-6.0_dp*(u(4,2,:)+Lq*u(1,2,:)+Ls(1,:))+15.0_dp*(u(3,2,:)+u(1,2,:))+Ls(2,:)) &
      +fy*(u(2,5,:)-6.0_dp*(u(2,4,:)+rbufS(2,3,:))+15.0_dp*(u(2,3,:)+u(2,1,:))+rbufS(2,2,:)))*deltaL(2)
    u(3,2,:)=Oomehd(qn)*u(3,2,:)+(rhs(3,2,:) &
!    u(3,2,:)=(rhs(3,2,:) &
      +fx*(u(6,2,:)-6.0_dp*(u(5,2,:)+u(1,2,:))+15.0_dp*(u(4,2,:)+u(2,2,:))+Lq*u(1,2,:)+Ls(1,:)) &
      +fy*(u(3,5,:)-6.0_dp*(u(3,4,:)+rbufS(3,3,:))+15.0_dp*(u(3,3,:)+u(3,1,:))+rbufS(3,2,:)))*delta
  else
    u(1,2,:)=Oomehd(qn)*u(1,2,:)+(rhs(1,2,:) &
!    u(1,2,:)=(rhs(1,2,:) &
      +fx*(u(4,2,:)-6.0_dp*(u(3,2,:)+rbufW(2,2,:))+15.0_dp*(u(2,2,:)+rbufW(3,2,:))+rbufW(1,2,:)) &
      +fy*(u(1,5,:)-6.0_dp*(u(1,4,:)+rbufS(1,3,:))+15.0_dp*(u(1,3,:)+u(1,1,:))+rbufS(1,2,:)))*delta
    u(2,2,:)=Oomehd(qn)*u(2,2,:)+(rhs(2,2,:) &
!    u(2,2,:)=(rhs(2,2,:) &
      +fx*(u(5,2,:)-6.0_dp*(u(4,2,:)+rbufW(3,2,:))+15.0_dp*(u(3,2,:)+u(1,2,:))+rbufW(2,2,:)) &
      +fy*(u(2,5,:)-6.0_dp*(u(2,4,:)+rbufS(2,3,:))+15.0_dp*(u(2,3,:)+u(2,1,:))+rbufS(2,2,:)))*delta
    u(3,2,:)=Oomehd(qn)*u(3,2,:)+(rhs(3,2,:) &
!    u(3,2,:)=(rhs(3,2,:) &
      +fx*(u(6,2,:)-6.0_dp*(u(5,2,:)+u(1,2,:))+15.0_dp*(u(4,2,:)+u(2,2,:))+rbufW(3,2,:)) &
      +fy*(u(3,5,:)-6.0_dp*(u(3,4,:)+rbufS(3,3,:))+15.0_dp*(u(3,3,:)+u(3,1,:))+rbufS(3,2,:)))*delta
  endif
  do i=4,nxf-3
    u(i,2,:)=Oomehd(qn)*u(i,2,:)+(rhs(i,2,:) &
!    u(i,2,:)=(rhs(i,2,:) &
      +fx*(u(i+3,2,:)-6.0_dp*(u(i+2,2,:)+u(i-2,2,:))+15.0_dp*(u(i+1,2,:)+u(i-1,2,:))+u(i-3,2,:)) &
      +fy*(u(i,5,:)-6.0_dp*(u(i,4,:)+rbufS(i,3,:))+15.0_dp*(u(i,3,:)+u(i,1,:))+rbufS(i,2,:)))*delta
  enddo
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-2,2,:)=Oomehd(qn)*u(nxf-2,2,:)+(rhs(nxf-2,2,:) &
!    u(nxf-2,2,:)=(rhs(nxf-2,2,:) &
      +fx*(Rq*u(nxf,2,:)+Rs(1,:)-6.0_dp*(u(nxf,2,:)+u(nxf-4,2,:))+15.0_dp*(u(nxf-1,2,:)+u(nxf-3,2,:))+u(nxf-5,2,:)) &
      +fy*(u(nxf-2,5,:)-6.0_dp*(u(nxf-2,4,:)+rbufS(nxf-2,3,:))+15.0_dp*(u(nxf-2,3,:)+u(nxf-2,1,:))+rbufS(nxf-2,2,:)))*delta
    u(nxf-1,2,:)=Oomehd(qn)*u(nxf-1,2,:)+(rhs(nxf-1,2,:) &
!    u(nxf-1,2,:)=(rhs(nxf-1,2,:) &
      +fx*(rs(2,:)-6.0_dp*(Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-3,2,:))+15.0_dp*(u(nxf,2,:)+u(nxf-2,2,:))+u(nxf-4,2,:)) &
      +fy*(u(nxf-1,5,:)-6.0_dp*(u(nxf-1,4,:)+rbufS(nxf-1,3,:))+15.0_dp*(u(nxf-1,3,:)+u(nxf-1,1,:))+rbufS(nxf-1,2,:)))*deltaR(2)
    u(nxf,2,:)=Oomehd(qn)*u(nxf,2,:)+(rhs(nxf,2,:) &
!    u(nxf,2,:)=(rhs(nxf,2,:) &
      +fx*(Rq*u(nxf-2,2,:)+Rs(3,:)-6.0_dp*(Rq*u(nxf-1,2,:)+Rs(2,:)+u(nxf-2,2,:))+15.0_dp*(rs(1,:)+u(nxf-1,2,:))+u(nxf-3,2,:)) &
      +fy*(u(nxf,5,:)-6.0_dp*(u(nxf,4,:)+rbufS(nxf,3,:))+15.0_dp*(u(nxf,3,:)+u(nxf,1,:))+rbufS(nxf,2,:)))*deltaR(1)
  else
    u(nxf-2,2,:)=Oomehd(qn)*u(nxf-2,2,:)+(rhs(nxf-2,2,:) &
!    u(nxf-2,2,:)=(rhs(nxf-2,2,:) &
      +fx*(rbufE(1,2,:)-6.0_dp*(u(nxf,2,:)+u(nxf-4,2,:))+15.0_dp*(u(nxf-1,2,:)+u(nxf-3,2,:))+u(nxf-5,2,:)) &
      +fy*(u(nxf-2,5,:)-6.0_dp*(u(nxf-2,4,:)+rbufS(nxf-2,3,:))+15.0_dp*(u(nxf-2,3,:)+u(nxf-2,1,:))+rbufS(nxf-2,2,:)))*delta
    u(nxf-1,2,:)=Oomehd(qn)*u(nxf-1,2,:)+(rhs(nxf-1,2,:) &
!    u(nxf-1,2,:)=(rhs(nxf-1,2,:) &
      +fx*(rbufE(2,2,:)-6.0_dp*(rbufE(1,2,:)+u(nxf-3,2,:))+15.0_dp*(u(nxf,2,:)+u(nxf-2,2,:))+u(nxf-4,2,:)) &
      +fy*(u(nxf-1,5,:)-6.0_dp*(u(nxf-1,4,:)+rbufS(nxf-1,3,:))+15.0_dp*(u(nxf-1,3,:)+u(nxf-1,1,:))+rbufS(nxf-1,2,:)))*delta
    u(nxf,2,:)=Oomehd(qn)*u(nxf,2,:)+(rhs(nxf,2,:) &
!    u(nxf,2,:)=(rhs(nxf,2,:) &
      +fx*(rbufE(3,2,:)-6.0_dp*(rbufE(2,2,:)+u(nxf-2,2,:))+15.0_dp*(rbufE(1,2,:)+u(nxf-1,2,:))+u(nxf-3,2,:)) &
      +fy*(u(nxf,5,:)-6.0_dp*(u(nxf,4,:)+rbufS(nxf,3,:))+15.0_dp*(u(nxf,3,:)+u(nxf,1,:))+rbufS(nxf,2,:)))*delta
  endif
!........................... THIRD ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,3,:)=Oomehd(qn)*u(1,3,:)+(rhs(1,3,:) &
!    u(1,3,:)=(rhs(1,3,:) &
      +fx*(u(4,3,:)-6.0_dp*(u(3,3,:)+Lq*u(2,3,:)+Ls(2,:))+15.0_dp*(u(2,3,:)+Ls(1,:))+Lq*u(3,3,:)+Ls(3,:)) &
      +fy*(u(1,6,:)-6.0_dp*(u(1,5,:)+u(1,1,:))+15.0_dp*(u(1,4,:)+u(1,2,:))+rbufS(1,3,:)))*deltaL(1)
    u(2,3,:)=Oomehd(qn)*u(2,3,:)+(rhs(2,3,:) &
!    u(2,3,:)=(rhs(2,3,:) &
      +fx*(u(5,3,:)-6.0_dp*(u(4,3,:)+Lq*u(1,3,:)+Ls(1,:))+15.0_dp*(u(3,3,:)+u(1,3,:))+Ls(2,:)) &
      +fy*(u(2,6,:)-6.0_dp*(u(2,5,:)+u(2,1,:))+15.0_dp*(u(2,4,:)+u(2,2,:))+rbufS(2,3,:)))*deltaL(2)
    u(3,3,:)=Oomehd(qn)*u(3,3,:)+(rhs(3,3,:) &
!    u(3,3,:)=(rhs(3,3,:) &
      +fx*(u(6,3,:)-6.0_dp*(u(5,3,:)+u(1,3,:))+15.0_dp*(u(4,3,:)+u(2,3,:))+Lq*u(1,3,:)+Ls(1,:)) &
      +fy*(u(3,6,:)-6.0_dp*(u(3,5,:)+u(3,1,:))+15.0_dp*(u(3,4,:)+u(3,2,:))+rbufS(3,3,:)))*delta
  else
    u(1,3,:)=Oomehd(qn)*u(1,3,:)+(rhs(1,3,:) &
!    u(1,3,:)=(rhs(1,3,:) &
      +fx*(u(4,3,:)-6.0_dp*(u(3,3,:)+rbufW(2,3,:))+15.0_dp*(u(2,3,:)+rbufW(3,3,:))+rbufW(1,3,:)) &
      +fy*(u(1,6,:)-6.0_dp*(u(1,5,:)+u(1,1,:))+15.0_dp*(u(1,4,:)+u(1,2,:))+rbufS(1,3,:)))*delta
    u(2,3,:)=Oomehd(qn)*u(2,3,:)+(rhs(2,3,:) &
!    u(2,3,:)=(rhs(2,3,:) &
      +fx*(u(5,3,:)-6.0_dp*(u(4,3,:)+rbufW(3,3,:))+15.0_dp*(u(3,3,:)+u(1,3,:))+rbufW(2,3,:)) &
      +fy*(u(2,6,:)-6.0_dp*(u(2,5,:)+u(2,1,:))+15.0_dp*(u(2,4,:)+u(2,2,:))+rbufS(2,3,:)))*delta
    u(3,3,:)=Oomehd(qn)*u(3,3,:)+(rhs(3,3,:) &
!    u(3,3,:)=(rhs(3,3,:) &
      +fx*(u(6,3,:)-6.0_dp*(u(5,3,:)+u(1,3,:))+15.0_dp*(u(4,3,:)+u(2,3,:))+rbufW(3,3,:)) &
      +fy*(u(3,6,:)-6.0_dp*(u(3,5,:)+u(3,1,:))+15.0_dp*(u(3,4,:)+u(3,2,:))+rbufS(3,3,:)))*delta
  endif
  do i=4,nxf-3
    u(i,3,:)=Oomehd(qn)*u(i,3,:)+(rhs(i,3,:) &
!    u(i,3,:)=(rhs(i,3,:) &
      +fx*(u(i+3,3,:)-6.0_dp*(u(i+2,3,:)+u(i-2,3,:))+15.0_dp*(u(i+1,3,:)+u(i-1,3,:))+u(i-3,3,:)) &
      +fy*(u(i,6,:)-6.0_dp*(u(i,5,:)+u(i,1,:))+15.0_dp*(u(i,4,:)+u(i,2,:))+rbufS(i,3,:)))*delta
  enddo
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-2,3,:)=Oomehd(qn)*u(nxf-2,3,:)+(rhs(nxf-2,3,:) &
!    u(nxf-2,3,:)=(rhs(nxf-2,3,:) &
      +fx*(Rq*u(nxf,3,:)+Rs(1,:)-6.0_dp*(u(nxf,3,:)+u(nxf-4,3,:))+15.0_dp*(u(nxf-1,3,:)+u(nxf-3,3,:))+u(nxf-5,3,:)) &
      +fy*(u(nxf-2,6,:)-6.0_dp*(u(nxf-2,5,:)+u(nxf-2,1,:))+15.0_dp*(u(nxf-2,4,:)+u(nxf-2,2,:))+rbufS(nxf-2,3,:)))*delta
    u(nxf-1,3,:)=Oomehd(qn)*u(nxf-1,3,:)+(rhs(nxf-1,3,:) &
!    u(nxf-1,3,:)=(rhs(nxf-1,3,:) &
      +fx*(rs(2,:)-6.0_dp*(Rq*u(nxf,3,:)+Rs(1,:)+u(nxf-3,3,:))+15.0_dp*(u(nxf,3,:)+u(nxf-2,3,:))+u(nxf-4,3,:)) &
      +fy*(u(nxf-1,6,:)-6.0_dp*(u(nxf-1,5,:)+u(nxf-1,1,:))+15.0_dp*(u(nxf-1,4,:)+u(nxf-1,2,:))+rbufS(nxf-1,3,:)))*deltaR(2)
    u(nxf,3,:)=Oomehd(qn)*u(nxf,3,:)+(rhs(nxf,3,:) &
!    u(nxf,3,:)=(rhs(nxf,3,:) &
      +fx*(Rq*u(nxf-2,3,:)+Rs(3,:)-6.0_dp*(Rq*u(nxf-1,3,:)+Rs(2,:)+u(nxf-2,3,:))+15.0_dp*(rs(1,:)+u(nxf-1,3,:))+u(nxf-3,3,:)) &
      +fy*(u(nxf,6,:)-6.0_dp*(u(nxf,5,:)+u(nxf,1,:))+15.0_dp*(u(nxf,4,:)+u(nxf,2,:))+rbufS(nxf,3,:)))*deltaR(1)
  else
    u(nxf-2,3,:)=Oomehd(qn)*u(nxf-2,3,:)+(rhs(nxf-2,3,:) &
!    u(nxf-2,3,:)=(rhs(nxf-2,3,:) &
      +fx*(rbufE(1,3,:)-6.0_dp*(u(nxf,3,:)+u(nxf-4,3,:))+15.0_dp*(u(nxf-1,3,:)+u(nxf-3,3,:))+u(nxf-5,3,:)) &
      +fy*(u(nxf-2,6,:)-6.0_dp*(u(nxf-2,5,:)+u(nxf-2,1,:))+15.0_dp*(u(nxf-2,4,:)+u(nxf-2,2,:))+rbufS(nxf-2,3,:)))*delta
    u(nxf-1,3,:)=Oomehd(qn)*u(nxf-1,3,:)+(rhs(nxf-1,3,:) &
!    u(nxf-1,3,:)=(rhs(nxf-1,3,:) &
      +fx*(rbufE(2,3,:)-6.0_dp*(rbufE(1,3,:)+u(nxf-3,3,:))+15.0_dp*(u(nxf,3,:)+u(nxf-2,3,:))+u(nxf-4,3,:)) &
      +fy*(u(nxf-1,6,:)-6.0_dp*(u(nxf-1,5,:)+u(nxf-1,1,:))+15.0_dp*(u(nxf-1,4,:)+u(nxf-1,2,:))+rbufS(nxf-1,3,:)))*delta
    u(nxf,3,:)=Oomehd(qn)*u(nxf,3,:)+(rhs(nxf,3,:) &
!    u(nxf,3,:)=(rhs(nxf,3,:) &
      +fx*(rbufE(3,3,:)-6.0_dp*(rbufE(2,3,:)+u(nxf-2,3,:))+15.0_dp*(rbufE(1,3,:)+u(nxf-1,3,:))+u(nxf-3,3,:)) &
      +fy*(u(nxf,6,:)-6.0_dp*(u(nxf,5,:)+u(nxf,1,:))+15.0_dp*(u(nxf,4,:)+u(nxf,2,:))+rbufS(nxf,3,:)))*delta
  endif
!................... AWAY FROM BOTTOM AND TOP BOUNDARIES ...................!
  do j=4,nyf-3
!...Left boundary
    if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
      u(1,j,:)=Oomehd(qn)*u(1,j,:)+(rhs(1,j,:) &
!      u(1,j,:)=(rhs(1,j,:) &
        +fx*(u(4,j,:)-6.0_dp*(u(3,j,:)+Lq*u(2,j,:)+Ls(2,:))+15.0_dp*(u(2,j,:)+Ls(1,:))+Lq*u(3,j,:)+Ls(3,:)) &
        +fy*(u(1,j+3,:)-6.0_dp*(u(1,j+2,:)+u(1,j-2,:))+15.0_dp*(u(1,j+1,:)+u(1,j-1,:))+u(1,j-3,:)))*deltaL(1)
      u(2,j,:)=Oomehd(qn)*u(2,j,:)+(rhs(2,j,:) &
!      u(2,j,:)=(rhs(2,j,:) &
        +fx*(u(5,j,:)-6.0_dp*(u(4,j,:)+Lq*u(1,j,:)+Ls(1,:))+15.0_dp*(u(3,j,:)+u(1,j,:))+Ls(2,:)) &
        +fy*(u(2,j+3,:)-6.0_dp*(u(2,j+2,:)+u(2,j-2,:))+15.0_dp*(u(2,j+1,:)+u(2,j-1,:))+u(2,j-3,:)))*deltaL(2)
      u(3,j,:)=Oomehd(qn)*u(3,j,:)+(rhs(3,j,:) &
!      u(3,j,:)=(rhs(3,j,:) &
        +fx*(u(6,j,:)-6.0_dp*(u(5,j,:)+u(1,j,:))+15.0_dp*(u(4,j,:)+u(2,j,:))+Lq*u(1,j,:)+Ls(1,:)) &
        +fy*(u(3,j+3,:)-6.0_dp*(u(3,j+2,:)+u(3,j-2,:))+15.0_dp*(u(3,j+1,:)+u(3,j-1,:))+u(3,j-3,:)))*delta
    else
      u(1,j,:)=Oomehd(qn)*u(1,j,:)+(rhs(1,j,:) &
!      u(1,j,:)=(rhs(1,j,:) &
        +fx*(u(4,j,:)-6.0_dp*(u(3,j,:)+rbufW(2,j,:))+15.0_dp*(u(2,j,:)+rbufW(3,j,:))+rbufW(1,j,:)) &
        +fy*(u(1,j+3,:)-6.0_dp*(u(1,j+2,:)+u(1,j-2,:))+15.0_dp*(u(1,j+1,:)+u(1,j-1,:))+u(1,j-3,:)))*delta
      u(2,j,:)=Oomehd(qn)*u(2,j,:)+(rhs(2,j,:) &
!      u(2,j,:)=(rhs(2,j,:) &
        +fx*(u(5,j,:)-6.0_dp*(u(4,j,:)+rbufW(3,j,:))+15.0_dp*(u(3,j,:)+u(1,j,:))+rbufW(2,j,:)) &
        +fy*(u(2,j+3,:)-6.0_dp*(u(2,j+2,:)+u(2,j-2,:))+15.0_dp*(u(2,j+1,:)+u(2,j-1,:))+u(2,j-3,:)))*delta
      u(3,j,:)=Oomehd(qn)*u(3,j,:)+(rhs(3,j,:) &
!      u(3,j,:)=(rhs(3,j,:) &
        +fx*(u(6,j,:)-6.0_dp*(u(5,j,:)+u(1,j,:))+15.0_dp*(u(4,j,:)+u(2,j,:))+rbufW(3,j,:)) &
        +fy*(u(3,j+3,:)-6.0_dp*(u(3,j+2,:)+u(3,j-2,:))+15.0_dp*(u(3,j+1,:)+u(3,j-1,:))+u(3,j-3,:)))*delta
    endif
    do i=4,nxf-3
      u(i,j,:)=Oomehd(qn)*u(i,j,:)+(rhs(i,j,:) &
!      u(i,j,:)=(rhs(i,j,:) &
        +fx*(u(i+3,j,:)-6.0_dp*(u(i+2,j,:)+u(i-2,j,:))+15.0_dp*(u(i+1,j,:)+u(i-1,j,:))+u(i-3,j,:)) &
        +fy*(u(i,j+3,:)-6.0_dp*(u(i,j+2,:)+u(i,j-2,:))+15.0_dp*(u(i,j+1,:)+u(i,j-1,:))+u(i,j-3,:)))*delta
    enddo
!...Right boundary
    if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
      u(nxf-2,j,:)=Oomehd(qn)*u(nxf-2,j,:)+(rhs(nxf-2,j,:) &
!      u(nxf-2,j,:)=(rhs(nxf-2,j,:) &
        +fx*(Rq*u(nxf,j,:)+Rs(1,:)-6.0_dp*(u(nxf,j,:)+u(nxf-4,j,:))+15.0_dp*(u(nxf-1,j,:)+u(nxf-3,j,:))+u(nxf-5,j,:)) &
        +fy*(u(nxf-2,j+3,:)-6.0_dp*(u(nxf-2,j+2,:)+u(nxf-2,j-2,:))+15.0_dp*(u(nxf-2,j+1,:)+u(nxf-2,j-1,:))+u(nxf-2,j-3,:)))*delta
      u(nxf-1,j,:)=Oomehd(qn)*u(nxf-1,j,:)+(rhs(nxf-1,j,:) &
!      u(nxf-1,j,:)=(rhs(nxf-1,j,:) &
        +fx*(rs(2,:)-6.0_dp*(Rq*u(nxf,j,:)+Rs(1,:)+u(nxf-3,j,:))+15.0_dp*(u(nxf,j,:)+u(nxf-2,j,:))+u(nxf-4,j,:)) &
        +fy*(u(nxf-1,j+3,:)-6.0_dp*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:))+15.0_dp*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:))+u(nxf-1,j-3,:)))*deltaR(2)
      u(nxf,j,:)=Oomehd(qn)*u(nxf,j,:)+(rhs(nxf,j,:) &
!      u(nxf,j,:)=(rhs(nxf,j,:) &
        +fx*(Rq*u(nxf-2,j,:)+Rs(3,:)-6.0_dp*(Rq*u(nxf-1,j,:)+Rs(2,:)+u(nxf-2,j,:))+15.0_dp*(rs(1,:)+u(nxf-1,j,:))+u(nxf-3,j,:)) &
        +fy*(u(nxf,j+3,:)-6.0_dp*(u(nxf,j+2,:)+u(nxf,j-2,:))+15.0_dp*(u(nxf,j+1,:)+u(nxf,j-1,:))+u(nxf,j-3,:)))*deltaR(1)
    else
      u(nxf-2,j,:)=Oomehd(qn)*u(nxf-2,j,:)+(rhs(nxf-2,j,:) &
!      u(nxf-2,j,:)=(rhs(nxf-2,j,:) &
        +fx*(rbufE(1,j,:)-6.0_dp*(u(nxf,j,:)+u(nxf-4,j,:))+15.0_dp*(u(nxf-1,j,:)+u(nxf-3,j,:))+u(nxf-5,j,:)) &
        +fy*(u(nxf-2,j+3,:)-6.0_dp*(u(nxf-2,j+2,:)+u(nxf-2,j-2,:))+15.0_dp*(u(nxf-2,j+1,:)+u(nxf-2,j-1,:))+u(nxf-2,j-3,:)))*delta
      u(nxf-1,j,:)=Oomehd(qn)*u(nxf-1,j,:)+(rhs(nxf-1,j,:) &
!      u(nxf-1,j,:)=(rhs(nxf-1,j,:) &
        +fx*(rbufE(2,j,:)-6.0_dp*(rbufE(1,j,:)+u(nxf-3,j,:))+15.0_dp*(u(nxf,j,:)+u(nxf-2,j,:))+u(nxf-4,j,:)) &
        +fy*(u(nxf-1,j+3,:)-6.0_dp*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:))+15.0_dp*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:))+u(nxf-1,j-3,:)))*delta
      u(nxf,j,:)=Oomehd(qn)*u(nxf,j,:)+(rhs(nxf,j,:) &
!      u(nxf,j,:)=(rhs(nxf,j,:) &
        +fx*(rbufE(3,j,:)-6.0_dp*(rbufE(2,j,:)+u(nxf-2,j,:))+15.0_dp*(rbufE(1,j,:)+u(nxf-1,j,:))+u(nxf-3,j,:)) &
        +fy*(u(nxf,j+3,:)-6.0_dp*(u(nxf,j+2,:)+u(nxf,j-2,:))+15.0_dp*(u(nxf,j+1,:)+u(nxf,j-1,:))+u(nxf,j-3,:)))*delta
    endif
  enddo
  call MPI_WAIT(req(8),stat(:,8),ierr)
!........................... THIRD TO TOP ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,nyf-2,:)=Oomehd(qn)*u(1,nyf-2,:)+(rhs(1,nyf-2,:) &
!    u(1,nyf-2,:)=(rhs(1,nyf-2,:) &
      +fx*(u(4,nyf-2,:)-6.0_dp*(u(3,nyf-2,:)+Lq*u(2,nyf-2,:)+Ls(2,:))+15.0_dp*(u(2,nyf-2,:)+Ls(1,:))+Lq*u(3,nyf-2,:)+Ls(3,:)) &
      +fy*(rbufN(1,1,:)-6.0_dp*(u(1,nyf,:)+u(1,nyf-4,:))+15.0_dp*(u(1,nyf-1,:)+u(1,nyf-3,:))+u(1,nyf-5,:)))*deltaL(1)
    u(2,nyf-2,:)=Oomehd(qn)*u(2,nyf-2,:)+(rhs(2,nyf-2,:) &
!    u(2,nyf-2,:)=(rhs(2,nyf-2,:) &
      +fx*(u(5,nyf-2,:)-6.0_dp*(u(4,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:))+15.0_dp*(u(3,nyf-2,:)+u(1,nyf-2,:))+Ls(2,:)) &
      +fy*(rbufN(2,1,:)-6.0_dp*(u(2,nyf,:)+u(2,nyf-4,:))+15.0_dp*(u(2,nyf-1,:)+u(2,nyf-3,:))+u(2,nyf-5,:)))*deltaL(2)
    u(3,nyf-2,:)=Oomehd(qn)*u(3,nyf-2,:)+(rhs(3,nyf-2,:) &
!    u(3,nyf-2,:)=(rhs(3,nyf-2,:) &
      +fx*(u(6,nyf-2,:)-6.0_dp*(u(5,nyf-2,:)+u(1,nyf-2,:))+15.0_dp*(u(4,nyf-2,:)+u(2,nyf-2,:))+Lq*u(1,nyf-2,:)+Ls(1,:)) &
      +fy*(rbufN(3,1,:)-6.0_dp*(u(3,nyf,:)+u(3,nyf-4,:))+15.0_dp*(u(3,nyf-1,:)+u(3,nyf-3,:))+u(3,nyf-5,:)))*delta
  else
    u(1,nyf-2,:)=Oomehd(qn)*u(1,nyf-2,:)+(rhs(1,nyf-2,:) &
!    u(1,nyf-2,:)=(rhs(1,nyf-2,:) &
      +fx*(u(4,nyf-2,:)-6.0_dp*(u(3,nyf-2,:)+rbufW(2,nyf-2,:))+15.0_dp*(u(2,nyf-2,:)+rbufW(3,nyf-2,:))+rbufW(1,nyf-2,:)) &
      +fy*(rbufN(1,1,:)-6.0_dp*(u(1,nyf,:)+u(1,nyf-4,:))+15.0_dp*(u(1,nyf-1,:)+u(1,nyf-3,:))+u(1,nyf-5,:)))*delta
    u(2,nyf-2,:)=Oomehd(qn)*u(2,nyf-2,:)+(rhs(2,nyf-2,:) &
!    u(2,nyf-2,:)=(rhs(2,nyf-2,:) &
      +fx*(u(5,nyf-2,:)-6.0_dp*(u(4,nyf-2,:)+rbufW(3,nyf-2,:))+15.0_dp*(u(3,nyf-2,:)+u(1,nyf-2,:))+rbufW(2,nyf-2,:)) &
      +fy*(rbufN(2,1,:)-6.0_dp*(u(2,nyf,:)+u(2,nyf-4,:))+15.0_dp*(u(2,nyf-1,:)+u(2,nyf-3,:))+u(2,nyf-5,:)))*delta
    u(3,nyf-2,:)=Oomehd(qn)*u(3,nyf-2,:)+(rhs(3,nyf-2,:) &
!    u(3,nyf-2,:)=(rhs(3,nyf-2,:) &
      +fx*(u(6,nyf-2,:)-6.0_dp*(u(5,nyf-2,:)+u(1,nyf-2,:))+15.0_dp*(u(4,nyf-2,:)+u(2,nyf-2,:))+rbufW(3,nyf-2,:)) &
      +fy*(rbufN(3,1,:)-6.0_dp*(u(3,nyf,:)+u(3,nyf-4,:))+15.0_dp*(u(3,nyf-1,:)+u(3,nyf-3,:))+u(3,nyf-5,:)))*delta
  endif
  do i=4,nxf-3
    u(i,nyf-2,:)=Oomehd(qn)*u(i,nyf-2,:)+(rhs(i,nyf-2,:) &
!    u(i,nyf-2,:)=(rhs(i,nyf-2,:) &
      +fx*(u(i+3,nyf-2,:)-6.0_dp*(u(i+2,nyf-2,:)+u(i-2,nyf-2,:))+15.0_dp*(u(i+1,nyf-2,:)+u(i-1,nyf-2,:))+u(i-3,nyf-2,:)) &
      +fy*(rbufN(i,1,:)-6.0_dp*(u(i,nyf,:)+u(i,nyf-4,:))+15.0_dp*(u(i,nyf-1,:)+u(i,nyf-3,:))+u(i,nyf-5,:)))*delta
  enddo
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-2,nyf-2,:)=Oomehd(qn)*u(nxf-2,nyf-2,:)+(rhs(nxf-2,nyf-2,:) &
!    u(nxf-2,nyf-2,:)=(rhs(nxf-2,nyf-2,:) &
      +fx*(Rq*u(nxf,nyf-2,:)+Rs(1,:)-6.0_dp*(u(nxf,nyf-2,:)+u(nxf-4,nyf-2,:))+15.0_dp*(u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:))+u(nxf-5,nyf-2,:)) &
      +fy*(rbufN(nxf-2,1,:)-6.0_dp*(u(nxf-2,nyf,:)+u(nxf-2,nyf-4,:))+15.0_dp*(u(nxf-2,nyf-1,:)+u(nxf-2,nyf-3,:))+u(nxf-2,nyf-5,:)))*delta
    u(nxf-1,nyf-2,:)=Oomehd(qn)*u(nxf-1,nyf-2,:)+(rhs(nxf-1,nyf-2,:) &
!    u(nxf-1,nyf-2,:)=(rhs(nxf-1,nyf-2,:) &
      +fx*(rs(2,:)-6.0_dp*(Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-3,nyf-2,:))+15.0_dp*(u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:))+u(nxf-4,nyf-2,:)) &
      +fy*(rbufN(nxf-1,1,:)-6.0_dp*(u(nxf-1,nyf,:)+u(nxf-1,nyf-4,:))+15.0_dp*(u(nxf-1,nyf-1,:)+u(nxf-1,nyf-3,:))+u(nxf-1,nyf-5,:)))*deltaR(2)
    u(nxf,nyf-2,:)=Oomehd(qn)*u(nxf,nyf-2,:)+(rhs(nxf,nyf-2,:) &
!    u(nxf,nyf-2,:)=(rhs(nxf,nyf-2,:) &
      +fx*(Rq*u(nxf-2,nyf-2,:)+Rs(3,:)-6.0_dp*(Rq*u(nxf-1,nyf-2,:)+Rs(2,:)+u(nxf-2,nyf-2,:))+15.0_dp*(rs(1,:)+u(nxf-1,nyf-2,:))+u(nxf-3,nyf-2,:)) &
      +fy*(rbufN(nxf,1,:)-6.0_dp*(u(nxf,nyf,:)+u(nxf,nyf-4,:))+15.0_dp*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:))+u(nxf,nyf-5,:)))*deltaR(1)
  else
    u(nxf-2,nyf-2,:)=Oomehd(qn)*u(nxf-2,nyf-2,:)+(rhs(nxf-2,nyf-2,:) &
!    u(nxf-2,nyf-2,:)=(rhs(nxf-2,nyf-2,:) &
      +fx*(rbufE(1,nyf-2,:)-6.0_dp*(u(nxf,nyf-2,:)+u(nxf-4,nyf-2,:))+15.0_dp*(u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:))+u(nxf-5,nyf-2,:)) &
      +fy*(rbufN(nxf-2,1,:)-6.0_dp*(u(nxf-2,nyf,:)+u(nxf-2,nyf-4,:))+15.0_dp*(u(nxf-2,nyf-1,:)+u(nxf-2,nyf-3,:))+u(nxf-2,nyf-5,:)))*delta
    u(nxf-1,nyf-2,:)=Oomehd(qn)*u(nxf-1,nyf-2,:)+(rhs(nxf-1,nyf-2,:) &
!    u(nxf-1,nyf-2,:)=(rhs(nxf-1,nyf-2,:) &
      +fx*(rbufE(2,nyf-2,:)-6.0_dp*(rbufE(1,nyf-2,:)+u(nxf-3,nyf-2,:))+15.0_dp*(u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:))+u(nxf-4,nyf-2,:)) &
      +fy*(rbufN(nxf-1,1,:)-6.0_dp*(u(nxf-1,nyf,:)+u(nxf-1,nyf-4,:))+15.0_dp*(u(nxf-1,nyf-1,:)+u(nxf-1,nyf-3,:))+u(nxf-1,nyf-5,:)))*delta
    u(nxf,nyf-2,:)=Oomehd(qn)*u(nxf,nyf-2,:)+(rhs(nxf,nyf-2,:) &
!    u(nxf,nyf-2,:)=(rhs(nxf,nyf-2,:) &
      +fx*(rbufE(3,nyf-2,:)-6.0_dp*(rbufE(2,nyf-2,:)+u(nxf-2,nyf-2,:))+15.0_dp*(rbufE(1,nyf-2,:)+u(nxf-1,nyf-2,:))+u(nxf-3,nyf-2,:)) &
      +fy*(rbufN(nxf,1,:)-6.0_dp*(u(nxf,nyf,:)+u(nxf,nyf-4,:))+15.0_dp*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:))+u(nxf,nyf-5,:)))*delta
  endif
!........................... SECOND TO TOP ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,nyf-1,:)=Oomehd(qn)*u(1,nyf-1,:)+(rhs(1,nyf-1,:) &
!    u(1,nyf-1,:)=(rhs(1,nyf-1,:) &
      +fx*(u(4,nyf-1,:)-6.0_dp*(u(3,nyf-1,:)+Lq*u(2,nyf-1,:)+Ls(2,:))+15.0_dp*(u(2,nyf-1,:)+Ls(1,:))+Lq*u(3,nyf-1,:)+Ls(3,:)) &
      +fy*(rbufN(1,2,:)-6.0_dp*(rbufN(1,1,:)+u(1,nyf-3,:))+15.0_dp*(u(1,nyf,:)+u(1,nyf-2,:))+u(1,nyf-4,:)))*deltaL(1)
    u(2,nyf-1,:)=Oomehd(qn)*u(2,nyf-1,:)+(rhs(2,nyf-1,:) &
!    u(2,nyf-1,:)=(rhs(2,nyf-1,:) &
      +fx*(u(5,nyf-1,:)-6.0_dp*(u(4,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:))+15.0_dp*(u(3,nyf-1,:)+u(1,nyf-1,:))+Ls(2,:)) &
      +fy*(rbufN(2,2,:)-6.0_dp*(rbufN(2,1,:)+u(2,nyf-3,:))+15.0_dp*(u(2,nyf,:)+u(2,nyf-2,:))+u(2,nyf-4,:)))*deltaL(2)
    u(3,nyf-1,:)=Oomehd(qn)*u(3,nyf-1,:)+(rhs(3,nyf-1,:) &
!    u(3,nyf-1,:)=(rhs(3,nyf-1,:) &
      +fx*(u(6,nyf-1,:)-6.0_dp*(u(5,nyf-1,:)+u(1,nyf-1,:))+15.0_dp*(u(4,nyf-1,:)+u(2,nyf-1,:))+Lq*u(1,nyf-1,:)+Ls(1,:)) &
      +fy*(rbufN(3,2,:)-6.0_dp*(rbufN(3,1,:)+u(3,nyf-3,:))+15.0_dp*(u(3,nyf,:)+u(3,nyf-2,:))+u(3,nyf-4,:)))*delta
  else
    u(1,nyf-1,:)=Oomehd(qn)*u(1,nyf-1,:)+(rhs(1,nyf-1,:) &
!    u(1,nyf-1,:)=(rhs(1,nyf-1,:) &
      +fx*(u(4,nyf-1,:)-6.0_dp*(u(3,nyf-1,:)+rbufW(2,nyf-1,:))+15.0_dp*(u(2,nyf-1,:)+rbufW(3,nyf-1,:))+rbufW(1,nyf-1,:)) &
      +fy*(rbufN(1,2,:)-6.0_dp*(rbufN(1,1,:)+u(1,nyf-3,:))+15.0_dp*(u(1,nyf,:)+u(1,nyf-2,:))+u(1,nyf-4,:)))*delta
    u(2,nyf-1,:)=Oomehd(qn)*u(2,nyf-1,:)+(rhs(2,nyf-1,:) &
!    u(2,nyf-1,:)=(rhs(2,nyf-1,:) &
      +fx*(u(5,nyf-1,:)-6.0_dp*(u(4,nyf-1,:)+rbufW(3,nyf-1,:))+15.0_dp*(u(3,nyf-1,:)+u(1,nyf-1,:))+rbufW(2,nyf-1,:)) &
      +fy*(rbufN(2,2,:)-6.0_dp*(rbufN(2,1,:)+u(2,nyf-3,:))+15.0_dp*(u(2,nyf,:)+u(2,nyf-2,:))+u(2,nyf-4,:)))*delta
    u(3,nyf-1,:)=Oomehd(qn)*u(3,nyf-1,:)+(rhs(3,nyf-1,:) &
!    u(3,nyf-1,:)=(rhs(3,nyf-1,:) &
      +fx*(u(6,nyf-1,:)-6.0_dp*(u(5,nyf-1,:)+u(1,nyf-1,:))+15.0_dp*(u(4,nyf-1,:)+u(2,nyf-1,:))+rbufW(3,nyf-1,:)) &
      +fy*(rbufN(3,2,:)-6.0_dp*(rbufN(3,1,:)+u(3,nyf-3,:))+15.0_dp*(u(3,nyf,:)+u(3,nyf-2,:))+u(3,nyf-4,:)))*delta
  endif
  do i=4,nxf-3
    u(i,nyf-1,:)=Oomehd(qn)*u(i,nyf-1,:)+(rhs(i,nyf-1,:) &
!    u(i,nyf-1,:)=(rhs(i,nyf-1,:) &
      +fx*(u(i+3,nyf-1,:)-6.0_dp*(u(i+2,nyf-1,:)+u(i-2,nyf-1,:))+15.0_dp*(u(i+1,nyf-1,:)+u(i-1,nyf-1,:))+u(i-3,nyf-1,:)) &
      +fy*(rbufN(i,2,:)-6.0_dp*(rbufN(i,1,:)+u(i,nyf-3,:))+15.0_dp*(u(i,nyf,:)+u(i,nyf-2,:))+u(i,nyf-4,:)))*delta
  enddo
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-2,nyf-1,:)=Oomehd(qn)*u(nxf-2,nyf-1,:)+(rhs(nxf-2,nyf-1,:) &
!    u(nxf-2,nyf-1,:)=(rhs(nxf-2,nyf-1,:) &
      +fx*(Rq*u(nxf,nyf-1,:)+Rs(1,:)-6.0_dp*(u(nxf,nyf-1,:)+u(nxf-4,nyf-1,:))+15.0_dp*(u(nxf-1,nyf-1,:)+u(nxf-3,nyf-1,:))+u(nxf-5,nyf-1,:)) &
      +fy*(rbufN(nxf-2,2,:)-6.0_dp*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-3,:))+15.0_dp*(u(nxf-2,nyf,:)+u(nxf-2,nyf-2,:))+u(nxf-2,nyf-4,:)))*delta
    u(nxf-1,nyf-1,:)=Oomehd(qn)*u(nxf-1,nyf-1,:)+(rhs(nxf-1,nyf-1,:) &
!    u(nxf-1,nyf-1,:)=(rhs(nxf-1,nyf-1,:) &
     +fx*(rs(2,:)-6.0_dp*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-3,nyf-1,:))+15.0_dp*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:))+u(nxf-4,nyf-1,:)) &
     +fy*(rbufN(nxf-1,2,:)-6.0_dp*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:))+15.0_dp*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:))+u(nxf-1,nyf-4,:)))*deltaR(2)
    u(nxf,nyf-1,:)=Oomehd(qn)*u(nxf,nyf-1,:)+(rhs(nxf,nyf-1,:) &
!    u(nxf,nyf-1,:)=(rhs(nxf,nyf-1,:) &
     +fx*(Rq*u(nxf-2,nyf-1,:)+Rs(3,:)-6.0_dp*(Rq*u(nxf-1,nyf-1,:)+Rs(2,:)+u(nxf-2,nyf-1,:))+15.0_dp*(rs(1,:)+u(nxf-1,nyf-1,:))+u(nxf-3,nyf-1,:))  &
     +fy*(rbufN(nxf,2,:)-6.0_dp*(rbufN(nxf,1,:)+u(nxf,nyf-3,:))+15.0_dp*(u(nxf,nyf,:)+u(nxf,nyf-2,:))+u(nxf,nyf-4,:)))*deltaR(1)
  else
    u(nxf-2,nyf-1,:)=Oomehd(qn)*u(nxf-2,nyf-1,:)+(rhs(nxf-2,nyf-1,:) &
!    u(nxf-2,nyf-1,:)=(rhs(nxf-2,nyf-1,:) &
      +fx*(rbufE(1,nyf-1,:)-6.0_dp*(u(nxf,nyf-1,:)+u(nxf-4,nyf-1,:))+15.0_dp*(u(nxf-1,nyf-1,:)+u(nxf-3,nyf-1,:))+u(nxf-5,nyf-1,:)) &
      +fy*(rbufN(nxf-2,2,:)-6.0_dp*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-3,:))+15.0_dp*(u(nxf-2,nyf,:)+u(nxf-2,nyf-2,:))+u(nxf-2,nyf-4,:)))*delta
    u(nxf-1,nyf-1,:)=Oomehd(qn)*u(nxf-1,nyf-1,:)+(rhs(nxf-1,nyf-1,:) &
!    u(nxf-1,nyf-1,:)=(rhs(nxf-1,nyf-1,:) &
      +fx*(rbufE(2,nyf-1,:)-6.0_dp*(rbufE(1,nyf-1,:)+u(nxf-3,nyf-1,:))+15.0_dp*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:))+u(nxf-4,nyf-1,:)) &
      +fy*(rbufN(nxf-1,2,:)-6.0_dp*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:))+15.0_dp*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:))+u(nxf-1,nyf-4,:)))*delta
    u(nxf,nyf-1,:)=Oomehd(qn)*u(nxf,nyf-1,:)+(rhs(nxf,nyf-1,:) &
!    u(nxf,nyf-1,:)=(rhs(nxf,nyf-1,:) &
      +fx*(rbufE(3,nyf-1,:)-6.0_dp*(rbufE(2,nyf-1,:)+u(nxf-2,nyf-1,:))+15.0_dp*(rbufE(1,nyf-1,:)+u(nxf-1,nyf-1,:))+u(nxf-3,nyf-1,:)) &
      +fy*(rbufN(nxf,2,:)-6.0_dp*(rbufN(nxf,1,:)+u(nxf,nyf-3,:))+15.0_dp*(u(nxf,nyf,:)+u(nxf,nyf-2,:))+u(nxf,nyf-4,:)))*delta
  endif
!........................... TOP ROW ...........................!
  if (iIDhd(qn)%a(ng)==0) then !.Left most process (along x)
    u(1,nyf,:)=Oomehd(qn)*u(1,nyf,:)+(rhs(1,nyf,:) &
!    u(1,nyf,:)=(rhs(1,nyf,:) &
      +fx*(u(4,nyf,:)-6.0_dp*(u(3,nyf,:)+Lq*u(2,nyf,:)+Ls(2,:))+15.0_dp*(u(2,nyf,:)+Ls(1,:))+Lq*u(3,nyf,:)+Ls(3,:)) &
      +fy*(rbufN(1,3,:)-6.0_dp*(rbufN(1,2,:)+u(1,nyf-2,:))+15.0_dp*(rbufN(1,1,:)+u(1,nyf-1,:))+u(1,nyf-3,:)))*deltaL(1)
    u(2,nyf,:)=Oomehd(qn)*u(2,nyf,:)+(rhs(2,nyf,:) &
!    u(2,nyf,:)=(rhs(2,nyf,:) &
      +fx*(u(5,nyf,:)-6.0_dp*(u(4,nyf,:)+Lq*u(1,nyf,:)+Ls(1,:))+15.0_dp*(u(3,nyf,:)+u(1,nyf,:))+Ls(2,:)) &
      +fy*(rbufN(2,3,:)-6.0_dp*(rbufN(2,2,:)+u(2,nyf-2,:))+15.0_dp*(rbufN(2,1,:)+u(2,nyf-1,:))+u(2,nyf-3,:)))*deltaL(2)
    u(3,nyf,:)=Oomehd(qn)*u(3,nyf,:)+(rhs(3,nyf,:) &
!    u(3,nyf,:)=(rhs(3,nyf,:) &
      +fx*(u(6,nyf,:)-6.0_dp*(u(5,nyf,:)+u(1,nyf,:))+15.0_dp*(u(4,nyf,:)+u(2,nyf,:))+Lq*u(1,nyf,:)+Ls(1,:)) &
      +fy*(rbufN(3,3,:)-6.0_dp*(rbufN(3,2,:)+u(3,nyf-2,:))+15.0_dp*(rbufN(3,1,:)+u(3,nyf-1,:))+u(3,nyf-3,:)))*delta
  else
    u(1,nyf,:)=Oomehd(qn)*u(1,nyf,:)+(rhs(1,nyf,:) &
!    u(1,nyf,:)=(rhs(1,nyf,:) &
      +fx*(u(4,nyf,:)-6.0_dp*(u(3,nyf,:)+rbufW(2,nyf,:))+15.0_dp*(u(2,nyf,:)+rbufW(3,nyf,:))+rbufW(1,nyf,:)) &
      +fy*(rbufN(1,3,:)-6.0_dp*(rbufN(1,2,:)+u(1,nyf-2,:))+15.0_dp*(rbufN(1,1,:)+u(1,nyf-1,:))+u(1,nyf-3,:)))*delta
    u(2,nyf,:)=Oomehd(qn)*u(2,nyf,:)+(rhs(2,nyf,:) &
!    u(2,nyf,:)=(rhs(2,nyf,:) &
      +fx*(u(5,nyf,:)-6.0_dp*(u(4,nyf,:)+rbufW(3,nyf,:))+15.0_dp*(u(3,nyf,:)+u(1,nyf,:))+rbufW(2,nyf,:)) &
      +fy*(rbufN(2,3,:)-6.0_dp*(rbufN(2,2,:)+u(2,nyf-2,:))+15.0_dp*(rbufN(2,1,:)+u(2,nyf-1,:))+u(2,nyf-3,:)))*delta
    u(3,nyf,:)=Oomehd(qn)*u(3,nyf,:)+(rhs(3,nyf,:) &
!    u(3,nyf,:)=(rhs(3,nyf,:) &
      +fx*(u(6,nyf,:)-6.0_dp*(u(5,nyf,:)+u(1,nyf,:))+15.0_dp*(u(4,nyf,:)+u(2,nyf,:))+rbufW(3,nyf,:)) &
      +fy*(rbufN(3,3,:)-6.0_dp*(rbufN(3,2,:)+u(3,nyf-2,:))+15.0_dp*(rbufN(3,1,:)+u(3,nyf-1,:))+u(3,nyf-3,:)))*delta
  endif
  do i=4,nxf-3
    u(i,nyf,:)=Oomehd(qn)*u(i,nyf,:)+(rhs(i,nyf,:) &
!    u(i,nyf,:)=(rhs(i,nyf,:) &
      +fx*(u(i+3,nyf,:)-6.0_dp*(u(i+2,nyf,:)+u(i-2,nyf,:))+15.0_dp*(u(i+1,nyf,:)+u(i-1,nyf,:))+u(i-3,nyf,:)) &
      +fy*(rbufN(i,3,:)-6.0_dp*(rbufN(i,2,:)+u(i,nyf-2,:))+15.0_dp*(rbufN(i,1,:)+u(i,nyf-1,:))+u(i,nyf-3,:)))*delta
  enddo
  if (iIDhd(qn)%a(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf-2,nyf,:)=Oomehd(qn)*u(nxf-2,nyf,:)+(rhs(nxf-2,nyf,:) &
!    u(nxf-2,nyf,:)=(rhs(nxf-2,nyf,:) &
      +fx*(Rq*u(nxf,nyf,:)+Rs(1,:)-6.0_dp*(u(nxf,nyf,:)+u(nxf-4,nyf,:))+15.0_dp*(u(nxf-1,nyf,:)+u(nxf-3,nyf,:))+u(nxf-5,nyf,:)) &
      +fy*(rbufN(nxf-2,3,:)-6.0_dp*(rbufN(nxf-2,2,:)+u(nxf-2,nyf-2,:))+15.0_dp*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-1,:))+u(nxf-2,nyf-3,:)))*delta
    u(nxf-1,nyf,:)=Oomehd(qn)*u(nxf-1,nyf,:)+(rhs(nxf-1,nyf,:) &
!    u(nxf-1,nyf,:)=(rhs(nxf-1,nyf,:) &
      +fx*(rs(2,:)-6.0_dp*(Rq*u(nxf,nyf,:)+Rs(1,:)+u(nxf-3,nyf,:))+15.0_dp*(u(nxf,nyf,:)+u(nxf-2,nyf,:))+u(nxf-4,nyf,:)) &
      +fy*(rbufN(nxf-1,3,:)-6.0_dp*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:))+15.0_dp*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:))+u(nxf-1,nyf-3,:)))*deltaR(2)
    u(nxf,nyf,:)=Oomehd(qn)*u(nxf,nyf,:)+(rhs(nxf,nyf,:) &
!    u(nxf,nyf,:)=(rhs(nxf,nyf,:) &
      +fx*(Rq*u(nxf-2,nyf,:)+Rs(3,:)-6.0_dp*(Rq*u(nxf-1,nyf,:)+Rs(2,:)+u(nxf-2,nyf,:))+15.0_dp*(rs(1,:)+u(nxf-1,nyf,:))+u(nxf-3,nyf,:)) &
      +fy*(rbufN(nxf,3,:)-6.0_dp*(rbufN(nxf,2,:)+u(nxf,nyf-2,:))+15.0_dp*(rbufN(nxf,1,:)+u(nxf,nyf-1,:))+u(nxf,nyf-3,:)))*deltaR(1)
  else
    u(nxf-2,nyf,:)=Oomehd(qn)*u(nxf-2,nyf,:)+(rhs(nxf-2,nyf,:) &
!    u(nxf-2,nyf,:)=(rhs(nxf-2,nyf,:) &
      +fx*(rbufE(1,nyf,:)-6.0_dp*(u(nxf,nyf,:)+u(nxf-4,nyf,:))+15.0_dp*(u(nxf-1,nyf,:)+u(nxf-3,nyf,:))+u(nxf-5,nyf,:)) &
      +fy*(rbufN(nxf-2,3,:)-6.0_dp*(rbufN(nxf-2,2,:)+u(nxf-2,nyf-2,:))+15.0_dp*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-1,:))+u(nxf-2,nyf-3,:)))*delta
    u(nxf-1,nyf,:)=Oomehd(qn)*u(nxf-1,nyf,:)+(rhs(nxf-1,nyf,:) &
!    u(nxf-1,nyf,:)=(rhs(nxf-1,nyf,:) &
      +fx*(rbufE(2,nyf,:)-6.0_dp*(rbufE(1,nyf,:)+u(nxf-3,nyf,:))+15.0_dp*(u(nxf,nyf,:)+u(nxf-2,nyf,:))+u(nxf-4,nyf,:)) &
      +fy*(rbufN(nxf-1,3,:)-6.0_dp*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:))+15.0_dp*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:))+u(nxf-1,nyf-3,:)))*delta
    u(nxf,nyf,:)=Oomehd(qn)*u(nxf,nyf,:)+(rhs(nxf,nyf,:) &
!    u(nxf,nyf,:)=(rhs(nxf,nyf,:) &
      +fx*(rbufE(3,nyf,:)-6.0_dp*(rbufE(2,nyf,:)+u(nxf-2,nyf,:))+15.0_dp*(rbufE(1,nyf,:)+u(nxf-1,nyf,:))+u(nxf-3,nyf,:)) &
      +fy*(rbufN(nxf,3,:)-6.0_dp*(rbufN(nxf,2,:)+u(nxf,nyf-2,:))+15.0_dp*(rbufN(nxf,1,:)+u(nxf,nyf-1,:))+u(nxf,nyf-3,:)))*delta
  endif
#ENDIF

  call MPI_WAITALL(4,(/req(1),req(2),req(3),req(4)/), &
    (/stat(:,1),stat(:,2),stat(:,3),stat(:,4)/),ierr)

  end subroutine relaxhd
!*********************************************************
  function reshd(ng,u,rhs,xBC,qn)
!.Returns minus the residual. Input quantities
!.are u, av, rhs, while the residual is returned in reshd.
!.reshd,u,rhs are (nxf,nyf) arrays
  implicit none
  include 'mpif.h'
  real(DP), allocatable :: reshd(:,:,:)
  real(DP), intent(in) :: u(:,:,:),rhs(:,:,:)
  integer(I4B), intent(in) :: ng,xBC(:),qn
  integer(I4B) :: i,j,k,Nrank,Srank,Erank,Wrank
  real(DP) :: fx,fy,fc,Lq,Rq
  real(DP), allocatable :: Ls(:,:),rs(:,:),loc2D(:,:)
  real(DP), allocatable, dimension(:,:,:) :: sbufW,rbufE,sbufE,rbufW, &
                                             sbufN,rbufS,sbufS,rbufN
  integer(I4B) :: req(8),stat(MPI_STATUS_SIZE,8),ierr

#IF (MGTIMER==1)
    call CPU_TIME(rest1)
#ENDIF

  Nrank=nborhd(qn)%a(ng,1)
  Erank=nborhd(qn)%a(ng,3)
  Srank=nborhd(qn)%a(ng,5)
  Wrank=nborhd(qn)%a(ng,7)

  allocate(reshd(nxf,nyf,nzL))
!.Boundary conditions are implemented with some amount of flexibility.
!.The right boundary for example, if u(nxf+1) is needed, it will be
!.written as
!.  u(nxf+1) = Rq*u(nxf) + Rr*u(nxf-1) + Rs
!.so that symmetry, extrapolation or Dirichlet conditions can be
!.easily imposed.
!......... BCs of u ................................!
  allocate(Ls(HDopd4,nzL),rs(HDopd4,nzL),loc2D(HDopd4,nzL))
!.LHS BC of u
  if (xBC(1)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Lq=-1.0_dp
    forall(i=1:HDopd4,k=1:nzL) loc2D(i,k)=sum(u(i,:,k))
    call MPI_ALLreduce(loc2D,ls,HDopd4*nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    ls=2.0_dp*ls/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC(1)==1 .OR. xBC(1)==-1) then ! symmetry
    Lq=dble(xBC(1))
    ls=0.0_dp
  endif
!.RHS BC of u
  if (xBC(2)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Rq=-1.0_dp
    forall(i=1:HDopd4,k=1:nzL) loc2D(i,k)=sum(u(nxf-i+1,:,k))
    call MPI_ALLreduce(loc2D,rs,HDopd4*nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMhd(qn)%a(ng,1),ierr)
    rs=2.0_dp*rs/dble(nyf*jprocshd(qn)%a(ng))
  else ! if (xBC(2)==1 .OR. xBC(2)==-1) then ! symmetry
    Rq=dble(xBC(2))
    rs=0.0_dp
  endif
!........................................................!
  fx=hDiff(1)*((dble(nxf*iprocshd(qn)%a(ng))/bLx)**HDopd2)
  fy=hDiff(2)*((dble(nyf*jprocshd(qn)%a(ng))/bLy)**HDopd2)

  allocate(sbufN(nxf,HDopd4,nzL),rbufS(nxf,HDopd4,nzL),sbufE(HDopd4,nyf,nzL),rbufW(HDopd4,nyf,nzL))
  allocate(sbufS(nxf,HDopd4,nzL),rbufN(nxf,HDopd4,nzL),sbufW(HDopd4,nyf,nzL),rbufE(HDopd4,nyf,nzL))
!.Send top 2 rows of u to North rank
  sbufN=u(:,nyf-HDopd4+1:nyf,:)
  call MPI_ISEND(sbufN,nxf*HDopd4*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(1),ierr)
!.Send bottom 2 rows of u to South rank
  sbufS=u(:,1:HDopd4,:)
  call MPI_ISEND(sbufS,nxf*HDopd4*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(2),ierr)
!.Send right 2 columns of u to East rank
  sbufE=u(nxf-HDopd4+1:nxf,:,:)
  call MPI_ISEND(sbufE,HDopd4*nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(3),ierr)
!.Send 2 left columns of u to West rank
  sbufW=u(1:HDopd4,:,:)
  call MPI_ISEND(sbufW,HDopd4*nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(4),ierr)
!.Receive u rows from South rank, for u(i,j-2) and u(i,j-1)
  call MPI_IRECV(rbufS,nxf*HDopd4*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMhd(qn)%a(ng,1),req(5),ierr)
!.Receive u rows from North rank, for u(i,j+2) and u(i,j+1)
  call MPI_IRECV(rbufN,nxf*HDopd4*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMhd(qn)%a(ng,1),req(6),ierr)
!.Receive u columns from West rank, for u(i-2,j) and u(i-1,j)
  call MPI_IRECV(rbufW,HDopd4*nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMhd(qn)%a(ng,2),req(7),ierr)
!.Receive u columns from East rank, for u(i+2,j) and u(i+1,j)
  call MPI_IRECV(rbufE,HDopd4*nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMhd(qn)%a(ng,2),req(8),ierr)

#IF (HDop==8)

  fc=6.0_dp*(fx+fy)

!.Interior points
  do j=3,nyf-2
    do i=3,nxf-2
      reshd(i,j,:)=rhs(i,j,:)-((1.0_dp+fc)*u(i,j,:) &
        +fx*(u(i+2,j,:)-4.0_dp*(u(i+1,j,:)+u(i-1,j,:))+u(i-2,j,:)) &
        +fy*(u(i,j+2,:)-4.0_dp*(u(i,j+1,:)+u(i,j-1,:))+u(i,j-2,:)))
    enddo
  enddo

!.Boundary points
  call MPI_WAIT(req(5),stat(:,5),ierr)
  do i=3,nxf-2 ! Bottom boundary
    reshd(i,1,:)=rhs(i,1,:)-((1.0_dp+fc)*u(i,1,:) &
      +fx*(u(i+2,1,:)-4.0_dp*(u(i+1,1,:)+u(i-1,1,:))+u(i-2,1,:)) &
      +fy*(u(i,3,:)-4.0_dp*(u(i,2,:)+rbufS(i,2,:))+rbufS(i,1,:)))
    reshd(i,2,:)=rhs(i,2,:)-((1.0_dp+fc)*u(i,2,:) &
      +fx*(u(i+2,2,:)-4.0_dp*(u(i+1,2,:)+u(i-1,2,:))+u(i-2,2,:)) &
      +fy*(u(i,4,:)-4.0_dp*(u(i,3,:)+u(i,1,:))+rbufS(i,2,:)))
  enddo
  call MPI_WAIT(req(6),stat(:,6),ierr)
  do i=3,nxf-2 ! Top boundary
    reshd(i,nyf-1,:)=rhs(i,nyf-1,:)-((1.0_dp+fc)*u(i,nyf-1,:) &
      +fx*(u(i+2,nyf-1,:)-4.0_dp*(u(i+1,nyf-1,:)+u(i-1,nyf-1,:))+u(i-2,nyf-1,:)) &
      +fy*(rbufN(i,1,:)-4.0_dp*(u(i,nyf,:)+u(i,nyf-2,:))+u(i,nyf-3,:)))
    reshd(i,nyf,:)=rhs(i,nyf,:)-((1.0_dp+fc)*u(i,nyf,:) &
      +fx*(u(i+2,nyf,:)-4.0_dp*(u(i+1,nyf,:)+u(i-1,nyf,:))+u(i-2,nyf,:)) &
      +fy*(rbufN(i,2,:)-4.0_dp*(rbufN(i,1,:)+u(i,nyf-1,:))+u(i,nyf-2,:)))
  enddo
!.Left boundary
  call MPI_WAIT(req(7),stat(:,7),ierr)
  if (iIDhd(qn)%a(ng) == 0) then !.Left most process (along x)
!...SW corner
    reshd(1,1,:)=rhs(1,1,:)-((1.0_dp+fc)*u(1,1,:) &
      +fx*(u(3,1,:)-4.0_dp*(u(2,1,:)+Lq*u(1,1,:)+Ls(1,:))+Lq*u(2,1,:)+Ls(2,:)) &
      +fy*(u(1,3,:)-4.0_dp*(u(1,2,:)+rbufS(1,2,:))+rbufS(1,1,:)))
    reshd(2,1,:)=rhs(2,1,:)-((1.0_dp+fc)*u(2,1,:) &
      +fx*(u(4,1,:)-4.0_dp*(u(3,1,:)+u(1,1,:))+Lq*u(1,1,:)+Ls(1,:)) &
      +fy*(u(2,3,:)-4.0_dp*(u(2,2,:)+rbufS(2,2,:))+rbufS(2,1,:)))
    reshd(1,2,:)=rhs(1,2,:)-((1.0_dp+fc)*u(1,2,:) &
      +fx*(u(3,2,:)-4.0_dp*(u(2,2,:)+Lq*u(1,2,:)+Ls(1,:))+Lq*u(2,2,:)+Ls(2,:)) &
      +fy*(u(1,4,:)-4.0_dp*(u(1,3,:)+u(1,1,:))+rbufS(1,2,:)))
    reshd(2,2,:)=rhs(2,2,:)-((1.0_dp+fc)*u(2,2,:) &
      +fx*(u(4,2,:)-4.0_dp*(u(3,2,:)+u(1,2,:))+Lq*u(1,2,:)+Ls(1,:)) &
      +fy*(u(2,4,:)-4.0_dp*(u(2,3,:)+u(2,1,:))+rbufS(2,2,:)))
    do j=3,nyf-2
      reshd(1,j,:)=rhs(1,j,:)-((1.0_dp+fc)*u(1,j,:) &
        +fx*(u(3,j,:)-4.0_dp*(u(2,j,:)+Lq*u(1,j,:)+Ls(1,:))+Lq*u(2,j,:)+Ls(2,:)) &
        +fy*(u(1,j+2,:)-4.0_dp*(u(1,j+1,:)+u(1,j-1,:))+u(1,j-2,:)))
      reshd(2,j,:)=rhs(2,j,:)-((1.0_dp+fc)*u(2,j,:) &
        +fx*(u(4,j,:)-4.0_dp*(u(3,j,:)+u(1,j,:))+Lq*u(1,j,:)+Ls(1,:)) &
        +fy*(u(2,j+2,:)-4.0_dp*(u(2,j+1,:)+u(2,j-1,:))+u(2,j-2,:)))
    enddo
!...NW corner
    reshd(1,nyf-1,:)=rhs(1,nyf-1,:)-((1.0_dp+fc)*u(1,nyf-1,:) &
      +fx*(u(3,nyf-1,:)-4.0_dp*(u(2,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:))+Lq*u(2,nyf-1,:)+Ls(2,:)) &
      +fy*(rbufN(1,1,:)-4.0_dp*(u(1,nyf,:)+u(1,nyf-2,:))+u(1,nyf-3,:)))
    reshd(2,nyf-1,:)=rhs(2,nyf-1,:)-((1.0_dp+fc)*u(2,nyf-1,:) &
      +fx*(u(4,nyf-1,:)-4.0_dp*(u(3,nyf-1,:)+u(1,nyf-1,:))+Lq*u(1,nyf-1,:)+Ls(1,:)) &
      +fy*(rbufN(2,1,:)-4.0_dp*(u(2,nyf,:)+u(2,nyf-2,:))+u(2,nyf-3,:)))
    reshd(1,nyf,:)=rhs(1,nyf,:)-((1.0_dp+fc)*u(1,nyf,:) &
      +fx*(u(3,nyf,:)-4.0_dp*(u(2,nyf,:)+Lq*u(1,nyf,:)+Ls(1,:))+Lq*u(2,nyf,:)+Ls(2,:)) &
      +fy*(rbufN(1,2,:)-4.0_dp*(rbufN(1,1,:)+u(1,nyf-1,:))+u(1,nyf-2,:)))
    reshd(2,nyf,:)=rhs(2,nyf,:)-((1.0_dp+fc)*u(2,nyf,:) &
      +fx*(u(4,nyf,:)-4.0_dp*(u(3,nyf,:)+u(1,nyf,:))+Lq*u(1,nyf,:)+Ls(1,:)) &
      +fy*(rbufN(2,2,:)-4.0_dp*(rbufN(2,1,:)+u(2,nyf-1,:))+u(2,nyf-2,:)))
  else
    reshd(1,1,:)=rhs(1,1,:)-((1.0_dp+fc)*u(1,1,:) &
      +fx*(u(3,1,:)-4.0_dp*(u(2,1,:)+rbufW(2,1,:))+rbufW(1,1,:)) &
      +fy*(u(1,3,:)-4.0_dp*(u(1,2,:)+rbufS(1,2,:))+rbufS(1,1,:)))
    reshd(2,1,:)=rhs(2,1,:)-((1.0_dp+fc)*u(2,1,:) &
      +fx*(u(4,1,:)-4.0_dp*(u(3,1,:)+u(1,1,:))+rbufW(2,1,:)) &
      +fy*(u(2,3,:)-4.0_dp*(u(2,2,:)+rbufS(2,2,:))+rbufS(2,1,:)))
    reshd(1,2,:)=rhs(1,2,:)-((1.0_dp+fc)*u(1,2,:) &
      +fx*(u(3,2,:)-4.0_dp*(u(2,2,:)+rbufW(2,2,:))+rbufW(1,2,:)) &
      +fy*(u(1,4,:)-4.0_dp*(u(1,3,:)+u(1,1,:))+rbufS(1,2,:)))
    reshd(2,2,:)=rhs(2,2,:)-((1.0_dp+fc)*u(2,2,:) &
      +fx*(u(4,2,:)-4.0_dp*(u(3,2,:)+u(1,2,:))+rbufW(2,2,:)) &
      +fy*(u(2,4,:)-4.0_dp*(u(2,3,:)+u(2,1,:))+rbufS(2,2,:)))
    do j=3,nyf-2
      reshd(1,j,:)=rhs(1,j,:)-((1.0_dp+fc)*u(1,j,:) &
        +fx*(u(3,j,:)-4.0_dp*(u(2,j,:)+rbufW(2,j,:))+rbufW(1,j,:)) &
        +fy*(u(1,j+2,:)-4.0_dp*(u(1,j+1,:)+u(1,j-1,:))+u(1,j-2,:)))
      reshd(2,j,:)=rhs(2,j,:)-((1.0_dp+fc)*u(2,j,:) &
        +fx*(u(4,j,:)-4.0_dp*(u(3,j,:)+u(1,j,:))+rbufW(2,j,:)) &
        +fy*(u(2,j+2,:)-4.0_dp*(u(2,j+1,:)+u(2,j-1,:))+u(2,j-2,:)))
    enddo
    reshd(1,nyf-1,:)=rhs(1,nyf-1,:)-((1.0_dp+fc)*u(1,nyf-1,:) &
      +fx*(u(3,nyf-1,:)-4.0_dp*(u(2,nyf-1,:)+rbufW(2,nyf-1,:))+rbufW(1,nyf-1,:)) &
      +fy*(rbufN(1,1,:)-4.0_dp*(u(1,nyf,:)+u(1,nyf-2,:))+u(1,nyf-3,:)))
    reshd(2,nyf-1,:)=rhs(2,nyf-1,:)-((1.0_dp+fc)*u(2,nyf-1,:) &
      +fx*(u(4,nyf-1,:)-4.0_dp*(u(3,nyf-1,:)+u(1,nyf-1,:))+rbufW(2,nyf-1,:)) &
      +fy*(rbufN(2,1,:)-4.0_dp*(u(2,nyf,:)+u(2,nyf-2,:))+u(2,nyf-3,:)))
    reshd(1,nyf,:)=rhs(1,nyf,:)-((1.0_dp+fc)*u(1,nyf,:) &
      +fx*(u(3,nyf,:)-4.0_dp*(u(2,nyf,:)+rbufW(2,nyf,:))+rbufW(1,nyf,:)) &
      +fy*(rbufN(1,2,:)-4.0_dp*(rbufN(1,1,:)+u(1,nyf-1,:))+u(1,nyf-2,:)))
    reshd(2,nyf,:)=rhs(2,nyf,:)-((1.0_dp+fc)*u(2,nyf,:) &
      +fx*(u(4,nyf,:)-4.0_dp*(u(3,nyf,:)+u(1,nyf,:))+rbufW(2,nyf,:)) &
      +fy*(rbufN(2,2,:)-4.0_dp*(rbufN(2,1,:)+u(2,nyf-1,:))+u(2,nyf-2,:)))
  endif
!.Right boundary
  call MPI_WAIT(req(8),stat(:,8),ierr)
  if (iIDhd(qn)%a(ng) == iprocsm1) then !.Right most process (along x)
!...SE corner
    reshd(nxf-1,1,:)=rhs(nxf-1,1,:)-((1.0_dp+fc)*u(nxf-1,1,:) &
      +fx*(Rq*u(nxf,1,:)+Rs(1,:)-4.0_dp*(u(nxf,1,:)+u(nxf-2,1,:))+u(nxf-3,1,:)) &
      +fy*(u(nxf-1,3,:)-4.0_dp*(u(nxf-1,2,:)+rbufS(nxf-1,2,:))+rbufS(nxf-1,1,:)))
    reshd(nxf,1,:)=rhs(nxf,1,:)-((1.0_dp+fc)*u(nxf,1,:) &
      +fx*(Rq*u(nxf-1,1,:)+Rs(2,:)-4.0_dp*(Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-1,1,:))+u(nxf-2,1,:)) &
      +fy*(u(nxf,3,:)-4.0_dp*(u(nxf,2,:)+rbufS(nxf,2,:))+rbufS(nxf,1,:)))
    reshd(nxf-1,2,:)=rhs(nxf-1,2,:)-((1.0_dp+fc)*u(nxf-1,2,:) &
      +fx*(Rq*u(nxf,2,:)+Rs(1,:)-4.0_dp*(u(nxf,2,:)+u(nxf-2,2,:))+u(nxf-3,2,:)) &
      +fy*(u(nxf-1,4,:)-4.0_dp*(u(nxf-1,3,:)+u(nxf-1,1,:))+rbufS(nxf-1,2,:)))
    reshd(nxf,2,:)=rhs(nxf,2,:)-((1.0_dp+fc)*u(nxf,2,:) &
      +fx*(Rq*u(nxf-1,2,:)+Rs(2,:)-4.0_dp*(Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-1,2,:))+u(nxf-2,2,:)) &
      +fy*(u(nxf,4,:)-4.0_dp*(u(nxf,3,:)+u(nxf,1,:))+rbufS(nxf,2,:)))
    do j=3,nyf-2
      reshd(nxf-1,j,:)=rhs(nxf-1,j,:)-((1.0_dp+fc)*u(nxf-1,j,:) &
        +fx*(Rq*u(nxf,j,:)+Rs(1,:)-4.0_dp*(u(nxf,j,:)+u(nxf-2,j,:))+u(nxf-3,j,:)) &
        +fy*(u(nxf-1,j+2,:)-4.0_dp*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:))+u(nxf-1,j-2,:)))
      reshd(nxf,j,:)=rhs(nxf,j,:)-((1.0_dp+fc)*u(nxf,j,:) &
        +fx*(Rq*u(nxf-1,j,:)+Rs(2,:)-4.0_dp*(Rq*u(nxf,j,:)+Rs(1,:)+u(nxf-1,j,:))+u(nxf-2,j,:)) &
        +fy*(u(nxf,j+2,:)-4.0_dp*(u(nxf,j+1,:)+u(nxf,j-1,:))+u(nxf,j-2,:)))
    enddo
!...NE corner
    reshd(nxf-1,nyf-1,:)=rhs(nxf-1,nyf-1,:)-((1.0_dp+fc)*u(nxf-1,nyf-1,:) &
      +fx*(Rq*u(nxf,nyf-1,:)+Rs(1,:)-4.0_dp*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:))+u(nxf-3,nyf-1,:)) &
      +fy*(rbufN(nxf-1,1,:)-4.0_dp*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:))+u(nxf-1,nyf-3,:)))
    reshd(nxf,nyf-1,:)=rhs(nxf,nyf-1,:)-((1.0_dp+fc)*u(nxf,nyf-1,:) &
      +fx*(Rq*u(nxf-1,nyf-1,:)+Rs(2,:)-4.0_dp*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-1,nyf-1,:))+u(nxf-2,nyf-1,:)) &
      +fy*(rbufN(nxf,1,:)-4.0_dp*(u(nxf,nyf,:)+u(nxf,nyf-2,:))+u(nxf,nyf-3,:)))
    reshd(nxf-1,nyf,:)=rhs(nxf-1,nyf,:)-((1.0_dp+fc)*u(nxf-1,nyf,:) &
      +fx*(Rq*u(nxf,nyf,:)+Rs(1,:)-4.0_dp*(u(nxf,nyf,:)+u(nxf-2,nyf,:))+u(nxf-3,nyf,:)) &
      +fy*(rbufN(nxf-1,2,:)-4.0_dp*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:))+u(nxf-1,nyf-2,:)))
    reshd(nxf,nyf,:)=rhs(nxf,nyf,:)-((1.0_dp+fc)*u(nxf,nyf,:) &
      +fx*(Rq*u(nxf-1,nyf,:)+Rs(2,:)-4.0_dp*(Rq*u(nxf,nyf,:)+Rs(1,:)+u(nxf-1,nyf,:))+u(nxf-2,nyf,:)) &
      +fy*(rbufN(nxf,2,:)-4.0_dp*(rbufN(nxf,1,:)+u(nxf,nyf-1,:))+u(nxf,nyf-2,:)))
  else
    reshd(nxf-1,1,:)=rhs(nxf-1,1,:)-((1.0_dp+fc)*u(nxf-1,1,:) &
      +fx*(rbufE(1,1,:)-4.0_dp*(u(nxf,1,:)+u(nxf-2,1,:))+u(nxf-3,1,:)) &
      +fy*(u(nxf-1,3,:)-4.0_dp*(u(nxf-1,2,:)+rbufS(nxf-1,2,:))+rbufS(nxf-1,1,:)))
    reshd(nxf,1,:)=rhs(nxf,1,:)-((1.0_dp+fc)*u(nxf,1,:) &
      +fx*(rbufE(2,1,:)-4.0_dp*(rbufE(1,1,:)+u(nxf-1,1,:))+u(nxf-2,1,:)) &
      +fy*(u(nxf,3,:)-4.0_dp*(u(nxf,2,:)+rbufS(nxf,2,:))+rbufS(nxf,1,:)))
    reshd(nxf-1,2,:)=rhs(nxf-1,2,:)-((1.0_dp+fc)*u(nxf-1,2,:) &
      +fx*(rbufE(1,2,:)-4.0_dp*(u(nxf,2,:)+u(nxf-2,2,:))+u(nxf-3,2,:)) &
      +fy*(u(nxf-1,4,:)-4.0_dp*(u(nxf-1,3,:)+u(nxf-1,1,:))+rbufS(nxf-1,2,:)))
    reshd(nxf,2,:)=rhs(nxf,2,:)-((1.0_dp+fc)*u(nxf,2,:) &
      +fx*(rbufE(2,2,:)-4.0_dp*(rbufE(1,2,:)+u(nxf-1,2,:))+u(nxf-2,2,:)) &
      +fy*(u(nxf,4,:)-4.0_dp*(u(nxf,3,:)+u(nxf,1,:))+rbufS(nxf,2,:)))
    do j=3,nyf-2
      reshd(nxf-1,j,:)=rhs(nxf-1,j,:)-((1.0_dp+fc)*u(nxf-1,j,:) &
        +fx*(rbufE(1,j,:)-4.0_dp*(u(nxf,j,:)+u(nxf-2,j,:))+u(nxf-3,j,:)) &
        +fy*(u(nxf-1,j+2,:)-4.0_dp*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:))+u(nxf-1,j-2,:)))
      reshd(nxf,j,:)=rhs(nxf,j,:)-((1.0_dp+fc)*u(nxf,j,:) &
        +fx*(rbufE(2,j,:)-4.0_dp*(rbufE(1,j,:)+u(nxf-1,j,:))+u(nxf-2,j,:)) &
        +fy*(u(nxf,j+2,:)-4.0_dp*(u(nxf,j+1,:)+u(nxf,j-1,:))+u(nxf,j-2,:)))
    enddo
    reshd(nxf-1,nyf-1,:)=rhs(nxf-1,nyf-1,:)-((1.0_dp+fc)*u(nxf-1,nyf-1,:) &
      +fx*(rbufE(1,nyf-1,:)-4.0_dp*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:))+u(nxf-3,nyf-1,:)) &
      +fy*(rbufN(nxf-1,1,:)-4.0_dp*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:))+u(nxf-1,nyf-3,:)))
    reshd(nxf,nyf-1,:)=rhs(nxf,nyf-1,:)-((1.0_dp+fc)*u(nxf,nyf-1,:) &
      +fx*(rbufE(2,nyf-1,:)-4.0_dp*(rbufE(1,nyf-1,:)+u(nxf-1,nyf-1,:))+u(nxf-2,nyf-1,:)) &
      +fy*(rbufN(nxf,1,:)-4.0_dp*(u(nxf,nyf,:)+u(nxf,nyf-2,:))+u(nxf,nyf-3,:)))
    reshd(nxf-1,nyf,:)=rhs(nxf-1,nyf,:)-((1.0_dp+fc)*u(nxf-1,nyf,:) &
      +fx*(rbufE(1,nyf,:)-4.0_dp*(u(nxf,nyf,:)+u(nxf-2,nyf,:))+u(nxf-3,nyf,:)) &
      +fy*(rbufN(nxf-1,2,:)-4.0_dp*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:))+u(nxf-1,nyf-2,:)))
    reshd(nxf,nyf,:)=rhs(nxf,nyf,:)-((1.0_dp+fc)*u(nxf,nyf,:) &
      +fx*(rbufE(2,nyf,:)-4.0_dp*(rbufE(1,nyf,:)+u(nxf-1,nyf,:))+u(nxf-2,nyf,:)) &
      +fy*(rbufN(nxf,2,:)-4.0_dp*(rbufN(nxf,1,:)+u(nxf,nyf-1,:))+u(nxf,nyf-2,:)))
  endif

#ELIF (HDop==12)

  fc=1.0_dp+20.0_dp*(fx+fy)

!.Interior points
  do j=4,nyf-3
    do i=4,nxf-3
      reshd(i,j,:)=rhs(i,j,:)-(fc*u(i,j,:) &
        -fx*(u(i+3,j,:)-6.0_dp*(u(i+2,j,:)+u(i-2,j,:))+15.0_dp*(u(i+1,j,:)+u(i-1,j,:))+u(i-3,j,:)) &
        -fy*(u(i,j+3,:)-6.0_dp*(u(i,j+2,:)+u(i,j-2,:))+15.0_dp*(u(i,j+1,:)+u(i,j-1,:))+u(i,j-3,:)))
    enddo
  enddo

!.Boundary points
  call MPI_WAIT(req(5),stat(:,5),ierr)
  do i=4,nxf-3 ! Bottom boundary
    reshd(i,1,:)=rhs(i,1,:)-(fc*u(i,1,:) &
      -fx*(u(i+3,1,:)-6.0_dp*(u(i+2,1,:)+u(i-2,1,:))+15.0_dp*(u(i+1,1,:)+u(i-1,1,:))+u(i-3,1,:)) &
      -fy*(u(i,4,:)-6.0_dp*(u(i,3,:)+rbufS(i,2,:))+15.0_dp*(u(i,2,:)+rbufS(i,3,:))+rbufS(i,1,:)))
    reshd(i,2,:)=rhs(i,2,:)-(fc*u(i,2,:) &
      -fx*(u(i+3,2,:)-6.0_dp*(u(i+2,2,:)+u(i-2,2,:))+15.0_dp*(u(i+1,2,:)+u(i-1,2,:))+u(i-3,2,:)) &
      -fy*(u(i,5,:)-6.0_dp*(u(i,4,:)+rbufS(i,3,:))+15.0_dp*(u(i,3,:)+u(i,1,:))+rbufS(i,2,:)))
    reshd(i,3,:)=rhs(i,3,:)-(fc*u(i,3,:) &
      -fx*(u(i+3,3,:)-6.0_dp*(u(i+2,3,:)+u(i-2,3,:))+15.0_dp*(u(i+1,3,:)+u(i-1,3,:))+u(i-3,3,:)) &
      -fy*(u(i,6,:)-6.0_dp*(u(i,5,:)+u(i,1,:))+15.0_dp*(u(i,4,:)+u(i,2,:))+rbufS(i,3,:)))
  enddo
  call MPI_WAIT(req(6),stat(:,6),ierr)
  do i=4,nxf-3 ! Top boundary
    reshd(i,nyf-2,:)=rhs(i,nyf-2,:)-(fc*u(i,nyf-2,:) &
      -fx*(u(i+3,nyf-2,:)-6.0_dp*(u(i+2,nyf-2,:)+u(i-2,nyf-2,:))+15.0_dp*(u(i+1,nyf-2,:)+u(i-1,nyf-2,:))+u(i-3,nyf-2,:)) &
      -fy*(rbufN(i,1,:)-6.0_dp*(u(i,nyf,:)+u(i,nyf-4,:))+15.0_dp*(u(i,nyf-1,:)+u(i,nyf-3,:))+u(i,nyf-5,:)))
    reshd(i,nyf-1,:)=rhs(i,nyf-1,:)-(fc*u(i,nyf-1,:) &
      -fx*(u(i+3,nyf-1,:)-6.0_dp*(u(i+2,nyf-1,:)+u(i-2,nyf-1,:))+15.0_dp*(u(i+1,nyf-1,:)+u(i-1,nyf-1,:))+u(i-3,nyf-1,:)) &
      -fy*(rbufN(i,2,:)-6.0_dp*(rbufN(i,1,:)+u(i,nyf-3,:))+15.0_dp*(u(i,nyf,:)+u(i,nyf-2,:))+u(i,nyf-4,:)))
    reshd(i,nyf,:)=rhs(i,nyf,:)-(fc*u(i,nyf,:) &
      -fx*(u(i+3,nyf,:)-6.0_dp*(u(i+2,nyf,:)+u(i-2,nyf,:))+15.0_dp*(u(i+1,nyf,:)+u(i-1,nyf,:))+u(i-3,nyf,:)) &
      -fy*(rbufN(i,3,:)-6.0_dp*(rbufN(i,2,:)+u(i,nyf-2,:))+15.0_dp*(rbufN(i,1,:)+u(i,nyf-1,:))+u(i,nyf-3,:)))
  enddo
!.Left boundary
  call MPI_WAIT(req(7),stat(:,7),ierr)
  if (iIDhd(qn)%a(ng) == 0) then !.Left most process (along x)
!...SW corner
    reshd(1,1,:)=rhs(1,1,:)-(fc*u(1,1,:) &
      -fx*(u(4,1,:)-6.0_dp*(u(3,1,:)+Lq*u(2,1,:)+Ls(2,:))+15.0_dp*(u(2,1,:)+Lq*u(1,1,:)+Ls(1,:))+Lq*u(3,1,:)+Ls(3,:)) &
      -fy*(u(1,4,:)-6.0_dp*(u(1,3,:)+rbufS(1,2,:))+15.0_dp*(u(1,2,:)+rbufS(1,3,:))+rbufS(1,1,:)))
    reshd(2,1,:)=rhs(2,1,:)-(fc*u(2,1,:) &
      -fx*(u(5,1,:)-6.0_dp*(u(4,1,:)+Lq*u(1,1,:)+Ls(1,:))+15.0_dp*(u(3,1,:)+u(1,1,:))+Lq*u(2,1,:)+Ls(2,:)) &
      -fy*(u(2,4,:)-6.0_dp*(u(2,3,:)+rbufS(2,2,:))+15.0_dp*(u(2,2,:)+rbufS(2,3,:))+rbufS(2,1,:)))
    reshd(3,1,:)=rhs(3,1,:)-(fc*u(3,1,:) &
      -fx*(u(6,1,:)-6.0_dp*(u(5,1,:)+u(1,1,:))+15.0_dp*(u(4,1,:)+u(2,1,:))+Lq*u(1,1,:)+Ls(1,:)) &
      -fy*(u(3,4,:)-6.0_dp*(u(3,3,:)+rbufS(3,2,:))+15.0_dp*(u(3,2,:)+rbufS(3,3,:))+rbufS(3,1,:)))
    reshd(1,2,:)=rhs(1,2,:)-(fc*u(1,2,:) &
      -fx*(u(4,2,:)-6.0_dp*(u(3,2,:)+Lq*u(2,2,:)+Ls(2,:))+15.0_dp*(u(2,2,:)+Lq*u(1,2,:)+Ls(1,:))+Lq*u(3,2,:)+Ls(3,:)) &
      -fy*(u(1,5,:)-6.0_dp*(u(1,4,:)+rbufS(1,3,:))+15.0_dp*(u(1,3,:)+u(1,1,:))+rbufS(1,2,:)))
    reshd(2,2,:)=rhs(2,2,:)-(fc*u(2,2,:) &
      -fx*(u(5,2,:)-6.0_dp*(u(4,2,:)+Lq*u(1,2,:)+Ls(1,:))+15.0_dp*(u(3,2,:)+u(1,2,:))+Lq*u(2,2,:)+Ls(2,:)) &
      -fy*(u(2,5,:)-6.0_dp*(u(2,4,:)+rbufS(2,3,:))+15.0_dp*(u(2,3,:)+u(2,1,:))+rbufS(2,2,:)))
    reshd(3,2,:)=rhs(3,2,:)-(fc*u(3,2,:) &
      -fx*(u(6,2,:)-6.0_dp*(u(5,2,:)+u(1,2,:))+15.0_dp*(u(4,2,:)+u(2,2,:))+Lq*u(1,2,:)+Ls(1,:)) &
      -fy*(u(3,5,:)-6.0_dp*(u(3,4,:)+rbufS(3,3,:))+15.0_dp*(u(3,3,:)+u(3,1,:))+rbufS(3,2,:)))
    reshd(1,3,:)=rhs(1,3,:)-(fc*u(1,3,:) &
      -fx*(u(4,3,:)-6.0_dp*(u(3,3,:)+Lq*u(2,3,:)+Ls(2,:))+15.0_dp*(u(2,3,:)+Lq*u(1,3,:)+Ls(1,:))+Lq*u(3,3,:)+Ls(3,:)) &
      -fy*(u(1,6,:)-6.0_dp*(u(1,5,:)+u(1,1,:))+15.0_dp*(u(1,4,:)+u(1,2,:))+rbufS(1,3,:)))
    reshd(2,3,:)=rhs(2,3,:)-(fc*u(2,3,:) &
      -fx*(u(5,3,:)-6.0_dp*(u(4,3,:)+Lq*u(1,3,:)+Ls(1,:))+15.0_dp*(u(3,3,:)+u(1,3,:))+Lq*u(2,3,:)+Ls(2,:)) &
      -fy*(u(2,6,:)-6.0_dp*(u(2,5,:)+u(2,1,:))+15.0_dp*(u(2,4,:)+u(2,2,:))+rbufS(2,3,:)))
    reshd(3,3,:)=rhs(3,3,:)-(fc*u(3,3,:) &
      -fx*(u(6,3,:)-6.0_dp*(u(5,3,:)+u(1,3,:))+15.0_dp*(u(4,3,:)+u(2,3,:))+Lq*u(1,3,:)+Ls(1,:)) &
      -fy*(u(3,6,:)-6.0_dp*(u(3,5,:)+u(3,1,:))+15.0_dp*(u(3,4,:)+u(3,2,:))+rbufS(3,3,:)))
    do j=4,nyf-3 ! Left boundary
      reshd(1,j,:)=rhs(1,j,:)-(fc*u(1,j,:) &
        -fx*(u(4,j,:)-6.0_dp*(u(3,j,:)+Lq*u(2,j,:)+Ls(2,:))+15.0_dp*(u(2,j,:)+Lq*u(1,j,:)+Ls(1,:))+Lq*u(3,j,:)+Ls(3,:)) &
        -fy*(u(1,j+3,:)-6.0_dp*(u(1,j+2,:)+u(1,j-2,:))+15.0_dp*(u(1,j+1,:)+u(1,j-1,:))+u(1,j-3,:)))
      reshd(2,j,:)=rhs(2,j,:)-(fc*u(2,j,:) &
        -fx*(u(5,j,:)-6.0_dp*(u(4,j,:)+Lq*u(1,j,:)+Ls(1,:))+15.0_dp*(u(3,j,:)+u(1,j,:))+Lq*u(2,j,:)+Ls(2,:)) &
        -fy*(u(2,j+3,:)-6.0_dp*(u(2,j+2,:)+u(2,j-2,:))+15.0_dp*(u(2,j+1,:)+u(2,j-1,:))+u(2,j-3,:)))
      reshd(3,j,:)=rhs(3,j,:)-(fc*u(3,j,:) &
        -fx*(u(6,j,:)-6.0_dp*(u(5,j,:)+u(1,j,:))+15.0_dp*(u(4,j,:)+u(2,j,:))+Lq*u(1,j,:)+Ls(1,:)) &
        -fy*(u(3,j+3,:)-6.0_dp*(u(3,j+2,:)+u(3,j-2,:))+15.0_dp*(u(3,j+1,:)+u(3,j-1,:))+u(3,j-3,:)))
    enddo
!...NW corner
    reshd(1,nyf-2,:)=rhs(1,nyf-2,:)-(fc*u(1,nyf-2,:) &
      -fx*(u(4,nyf-2,:)-6.0_dp*(u(3,nyf-2,:)+Lq*u(2,nyf-2,:)+Ls(2,:))+15.0_dp*(u(2,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:))+Lq*u(3,nyf-2,:)+Ls(3,:)) &
      -fy*(rbufN(1,1,:)-6.0_dp*(u(1,nyf,:)+u(1,nyf-4,:))+15.0_dp*(u(1,nyf-1,:)+u(1,nyf-3,:))+u(1,nyf-5,:)))
    reshd(2,nyf-2,:)=rhs(2,nyf-2,:)-(fc*u(2,nyf-2,:) &
      -fx*(u(5,nyf-2,:)-6.0_dp*(u(4,nyf-2,:)+Lq*u(1,nyf-2,:)+Ls(1,:))+15.0_dp*(u(3,nyf-2,:)+u(1,nyf-2,:))+Lq*u(2,nyf-2,:)+Ls(2,:)) &
      -fy*(rbufN(2,1,:)-6.0_dp*(u(2,nyf,:)+u(2,nyf-4,:))+15.0_dp*(u(2,nyf-1,:)+u(2,nyf-3,:))+u(2,nyf-5,:)))
    reshd(3,nyf-2,:)=rhs(3,nyf-2,:)-(fc*u(3,nyf-2,:) &
      -fx*(u(6,nyf-2,:)-6.0_dp*(u(5,nyf-2,:)+u(1,nyf-2,:))+15.0_dp*(u(4,nyf-2,:)+u(2,nyf-2,:))+Lq*u(1,nyf-2,:)+Ls(1,:)) &
      -fy*(rbufN(3,1,:)-6.0_dp*(u(3,nyf,:)+u(3,nyf-4,:))+15.0_dp*(u(3,nyf-1,:)+u(3,nyf-3,:))+u(3,nyf-5,:)))
    reshd(1,nyf-1,:)=rhs(1,nyf-1,:)-(fc*u(1,nyf-1,:) &
      -fx*(u(4,nyf-1,:)-6.0_dp*(u(3,nyf-1,:)+Lq*u(2,nyf-1,:)+Ls(2,:))+15.0_dp*(u(2,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:))+Lq*u(3,nyf-1,:)+Ls(3,:)) &
      -fy*(rbufN(1,2,:)-6.0_dp*(rbufN(1,1,:)+u(1,nyf-3,:))+15.0_dp*(u(1,nyf,:)+u(1,nyf-2,:))+u(1,nyf-4,:)))
    reshd(2,nyf-1,:)=rhs(2,nyf-1,:)-(fc*u(2,nyf-1,:) &
      -fx*(u(5,nyf-1,:)-6.0_dp*(u(4,nyf-1,:)+Lq*u(1,nyf-1,:)+Ls(1,:))+15.0_dp*(u(3,nyf-1,:)+u(1,nyf-1,:))+Lq*u(2,nyf-1,:)+Ls(2,:)) &
      -fy*(rbufN(2,2,:)-6.0_dp*(rbufN(2,1,:)+u(2,nyf-3,:))+15.0_dp*(u(2,nyf,:)+u(2,nyf-2,:))+u(2,nyf-4,:)))
    reshd(3,nyf-1,:)=rhs(3,nyf-1,:)-(fc*u(3,nyf-1,:) &
      -fx*(u(6,nyf-1,:)-6.0_dp*(u(5,nyf-1,:)+u(1,nyf-1,:))+15.0_dp*(u(4,nyf-1,:)+u(2,nyf-1,:))+Lq*u(1,nyf-1,:)+Ls(1,:)) &
      -fy*(rbufN(3,2,:)-6.0_dp*(rbufN(3,1,:)+u(3,nyf-3,:))+15.0_dp*(u(3,nyf,:)+u(3,nyf-2,:))+u(3,nyf-4,:)))
    reshd(1,nyf,:)=rhs(1,nyf,:)-(fc*u(1,nyf,:) &
      -fx*(u(4,nyf,:)-6.0_dp*(u(3,nyf,:)+Lq*u(2,nyf,:)+Ls(2,:))+15.0_dp*(u(2,nyf,:)+Lq*u(1,nyf,:)+Ls(1,:))+Lq*u(3,nyf,:)+Ls(3,:)) &
      -fy*(rbufN(1,3,:)-6.0_dp*(rbufN(1,2,:)+u(1,nyf-2,:))+15.0_dp*(rbufN(1,1,:)+u(1,nyf-1,:))+u(1,nyf-3,:)))
    reshd(2,nyf,:)=rhs(2,nyf,:)-(fc*u(2,nyf,:) &
      -fx*(u(5,nyf,:)-6.0_dp*(u(4,nyf,:)+Lq*u(1,nyf,:)+Ls(1,:))+15.0_dp*(u(3,nyf,:)+u(1,nyf,:))+Lq*u(2,nyf,:)+Ls(2,:)) &
      -fy*(rbufN(2,3,:)-6.0_dp*(rbufN(2,2,:)+u(2,nyf-2,:))+15.0_dp*(rbufN(2,1,:)+u(2,nyf-1,:))+u(2,nyf-3,:)))
    reshd(3,nyf,:)=rhs(3,nyf,:)-(fc*u(3,nyf,:) &
      -fx*(u(6,nyf,:)-6.0_dp*(u(5,nyf,:)+u(1,nyf,:))+15.0_dp*(u(4,nyf,:)+u(2,nyf,:))+Lq*u(1,nyf,:)+Ls(1,:)) &
      -fy*(rbufN(3,3,:)-6.0_dp*(rbufN(3,2,:)+u(3,nyf-2,:))+15.0_dp*(rbufN(3,1,:)+u(3,nyf-1,:))+u(3,nyf-3,:)))
  else
    reshd(1,1,:)=rhs(1,1,:)-(fc*u(1,1,:) &
      -fx*(u(4,1,:)-6.0_dp*(u(3,1,:)+rbufW(2,1,:))+15.0_dp*(u(2,1,:)+rbufW(3,1,:))+rbufW(1,1,:)) &
      -fy*(u(1,4,:)-6.0_dp*(u(1,3,:)+rbufS(1,2,:))+15.0_dp*(u(1,2,:)+rbufS(1,3,:))+rbufS(1,1,:)))
    reshd(2,1,:)=rhs(2,1,:)-(fc*u(2,1,:) &
      -fx*(u(5,1,:)-6.0_dp*(u(4,1,:)+rbufW(3,1,:))+15.0_dp*(u(3,1,:)+u(1,1,:))+rbufW(2,1,:)) &
      -fy*(u(2,4,:)-6.0_dp*(u(2,3,:)+rbufS(2,2,:))+15.0_dp*(u(2,2,:)+rbufS(2,3,:))+rbufS(2,1,:)))
    reshd(3,1,:)=rhs(3,1,:)-(fc*u(3,1,:) &
      -fx*(u(6,1,:)-6.0_dp*(u(5,1,:)+u(1,1,:))+15.0_dp*(u(4,1,:)+u(2,1,:))+rbufW(3,1,:)) &
      -fy*(u(3,4,:)-6.0_dp*(u(3,3,:)+rbufS(3,2,:))+15.0_dp*(u(3,2,:)+rbufS(3,3,:))+rbufS(3,1,:)))
    reshd(1,2,:)=rhs(1,2,:)-(fc*u(1,2,:) &
      -fx*(u(4,2,:)-6.0_dp*(u(3,2,:)+rbufW(2,2,:))+15.0_dp*(u(2,2,:)+rbufW(3,2,:))+rbufW(1,2,:)) &
      -fy*(u(1,5,:)-6.0_dp*(u(1,4,:)+rbufS(1,3,:))+15.0_dp*(u(1,3,:)+u(1,1,:))+rbufS(1,2,:)))
    reshd(2,2,:)=rhs(2,2,:)-(fc*u(2,2,:) &
      -fx*(u(5,2,:)-6.0_dp*(u(4,2,:)+rbufW(3,2,:))+15.0_dp*(u(3,2,:)+u(1,2,:))+rbufW(2,2,:)) &
      -fy*(u(2,5,:)-6.0_dp*(u(2,4,:)+rbufS(2,3,:))+15.0_dp*(u(2,3,:)+u(2,1,:))+rbufS(2,2,:)))
    reshd(3,2,:)=rhs(3,2,:)-(fc*u(3,2,:) &
      -fx*(u(6,2,:)-6.0_dp*(u(5,2,:)+u(1,2,:))+15.0_dp*(u(4,2,:)+u(2,2,:))+rbufW(3,2,:)) &
      -fy*(u(3,5,:)-6.0_dp*(u(3,4,:)+rbufS(3,3,:))+15.0_dp*(u(3,3,:)+u(3,1,:))+rbufS(3,2,:)))
    reshd(1,3,:)=rhs(1,3,:)-(fc*u(1,3,:) &
      -fx*(u(4,3,:)-6.0_dp*(u(3,3,:)+rbufW(2,3,:))+15.0_dp*(u(2,3,:)+rbufW(3,3,:))+rbufW(1,3,:)) &
      -fy*(u(1,6,:)-6.0_dp*(u(1,5,:)+u(1,1,:))+15.0_dp*(u(1,4,:)+u(1,2,:))+rbufS(1,3,:)))
    reshd(2,3,:)=rhs(2,3,:)-(fc*u(2,3,:) &
      -fx*(u(5,3,:)-6.0_dp*(u(4,3,:)+rbufW(3,3,:))+15.0_dp*(u(3,3,:)+u(1,3,:))+rbufW(2,3,:)) &
      -fy*(u(2,6,:)-6.0_dp*(u(2,5,:)+u(2,1,:))+15.0_dp*(u(2,4,:)+u(2,2,:))+rbufS(2,3,:)))
    reshd(3,3,:)=rhs(3,3,:)-(fc*u(3,3,:) &
      -fx*(u(6,3,:)-6.0_dp*(u(5,3,:)+u(1,3,:))+15.0_dp*(u(4,3,:)+u(2,3,:))+rbufW(3,3,:)) &
      -fy*(u(3,6,:)-6.0_dp*(u(3,5,:)+u(3,1,:))+15.0_dp*(u(3,4,:)+u(3,2,:))+rbufS(3,3,:)))
    do j=4,nyf-3
      reshd(1,j,:)=rhs(1,j,:)-(fc*u(1,j,:) &
        -fx*(u(4,j,:)-6.0_dp*(u(3,j,:)+rbufW(2,j,:))+15.0_dp*(u(2,j,:)+rbufW(3,j,:))+rbufW(1,j,:)) &
        -fy*(u(1,j+3,:)-6.0_dp*(u(1,j+2,:)+u(1,j-2,:))+15.0_dp*(u(1,j+1,:)+u(1,j-1,:))+u(1,j-3,:)))
      reshd(2,j,:)=rhs(2,j,:)-(fc*u(2,j,:) &
        -fx*(u(5,j,:)-6.0_dp*(u(4,j,:)+rbufW(3,j,:))+15.0_dp*(u(3,j,:)+u(1,j,:))+rbufW(2,j,:)) &
        -fy*(u(2,j+3,:)-6.0_dp*(u(2,j+2,:)+u(2,j-2,:))+15.0_dp*(u(2,j+1,:)+u(2,j-1,:))+u(2,j-3,:)))
      reshd(3,j,:)=rhs(3,j,:)-(fc*u(3,j,:) &
        -fx*(u(6,j,:)-6.0_dp*(u(5,j,:)+u(1,j,:))+15.0_dp*(u(4,j,:)+u(2,j,:))+rbufW(3,j,:)) &
        -fy*(u(3,j+3,:)-6.0_dp*(u(3,j+2,:)+u(3,j-2,:))+15.0_dp*(u(3,j+1,:)+u(3,j-1,:))+u(3,j-3,:)))
    enddo
    reshd(1,nyf-2,:)=rhs(1,nyf-2,:)-(fc*u(1,nyf-2,:) &
      -fx*(u(4,nyf-2,:)-6.0_dp*(u(3,nyf-2,:)+rbufW(2,nyf-2,:))+15.0_dp*(u(2,nyf-2,:)+rbufW(3,nyf-2,:))+rbufW(1,nyf-2,:)) &
      -fy*(rbufN(1,1,:)-6.0_dp*(u(1,nyf,:)+u(1,nyf-4,:))+15.0_dp*(u(1,nyf-1,:)+u(1,nyf-3,:))+u(1,nyf-5,:)))
    reshd(2,nyf-2,:)=rhs(2,nyf-2,:)-(fc*u(2,nyf-2,:) &
      -fx*(u(5,nyf-2,:)-6.0_dp*(u(4,nyf-2,:)+rbufW(3,nyf-2,:))+15.0_dp*(u(3,nyf-2,:)+u(1,nyf-2,:))+rbufW(2,nyf-2,:)) &
      -fy*(rbufN(2,1,:)-6.0_dp*(u(2,nyf,:)+u(2,nyf-4,:))+15.0_dp*(u(2,nyf-1,:)+u(2,nyf-3,:))+u(2,nyf-5,:)))
    reshd(3,nyf-2,:)=rhs(3,nyf-2,:)-(fc*u(3,nyf-2,:) &
      -fx*(u(6,nyf-2,:)-6.0_dp*(u(5,nyf-2,:)+u(1,nyf-2,:))+15.0_dp*(u(4,nyf-2,:)+u(2,nyf-2,:))+rbufW(3,nyf-2,:)) &
      -fy*(rbufN(3,1,:)-6.0_dp*(u(3,nyf,:)+u(3,nyf-4,:))+15.0_dp*(u(3,nyf-1,:)+u(3,nyf-3,:))+u(3,nyf-5,:)))
    reshd(1,nyf-1,:)=rhs(1,nyf-1,:)-(fc*u(1,nyf-1,:) &
      -fx*(u(4,nyf-1,:)-6.0_dp*(u(3,nyf-1,:)+rbufW(2,nyf-1,:))+15.0_dp*(u(2,nyf-1,:)+rbufW(3,nyf-1,:))+rbufW(1,nyf-1,:)) &
      -fy*(rbufN(1,2,:)-6.0_dp*(rbufN(1,1,:)+u(1,nyf-3,:))+15.0_dp*(u(1,nyf,:)+u(1,nyf-2,:))+u(1,nyf-4,:)))
    reshd(2,nyf-1,:)=rhs(2,nyf-1,:)-(fc*u(2,nyf-1,:) &
      -fx*(u(5,nyf-1,:)-6.0_dp*(u(4,nyf-1,:)+rbufW(3,nyf-1,:))+15.0_dp*(u(3,nyf-1,:)+u(1,nyf-1,:))+rbufW(2,nyf-1,:)) &
      -fy*(rbufN(2,2,:)-6.0_dp*(rbufN(2,1,:)+u(2,nyf-3,:))+15.0_dp*(u(2,nyf,:)+u(2,nyf-2,:))+u(2,nyf-4,:)))
    reshd(3,nyf-1,:)=rhs(3,nyf-1,:)-(fc*u(3,nyf-1,:) &
      -fx*(u(6,nyf-1,:)-6.0_dp*(u(5,nyf-1,:)+u(1,nyf-1,:))+15.0_dp*(u(4,nyf-1,:)+u(2,nyf-1,:))+rbufW(3,nyf-1,:)) &
      -fy*(rbufN(3,2,:)-6.0_dp*(rbufN(3,1,:)+u(3,nyf-3,:))+15.0_dp*(u(3,nyf,:)+u(3,nyf-2,:))+u(3,nyf-4,:)))
    reshd(1,nyf,:)=rhs(1,nyf,:)-(fc*u(1,nyf,:) &
      -fx*(u(4,nyf,:)-6.0_dp*(u(3,nyf,:)+rbufW(2,nyf,:))+15.0_dp*(u(2,nyf,:)+rbufW(3,nyf,:))+rbufW(1,nyf,:)) &
      -fy*(rbufN(1,3,:)-6.0_dp*(rbufN(1,2,:)+u(1,nyf-2,:))+15.0_dp*(rbufN(1,1,:)+u(1,nyf-1,:))+u(1,nyf-3,:)))
    reshd(2,nyf,:)=rhs(2,nyf,:)-(fc*u(2,nyf,:) &
      -fx*(u(5,nyf,:)-6.0_dp*(u(4,nyf,:)+rbufW(3,nyf,:))+15.0_dp*(u(3,nyf,:)+u(1,nyf,:))+rbufW(2,nyf,:)) &
      -fy*(rbufN(2,3,:)-6.0_dp*(rbufN(2,2,:)+u(2,nyf-2,:))+15.0_dp*(rbufN(2,1,:)+u(2,nyf-1,:))+u(2,nyf-3,:)))
    reshd(3,nyf,:)=rhs(3,nyf,:)-(fc*u(3,nyf,:) &
      -fx*(u(6,nyf,:)-6.0_dp*(u(5,nyf,:)+u(1,nyf,:))+15.0_dp*(u(4,nyf,:)+u(2,nyf,:))+rbufW(3,nyf,:)) &
      -fy*(rbufN(3,3,:)-6.0_dp*(rbufN(3,2,:)+u(3,nyf-2,:))+15.0_dp*(rbufN(3,1,:)+u(3,nyf-1,:))+u(3,nyf-3,:)))
  endif
!.Right boundary
  call MPI_WAIT(req(8),stat(:,8),ierr)
  if (iIDhd(qn)%a(ng) == iprocsm1) then !.Right most process (along x)
!...SE corner
    reshd(nxf-2,1,:)=rhs(nxf-2,1,:)-(fc*u(nxf-2,1,:) &
      -fx*(Rq*u(nxf,1,:)+Rs(1,:)-6.0_dp*(u(nxf,1,:)+u(nxf-4,1,:))+15.0_dp*(u(nxf-1,1,:)+u(nxf-3,1,:))+u(nxf-5,1,:)) &
      -fy*(u(nxf-2,4,:)-6.0_dp*(u(nxf-2,3,:)+rbufS(nxf-2,2,:))+15.0_dp*(u(nxf-2,2,:)+rbufS(nxf-2,3,:))+rbufS(nxf-2,1,:)))
    reshd(nxf-1,1,:)=rhs(nxf-1,1,:)-(fc*u(nxf-1,1,:) &
      -fx*(Rq*u(nxf-1,1,:)+Rs(2,:)-6.0_dp*(Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-3,1,:))+15.0_dp*(u(nxf,1,:)+u(nxf-2,1,:))+u(nxf-4,1,:)) &
      -fy*(u(nxf-1,4,:)-6.0_dp*(u(nxf-1,3,:)+rbufS(nxf-1,2,:))+15.0_dp*(u(nxf-1,2,:)+rbufS(nxf-1,3,:))+rbufS(nxf-1,1,:)))
    reshd(nxf,1,:)=rhs(nxf,1,:)-(fc*u(nxf,1,:) &
      -fx*(Rq*u(nxf-2,1,:)+Rs(3,:)-6.0_dp*(Rq*u(nxf-1,1,:)+Rs(2,:)+u(nxf-2,1,:))+15.0_dp*(Rq*u(nxf,1,:)+Rs(1,:)+u(nxf-1,1,:))+u(nxf-3,1,:)) &
      -fy*(u(nxf,4,:)-6.0_dp*(u(nxf,3,:)+rbufS(nxf,2,:))+15.0_dp*(u(nxf,2,:)+rbufS(nxf,3,:))+rbufS(nxf,1,:)))
    reshd(nxf-2,2,:)=rhs(nxf-2,2,:)-(fc*u(nxf-2,2,:) &
      -fx*(Rq*u(nxf,2,:)+Rs(1,:)-6.0_dp*(u(nxf,2,:)+u(nxf-4,2,:))+15.0_dp*(u(nxf-1,2,:)+u(nxf-3,2,:))+u(nxf-5,2,:)) &
      -fy*(u(nxf-2,5,:)-6.0_dp*(u(nxf-2,4,:)+rbufS(nxf-2,3,:))+15.0_dp*(u(nxf-2,3,:)+u(nxf-2,1,:))+rbufS(nxf-2,2,:)))
    reshd(nxf-1,2,:)=rhs(nxf-1,2,:)-(fc*u(nxf-1,2,:) &
      -fx*(Rq*u(nxf-1,2,:)+Rs(2,:)-6.0_dp*(Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-3,2,:))+15.0_dp*(u(nxf,2,:)+u(nxf-2,2,:))+u(nxf-4,2,:)) &
      -fy*(u(nxf-1,5,:)-6.0_dp*(u(nxf-1,4,:)+rbufS(nxf-1,3,:))+15.0_dp*(u(nxf-1,3,:)+u(nxf-1,1,:))+rbufS(nxf-1,2,:)))
    reshd(nxf,2,:)=rhs(nxf,2,:)-(fc*u(nxf,2,:) &
      -fx*(Rq*u(nxf-2,2,:)+Rs(3,:)-6.0_dp*(Rq*u(nxf-1,2,:)+Rs(2,:)+u(nxf-2,2,:))+15.0_dp*(Rq*u(nxf,2,:)+Rs(1,:)+u(nxf-1,2,:))+u(nxf-3,2,:)) &
      -fy*(u(nxf,5,:)-6.0_dp*(u(nxf,4,:)+rbufS(nxf,3,:))+15.0_dp*(u(nxf,3,:)+u(nxf,1,:))+rbufS(nxf,2,:)))
    reshd(nxf-2,3,:)=rhs(nxf-2,3,:)-(fc*u(nxf-2,3,:) &
      -fx*(Rq*u(nxf,3,:)+Rs(1,:)-6.0_dp*(u(nxf,3,:)+u(nxf-4,3,:))+15.0_dp*(u(nxf-1,3,:)+u(nxf-3,3,:))+u(nxf-5,3,:)) &
      -fy*(u(nxf-2,6,:)-6.0_dp*(u(nxf-2,5,:)+u(nxf-2,1,:))+15.0_dp*(u(nxf-2,4,:)+u(nxf-2,2,:))+rbufS(nxf-2,3,:)))
    reshd(nxf-1,3,:)=rhs(nxf-1,3,:)-(fc*u(nxf-1,3,:) &
      -fx*(Rq*u(nxf-1,3,:)+Rs(2,:)-6.0_dp*(Rq*u(nxf,3,:)+Rs(1,:)+u(nxf-3,3,:))+15.0_dp*(u(nxf,3,:)+u(nxf-2,3,:))+u(nxf-4,3,:)) &
      -fy*(u(nxf-1,6,:)-6.0_dp*(u(nxf-1,5,:)+u(nxf-1,1,:))+15.0_dp*(u(nxf-1,4,:)+u(nxf-1,2,:))+rbufS(nxf-1,3,:)))
    reshd(nxf,3,:)=rhs(nxf,3,:)-(fc*u(nxf,3,:) &
      -fx*(Rq*u(nxf-2,3,:)+Rs(3,:)-6.0_dp*(Rq*u(nxf-1,3,:)+Rs(2,:)+u(nxf-2,3,:))+15.0_dp*(Rq*u(nxf,3,:)+Rs(1,:)+u(nxf-1,3,:))+u(nxf-3,3,:)) &
      -fy*(u(nxf,6,:)-6.0_dp*(u(nxf,5,:)+u(nxf,1,:))+15.0_dp*(u(nxf,4,:)+u(nxf,2,:))+rbufS(nxf,3,:)))
    do j=4,nyf-3 ! Right boundary
      reshd(nxf-2,j,:)=rhs(nxf-2,j,:)-(fc*u(nxf-2,j,:) &
        -fx*(Rq*u(nxf,j,:)+Rs(1,:)-6.0_dp*(u(nxf,j,:)+u(nxf-4,j,:))+15.0_dp*(u(nxf-1,j,:)+u(nxf-3,j,:))+u(nxf-5,j,:)) &
        -fy*(u(nxf-2,j+3,:)-6.0_dp*(u(nxf-2,j+2,:)+u(nxf-2,j-2,:))+15.0_dp*(u(nxf-2,j+1,:)+u(nxf-2,j-1,:))+u(nxf-2,j-3,:)))
      reshd(nxf-1,j,:)=rhs(nxf-1,j,:)-(fc*u(nxf-1,j,:) &
        -fx*(Rq*u(nxf-1,j,:)+Rs(2,:)-6.0_dp*(Rq*u(nxf,j,:)+Rs(1,:)+u(nxf-3,j,:))+15.0_dp*(u(nxf,j,:)+u(nxf-2,j,:))+u(nxf-4,j,:)) &
        -fy*(u(nxf-1,j+3,:)-6.0_dp*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:))+15.0_dp*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:))+u(nxf-1,j-3,:)))
      reshd(nxf,j,:)=rhs(nxf,j,:)-(fc*u(nxf,j,:) &
        -fx*(Rq*u(nxf-2,j,:)+Rs(3,:)-6.0_dp*(Rq*u(nxf-1,j,:)+Rs(2,:)+u(nxf-2,j,:))+15.0_dp*(Rq*u(nxf,j,:)+Rs(1,:)+u(nxf-1,j,:))+u(nxf-3,j,:)) &
        -fy*(u(nxf,j+3,:)-6.0_dp*(u(nxf,j+2,:)+u(nxf,j-2,:))+15.0_dp*(u(nxf,j+1,:)+u(nxf,j-1,:))+u(nxf,j-3,:)))
    enddo
!...NE corner
    reshd(nxf-2,nyf-2,:)=rhs(nxf-2,nyf-2,:)-(fc*u(nxf-2,nyf-2,:) &
      -fx*(Rq*u(nxf,nyf-2,:)+Rs(1,:)-6.0_dp*(u(nxf,nyf-2,:)+u(nxf-4,nyf-2,:))+15.0_dp*(u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:))+u(nxf-5,nyf-2,:)) &
      -fy*(rbufN(nxf-2,1,:)-6.0_dp*(u(nxf-2,nyf,:)+u(nxf-2,nyf-4,:))+15.0_dp*(u(nxf-2,nyf-1,:)+u(nxf-2,nyf-3,:))+u(nxf-2,nyf-5,:)))
    reshd(nxf-1,nyf-2,:)=rhs(nxf-1,nyf-2,:)-(fc*u(nxf-1,nyf-2,:) &
      -fx*(Rq*u(nxf-1,nyf-2,:)+Rs(2,:)-6.0_dp*(Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-3,nyf-2,:))+15.0_dp*(u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:))+u(nxf-4,nyf-2,:)) &
      -fy*(rbufN(nxf-1,1,:)-6.0_dp*(u(nxf-1,nyf,:)+u(nxf-1,nyf-4,:))+15.0_dp*(u(nxf-1,nyf-1,:)+u(nxf-1,nyf-3,:))+u(nxf-1,nyf-5,:)))
    reshd(nxf,nyf-2,:)=rhs(nxf,nyf-2,:)-(fc*u(nxf,nyf-2,:) &
      -fx*(Rq*u(nxf-2,nyf-2,:)+Rs(3,:)-6.0_dp*(Rq*u(nxf-1,nyf-2,:)+Rs(2,:)+u(nxf-2,nyf-2,:))+15.0_dp*(Rq*u(nxf,nyf-2,:)+Rs(1,:)+u(nxf-1,nyf-2,:))+u(nxf-3,nyf-2,:)) &
      -fy*(rbufN(nxf,1,:)-6.0_dp*(u(nxf,nyf,:)+u(nxf,nyf-4,:))+15.0_dp*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:))+u(nxf,nyf-5,:)))
    reshd(nxf-2,nyf-1,:)=rhs(nxf-2,nyf-1,:)-(fc*u(nxf-2,nyf-1,:) &
      -fx*(Rq*u(nxf,nyf-1,:)+Rs(1,:)-6.0_dp*(u(nxf,nyf-1,:)+u(nxf-4,nyf-1,:))+15.0_dp*(u(nxf-1,nyf-1,:)+u(nxf-3,nyf-1,:))+u(nxf-5,nyf-1,:)) &
      -fy*(rbufN(nxf-2,2,:)-6.0_dp*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-3,:))+15.0_dp*(u(nxf-2,nyf,:)+u(nxf-2,nyf-2,:))+u(nxf-2,nyf-4,:)))
    reshd(nxf-1,nyf-1,:)=rhs(nxf-1,nyf-1,:)-(fc*u(nxf-1,nyf-1,:) &
      -fx*(Rq*u(nxf-1,nyf-1,:)+Rs(2,:)-6.0_dp*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-3,nyf-1,:))+15.0_dp*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:))+u(nxf-4,nyf-1,:)) &
      -fy*(rbufN(nxf-1,2,:)-6.0_dp*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:))+15.0_dp*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:))+u(nxf-1,nyf-4,:)))
    reshd(nxf,nyf-1,:)=rhs(nxf,nyf-1,:)-(fc*u(nxf,nyf-1,:) &
      -fx*(Rq*u(nxf-2,nyf-1,:)+Rs(3,:)-6.0_dp*(Rq*u(nxf-1,nyf-1,:)+Rs(2,:)+u(nxf-2,nyf-1,:))+15.0_dp*(Rq*u(nxf,nyf-1,:)+Rs(1,:)+u(nxf-1,nyf-1,:))+u(nxf-3,nyf-1,:)) &
      -fy*(rbufN(nxf,2,:)-6.0_dp*(rbufN(nxf,1,:)+u(nxf,nyf-3,:))+15.0_dp*(u(nxf,nyf,:)+u(nxf,nyf-2,:))+u(nxf,nyf-4,:)))
    reshd(nxf-2,nyf,:)=rhs(nxf-2,nyf,:)-(fc*u(nxf-2,nyf,:) &
      -fx*(Rq*u(nxf,nyf,:)+Rs(1,:)-6.0_dp*(u(nxf,nyf,:)+u(nxf-4,nyf,:))+15.0_dp*(u(nxf-1,nyf,:)+u(nxf-3,nyf,:))+u(nxf-5,nyf,:)) &
      -fy*(rbufN(nxf-2,3,:)-6.0_dp*(rbufN(nxf-2,2,:)+u(nxf-2,nyf-2,:))+15.0_dp*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-1,:))+u(nxf-2,nyf-3,:)))
    reshd(nxf-1,nyf,:)=rhs(nxf-1,nyf,:)-(fc*u(nxf-1,nyf,:) &
      -fx*(Rq*u(nxf-1,nyf,:)+Rs(2,:)-6.0_dp*(Rq*u(nxf,nyf,:)+Rs(1,:)+u(nxf-3,nyf,:))+15.0_dp*(u(nxf,nyf,:)+u(nxf-2,nyf,:))+u(nxf-4,nyf,:)) &
      -fy*(rbufN(nxf-1,3,:)-6.0_dp*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:))+15.0_dp*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:))+u(nxf-1,nyf-3,:)))
    reshd(nxf,nyf,:)=rhs(nxf,nyf,:)-(fc*u(nxf,nyf,:) &
      -fx*(Rq*u(nxf-2,nyf,:)+Rs(3,:)-6.0_dp*(Rq*u(nxf-1,nyf,:)+Rs(2,:)+u(nxf-2,nyf,:))+15.0_dp*(Rq*u(nxf,nyf,:)+Rs(1,:)+u(nxf-1,nyf,:))+u(nxf-3,nyf,:)) &
      -fy*(rbufN(nxf,3,:)-6.0_dp*(rbufN(nxf,2,:)+u(nxf,nyf-2,:))+15.0_dp*(rbufN(nxf,1,:)+u(nxf,nyf-1,:))+u(nxf,nyf-3,:)))
  else
    reshd(nxf-2,1,:)=rhs(nxf-2,1,:)-(fc*u(nxf-2,1,:) &
      -fx*(rbufE(1,1,:)-6.0_dp*(u(nxf,1,:)+u(nxf-4,1,:))+15.0_dp*(u(nxf-1,1,:)+u(nxf-3,1,:))+u(nxf-5,1,:)) &
      -fy*(u(nxf-2,4,:)-6.0_dp*(u(nxf-2,3,:)+rbufS(nxf-2,2,:))+15.0_dp*(u(nxf-2,2,:)+rbufS(nxf-2,3,:))+rbufS(nxf-2,1,:)))
    reshd(nxf-1,1,:)=rhs(nxf-1,1,:)-(fc*u(nxf-1,1,:) &
      -fx*(rbufE(2,1,:)-6.0_dp*(rbufE(1,1,:)+u(nxf-3,1,:))+15.0_dp*(u(nxf,1,:)+u(nxf-2,1,:))+u(nxf-4,1,:)) &
      -fy*(u(nxf-1,4,:)-6.0_dp*(u(nxf-1,3,:)+rbufS(nxf-1,2,:))+15.0_dp*(u(nxf-1,2,:)+rbufS(nxf-1,3,:))+rbufS(nxf-1,1,:)))
    reshd(nxf,1,:)=rhs(nxf,1,:)-(fc*u(nxf,1,:) &
      -fx*(rbufE(3,1,:)-6.0_dp*(rbufE(2,1,:)+u(nxf-2,1,:))+15.0_dp*(rbufE(1,1,:)+u(nxf-1,1,:))+u(nxf-3,1,:)) &
      -fy*(u(nxf,4,:)-6.0_dp*(u(nxf,3,:)+rbufS(nxf,2,:))+15.0_dp*(u(nxf,2,:)+rbufS(nxf,3,:))+rbufS(nxf,1,:)))
    reshd(nxf-2,2,:)=rhs(nxf-2,2,:)-(fc*u(nxf-2,2,:) &
      -fx*(rbufE(1,2,:)-6.0_dp*(u(nxf,2,:)+u(nxf-4,2,:))+15.0_dp*(u(nxf-1,2,:)+u(nxf-3,2,:))+u(nxf-5,2,:)) &
      -fy*(u(nxf-2,5,:)-6.0_dp*(u(nxf-2,4,:)+rbufS(nxf-2,3,:))+15.0_dp*(u(nxf-2,3,:)+u(nxf-2,1,:))+rbufS(nxf-2,2,:)))
    reshd(nxf-1,2,:)=rhs(nxf-1,2,:)-(fc*u(nxf-1,2,:) &
      -fx*(rbufE(2,2,:)-6.0_dp*(rbufE(1,2,:)+u(nxf-3,2,:))+15.0_dp*(u(nxf,2,:)+u(nxf-2,2,:))+u(nxf-4,2,:)) &
      -fy*(u(nxf-1,5,:)-6.0_dp*(u(nxf-1,4,:)+rbufS(nxf-1,3,:))+15.0_dp*(u(nxf-1,3,:)+u(nxf-1,1,:))+rbufS(nxf-1,2,:)))
    reshd(nxf,2,:)=rhs(nxf,2,:)-(fc*u(nxf,2,:) &
      -fx*(rbufE(3,2,:)-6.0_dp*(rbufE(2,2,:)+u(nxf-2,2,:))+15.0_dp*(rbufE(1,2,:)+u(nxf-1,2,:))+u(nxf-3,2,:)) &
      -fy*(u(nxf,5,:)-6.0_dp*(u(nxf,4,:)+rbufS(nxf,3,:))+15.0_dp*(u(nxf,3,:)+u(nxf,1,:))+rbufS(nxf,2,:)))
    reshd(nxf-2,3,:)=rhs(nxf-2,3,:)-(fc*u(nxf-2,3,:) &
      -fx*(rbufE(1,3,:)-6.0_dp*(u(nxf,3,:)+u(nxf-4,3,:))+15.0_dp*(u(nxf-1,3,:)+u(nxf-3,3,:))+u(nxf-5,3,:)) &
      -fy*(u(nxf-2,6,:)-6.0_dp*(u(nxf-2,5,:)+u(nxf-2,1,:))+15.0_dp*(u(nxf-2,4,:)+u(nxf-2,2,:))+rbufS(nxf-2,3,:)))
    reshd(nxf-1,3,:)=rhs(nxf-1,3,:)-(fc*u(nxf-1,3,:) &
      -fx*(rbufE(2,3,:)-6.0_dp*(rbufE(1,3,:)+u(nxf-3,3,:))+15.0_dp*(u(nxf,3,:)+u(nxf-2,3,:))+u(nxf-4,3,:)) &
      -fy*(u(nxf-1,6,:)-6.0_dp*(u(nxf-1,5,:)+u(nxf-1,1,:))+15.0_dp*(u(nxf-1,4,:)+u(nxf-1,2,:))+rbufS(nxf-1,3,:)))
    reshd(nxf,3,:)=rhs(nxf,3,:)-(fc*u(nxf,3,:) &
      -fx*(rbufE(3,3,:)-6.0_dp*(rbufE(2,3,:)+u(nxf-2,3,:))+15.0_dp*(rbufE(1,3,:)+u(nxf-1,3,:))+u(nxf-3,3,:)) &
      -fy*(u(nxf,6,:)-6.0_dp*(u(nxf,5,:)+u(nxf,1,:))+15.0_dp*(u(nxf,4,:)+u(nxf,2,:))+rbufS(nxf,3,:)))
    do j=4,nyf-3
      reshd(nxf-2,j,:)=rhs(nxf-2,j,:)-(fc*u(nxf-2,j,:) &
        -fx*(rbufE(1,j,:)-6.0_dp*(u(nxf,j,:)+u(nxf-4,j,:))+15.0_dp*(u(nxf-1,j,:)+u(nxf-3,j,:))+u(nxf-5,j,:)) &
        -fy*(u(nxf-2,j+3,:)-6.0_dp*(u(nxf-2,j+2,:)+u(nxf-2,j-2,:))+15.0_dp*(u(nxf-2,j+1,:)+u(nxf-2,j-1,:))+u(nxf-2,j-3,:)))
      reshd(nxf-1,j,:)=rhs(nxf-1,j,:)-(fc*u(nxf-1,j,:) &
        -fx*(rbufE(2,j,:)-6.0_dp*(rbufE(1,j,:)+u(nxf-3,j,:))+15.0_dp*(u(nxf,j,:)+u(nxf-2,j,:))+u(nxf-4,j,:)) &
        -fy*(u(nxf-1,j+3,:)-6.0_dp*(u(nxf-1,j+2,:)+u(nxf-1,j-2,:))+15.0_dp*(u(nxf-1,j+1,:)+u(nxf-1,j-1,:))+u(nxf-1,j-3,:)))
      reshd(nxf,j,:)=rhs(nxf,j,:)-(fc*u(nxf,j,:) &
        -fx*(rbufE(3,j,:)-6.0_dp*(rbufE(2,j,:)+u(nxf-2,j,:))+15.0_dp*(rbufE(1,j,:)+u(nxf-1,j,:))+u(nxf-3,j,:)) &
        -fy*(u(nxf,j+3,:)-6.0_dp*(u(nxf,j+2,:)+u(nxf,j-2,:))+15.0_dp*(u(nxf,j+1,:)+u(nxf,j-1,:))+u(nxf,j-3,:)))
    enddo
    reshd(nxf-2,nyf-2,:)=rhs(nxf-2,nyf-2,:)-(fc*u(nxf-2,nyf-2,:) &
      -fx*(rbufE(1,nyf-2,:)-6.0_dp*(u(nxf,nyf-2,:)+u(nxf-4,nyf-2,:))+15.0_dp*(u(nxf-1,nyf-2,:)+u(nxf-3,nyf-2,:))+u(nxf-5,nyf-2,:)) &
      -fy*(rbufN(nxf-2,1,:)-6.0_dp*(u(nxf-2,nyf,:)+u(nxf-2,nyf-4,:))+15.0_dp*(u(nxf-2,nyf-1,:)+u(nxf-2,nyf-3,:))+u(nxf-2,nyf-5,:)))
    reshd(nxf-1,nyf-2,:)=rhs(nxf-1,nyf-2,:)-(fc*u(nxf-1,nyf-2,:) &
      -fx*(rbufE(2,nyf-2,:)-6.0_dp*(rbufE(1,nyf-2,:)+u(nxf-3,nyf-2,:))+15.0_dp*(u(nxf,nyf-2,:)+u(nxf-2,nyf-2,:))+u(nxf-4,nyf-2,:)) &
      -fy*(rbufN(nxf-1,1,:)-6.0_dp*(u(nxf-1,nyf,:)+u(nxf-1,nyf-4,:))+15.0_dp*(u(nxf-1,nyf-1,:)+u(nxf-1,nyf-3,:))+u(nxf-1,nyf-5,:)))
    reshd(nxf,nyf-2,:)=rhs(nxf,nyf-2,:)-(fc*u(nxf,nyf-2,:) &
      -fx*(rbufE(3,nyf-2,:)-6.0_dp*(rbufE(2,nyf-2,:)+u(nxf-2,nyf-2,:))+15.0_dp*(rbufE(1,nyf-2,:)+u(nxf-1,nyf-2,:))+u(nxf-3,nyf-2,:)) &
      -fy*(rbufN(nxf,1,:)-6.0_dp*(u(nxf,nyf,:)+u(nxf,nyf-4,:))+15.0_dp*(u(nxf,nyf-1,:)+u(nxf,nyf-3,:))+u(nxf,nyf-5,:)))
    reshd(nxf-2,nyf-1,:)=rhs(nxf-2,nyf-1,:)-(fc*u(nxf-2,nyf-1,:) &
      -fx*(rbufE(1,nyf-1,:)-6.0_dp*(u(nxf,nyf-1,:)+u(nxf-4,nyf-1,:))+15.0_dp*(u(nxf-1,nyf-1,:)+u(nxf-3,nyf-1,:))+u(nxf-5,nyf-1,:)) &
      -fy*(rbufN(nxf-2,2,:)-6.0_dp*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-3,:))+15.0_dp*(u(nxf-2,nyf,:)+u(nxf-2,nyf-2,:))+u(nxf-2,nyf-4,:)))
    reshd(nxf-1,nyf-1,:)=rhs(nxf-1,nyf-1,:)-(fc*u(nxf-1,nyf-1,:) &
      -fx*(rbufE(2,nyf-1,:)-6.0_dp*(rbufE(1,nyf-1,:)+u(nxf-3,nyf-1,:))+15.0_dp*(u(nxf,nyf-1,:)+u(nxf-2,nyf-1,:))+u(nxf-4,nyf-1,:)) &
      -fy*(rbufN(nxf-1,2,:)-6.0_dp*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-3,:))+15.0_dp*(u(nxf-1,nyf,:)+u(nxf-1,nyf-2,:))+u(nxf-1,nyf-4,:)))
    reshd(nxf,nyf-1,:)=rhs(nxf,nyf-1,:)-(fc*u(nxf,nyf-1,:) &
      -fx*(rbufE(3,nyf-1,:)-6.0_dp*(rbufE(2,nyf-1,:)+u(nxf-2,nyf-1,:))+15.0_dp*(rbufE(1,nyf-1,:)+u(nxf-1,nyf-1,:))+u(nxf-3,nyf-1,:)) &
      -fy*(rbufN(nxf,2,:)-6.0_dp*(rbufN(nxf,1,:)+u(nxf,nyf-3,:))+15.0_dp*(u(nxf,nyf,:)+u(nxf,nyf-2,:))+u(nxf,nyf-4,:)))
    reshd(nxf-2,nyf,:)=rhs(nxf-2,nyf,:)-(fc*u(nxf-2,nyf,:) &
      -fx*(rbufE(1,nyf,:)-6.0_dp*(u(nxf,nyf,:)+u(nxf-4,nyf,:))+15.0_dp*(u(nxf-1,nyf,:)+u(nxf-3,nyf,:))+u(nxf-5,nyf,:)) &
      -fy*(rbufN(nxf-2,3,:)-6.0_dp*(rbufN(nxf-2,2,:)+u(nxf-2,nyf-2,:))+15.0_dp*(rbufN(nxf-2,1,:)+u(nxf-2,nyf-1,:))+u(nxf-2,nyf-3,:)))
    reshd(nxf-1,nyf,:)=rhs(nxf-1,nyf,:)-(fc*u(nxf-1,nyf,:) &
      -fx*(rbufE(2,nyf,:)-6.0_dp*(rbufE(1,nyf,:)+u(nxf-3,nyf,:))+15.0_dp*(u(nxf,nyf,:)+u(nxf-2,nyf,:))+u(nxf-4,nyf,:)) &
      -fy*(rbufN(nxf-1,3,:)-6.0_dp*(rbufN(nxf-1,2,:)+u(nxf-1,nyf-2,:))+15.0_dp*(rbufN(nxf-1,1,:)+u(nxf-1,nyf-1,:))+u(nxf-1,nyf-3,:)))
    reshd(nxf,nyf,:)=rhs(nxf,nyf,:)-(fc*u(nxf,nyf,:) &
      -fx*(rbufE(3,nyf,:)-6.0_dp*(rbufE(2,nyf,:)+u(nxf-2,nyf,:))+15.0_dp*(rbufE(1,nyf,:)+u(nxf-1,nyf,:))+u(nxf-3,nyf,:)) &
      -fy*(rbufN(nxf,3,:)-6.0_dp*(rbufN(nxf,2,:)+u(nxf,nyf-2,:))+15.0_dp*(rbufN(nxf,1,:)+u(nxf,nyf-1,:))+u(nxf,nyf-3,:)))
  endif

#ENDIF

  call MPI_WAITALL(4,(/req(1),req(2),req(3),req(4)/), &
    (/stat(:,1),stat(:,2),stat(:,3),stat(:,4)/),ierr)

#IF (MGTIMER==1)
    call CPU_TIME(rest2)
    resthd=resthd+rest2-rest1
#ENDIF
  end function reshd
!*********************************************************
#ENDIF !.(HDop==2,3) else (HDop==8,12)
#ENDIF !.(HDop==4,6) else (HDop==2,3,8,12)
!*********************************************************
  recursive subroutine Gcycpoi(ng,u,rhs,xBC)
! Recursive multigrid iteration. On input, ng is the current level, u is
! the current value of the solution, and rhs is the right-hand side.
! On output u contains the improved solution at the current level.
  implicit none
  integer(I4B), intent(in) :: ng
  real(DP), intent(inout) :: u(:,:,:)
  real(DP), allocatable, intent(in) :: rhs(:,:,:)
  integer(I4B), intent(in) :: xBC(:)
  integer(I4B) :: jpost, jpre, i,j
  real(DP), allocatable :: res(:,:,:),v(:,:,:)

  iprocsm1=iprocsph(ng)-1 !.ID of last rank in x COMM
  jprocsm1=jprocsph(ng)-1 !.ID of last rank in y COMM

  if (ng==1) then        !.Bottom of gamma cycle: Solve on coarsest grid
#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do i=1,NU2ph
      call relaxpoi(ng,u,rhs,xBC)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    postXtpo=postXtpo+count2-count1
#ENDIF
  else                  !.On downward stroke towards coarser grids
#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpre=1,NU1ph      !.Pre-smoothing
      call relaxpoi(ng,u,rhs,xBC)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    preXtpo=preXtpo+count2-count1
#ENDIF

!...Number of grid points in coarse grid
    nxc=nxCph(ng)
    nyc=nyCph(ng)
!...Number of grid points in coarse grid after subdomain unification
    nxu=nxFph(ng-1)
    nyu=nyFph(ng-1)
    allocate(res(nxu,nyu,nzL),v(nxu,nyu,nzL))
#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
!...Restricted residual is the next RHS
    res=rstrct(respoi(ng,u,rhs,xBC),wdcph(ng),xBC, &
      acuph(ng-1,:),urnkph(ng-1)%a,nurnkph(ng-1,:))
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    rNrtpo=rNrtpo+count2-count1
#ENDIF

!...Next fine grid is the coarse, or unified subdomain, grid
    nxf=nxu
    nyf=nyu

    v=0.0_dp               !.Zero for initial guess in next relaxation

    if (acuph(ng-1,1)) then
!.....Only active processes proceed to next coarser grid
      do i=1,GAMph
!.......Recursive call for the coarse grid correction
        call Gcycpoi(ng-1,v,res,xBC)
      enddo
    endif

!...Redefine quantities that were modified in coarser grids
    iprocsm1=iprocsph(ng)-1
    jprocsm1=jprocsph(ng)-1
    nxc=nxCph(ng)
    nyc=nyCph(ng)
    nxf=nxFph(ng)
    nyf=nyFph(ng)

#IF (MGTIMER==1)
    call CPU_TIME(prolt1)
#ENDIF
!...Upward stroke -> finer grid: prolong coarse grid error + correct solution
    u=u+prolng(v,wdcph(ng),xBC,COMMph(ng,:),iIDph(ng),jprocsph(ng), &
      nborph(ng,:),acuph(ng-1,:),urnkph(ng-1)%a,nurnkph(ng-1,:))
#IF (MGTIMER==1)
    call CPU_TIME(prolt2)
    proltpo=proltpo+prolt2-prolt1
#ENDIF

#IF (MGTIMER==1)
    call CPU_TIME(count1)
#ENDIF
    do jpost=1,NU2ph      !.Post-smoothing
      call relaxpoi(ng,u,rhs,xBC)
    enddo
#IF (MGTIMER==1)
    call CPU_TIME(count2)
    postXtpo=postXtpo+count2-count1
#ENDIF
  endif

  end subroutine Gcycpoi
!*********************************************************
  subroutine relaxpoi(ng,u,rhs,xBC)
!.Red-black Gauss-Seidel relaxation of Poisson equation
!.and calculation of the residue res=rhs-div(lam*grad(u)).
  implicit none
  include 'mpif.h'
  real(DP), intent(inout) :: u(:,:,:)
  real(DP), allocatable, intent(in) :: rhs(:,:,:)
  integer(I4B), intent(in) :: ng,xBC(:)
  real(DP), allocatable :: Ls(:),rs(:),loc1D(:)
  integer(I4B) :: i,j,k,il,nxfd2,nyfd2,Nrank,Srank,Erank,Wrank
  real(DP) :: fx,fy,fc,delta,deltaL,deltaR
  real(DP) :: Lq,Rq,Lr,rr
  real(DP), allocatable, dimension(:,:) :: sbufW,rbufE,sbufE,rbufW, &
                                           sbufN,rbufS,sbufS,rbufN
  integer(I4B) :: req(4),stat(MPI_STATUS_SIZE,4),ierr

  nxfd2=nxf/2
  nyfd2=nyf/2

  Nrank=nborph(ng,1)
  Erank=nborph(ng,3)
  Srank=nborph(ng,5)
  Wrank=nborph(ng,7)

!.Boundary conditions are implemented with some amount of flexibility.
!.The right boundary for example, if u(nxf+1) is needed, it will be
!.written as
!.  u(nxf+1) = Rq*u(nxf) + Rr*u(nxf-1) + Rs
!.so that symmetry, extrapolation or Dirichlet conditions can be
!.easily imposed.
!......... BCs of u ................................!
!.LHS BC of u
  allocate(Ls(nzL),loc1D(nzL))
  if (xBC(1)==0) then
!...linear extrapolation w/o boundary value
    Lq=2.0_dp
    Lr=-1.0_dp
    ls=0.0_dp
  elseif (xBC(1)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Lq=-1.0_dp
    Lr=0.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u(1,:,k))
    call MPI_ALLreduce(loc1D,ls,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMph(ng,1),ierr)
    ls=2.0_dp*ls/dble(nyf*jprocsph(ng))
  else ! if (xBC(1)==1 .OR. xBC(1)==-1) then ! symmetry
    Lq=dble(xBC(1))
    Lr=0.0_dp
    ls=0.0_dp
  endif
!.RHS BC of u
  allocate(rs(nzL))
  if (xBC(2)==0) then
!...linear extrapolation w/o boundary value
    Rq=2.0_dp
    rr=-1.0_dp
    rs=0.0_dp
  elseif (xBC(2)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Rq=-1.0_dp
    rr=0.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u(nxf,:,k))
    call MPI_ALLreduce(loc1D,rs,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMph(ng,1),ierr)
    rs=2.0_dp*rs/dble(nyf*jprocsph(ng))
  else ! if (xBC(2)==1 .OR. xBC(2)==-1) then ! symmetry
    Rq=dble(xBC(2))
    rr=0.0_dp
    rs=0.0_dp
  endif
!........................................................!
  fx=(dble(nxf*iprocsph(ng))/bLx)**2
  fy=(dble(nyf*jprocsph(ng))/bLy)**2
  fc=2.0_dp*(fx+fy)
  delta=omeph/fc
  deltaL=omeph/(fc-fx*Lq)
  deltaR=omeph/(fc-fx*Rq)

!.Send buffer names indicate direction data is sent to
!.and receive buffers direction data is received from.
  allocate(sbufS(nxfd2,nzL),rbufN(nxfd2,nzL),sbufW(nyfd2,nzL),rbufE(nyfd2,nzL))
!.Comunicate black grid points needed
!.Send bottom row of u to South rank
  sbufS=u(2:nxf:2,1,:)
  call MPI_ISEND(sbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(1),ierr)
!.Receive u row from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(2),ierr)
!.Send first column of u to West rank
  sbufW=u(1,2:nyf:2,:)
  call MPI_ISEND(sbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(3),ierr)
!.Receive u column from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(4),ierr)

!.First do the even-even and odd-odd squares of the grid,
!.i.e., the red squares of the checker-board
  do j=2,nyf-2,2
    do i=2,nxf-2,2
      u(i,j,:)=Oomeph*u(i,j,:)+(fx*(u(i+1,j,:)+u(i-1,j,:)) &
        +fy*(u(i,j+1,:)+u(i,j-1,:))-rhs(i,j,:))*delta
    enddo
  enddo
!.Top, even i
  call MPI_WAIT(req(2),stat(:,2),ierr)
  do i=2,nxf-2,2
    u(i,nyf,:)=Oomeph*u(i,nyf,:)+(fx*(u(i+1,nyf,:)+u(i-1,nyf,:)) &
      +fy*(rbufN(i/2,:)+u(i,nyf-1,:))-rhs(i,nyf,:))*delta
  enddo
!.Right boundary, even j
    call MPI_WAIT(req(4),stat(:,4),ierr)
  if (iIDph(ng) .eq. iprocsm1) then !.Right most process (along x)
    do j=2,nyf-2,2
      u(nxf,j,:)=Oomeph*u(nxf,j,:)+(fx*(rs+u(nxf-1,j,:)) &
        +fy*(u(nxf,j+1,:)+u(nxf,j-1,:))-rhs(nxf,j,:))*deltaR
    enddo
    u(nxf,nyf,:)=Oomeph*u(nxf,nyf,:)+(fx*(rs+u(nxf-1,nyf,:)) &
      +fy*(rbufN(nxfd2,:)+u(nxf,nyf-1,:))-rhs(nxf,nyf,:))*deltaR
  else
    do j=2,nyf-2,2
      u(nxf,j,:)=Oomeph*u(nxf,j,:)+(fx*(rbufE(j/2,:)+u(nxf-1,j,:)) &
        +fy*(u(nxf,j+1,:)+u(nxf,j-1,:))-rhs(nxf,j,:))*delta
    enddo
    u(nxf,nyf,:)=Oomeph*u(nxf,nyf,:)+(fx*(rbufE(nyfd2,:)+u(nxf-1,nyf,:)) &
      +fy*(rbufN(nxfd2,:)+u(nxf,nyf-1,:))-rhs(nxf,nyf,:))*delta
  endif

    call MPI_WAITALL(2,(/req(1),req(3)/),(/stat(:,1),stat(:,3)/),ierr)
  allocate(sbufN(nxfd2,nzL),rbufS(nxfd2,nzL),sbufE(nyfd2,nzL),rbufW(nyfd2,nzL))
!.Send top row of array to North rank
  sbufN=u(1:nxf-1:2,nyf,:)
  call MPI_ISEND(sbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(1),ierr)
!.Receive u row from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(2),ierr)
!.Send last column of array to East rank
  sbufE=u(nxf,1:nyf-1:2,:)
  call MPI_ISEND(sbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(3),ierr)
!.Receive u column from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(4),ierr)
  do j=3,nyf-1,2
    do i=3,nxf-1,2
      u(i,j,:)=Oomeph*u(i,j,:)+(fx*(u(i+1,j,:)+u(i-1,j,:)) &
        +fy*(u(i,j+1,:)+u(i,j-1,:))-rhs(i,j,:))*delta
    enddo
  enddo
!.Bottom, odd i
  call MPI_WAIT(req(2),stat(:,2),ierr)
  do i=3,nxf-1,2
    u(i,1,:)=Oomeph*u(i,1,:)+(fx*(u(i+1,1,:)+u(i-1,1,:)) &
      +fy*(u(i,2,:)+rbufS((i+1)/2,:))-rhs(i,1,:))*delta
  enddo
!.Left boundary, odd j
  call MPI_WAIT(req(4),stat(:,4),ierr)
  if (iIDph(ng)==0) then !.Left most process (along x)
    u(1,1,:)=Oomeph*u(1,1,:)+(fx*(u(2,1,:)+Ls) &
      +fy*(u(1,2,:)+rbufS(1,:))-rhs(1,1,:))*deltaL ! SW corner
    do j=3,nyf-1,2
      u(1,j,:)=Oomeph*u(1,j,:)+(fx*(u(2,j,:)+Ls) &
        +fy*(u(1,j+1,:)+u(1,j-1,:))-rhs(1,j,:))*deltaL
    enddo
  else
    u(1,1,:)=Oomeph*u(1,1,:)+(fx*(u(2,1,:)+rbufW(1,:)) &
      +fy*(u(1,2,:)+rbufS(1,:))-rhs(1,1,:))*delta
    do j=3,nyf-1,2
      u(1,j,:)=Oomeph*u(1,j,:)+(fx*(u(2,j,:)+rbufW((j+1)/2,:)) &
        +fy*(u(1,j+1,:)+u(1,j-1,:))-rhs(1,j,:))*delta
    enddo
  endif

    call MPI_WAITALL(2,(/req(1),req(3)/),(/stat(:,1),stat(:,3)/),ierr)
!.Now do even-odd and odd-even squares of the grid, i.e, the black
!.squares of the checker-board
!.Comunicate red grid points needed
!.Send bottom row of u to South rank
  sbufS=u(1:nxf-1:2,1,:)
  call MPI_ISEND(sbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(1),ierr)
!.Receive u row from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(2),ierr)
!.Send last column of array to East rank
  sbufE=u(nxf,2:nyf:2,:)
  call MPI_ISEND(sbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(3),ierr)
!.Receive u column from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(4),ierr)
  do j=2,nyf-2,2
    do i=3,nxf-1,2
      u(i,j,:)=Oomeph*u(i,j,:)+(fx*(u(i+1,j,:)+u(i-1,j,:)) &
        +fy*(u(i,j+1,:)+u(i,j-1,:))-rhs(i,j,:))*delta
    enddo
  enddo
!.Top, odd i
  call MPI_WAIT(req(2),stat(:,2),ierr)
  do i=3,nxf-1,2
    u(i,nyf,:)=Oomeph*u(i,nyf,:)+(fx*(u(i+1,nyf,:)+u(i-1,nyf,:)) &
      +fy*(rbufN((i+1)/2,:)+u(i,nyf-1,:))-rhs(i,nyf,:))*delta
  enddo
!.Left boundary, even j
  call MPI_WAIT(req(4),stat(:,4),ierr)
  if (iIDph(ng)==0) then !.Left most process (along x)
    do j=2,nyf-2,2
      u(1,j,:)=Oomeph*u(1,j,:)+(fx*(u(2,j,:)+Ls) &
        +fy*(u(1,j+1,:)+u(1,j-1,:))-rhs(1,j,:))*deltaL
    enddo
    u(1,nyf,:)=Oomeph*u(1,nyf,:)+(fx*(u(2,nyf,:)+Ls) &
      +fy*(rbufN(1,:)+u(1,nyf-1,:))-rhs(1,nyf,:))*deltaL !.NW corner
  else
    do j=2,nyf-2,2
      u(1,j,:)=Oomeph*u(1,j,:)+(fx*(u(2,j,:)+rbufW(j/2,:)) &
        +fy*(u(1,j+1,:)+u(1,j-1,:))-rhs(1,j,:))*delta
    enddo
    u(1,nyf,:)=Oomeph*u(1,nyf,:)+(fx*(u(2,nyf,:)+rbufW(nyfd2,:)) &
      +fy*(rbufN(1,:)+u(1,nyf-1,:))-rhs(1,nyf,:))*delta
  endif

    call MPI_WAITALL(2,(/req(1),req(3)/),(/stat(:,1),stat(:,3)/),ierr)
!.Send top row of array to North rank
  sbufN=u(2:nxf:2,nyf,:)
  call MPI_ISEND(sbufN,nxfd2*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(1),ierr)
!.Receive u row from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxfd2*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(2),ierr)
!.Send first column of u to West rank
  sbufW=u(1,1:nyf-1:2,:)
  call MPI_ISEND(sbufW,nyfd2*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(3),ierr)
!.Receive u column from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,nyfd2*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(4),ierr)
  do j=3,nyf-1,2
    do i=2,nxf-2,2
      u(i,j,:)=Oomeph*u(i,j,:)+(fx*(u(i+1,j,:)+u(i-1,j,:)) &
        +fy*(u(i,j+1,:)+u(i,j-1,:))-rhs(i,j,:))*delta
    enddo
  enddo
!.Bottom, even i
  call MPI_WAIT(req(2),stat(:,2),ierr)
  do i=2,nxf-2,2
    u(i,1,:)=Oomeph*u(i,1,:)+(fx*(u(i+1,1,:)+u(i-1,1,:)) &
      +fy*(u(i,2,:)+rbufS(i/2,:))-rhs(i,1,:))*delta
  enddo
!.Right boundary, odd j
  call MPI_WAIT(req(4),stat(:,4),ierr)
  if (iIDph(ng)==iprocsm1) then !.Right most process (along x)
    u(nxf,1,:)=Oomeph*u(nxf,1,:)+(fx*(rs+u(nxf-1,1,:)) &
      +fy*(u(nxf,2,:)+rbufS(nxfd2,:))-rhs(nxf,1,:))*deltaR !.SE corner
    do j=3,nyf-1,2
      u(nxf,j,:)=Oomeph*u(nxf,j,:)+(fx*(rs+u(nxf-1,j,:)) &
        +fy*(u(nxf,j+1,:)+u(nxf,j-1,:))-rhs(nxf,j,:))*deltaR
    enddo
  else
    u(nxf,1,:)=Oomeph*u(nxf,1,:)+(fx*(rbufE(1,:)+u(nxf-1,1,:)) &
      +fy*(u(nxf,2,:)+rbufS(nxfd2,:))-rhs(nxf,1,:))*delta
    do j=3,nyf-1,2
      u(nxf,j,:)=Oomeph*u(nxf,j,:)+(fx*(rbufE((j+1)/2,:)+u(nxf-1,j,:)) &
        +fy*(u(nxf,j+1,:)+u(nxf,j-1,:))-rhs(nxf,j,:))*delta
    enddo
  endif

    call MPI_WAITALL(2,(/req(1),req(3)/),(/stat(:,1),stat(:,3)/),ierr)
  end subroutine relaxpoi
!*********************************************************
  function respoi(ng,u,rhs,xBC)
! Returns minus the residual. Input quantities are u and rhs, while
! the residual is returned in res.
  implicit none
  include 'mpif.h'
  real(DP), dimension(:,:,:), intent(in) :: u, rhs
  integer(I4B), intent(in) :: ng,xBC(:)
  real(DP), allocatable :: respoi(:,:,:)
  real(DP), allocatable :: Ls(:),rs(:),loc1D(:)
  integer(I4B) :: i,j,k,il,Nrank,Srank,Erank,Wrank
  real(DP) :: fx,fy,fc,Lq,Lr,Rq,Rr
  real(DP), allocatable, dimension(:,:) :: sbufW,rbufE,sbufE,rbufW, &
                                           sbufN,rbufS,sbufS,rbufN
  integer(I4B) :: req(8),stat(MPI_STATUS_SIZE,8)
  integer(I4B) :: ierr

#IF (MGTIMER==1)
    call CPU_TIME(rest1)
#ENDIF
  Nrank=nborph(ng,1)
  Erank=nborph(ng,3)
  Srank=nborph(ng,5)
  Wrank=nborph(ng,7)

  allocate(respoi(nxf,nyf,nzL))
!.Boundary conditions are implemented with some amount of flexibility.
!.The right boundary for example, if u(nxf+1) is needed, it will be
!.written as
!.  u(nxf+1) = Rq*u(nxf) + Rr*u(nxf-1) + Rs
!.so that symmetry, extrapolation or Dirichlet conditions can be
!.easily imposed.
!.LHS BC of u
!  allocate(Lsu(nyf,nzL))
  allocate(Ls(nzL),loc1D(nzL))
  if (xBC(1)==0) then
!...linear extrapolation w/o boundary value
    Lq=2.0_dp
    Lr=-1.0_dp
    ls=0.0_dp
  elseif (xBC(1)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Lq=-1.0_dp
    Lr=0.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u(1,:,k))
    call MPI_ALLreduce(loc1D,ls,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMph(ng,1),ierr)
    ls=2.0_dp*ls/dble(nyf*jprocsph(ng))
  else ! if (xBC(1)==1 .OR. xBC(1)==-1) then ! symmetry
    Lq=dble(xBC(1))
    Lr=0.0_dp
    ls=0.0_dp
  endif
!.RHS BC of u
  allocate(Rs(nzL))
  if (xBC(2)==0) then
!...linear extrapolation w/o boundary value
    Rq=2.0_dp
    Rr=-1.0_dp
    rs=0.0_dp
  elseif (xBC(2)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Rq=-1.0_dp
    Rr=0.0_dp
    forall(k=1:nzL) loc1D(k)=sum(u(nxf,:,k))
    call MPI_ALLreduce(loc1D,rs,nzL,MPI_DOUBLE_PRECISION, &
                       MPI_SUM,COMMph(ng,1),ierr)
    rs=2.0_dp*rs/dble(nyf*jprocsph(ng))
  else ! if (xBC(2)==1 .OR. xBC(2)==-1) then ! symmetry
    Rq=dble(xBC(2))
    Rr=0.0_dp
    rs=0.0_dp
  endif
  fx=(dble(nxf*iprocsph(ng))/bLx)**2         !.1/(hx^2)
  fy=(dble(nyf*jprocsph(ng))/bLy)**2         !.1/(hy^2)
  fc=2.0_dp*(fx+fy)

  allocate(sbufS(nxf,nzL),rbufN(nxf,nzL),sbufW(nyf,nzL),rbufE(nyf,nzL))
  allocate(sbufN(nxf,nzL),rbufS(nxf,nzL),sbufE(nyf,nzL),rbufW(nyf,nzL))
!.Send top row of u to North rank
  sbufN=u(:,nyf,:)
  call MPI_ISEND(sbufN,nxf*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(1),ierr)
!.Send right column of u to East rank
  sbufE=u(nxf,:,:)
  call MPI_ISEND(sbufE,nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(2),ierr)
!.Send bottom row of u to South rank
  sbufS=u(:,1,:)
  call MPI_ISEND(sbufS,nxf*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(3),ierr)
!.Send left column of u to West rank
  sbufW=u(1,:,:)
  call MPI_ISEND(sbufW,nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(4),ierr)
!.Receive u row from South rank, for u(i,j-1)
  call MPI_IRECV(rbufS,nxf*nzL,MPI_DOUBLE_PRECISION,Srank, &
                 0,COMMph(ng,1),req(5),ierr)
!.Receive u column from West rank, for u(i-1,j)
  call MPI_IRECV(rbufW,nyf*nzL,MPI_DOUBLE_PRECISION,Wrank, &
                 0,COMMph(ng,2),req(6),ierr)
!.Receive u row from North rank, for u(i,j+1)
  call MPI_IRECV(rbufN,nxf*nzL,MPI_DOUBLE_PRECISION,Nrank, &
                 0,COMMph(ng,1),req(7),ierr)
!.Receive u column from East rank, for u(i+1,j)
  call MPI_IRECV(rbufE,nyf*nzL,MPI_DOUBLE_PRECISION,Erank, &
                 0,COMMph(ng,2),req(8),ierr)
!.Interior points
  do j=2,nyf-1
    do i=2,nxf-1
      respoi(i,j,:)=rhs(i,j,:)-(fx*(u(i+1,j,:)+u(i-1,j,:)) &
        +fy*(u(i,j+1,:)+u(i,j-1,:))-fc*u(i,j,:))
    enddo
  enddo

!.Boundary points
  call MPI_WAIT(req(5),stat(:,5),ierr)
  do i=2,nxf-1 !.Bottom boundary
    respoi(i,1,:)=rhs(i,1,:)-(fx*(u(i+1,1,:)+u(i-1,1,:)) &
      +fy*(u(i,2,:)+rbufS(i,:))-fc*u(i,1,:))
  enddo
  call MPI_WAIT(req(6),stat(:,6),ierr)
  if (iIDph(ng) == 0) then !.Left most process (along x)
    respoi(1,1,:)=rhs(1,1,:)-(fx*(u(2,1,:)+Lq*u(1,1,:)+Lr*u(2,1,:)+Ls) &
      +fy*(u(1,2,:)+rbufS(1,:))-fc*u(1,1,:)) !.SW corner
    do j=2,nyf-1 !.Left boundary
      respoi(1,j,:)=rhs(1,j,:)-(fx*(u(2,j,:)+Lq*u(1,j,:)+Lr*u(2,j,:)+Ls) &
        +fy*(u(1,j+1,:)+u(1,j-1,:))-fc*u(1,j,:))
    enddo
  else
    respoi(1,1,:)=rhs(1,1,:)-(fx*(u(2,1,:)+rbufW(1,:)) &
      +fy*(u(1,2,:)+rbufS(1,:))-fc*u(1,1,:))
    do j=2,nyf-1
      respoi(1,j,:)=rhs(1,j,:)-(fx*(u(2,j,:)+rbufW(j,:)) &
        +fy*(u(1,j+1,:)+u(1,j-1,:))-fc*u(1,j,:))
    enddo
  endif

  call MPI_WAIT(req(7),stat(:,7),ierr)
  if (iIDph(ng)==0) then !.Left most process (along x)
!...NW corner
    respoi(1,nyf,:)=rhs(1,nyf,:)-(fx*(u(2,nyf,:)+Lq*u(1,nyf,:)+Lr*u(2,nyf,:)+Ls) &
      +fy*(rbufN(1,:)+u(1,nyf-1,:))-fc*u(1,nyf,:))
  else
    respoi(1,nyf,:)=rhs(1,nyf,:)-(fx*(u(2,nyf,:)+rbufW(nyf,:)) &
      +fy*(rbufN(1,:)+u(1,nyf-1,:))-fc*u(1,nyf,:))
  endif
  do i=2,nxf-1 !.Top boundary
    respoi(i,nyf,:)=rhs(i,nyf,:)-(fx*(u(i+1,nyf,:)+u(i-1,nyf,:)) &
      +fy*(rbufN(i,:)+u(i,nyf-1,:))-fc*u(i,nyf,:))
  enddo

  call MPI_WAIT(req(8),stat(:,8),ierr)
  if (iIDph(ng) == iprocsm1) then !.Right most process (along x)
    respoi(nxf,1,:)=rhs(nxf,1,:)-(fx*(Rq*u(nxf,1,:)+Rr*u(nxf-1,1,:)+Rs+u(nxf-1,1,:)) &
      +fy*(u(nxf,2,:)+rbufS(nxf,:))-fc*u(nxf,1,:)) !.SE corner
    do j=2,nyf-1 !.Right boundary
      respoi(nxf,j,:)=rhs(nxf,j,:)-(fx*(Rq*u(nxf,j,:)+Rr*u(nxf-1,j,:)+Rs+u(nxf-1,j,:)) &
        +fy*(u(nxf,j+1,:)+u(nxf,j-1,:))-fc*u(nxf,j,:))
    enddo
    respoi(nxf,nyf,:)=rhs(nxf,nyf,:)-(fx*(Rq*u(nxf,nyf,:)+Rr*u(nxf-1,nyf,:)+Rs+u(nxf-1,nyf,:)) &
      +fy*(rbufN(nxf,:)+u(nxf,nyf-1,:))-fc*u(nxf,nyf,:)) !.NE corner
  else
    respoi(nxf,1,:)=rhs(nxf,1,:)-(fx*(rbufE(1,:)+u(nxf-1,1,:)) &
      +fy*(u(nxf,2,:)+rbufS(nxf,:))-fc*u(nxf,1,:))
    do j=2,nyf-1 !.Right boundary
      respoi(nxf,j,:)=rhs(nxf,j,:)-(fx*(rbufE(j,:)+u(nxf-1,j,:)) &
        +fy*(u(nxf,j+1,:)+u(nxf,j-1,:))-fc*u(nxf,j,:))
    enddo
    respoi(nxf,nyf,:)=rhs(nxf,nyf,:)-(fx*(rbufE(nyf,:)+u(nxf-1,nyf,:)) &
      +fy*(rbufN(nxf,:)+u(nxf,nyf-1,:))-fc*u(nxf,nyf,:)) !.NE corner
  endif
    call MPI_WAITALL(4,(/req(1),req(2),req(3),req(4)/), &
                       (/stat(:,1),stat(:,2),stat(:,3),stat(:,4)/),ierr)

#IF (MGTIMER==1)
    call CPU_TIME(rest2)
    restph=restph+rest2-rest1
#ENDIF
  end function respoi
!*********************************************************
  subroutine MGphi(u,lambda,rho,xBC,ur)
! Multigrid (GAMMA-cycle) algorithm for solution of linear elliptic equation
!                div(lambda*grad(u))=rho.
! Input u(inx,iny,nzL) contains the initial guess to u, while on output it
! contains the final approximation. 
! Right hand side Dirichlet values inputted in urb
  implicit none
  real(DP), intent(inout) :: u(:,:,:)
  real(DP), allocatable, intent(in) :: rho(:,:,:),lambda(:,:,:),ur(:,:)
  integer(I4B), intent(in) :: xBC(:) !.Left & right BCs
  integer(I4B) :: jcycle,ng
#IF (DispRes==1)
  include 'mpif.h'
  real(DP), allocatable :: res(:,:,:)
  real(DP) :: nresL,nres
  integer(I4B) :: i,j,ierr
#ENDIF

  nxf=nxFph(lgph)
  nyf=nyFph(lgph)

#IF (MGTIMER==1)
  resDtph=0.0_dp
  proltph=0.0_dp
  preXtph=0.0_dp
  postXtph=0.0_dp
  rNrtph=0.0_dp
  restph=0.0_dp
  call CPU_TIME(mg1)
#ENDIF

#IF (MGTIMER==1)
  call CPU_TIME(rstdiff1)
#ENDIF
!.Also need to compute coarsened lambda at each grid (Diff)
  Diff(lgph)%a=lambda
  do ng=lgph,2,-1
    nxf=nxFph(ng)
    nyf=nyFph(ng)
    nxc=nxCph(ng)
    nyc=nyCph(ng)
    nxu=nxFph(ng-1)
    nyu=nyFph(ng-1)
    if (acuph(ng,1)) Diff(ng-1)%a=rstrct(Diff(ng)%a,wdcph(ng),(/0,0/), &
        acuph(ng-1,:),urnkph(ng-1)%a,nurnkph(ng-1,:))
  enddo
  call MPI_BARRIER(XYcomm0,jcycle)
#IF (MGTIMER==1)
  call CPU_TIME(rstdiff2)
  resDtph=resDtph+rstdiff2-rstdiff1
#ENDIF

  nxf=nxFph(lgph)
  nyf=nyFph(lgph)

  do jcycle=1,BETph

    call Gcycphi(lgph,u,Diff,rho,xBC,ur) !.G-cycle

#IF (DispRes==1)
    allocate(res(nxf,nyf,nzL))
    iprocsm1=iprocsph(lgph)-1
    jprocsm1=jprocsph(lgph)-1
    res=0.0_dp
    nresL=0.0_dp
    nres=0.0_dp
    res=resphi(lgph,u,Diff(lgph)%a,rho,xBC,ur)
    do i=1,nxf; do j=1,nyf
      nresL=nresL+(res(i,j,1))**2
    enddo; enddo
    call MPI_REDUCE(nresL,nres,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,XYcomm0,ierr)
    nres=nres/dble(nxf*iprocsph(lgph)*nyf*jprocsph(lgph))
    if (wrnk0==0) write(*,"('jcyc:',I4,'  |          nres=',E16.8E3)") jcycle,nres
    deallocate(res)
#ENDIF

  enddo
#IF (MGTIMER==1)
  call CPU_TIME(mg2)
  mgtph=mg2-mg1

  if (wrnk0==0) then
    open(unit=10,file='mgtimes.out',status='unknown',form='unformatted')
    write(10) (/initPMG,mgtph,resDtph,preXtph,postXtph,proltph,rNrtph,restph/)
    close(unit=10)
  endif
#ENDIF

  end subroutine MGphi
!*********************************************************
  subroutine MGpsi(u,rho,iksq,xBC)
!.Multigrid (V-cycle) algorithm for solution of linear elliptic equation
!                u - ksq*(L)u=rho
!.where L is the Laplacian operator. Input u contains the
!.initial guess to the solution u, rho the source RHS. On output u
!.contains the final approximation.
  implicit none
  real(DP), intent(inout) :: u(:,:,:),rho(:,:,:)
  real(DP), intent(in) :: iksq
  integer(I4B), intent(in) :: xBC(:)
  integer(I4B) :: jcycle
#IF (DispRes==1)
  include 'mpif.h'
  real(DP), allocatable :: res(:,:,:)
  real(DP) :: nresL,nres
  integer(I4B) :: i,j,ierr
#ENDIF

  ksq=iksq

  nxf=nxFps(lgps)
  nyf=nyFps(lgps)

#IF (MGTIMER==1)
  resDtps=0.0_dp
  proltps=0.0_dp
  preXtps=0.0_dp
  postXtps=0.0_dp
  rNrtps=0.0_dp
  restps=0.0_dp
  call CPU_TIME(mg1)
#ENDIF

  do jcycle=1,BETps

    call Gcycpsi(lgps,u,rho,xBC) !.G-cycle loop

!..............................................................
#IF (DispRes == 1)
!...PRINT NORM OF RESIDUE FOR TESTING
    allocate(res(nxf,nyf,nzL))
    iprocsm1=iprocsps(lgps)-1
    jprocsm1=jprocsps(lgps)-1
    res=0.0_dp
    nresL=0.0_dp
    nres=0.0_dp
    res=respsi(lgps,u,rho,xBC)
    do i=1,nxf; do j=1,nyf
      nresL=nresL+(res(i,j,1))**2
    enddo; enddo
    call MPI_REDUCE(nresL,nres,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,XYcomm0,ierr)
    nres=nres/dble(nxf*iprocsps(lgps)*nyf*jprocsps(lgps))
    if (wrnk0==0) write(*,"('jcyc:',I4,'  |          nres=',E16.8E3)") jcycle,nres
    deallocate(res)
#ENDIF
!..............................................................

  enddo

#IF (MGTIMER==1)
  call CPU_TIME(mg2)
  mgtps=mg2-mg1

  if (wrnk0==0) then
    open(unit=10,file='mgtimes.out',status='unknown',form='unformatted')
    write(10) (/initPMG,mgtps,resDtps,preXtps,postXtps,proltps,rNrtps,restps/)
    close(unit=10)
  endif
#ENDIF

  end subroutine MGpsi
!*********************************************************
  subroutine MGhd(u,rho,iDiff,xBC,qn)
!.Multigrid (V-cycle) algorithm for solution of linear elliptic equation
!                u + Diff*(L^pd)u=rho
!.where Diff is a constant diffusion coefficient, L the Laplacian
!.operator and pd the number of Laplacians whose composite make the
!.hyper-diffusion operator. Input u(inx,iny) contains the initial guess
!.to the solution u, av the initial guess to auxiliary variable values,
!.rho the source RHS. On output u contains the final approximation.
!.Domain and grid defined by aLx, aLy, inx and iny
  implicit none
  real(DP), intent(inout) :: u(:,:,:)
  real(DP), intent(in) :: iDiff(:),rho(:,:,:)
  integer(I4B), intent(in) :: xBC(:),qn
  integer(I4B) :: jcycle
#IF (HDop==4 .OR. HDop==6)
  real(DP), allocatable :: rhs2(:,:,:), au1(:,:,:) !.au=auxiliary variable
#IF (HDop==6)
  real(DP), allocatable :: au2(:,:,:) !.second auxiliary variable
#ENDIF
#ENDIF
#IF (DispRes == 1)
  include 'mpif.h'
  real(DP), allocatable :: res1(:,:,:)
#IF (HDop==4)
  real(DP), allocatable :: res2(:,:,:)
#ELIF (HDop==6)
  real(DP), allocatable :: res2(:,:,:),res3(:,:,:)
#ENDIF
  real(DP) :: nresL,nres
  integer(I4B) :: i,j,ierr
#ENDIF

  hDiff=iDiff

  nxf=nxFhd(qn)%a(lghd(qn))
  nyf=nyFhd(qn)%a(lghd(qn))

#IF (HDop==4 .OR. HDop==6)
  allocate(rhs2(nxf,nyf,nzL),au1(nxf,nyf,nzL))
  rhs2=0.0_dp
  au1=0.0_dp
#IF (HDop==6)
  allocate(au2(nxf,nyf,nzL))
  au2=0.0_dp
#ENDIF
#ENDIF

#IF (MGTIMER==1)
  resDthd=0.0_dp
  prolthd=0.0_dp
  preXthd=0.0_dp
  postXthd=0.0_dp
  rNrthd=0.0_dp
  resthd=0.0_dp
  call CPU_TIME(mg1)
#ENDIF
  do jcycle=1,BEThd(qn) !.G-cycle(s)

#IF (HDop==4)
    call Gcychd(lghd(qn),au1,u,rho,rhs2,xBC,qn) !.Split Gamma-cycle loop
#ELIF (HDop==6)
    call Gcychd(lghd(qn),au1,au2,u,rho,rhs2,rhs2,xBC,qn) !.Split Gamma-cycle loop
#ELSE
!...Bi/TriHelmholtz (not split) Gamma-cycle or Nonisotropic diffusion
    call Gcychd(lghd(qn),u,rho,xBC,qn)
#ENDIF

!..............................................................
#IF (DispRes == 1)
!...PRINT NORM OF RESIDUE FOR TESTING
    nxf=nxFhd(qn)%a(lghd(qn))
    nyf=nyFhd(qn)%a(lghd(qn))
    allocate(res1(nxf,nyf,nzL))
    iprocsm1=iprocshd(qn)%a(lghd(qn))-1
    jprocsm1=jprocshd(qn)%a(lghd(qn))-1
    res1=0.0_dp
    nresL=0.0_dp
    nres=0.0_dp
#IF (HDop==4)
    allocate(res2(nxf,nyf,nzL))
    res2=0.0_dp
    call reshd(lghd(qn),res1,res2,au1,u,rho,rhs2,xBC,qn)
#ELIF (HDop==6)
    allocate(res2(nxf,nyf,nzL),res3(nxf,nyf,nzL))
    res2=0.0_dp
    res3=0.0_dp
    call reshd(lghd(qn),res1,res2,res3,au1,au2,u,rho,rhs2,rhs2,xBC,qn)
#ELSE
    res1=reshd(lghd(qn),u,rho,xBC,qn)
#ENDIF
    do i=1,nxf; do j=1,nyf
      nresL=nresL+(res1(i,j,1))**2
    enddo; enddo
    call MPI_REDUCE(nresL,nres,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,XYcomm0,ierr)
    nres=nres/dble(nxf*iprocshd(qn)%a(lghd(qn))*nyf*jprocshd(qn)%a(lghd(qn)))
    if (wrnk0==0) write(*,"('jcyc:',I4,'  |          nres=',E16.8E3)") jcycle,nres
    deallocate(res1)
#IF (HDop==4)
    deallocate(res2)
#ELIF (HDop==6)
    deallocate(res2,res3)
#ENDIF
#ENDIF
!..............................................................

  enddo

#IF (MGTIMER==1)
  call CPU_TIME(mg2)
  mgthd=mg2-mg1

  if (wrnk0==0) then
    open(unit=10,file='mgtimes.out',status='unknown',form='unformatted')
    write(10) (/initPMG,mgthd,resDthd,preXthd,postXthd,prolthd,rNrthd,resthd/)
    close(unit=10)
  endif
#ENDIF

  end subroutine MGhd
!*********************************************************
  subroutine FMGphi(u,lambda,rhs,xBC,urb)
! Full Multigrid Algorithm for solution of linear elliptic equation
!                div(lambda*grad(u))=rhs
! u(inx,iny,nzL) is filled on output with the final approximation.
! Right hand side Dirichlet values inputted in urb
  implicit none
  real(DP), intent(inout) :: u(:,:,:)
  real(DP), intent(in) :: rhs(:,:,:),urb(:,:),lambda(:,:,:)
!.Left and right boundary conditions
  integer(I4B), intent(in) :: xBC(:)
  integer(I4B) :: jcycle,ng
  type(allo3d), allocatable :: rho(:)
!.Array for RHS Dirichlet values at each grid.
  type(allo2d), allocatable :: ur(:)
!.These are pointers to current (u_n) and previous (u_{n-1}) estimate
  real(DP), pointer :: un(:,:,:),un1(:,:,:)
  real(DP), allocatable :: tmp(:,:,:),tmpc(:,:,:)
#IF (DispRes == 1)
  include 'mpif.h'
  real(DP), allocatable :: res(:,:,:)
  real(DP) :: nresL,nres
  integer(I4B) :: i,j,ierr
#ENDIF

  nxf=nxFph(lgph)
  nyf=nyFph(lgph)
  allocate(rho(lgph))
  allocate(rho(lgph)%a(nxf,nyf,nzL))  !.Allocate storage for RHS on grid ng
  rho(lgph)%a=rhs        !.fill it with the input RHS

  allocate(ur(lgph))
  allocate(ur(lgph)%a(nyf,nzL))  !.Dirichlet boundary values on grid ng
  ur(lgph)%a=urb        !.fill it with u(x=Lx) on finest grid

  do ng=lgph,2,-1
!...These two are not needed with simplest coarsening stencil
!    iprocsm1=iprocsph(ng)-1
!    jprocsm1=jprocsph(ng)-1
!...Extract fine and coarse grids dimensions from nxv, nyv
    nxf=nxFph(ng)
    nyf=nyFph(ng)
    nxc=nxCph(ng)
    nyc=nyCph(ng)
    nxu=nxFph(ng-1)
    nyu=nyFph(ng-1)
!...Full Multi-grid needs the source RHS at every grid. Restrict and save
!...IMPORTANT: check the symmetries here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(rho(ng-1)%a(nxu,nyu,nzL))
    if (acuph(ng,1)) rho(ng-1)%a=rstrct(rho(ng)%a,wdcph(ng),(/0,0/), &
      acuph(ng-1,:),urnkph(ng-1)%a,nurnkph(ng-1,:))

    allocate(ur(ng-1)%a(nyu,nzL)) !.Only restrict ur during y-coarsening
!...THIS WAY OF DOING IT WASTES MEMORY. Saving copies of ur.
    if (wdcph(ng)==0 .OR. wdcph(ng)==2) then
      allocate(tmp(nxf,nyf,nzL),tmpc(nxu,nyu,nzL))
      tmp=0.0_dp
      tmp(1,:,:)=ur(ng)%a
      tmp(2,:,:)=ur(ng)%a
      if (acuph(ng,1)) then
!.......This one should use 2 instead of wdcph(ng), but current domain
!.......unification code doesn't allow it
        tmpc=rstrct(tmp,wdcph(ng),(/0,0/),acuph(ng-1,:), &
          urnkph(ng-1)%a,nurnkph(ng-1,:),ur(ng-1)%a)
        ur(ng-1)%a=tmpc(1,:,:)
      endif
      deallocate(tmp,tmpc)
    else
      if (acuph(ng,1)) ur(ng-1)%a=ur(ng)%a
    endif
  enddo

!.Also need to compute coarsened lambda at each grid (Diff)
  Diff(lgph)%a=lambda
  do ng=lgph,2,-1
!...These two are not needed with simplest coarsening stencil
!    iprocsm1=iprocsph(ng)-1
!    jprocsm1=jprocsph(ng)-1
    nxf=nxFph(ng)
    nyf=nyFph(ng)
    nxc=nxCph(ng)
    nyc=nyCph(ng)
    nxu=nxFph(ng-1)
    nyu=nyFph(ng-1)
    if (acuph(ng,1)) Diff(ng-1)%a=rstrct(Diff(ng)%a,wdcph(ng),(/0,0/), &
        acuph(ng-1,:),urnkph(ng-1)%a,nurnkph(ng-1,:))
  enddo
  call MPI_BARRIER(XYcomm0,ng)

  iprocsm1=iprocsph(1)-1
  jprocsm1=jprocsph(1)-1
  nxf=nxFph(1)
  nyf=nyFph(1)
  allocate(un(nxf,nyf,nzL))
  un=0.0_dp !.Initial guess

  if (acuph(1,1)) then
    do jcycle=1,10 !NU2ph
      call relaxphi(1,un,Diff(1)%a,rho(1)%a,xBC,ur(1)%a)
    enddo
  endif

  do ng=2,lgph

    iprocsm1=iprocsph(ng)-1
    jprocsm1=jprocsph(ng)-1
!...New fine grid is twice of last coarse grid (minus 1 for x due to wall point)
    nxf=nxFph(ng)
    nyf=nyFph(ng)
    nxc=nxCph(ng)
    nyc=nyCph(ng)
!...Associate pointer of previous estimate (u_{n-1}) to current estimate (u_n)
    un1=>un
    allocate(un(nxf,nyf,nzL))
!...from grid ng-1 to next finer grid
    if (acuph(ng,1)) & 
      un=prolng(un1,wdcph(ng),xBC,COMMph(ng,:),iIDph(ng),jprocsph(ng), &
        nborph(ng,:),acuph(ng-1,:),urnkph(ng-1)%a,nurnkph(ng-1,:),ur(ng-1)%a)
    deallocate(un1)

    if (acuph(ng,1)) then
      do jcycle=1,BETph

        call Gcycphi(ng,un,Diff,rho(ng)%a,xBC,ur(ng)%a) !.G-cycle

!..............................................................
#IF (DispRes == 1)
!.......PRINT NORM OF RESIDUE FOR TESTING
        if (ng==lgph) then
          allocate(res(nxf,nyf,nzL))
          iprocsm1=iprocsph(lgph)-1
          jprocsm1=jprocsph(lgph)-1
          res=0.0_dp
          nresL=0.0_dp
          nres=0.0_dp
          res=resphi(ng,un,Diff(ng)%a,rho(ng)%a,xBC,ur(ng)%a)
          do i=1,nxf; do j=1,nyf
            nresL=nresL+(res(i,j,1))**2
          enddo; enddo
          call MPI_REDUCE(nresL,nres,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,XYcomm0,ierr)
          nres=nres/dble(nxf*iprocsph(lgph)*nyf*jprocsph(lgph))
          if (wrnk0==0) write(*,"('jcyc:',I4,'  |          nres=',E16.8E3)") jcycle,nres
          deallocate(res)
        endif
#ENDIF
!..............................................................

      enddo
    endif

  enddo
  u=un

  deallocate(un)
  do ng=1,lgph
    deallocate(rho(ng)%a,ur(ng)%a)
  enddo
  deallocate(rho,ur)

  end subroutine FMGphi
!*********************************************************
  subroutine MGpoi(u,rho,xBC)
!.Multigrid (GAMMA-cycle) algorithm for solution of Poisson equation
!.               nabla^2(u)=rho.
!.Input u(inx,iny,nzL) contains the initial guess to u, while on
!.output it contains the final approximation.
  implicit none
  real(DP), intent(inout) :: u(:,:,:)
  real(DP), allocatable, intent(in) :: rho(:,:,:)
  integer(I4B), intent(in) :: xBC(:) !.Left & right BCs
  integer(I4B) :: jcycle,ng
#IF (DispRes==1)
  include 'mpif.h'
  real(DP), allocatable :: res(:,:,:)
  real(DP) :: nresL,nres
  integer(I4B) :: i,j,ierr
#ENDIF

  nxf=nxFph(lgph)
  nyf=nyFph(lgph)

#IF (MGTIMER==1)
  resDtpo=0.0_dp
  proltpo=0.0_dp
  preXtpo=0.0_dp
  postXtpo=0.0_dp
  rNrtpo=0.0_dp
  restpo=0.0_dp
  call CPU_TIME(mg1)
#ENDIF

  nxf=nxFph(lgph)
  nyf=nyFph(lgph)

  do jcycle=1,BETph

    call Gcycpoi(lgph,u,rho,xBC) !.G-cycle

#IF (DispRes==1)
    allocate(res(nxf,nyf,nzL))
    iprocsm1=iprocsph(lgph)-1
    jprocsm1=jprocsph(lgph)-1
    res=0.0_dp
    nresL=0.0_dp
    nres=0.0_dp
    res=respoi(lgph,u,rho,xBC)
    do i=1,nxf; do j=1,nyf
      nresL=nresL+(res(i,j,1))**2
    enddo; enddo
    call MPI_REDUCE(nresL,nres,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,XYcomm0,ierr)
    nres=nres/dble(nxf*iprocsph(lgph)*nyf*jprocsph(lgph))
    if (wrnk0==0) write(*,"('jcyc:',I4,'  |          nres=',E16.8E3)") jcycle,nres
    deallocate(res)
#ENDIF

  enddo
#IF (MGTIMER==1)
  call CPU_TIME(mg2)
  mgtpo=mg2-mg1

  if (wrnk0==0) then
    open(unit=10,file='mgtimes.out',status='unknown',form='unformatted')
    write(10) (/initPMG,mgtpo,resDtpo,preXtpo,postXtpo,proltpo,rNrtpo,restpo/)
    close(unit=10)
  endif
#ENDIF

  end subroutine MGpoi
!*********************************************************
  subroutine FMGpoi(u,rhs,xBC)
! Full Multigrid Algorithm for solution of Poisson equation
!                nabla^2(u)=rhs
! u(inx,iny,nzL) is filled on output with the final approximation.
  implicit none
  real(DP), intent(inout) :: u(:,:,:)
  real(DP), intent(in) :: rhs(:,:,:)
!.Left and right boundary conditions
  integer(I4B), intent(in) :: xBC(:)
  integer(I4B) :: jcycle,ng
  type(allo3d), allocatable :: rho(:)
!.These are pointers to current (u_n) and previous (u_{n-1}) estimate
  real(DP), pointer :: un(:,:,:),un1(:,:,:)
#IF (DispRes == 1)
  include 'mpif.h'
  real(DP), allocatable :: res(:,:,:)
  real(DP) :: nres,nresL
  integer(I4B) :: i,j,ierr
#ENDIF

  nxf=nxFph(lgph)
  nyf=nyFph(lgph)
  allocate(rho(lgph))
  allocate(rho(lgph)%a(nxf,nyf,nzL))  !.Allocate storage for RHS on grid ng
  rho(lgph)%a=rhs        !.fill it with the input RHS

  do ng=lgph,2,-1
!...These two are not needed with simplest coarsening stencil
!    iprocsm1=iprocsph(ng)-1
!    jprocsm1=jprocsph(ng)-1
!...Extract fine and coarse grids dimensions from nxv, nyv
    nxf=nxFph(ng)
    nyf=nyFph(ng)
    nxc=nxCph(ng)
    nyc=nyCph(ng)
    nxu=nxFph(ng-1)
    nyu=nyFph(ng-1)
!...Full Multi-grid needs the source RHS at every grid. Restrict and save
!...IMPORTANT: check the symmetries here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate(rho(ng-1)%a(nxu,nyu,nzL))
    if (acuph(ng,1)) rho(ng-1)%a=rstrct(rho(ng)%a,wdcph(ng),(/0,0/), &
      acuph(ng-1,:),urnkph(ng-1)%a,nurnkph(ng-1,:))
  enddo
  call MPI_BARRIER(XYcomm0,ng)

  iprocsm1=iprocsph(1)-1
  jprocsm1=jprocsph(1)-1
  nxf=nxFph(1)
  nyf=nyFph(1)
  allocate(un(nxf,nyf,nzL))
  un=0.0_dp !.Initial guess

  if (acuph(1,1)) then
    do jcycle=1,10 !NU2ph
      call relaxpoi(1,un,rho(1)%a,xBC)
    enddo
  endif

  do ng=2,lgph

    iprocsm1=iprocsph(ng)-1
    jprocsm1=jprocsph(ng)-1
!...New fine grid is twice of last coarse grid (minus 1 for x due to wall point)
    nxf=nxFph(ng)
    nyf=nyFph(ng)
    nxc=nxCph(ng)
    nyc=nyCph(ng)
!...Associate pointer of previous estimate (u_{n-1}) to current estimate (u_n)
    un1=>un
    allocate(un(nxf,nyf,nzL))
!...from grid ng-1 to next finer grid
    if (acuph(ng,1)) &
      un=prolng(un1,wdcph(ng),xBC,COMMph(ng,:),iIDph(ng),jprocsph(ng), &
        nborph(ng,:),acuph(ng-1,:),urnkph(ng-1)%a,nurnkph(ng-1,:))
    deallocate(un1)

    if (acuph(ng,1)) then
      do jcycle=1,BETph

        call Gcycpoi(ng,un,rho(ng)%a,xBC) !.G-cycle

!..............................................................
#IF (DispRes == 1)
!.......PRINT NORM OF RESIDUE FOR TESTING
        if (ng==lgph) then
          allocate(res(nxf,nyf,nzL))
          iprocsm1=iprocsph(lgph)-1
          jprocsm1=jprocsph(lgph)-1
          res=0.0_dp
          nresL=0.0_dp
          nres=0.0_dp
          res=respoi(ng,un,rho(ng)%a,xBC)
          do i=1,nxf; do j=1,nyf
            nresL=nresL+(res(i,j,1))**2
          enddo; enddo
          call MPI_REDUCE(nresL,nres,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,XYcomm0,ierr)
          nres=nres/dble(nxf*iprocsph(lgph)*nyf*jprocsph(lgph))
          if (wrnk0==0) write(*,"('jcyc:',I4,'  |          nres=',E16.8E3)") jcycle,nres
          deallocate(res)
        endif
#ENDIF
!..............................................................

      enddo
    endif

  enddo
  u=un

  deallocate(un)
  do ng=1,lgph
    deallocate(rho(ng)%a)
  enddo
  deallocate(rho)

  end subroutine FMGpoi
!*********************************************************
  function rstrct(uf,ldc,xBC,acu,urnk,nurnk,ufr)
!.Restrict the fine-grid solution uf(nxf,nyf,nzL) and output
!.the coarse-grid solution rstrct(nxc,nyc,nzL).
!.Unify SubDomains because nxc<nxmin or nyc<nymin
!.so that the volume to surface area of each process remains
!.sufficiently high for good parallel performance
!.Other inputs:
! ldc: last dimension coarsened. If =0 restrict along both dimensions,
!      if =1 restrict along x, if =2 restrict along y
! (l/r)sym: left/right BC for uf, =1 even, =-1 odd, =0 linear extrapolation
! ufr: u(x=Lx) is the Dirichlet value at right boundary
  implicit none
  include 'mpif.h'
  real(DP), allocatable, intent(in) :: uf(:,:,:)
  real(DP), intent(in), optional :: ufr(:,:)
  integer(I4B), intent(in) :: ldc,xBC(:),urnk(:,:),nurnk(:)
  logical, intent(in) :: acu(:)
  real(DP), allocatable :: rstrct(:,:,:),rst(:,:,:)
  real(DP), allocatable :: sbuf(:,:,:),rbuf(:,:,:)
  integer(I4B) :: ic,iff,jc,jf,k,iL,iU,jL,jU,req,ierr,cof
  integer(I4B) :: stat(MPI_STATUS_SIZE)

  allocate(rstrct(nxu,nyu,nzL),rst(nxc,nyc,nzL))
  cof=maxval((/ nxf/nxc,nyf/nyc /))

!.Restrict interior points
  if (ldc==0) then
!.Restrict along both directions

#IF (IR==0)
!...4-point averaging
    do jc=2,nyc-1
      jf=2*jc
      do ic=2,nxc-1
        iff=2*ic
        rst(ic,jc,:)=0.25_dp*(uf(iff,jf,:)+uf(iff,jf-1,:) &
                             +uf(iff-1,jf,:)+uf(iff-1,jf-1,:))
      enddo
    enddo
#ELIF (IR==1)
!...Kwak restriction
#ELIF (IR==2)
!...Bilinear averaging (1/64=.015625)
#ELIF (IR==3)
!...Cubic interpolation (1/256=.00390625)
#ENDIF

  elseif (ldc==1) then
!.Restrict along x (nyc=nyf) with linear weighting

#IF (ISR==1) 
!...Linear weighting
    if (cof==2) then
      do ic=2,nxc-1
        iff=2*ic
        rst(ic,:,:)=0.5_dp*(uf(iff,:,:)+uf(iff-1,:,:))
      enddo
    else !if (cof==4) then
      do ic=2,nxc-1
        iff=4*ic
        rst(ic,:,:)=0.25_dp*(uf(iff,:,:)+uf(iff-1,:,:) &
                            +uf(iff-2,:,:)+uf(iff-3,:,:))
!        rst(ic,:,:)=(uf(iff,:,:)+2.0_dp*(uf(iff-1,:,:) &
!                    +uf(iff-2,:,:))+uf(iff-3,:,:))/6.0_dp
!        rst(ic,:,:)=(-uf(iff,:,:)+9.0_dp*(uf(iff-1,:,:) &
!                     +uf(iff-2,:,:))-uf(iff-3,:,:))/16.0_dp
      enddo
    endif
#ELIF (ISR==2) 
!...Cubic interpolation
#ENDIF

  else !if (ldc==2) then
!.Restrict along y (nxc=nxf)

#IF (ISR==1) 
!...Linear weighting
    if (cof==2) then
      do jc=2,nyc-1
        jf=2*jc
        rst(:,jc,:)=0.5_dp*(uf(:,jf,:)+uf(:,jf-1,:))
      enddo
    else !if (cof==4) then
      do jc=2,nyc-1
        jf=4*jc
        rst(:,jc,:)=0.25_dp*(uf(:,jf,:)+uf(:,jf-1,:) &
                            +uf(:,jf-2,:)+uf(:,jf-3,:))
!        rst(:,jc,:)=(uf(:,jf,:)+2.0_dp*(uf(:,jf-1,:) &
!                    +uf(:,jf-2,:))+uf(:,jf-3,:))/6.0_dp
!        rst(:,jc,:)=(-uf(:,jf,:)+9.0_dp*(uf(:,jf-1,:) &
!                     +uf(:,jf-2,:))-uf(:,jf-3,:))/16.0_dp
      enddo
    endif
#ELIF (ISR==2) 
!...Cubic interpolation
#ENDIF
  endif
  
!.Restrict points closest to x boundaries separately
  if (ldc==0) then
!.Restrict along both directions

#IF (IRES==0)
!...4-point averaging
    do ic=2,nxc-1
      iff=2*ic
      rst(ic,1,:)=0.25_dp*(uf(iff,2,:)+uf(iff,1,:) &
                          +uf(iff-1,2,:)+uf(iff-1,1,:))
      rst(ic,nyc,:)=0.25_dp*(uf(iff,nyf,:)+uf(iff,nyf-1,:) &
                            +uf(iff-1,nyf,:)+uf(iff-1,nyf-1,:))
    enddo
    do jc=1,nyc
      jf=2*jc
      rst(1,jc,:)=0.25_dp*(uf(2,jf,:)+uf(2,jf-1,:) &
                          +uf(1,jf,:)+uf(1,jf-1,:))
      rst(nxc,jc,:)=0.25_dp*(uf(nxf,jf,:)+uf(nxf,jf-1,:) &
                            +uf(nxf-1,jf,:)+uf(nxf-1,jf-1,:))
    enddo
#ELIF (IRES==1)
!...Kwak restriction
#ELIF (IR==2)
!...Bilinear averaging (1/64=.015625)
#ELIF (IR==3)
!...Cubic interpolation (1/256=.00390625)
#ENDIF

  elseif (ldc==1) then
!.Restrict along x (nyc=nyf)

#IF (ISR==1)
!...Linear weighting
    if (cof==2) then
      rst(1,:,:)=0.5_dp*(uf(2,:,:)+uf(1,:,:))
      rst(nxc,:,:)=0.5_dp*(uf(nxf,:,:)+uf(nxf-1,:,:))
    else !if (cof==4) then
      rst(1,:,:)=0.25_dp*(uf(4,:,:)+uf(3,:,:)+uf(2,:,:)+uf(1,:,:))
      rst(nxc,:,:)=0.25_dp*(uf(nxf,:,:)+uf(nxf-1,:,:)+uf(nxf-2,:,:)+uf(nxf-3,:,:))
!      rst(1,:,:)=(uf(4,:,:)+2.0_dp*(uf(3,:,:)+uf(2,:,:))+uf(1,:,:))/6.0_dp
!      rst(nxc,:,:)=(uf(nxf,:,:)+2.0_dp*(uf(nxf-1,:,:) &
!                   +uf(nxf-2,:,:))+uf(nxf-3,:,:))/6.0_dp
!      rst(1,:,:)=(-uf(4,:,:)+9.0_dp*(uf(3,:,:)+uf(2,:,:))-uf(1,:,:))/16.0_dp
!      rst(nxc,:,:)=(-uf(nxf,:,:)+9.0_dp*(uf(nxf-1,:,:) &
!                    +uf(nxf-2,:,:))-uf(nxf-3,:,:))/16.0_dp
    endif
#ELIF (ISR==2)
!...Cubic interpolation
#ENDIF

  else !if (ldc==2) then
!.Restrict along y (nxc=nxf)

#IF (ISR==1)
!...Linear weighting
    if (cof==2) then
      rst(:,1,:)=0.5_dp*(uf(:,2,:)+uf(:,1,:))
      rst(:,nyc,:)=0.5_dp*(uf(:,nyf,:)+uf(:,nyf-1,:))
    else !if (cof==4) then
      rst(:,1,:)=0.25_dp*(uf(:,4,:)+uf(:,3,:)+uf(:,2,:)+uf(:,1,:))
      rst(:,nyc,:)=0.25_dp*(uf(:,nyf,:)+uf(:,nyf-1,:) &
                           +uf(:,nyf-2,:)+uf(:,nyf-3,:))
!      rst(:,1,:)=(uf(:,4,:)+2.0_dp*(uf(:,3,:)+uf(:,2,:))+uf(:,1,:))/6.0_dp
!      rst(:,nyc,:)=(uf(:,nyf,:)+2.0_dp*(uf(:,nyf-1,:) &
!                   +uf(:,nyf-2,:))+uf(:,nyf-3,:))/6.0_dp
!      rst(:,1,:)=(-uf(:,4,:)+9.0_dp*(uf(:,3,:) &
!                  +uf(:,2,:))-uf(:,1,:))/16.0_dp
!      rst(:,nyc,:)=(-uf(:,nyf,:)+9.0_dp*(uf(:,nyf-1,:) &
!                    +uf(:,nyf-2,:))-uf(:,nyf-3,:))/16.0_dp
    endif
#ELIF (ISR==2)
!...Cubic interpolation
#ENDIF

  endif

  if (acu(2)) then
    if (acu(1)) then
!.....These processes receive data from others
      allocate(rbuf(nxc,nyc,nzL))
 
      do jc=1,nurnk(2)
        do ic=1,nurnk(1)
!.........location (in uc) of received data from other procs in unification
          iL=(ic-1)*nxc+1
          iU=ic*nxc
          jL=(jc-1)*nyc+1
          jU=jc*nyc
          if (urnk(ic,jc)==myrank0) then 
            rstrct(iL:iU,jL:jU,:)=rst
          else
!            call MPI_IRECV(rbuf,nxc*nyc*nzL,MPI_DOUBLE_PRECISION, &
!                           urnk(ic,jc),0,XYcomm0,req,ierr)
!            call MPI_WAIT(req,stat,ierr)
            call MPI_RECV(rbuf,nxc*nyc*nzL,MPI_DOUBLE_PRECISION, &
                          urnk(ic,jc),0,XYcomm0,stat,ierr)
            rstrct(iL:iU,jL:jU,:)=rbuf
          endif
        enddo
      enddo
    else
!.....Send data to process closer to the center of global domain
      allocate(sbuf(nxc,nyc,nzL))
      sbuf=rst
!      call MPI_ISEND(sbuf,nxc*nyc*nzL,MPI_DOUBLE_PRECISION, &
!                     urnk(1,1),0,XYcomm0,req,ierr)
!      call MPI_WAIT(req,stat,ierr)
      call MPI_SEND(sbuf,nxc*nyc*nzL,MPI_DOUBLE_PRECISION, &
                    urnk(1,1),0,XYcomm0,ierr)
    endif
  else
    rstrct=rst
  endif

  end function rstrct
!*********************************************************
  function prolng(ucU,ldc,xBC,COMMs,iID,jprocs,nbors,acu,urnk,nurnk,ucr)
!.Coarse-to-fine prolongation. The coarse-grid quantity is
!.input as uc, the fine-grid equivalent is returned as prolng.
  implicit none
  include 'mpif.h'
  real(DP), intent(in) :: ucU(:,:,:)
  integer(I4B), intent(in) :: ldc,xBC(:),iID,jprocs,urnk(:,:), &
                              nurnk(:),COMMs(:),nbors(:)
  logical, intent(in) :: acu(:)
  real(DP), intent(in), optional :: ucr(:,:)
  real(DP), allocatable :: prolng(:,:,:),uc(:,:,:)
  real(DP), allocatable :: sbuf(:,:,:),rbuf(:,:,:)
  integer(I4B) :: ic,iff,jc,jf,k,k1,k2,iL,iU,jL,jU,ierr,cof,p,q
  integer(I4B) :: Nrank,Srank,Erank,Wrank
  integer(I4B), allocatable :: Drank(:)
  real(DP), allocatable :: Lsu(:),Rsu(:),loc1D(:)
  real(DP), allocatable, dimension(:,:) :: rbufN,sbufN,rbufE,sbufE, &
                                           rbufS,sbufS,rbufW,sbufW
  real(DP), allocatable, dimension(:) :: rNE,rSE,rSW,rNW,sNE,sSE,sSW,sNW
  real(DP) :: Lqu,Rqu,Lru,Rru
  integer(I4B) :: req(16),stat(MPI_STATUS_SIZE,16)

  Nrank=nbors(1)
  Erank=nbors(3)
  Srank=nbors(5)
  Wrank=nbors(7)
  allocate(Drank(4))
  Drank=(/nbors(2),nbors(4),nbors(6),nbors(8)/)

  allocate(uc(nxc,nyc,nzL),prolng(nxf,nyf,nzL))
  cof=maxval((/ nxf/nxc,nyf/nyc /))

!.Re-distribute subdomains if the prolongated grid has
!.nxf>=nxmin or nyf>=nymin
  if (acu(2)) then
    if (acu(1)) then
!.....These processes send part of their data to others
      allocate(sbuf(nxc,nyc,nzL))
      do jc=1,nurnk(2)
        do ic=1,nurnk(1)
!.........location (in uc) of data to be distributed to other procs
          iL=(ic-1)*nxc+1
          iU=ic*nxc
          jL=(jc-1)*nyc+1
          jU=jc*nyc
          if (urnk(ic,jc)==myrank0) then 
            uc=ucU(iL:iU,jL:jU,:)
          else
            sbuf=ucU(iL:iU,jL:jU,:)
!.For some reason SEND/RECV interfere with MPI_WAITs below
!            call MPI_SEND(sbuf,nxc*nyc*nzL,MPI_DOUBLE_PRECISION, &
!                          urnk(ic,jc),0,XYcomm0,ierr)
            call MPI_ISEND(sbuf,nxc*nyc*nzL,MPI_DOUBLE_PRECISION, &
                           urnk(ic,jc),0,XYcomm0,req(1),ierr)
            call MPI_WAIT(req(1),stat(:,1),ierr)
          endif
        enddo
      enddo
    else
!.....Receive data from process closer to the center of global domain
      allocate(rbuf(nxc,nyc,nzL))
!.For some reason SEND/RECV interfere with MPI_WAITs below
!      call MPI_RECV(rbuf,nxc*nyc*nzL,MPI_DOUBLE_PRECISION, &
!                    urnk(1,1),0,XYcomm0,stat(:,1),ierr)
      call MPI_IRECV(rbuf,nxc*nyc*nzL,MPI_DOUBLE_PRECISION, &
                     urnk(1,1),0,XYcomm0,req(1),ierr)
      call MPI_WAIT(req(1),stat(:,1),ierr)
      uc=rbuf
    endif
  else
    uc=ucU
  endif

!.First choose the value of certain factors in the stencil/matrix
!.based on the input boundary conditions
!......... BCs of u ................................!
!.LHS BC of u
  allocate(Lsu(nzL),loc1D(nzL))
  if (xBC(1)==0) then
!...linear extrapolation w/o boundary value
    Lqu=2.0_dp
    Lru=-1.0_dp
    Lsu=0.0_dp
  elseif (xBC(1)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Lqu=-1.0_dp
    Lru=0.0_dp
    forall(k=1:nzL) loc1D(k)=sum(uc(1,:,k))
    call MPI_ALLreduce(loc1D,Lsu,nzL,MPI_DOUBLE_PRECISION,MPI_SUM,COMMs(1),ierr)
    Lsu=2.0_dp*Lsu/dble(nyc*jprocs)
  else ! if (xBC(1)==1 .OR. xBC(1)==-1) then ! symmetry
    Lqu=dble(xBC(1))
    Lru=0.0_dp
    Lsu=0.0_dp
  endif
!.RHS BC of u
  allocate(Rsu(nzL))
  if (xBC(2)==0) then
!...linear extrapolation w/o boundary value
    Rqu=2.0_dp
    Rru=-1.0_dp
    Rsu=0.0_dp
  elseif (xBC(2)==2) then
!...ky=0 component has even symmetry, ky.neq.0 component is odd
    Rqu=-1.0_dp
    Rru=0.0_dp
    forall(k=1:nzL) loc1D(k)=sum(uc(nxc,:,k))
    call MPI_ALLreduce(loc1D,Rsu,nzL,MPI_DOUBLE_PRECISION,MPI_SUM,COMMs(1),ierr)
    Rsu=2.0_dp*Rsu/dble(nyc*jprocs)
  else ! if (xBC(2)==1 .OR. xBC(2)==-1) then ! symmetry
    Rqu=dble(xBC(2))
    Rru=0.0_dp
    Rsu=0.0_dp
  endif

  if (ldc==0) then
!.Restrict along both directions

!..Need to communicate with eight neighboring procs
   allocate(rbufN(nxc,nzL),sbufN(nxc,nzL),rbufS(nxc,nzL),sbufS(nxc,nzL), &
            rbufE(nYc,nzL),sbufE(nyc,nzL),rbufW(nyc,nzL),sbufW(nyc,nzL))
   allocate(sNE(nzL),sSE(nzL),sSW(nzL),sNW(nzL),rNE(nzL),rSE(nzL),rSW(nzL),rNW(nzL))
!..Send top row to North rank
   sbufN=uc(:,nyc,:)
   call MPI_ISEND(sbufN,nxc*nzL,MPI_DOUBLE_PRECISION,Nrank,0,COMMs(1),req(1),ierr)
!..Receive uc row from South rank, for uc(i,j-1)
   call MPI_IRECV(rbufS,nxc*nzL,MPI_DOUBLE_PRECISION,Srank,0,COMMs(1),req(2),ierr)
!..Send bottom row to South rank
   sbufS=uc(:,1,:)
   call MPI_ISEND(sbufS,nxc*nzL,MPI_DOUBLE_PRECISION,Srank,1,COMMs(1),req(3),ierr)
!..Receive uc row from North rank, for uc(i,j+1)
   call MPI_IRECV(rbufN,nxc*nzL,MPI_DOUBLE_PRECISION,Nrank,1,COMMs(1),req(4),ierr)
!..Send right column to East rank
   sbufE=uc(nxc,:,:)
   call MPI_ISEND(sbufE,nyc*nzL,MPI_DOUBLE_PRECISION,Erank,2,COMMs(2),req(5),ierr)
!..Receive uc column from West rank, for uc(i-1,j)
   call MPI_IRECV(rbufW,nyc*nzL,MPI_DOUBLE_PRECISION,Wrank,2,COMMs(2),req(6),ierr)
!..Send top-right corner value to NE rank
   sNE=uc(nxc,nyc,:)
   call MPI_ISEND(sNE,nzL,MPI_DOUBLE_PRECISION,Drank(1),4,COMMs(3),req(7),ierr)
!..Receive uc value from SW rank
   call MPI_IRECV(rSW,nzL,MPI_DOUBLE_PRECISION,Drank(3),4,COMMs(3),req(8),ierr)
!..Send bottom-right corner value to SE rank
   sSE=uc(nxc,1,:)
   call MPI_ISEND(sSE,nzL,MPI_DOUBLE_PRECISION,Drank(2),5,COMMs(3),req(9),ierr)
!..Receive uc value from NW rank
   call MPI_IRECV(rNW,nzL,MPI_DOUBLE_PRECISION,Drank(4),5,COMMs(3),req(10),ierr)
!..Send left column to West rank
   sbufW=uc(1,:,:)
   call MPI_ISEND(sbufW,nyc*nzL,MPI_DOUBLE_PRECISION,Wrank,3,COMMs(2),req(11),ierr)
!..Receive uc column from East rank, for uc(i+1,j)
   call MPI_IRECV(rbufE,nyc*nzL,MPI_DOUBLE_PRECISION,Erank,3,COMMs(2),req(12),ierr)
!..Send top-left corner value to NW rank
   sNW=uc(1,nyc,:)
   call MPI_ISEND(sNW,nzL,MPI_DOUBLE_PRECISION,Drank(4),7,COMMs(3),req(13),ierr)
!..Receive uc value from SE rank
   call MPI_IRECV(rSE,nzL,MPI_DOUBLE_PRECISION,Drank(2),7,COMMs(3),req(14),ierr)
!..Send bottom-left corner value to SW rank
   sSW=uc(1,1,:)
   call MPI_ISEND(sSW,nzL,MPI_DOUBLE_PRECISION,Drank(3),6,COMMs(3),req(15),ierr)
!..Receive uc value from SW rank
   call MPI_IRECV(rNE,nzL,MPI_DOUBLE_PRECISION,Drank(1),6,COMMs(3),req(16),ierr)

!.....Bilinear interpolation..................................!
!      do jf=2,nyf-1
!        jc=(jf+1)/2
!        k2=(-1)**IAND(jf,1) ! k2=1 even jf, k2=-1 odd jf
!        do iff=2,nxf-1
!          ic=(iff+1)/2
!          k1=(-1)**IAND(iff,1) ! k=1 even iff, k=-1 odd iff
!          prolng(iff,jf,:)=0.0625_dp*(9.0_dp*uc(ic,jc,:) &
!            +3.0_dp*(uc(ic+k1,jc,:)+uc(ic,jc+k2,:))+uc(ic+k1,jc+k2,:))
!        enddo
!      enddo
      do jf=2,nyf-2,2
        jc=(jf+1)/2
!.......Grid points with iff,jf=even
        do iff=2,nxf-2,2
          ic=(iff+1)/2
          prolng(iff,jf,:)=0.0625_dp*(9.0_dp*uc(ic,jc,:) &
                           +3.0_dp*(uc(ic+1,jc,:)+uc(ic,jc+1,:))+uc(ic+1,jc+1,:))
        enddo
!.......Grid points with iff=odd,jf=even
        do iff=3,nxf-1,2
          ic=(iff+1)/2
          prolng(iff,jf,:)=0.0625_dp*(9.0_dp*uc(ic,jc,:) &
                           +3.0_dp*(uc(ic-1,jc,:)+uc(ic,jc+1,:))+uc(ic-1,jc+1,:))
        enddo
      enddo
      do jf=3,nyf-1,2
        jc=(jf+1)/2
!.......Grid points with iff=even,jf=odd
        do iff=2,nxf-2,2
          ic=(iff+1)/2
          prolng(iff,jf,:)=0.0625_dp*(9.0_dp*uc(ic,jc,:) &
                           +3.0_dp*(uc(ic+1,jc,:)+uc(ic,jc-1,:))+uc(ic+1,jc-1,:))
        enddo
!.......Grid points with iff,jf=odd
        do iff=3,nxf-1,2
          ic=(iff+1)/2
          prolng(iff,jf,:)=0.0625_dp*(9.0_dp*uc(ic,jc,:) &
                           +3.0_dp*(uc(ic-1,jc,:)+uc(ic,jc-1,:))+uc(ic-1,jc-1,:))
        enddo
      enddo

    call MPI_WAIT(req(2),stat(:,2),ierr)
!      do iff=2,nxf-1 ! Bottom row
!        ic=(iff+1)/2
!        k1=(-1)**IAND(iff,1) ! k=1 even iff, k=-1 odd iff
!        prolng(iff,1,:)=0.0625_dp*(9.0_dp*uc(ic,1,:) &
!          +3.0_dp*(uc(ic+k1,1,:)+rbufS(ic,:))+rbufS(ic+k1,:))
!      enddo
      do iff=2,nxf-2,2 ! Grid points with iff=even
        ic=iff/2
        prolng(iff,1,:)=0.0625_dp*(9.0_dp*uc(ic,1,:) &
                        +3.0_dp*(uc(ic+1,1,:)+rbufS(ic,:))+rbufS(ic+1,:))
      enddo
      do iff=3,nxf-1,2 ! Grid points with iff=odd
        ic=(iff+1)/2
        prolng(iff,1,:)=0.0625_dp*(9.0_dp*uc(ic,1,:) &
                        +3.0_dp*(uc(ic-1,1,:)+rbufS(ic,:))+rbufS(ic-1,:))
      enddo
    call MPI_WAIT(req(4),stat(:,4),ierr)
!      do iff=2,nxf-1 ! Top row
!        ic=(iff+1)/2
!        k1=(-1)**IAND(iff,1) ! k=1 even iff, k=-1 odd iff
!        prolng(iff,nyf,:)=0.0625_dp*(9.0_dp*uc(ic,nyc,:) &
!          +3.0_dp*(uc(ic+k1,nyc,:)+rbufN(ic,:))+rbufN(ic+k1,:))
!      enddo
      do iff=2,nxf-2,2 ! Grid points with iff=even
        ic=iff/2
        prolng(iff,nyf,:)=0.0625_dp*(9.0_dp*uc(ic,nyc,:) &
                          +3.0_dp*(uc(ic+1,nyc,:)+rbufN(ic,:))+rbufN(ic+1,:))
      enddo
      do iff=3,nxf-1,2 ! Grid points with iff=odd
        ic=(iff+1)/2
        prolng(iff,nyf,:)=0.0625_dp*(9.0_dp*uc(ic,nyc,:) &
                          +3.0_dp*(uc(ic-1,nyc,:)+rbufN(ic,:))+rbufN(ic-1,:))
      enddo
!.....Left boundary
!      call MPI_WAIT(req(6),stat(:,6),ierr)
      call MPI_WAITALL(3,(/req(6),req(8),req(10)/), &
        (/stat(:,6),stat(:,8),stat(:,10)/),ierr)
    if (iID==0) then
      prolng(1,1,:)=0.0625_dp*(9.0_dp*uc(1,1,:)+3.0_dp*(Lqu*uc(1,1,:) &
        +Lru*uc(2,1,:)+Lsu+rbufS(1,:))+Lqu*rbufS(1,:)+Lru*rbufS(2,:)+Lsu)
!.For some reason this loop causes issues with -O3 on Discovery
!      do jf=2,nyf-1
!        jc=(jf+1)/2
!        k2=(-1)**IAND(jf,1) ! k=1 even jf, k=-1 odd jf
!        prolng(1,jf,:)=0.0625_dp*(9.0_dp*uc(1,jc,:)+3.0_dp*(Lqu*uc(1,jc,:) &
!          +Lru*uc(2,jc,:)+Lsu+uc(1,jc+k2,:))+Lqu*uc(1,jc+k2,:)+Lru*uc(2,jc+k2,:)+Lsu)
!      enddo
      do jf=2,nyf-2,2 ! Grid points with jf=even
        jc=jf/2
        prolng(1,jf,:)=0.0625_dp*(9.0_dp*uc(1,jc,:)+3.0_dp*(Lqu*uc(1,jc,:) &
          +Lru*uc(2,jc,:)+Lsu+uc(1,jc+1,:))+Lqu*uc(1,jc+1,:)+Lru*uc(2,jc+1,:)+Lsu)
      enddo
      do jf=3,nyf-1,2 ! Grid points with jf=odd
        jc=(jf+1)/2
        prolng(1,jf,:)=0.0625_dp*(9.0_dp*uc(1,jc,:)+3.0_dp*(Lqu*uc(1,jc,:) &
          +Lru*uc(2,jc,:)+Lsu+uc(1,jc-1,:))+Lqu*uc(1,jc-1,:)+Lru*uc(2,jc-1,:)+Lsu)
      enddo
      prolng(1,nyf,:)=0.0625_dp*(9.0_dp*uc(1,nyc,:)+3.0_dp*(Lqu*uc(1,nyc,:) &
        +Lru*uc(2,nyc,:)+Lsu+rbufN(1,:))+Lqu*rbufN(1,:)+Lru*rbufN(2,:)+Lsu)
    else
      prolng(1,1,:)=0.0625_dp*(9.0_dp*uc(1,1,:) &
        +3.0_dp*(rbufW(1,:)+rbufS(1,:))+rSW)
!.This loop my cause issues with -O3 on Discovery
!      do jf=2,nyf-1
!        jc=(jf+1)/2
!        k2=(-1)**IAND(jf,1) ! k2=1 even jf, k2=-1 odd jf
!        prolng(1,jf,:)=0.0625_dp*(9.0_dp*uc(1,jc,:) &
!          +3.0_dp*(rbufW(jc,:)+uc(1,jc+k2,:))+rbufW(jc+k2,:))
!      enddo
      do jf=2,nyf-2,2
        jc=(jf+1)/2
!.......Grid points with iff=odd,jf=even
        prolng(1,jf,:)=0.0625_dp*(9.0_dp*uc(1,jc,:) &
                       +3.0_dp*(rbufW(jc,:)+uc(1,jc+1,:))+rbufW(jc+1,:))
      enddo
      do jf=3,nyf-1,2
        jc=(jf+1)/2
!.......Grid points with iff,jf=odd
        prolng(1,jf,:)=0.0625_dp*(9.0_dp*uc(1,jc,:) &
                       +3.0_dp*(rbufW(jc,:)+uc(1,jc-1,:))+rbufW(jc-1,:))
      enddo
      prolng(1,nyf,:)=0.0625_dp*(9.0_dp*uc(1,nyc,:) &
        +3.0_dp*(rbufW(nyc,:)+rbufN(1,:))+rNW)
    endif
      call MPI_WAITALL(3,(/req(12),req(14),req(16)/), &
        (/stat(:,12),stat(:,14),stat(:,16)/),ierr)
    if (iID==iprocsm1) then
!.....Right boundary
      prolng(nxf,1,:)=0.0625_dp*(9.0_dp*uc(nxc,1,:) &
        +3.0_dp*(Rqu*uc(nxc,1,:)+Rru*uc(nxc-1,1,:)+Rsu(:) &
        +rbufS(nxc,:))+Rqu*rbufS(nxc,:)+Rru*rbufS(nxc-1,:)+Rsu(:))
!      do jf=2,nyf-1
!        jc=(jf+1)/2
!        k2=(-1)**IAND(jf,1) ! k=1 even jf, k=-1 odd jf
!        prolng(nxf,jf,:)=0.0625_dp*(9.0_dp*uc(nxc,jc,:) &
!          +3.0_dp*(Rqu*uc(nxc,jc,:)+Rru*uc(nxc-1,jc,:)+Rsu(:) &
!          +uc(nxc,jc+k2,:))+Rqu*uc(nxc,jc+k2,:)+Rru*uc(nxc-1,jc+k2,:)+Rsu(:))
!      enddo
      do jf=2,nyf-2,2 ! Grid points with jf=even
        jc=jf/2
        prolng(nxf,jf,:)=0.0625_dp*(9.0_dp*uc(nxc,jc,:) &
          +3.0_dp*(Rqu*uc(nxc,jc,:)+Rru*uc(nxc-1,jc,:)+Rsu &
          +uc(nxc,jc+1,:))+Rqu*uc(nxc,jc+1,:)+Rru*uc(nxc-1,jc+1,:)+Rsu)
      enddo
      do jf=3,nyf-1,2 ! Grid points with jf=odd
        jc=(jf+1)/2
        prolng(nxf,jf,:)=0.0625_dp*(9.0_dp*uc(nxc,jc,:) &
          +3.0_dp*(Rqu*uc(nxc,jc,:)+Rru*uc(nxc-1,jc,:)+Rsu+uc(nxc,jc-1,:)) &
          +Rqu*uc(nxc,jc-1,:)+Rru*uc(nxc-1,jc-1,:)+Rsu)
      enddo
      prolng(nxf,nyf,:)=0.0625_dp*(9.0_dp*uc(nxc,nyc,:) &
        +3.0_dp*(Rqu*uc(nxc,nyc,:)+Rru*uc(nxc-1,nyc,:)+Rsu(:) &
        +rbufN(nxc,:))+Rqu*rbufN(nxc,:)+Rru*rbufN(nxc-1,:)+Rsu(:))
    else
      prolng(nxf,1,:)=0.0625_dp*(9.0_dp*uc(nxc,1,:) &
        +3.0_dp*(rbufE(1,:)+rbufS(nxc,:))+rSE)
!      do jf=2,nyf-1
!        jc=(jf+1)/2
!        k2=(-1)**IAND(jf,1) ! k2=1 even jf, k2=-1 odd jf
!        prolng(nxf,jf,:)=0.0625_dp*(9.0_dp*uc(nxc,jc,:) &
!          +3.0_dp*(rbufE(jc,:)+uc(nxc,jc+k2,:))+rbufE(jc+k2,:))
!      enddo
      do jf=2,nyf-2,2
        jc=(jf+1)/2
!.......Grid points with iff,jf=even
        prolng(nxf,jf,:)=0.0625_dp*(9.0_dp*uc(nxc,jc,:) &
                         +3.0_dp*(rbufE(jc,:)+uc(nxc,jc+1,:))+rbufE(jc+1,:))
      enddo
      do jf=3,nyf-1,2
        jc=(jf+1)/2
!.......Grid points with iff=even,jf=odd
        prolng(nxf,jf,:)=0.0625_dp*(9.0_dp*uc(nxc,jc,:) &
                         +3.0_dp*(rbufE(jc,:)+uc(nxc,jc-1,:))+rbufE(jc-1,:))
      enddo
      prolng(nxf,nyf,:)=0.0625_dp*(9.0_dp*uc(nxc,nyc,:) &
        +3.0_dp*(rbufE(nyc,:)+rbufN(nxc,:))+rNE)
    endif
      call MPI_WAITALL(8,(/req(1),req(3),req(5),req(7),req(9),req(11), &
        req(13),req(15)/),(/stat(:,1),stat(:,3),stat(:,5),stat(:,7), &
        stat(:,9),stat(:,11),stat(:,13),stat(:,15)/),ierr)

  elseif (ldc==1) then
!.Restrict along x (nyc=nyf)

!...Need to communicate with East and West neighboring procs
    allocate(rbufE(nyc,nzL),sbufE(nyc,nzL),rbufW(nyc,nzL),sbufW(nyc,nzL))
!...Send right column to East rank
    sbufE=uc(nxc,:,:)
    call MPI_ISEND(sbufE,nyc*nzL,MPI_DOUBLE_PRECISION,Erank,1,COMMs(2),req(1),ierr)
!...Receive uc column from West rank, for uc(i-1,j)
    call MPI_IRECV(rbufW,nyc*nzL,MPI_DOUBLE_PRECISION,Wrank,1,COMMs(2),req(2),ierr)
!...Send left column to West rank
    sbufW=uc(1,:,:)
    call MPI_ISEND(sbufW,nyc*nzL,MPI_DOUBLE_PRECISION,Wrank,2,COMMs(2),req(3),ierr)
!...Receive uc column from East rank, for uc(i+1,j)
    call MPI_IRECV(rbufE,nyc*nzL,MPI_DOUBLE_PRECISION,Erank,2,COMMs(2),req(4),ierr)

    if (cof==2) then
      do iff=2,nxf-1
        ic=(iff+1)/2
        p=(-1)**IAND(iff,1) ! p=1 for even iff, p=-1 for odd iff
        prolng(iff,:,:)=0.25_dp*(3.0_dp*uc(ic,:,:)+uc(ic+p,:,:))
      enddo

      call MPI_WAIT(req(2),stat(:,2),ierr)
      if (iID==0) then
        forall(k=1:nzL) prolng(1,:,k)=0.25_dp*(Lqu*uc(1,:,k) &
          +Lru*uc(2,:,k)+Lsu(k))+0.75_dp*uc(1,:,k)
      else
        prolng(1,:,:)=0.25_dp*rbufW+0.75_dp*uc(1,:,:)
      endif
      call MPI_WAIT(req(4),stat(:,4),ierr)
      if (iID==iprocsm1) then
        forall(k=1:nzL) prolng(nxf,:,k)=0.75_dp*uc(nxc,:,k)+0.25_dp*(Rqu*uc(nxc,:,k) &
          +Rru*uc(nxc-1,:,k)+Rsu(k))
      else
        prolng(nxf,:,:)=0.75_dp*uc(nxc,:,:)+0.25_dp*rbufE
      endif

    else
      do iff=4,nxf-4,4
        ic=(iff+1)/4
        prolng(iff-1,:,:)=0.875_dp*uc(ic,:,:)+0.125_dp*uc(ic+1,:,:)
        prolng(iff,:,:)=0.625_dp*uc(ic,:,:)+0.375_dp*uc(ic+1,:,:)
        prolng(iff+1,:,:)=0.375_dp*uc(ic,:,:)+0.625_dp*uc(ic+1,:,:)
        prolng(iff+2,:,:)=0.125_dp*uc(ic,:,:)+0.875_dp*uc(ic+1,:,:)
      enddo

      call MPI_WAIT(req(2),stat(:,2),ierr) 
      if (iID==0) then
        forall(k=1:nzL) prolng(1,:,k)=0.375_dp*(Lqu*uc(1,:,k) &
          +Lru*uc(2,:,k)+Lsu(k))+0.625_dp*uc(1,:,k)
        forall(k=1:nzL) prolng(2,:,k)=0.125_dp*(Lqu*uc(1,:,k) &
          +Lru*uc(2,:,k)+Lsu(k))+0.875_dp*uc(1,:,k)
      else
        prolng(1,:,:)=0.375_dp*rbufW+0.625_dp*uc(1,:,:)
        prolng(2,:,:)=0.125_dp*rbufW+0.875_dp*uc(1,:,:)
      endif
        call MPI_WAIT(req(4),stat(:,4),ierr) 
      if (iID==iprocsm1) then
        forall(k=1:nzL) prolng(nxf-1,:,k)=0.875_dp*uc(nxc,:,k)+0.125_dp*(Rqu*uc(nxc,:,k) &
          +Rru*uc(nxc-1,:,k)+Rsu(k))
        forall(k=1:nzL) prolng(nxf,:,k)=0.625_dp*uc(nxc,:,k)+0.375_dp*(Rqu*uc(nxc,:,k) &
          +Rru*uc(nxc-1,:,k)+Rsu(k))
      else
        prolng(nxf-1,:,:)=0.875_dp*uc(nxc,:,:)+0.125_dp*rbufE
        prolng(nxf,:,:)=0.625_dp*uc(nxc,:,:)+0.375_dp*rbufE
      endif

    endif
    
    call MPI_WAITALL(2,(/req(1),req(3)/),(/stat(:,1),stat(:,3)/),ierr)

  else ! if (ldc==2) then
!.Restrict along y (nxc=nxf)

!...Need to communicate with North and South neighboring procs
    allocate(rbufN(nxc,nzL),sbufN(nxc,nzL),rbufS(nxc,nzL),sbufS(nxc,nzL))
!...Send top row to North rank
    sbufN=uc(:,nyc,:)
    call MPI_ISEND(sbufN,nxc*nzL,MPI_DOUBLE_PRECISION,Nrank,0,COMMs(1),req(1),ierr)
!...Receive uc row from South rank, for uc(i,j-1)
    call MPI_IRECV(rbufS,nxc*nzL,MPI_DOUBLE_PRECISION,Srank,0,COMMs(1),req(2),ierr)
!...Send bottom row to South rank
    sbufS=uc(:,1,:)
    call MPI_ISEND(sbufS,nxc*nzL,MPI_DOUBLE_PRECISION,Srank,1,COMMs(1),req(3),ierr)
!...Receive uc row from North rank, for uc(i,j+1)
    call MPI_IRECV(rbufN,nxc*nzL,MPI_DOUBLE_PRECISION,Nrank,1,COMMs(1),req(4),ierr)

    if (cof==2) then
      do jf=2,nyf-1
        jc=(jf+1)/2
        q=(-1)**IAND(jf,1) ! q=1 for even jf, q=-1 for odd jf
        prolng(:,jf,:)=0.25_dp*(3.0_dp*uc(:,jc,:)+uc(:,jc+q,:))
      enddo

      call MPI_WAIT(req(2),stat(:,2),ierr)
      prolng(:,1,:)=0.25_dp*rbufS+0.75_dp*uc(:,1,:)
      call MPI_WAIT(req(4),stat(:,4),ierr)
      prolng(:,nyf,:)=0.75_dp*uc(:,nyc,:)+0.25_dp*rbufN

    else
      do jf=4,nyf-4,4
        jc=(jf+1)/4
        prolng(:,jf-1,:)=0.875_dp*uc(:,jc,:)+0.125_dp*uc(:,jc+1,:)
        prolng(:,jf,:)=0.625_dp*uc(:,jc,:)+0.375_dp*uc(:,jc+1,:)
        prolng(:,jf+1,:)=0.375_dp*uc(:,jc,:)+0.625_dp*uc(:,jc+1,:)
        prolng(:,jf+2,:)=0.125_dp*uc(:,jc,:)+0.875_dp*uc(:,jc+1,:)
      enddo

      call MPI_WAIT(req(2),stat(:,2),ierr) 
      prolng(:,1,:)=0.375_dp*rbufS+0.625_dp*uc(:,1,:)
      prolng(:,2,:)=0.125_dp*rbufS+0.875_dp*uc(:,1,:)
      call MPI_WAIT(req(4),stat(:,4),ierr) 
      prolng(:,nyf-1,:)=0.875_dp*uc(:,nyc,:)+0.125_dp*rbufN
      prolng(:,nyf,:)=0.625_dp*uc(:,nyc,:)+0.375_dp*rbufN
    endif
    call MPI_WAITALL(2,(/req(1),req(3)/),(/stat(:,1),stat(:,3)/),ierr)
  endif

  end function prolng
!*********************************************************
end module PMGsolver
!*********************************************************
!*********************************************************
