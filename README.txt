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
! 4. A simple Poisson solver
!                 div(grad(u))=rho
!
! In the routines below we use the suffix 'ph' to indicate
! routines to obtain phi, 'ps' for psi and 'hd' for hyperdiffusion.
!
! Though these are 2D solvers, we allow for a 3rd dimension so that
! PMGsolver effetively supports solving many xy planes of
! independent 2D problems.
!
! A program, pmg4gdb.F90, has been included to show the use of some
! of these solvers with a test function.
!
! Manaure Francisquez
! mana@dartmouth.edu
! 2017
!
!********************************************************************
! IMPORTANT NOTES
!   i) Not all solvers support the same boundary conditions
!  ii) Modified Helmholtz relaxation has damped relaxation commented out
! iii) Some restriction and prolongation operators have been removed
!      for simplicity. I believe IR=0, ISR=1, IPRO=2 are the only options
!      available in this file.
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
!   6. Use boundary linear extrapolation when restricting lambda
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
