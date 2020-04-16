# 1 "main.f90"
# 1 "<built-in>"
# 1 "<command-line>"
# 31 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 32 "<command-line>" 2
# 1 "main.f90"
!
! CCCCCCCCCCCCC NNNNNNNN NNNNNNNN SSSSSSSSSSSSSSS
! CCC::::::::::::C N:::::::N N::::::N SS:::::::::::::::S
! CC:::::::::::::::C N::::::::N N::::::N S:::::SSSSSS::::::S
! C:::::CCCCCCCC::::C N:::::::::N N::::::N S:::::S SSSSSSS
! C:::::C CCCCCC aaaaaaaaaaaaa N::::::::::N N::::::N S:::::S
!C:::::C a::::::::::::a N:::::::::::N N::::::N S:::::S
!C:::::C aaaaaaaaa:::::a N:::::::N::::N N::::::N S::::SSSS
!C:::::C a::::a N::::::N N::::N N::::::N SS::::::SSSSS
!C:::::C aaaaaaa:::::a N::::::N N::::N:::::::N SSS::::::::SS
!C:::::C aa::::::::::::a N::::::N N:::::::::::N SSSSSS::::S
!C:::::C a::::aaaa::::::a N::::::N N::::::::::N S:::::S
! C:::::C CCCCCC a::::a a:::::a N::::::N N:::::::::N S:::::S
! C:::::CCCCCCCC::::C a::::a a:::::a N::::::N N::::::::N SSSSSSS S:::::S
! CC:::::::::::::::C a:::::aaaa::::::a N::::::N N:::::::N S::::::SSSSSS:::::S
! CCC::::::::::::C a::::::::::aa:::a N::::::N N::::::N S:::::::::::::::SS
! CCCCCCCCCCCCC aaaaaaaaaa aaaa NNNNNNNN NNNNNNN SSSSSSSSSSSSSSS
!-------------------------------------------------------------------------------------
! CaNS -- Canonical Navier-Stokes Solver
! Pedro Costa (p.simoes.costa@gmail.com)
!-------------------------------------------------------------------------------------
program cans
  use iso_c_binding , only: C_PTR
  use mpi
  use decomp_2d
  use mod_bound , only: boundp,bounduvw,updt_rhs_b
  use mod_chkdiv , only: chkdiv
  use mod_chkdt , only: chkdt
  use mod_common_mpi , only: myid,ierr
  use mod_correc , only: correc
  use mod_debug , only: chkmean
  use mod_fft , only: fftini,fftend
  use mod_fillps , only: fillps
  use mod_initflow , only: initflow
  use mod_initgrid , only: initgrid
  use mod_initmpi , only: initmpi
  use mod_initsolver , only: initsolver
  use mod_load , only: load
  use mod_rk , only: rk,rk_id
  use mod_output , only: out0d,out1d,out1d_2,out2d,out3d
  use mod_param , only: itot,jtot,ktot,lx,ly,lz,dx,dy,dz,dxi,dyi,dzi,uref,lref,rey,visc,small, &
                             cbcvel,bcvel,cbcpre,bcpre, &
                             icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                             nstep,time_max,tw_max,stop_type,restart, &
                             rkcoeff, &
                             datadir, &
                             cfl,dtmin, &
                             inivel, &
                             imax,jmax,dims, &
                             nthreadsmax, &
                             gr, &
                             is_outflow,no_outflow,is_forced,bforce, &
                             n,ng,l,dl,dli, &
                             read_input
  use mod_sanity , only: test_sanity
  use mod_solver , only: solver
  use mod_types
  !$ use omp_lib
  implicit none
  real(rp), allocatable, dimension(:,:,:) :: u,v,w,p,up,vp,wp,pp
  real(rp), allocatable, dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
  real(rp), dimension(3) :: tauxo,tauyo,tauzo
  real(rp), dimension(3) :: f
  type(C_PTR), dimension(2,2) :: arrplanp
  real(rp), allocatable, dimension(:,:) :: lambdaxyp
  real(rp), allocatable, dimension(:) :: ap,bp,cp
  real(rp) :: normfftp
  type rhs_bound
    real(rp), allocatable, dimension(:,:,:) :: x
    real(rp), allocatable, dimension(:,:,:) :: y
    real(rp), allocatable, dimension(:,:,:) :: z
  end type rhs_bound
  type(rhs_bound) :: rhsbp
# 83 "main.f90"
  real(rp) :: ristep
  real(rp) :: dt,dti,dtmax,time,dtrk,dtrki,divtot,divmax
  integer :: irk,istep
  real(rp), allocatable, dimension(:) :: dzc,dzf,zc,zf,dzci,dzfi
  real(rp) :: meanvelu,meanvelv,meanvelw
  real(rp), dimension(3) :: dpdl
  !real(rp), allocatable, dimension(:) :: var
  real(rp), dimension(10) :: var



  real(rp) :: twi,tw
  character(len=7) :: fldnum
  integer :: kk
  logical :: is_done,kill
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  !
  ! read parameter file
  !
  call read_input(myid)
  !
  ! initialize MPI/OpenMP
  !
  !$call omp_set_num_threads(nthreadsmax)
  call initmpi(ng,cbcpre)
  twi = MPI_WTIME()
  !
  ! allocate variables
  !
  allocate(u( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           v( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           w( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           p( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           up(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           vp(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           wp(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
           pp(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(dudtrko(n(1),n(2),n(3)), &
           dvdtrko(n(1),n(2),n(3)), &
           dwdtrko(n(1),n(2),n(3)))
  allocate(lambdaxyp(n(1),n(2)))
  allocate(ap(n(3)),bp(n(3)),cp(n(3)))
  allocate(dzc( 0:n(3)+1), &
           dzf( 0:n(3)+1), &
           zc( 0:n(3)+1), &
           zf( 0:n(3)+1), &
           dzci(0:n(3)+1), &
           dzfi(0:n(3)+1))
  allocate(rhsbp%x(n(2),n(3),0:1), &
           rhsbp%y(n(1),n(3),0:1), &
           rhsbp%z(n(1),n(2),0:1))
# 154 "main.f90"
  if(myid.eq.0) print*, '*******************************'
  if(myid.eq.0) print*, '*** Beginning of simulation ***'
  if(myid.eq.0) print*, '*******************************'
  if(myid.eq.0) print*, ''
  call initgrid(inivel,n(3),gr,lz,dzc,dzf,zc,zf)
  dzci(:) = dzc(:)**(-1)
  dzfi(:) = dzf(:)**(-1)
  if(myid.eq.0) then
    open(99,file=trim(datadir)
    write(99,rec=1) dzc(1:n(3)),dzf(1:n(3)),zc(1:n(3)),zf(1:n(3))
    close(99)
    open(99,file=trim(datadir)
    do kk=0,ktot+1
      write(99,'(5E15.7)') 0.,zf(kk),zc(kk),dzf(kk),dzc(kk)
    enddo
    close(99)
  endif
  !
  ! test input files before proceeding with the calculation
  !
  call test_sanity(ng,n,dims,stop_type,cbcvel,cbcpre,bcvel,bcpre,is_outflow,is_forced, &
                   dli,dzci,dzfi)
  !
  if(.not.restart) then
    istep = 0
    time = 0.
    call initflow(inivel,n,zc/lz,dzc/lz,dzf/lz,visc,u,v,w,p)
    if(myid.eq.0) print*, '*** Initial condition succesfully set ***'
  else
    call load('r',trim(datadir)
                                             v(1:n(1),1:n(2),1:n(3)), &
                                             w(1:n(1),1:n(2),1:n(3)), &
                                             p(1:n(1),1:n(2),1:n(3)), &
                                             time,ristep)
    istep = nint(ristep)
    if(myid.eq.0) print*, '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
  endif
  call bounduvw(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,u,v,w)
  call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
  !
  ! post-process and write initial condition
  !
  write(fldnum,'(i7.7)') istep
  include 'out1d.h90'
  include 'out2d.h90'
  include 'out3d.h90'
  !
  dudtrko(:,:,:) = 0.
  dvdtrko(:,:,:) = 0.
  dwdtrko(:,:,:) = 0.
  call chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
  dt = min(cfl*dtmax,dtmin)
  if(myid.eq.0) print*, 'dtmax = ', dtmax, 'dt = ',dt
  dti = 1./dt
  kill = .false.
  !
  ! initialize Poisson solver
  !
  call initsolver(n,dli,dzci,dzfi,cbcpre,bcpre(:,:),lambdaxyp,(/'c','c','c'/),ap,bp,cp,arrplanp,normfftp,rhsbp%x,rhsbp%y,rhsbp%z)
# 221 "main.f90"
  !
  ! main loop
  !
  if(myid.eq.0) print*, '*** Calculation loop starts now ***'
  is_done = .false.
  do while(.not.is_done)



    istep = istep + 1
    time = time + dt
    if(myid.eq.0) print*, 'Timestep #', istep, 'Time = ', time
    dpdl(:) = 0.
    tauxo(:) = 0.
    tauyo(:) = 0.
    tauzo(:) = 0.
    do irk=1,3
      dtrk = sum(rkcoeff(:,irk))*dt
      dtrki = dtrk**(-1)

      call rk(rkcoeff(:,irk),n,dli,dzci,dzfi,dzf/lz,dzc/lz,visc,dt,l, &
                 u,v,w,p,dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)




      if(is_forced(1)) up(1:n(1),1:n(2),1:n(3)) = up(1:n(1),1:n(2),1:n(3)) + f(1)
      if(is_forced(2)) vp(1:n(1),1:n(2),1:n(3)) = vp(1:n(1),1:n(2),1:n(3)) + f(2)
      if(is_forced(3)) wp(1:n(1),1:n(2),1:n(3)) = wp(1:n(1),1:n(2),1:n(3)) + f(3)
# 271 "main.f90"
      dpdl(:) = dpdl(:) + f(:)
      call bounduvw(cbcvel,n,bcvel,no_outflow,dl,dzc,dzf,up,vp,wp) ! outflow BC only at final velocity
# 287 "main.f90"
      call fillps(n,dli,dzfi,dtrki,up,vp,wp,pp)
      call updt_rhs_b((/'c','c','c'/),cbcpre,n,rhsbp%x,rhsbp%y,rhsbp%z,pp)
      call solver(n,arrplanp,normfftp,lambdaxyp,ap,bp,cp,cbcpre(:,3),(/'c','c','c'/),pp)
      call boundp(cbcpre,n,bcpre,dl,dzc,dzf,pp)
      call correc(n,dli,dzci,dtrk,pp,up,vp,wp,u,v,w)
      call bounduvw(cbcvel,n,bcvel,is_outflow,dl,dzc,dzf,u,v,w)
# 317 "main.f90"
      !$OMP WORKSHARE
      p(1:n(1),1:n(2),1:n(3)) = p(1:n(1),1:n(2),1:n(3)) + pp(1:n(1),1:n(2),1:n(3))
      !$OMP END WORKSHARE

      call boundp(cbcpre,n,bcpre,dl,dzc,dzf,p)
    enddo
    dpdl(:) = -dpdl(:)*dti
    !
    ! check simulation stopping criteria
    !
    if(stop_type(1)) then ! maximum number of time steps reached
      if(istep.ge.nstep ) is_done = is_done.or..true.
    endif
    if(stop_type(2)) then ! maximum simulation time reached
      if(time .ge.time_max) is_done = is_done.or..true.
    endif
    if(stop_type(3)) then ! maximum wall-clock time reached
      tw = (MPI_WTIME()-twi)/3600.
      if(tw .ge.tw_max ) is_done = is_done.or..true.
    endif
    if(mod(istep,icheck).eq.0) then
      if(myid.eq.0) print*, 'Checking stability and divergence...'
      call chkdt(n,dl,dzci,dzfi,visc,u,v,w,dtmax)
      dt = min(cfl*dtmax,dtmin)
      if(myid.eq.0) print*, 'dtmax = ', dtmax, 'dt = ',dt
      if(dtmax.lt.small) then
        if(myid.eq.0) print*, 'ERROR: timestep is too small.'
        if(myid.eq.0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      endif
      dti = 1./dt
      call chkdiv(n,dli,dzfi,u,v,w,divtot,divmax)
      if(divmax.gt.small.or.divtot.ne.divtot) then
        if(myid.eq.0) print*, 'ERROR: maximum divergence is too large.'
        if(myid.eq.0) print*, 'Aborting...'
        is_done = .true.
        kill = .true.
      endif
    endif
    !
    ! output routines below
    !
    if(mod(istep,iout0d).eq.0) then
      !allocate(var(4))
      var(1) = 1.*istep
      var(2) = dt
      var(3) = time
      call out0d(trim(datadir)
      !
      if(any(is_forced(:)).or.any(abs(bforce(:)).gt.0.)) then
        meanvelu = 0.
        meanvelv = 0.
        meanvelw = 0.
        if(is_forced(1).or.abs(bforce(1)).gt.0.) then
          call chkmean(n,dzf/lz,up,meanvelu)
        endif
        if(is_forced(2).or.abs(bforce(2)).gt.0.) then
          call chkmean(n,dzf/lz,vp,meanvelv)
        endif
        if(is_forced(3).or.abs(bforce(3)).gt.0.) then
          call chkmean(n,dzc/lz,wp,meanvelw)
        endif
        if(.not.any(is_forced(:))) dpdl(:) = -bforce(:) ! constant pressure gradient
        var(1) = time
        var(2:4) = dpdl(1:3)
        var(5:7) = (/meanvelu,meanvelv,meanvelw/)
        call out0d(trim(datadir)
      endif
      !deallocate(var)
    endif
    write(fldnum,'(i7.7)') istep
    if(mod(istep,iout1d).eq.0) then
      include 'out1d.h90'
    endif
    if(mod(istep,iout2d).eq.0) then
      include 'out2d.h90'
    endif
    if(mod(istep,iout3d).eq.0) then
      include 'out3d.h90'
    endif
    if(mod(istep,isave ).eq.0.or.(is_done.and..not.kill)) then
      ristep = 1.*istep
      call load('w',trim(datadir)
                                               v(1:n(1),1:n(2),1:n(3)), &
                                               w(1:n(1),1:n(2),1:n(3)), &
                                               p(1:n(1),1:n(2),1:n(3)), &
                                               time,ristep)
      if(myid.eq.0) print*, '*** Checkpoint saved at time = ', time, 'time step = ', istep, '. ***'
    endif
# 415 "main.f90"
  enddo
  !
  ! clear ffts
  !
  call fftend(arrplanp)





  !
  ! deallocate variables
  !
  deallocate(u,v,w,p,up,vp,wp,pp)
  deallocate(dudtrko,dvdtrko,dwdtrko)
  deallocate(lambdaxyp)
  deallocate(ap,bp,cp)
  deallocate(dzc,dzf,zc,zf,dzci,dzfi)
  deallocate(rhsbp%x,rhsbp%y,rhsbp%z)







  if(myid.eq.0.and.(.not.kill)) print*, '*** Fim ***'
  call decomp_2d_finalize
  call MPI_FINALIZE(ierr)
  call exit
end program cans
