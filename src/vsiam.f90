Module MOD_MOM_VSIAM
  Use mpi
  use mod_param     , only: dims, bforce
  use mod_common_mpi, only: ierr
  Implicit None
  private
  save
  public

  Real(rp) :: tvb_m = 0._rp    ! totoal varation bounded for CIP-CSL3
  Real(rp) :: tec_m = 1.0_rp  ! modification for time evolution compute
  Real(rp) :: grad_m = 1._rp  ! modification for time evolution compute
  Real(rp), dimension(:,:,:), allocatable :: ux, uy, yz, vx, vy, vz, wx, wy, wz
  Real(rp), dimension(:), allocatable :: dx, dy, dz, dxc, dyc, dzc
Contains

  Subroutine Fake_VSIAM_init(n,dli)
    Implicit None
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in   ), dimension(3) :: dli
    ! CIP variables
    allocate(ucip( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
        vcip(0:n(1)+1,0:n(2)+1,0:n(3)+1), &
        wcip(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    allocate(ux( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
        uy( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
        uz( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
        vx( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
        vy( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
        vz( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
        wx( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
        wy( 0:n(1)+1,0:n(2)+1,0:n(3)+1), &
        wz(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    allocate(dx(0:nx+1))
    allocate(dy(0:ny+1))
    allocate(dz(0:nz+1))
    allocate(dxc(0:nx+1))
    allocate(dyc(0:ny+1))
    allocate(dzc(0:nz+1))
    dx  = dli(1)
    dxc = dli(1)
    dy  = dli(2)
    dyc = dli(2)
    dz  = dli(3)
    dzc = dli(3)

  end Subroutine Fake_VSIAM_init

  subroutine Fake_VSIAM(rkpar,n,dli,dzci,dzfi,dzflzi,dzclzi,visc,dt,l,u,v,w,p,dudtrko,dvdtrko,dwdtrko,tauxo,tauyo,tauzo,up,vp,wp,f)
    !
    ! A very simplified VSIAM...
    ! Cell-centered variables are obtained from the interpolation of the Face-centered variables
    !
    implicit none
    real(rp), intent(in), dimension(2) :: rkpar
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in) :: visc,dt
    real(rp), intent(in   ), dimension(3) :: dli,l 
    real(rp), intent(in   ), dimension(0:) :: dzci,dzfi,dzflzi,dzclzi
    real(rp), intent(in   ), dimension(0:,0:,0:) :: u ,v ,w,p
    real(rp), intent(inout), dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
    real(rp), intent(inout), dimension(3) :: tauxo,tauyo,tauzo
    real(rp), intent(out), dimension(0:,0:,0:) :: up,vp,wp
    real(rp), intent(out), dimension(3) :: f
    real(rp),              dimension(n(1),n(2),n(3)) :: dudtrk,dvdtrk,dwdtrk
    real(rp),              dimension(0:n(1),0:n(2),0:n(3)) :: dudtrk,dvdtrk,dwdtrk
    real(rp), dimension(0:,0:,0:) :: uc,vc,wc
    real(rp) :: factor1,factor2,factor12
    real(rp), dimension(3) :: taux,tauy,tauz
    integer :: i,j,k
    real(rp) :: mean

    uc = 0.0_rp
    vc = 0.0_rp
    wc = 0.0_rp


    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    !

    call momxp(n(1),n(2),n(3),dli(1),p,dudtrk)
    call momyp(n(1),n(2),n(3),dli(2),p,dvdtrk)
    call momzp(n(1),n(2),n(3),dzci  ,p,dwdtrk)


    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor12,u,v,w,up,vp,wp,dudtrk,dvdtrk,dwdtrk)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          up(i,j,k) = u(i,j,k) + factor12*dudtrk(i,j,k)
          vp(i,j,k) = v(i,j,k) + factor12*dvdtrk(i,j,k)
          wp(i,j,k) = w(i,j,k) + factor12*dwdtrk(i,j,k)
        enddo
      enddo
    enddo

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor12,u,v,w,up,vp,wp,dudtrk,dvdtrk,dwdtrk)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          UN(i,j,j) = UP(i,j,k)
          VN(i,j,j) = WP(i,j,k)
          WN(i,j,j) = WP(i,j,k)
          ucip(i,j,j) = UP(i,j,k)
          vcip(i,j,j) = VP(i,j,k)
          wcip(i,j,j) = WP(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO

    itp = 1
    Call CIP3D(itp,UCIP,UX,UY,UZ,VCIP,VX,VY,VZ,WCIP,WX,WY,WZ &
        &		, DX,DXC,DY,DYC,DZ,DZC , nx,ny,nz,dt )
    itp = 2
    Call CIP3D(itp,UCIP,UX,UY,UZ,VCIP,VX,VY,VZ,WCIP,WX,WY,WZ &
        &		, DX,DXC,DY,DYC,DZ,DZC , nx,ny,nz,dt )
    itp = 3
    Call CIP3D(itp,UCIP,UX,UY,UZ,VCIP,VX,VY,VZ,WCIP,WX,WY,WZ &
        &		, DX,DXC,DY,DYC,DZ,DZC , nx,ny,nz,dt )

    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          dudtrk(i,j,j) = ucip(i,j,k) - wn(i,j,k)
          dvdtrk(i,j,j) = vcip(i,j,k) - vn(i,j,k)
          dwdtrk(i,j,j) = wcip(i,j,k) - wn(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call momxad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,dzflzi,visc,u,v,w,dudtrk,taux)
    call momyad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,dzflzi,visc,u,v,w,dvdtrk,tauy)
    call momzad(n(1),n(2),n(3),dli(1),dli(2),dli(3),dzci,dzfi,dzflzi,visc,u,v,w,dwdtrk,tauz)

    f(1) = (factor1*sum(taux(:)/l(:)) + factor2*sum(tauxo(:)/l(:)))
    f(2) = (factor1*sum(tauy(:)/l(:)) + factor2*sum(tauyo(:)/l(:)))
    f(3) = (factor1*sum(tauz(:)/l(:)) + factor2*sum(tauzo(:)/l(:)))
    tauxo(:) = taux(:)
    tauyo(:) = tauy(:)
    tauzo(:) = tauz(:)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(n,factor1,factor2,up,vp,wp,dudtrk,dvdtrk,dwdtrk,dudtrko,dvdtrko,dwdtrko)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          ! could be split in two loops, because factor2=0 for istep=1, but like this reads nicer
          up(i,j,k) = up(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k)
          vp(i,j,k) = vp(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k)
          wp(i,j,k) = wp(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k)
          dudtrko(i,j,k) = dudtrk(i,j,k)
          dvdtrko(i,j,k) = dvdtrk(i,j,k)
          dwdtrko(i,j,k) = dwdtrk(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !
    ! bulk velocity forcing
    !
    f(:) = 0.
    if(is_forced(1)) then
      call chkmean(n,dzflzi,up,mean)
      f(1) = velf(1) - mean
    endif
    if(is_forced(2)) then
      call chkmean(n,dzflzi,vp,mean)
      f(2) = velf(2) - mean
    endif
    if(is_forced(3)) then
      call chkmean(n,dzclzi,wp,mean)
      f(3) = velf(3) - mean
    endif
    return
  end subroutine Fake_VSIAM

subroutine CIP3D(itp,U,UX,UY,UZ,V,VX,VY,VZ,W,WX,WY,WZ &
			&		, DX,DXC,DY,DYC,DZ,DZC , nx,ny,nz,dt )

	implicit none

	real(8) , intent(in)	::	itp
    
	real(8) , intent(in)	::	dt
	integer , intent(in)	::	nx
	integer , intent(in)	::	ny
	integer , intent(in)	::	nz
    

	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1) , intent(inout)	::	U
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1) , intent(inout)	::	UX
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1) , intent(inout)	::	UY
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1) , intent(inout)	::	UZ
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1) , intent(inout)	::	V
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1) , intent(inout)	::	VX
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1) , intent(inout)	::	VY
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1) , intent(inout)	::	VZ
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1) , intent(inout)	::	W
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1) , intent(inout)	::	WX
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1) , intent(inout)	::	WY
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1) , intent(inout)	::	WZ

	real(8) , dimension(0:nx+1) , intent(in)	::	DX
	real(8) , dimension(0:nx+1) , intent(in)	::	DXC
	real(8) , dimension(0:ny+1) , intent(in)	::	DY
	real(8) , dimension(0:ny+1) , intent(in)	::	DYC
	real(8) , dimension(0:nz+1) , intent(in)	::	DZ
	real(8) , dimension(0:nz+1) , intent(in)	::	DZC
    
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1)	::	F   
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1)	::	GX
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1)	::	GY
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1)	::	GZ

	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1)	::	FA
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1)	::	UC
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1)	::	VC
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1)	::	WC
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1)	::	GXA
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1)	::	GYA
	real(8) , dimension(0:nx+1,0:ny+1,0:nz+1)	::	GZA

	real(8) , dimension(0:nx+1)	::	DXA
	real(8) , dimension(0:ny+1)	::	DYA
	real(8) , dimension(0:nz+1)	::	DZA

	real(8)	::	csmt,csmt1
	real(8)	::	xx,yy,zz
	real(8)	::	dx1,dy1,dz1 , dx2,dy2,dz2 , dx3,dy3,dz3
	real(8)	::	a300,a030,a003,a210,a120,a201,a102,a021,a012,a111
	real(8)	::	a200,a020,a002,a110,a101,a011,a100,a010,a001,a000
	real(8)	::	c00,c10,c01,c11
	real(8)	::	tmp1,tmp2,tmp3
	real(8)	::	tmp210,tmp021,tmp201
	real(8)	::	tmp120,tmp012,tmp102
	real(8)	::	chck
	real(8)	::	phi20
	real(8)	::	fa0,gxa0,gya0
    
	integer	::	im1,jm1,km1 , isn,jsn,ksn , ima1,jma1,kma1 
	integer	::	i,j,k

	if (int(itp).EQ.1) then
        F=U
        GX=UX
        GY=UY
        GZ=UZ
    elseif (int(itp).EQ.2) then
        F=V
        GX=VX
        GY=VY
        GZ=VZ
    elseif (int(itp).EQ.3) then
    	F=W
    	GX=WX
    	GY=WY
    	GZ=WZ
    end if

	csmt  = 1.0

	DXA(0:nx+1) =  DX(0:nx+1)
	DYA(0:ny+1) =  DY(0:ny+1)
	DZA(0:nz+1) =  DZ(0:nz+1)

	if ( int(itp).EQ.1 ) then
 
		UC(1:nx,1:ny,1:nz) = U(1:nx,1:ny,1:nz)
		do j=1,ny
		do k=1,nz
			VC(1:nx,j,k)=0.25*( ( V(1:nx,j,k)   + V(1:nx,j-1,k) )*DX(2:nx+1) +  &
				&			  ( V(2:nx+1,j,k) + V(2:nx+1,j-1,k) )*DX(1:nx)   )  &
                &           / DXC(1:nx)
            WC(1:nx,j,k)=0.25*( ( W(1:nx,j,k)   + W(1:nx,j,k-1) )*DX(2:nx+1) +  &
				&			  ( W(2:nx+1,j,k) + W(2:nx+1,j,k-1) )*DX(1:nx)   )  &
                &           / DXC(1:nx)
		end do
		end do
 
		DXA(1:nx) 	=  DX(1:nx)
		DYA(1:ny)	=  DYC(0:ny-1)
		DZA(1:nz)	=  DZC(0:nz-1)

		UC(    0, 1:ny, 1:nz )	=  0.d0
		VC(    0, 1:ny, 1:nz )	=  0.d0
		WC(    0, 1:ny, 1:nz )	=  0.d0
		UC(   nx, 1:ny, 1:nz )	=  0.d0                                 
		VC(   nx, 1:ny, 1:nz )	=  0.d0
		WC(   nx, 1:ny, 1:nz )	=  0.d0
 
		UC( 1:nx,    0, 1:nz )	=  UC( 1:nx,    1, 1:nz )
		VC( 1:nx,    0, 1:nz )	= -VC( 1:nx,    1, 1:nz )
		WC( 1:nx,    0, 1:nz )	=  WC( 1:nx,    1, 1:nz )
		UC( 1:nx, ny+1, 1:nz )	= -UC( 1:nx,   ny, 1:nz )
		VC( 1:nx, ny+1, 1:nz )	=  VC( 1:nx,   ny, 1:nz )
		WC( 1:nx, ny+1, 1:nz )	= -WC( 1:nx,   ny, 1:nz )

		UC( 1:nx, 1:ny,    0 )	= -UC( 1:nx, 1:ny,    1 )
		VC( 1:nx, 1:ny,    0 )	= -VC( 1:nx, 1:ny,    1 )
		WC( 1:nx, 1:ny,    0 )	=  WC( 1:nx, 1:ny,    1 )
		UC( 1:nx, 1:ny, nz+1 )	= -UC( 1:nx, 1:ny,   nz )
		VC( 1:nx, 1:ny, nz+1 )	= -VC( 1:nx, 1:ny,   nz )
		WC( 1:nx, 1:ny, nz+1 )	=  WC( 1:nx, 1:ny,   nz )

	elseif ( int(itp).EQ.2 ) then

		VC(1:nx,1:ny,1:nz) = V(1:nx,1:ny,1:nz)
		do k = 1,nz
		do i = 1,nx
			UC(i,1:ny,k)=0.25*( (U(i,1:ny,k)   + U(i-1,1:ny,k)  )*DY(2:ny+1) + &
				&			  (U(i,2:ny+1,k) + U(i-1,2:ny+1,k))*DY(1:ny)  )  &
				&		/DYC(1:ny)
			WC(i,1:ny,k)=0.25*( (W(i,1:ny,k)   + W(i,1:ny,k-1)  )*DY(2:ny+1) + &
				&			  (W(i,2:ny+1,k) + W(i,2:ny+1,k-1))*DY(1:ny)  )  &
				&		/DYC(1:ny)
		end do
		end do
 
		DXA(1:nx) 	=  DXC(0:nx-1)
		DYA(1:ny)	=  DY(1:ny)
		DZA(1:nz)	=  DZC(0:nz-1)

		UC(    0, 1:ny, 1:nz )	= -UC(    1, 1:ny, 1:nz )
		VC(    0, 1:ny, 1:nz )	=  VC(    1, 1:ny, 1:nz )
		WC(    0, 1:ny, 1:nz )	=  WC(    1, 1:ny, 1:nz )
		UC( nx+1, 1:ny, 1:nz )	=  UC(   nx, 1:ny, 1:nz )                              
		VC( nx+1, 1:ny, 1:nz )	= -VC(   nx, 1:ny, 1:nz )
		WC( nx+1, 1:ny, 1:nz )	= -WC(   nx, 1:ny, 1:nz )

		UC( 1:nx,    0, 1:nz )	=  0.d0
		VC( 1:nx,    0, 1:nz )	=  0.d0
		WC( 1:nx,    0, 1:nz )	=  0.d0
		UC( 1:nx,   ny, 1:nz )	=  0.d0
		VC( 1:nx,   ny, 1:nz )	=  0.d0
		WC( 1:nx,   ny, 1:nz )	=  0.d0

		UC( 1:nx, 1:ny,    0 )	=  UC( 1:nx, 1:ny,    1 )
		VC( 1:nx, 1:ny,    0 )	=  VC( 1:nx, 1:ny,    1 )
		WC( 1:nx, 1:ny,    0 )	= -WC( 1:nx, 1:ny,    1 )
		UC( 1:nx, 1:ny, nz+1 )	= -UC( 1:nx, 1:ny,   nz )
		VC( 1:nx, 1:ny, nz+1 )	= -VC( 1:nx, 1:ny,   nz )
		WC( 1:nx, 1:ny, nz+1 )	=  WC( 1:nx, 1:ny,   nz )

	elseif ( int(itp).EQ.3 ) then
 
		WC(1:nx,1:ny,1:nz) = W(1:nx,1:ny,1:nz)
		do j = 1,ny
		do i = 1,nx
			UC(i,j,1:nz)=0.25*( (U(i,j,1:nz)   + U(i-1,j,1:nz)  )*DZ(2:nz+1) + &
				&			  (U(i,j,2:nz+1) + U(i-1,j,2:nz+1))*DZ(1:nz)  )  &
				&		/DZC(1:nz)
			VC(i,j,1:nz)=0.25*( (V(i,j,1:nz)   + V(i,j-1,1:nz)  )*DZ(2:nz+1) + &
				&			  (V(i,j,2:nz+1) + V(i,j-1,2:nz+1))*DZ(1:nz)  )  &
				&		/DZC(1:nz)
		end do
		end do
 
		DXA(1:nx) 	=  DXC(0:nx-1)
		DYA(1:ny)	=  DYC(0:ny-1)
		DZA(1:nz)	=  DZ(1:nz)

		UC(    0, 1:ny, 1:nz )	= -UC(    1, 1:ny, 1:nz )
 		VC(    0, 1:ny, 1:nz )	=  VC(    1, 1:ny, 1:nz )
 		WC(    0, 1:ny, 1:nz )	=  WC(    1, 1:ny, 1:nz )
 		UC( nx+1, 1:ny, 1:nz )	=  UC(   nx, 1:ny, 1:nz )                              
		VC( nx+1, 1:ny, 1:nz )	= -VC(   nx, 1:ny, 1:nz )
		WC( nx+1, 1:ny, 1:nz )	= -WC(   nx, 1:ny, 1:nz )
 
		UC( 1:nx,    0, 1:nz )	=  UC( 1:nx,    1, 1:nz )
		VC( 1:nx,    0, 1:nz )	= -VC( 1:nx,    1, 1:nz )
		WC( 1:nx,    0, 1:nz )	=  WC( 1:nx,    1, 1:nz )
		UC( 1:nx, ny+1, 1:nz )	= -UC( 1:nx,   ny, 1:nz )
		VC( 1:nx, ny+1, 1:nz )	=  VC( 1:nx,   ny, 1:nz )
		WC( 1:nx, ny+1, 1:nz )	= -WC( 1:nx,   ny, 1:nz )

		UC( 1:nx, 1:ny,    0 )	=  0.d0
		VC( 1:nx, 1:ny,    0 )	=  0.d0
		WC( 1:nx, 1:ny,    0 )	=  0.d0
		UC( 1:nx, 1:ny,   nz )	=  0.d0
		VC( 1:nx, 1:ny,   nz )	=  0.d0
		WC( 1:nx, 1:ny,   nz )	=  0.d0
 
	else
 
		UC(1:nx,1:ny,1:nz) = 0.5 * ( U(1:nx,1:ny,1:nz)+U(0:nx-1,1:ny,1:nz) )
        VC(1:nx,1:ny,1:nz) = 0.5 * ( V(1:nx,1:ny,1:nz)+V(1:nx,0:ny-1,1:nz) )
        WC(1:nx,1:ny,1:nz) = 0.5 * ( W(1:nx,1:ny,1:nz)+W(1:nx,1:ny,0:nz-1) )

        DXA(1:nx)=DX(0:nx-1)
		DYA(1:ny)=DY(0:ny-1)
		DZA(1:nz)=DZ(0:nz-1)
 
		UC(    0, 1:ny, 1:nz )	=  UC(    1, 1:ny, 1:nz )
    	VC(    0, 1:ny, 1:nz )	= -VC(    1, 1:ny, 1:nz )
    	WC(    0, 1:ny, 1:nz )	= -WC(    1, 1:ny, 1:nz )
    	UC( nx+1, 1:ny, 1:nz )	=  UC(   nx, 1:ny, 1:nz )                              
    	VC( nx+1, 1:ny, 1:nz )	= -VC(   nx, 1:ny, 1:nz )
    	WC( nx+1, 1:ny, 1:nz )	= -WC(   nx, 1:ny, 1:nz )
 
    	UC( 1:nx,    0, 1:nz )	= -UC( 1:nx,    1, 1:nz )
    	VC( 1:nx,    0, 1:nz )	=  VC( 1:nx,    1, 1:nz )
    	WC( 1:nx,    0, 1:nz )	= -WC( 1:nx,    1, 1:nz )
    	UC( 1:nx, ny+1, 1:nz )	= -UC( 1:nx,   ny, 1:nz )
    	VC( 1:nx, ny+1, 1:nz )	=  VC( 1:nx,   ny, 1:nz )
    	WC( 1:nx, ny+1, 1:nz )	= -WC( 1:nx,   ny, 1:nz )

    	UC( 1:nx, 1:ny,    0 )	= -UC( 1:nx, 1:ny,    1 )
    	VC( 1:nx, 1:ny,    0 )	= -VC( 1:nx, 1:ny,    1 )
    	WC( 1:nx, 1:ny,    0 )	=  WC( 1:nx, 1:ny,    1 )
    	UC( 1:nx, 1:ny, nz+1 )	= -UC( 1:nx, 1:ny,   nz )
    	VC( 1:nx, 1:ny, nz+1 )	= -VC( 1:nx, 1:ny,   nz )
    	WC( 1:nx, 1:ny, nz+1 )	=  WC( 1:nx, 1:ny,   nz )

 	end if
 
 	DXA(0)=DXA(1)
 	DYA(0)=DYA(1)
 	DZA(0)=DZA(1)
	DXA(nx+1)=DXA(nx)
 	DYA(ny+1)=DYA(ny)
 	DZA(nz+1)=DZA(nz)
 
 	do k = 1,nz
 	do j = 1,ny
 	do i = 1,nx
 
 		IF( (int(itp).EQ.1).AND.(i.EQ.nx) ) CYCLE
    	IF( (int(itp).EQ.2).AND.(j.EQ.ny) ) CYCLE
    	IF( (int(itp).EQ.3).AND.(k.EQ.nz) ) CYCLE
 
    	xx   = -UC(i,j,k)*dt
		yy   = -VC(i,j,k)*dt
		zz   = -WC(i,j,k)*dt
		isn  =  dsign(1.d0,UC(i,j,k))
		jsn  =  dsign(1.d0,VC(i,j,k))
		ksn  =  dsign(1.d0,WC(i,j,k))
        if (UC(i,j,k).EQ.0) isn=0
        if (VC(i,j,k).EQ.0) jsn=0
        if (WC(i,j,k).EQ.0) ksn=0
		im1  =  i-isn
		jm1  =  j-jsn
		km1  =  k-ksn

		ima1 =  MAX(im1,I)
		jma1 =  MAX(jm1,J)
		kma1 =  MAX(km1,k)
        
		dx1  =  1.d0/DXA(ima1)
		dy1  =  1.d0/DYA(jma1)
		dz1  =  1.d0/DZA(kma1)
		dx2  =  dx1*dx1
		dx3  =  dx2*dx1
		dy2  =  dy1*dy1
		dy3  =  dy2*dy1
		dz2  =  dz1*dz1
		dz3  =  dz2*dz1

		a300 =  ((GX(im1,j,k)+GX(i,j,k))*DXA(ima1)*isn-2.0*(F(i,j,k)-F(im1,j,k)))*dx3*isn  
		a030 =  ((GY(i,jm1,k)+GY(i,j,k))*DYA(jma1)*jsn-2.0*(F(i,j,k)-F(i,jm1,k)))*dy3*jsn
		a003 =  ((GZ(i,j,km1)+GZ(i,j,k))*DZA(kma1)*ksn-2.0*(F(i,j,k)-F(i,j,km1)))*dz3*ksn

		a200 =  (3.0*(F(im1,j,k)-F(i,j,k))+(GX(im1,j,k)+2.0*GX(i,j,k))*DXA(ima1)*isn)*dx2
		a020 =  (3.0*(F(i,jm1,k)-F(i,j,k))+(GY(i,jm1,k)+2.0*GY(i,j,k))*DYA(jma1)*jsn)*dy2
		a002 =  (3.0*(F(i,j,km1)-F(i,j,k))+(GZ(i,j,km1)+2.0*GZ(i,j,k))*DZA(kma1)*ksn)*dz2

		tmp1 =  F(i,j,k)-F(i,jm1,k)-F(im1,j,k)+F(im1,jm1,k)
		tmp2 =  F(i,j,k)-F(im1,j,k)-F(i,j,km1)+F(im1,j,km1)
		tmp3 =  F(i,j,k)-F(i,jm1,k)-F(i,j,km1)+F(i,jm1,km1)
        
        tmp210 =  GX(i,jm1,k)-GX(i,j,k)
        tmp120 =  GY(im1,j,k)-GY(i,j,k)
    	tmp201 =  GX(i,j,km1)-GX(i,j,k)
    	tmp102 =  GZ(im1,j,k)-GZ(i,j,k)
    	tmp021 =  GY(i,j,km1)-GY(i,j,k)
    	tmp012 =  GZ(i,jm1,k)-GZ(i,j,k)

    	a210 =  (-tmp1 -tmp210*DXA(ima1)*isn)*dx2*dy1*jsn      
    	a120 =  (-tmp1 -tmp120*DYA(jma1)*jsn)*dx1*dy2*isn  
    	a201 =  (-tmp2 -tmp201*DXA(ima1)*isn)*dx2*dz1*ksn
    	a102 =  (-tmp2 -tmp102*DZA(kma1)*ksn)*dx1*dz2*isn  
    	a021 =  (-tmp3 -tmp021*DYA(jma1)*jsn)*dy2*dz1*ksn  
    	a012 =  (-tmp3 -tmp012*DZA(kma1)*ksn)*dy1*dz2*jsn 

    	a110 =  (-tmp1-tmp210*DXA(ima1)*isn-tmp120*DYA(jma1)*jsn)*dx1*dy1*isn*jsn
    	a101 =  (-tmp2-tmp201*DXA(ima1)*isn-tmp102*DZA(kma1)*ksn)*dx1*dz1*isn*ksn 
    	a011 =  (-tmp3-tmp021*DYA(jma1)*jsn-tmp012*DZA(kma1)*ksn)*dy1*dz1*jsn*ksn
        
		
		a111 =  ( F(i,j,k)-F(im1,j,k)-F(i,jm1,k)-F(i,j,km1)	&
				+ F(im1,jm1,k)+F(im1,j,km1)+F(i,jm1,km1)-F(im1,jm1,km1))*dx1*dy1*dz1*isn*jsn*ksn
!------------------------------------------------------------------------
        
        FA(i,j,k) =  ( (a300*xx + a210*yy + a201*zz + a200)*xx + a110*yy + GX(i,j,k) )*xx  &
                  +  ( (a030*yy + a120*xx + a021*zz + a020)*yy + a011*zz + GY(i,j,k) )*yy  &
                  +  ( (a003*zz + a102*xx + a012*yy + a002)*zz + a101*xx + GZ(i,j,k) )*zz  &
                  +   a111*xx*yy*zz + F(i,j,k)
		GXA(i,j,k)=  (3.d0*a300*xx + 2.d0*(a210*yy+a201*zz+a200))*xx + (a120*yy+a111*zz+a110)*yy + (a102*zz+a101)*zz + GX(i,j,k)
        GYA(i,j,k)=  (3.d0*a030*yy + 2.d0*(a120*xx+a021*zz+a020))*yy + (a210*xx+a111*zz+a110)*xx + (a012*zz+a011)*zz + GY(i,j,k)
		GZA(i,j,k)=  (3.d0*a003*zz + 2.d0*(a102*xx+a012*yy+a002))*zz + (a201*xx+a111*yy+a101)*xx + (a021*xx+a011)*yy + GZ(i,j,k)
 
  
  	end do
  	end do
  	end do
  
  	do k=1,nz
    do j=1,ny
    do i=1,nx
    	if( (int(itp).EQ.1) .AND. (i.EQ.nx) ) cycle
    	if( (int(itp).EQ.2) .AND. (j.EQ.ny) ) cycle
    	if( (int(itp).EQ.3) .AND. (k.EQ.nz) ) cycle
      
        F(i,j,k)=FA(i,j,k)

        GX(i,j,k)=GXA(i,j,k)-                                   &
     &          (GXA(i,j,k)*(UC(i+1,j,k)-UC(i-1,J,k))           &
     &          +GYA(i,j,k)*(VC(i+1,j,k)-VC(i-1,J,k))           &
     &          +GZA(i,j,k)*(WC(i+1,j,k)-WC(i-1,J,k)))*dt       &
     &          /( 3.d0/2.d0*(DXA(I)+DXA(I+1)) )

        GY(i,j,k)=GYA(i,j,k)-                                   &
     &          (GXA(i,j,k)*(UC(i,j+1,K)-UC(i,j-1,k))           &
     &          +GYA(i,j,k)*(VC(i,j+1,K)-VC(i,j-1,k))           &
     &          +GZA(i,j,k)*(WC(i,j+1,K)-WC(i,j-1,k)))*dt       &
     &          /( 3.d0/2.d0*(DYA(j)+DYA(j+1)) )

        GZ(i,j,k)=GZA(i,j,k)-                                   &
     &          (GXA(i,j,k)*(UC(i,j,k+1)-UC(i,j,k-1))           &
     &          +GYA(i,j,k)*(VC(i,j,k+1)-VC(i,j,k-1))           &
     &          +GZA(i,j,k)*(WC(i,j,k+1)-WC(i,j,k-1)))*dt       &
     &          /( 3.d0/2.d0*(DZA(k)+DZA(k+1)) )
  
	end do
    end do
    end do
    
    if (int(itp).EQ.1) then
        U=F
        UX=GX
        UY=GY
        UZ=GZ
    elseif (int(itp).EQ.2) then
        V=F
        VX=GX
        VY=GY
        VZ=GZ
    elseif (int(itp).EQ.3) then
        W=F
        WX=GX
        WY=GY
        WZ=GZ
    end if

end subroutine


End Module MOD_MOM_VSIAM
