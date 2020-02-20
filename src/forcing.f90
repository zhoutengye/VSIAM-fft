module mod_forcing
  use mod_types
  use mod_common_mpi, only: ierr
  implicit none
  private
  public force_vel,force_bulk_vel,chkmean_ibm
  contains
  subroutine force_vel(n,dl,dzc,dzf,l,psi_u,psi_v,psi_w,u,v,w,f)
    !
    ! force velocity field using the volume-penalisation IBM:
    ! the force is proportional to the volume fraction of
    ! solid in a grid cell i,j,k.
    ! the volume fraction is defined in the cell centers and
    ! is interpolated to the cell faces where the velocity
    ! components are defined. This results in a smoothing of the
    ! local volume fraction field.
    ! 
    ! note: for some systems it may be convenient to save the force
    !       distribution
    !
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in   ) , dimension(3) :: dl,l
    real(rp), intent(in   ) , dimension(0:) :: dzc,dzf
    real(rp), intent(in   ) , dimension(0:,0:,0:) :: psi_u,psi_v,psi_w
    real(rp), intent(inout) , dimension(0:,0:,0:) :: u,v,w
    real(rp), intent(inout), dimension(3) :: f
    real(rp) :: psix,psiy,psiz,fx,fy,fz,fxtot,fytot,fztot
    integer :: i,j,k!,ip,jp,kp
    !
    fxtot = 0.
    fytot = 0.
    fztot = 0.
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,u,v,w,psi_u,psi_v,psi_w,dl,dzc,dzf) &
    !$OMP PRIVATE(i,j,k,fx,fy,fz,psix,psiy,psiz) &
    !$OMP REDUCTION(+:fxtot,fytot,fztot)
    do k=1,n(3)
      !kp = k+1
      do j=1,n(2)
        !jp = j+1
        do i=1,n(1)
          !ip = i+1
          !
          ! note: these loops should be splitted into three for more efficient
          !       memory access
          !
          psix  = psi_u(i,j,k)!0.5*(psi(ip,j ,k )+psi(i,j,k)) ! liquid volume fraction
          psiy  = psi_v(i,j,k)!0.5*(psi(i ,jp,k )+psi(i,j,k)) ! liquid volume fraction
          psiz  = psi_w(i,j,k)!0.5*(psi(i ,j ,kp)+psi(i,j,k)) ! liquid volume fraction
          fx    = - u(i,j,k)*psix ! (u(i,j,k)*(1.-psix)-u(i,j,k))*dti
          fy    = - v(i,j,k)*psiy ! (v(i,j,k)*(1.-psiy)-v(i,j,k))*dti
          fz    = - w(i,j,k)*psiz ! (w(i,j,k)*(1.-psiz)-w(i,j,k))*dti
          u(i,j,k) = u(i,j,k) + fx
          v(i,j,k) = v(i,j,k) + fy
          w(i,j,k) = w(i,j,k) + fz
          fxtot = fxtot + fx*dl(1)*dl(2)*dzf(k)
          fytot = fytot + fy*dl(1)*dl(2)*dzf(k)
          fztot = fztot + fz*dl(1)*dl(2)*dzc(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,f(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    f(1) = fxtot/(l(1)*l(2)*l(3))
    f(2) = fytot/(l(1)*l(2)*l(3))
    f(3) = fztot/(l(1)*l(2)*l(3))
    return
  end subroutine force_vel
  !
  subroutine force_bulk_vel(n,idir,dl,dzc,dzf,l,psi,p,velf,f)
    !
    ! bulk velocity forcing only in a region of the domain
    ! where psi is non-zero
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   ) :: idir
    real(rp), intent(in   ) , dimension(3) :: dl,l
    real(rp), intent(in   ) , dimension(0:) :: dzc,dzf
    real(rp), intent(in   ), dimension(0:,0:,0:) :: psi
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp), intent(in   ) :: velf
    real(rp), intent(out  ) :: f
    real(rp)              :: mean_val,mean_psi
    integer :: i,j,k
    integer, dimension(3) :: q
    real(rp) :: psis
    real(rp), dimension(n(3)) :: dz
    !
    select case(idir)
    case(0)
      q = (/0,0,0/)
      dz(:) = dzf(1:n(3))
    case(1)
      q = (/1,0,0/)
      dz(:) = dzf(1:n(3))
    case(2)
      q = (/0,1,0/)
      dz(:) = dzf(1:n(3))
    case(3)
      q = (/0,0,1/)
      dz(:) = dzc(1:n(3))
    end select
    mean_val = 0.
    mean_psi = 0.
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,q,psi,p,dl,dz) &
    !$OMP PRIVATE(i,j,k,psis) &
    !$OMP REDUCTION(+:mean_val,mean_psi)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !psis = 1.-0.5*(psi(i,j,k)+psi(i+q(1),j+q(2),k+q(3)))
          psis = 1.-psi(i,j,k) ! psi already shifted according to staggered position
          mean_val = mean_val + p(i,j,k)*psis*dl(1)*dl(2)*dz(k)
          mean_psi = mean_psi + psis*dl(1)*dl(2)*dz(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,mean_val,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,mean_psi,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    mean_val = mean_val/mean_psi
    f = 0.
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,q,psi,p,mean_val,velf,dl,dz) &
    !$OMP PRIVATE(i,j,k,psis) &
    !$OMP REDUCTION(+:f)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !psis = 1.-0.5*(psi(i,j,k)+psi(i+q(1),j+q(2),k+q(3)))
          !psis = 1.-psi(i,j,k) ! (if bulk velocity forced only inside the fluid)
          psis = 1.
          p(i,j,k) = p(i,j,k) + (velf-mean_val)*psis
          f = f + (velf-mean_val)*psis*dl(1)*dl(2)*dz(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,f,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    f = f/(l(1)*l(2)*l(3))
    !f = f/(l(1)*l(2)*l(3)-mean_psi) ! (if bulk velocity forced only inside the fluid)
    return
  end subroutine force_bulk_vel
  subroutine chkmean_ibm(n,idir,dl,dzc,dzf,l,psi,p,mean)
    !
    ! bulk velocity forcing only in a region of the domain
    ! where psi is non-zero
    !
    implicit none
    integer , intent(in   ), dimension(3) :: n
    integer , intent(in   ) :: idir
    real(rp), intent(in   ) , dimension(3) :: dl,l ! l not used
    real(rp), intent(in   ) , dimension(0:) :: dzc,dzf
    real(rp), intent(in   ), dimension(0:,0:,0:) :: psi
    real(rp), intent(inout), dimension(0:,0:,0:) :: p
    real(rp), intent(out  ) :: mean
    real(rp)              :: mean_val,mean_psi
    integer :: i,j,k
    integer, dimension(3) :: q
    real(rp) :: psis
    real(rp), dimension(n(3)) :: dz
    !
    select case(idir)
    case(0)
      q = (/0,0,0/)
      dz(:) = dzf(1:n(3))
    case(1)
      q = (/1,0,0/)
      dz(:) = dzf(1:n(3))
    case(2)
      q = (/0,1,0/)
    case(3)
      q = (/0,0,1/)
      dz(:) = dzc(1:n(3))
    end select
    mean_val = 0.
    mean_psi = 0.
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(n,q,psi,p,dl,dz) &
    !$OMP PRIVATE(i,j,k,psis) &
    !$OMP REDUCTION(+:mean_val,mean_psi)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          !psis = 1.-0.5*(psi(i,j,k)+psi(i+q(1),j+q(2),k+q(3)))
          psis = 1.-psi(i,j,k)
          mean_val = mean_val + p(i,j,k)*psis*dl(1)*dl(2)*dz(k)
          mean_psi = mean_psi + psis*dl(1)*dl(2)*dz(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,mean_val,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    call mpi_allreduce(MPI_IN_PLACE,mean_psi,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    mean = mean_val/mean_psi
    return
  end subroutine chkmean_ibm
end module mod_forcing
