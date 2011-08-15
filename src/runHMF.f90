program runHMF
  use HMF_module
  use ParseText
  implicit none
  
  type(HMF) :: H
  
  type(PTo) :: HCF
  type(datafile_h5) :: h5hmf

  integer :: Nx, Nv
  double precision :: width, bag, e0, m0, epsilon
  double precision :: vmax, DT
  integer :: i,m,t,t_top, n_images, t_images
  integer :: n_steps, n_top
  
  double precision :: norme, masse
  integer, parameter :: nc=3
  real*4 :: perim(nc), z(nc)
  real*4, allocatable :: data(:,:), x(:), y(:)
  logical :: do_conrec
  double precision :: com
  integer, parameter :: time_number = 11
  double precision :: vals(time_number)
  character(len=32) :: time_names(time_number), IC
  character(len=13) :: fname

! declaration for the routines allowing the inversion of m = m(width)
  external zbrent, ftbz
  double precision zbrent

  interface 
     subroutine conrec(d,ilb,iub,jlb,jub,x,y,nc,z,funit,s)
       real*4, intent(in) :: d(ilb:iub,jlb:jub)
       integer, intent(in) :: ilb, iub, jlb, jub
       real*4, intent(in) :: x(ilb:iub), y(jlb:jub)
       integer, intent(in) :: nc
       real*4, intent(in) :: z(nc)
       integer, intent(in) :: funit
       real*4, intent(out) :: s(nc)
     end subroutine conrec
  end interface
! declaration for the contour routine
!     subroutine conrec(d,ilb,iub,jlb,jub,x,y,nc,z, s)
!     real*8 d(ilb:iub,jlb:jub)  ! matrix of data to contour
!     integer ilb,iub,jlb,jub    ! index bounds of data matrix
!     real*8 x(ilb:iub)          ! data matrix column coordinates
!     real*8 y(jlb,jub)          ! data matrix row coordinates
!     integer nc                 ! number of contour levels
!     real*8 z(1:nc)             ! contour levels in increasing order
!     real*8 s(nc)               ! length of the contour, one length by contour
!  external conrec

! Parse config file and create HMF and datafile data

  call PTparse(HCF,'HMF_in',7)
  Nx = PTread_i(HCF,'Nx')
  Nv = PTread_i(HCF,'Nv')
  vmax = PTread_d(HCF,'vmax')
  n_steps = PTread_i(HCF,'n_steps')
  n_top = PTread_i(HCF,'n_top')
  DT = PTread_d(HCF,'DT')
  IC = PTread_s(HCF,'IC')
  n_images = PTread_i(HCF, 'n_images')
  do_conrec = PTread_l(HCF, 'perimeter')
  t_images = 1

  call newHMF(H,Nx,Nv,vmax, Nedf=PTread_i(HCF,'Nedf'), model=PTread_s(HCF, 'model'), epsilon=PTread_d(HCF, 'epsilon'), Hfield=PTread_d(HCF,'Hfield'))
  H%V%DT = DT
  call create_h5(h5hmf, H%V, PTread_s(HCF, 'out_file'), n_top, time_number)

  write(*,'(3f30.26)') H%V%xmin, get_x(H%V,H%V%Nx), H%V%xmax
  write(*,'(2f30.26)') sin(H%V%xmin), sin(H%V%xmax)

  time_names = (/'time    ', 'mass    ', 'energy  ', 'int     ', 'kin     ', 'momentum', 'Mx      ', 'My      ', 'I2      ', 'I3      ', 'entropy '/)
  call write_info_string_array_h5(h5hmf, 'time_names', time_names)


  write(*,*) 'launching HMF with ', Nx, 'x', Nv,  'points'
  write(*,*) 'is_ext', H%is_ext

  if (allocated(H%edf)) then
     call write_info_double_h5(h5hmf, 'e_min', H%e_min)
     call write_info_double_h5(h5hmf, 'e_max', H%e_max)
     call write_info_double_h5(h5hmf, 'de', H%de)
  end if

! Handle initial condition

  if (IC.eq.'waterbag') then
     width = PTread_d(HCF, 'width')
     bag = PTread_d(HCF, 'bag')
     call init_carre(H%V, width, bag)
     H%f0 = 1d0/(4*width*bag)
  else if (IC.eq.'wb_eps') then
     width = PTread_d(HCF, 'width')
     bag = PTread_d(HCF, 'bag')
     epsilon = PTread_d(HCF,'epsilon')
     call init_carre(H%V, width, bag, epsilon=epsilon)
     H%f0 = 1d0/(4*width*bag)
  else if (IC.eq.'wb_m0') then
     e0 = PTread_d(HCF, 'e0')
     m0 = PTread_d(HCF, 'm0')
     if (m0.lt.0.001d0) then
        width=PI
     else if (m0.gt.0.99d0) then
        stop 'too singular (close to 0) width for the waterbag'
     else
        width = zbrent(ftbz, 1d-8,4.d0*atan(1.d0), m0, 1d-6)
     end if
     bag = sqrt(6.d0*(e0-0.5d0+m0**2*0.5d0))
     call init_carre(H%V, width, bag)
     H%f0 = 1d0/(4*width*bag)
  else if (IC.eq.'wb_m0_eps') then
     e0 = PTread_d(HCF, 'e0')
     m0 = PTread_d(HCF, 'm0')
     epsilon = PTread_d(HCF,'epsilon')
     if (m0.lt.0.001d0) then
        width=PI
     else if (m0.gt.0.99d0) then
        stop 'too singular (close to 0) width for the waterbag'
     else
        width = zbrent(ftbz, 1d-8,4.d0*atan(1.d0), m0, 1d-6)
     end if
     bag = sqrt(6.d0*(e0-0.5d0+m0**2*0.5d0))
     call init_carre(H%V, width, bag,epsilon=epsilon)
     H%f0 = 1d0/(4*width*bag)
     
     
  else if (IC.eq.'gaussian') then
     e0 = PTread_d(HCF, 'e0')
     if (e0.le.0.5d0) then
        stop 'e0 lt 0.5 in gaussian IC'
     end if
     call init_gaussian(H%V, 0.5d0/(e0-0.5d0))
     write(*,*) sum(H%V%f)*H%V%dx*H%V%dv
  else if (IC.eq.'gauss_beta_eps') then
     call init_gaussian(H%V, beta=PTread_d(HCF, 'beta'), epsilon=PTread_d(HCF,'epsilon'))
     write(*,*) sum(H%V%f(1:H%V%Nx,:))*H%V%dx*H%V%dv-1.d0
  else if (IC.eq.'gaussian_eps') then
     e0 = PTread_d(HCF, 'e0')
     epsilon = PTread_d(HCF,'epsilon')
     if (e0.le.0.5d0) then
        stop 'e0 lt 0.5 in gaussian IC'
     end if
     call init_gaussian(H%V, 0.5d0/(e0-0.5d0),epsilon=epsilon)
     write(*,*) sum(H%V%f)*H%V%dx*H%V%dv
  else if (IC.eq.'fermi_eps') then
     width = PTread_d(HCF, 'beta')
     e0 = PTread_d(HCF, 'a')
     epsilon = PTread_d(HCF, 'epsilon')
     do i=1,H%V%Nx
        do m=1,H%V%Nv
           H%V%f(i,m) = 1.d0/(1.d0+exp(width * (get_v(H%V,m)**2 - e0)/2.d0)) * (1.d0 + epsilon*sin(get_x(H%V,i)))
        end do
     end do

     norme = sum(H%V%f(1:H%V%Nx,:))* H%V%dx * H%V%dv
     H%V%f = H%V%f / norme
     write(*,*) sum(H%V%f(1:H%V%Nx,:))*H%V%dx*H%V%dv-1.d0
  else if (IC.eq.'two_streams') then
     call init_two_streams(H%V, PTread_d(HCF,'v_min'), PTread_d(HCF, 'v_max'))
  else if (IC.eq.'two_streams_width') then
     call init_two_streams(H%V, PTread_d(HCF,'v_min'), PTread_d(HCF, 'v_max'), width=PTread_d(HCF,'width'))
  else
     stop 'unknown IC'
  end if

  call PTkill(HCF)

  call write_grid_info(h5hmf, H%V, 'HMF')

! fin de la condition initiale

!
! data is a real*4 array that will be passed to conrec, the contouring subroutine
! x and y are real*4 array for the x and v position arrays to be passed to conrec
! z(:) is a list of values at which the perimeter is to be computed
!

  if (do_conrec) then
     allocate(data(0:H%V%Nx-1,0:H%V%Nv-1))
     allocate(x(0:H%V%Nx-1)) ; allocate(y(0:H%V%Nv-1))
     data = 0.
     x = (/ (real(H%V%xmin)+real(i)*real(H%V%dx), i=0, H%V%Nx-1) /)
     y = (/ (real(H%V%vmin)+real(m)*real(H%V%dv), m=0, H%V%Nv-1) /)
     z(1) = real(H%f0*0.25d0)
     z(2) = real(H%f0*0.5d0)
     z(3) = real(H%f0*0.75d0)

     open(17, file='t_perim')
  end if
  
! Simulation

!
! Computing of the macroscopic quantities of the initial condition, before entering the simulation loop
!

  t_top=0
  call compute_rho(H%V)
  call compute_M(H)
  call compute_phi(H%V)
  call write_data_group_h5(h5hmf, H%V, 0)
  call compute_phys(H, vals, t_top*DT)
  call write_time_slice_h5(h5hmf, t_top, vals)
  if (allocated(H%edf)) then
     call compute_edf(H)
     call write_edf_h5(h5hmf, H, t_top)
  end if

  if (do_conrec) then
     write(fname, '(a8,i5.5)') 'images/c', t_top
     open(18, file=fname)
     data(0:H%V%Nx-1,0:H%V%Nv-1) = real(H%V%f(1:H%V%Nx,1:H%V%Nv))
     call conrec(data, 0, H%V%Nx-1, 0,H%V%Nv-1, x, y, nc, z, 18, perim)
     close(18)
     write(17,'(4e20.12)') t_top*n_steps*DT, perim(1), perim(2), perim(3)
  end if


!
! Simulation loop
!

  do t_top = 1,n_top

     call spline_x(H%V)
     call advection_x_demi(H%V)
     H%V%f = H%V%g
     do t=1,n_steps-1
        
        call compute_rho(H%V)
        call compute_force(H)
        call spline_v(H%V)
        call advection_v(H%V)
        H%V%f = H%V%g
     
        call spline_x(H%V)
        call advection_x(H%V)
        H%V%f = H%V%g
     end do

     call compute_rho(H%V)
     call compute_force(H)
     call spline_v(H%V)
     call advection_v(H%V)
     H%V%f = H%V%g
     call spline_x(H%V)
     call advection_x_demi(H%V)
     H%V%f = H%V%g

     call compute_rho(H%V)
     call compute_M(H)
     call compute_phys(H, vals, t_top*n_steps*DT)
     call compute_phi(H%V)
     if (t_top*n_images/n_top.ge.t_images) then
        call write_data_group_h5(h5hmf, H%V, t_top)
        t_images = t_images + 1
        if (allocated(H%edf)) then
           call compute_edf(H)
           call write_edf_h5(h5hmf, H, t_top)
        end if
        if (do_conrec) then
           write(fname, '(a8,i5.5)') 'images/c', t_top
           open(18, file=fname)
           data(0:H%V%Nx-1,0:H%V%Nv-1) = real(H%V%f(1:H%V%Nx,1:H%V%Nv))
           call conrec(data, 0, H%V%Nx-1, 0,H%V%Nv-1, x, y, nc, z, 18, perim)
           close(18)
        end if
     end if
     call write_time_slice_h5(h5hmf, t_top, vals)

     if (do_conrec) then
        data(0:H%V%Nx-1,0:H%V%Nv-1) = real(H%V%f(1:H%V%Nx,1:H%V%Nv))
        call conrec(data, 0, H%V%Nx-1, 0,H%V%Nv-1, x, y, nc, z, 0, perim)

        write(17,'(4e20.12)') t_top*n_steps*DT, perim(1), perim(2), perim(3)
     end if

     
!     call h5fflush_f(h5hmf%file_id, H5F_SCOPE_GLOBAL_F, h5hmf%error)

  end do

!
! End of the simulation loop
!

  call close_h5(h5hmf)
  if (do_conrec) then
     close(17)
  end if

end program runHMF

function ftbz(s_theta,bunching)
  implicit none
  double precision s_theta, bunching, ftbz

  ftbz = sin(s_theta)/s_theta - bunching

end function ftbz
