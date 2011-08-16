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
  double precision :: com
  integer, parameter :: time_number = 11
  double precision :: vals(time_number)
  character(len=32) :: time_names(time_number), IC
  character(len=13) :: fname

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
  else
     stop 'unknown IC'
  end if

  call PTkill(HCF)

  call write_grid_info(h5hmf, H%V, 'HMF')

! fin de la condition initiale

  
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
     end if
     call write_time_slice_h5(h5hmf, t_top, vals)

!     call h5fflush_f(h5hmf%file_id, H5F_SCOPE_GLOBAL_F, h5hmf%error)

  end do

!
! End of the simulation loop
!

  call close_h5(h5hmf)

end program runHMF
