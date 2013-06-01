! Copyright (C) 2009-2011 Pierre de Buyl

! This file is part of vmf90

! vmf90 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! vmf90 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with vmf90.  If not, see <http://www.gnu.org/licenses/>.

!> This program runs a Vlasov simulation of the Hamiltonian Mean-Field model, using
!! the vmf90 code.
program runHMF
  use HMF_module
  use ParseText
  use vmf90
  implicit none
  
  type(HMF) :: H
  
  type(PTo) :: HCF
  type(datafile_h5) :: h5hmf

  integer :: Nx, Nv
  double precision :: width, bag, e0, m0, epsilon
  double precision :: beta, fzero, ics
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

  call newHMF(H,Nx,Nv,vmax, Nedf=PTread_i(HCF,'Nedf'), model=PTread_s(HCF, 'model'), Hfield=PTread_d(HCF,'Hfield') )
  H%V%DT = DT
  call create_h5(h5hmf, H%V, PTread_s(HCF, 'out_file'), n_top, time_number)


  time_names = (/'time    ', 'mass    ', 'energy  ', 'int     ', &
                 'kin     ', 'momentum', 'Mx      ', 'My      ', &
                 'I2      ', 'I3      ', 'entropy '/)
  call write_info_string_array_h5(h5hmf, 'time_names', time_names)

  call vmf90_info()
  write(*,*) 'Running HMF with ', Nx, 'x', Nv,  'points'

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
  else if (IC.eq.'wb_selfc') then
     call init_wb_selfc(H%V, PTread_d(HCF, 'wb_y0'), PTread_d(HCF, 'wb_mx'))
  else if (IC.eq.'gaussian') then
     call init_gaussian(H%V, beta=PTread_d(HCF, 'gaussian_beta'))
  else if (IC.eq.'gaussian_eps') then
     call init_gaussian(H%V, beta=PTread_d(HCF, 'gaussian_beta'), &
          epsilon=PTread_d(HCF,'gaussian_epsilon'))
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
  else if (IC.eq.'LB') then
     beta = PTread_d(HCF, 'beta')
     fzero = PTread_d(HCF, 'fzero')
     ics = PTread_d(HCF, 'ics')
     m0 = PTread_d(HCF, 'm0')
     do i=1,H%V%Nx
        do m=1,H%V%Nv
           H%V%f(i,m) = fzero/(1.d0+ &
           exp( beta*(0.5d0*get_v(H%V,m)**2 - m0*cos(get_x(H%V,i)) )) / ics &
                )
        end do
     end do
     norme = sum(H%V%f(1:H%V%Nx,:))* H%V%dx * H%V%dv
     H%V%f = H%V%f / norme
  else if (IC.eq.'squared_lorentzian') then
     call init_squared_lorentzian(H%V, PTread_d(HCF,'lorentz_gamma'), &
          epsilon=PTread_d(HCF,'lorentz_epsilon'))
  else if (IC.eq.'from_file') then
     call load_data_from_h5(H%V,PTread_s(HCF,'IC_file'),trim(PTread_s(HCF,'IC_position')))
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

  end do

!
! End of the simulation loop
!

  call close_h5(h5hmf)

end program runHMF
