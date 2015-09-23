! Copyright (C) 2009-2015 Pierre de Buyl

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
  use h5md
  implicit none
  
  type(HMF) :: H
  
  type(PTo) :: HCF

  integer :: Nx, Nv
  double precision :: width, bag, e0, m0, epsilon
  double precision :: beta, fzero, ics
  double precision :: vmax, DT
  integer :: i,m,t,t_top, n_images, t_images
  integer :: n_pmoments, n_hmoments
  integer :: n_steps, n_top
  
  double precision :: norme
  double precision :: realtime
  character(len=32) :: IC

  integer(HID_T) :: file_ID
  type(h5md_t) :: mass_ID, energy_ID, int_ID, kin_ID, momentum_ID, Mx_ID, My_ID, I2_ID, I3_ID, pn_ID, hn_ID
  type(h5md_t) :: edf_ID
  type(h5md_t) :: f_ID, rho_ID, phi_ID

  call h5open_f(h5_error)

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
  n_pmoments = PTread_i(HCF, 'n_pmoments')
  n_hmoments = PTread_i(HCF, 'n_hmoments')

  call newHMF(H,Nx,Nv,vmax, Nedf=PTread_i(HCF,'Nedf'), model=PTread_s(HCF, 'model'), &
       Hfield=PTread_d(HCF,'Hfield'), n_pmoments=n_pmoments, n_hmoments=n_hmoments)
  H%V%DT = DT

  call h5md_create_file(file_ID, 'hmf.h5', 'Pierre de Buyl <pdebuyl@ulb.ac.be>', 'vmf90_hmf', trim(vmf90_version()))

  call vmf90_info()
  write(*,*) 'Running HMF with ', Nx, 'x', Nv,  'points'

  if (allocated(H%edf)) then
     call h5md_write_par(file_ID, 'e_min', H%e_min)  
     call h5md_write_par(file_ID, 'e_max', H%e_max)  
     call h5md_write_par(file_ID, 'de', H%de)  
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
     !call load_data_from_h5(H%V,PTread_s(HCF,'IC_file'),trim(PTread_s(HCF,'IC_position')))
     stop 'no reload allowed yet'
  else
     stop 'unknown IC'
  end if

  call PTkill(HCF)

  call write_grid_info(file_ID, H%V, 'HMF')

! fin de la condition initiale

  
! Simulation

!
! Computing of the macroscopic quantities of the initial condition, before entering the simulation loop
!

  call begin_h5md
  t_top=0
  realtime=0.d0
  call compute_rho(H%V)
  call compute_M(H)
  call compute_phi(H%V)
  call write_fields
  call compute_phys(H)
  call write_obs
  if (allocated(H%edf)) then
     call compute_edf(H)
     call h5md_write_obs(edf_id, H%edf, t_top, realtime)
  end if


!
! Simulation loop
!

  do t_top = 1,n_top

     call advance_x(H%V, 0.5d0)
     do t=1,n_steps-1
        
        call compute_rho(H%V)
        call compute_force(H)
        call advance_v(H%V, 1.d0)
     
        call advance_x(H%V, 1.d0)

        realtime = realtime + DT
     end do

     call compute_rho(H%V)
     call compute_force(H)
     call advance_v(H%V, 1.d0)
     call advance_x(H%V, 0.5d0)
     realtime = realtime + DT

     call compute_rho(H%V)
     call compute_M(H)
     call compute_phys(H)
     call compute_phi(H%V)
     if (t_top*n_images/n_top.ge.t_images) then
        call write_fields
        t_images = t_images + 1
        if (allocated(H%edf)) then
           call compute_edf(H)
           call h5md_write_obs(edf_id, H%edf, t_top, realtime)
        end if
     end if
     call write_obs

  end do

!
! End of the simulation loop
!

  call h5close_f(h5_error)

contains

  subroutine begin_h5md
    integer :: i
    call h5md_create_obs(file_ID, 'mass', mass_ID, H%V%masse)
    call h5md_create_obs(file_ID, 'energy', energy_ID, H%V%energie, link_from='mass')
    call h5md_create_obs(file_ID, 'en_int', int_ID, H%V%en_int, link_from='mass')
    call h5md_create_obs(file_ID, 'en_kin', kin_ID, H%V%en_kin, link_from='mass')
    call h5md_create_obs(file_ID, 'momentum', momentum_ID, H%V%momentum, link_from='mass')
    call h5md_create_obs(file_ID, 'Mx', Mx_ID, H%Mx, link_from='mass')
    call h5md_create_obs(file_ID, 'My', My_ID, H%My, link_from='mass')
    call h5md_create_obs(file_ID, 'I2', I2_ID, H%I2, link_from='mass')
    call h5md_create_obs(file_ID, 'I3', I3_ID, H%I3, link_from='mass')
    call h5md_create_obs(file_ID, 'edf', edf_ID, H%edf)
    if (allocated(H%pn)) then
       call h5md_create_obs(file_ID, 'pn', pn_ID, H%pn, link_from='mass')
    end if
    if (allocated(H%hn)) then
       call h5md_create_obs(file_ID, 'hn', hn_ID, H%hn, link_from='mass')
    end if

    call create_fields_group(file_ID)
    call h5md_create_obs(file_ID, 'f', f_ID, H%V%f, override_obs='fields')
    call h5md_create_obs(file_ID, 'rho', rho_ID, H%V%rho, link_from='f', override_obs='fields')
    call h5md_create_obs(file_ID, 'phi', phi_ID, H%V%phi, link_from='f', override_obs='fields')
  end subroutine begin_h5md

  subroutine write_obs
    call h5md_write_obs(mass_ID, H%V%masse, t_top, realtime)
    call h5md_write_obs(energy_ID, H%V%energie, t_top, realtime)
    call h5md_write_obs(int_ID, H%V%en_int, t_top, realtime)
    call h5md_write_obs(kin_ID, H%V%en_kin, t_top, realtime)
    call h5md_write_obs(momentum_ID, H%V%momentum, t_top, realtime)
    call h5md_write_obs(Mx_ID, H%Mx, t_top, realtime)
    call h5md_write_obs(My_ID, H%My, t_top, realtime)
    call h5md_write_obs(I2_ID, H%I2, t_top, realtime)
    call h5md_write_obs(I3_ID, H%I3, t_top, realtime)
    if (allocated(H%pn)) then
       call h5md_write_obs(pn_ID, H%pn, t_top, realtime)
    end if
    if (allocated(H%hn)) then
       call h5md_write_obs(hn_ID, H%hn, t_top, realtime)
    end if
  end subroutine write_obs

  subroutine write_fields
    call h5md_write_obs(f_ID, H%V%f, t_top, realtime)
    call h5md_write_obs(rho_ID, H%V%rho, t_top, realtime)
    call h5md_write_obs(phi_ID, H%V%phi, t_top, realtime)
  end subroutine write_fields

end program runHMF
