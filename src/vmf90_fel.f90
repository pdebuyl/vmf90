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

!> This program runs a Vlasov simulation of the Colson-Bonifacio model for the free electron
!! laser. It makes use of the vmf90 code.
program runFEL
  use FEL_module
  use ParseText
  use vmf90
  use h5md
  implicit none
  
  type(FEL) :: F
  
  type(PTo) :: HCF
  integer(HID_T) :: file_ID
  type(h5md_t) :: mass_ID, energy_ID, int_ID, kin_ID, momentum_ID, I_ID, varphi_ID, &
       entropy_ID, Ax_ID, Ay_ID, L2_ID
  type(h5md_t) :: f_ID, rho_ID, phi_ID

  integer :: Nx, Nv
  double precision :: width, bag, vmin, vmax, DT, Mx, My, p0
  integer :: t,t_top, n_images, t_images
  integer :: n_steps, n_top
  double precision :: realtime
  
  character(len=12) :: IC
  
! Parse config file and create HMF and datafile data

  call PTparse(HCF,'FEL_in',7)
  Nx = PTread_i(HCF,'Nx')
  Nv = PTread_i(HCF,'Nv')
  vmax = PTread_d(HCF,'vmax')
  vmin = PTread_d(HCF,'vmin')
  n_steps = PTread_i(HCF,'n_steps')
  n_top = PTread_i(HCF,'n_top')
  DT = PTread_d(HCF,'DT')
  IC = PTread_s(HCF, 'IC')
  n_images = PTread_i(HCF, 'n_images')
  t_images = 1

  call newFEL(F,Nx,Nv,vmax,vmin, delta=PTread_d(HCF,'delta'))
  F%V%DT = DT
  F%phi = 0.d0
  F%I = 0.0001d0
  
  call vmf90_info()
  write(*,*) 'launching FEL with ', Nx, 'x', Nv,  'points'
  F%Ax = 0.0001d0;
  F%Ay = 0.d0;

! Handle initial condition

  if (IC.eq.'waterbag') then
     width = PTread_d(HCF, 'width')
     bag = PTread_d(HCF, 'bag')
     call init_carre(F%V, width, bag)
     F%f0 = .25d0/(width*bag)
  else if (IC.eq.'wb_p0') then
     width = PTread_d(HCF, 'width')
     bag = PTread_d(HCF, 'bag')
     p0 = PTread_d(HCF, 'p0')
     call init_carre(F%V, width, bag, p0 = p0)
     F%f0 = .25d0/(width*bag)
  else if (IC.eq.'gaussian_p0') then
     call init_gaussian(F%V, PTread_d(HCF, 'beta'), p0=PTread_d(HCF,'p0'))
  else
     stop 'unknown IC'
  end if

  call h5open_f(h5_error)

  call h5md_create_file(file_ID, 'fel.h5', 'Pierre de Buyl <pdebuyl@ulb.ac.be>', 'vmf90_fel', 'No version information yet')

  call write_grid_info(file_ID, F%V, 'FEL')
  call h5md_write_par(file_ID,'delta',F%delta)


  call PTkill(HCF)

  call begin_h5md

  t_top=0

  call compute_rho(F%V)
  call compute_M(F)
  call compute_phys(F)
  realtime = 0.d0
  call write_obs
  call write_fields

  do t_top = 1,n_top

     call advance_x(F%V, 0.5d0)

     call compute_rho(F%V)
     call compute_M(F)

     do t=1,n_steps-1
        
        Mx = F%Mx
        My = F%My

        call compute_force_A(F)
        call advance_v(F%V, 0.5d0)
        
        call advance_x(F%V, 1.d0)

        call compute_rho(F%V)
        call compute_M(F)

        F%Ax(1) = F%Ax(1) + 0.5d0 * (Mx+F%Mx) * F%V%DT
        F%Ay(1) = F%Ay(1) - 0.5d0 * (My+F%My) * F%V%DT
     
     end do

     Mx = F%Mx
     My = F%My

     call compute_force_A(F)
     call advance_v(F%V,1.d0)

     call advance_x(F%V, 0.5d0)

     realtime = n_steps*DT*t_top

     call compute_rho(F%V)
     call compute_M(F)

     F%Ax(1) = F%Ax(1) + F%Mx * F%V%DT 
     F%Ay(1) = F%Ay(1) - F%My * F%V%DT

     F%I = F%Ax(1)**2+F%Ay(1)**2
     F%phi = atan2(F%Ay(1),F%Ax(1))

     call compute_phys(F)
     call write_obs

     if (t_top*n_images/n_top.ge.t_images) then
        call compute_phi(F%V)
        call write_fields
        t_images = t_images + 1
     end if
     
  end do

  call h5close_f(h5_error)

contains

  subroutine begin_h5md
    call h5md_create_obs(file_ID, 'mass', mass_ID, F%V%masse)
    call h5md_create_obs(file_ID, 'energy', energy_ID, F%V%energie, link_from='mass')
    call h5md_create_obs(file_ID, 'en_int', int_ID, F%V%en_int, link_from='mass')
    call h5md_create_obs(file_ID, 'en_kin', kin_ID, F%V%en_kin, link_from='mass')
    call h5md_create_obs(file_ID, 'momentum', momentum_ID, F%V%momentum, link_from='mass')
    call h5md_create_obs(file_ID, 'I', I_ID, F%I, link_from='mass')
    call h5md_create_obs(file_ID, 'phi', varphi_ID, F%phi, link_from='mass')
    call h5md_create_obs(file_ID, 'entropy', entropy_ID, F%entropy, link_from='mass')
    call h5md_create_obs(file_ID, 'Ax', Ax_ID, F%Ax(1), link_from='mass')
    call h5md_create_obs(file_ID, 'Ay', Ay_ID, F%Ay(1), link_from='mass')
    call h5md_create_obs(file_ID, 'L2', L2_ID, F%L2, link_from='mass')

    call create_fields_group(file_ID)
    call h5md_create_obs(file_ID, 'f', f_ID, F%V%f, override_obs='fields')
    call h5md_create_obs(file_ID, 'rho', rho_ID, F%V%rho, link_from='f', override_obs='fields')
    call h5md_create_obs(file_ID, 'phi', phi_ID, F%V%phi, link_from='f', override_obs='fields')
  end subroutine begin_h5md

  subroutine write_obs
    call h5md_write_obs(mass_ID, F%V%masse, t_top, realtime)
    call h5md_write_obs(energy_ID, F%V%energie, t_top, realtime)
    call h5md_write_obs(int_ID, F%V%en_int, t_top, realtime)
    call h5md_write_obs(kin_ID, F%V%en_kin, t_top, realtime)
    call h5md_write_obs(momentum_ID, F%V%momentum, t_top, realtime)
    call h5md_write_obs(I_ID, F%I, t_top, realtime)
    call h5md_write_obs(varphi_ID, F%phi, t_top, realtime)
    call h5md_write_obs(entropy_ID, F%entropy, t_top, realtime)
    call h5md_write_obs(Ax_ID, F%Ax(1), t_top, realtime)
    call h5md_write_obs(Ay_ID, F%Ay(1), t_top, realtime)
    call h5md_write_obs(L2_ID, F%L2, t_top, realtime)
  end subroutine write_obs

  subroutine write_fields
    call h5md_write_obs(f_ID, F%V%f, t_top, realtime)
    call h5md_write_obs(rho_ID, F%V%rho, t_top, realtime)
    call h5md_write_obs(phi_ID, F%V%phi, t_top, realtime)
  end subroutine write_fields

end program runFEL
