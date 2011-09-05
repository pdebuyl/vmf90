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
  implicit none
  
  type(FEL) :: F
  
  type(PTo) :: HCF
  type(datafile_h5) :: h5fel

  integer :: Nx, Nv
  double precision :: width, bag, vmin, vmax, DT, x1, x2, x3, Mx, My, p0
  integer :: i,m,t,t_top, n_images, t_images
  integer :: n_steps, n_top
  
  integer, parameter :: time_number = 12
  double precision :: vals(time_number)
  character(len=12) :: time_names(time_number), IC
  
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
  else
     stop 'unknown IC'
  end if

  call create_h5(h5fel, F%V, PTread_s(HCF,'out_file'), n_top, time_number)
  call write_grid_info(h5fel, F%V, 'FEL')
  call write_info_double_h5(h5fel,'delta',F%delta)

  time_names = (/'time', 'mass', 'energy', 'int', 'kin', 'momentum', 'I', 'phi', 'entropy','Ax','Ay', 'L2'/)
  call write_info_string_array_h5(h5fel, 'time_names', time_names)

  call PTkill(HCF)

  t_top=0

  call compute_rho(F%V)
  call compute_M(F)
  call compute_phys(F, vals, t_top*n_steps*DT)
  call write_data_group_h5(h5fel, F%V, 0)
  call write_time_slice_h5(h5fel, t_top, vals)


  do t_top = 1,n_top

     call spline_x(F%V)
     call advection_x_demi(F%V)
     F%V%f = F%V%g

     call compute_rho(F%V)
     call compute_M(F)

     do t=1,n_steps-1
        
        Mx = F%Mx
        My = F%My

        call compute_force_A(F)
        call spline_v(F%V)
        call advection_v(F%V)
        F%V%f = F%V%g
        
        call spline_x(F%V)
        call advection_x(F%V)
        F%V%f = F%V%g

        call compute_rho(F%V)
        call compute_M(F)

        F%Ax(1) = F%Ax(1) + 0.5d0 * (Mx+F%Mx) * F%V%DT
        F%Ay(1) = F%Ay(1) - 0.5d0 * (My+F%My) * F%V%DT
     
     end do

     Mx = F%Mx
     My = F%My

     call compute_force_A(F)
     call spline_v(F%V)
     call advection_v(F%V)
     F%V%f = F%V%g

     call spline_x(F%V)
     call advection_x_demi(F%V)
     F%V%f = F%V%g

     call compute_rho(F%V)
     call compute_M(F)

     F%Ax(1) = F%Ax(1) + F%Mx * F%V%DT 
     F%Ay(1) = F%Ay(1) - F%My * F%V%DT

     F%I = F%Ax(1)**2+F%Ay(1)**2
     F%phi = atan2(F%Ay(1),F%Ax(1))

     call compute_phys(F, vals, t_top*n_steps*DT)
     call write_time_slice_h5(h5fel, t_top, vals)

     if (t_top*n_images/n_top.ge.t_images) then
        call write_data_group_h5(h5fel, F%V, t_top)
        t_images = t_images + 1
     end if
     
  end do

  call close_h5(h5fel)

end program runFEL
