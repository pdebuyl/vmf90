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

!> This module defines a Vlasov setup for a variation on the Hamiltonian Mean-Field model in which
!! the system is divided in two subsystem of variable importance.
!!
!! Routines are provided to compute the force field and the thermodynamical observables.

module AB_module
  use Vlasov_module
  implicit none

  type AB
     type(grid) :: A,B
     double precision :: gamma
     double precision :: MAx, MAy, MBx, MBy, Mx, My
     double precision :: fA0, fB0
  end type AB

contains
  subroutine newAB(this, NAx, NAv, NBx, NBv, vmax, vmin)
    type(AB), intent(out) :: this
    integer, intent(in) :: NAx, NAv, NBx, NBv
    double precision, intent(in) :: vmax
    double precision, intent(in), optional :: vmin

    double precision :: vmin_final

    if (present(vmin)) then
       vmin_final = vmin
    else
       vmin_final = -vmax
    end if

    call new(this%A, NAx, NAv, PI, vmax, vmin=vmin_final, is_periodic=.true.)
    call new(this%B, NBx, NBv, PI, vmax, vmin=vmin_final, is_periodic=.true.)
  end subroutine newAB

  subroutine create_AB_h5(this, thisAB, filename, ntime, nvals)
    type(datafile_h5), intent(out) :: this
    type(AB), intent(in) :: thisAB
    character(len=*), intent(in) :: filename
    integer, intent(in) :: ntime, nvals

    integer(HSIZE_T) :: ts_one_dims(1)

    call h5open_f(this%error)
    this%filename=filename
    call h5fcreate_f(this%filename, H5F_ACC_TRUNC_F, this%file_id, this%error)
    
    call h5gcreate_f(this%file_id, this%info_g_name, this%info_g_id, this%error)
    call h5gcreate_f(this%file_id, this%data_g_name, this%data_g_id, this%error)
    call h5gcreate_f(this%file_id, this%data_g2_name, this%data_g2_id, this%error)

    this%ntime = ntime
    this%nvals = nvals
    this%time_dims = (/ntime+1,nvals/)
    call h5screate_simple_f(this%time_rank, this%time_dims, this%time_s_id, this%error)
    call h5dcreate_f(this%file_id, this%time_d_name, H5T_NATIVE_DOUBLE, this%time_s_id, this%time_d_id, this%error)

    this%f_dims = (/thisAB%A%Nx,thisAB%A%Nv/)
    this%rho_dims = (/thisAB%A%Nx/)
    this%phi_dims = (/thisAB%A%Nv/)
    
    ts_one_dims(1) = this%nvals
    call h5screate_simple_f(1, ts_one_dims, this%ts_s_id, this%error)

  end subroutine create_AB_h5
  subroutine close_AB_h5(this)
    type(datafile_h5), intent(inout) :: this
    
    call h5dclose_f(this%time_d_id, this%error)
    call h5sclose_f(this%time_s_id, this%error)
    
    call h5gclose_f(this%data_g_id, this%error)
    call h5gclose_f(this%data_g2_id, this%error)
    
    call h5gclose_f(this%info_g_id, this%error)
    
    call h5fclose_f(this%file_id, this%error)
    
    call h5close_f(this%error)
  end subroutine close_AB_h5

  subroutine compute_M(this)
    type(AB), intent(inout) :: this

    integer :: i

    this%MAx=0.d0 ; this%MAy=0.d0 ; this%MBx=0.d0 ; this%MBY=0.d0 ; 

    do i=1,this%A%Nx
       this%MAx = this%MAx + cos(get_x(this%A,i)) * this%A%rho(i)
       this%MAy = this%MAy + sin(get_x(this%A,i)) * this%A%rho(i)
    end do
    do i=1,this%B%Nx
       this%MBx = this%MBx + cos(get_x(this%B,i)) * this%B%rho(i)
       this%MBy = this%MBy + sin(get_x(this%B,i)) * this%B%rho(i)
    end do

    this%MAx = this%MAx * this%A%dx
    this%MAy = this%MAy * this%A%dx
    this%MBx = this%MBx * this%B%dx
    this%MBy = this%MBy * this%B%dx

    this%Mx = (1.d0 - this%gamma) * this%MAx + this%gamma * this%MBx
    this%My = (1.d0 - this%gamma) * this%MAy + this%gamma * this%MBy

  end subroutine compute_M

  subroutine compute_force(this)
    type(AB), intent(inout) :: this

    integer :: i

    do i=1,this%A%Nx
       this%A%force(i) = - this%Mx  * sin(get_x(this%A,i)) + this%My * cos(get_x(this%A,i))
    end do
    do i=1,this%B%Nx
       this%B%force(i) = - this%Mx * sin(get_x(this%B,i)) + this%My * cos(get_x(this%B,i))
    end do

  end subroutine compute_force

  subroutine compute_phys(this, phys, time)
    type(AB), intent(inout) :: this
    double precision, intent(out) :: phys(:)
    double precision, intent(in) :: time

    integer :: i,m
    double precision :: masseA, masseB, energy, intA, intB, intAB, kinA, kinB
    double precision :: fA, fB
    double precision :: dAx, dAy, dBx, dBy, pFA, pFB
    
    masseA = 0.d0 ; masseB = 0.d0 ; kinA = 0.d0 ; kinB = 0.d0
    dAx=0.d0 ; dAy=0.d0 ; dBx=0.d0 ; dBy=0.d0 ; pFA = 0.d0 ; pFB = 0.d0 ;
    do i=1,this%A%Nx
       do m=1,this%A%Nv
          fA = this%A%f(i,m)
          masseA = masseA + fA
          kinA = kinA + fA * get_v(this%A,m)**2
          dAx = dAx - fA*sin(get_x(this%A,i))*get_v(this%A,m)
          dAy = dAy + fA*cos(get_x(this%A,i))*get_v(this%A,m)
          pFA = pFA + fA*get_v(this%A,m)*(-this%Mx*sin(get_x(this%A,i))+this%My*cos(get_x(this%A,i)) )
       end do
    end do

    do i=1,this%B%Nx
       do m=1,this%B%Nv
          fB = this%B%f(i,m)
          masseB = masseB + fB
          kinB = kinB + fB * get_v(this%B,m)**2
          dBx = dBx - fB*sin(get_x(this%B,i))*get_v(this%B,m)
          dBy = dBy + fB*cos(get_x(this%B,i))*get_v(this%B,m)
          pFB = pFB + fB*get_v(this%B,m)*(-this%Mx*sin(get_x(this%B,i))+this%My*cos(get_x(this%B,i)) )
       end do
    end do

    masseA = masseA * this%A%dx * this%A%dv
    masseB = masseB * this%B%dx * this%B%dv
    kinA = kinA * .5d0 * this%A%dx * this%A%dv
    kinB = kinB * .5d0 * this%B%dx * this%B%dv

    dAx = dAx * this%A%dx * this%A%dv
    dAy = dAy * this%A%dx * this%A%dv
    pFA = pFA * this%A%dx * this%A%dv

    dBx = dBx * this%B%dx * this%B%dv
    dBy = dBy * this%B%dx * this%B%dv
    pFB = pFB * this%B%dx * this%B%dv

    intAB = 0.5d0 * ( 1.d0 - this % Mx**2 - this % My**2 )

    energy = (1.d0 - this % gamma) * kinA + this % gamma * kinB + intAB

    phys = (/time, masseA, masseB, energy, intAB, kinA, kinB, this%MAx, this%MAy, this%MBx, this%MBy, this%Mx, this%My , dAx, dAy, dBx, dBy , pFA , pFB /)

  end subroutine compute_phys


end module AB_module
