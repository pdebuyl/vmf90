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

!> This module provides facilities to describe a (x,p) phase space grid and to perform
!! operations on it to solve the Vlasov equation with the second order Cheng-Knorr splitting,
!! for now only with cubic spline interpolation.
!!
!! Utility routines are provided for setting up the system and to obtain grid information such
!! as position and velocity, on the basis of the indices, or computation of the position- and 
!! velocity marginals.
!!
!! The storage of the one-particle phase space distribution function is found in the array f.
!! Storage of a copy of f is necessary for the algorithm, as it is written at the moment.
!! 
!! The semi-Lagrangian methodology that is used in this code is found in 
!! E. Sonnendrucker, J. Roche, P. Bertrand and A. Ghizzo. The Semi-Lagrangian Method for the 
!! Numerical Resolution of the Vlasov Equation, J. Comp. Phys. vol. 149, pp. 201-220 (1999).
!!

module Vlasov_module
  use spline_module
  use HDF5
  implicit none

  !> Value of Pi computed via N[Pi,35] in Mathematica
  double precision, parameter :: PI = 3.1415926535897932384626433832795029d0

  !> The type grid allows to describe (position and velocity coordinates) and store (array f) a 
  !! numerical distribution function in the one-particle phase space.
  type grid
     integer :: Nx, Nv
     double precision :: xmin, xmax, vmin, vmax
     double precision :: dx, dv
     double precision :: DT
     double precision :: en_int, en_kin, energie, masse, momentum
     logical :: is_periodic

     double precision, allocatable :: f(:,:)
     double precision, allocatable :: f2(:,:), g(:,:)
     double precision, allocatable :: rho(:), phi(:), force(:)
  end type grid

  !> The type datafile_h5 is used to store HDF5 variables for a Vlasov simulation.
  type datafile_h5
     character(len=64) :: filename
     character(len=4) :: info_g_name="info"
     character(len=4) :: data_g_name="data"
     character(len=5) :: data_g2_name="data2"
     character(len=1) :: f_d_name = "f"
     character(len=3) :: rho_d_name = "rho"
     character(len=3) :: phi_d_name = "phi"
     character(len=4) :: time_d_name = "time"
     integer(HID_T) :: file_id
     integer(HID_T) :: info_g_id, data_g_id, data_g2_id
     integer(HID_T) :: f_s_id, rho_s_id, phi_s_id, time_s_id, ts_s_id
     integer(HID_T) :: f_d_id, rho_d_id, phi_d_id, time_d_id
     integer :: f_rank=2, rho_rank=1, phi_rank=1, time_rank=2
     integer(HSIZE_T) :: f_dims(2), rho_dims(1), phi_dims(1), time_dims(2)
     integer :: ntime
     integer :: error, nvals
  end type datafile_h5

  contains

    !> Creates a new type(grid) variable with given dimensions.
    !! 
    !! The routine new allocates data arrays and sets all the descriptive variables of the type(grid)
    !! variable.
    !! @param this A type(grid) variable.
    !! @param Nx The number of grid points in the x-direction.
    !! @param Nv The number of grid points in the v-direction.
    !! @param xmax The maximum x coordinate.
    !! @param xmin The minimum x coordinate.
    !! @param vmax The maximum v coordinate.
    !! @param vmin The minimum v coordinate.
    !! @param is_periodic Sets the type(grid) variable to be periodic or not.
    subroutine new(this,Nx,Nv,xmax,vmax,xmin,vmin,is_periodic)
      type(grid), intent(out) :: this
      integer, intent(in) :: Nx, Nv
      double precision, intent(in) :: xmax, vmax
      double precision, intent(in), optional :: xmin, vmin
      logical, intent(in), optional :: is_periodic
      integer :: i,m, Nmax, bonus_gridpoint

      this%Nx = Nx
      this%Nv = Nv

      if (present(is_periodic)) then
         this%is_periodic=is_periodic
      else
         this%is_periodic=.true.
      end if

      if (this%is_periodic) then
         bonus_gridpoint=1
      else
         bonus_gridpoint=0
      end if
      
      allocate(this%f(Nx+bonus_gridpoint,Nv))
      allocate(this%f2(Nx+bonus_gridpoint,Nv))
      allocate(this%g(Nx+bonus_gridpoint,Nv))
      allocate(this%rho(Nx))
      allocate(this%force(Nx))
      allocate(this%phi(Nv))

      this%xmax=xmax
      this%vmax=vmax
      if (present(xmin)) then
         this%xmin=xmin
      else
         this%xmin=-xmax
      end if
      if (present(vmin)) then
         this%vmin=vmin
      else
         this%vmin=-vmax
      end if

      if (this%is_periodic) then
         this%dx = (this%xmax-this%xmin)/dble(Nx)
         this%dv = (this%vmax-this%vmin)/dble(Nv-1)
      else
         this%dx = (this%xmax-this%xmin)/dble(Nx-1)
         this%dv = (this%vmax-this%vmin)/dble(Nv-1)
      end if      

    end subroutine new

    !> Sets the timestep for the grid.
    !! @param this A type(grid) variable.
    !! @param DT The timestep.
    subroutine set_DT(this,DT)
      type(grid), intent(inout) :: this
      double precision, intent(in) :: DT
      this%DT=DT
    end subroutine set_DT
    !> Returns the time step for the type(grid) variable.
    !! @param this A type(grid) variable.
    !! @return The timestep.
    double precision function get_DT(this)
      type(grid), intent(in) :: this
      get_DT = this%DT
    end function get_DT
    !> Returns the maximum x coordinate for the type(grid) variable.
    !! @param this A type(grid) variable.
    !! @return xmax.
    double precision function get_xmax(this)
      type(grid), intent(in) :: this
      get_xmax = this%xmax
    end function get_xmax
    !> Returns the minimum x coordinate for the type(grid) variable.
    !! @param this A type(grid) variable.
    !! @return xmin.
    double precision function get_xmin(this)
      type(grid), intent(in) :: this
      get_xmin = this%xmin
    end function get_xmin
    !> Returns the maximum v coordinate for the type(grid) variable.
    !! @param this A type(grid) variable.
    !! @return vmax.
    double precision function get_vmax(this)
      type(grid), intent(in) :: this
      get_vmax = this%vmax
    end function get_vmax
    !> Returns the minimum v coordinate for the type(grid) variable.
    !! @param this A type(grid) variable.
    !! @return vmin.
    double precision function get_vmin(this)
      type(grid), intent(in) :: this
      get_vmin = this%vmin
    end function get_vmin
    !> Returns the x coordinate for the type(grid) variable for the column index ix.
    !! @param this A type(grid) variable.
    !! @param ix A column index for the position.
    !! @return x.
    double precision function get_x(this,ix)
      type(grid), intent(in) :: this
      integer, intent(in) :: ix
      get_x = this%xmin + (ix-1)*this%dx
    end function get_x
    !> Returns the v coordinate for the type(grid) variable for the row index mv.
    !! @param this A type(grid) variable.
    !! @param mv A row index for the velocity.
    !! @return v.
    double precision function get_v(this,mv)
      type(grid), intent(in) :: this
      integer, intent(in) :: mv
      get_v = this%vmin + (mv-1)*this%dv
    end function get_v


    !> Computes the second derivatives for spline interpolation in the x-direction.
    !!
    !! @param this A type(grid) variable.
    subroutine spline_x(this)
      type(grid), intent(inout) :: this

      integer m
      
      ! the last column is set to be equal to the first to allow interpolation beyond xmin+Nx*dx (up to xmax!)
      this%f(this%Nx+1,:) = this%f(1,:)

      do m=1,this%Nv
         call spline_periodic(this%f(1:this%Nx,m), this%dx, this%f2(1:this%Nx,m))
         this% f2(this%Nx+1,m) = this% f2(1,m)
      end do

    end subroutine spline_x

    !> Returns the interpolated value on the m-th row of the grid.
    !!
    !! If the grid is periodic, takes it into account. Else, returns 0 outside of the box.
    !!
    !! @param this A type(grid) variable.
    !! @param x_in The point at which to interpolate.
    !! @param m The velocity row index.
    !! @return Interpolated value of f.
    double precision function splint_x(this,x_in,m)
      type(grid), intent(in) :: this
      double precision, intent(in) :: x_in
      integer, intent(in) :: m
      double precision :: x

      x = x_in

      if (x.le.this%xmin .or. x.ge.this%xmax) then
         if (this%is_periodic) then
            if (x.le.this%xmin) then
               x = x + this%xmax-this%xmin
            end if
            if (x.ge.this%xmax) then 
               x = x - (this%xmax-this%xmin)
            end if
         else
            splint_x = 0.d0
            return
         end if
      end if

      splint_x = spline_2(this% f(:,m), this% dx, this% f2(:,m), x - this% xmin)
      
    end function splint_x
    
    !> Computes the second derivatives for spline interpolation in the v-direction.
    !!
    !! @param this A type(grid) variable.
    subroutine spline_v(this)
      type(grid), intent(inout) :: this

      integer i
      
      do i=1,this%Nx
         call spline_natural(this% f(i,:), this% dv, this% f2(i,:))
      end do

    end subroutine spline_v

    !> Returns the interpolated value on the i-th column of the grid.
    !!
    !! Returns 0 outside of the box.
    !!
    !! @param this A type(grid) variable.
    !! @param v The point at which to interpolate.
    !! @param i The position column index.
    !! @return Interpolated value of f.
    double precision function splint_v(this,v,i)
      type(grid), intent(in) :: this
      double precision, intent(in) :: v
      integer, intent(in) :: i

      if (v.lt.this%vmin .or. v.ge.this%vmax) then
         splint_v = 0.d0
         return
      end if

      splint_v = spline_2(this% f(i,:), this% dv, this% f2(i,:), v - this% vmin)
      
    end function splint_v


    !> Performs a half timestep advection in the x-direction.
    !!
    !! @param this A type(grid) variable.
    subroutine advection_x_demi(this)
      type(grid), intent(inout) :: this
      
      integer :: i,m

      do i=1,this%Nx
         do m=1,this%Nv
            this%g(i,m) = splint_x(this,get_x(this,i)-this%DT*get_v(this,m)*0.5d0,m)
         end do
      end do

    end subroutine advection_x_demi

    !> Performs a full timestep advection in the x-direction.
    !!
    !! @param this A type(grid) variable.
    subroutine advection_x(this)
      type(grid), intent(inout) :: this
      
      integer :: i,m

      do i=1,this%Nx
         do m=1,this%Nv
            this%g(i,m) = splint_x(this,get_x(this,i)-this%DT*get_v(this,m),m)
         end do
      end do

    end subroutine advection_x

    !> Performs a full timestep advection in the v-direction.
    !!
    !! @param this A type(grid) variable.
    subroutine advection_v(this)
      type(grid), intent(inout) :: this

      integer :: i,m

      do i=1,this%Nx
         do m=1,this%Nv
            this%g(i,m) = splint_v(this,get_v(this,m)-this%DT*this%force(i),i)
         end do
      end do
         
    end subroutine advection_v
    
    !> Writes the distribution function f in a file of name "xvf.iiiii" in the directory
    !! "images" that should exist beforehand.
    !!
    !! The optional "set_name" argument allows to postfix the filename with one character.
    !! @param this A type(grid) variable.
    !! @param itime The integer time on which the filename is based.
    !! @param unit The Fortran file unit to use.
    !! @param set_name An optional one character postfix.
    subroutine xvf_write(this,itime,unit,set_name)
      type(grid), intent(in) :: this
      integer, intent(in) :: itime,unit
      character(len=1), intent(in), optional :: set_name

      integer :: i,m
      character(len=20) :: filename

      if (present(set_name)) then
         write(filename,'(a11,a1,a1,i5.5)') 'images/xvf.',set_name,'.', itime
      else
         write(filename,'(a11,i5.5)') 'images/xvf.', itime
      end if
      

      open(unit, file=filename)

      do i=1,this%Nx
         do m=1,this%Nv
            write(unit, '(3e20.10e3)') get_x(this,i), get_v(this,m), this%f(i,m)
         end do
         write(unit, *) ''
      end do
      close(unit)
    end subroutine xvf_write

    !> Computes the position-wise marginal distribution.
    !!
    !! \f$ \rho_i = \int_{-v_{max}}^{v_{max}} f(x,v) dv \f$
    !! @param this A type(grid) variable.
    subroutine compute_rho(this)
      type(grid), intent(inout) :: this

      integer i

      do i=1,this%Nx
         this%rho(i) = sum(this%f(i,:))
         this%rho(i) = this%rho(i) * this%dv
      end do
    
    end subroutine compute_rho

    !> Computes the velocity-wise marginal distribution.
    !!
    !! \f$ \phi_i = \int_{-x_{max}}^{x_{max}} f(x,v) dx \f$
    !! @param this A type(grid) variable.
    subroutine compute_phi(this)
      type(grid), intent(inout) :: this
      
      integer m

      do m=1,this%Nv
         this%phi(m) = sum(this%f(1:this%Nx,m))
         this%phi(m) = this%phi(m) * this%dx
      end do

    end subroutine compute_phi

    !> Initializes the distribution function f with a waterbag distribution.
    !!
    !! The waterbag is delimited by [-width:width] in position and [-bag:bag] in velocity.
    !! The optional argument epsilon allows to apply a perturbation in the form \f$ f \propto 1+\epsilon\cos x\f$
    !! to the distribution. The parameter p0 allows to shift the c.o.m. velocity by p0.
    !! @param this A type(grid) variable.
    !! @param width Half-width of the waterbag distribution.
    !! @param bag Half-width in velocity of the waterbag distribution.
    !! @param epsilon Amplitude of the cosine perturbation.
    !! @param p0 Shift of the center of mass velocity.
    subroutine init_carre(this, width, bag, epsilon, p0)
      type(grid), intent(inout) :: this
      double precision, intent(in) :: width, bag
      double precision, intent(in), optional :: epsilon, p0
      
      integer :: i, m
      double precision :: norme, e, p0var
      
      if (present(epsilon)) then
         e = epsilon
      else
         e = 0.d0
      end if

      if (present(p0)) then
         p0var = p0
      else
         p0var = 0.d0
      end if

      do i=1,this%Nx
         do m=1,this%Nv
            if ( (abs(get_x(this,i)).le.width) .and. (abs(get_v(this,m)-p0var).le.bag) ) then
               this%f(i,m) = 1.d0 + e*cos(get_x(this,i))
            else
               this%f(i,m) = 0.d0
            end if
         end do
      end do
      
      norme = sum(this%f(1:this%Nx,:)) * this%dx * this%dv
      this%f = this%f / norme

    end subroutine init_carre

    subroutine init_two_streams(this, v_min, v_max, width)
      type(grid), intent(inout) :: this
      double precision, intent(in) :: v_min, v_max
      double precision, intent(in), optional :: width

      integer :: i,m
      double precision :: norme, v, w

      if (present(width)) then
         w = width
      else
         w = atan(1.d0)/4.d0
      end if

      do i=1,this%Nx
         do m=1,this%Nv
            v = abs(get_v(this,m))
            if (v.ge.v_min .and. v.le.v_max .and. abs(get_x(this,i)).le.w) then
               this%f(i,m) = 1.d0
            else
               this%f(i,m) = 0.d0
            end if
         end do
      end do
      norme = sum(this%f(1:this%Nx,:)) * this%dx * this%dv
      this%f = this%f / norme

    end subroutine init_two_streams

    !> Initializes the distribution function with an homogeneous gaussian distribution.
    !!
    !! An optional perturbation epsilon can be applied.
    !! @param this A type(grid) variable.
    !! @param beta The inverse temperature of the distribution.
    !! @param epsilon The amplitude of the sinusoidal perturbation.
   subroutine init_gaussian(this, beta, epsilon)
      type(grid), intent(inout) :: this
      double precision, intent(in) :: beta
      double precision, intent(in), optional :: epsilon

      integer :: i,m
      double precision :: norme, e

      if (present(epsilon)) then
         e = epsilon
      else
         e = 0.d0
      end if

      do i=1,this%Nx
         do m=1,this%Nv
            this%f(i,m) = sqrt(beta/(2.d0*PI)) * exp(-beta * get_v(this,m)**2 /2.d0) * (1.d0 + e*sin(get_x(this,i)))
         end do
      end do

      norme = sum(this%f(1:this%Nx,:))* this%dx * this%dv
      this%f = this%f / norme
         
    end subroutine init_gaussian
    
    subroutine create_h5(this, thisgrid, filename, ntime, nvals)
      type(datafile_h5), intent(out) :: this
      type(grid), intent(in) :: thisgrid
      character(len=*), intent(in) :: filename
      integer, intent(in) :: ntime, nvals

      integer(HSIZE_T) :: ts_one_dims(1)

      call h5open_f(this%error)
      this%filename = filename
      call h5fcreate_f(this%filename, H5F_ACC_TRUNC_F, this%file_id, this%error)

      call h5gcreate_f(this%file_id, this%info_g_name, this%info_g_id, this%error)
      call h5gcreate_f(this%file_id, this%data_g_name, this%data_g_id, this%error)

      this%ntime = ntime
      this%nvals = nvals
      this%time_dims = (/ntime+1,nvals/)
      call h5screate_simple_f(this%time_rank, this%time_dims, this%time_s_id, this%error)
      call h5dcreate_f(this%file_id, this%time_d_name, H5T_NATIVE_DOUBLE, this%time_s_id, this%time_d_id, this%error)


      ts_one_dims(1) = this%nvals
      call h5screate_simple_f(1, ts_one_dims, this%ts_s_id, this%error)

    end subroutine create_h5

    subroutine write_time_slice_h5(this, time, vals)
      type(datafile_h5), intent(inout) :: this
      integer, intent(in) :: time
      double precision, intent(in) :: vals(:)

      integer :: l
      integer(HSIZE_T) :: ts_dims(2), offset(2)

      if (time.gt.this%ntime) then
         stop 'time of the time table exceeded in write_time_slice_h5'
      end if
      l = size(vals)
      if (l.ne.this%nvals) then
         stop 'wrong vals size in write_time_slice_h5'
      end if

      ts_dims = (/1,l/)
      offset = (/time, 0/)

      call h5sselect_hyperslab_f(this%time_s_id, H5S_SELECT_SET_F, offset, ts_dims, this%error)
      call h5dwrite_f(this%time_d_id, H5T_NATIVE_DOUBLE, vals, ts_dims, this%error, this%ts_s_id, this%time_s_id)

    end subroutine write_time_slice_h5

    subroutine close_h5(this)
      type(datafile_h5), intent(inout) :: this

      call h5dclose_f(this%time_d_id, this%error)
      call h5sclose_f(this%time_s_id, this%error)

      call h5gclose_f(this%data_g_id, this%error)

      call h5gclose_f(this%info_g_id, this%error)
      
      call h5fclose_f(this%file_id, this%error)

      call h5close_f(this%error)
    end subroutine close_h5

    subroutine write_data_group_h5(this, thisgrid, time, group2)
      type(datafile_h5), intent(inout) :: this
      type(grid), intent(in) :: thisgrid
      integer, intent(in) :: time
      logical, intent(in), optional :: group2

      character(len=32) g_name, f_name, rho_name, phi_name
      integer(HID_T) :: g_id

      write(g_name, '(a1,i5.5)') 't', time

      if (present(group2)) then
         if (group2) then
            call h5gcreate_f(this%data_g2_id, g_name, g_id, this%error)
            else
               stop 'you should not set group2 to something else than true in write_data_group_h5'
         end if
      else
         call h5gcreate_f(this%data_g_id, g_name, g_id, this%error)
      end if

      this%f_dims = (/thisgrid%Nx,thisgrid%Nv/)
      this%rho_dims = (/thisgrid%Nx/)
      this%phi_dims = (/thisgrid%Nv/)

      call h5screate_simple_f(this%f_rank, this%f_dims, this%f_s_id, this%error)
      call h5screate_simple_f(this%rho_rank, this%rho_dims, this%rho_s_id, this%error)
      call h5screate_simple_f(this%phi_rank, this%phi_dims, this%phi_s_id, this%error)

      write(f_name, '(a6,a2)') g_name, '/f'
      write(rho_name, '(a6,a4)') g_name, '/rho'
      write(phi_name, '(a6,a4)') g_name, '/phi'
      call h5dcreate_f(g_id, this%f_d_name, H5T_NATIVE_DOUBLE, this%f_s_id, this%f_d_id, this%error)
      call h5dcreate_f(g_id, this%rho_d_name, H5T_NATIVE_DOUBLE, this%rho_s_id, this%rho_d_id, this%error)
      call h5dcreate_f(g_id, this%phi_d_name, H5T_NATIVE_DOUBLE, this%phi_s_id, this%phi_d_id, this%error)

      call h5dwrite_f(this%f_d_id, H5T_NATIVE_DOUBLE, thisgrid%f(1:thisgrid%Nx,:), this%f_dims, this%error)
      call h5dwrite_f(this%rho_d_id, H5T_NATIVE_DOUBLE, thisgrid%rho, this%rho_dims, this%error)
      call h5dwrite_f(this%phi_d_id, H5T_NATIVE_DOUBLE, thisgrid%phi, this%phi_dims, this%error)

      call h5dclose_f(this%f_d_id, this%error)
      call h5sclose_f(this%f_s_id, this%error)
      call h5dclose_f(this%rho_d_id, this%error)
      call h5sclose_f(this%rho_s_id, this%error)
      call h5dclose_f(this%phi_d_id, this%error)
      call h5sclose_f(this%phi_s_id, this%error)


      call h5gclose_f(g_id, this%error)

    end subroutine write_data_group_h5

    subroutine write_grid_info(this, thisgrid, model)
      type(datafile_h5), intent(inout) :: this
      type(grid), intent(in) :: thisgrid
      character(len=*), intent(in) :: model
      
      integer(HID_T) :: s_id, d_id
      integer(HSIZE_T) :: dims(1)
      
      call write_info_string_h5(this, 'model', 'FEL')
      call write_info_integer_h5(this, 'Nx', thisgrid%Nx)
      call write_info_integer_h5(this, 'Nv', thisgrid%Nv)
      call write_info_double_h5(this, 'DT', thisgrid%DT)
      call write_info_double_h5(this, 'xmin', thisgrid%xmin)
      call write_info_double_h5(this, 'xmax', thisgrid%xmax)
      call write_info_double_h5(this, 'vmin', thisgrid%vmin)
      call write_info_double_h5(this, 'vmax', thisgrid%vmax)
      call write_info_double_h5(this, 'dx', thisgrid%dx)
      call write_info_double_h5(this, 'dv', thisgrid%dv)
      call write_info_integer_h5(this, 'ntime', this%ntime)
     
    end subroutine write_grid_info

    subroutine write_info_string_h5(this, info_name, info_element)
      type(datafile_h5), intent(inout) :: this
      character(len=*), intent(in) :: info_name
      character(len=*), intent(in) :: info_element
      
      character(len=32) :: string, formatstring
      character(len=64) :: datastring
      integer(HID_T) :: DSP_id, DSET_id, STR_T
      integer(HSIZE_T) :: dims(1)
      integer(SIZE_T) :: l
      dims(1) = 1
      
      write(formatstring, '(a5,i2.2,a1)') "(a5,a", len(info_name), ")"
      write(string, formatstring) 'info/', info_name
      
      l = len(trim(info_element))+1
      call h5tcopy_f(H5T_NATIVE_CHARACTER, STR_T, this%error)
      call h5tset_size_f(STR_T, l, this%error)

      call h5screate_f(H5S_SCALAR_F, DSP_id, this%error)
      call h5dcreate_f(this%file_id, string, STR_T, DSP_id, DSET_id, this%error)
      call h5dwrite_f(DSET_id, STR_T, info_element, dims, this%error)
      call h5dclose_f(DSET_id, this%error)
      call h5sclose_f(DSP_id, this%error)

    end subroutine write_info_string_h5

    subroutine write_info_string_array_h5(this, info_name, info_element)
      type(datafile_h5), intent(inout) :: this
      character(len=*), intent(in) :: info_name
      character(len=*), intent(in) :: info_element(:)
      
      character(len=32) :: string, formatstring
      character(len=64) :: datastring
      integer(HID_T) :: DSP_id, DSET_id, STR_T
      integer(HSIZE_T) :: dims(1)
      integer(SIZE_T) :: l
      dims(1) = size(info_element)
      
      write(formatstring, '(a5,i2.2,a1)') "(a5,a", len(info_name), ")"
      write(string, formatstring) 'info/', info_name

      l=len(info_element(1))
      call h5tcopy_f(H5T_NATIVE_CHARACTER, STR_T, this%error)
      call h5tset_size_f(STR_T, l, this%error)

      call h5screate_simple_f(1, dims, DSP_id, this%error)
      call h5dcreate_f(this%file_id, string, STR_T, DSP_id, DSET_id, this%error)
      call h5dwrite_f(DSET_id, STR_T, info_element, dims, this%error)
      call h5dclose_f(DSET_id, this%error)
      call h5sclose_f(DSP_id, this%error)

    end subroutine write_info_string_array_h5


    subroutine write_info_integer_h5(this, info_name, info_element)
      type(datafile_h5), intent(inout) :: this
      character(len=*), intent(in) :: info_name
      integer, intent(in) :: info_element
      
      character(len=32) :: string, formatstring
      integer(HID_T) :: DSP_id, DSET_id
      integer(HSIZE_T) :: dims(1)
      dims(1) = 1
      
      write(formatstring, '(a5,i2.2,a1)') "(a5,a", len(info_name), ")"
      write(string, formatstring) 'info/', info_name
      
      call h5screate_f(H5S_SCALAR_F, DSP_id, this%error)
      call h5dcreate_f(this%file_id, string, H5T_NATIVE_INTEGER, DSP_id, DSET_id, this%error)
      call h5dwrite_f(DSET_id, H5T_NATIVE_INTEGER, info_element, dims, this%error)
      call h5dclose_f(DSET_id, this%error)
      call h5sclose_f(DSP_id, this%error)

    end subroutine write_info_integer_h5


    subroutine write_info_double_h5(this, info_name, info_element)
      type(datafile_h5), intent(inout) :: this
      character(len=*), intent(in) :: info_name
      double precision, intent(in) :: info_element
      
      character(len=32) :: string, formatstring
      integer(HID_T) :: DSP_id, DSET_id
      integer(HSIZE_T) :: dims(1)
      dims(1) = 1
      
      write(formatstring, '(a5,i2.2,a1)') "(a", len(info_name), ")"
      write(string, formatstring) info_name
      
      call h5screate_f(H5S_SCALAR_F, DSP_id, this%error)
      call h5dcreate_f(this%info_g_id, string, H5T_NATIVE_DOUBLE, DSP_id, DSET_id, this%error)
      call h5dwrite_f(DSET_id, H5T_NATIVE_DOUBLE, info_element, dims, this%error)
      call h5dclose_f(DSET_id, this%error)
      call h5sclose_f(DSP_id, this%error)

    end subroutine write_info_double_h5


    subroutine read_time_slice_h5(this, time, vals)
      type(datafile_h5), intent(inout) :: this
      integer, intent(in) :: time
      double precision, intent(inout) :: vals(:)

      integer :: l
      integer(HSIZE_T) :: ts_dims(2), offset(2)

      if (time.gt.this%ntime) then
         stop 'time of the time table exceeded in write_time_slice_h5'
      end if
      l = size(vals)
      if (l.ne.this%nvals) then
         stop 'wrong vals size in write_time_slice_h5'
      end if

      ts_dims = (/1,l/)
      offset = (/time, 0/)

      call h5sselect_hyperslab_f(this%time_s_id, H5S_SELECT_SET_F, offset, ts_dims, this%error)
      call h5dread_f(this%time_d_id, H5T_NATIVE_DOUBLE, vals, ts_dims, this%error, this%ts_s_id, this%time_s_id)

    end subroutine read_time_slice_h5

    subroutine load_data_from_h5(thisgrid, f_name, d_name)
      type(grid), intent(inout) :: thisgrid
      character(len=*), intent(in) :: f_name, d_name

      integer(HID_T) :: f_id, d_id, s_id
      integer :: ndims
      integer(HSIZE_T) :: dims(2), dumb(2)
      integer :: err

      call h5fopen_f(f_name, H5F_ACC_RDONLY_F, f_id, err)
      call h5dopen_f(f_id, d_name, d_id, err)
      call h5dget_space_f(d_id, s_id, err)
      
      call h5sget_simple_extent_ndims_f(s_id, ndims, err)
      if (ndims.ne.2) stop 'wrong number of dims in load_data_from_h5'
      call h5sget_simple_extent_dims_f(s_id,dims, dumb, err)
      if ( (dims(1).ne.thisgrid%Nx) .or. (dims(2).ne.thisgrid%Nv) ) stop 'mismatch in the actual dimensions and in the stored dimensions'

      call h5dread_f(d_id, H5T_NATIVE_DOUBLE, thisgrid%f(1:thisgrid%Nx,1:thisgrid%Nv), dims, err)

      call h5sclose_f(s_id, err)
      call h5dclose_f(d_id, err)
      call h5fclose_f(f_id, err)

    end subroutine load_data_from_h5


end module Vlasov_module
