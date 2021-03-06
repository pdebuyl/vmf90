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
  use h5md
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
     double precision, allocatable :: d2(:), copy(:)
     double precision, allocatable :: rho(:), phi(:), force(:)
  end type grid

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
      allocate(this%d2(max(Nx+bonus_gridpoint,Nv)))
      allocate(this%copy(max(Nx+bonus_gridpoint,Nv)))
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


    !> Advances the solution f by spline interpolation in the x-direction.
    !!
    !! @param this A type(grid) variable.
    !! @param h fraction of timestep this%DT to use.
    subroutine advance_x(this, h)
      type(grid), intent(inout) :: this
      double precision, intent(in) :: h

      integer :: i, m, l
      double precision :: x

      ! the last column is set to be equal to the first to allow interpolation beyond xmin+Nx*dx (up to xmax!)
      this%f(this%Nx+1,:) = this%f(1,:)

      l = size(this% f, dim=1)

      do m=1,this%Nv
         this% copy(1:l) = this% f(:,m)
         call spline_periodic(this%copy(1:this%Nx), this%dx, this%d2(1:this%Nx))
         this% d2(this%Nx+1) = this% d2(1)
         do i=1,this%Nx
            x = get_x(this,i)-this%DT*h*get_v(this,m)
            if (x.le.this%xmin .or. x.ge.this%xmax) then
               if (this%is_periodic) then
                  if (x.le.this%xmin) then
                     x = x + this%xmax-this%xmin
                  end if
                  if (x.ge.this%xmax) then
                     x = x - (this%xmax-this%xmin)
                  end if
               else
                  this%f(i,m) = 0.d0
               end if
            end if

            this%f(i,m) = spline_2(this% copy(1:l), this% dx, this% d2(1:l), x - this% xmin)

         end do
      end do

    end subroutine advance_x

    !> Advances the solution f by spline interpolation in the v-direction.
    !!
    !! @param this A type(grid) variable.
    !! @param h fraction of timestep this%DT to use.
    subroutine advance_v(this, h)
      type(grid), intent(inout) :: this
      double precision, intent(in) :: h

      integer :: i, m, l
      double precision :: v

      l = size(this% f, dim=2)

      do i=1,this%Nx
         this% copy(1:l) = this% f(i,:)
         call spline_natural(this% copy(1:l), this% dv, this% d2(1:l))
         do m=1,this%Nv
            v = get_v(this,m)-this%DT*h*this%force(i)
            if (v-this% vmin.le.0.d0 .or. v-this%vmin.ge.this%vmax-this%vmin) then
               this%f(i,m) = 0.d0
            else
               this%f(i,m) = spline_2(this% copy(1:l), this% dv, this% d2(1:l), v - this% vmin)
            end if
         end do
      end do

    end subroutine advance_v

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

    !> Initializes the distribution function f with an self-consistent waterbag distribution.
    !!
    !! The waterbag is delimited in velocity by \f$ \sqrt{ 2 \left( y0 + m_x \cos\theta \right) } \f$.
    !! The optional argument epsilon allows to apply a perturbation in the form \f$ f \propto 1+\epsilon\cos x\f$
    !! to the distribution.
    !! @param this A type(grid) variable.
    !! @param y0 parameter for the self-consistent waterbag distribution.
    !! @param magnetization Self-consistent value of the magnetization.
    !! @param epsilon Amplitude of the cosine perturbation.
    subroutine init_wb_selfc(this, y0, magnetization, epsilon)
      type(grid), intent(inout) :: this
      double precision, intent(in) :: y0, magnetization
      double precision, intent(in), optional :: epsilon

      integer :: i, m
      double precision :: norme, e, p0var

      if (present(epsilon)) then
         e = epsilon
      else
         e = 0.d0
      end if

      do i=1,this%Nx
         do m=1,this%Nv
            if (abs(get_v(this,m)).le. &
                 sqrt(2.d0*(y0+magnetization*cos(get_x(this,i)))) &
                 ) then
               this%f(i,m) = 1.d0 + e*cos(get_x(this,i))
            else
               this%f(i,m) = 0.d0
            end if
         end do
      end do

      norme = sum(this%f(1:this%Nx,:)) * this%dx * this%dv
      this%f = this%f / norme

    end subroutine init_wb_selfc

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
    !! The function is \f$ f(\theta, v) = \sqrt{\frac{\beta}{2\pi}}
    !! e^{-\beta \frac{(v-p_0)^2}{2}} \left(1+\epsilon\cos\theta\right) \f$ where
    !! \f$ \epsilon = 0 \f$ if unspecified. It represents a small pertubation to
    !! test stability properties.
    !!
    !! @param this A type(grid) variable.
    !! @param beta The inverse temperature of the distribution.
    !! @param epsilon The amplitude of the sinusoidal perturbation.
    !! @param p0 Shift of the gaussian in the v-direction.
   subroutine init_gaussian(this, beta, epsilon, p0)
      type(grid), intent(inout) :: this
      double precision, intent(in) :: beta
      double precision, intent(in), optional :: epsilon, p0

      integer :: i,m
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
            this%f(i,m) = sqrt(beta/(2.d0*PI)) * exp(-beta * (get_v(this,m)-p0var)**2 /2.d0) * (1.d0 + e*cos(get_x(this,i)))
         end do
      end do

      norme = sum(this%f(1:this%Nx,:))* this%dx * this%dv
      this%f = this%f / norme
         
    end subroutine init_gaussian

    !> Initializes the distribution function with an homogeneous square Lorentzian distribution.
    !!
    !! The shape of this distribution is
    !! \f$ f(v) = \gamma^3 \pi^{-2} \frac{1}{\left( v^2 + \gamma^2 \right)^2} \f$.
    !!
    !! An optional perturbation epsilon can be applied.
    !! @param this A type(grid) variable.
    !! @param gamma So-called height-parameter.
    !! @param epsilon The amplitude of the sinusoidal perturbation.
   subroutine init_squared_lorentzian(this, gamma, epsilon)
      type(grid), intent(inout) :: this
      double precision, intent(in) :: gamma
      double precision, intent(in), optional :: epsilon

      integer :: i,m
      double precision :: e

      if (present(epsilon)) then
         e = epsilon
      else
         e = 0.d0
      end if

      do i=1,this%Nx
         do m=1,this%Nv
            this%f(i,m) = ( gamma**3/ ( PI**2 * (get_v(this,m)**2+gamma**2)**2 ) ) * &
                 (1.d0 + e*cos(get_x(this,i)))
         end do
      end do

    end subroutine init_squared_lorentzian
    
    subroutine write_grid_info(this, thisgrid, model)
      integer(HID_T), intent(inout) :: this
      type(grid), intent(in) :: thisgrid
      character(len=*), intent(in) :: model

      call h5md_write_par(this, 'model', model)
      call h5md_write_par(this, 'Nx', thisgrid%Nx)
      call h5md_write_par(this, 'Nv', thisgrid%Nv)
      call h5md_write_par(this, 'DT', thisgrid%DT)
      call h5md_write_par(this, 'xmin', thisgrid%xmin)
      call h5md_write_par(this, 'xmax', thisgrid%xmax)
      call h5md_write_par(this, 'vmin', thisgrid%vmin)
      call h5md_write_par(this, 'vmax', thisgrid%vmax)
      call h5md_write_par(this, 'dx', thisgrid%dx)
      call h5md_write_par(this, 'dv', thisgrid%dv)
     
    end subroutine write_grid_info

  !> Adds a fields group in a h5md file
  !! @param file_id the HDF5 ID of the file.
  subroutine create_fields_group(file_id)
    integer(HID_T), intent(inout) :: file_id

    integer(HID_T) :: g_id

    call h5gcreate_f(file_id, 'fields', g_id, h5_error)
    call h5gclose_f(g_id, h5_error)

  end subroutine create_fields_group

end module Vlasov_module
