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

!> This module defines a Vlasov setup for the Hamiltonian Mean-Field model.
!!
!! Routines are provided to compute the force field, the thermodynamical observables and the energy
!! distribution function.
!!
!! The HMF model has been introduced in
!! M. Antoni & S. Ruffo, Clustering and relaxation in Hamiltonian long-range dynamics, Phys. Rev. E
!! vol. 52, pp. 2361-2374 (1995).
!!
!! Results of Vlasov simulations for the HMF model are found in
!! P. de Buyl, Numerical resolution of the Vlasov equation for the Hamiltonian Mean-Field model, 
!! Commun. Nonlinear Sci. Numer. Simulat. vol. 15, pp. 2133-2139 (2010).
!!

module HMF_module
  use Vlasov_module
  implicit none

  !> This type holds the data to describe a HMF system.
  type HMF
     !> The grid holding the distribution function.
     type(grid) :: V
     !> The x component of the magnetization.
     double precision :: Mx
     !> The y component of the magnetization.
     double precision :: My
     !> The height of the waterbag.
     double precision :: f0
     !> The energy distribution function.
     double precision, allocatable :: edf(:)
     !> Minimum for the energy DF.
     double precision :: e_min
     !> Maximum for the energy DF.
     double precision :: e_max
     !> Step for the energy DF.
     double precision :: de
     !> If true, an external field is applied.
     logical :: is_ext
     !> Optional coupling parameter.
     double precision :: epsilon
     !> External field.
     double precision :: Hfield
     !> Entropy
     double precision :: entropy
     !> I2
     double precision :: I2
     !> I3
     double precision :: I3
     !> Moments of p
     double precision, allocatable :: pn(:)
  end type HMF
  
contains
  
  !> Initializes a HMF variable.
  !!
  !! @param this A type(HMF) variable.
  !! @param Nx The number of grid points in the x direction.
  !! @param Nv The number of grid points in the v direction.
  !! @param vmax Upper bound of the box in velocity.
  !! @param vmin Lower bound of the box in velocity.
  !! @param Nedf Number of points for the energy DF.
  !! @param model Textual declaration of the model.
  !! @param epsilon Coupling parameter.
  !! @param Hfield External field.
  !! @param n_moments number of moments of p to compute
  subroutine newHMF(this,Nx,Nv,vmax,vmin, Nedf, model, epsilon, Hfield, n_moments)
    type(HMF), intent(out) :: this
    integer, intent(in) :: Nx, Nv
    double precision, intent(in) :: vmax
    double precision, intent(in), optional :: vmin
    integer, intent(in), optional :: Nedf
    character(len=12), optional, intent(in) :: model
    double precision, optional, intent(in) :: epsilon, Hfield
    integer, optional, intent(in) :: n_moments
    double precision :: vmin_final

    if (present(vmin)) then
       vmin_final = vmin
    else
       vmin_final = -vmax
    end if

    call new(this%V,Nx,Nv,PI,vmax,vmin=vmin_final,is_periodic=.true.)

    if (present(Nedf).and.Nedf.gt.1) then
       allocate(this%edf(Nedf))
       this%e_min = 0.d0
       this%e_max = max(vmin_final,vmax)**2*0.5d0 + 1.d0
       this%de = (this%e_max - this%e_min)/(Nedf-1)
    end if

    if (present(model)) then
       if (model.eq.'HMFext') then
          this%is_ext = .true.
          if (present(epsilon)) then
             this%epsilon = epsilon
          else
             this%epsilon = 1.d0
          end if
          if (present(Hfield)) then
             this%Hfield = Hfield
          else
             this%Hfield = 0.d0
          end if
       else
          this%is_ext = .false.
          this%epsilon = 1.d0
          this%Hfield = 0.d0
       end if

    end if

    if (present(n_moments)) then
       if (n_moments .gt. 0) then
          allocate(this%pn(n_moments))
          this%pn = 0.d0
       end if
    end if

  end subroutine newHMF
  
  !> Computes the magnetization.
  !!
  !! @param this A type(HMF) variable.
  subroutine compute_M(this)
    type(HMF), intent(inout) :: this

    integer :: i

    this%Mx = 0.d0 ; this%My = 0.d0 ;
    do i=1,this%V%Nx
       this%Mx = this%Mx + cos(get_x(this%V,i)) * this%V%rho(i)
       this%My = this%My + sin(get_x(this%V,i)) * this%V%rho(i)
    end do

    this%Mx = this%Mx * this%V%dx
    this%My = this%My * this%V%dx

  end subroutine compute_M

  !> Computes the force field.
  !!
  !! @param this A type(HMF) variable.
  subroutine compute_force(this)
    type(HMF), intent(inout) :: this
    
    integer :: i

    call compute_M(this)

    if (this%is_ext) then
       do i=1,this%V%Nx
          this%V%force(i) = cos(get_x(this%V,i))*(this%epsilon*this%My) &
               - sin(get_x(this%V,i)) * (this%Mx*this%epsilon + this%Hfield)
       end do
    else
       do i=1,this%V%Nx
          this%V%force(i) = cos(get_x(this%V,i))*this%My - sin(get_x(this%V,i)) * this%Mx
       end do
    end if

  end subroutine compute_force

  !> Computes the macroscopic observables.
  !!
  !! @param this A type(HMF) variable.
  subroutine compute_phys(this)
    type(HMF), intent(inout) :: this

    integer :: i,m
    double precision :: loopf
    integer :: n, n_moments
    logical :: do_pn

    if (allocated(this%pn)) then
       this%pn = 0.d0
       do_pn = .true.
       n_moments = size(this%pn)
    else
       do_pn = .false.
    end if

    if (this%is_ext) then
       this%V%en_int= this%epsilon * 0.5d0 * (1.d0 - this%Mx**2 - this%My**2) + this%Hfield*(1-this%Mx)
    else
       this%V%en_int= 0.5d0 * (1.d0 - this%Mx**2 - this%My**2)
    end if
    this%V%en_kin = 0.d0
    this%V%momentum = 0.d0
    this%V%masse = 0.d0
    this%I2 = 0.d0 ; this%I3 = 0.d0 ; this%entropy = 0.d0
    do i=1,this%V%Nx
       do m=1,this%V%Nv
          loopf = this%V%f(i,m)
          this%V%en_kin = this%V%en_kin + get_v(this%V,m)**2 * loopf
          this%V%momentum = this%V%momentum + get_v(this%V,m) * loopf
          this%I2 = this%I2 + loopf**2
          this%I3 = this%I3 + loopf**3
          if (do_pn) then
             do n=1,n_moments
                this%pn(n) = this%pn(n) + (get_v(this%V, m)**n) * loopf
             end do
          end if
          if (loopf.gt.0.d0 .and. loopf.lt.this%f0) then
             loopf = loopf/this%f0
             this%entropy = this%entropy + loopf*log(loopf) + (1.d0-loopf)*log(1.d0-loopf)
          end if
       end do
    end do
    this%V%masse = sum(this%V%f(1:this%V%Nx,:))
    this%V%en_kin = this%V%en_kin*0.5d0  * this%V%dx * this%V%dv
    this%V%momentum = this%V%momentum    * this%V%dx * this%V%dv
    this%V%masse    = this%V%masse       * this%V%dx * this%V%dv
    this%I2  = this%I2                   * this%V%dx * this%V%dv
    this%I3  = this%I3                   * this%V%dx * this%V%dv
    this%entropy  = this%entropy         * this%V%dx * this%V%dv
    this%V%energie = this%V%en_int + this%V%en_kin
    if (do_pn) then
       do n=1,n_moments
          this%pn(n) = this%pn(n) * this%V%dx * this%V%dv
       end do
    end if

  end subroutine compute_phys

  !> Computes the energy distribution function.
  !!
  !! @param this A type(HMF) variable.
  subroutine compute_edf(this)
    type(HMF), intent(inout) :: this

    integer :: i, m, n, pos
    double precision :: h

    n = size(this%edf)
    this%edf=0.d0
    do i=1,this%V%Nx
       do m=1,this%V%Nv
          h = get_v(this%V,m)**2*0.5d0 &
               + this%epsilon*(1.d0-(this%Mx*cos(get_x(this%V,i))+this%My*sin(get_x(this%V,i)))) &
               + this%Hfield*(1.d0-cos(get_x(this%V,i)))
          pos = max(min(floor((h-this%e_min)/this%de),size(this%edf)-1)+1, 1)
          this%edf(pos) = this%edf(pos) + this%V%f(i,m)*this%V%dx*this%V%dv/this%de
       end do
    end do
     
  end subroutine compute_edf

end module HMF_module
