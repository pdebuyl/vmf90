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

!> This module defines a Vlasov setup for the Colson-Bonifacio model for the free electron
!! laser.
!!
!! Routines are provided to compute the force field and the thermodynamical observables.
!!
!! The Vlasov equation for this model has been introduced in
!! J. Barré, T. Dauxois, G. De Ninno, D. Fanelli and S. Ruffo, Statistical theory of high-gain
!! free-electron laser saturation, Phys. Rev. E vol. 69, 045501 (2004).
!!
!! Results of Vlasov simulations for this model are found in
!! P. de Buyl, D. Fanelli, R. Bachelard and G. De Ninno, Out-of-equilibrium mean-field dynamics
!! of a model for wave-particle interaction, Phys. Rev. ST Accel. Beams vol. 12, 060704 (2009).
!!

module FEL_module
  use Vlasov_module
  implicit none

  !> This type holds the data to describe a HMF system.
  type FEL
     !> The grid holding the distribution function.
     type(grid) :: V
     !> The intensity of the field.
     double precision :: I
     !> The phase of the field.
     double precision :: phi
     !> The x component of the instantaneous magnetization.
     double precision :: Mx
     !> The y component of the instantaneous magnetization.
     double precision :: My
     !> The detuning.
     double precision :: delta
     !> The height of the waterbag.
     double precision :: f0
     !> The x component of the field. Ax(1) is the present step and the three other components
     !! are previous steps for the predictor-corrector algorithm.
     double precision :: Ax(4)
     !> The y component of the field. Ay(1) is the present step and the three other components
     !! are previous steps for the predictor-corrector algorithm.
     double precision :: Ay(4)
     !> The y component of the field derivative. Bx(1) is the present step and the three
     !! other components are previous steps for the predictor-corrector algorithm.
     double precision :: Bx(4)
     !> The y component of the field derivative. By(1) is the present step and the three
     !! other components are previous steps for the predictor-corrector algorithm.
     double precision :: By(4)
     !> Entropy
     double precision :: entropy
     !> L2 norm
     double precision :: L2
  end type FEL
  
contains
  
  !> Initializes a FEL variable.
  !!
  !! @param this A type(FEL) variable.
  !! @param Nx The number of grid points in the x direction.
  !! @param Nv The number of grid points in the v direction.
  !! @param vmax Upper bound of the box in velocity.
  !! @param vmin Lower bound of the box in velocity.
  !! @param delta Detuning.
  subroutine newFEL(this,Nx,Nv,vmax,vmin, delta)
    type(FEL), intent(out) :: this
    integer, intent(in) :: Nx, Nv
    double precision, intent(in) :: vmax
    double precision, intent(in), optional :: vmin
    double precision, intent(in), optional :: delta
    
    double precision :: vmin_final

    if (present(vmin)) then
       vmin_final = vmin
    else
       vmin_final = -vmax
    end if

    if (present(delta)) then
       this%delta = delta
    else
       this%delta = 0.d0
    end if

    this%Ax = 0. ; this%Ay = 0. ; this%Bx = 0. ; this%By = 0. ;

    call new(this%V,Nx,Nv,PI,vmax,vmin=vmin_final,is_periodic=.true.)

  end subroutine newFEL
  
  !> Computes the magnetization.
  !!
  !! @param this A type(HMF) variable.
  subroutine compute_M(this)
    type(FEL), intent(inout) :: this

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
  subroutine compute_force_A(this)
    type(FEL), intent(inout) :: this

    integer :: i

    do i=1, this%V%Nx
       this%V%force(i) = -2.d0*(this%Ax(1)*cos(get_x(this%V,i))-this%Ay(1)*sin(get_x(this%V,i)))
    end do
  end subroutine compute_force_A

  !> Predictor step for A.
  !!
  !! @param this A type(HMF) variable.
  subroutine compute_A(this)
    type(FEL), intent(inout) :: this

    integer :: i
    double precision :: temp_Bx, temp_By

    do i=2,4
       this%Ax(i) = this%Ax(i-1) ; this%Ay(i) = this%Ay(i-1) ;
       this%Bx(i) = this%Bx(i-1) ; this%By(i) = this%By(i-1) ;
    end do

    temp_Bx = this%Mx ; temp_By = this%My ;

    this%Bx(2) = - this%delta*this%Ay(2) + temp_Bx
    this%By(2) =   this%delta*this%Ax(2) - temp_By

    this%Ax(1) = this%Ax(2)  + this%V%DT/12.d0 * (23.d0*this%Bx(2) - 16.d0*this%Bx(3) + 5.d0*this%Bx(4))
    this%Ay(1) = this%Ay(2)  + this%V%DT/12.d0 * (23.d0*this%By(2) - 16.d0*this%By(3) + 5.d0*this%By(4))

  end subroutine compute_A

  !> Corrector step for A.
  !!
  !! @param this A type(HMF) variable.
  subroutine correct_A(this)
    type(FEL), intent(inout) :: this

    this%Bx(1) = (this%Ax(1)-this%Ax(2))/this%V%DT
    this%By(1) = (this%Ay(1)-this%Ay(2))/this%V%DT
    
    this%Ax(1) = this%Ax(2) + this%V%DT/12.d0 * (5.d0*this%Bx(1) + 8.d0*this%Bx(2) - 1.d0*this%Bx(3))
    this%Ay(1) = this%Ay(2) + this%V%DT/12.d0 * (5.d0*this%By(1) + 8.d0*this%By(2) - 1.d0*this%By(3))

  end subroutine correct_A

  !> Computes the macroscopic observables.
  !!
  !! @param this A type(FEL) variable.
  !! @param phys An array holding the observables.
  !! @param time The real-valued time. Is inserted with the observables in phys.
  subroutine compute_phys(this)
    type(FEL), intent(inout) :: this

    integer :: i,m
    double precision :: loopf

    this%V%masse = sum(this%V%rho)*this%V%dx
    this%V%en_int = 2.d0*(this%Ax(1)*this%My + this%Ay(1)*this%Mx)
    this%V%en_kin = 0.d0
    this%V%momentum = 0.d0
    this%entropy=0.d0
    this%L2=0.d0
    do i=1,this%V%Nx
       do m=1,this%V%Nv
          loopf = this%V%f(i,m)
          this%V%en_kin = this%V%en_kin + get_v(this%V,m)**2 * loopf
          this%V%momentum = this%V%momentum + get_v(this%V,m) * loopf
          this%L2 = this%L2+loopf**2
          loopf=loopf/this%f0
          if (loopf.gt.0.d0) then
             this%entropy = this%entropy + loopf*log(loopf)
          end if
          if (loopf.lt.1.d0) then
             this%entropy = this%entropy + (1.d0-loopf)*log(1.d0-loopf)
          end if
       end do
    end do
    this%V%en_kin   = this%V%en_kin * 0.5d0  * this%V%dx * this%V%dv
    this%V%momentum = this%V%momentum        * this%V%dx * this%V%dv
    this%entropy    = -this%entropy          * this%V%dx * this%V%dv
    this%L2         = this%L2                * this%V%dx * this%V%dv
    this%V%energie  = this%V%en_int + this%V%en_kin

  end subroutine compute_phys

end module FEL_module
