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

!  double precision, parameter :: PI = 3.141592653589793115997963468544185161590576171875 !atan(1.d0)*4.d0

  type HMF
     type(grid) :: V
     double precision :: Mx, My
     double precision :: f0
     double precision, allocatable :: edf(:)
     double precision :: e_min, e_max, de
     logical :: is_ext
     double precision :: epsilon, Hfield
  end type HMF
  
contains
  
  subroutine newHMF(this,Nx,Nv,vmax,vmin, Nedf, model, epsilon, Hfield)
    type(HMF), intent(out) :: this
    integer, intent(in) :: Nx, Nv
    double precision, intent(in) :: vmax
    double precision, intent(in), optional :: vmin
    integer, intent(in), optional :: Nedf
    character(len=12), optional, intent(in) :: model
    double precision, optional, intent(in) :: epsilon, Hfield
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
          this%epsilon = epsilon
          this%Hfield = Hfield
       else
          this%is_ext = .false.
          this%epsilon = 1.d0
          this%Hfield = 0.d0
       end if

    end if

  end subroutine newHMF
  
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

  subroutine compute_force(this)
    type(HMF), intent(inout) :: this
    
    integer :: i

    call compute_M(this)

    if (this%is_ext) then
       do i=1,this%V%Nx
          this%V%force(i) = cos(get_x(this%V,i))*(this%epsilon*this%My) - sin(get_x(this%V,i)) * (this%Mx*this%epsilon + this%Hfield)
       end do
    else
       do i=1,this%V%Nx
          this%V%force(i) = cos(get_x(this%V,i))*this%My - sin(get_x(this%V,i)) * this%Mx
       end do
    end if

  end subroutine compute_force

  subroutine compute_phys(this, phys, time)
    type(HMF), intent(inout) :: this
    double precision, intent(out) :: phys(:)
    double precision, intent(in) :: time

    integer :: i,m
    double precision :: masse, I2, I3, entropy, loopf

    if (this%is_ext) then
       this%V%en_int= this%epsilon * 0.5d0 * (1.d0 - this%Mx**2 - this%My**2) + this%Hfield*(1-this%Mx)
    else
       this%V%en_int= 0.5d0 * (1.d0 - this%Mx**2 - this%My**2)
    end if
    this%V%en_kin = 0.d0
    this%V%momentum = 0.d0
    masse = 0.d0
    I2 = 0.d0 ; I3 = 0.d0 ; entropy = 0.d0
    do i=1,this%V%Nx
       do m=1,this%V%Nv
          loopf = this%V%f(i,m)
          this%V%en_kin = this%V%en_kin + get_v(this%V,m)**2 * loopf
          this%V%momentum = this%V%momentum + get_v(this%V,m) * loopf
          I2 = I2 + loopf**2
          I3 = I3 + loopf**3
          if (loopf.gt.0.d0 .and. loopf.lt.this%f0) then
             loopf = loopf/this%f0
             entropy = entropy + loopf*log(loopf) + (1.d0-loopf)*log(1.d0-loopf)
          end if
       end do
    end do
    masse = sum(this%V%f(1:this%V%Nx,:))
    this%V%en_kin = this%V%en_kin*0.5d0  * this%V%dx * this%V%dv
    this%V%momentum = this%V%momentum    * this%V%dx * this%V%dv
    masse    = masse                     * this%V%dx * this%V%dv
    I2       = I2                        * this%V%dx * this%V%dv
    I3       = I3                        * this%V%dx * this%V%dv
    entropy  = entropy                   * this%V%dx * this%V%dv
    this%V%energie = this%V%en_int + this%V%en_kin

    phys = (/time, masse, this%V%energie, this%V%en_int, this%V%en_kin, this%V%momentum, this%Mx, this%My, I2, I3, entropy/)


  end subroutine compute_phys

  subroutine compute_edf(this)
    type(HMF), intent(inout) :: this

    integer :: i, m, n, pos
    double precision :: h

    n = size(this%edf)
    this%edf=0.d0
    do i=1,this%V%Nx
       do m=1,this%V%Nv
          h = get_v(this%V,m)**2*0.5d0 + this%epsilon*(1.d0-(this%Mx*cos(get_x(this%V,i))+this%My*sin(get_x(this%V,i)))) + this%Hfield*(1.d0-cos(get_x(this%V,i)))
          pos = max(min(floor((h-this%e_min)/this%de),size(this%edf)-1)+1, 1)
          this%edf(pos) = this%edf(pos) + this%V%f(i,m)*this%V%dx*this%V%dv/this%de
       end do
    end do
     
  end subroutine compute_edf

  subroutine write_edf_h5(this, thisHMF, time)
    type(datafile_h5), intent(inout) :: this
    type(HMF), intent(in) :: thisHMF
    integer, intent(in) :: time

    character(len=32) g_name, edf_name
    integer(HID_T) :: g_id, DSP_id, DSET_id
    integer(HSIZE_T) :: dims(1)

    write(g_name, '(a1,i5.5)') 't', time

    dims(1) = size(thisHMF%edf)

    call h5gopen_f(this%data_g_id, g_name, g_id, this%error)
    call h5screate_simple_f(1, dims, DSP_id, this%error)
    call h5dcreate_f(g_id,'edf', H5T_NATIVE_DOUBLE, DSP_id, DSET_id, this%error)

    call h5dwrite_f(DSET_id, H5T_NATIVE_DOUBLE, thisHMF%edf, dims, this%error)

    call h5dclose_f(DSET_id, this%error)
    call h5sclose_f(DSP_id, this%error)

    call h5gclose_f(g_id, this%error)

  end subroutine write_edf_h5

end module HMF_module
