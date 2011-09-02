module FEL_module
  use Vlasov_module
  implicit none

!  double precision, parameter :: PI = 3.141592653589793115997963468544185161590576171875 !atan(1.d0)*4.d0

  type FEL
     type(grid) :: V
     double precision :: I, phi
     double precision :: Mx, My
     double precision :: delta
     double precision :: sommesin, sommecos
     double precision :: Ax(4), Ay(4), Bx(4), By(4)
  end type FEL
  
contains
  
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
  
  subroutine compute_M(this)
    type(FEL), intent(inout) :: this

    integer :: i

    this%Mx = 0.d0 ; this%My = 0.d0 ;
    do i=1,this%V%Nx
       this%Mx = this%Mx + cos(get_x(this%V,i)) * this%V%rho(i)
       this%My = this%My + sin(get_x(this%V,i)) * this%V%rho(i)
    end do
!    this%Mx = this%Mx - 0.5d0 * (cos(get_x(this%V,1))*this%V%rho(1) + cos(get_x(this%V,this%V%Nx))*this%V%rho(this%V%Nx))
!    this%My = this%My - 0.5d0 * (sin(get_x(this%V,1))*this%V%rho(1) + sin(get_x(this%V,this%V%Nx))*this%V%rho(this%V%Nx))

    this%Mx = this%Mx * this%V%dx
    this%My = this%My * this%V%dx

  end subroutine compute_M

  subroutine compute_force(this)
    type(FEL), intent(inout) :: this
    
    integer :: i

    do i=1,this%V%Nx
       this%V%force(i) = - 2.d0 * sqrt(this%I)*(cos(get_x(this%V,i)-this%phi))
    end do
    
  end subroutine compute_force

  subroutine compute_force_A(this)
    type(FEL), intent(inout) :: this

    integer :: i

    do i=1, this%V%Nx
       this%V%force(i) = -2.d0*(this%Ax(1)*cos(get_x(this%V,i))-this%Ay(1)*sin(get_x(this%V,i)))
    end do
  end subroutine compute_force_A

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

  subroutine correct_A(this)
    type(FEL), intent(inout) :: this

    this%Bx(1) = (this%Ax(1)-this%Ax(2))/this%V%DT
    this%By(1) = (this%Ay(1)-this%Ay(2))/this%V%DT
    
    this%Ax(1) = this%Ax(2) + this%V%DT/12.d0 * (5.d0*this%Bx(1) + 8.d0*this%Bx(2) - 1.d0*this%Bx(3))
    this%Ay(1) = this%Ay(2) + this%V%DT/12.d0 * (5.d0*this%By(1) + 8.d0*this%By(2) - 1.d0*this%By(3))

  end subroutine correct_A
  
  subroutine c_somme(this)
    type(FEL), intent(inout) :: this
    
    integer :: i

    this%sommesin = 0.d0
    this%sommecos = 0.d0
    do i=1,this%V%Nx
       this%sommesin = this%sommesin + sin(get_x(this%V,i)-this%phi)*this%V%rho(i)
       this%sommecos = this%sommecos + cos(get_x(this%V,i)-this%phi)*this%V%rho(i)
    end do
    this%sommesin = this%sommesin * this%V%dx
    this%sommecos = this%sommecos * this%V%dx
    
  end subroutine c_somme

end module FEL_module
