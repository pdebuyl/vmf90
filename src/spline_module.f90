!> Module for the computation of natural and periodic spline interpolation.
!!
!! The natural spline has zero-valued second derivatives at the boundaries. The periodic splines
!! are assuming a periodic data series. The interpolation step is done with spline_2 for both 
!! kinds of splines.

module spline_module
  implicit none

contains

  !> Solves a tridiagonal equation.
  !!
  !! The problem is in the form \f$ a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i \f$.
  !!
  !! @param a Lower off-diagonal coefficients.
  !! @param b Diagonal coefficients.
  !! @param c Upper off-diagonal coefficients.
  !! @param d Right-hand side vector.
  !! @param x Left hand-side vector. This is the unknown of the problem.
  subroutine tridiag(a,b,c,x,d)
    implicit none
    double precision, dimension(:), intent(inout) :: a,b,c
    double precision, dimension(:), intent(inout) :: d
    double precision, dimension(size(d)), intent(out) :: x

    integer :: i
    
    c(1) = c(1) / b(1)
    do i=2,size(d)-1
       c(i) = c(i)/(b(i)-c(i-1)*a(i))
    end do
    d(1) = d(1)/b(1)
    do i=2,size(d)
       d(i) = (d(i)-d(i-1)*a(i))/(b(i)-c(i-1)*a(i))
    end do
    x(size(d)) = d(size(d))
    do i=size(d)-1,1,-1
       x(i) = d(i) - c(i)*x(i+1)
    end do
  end subroutine tridiag


  !> Computes the second derivative of the data array y.
  !! This subroutine is suited for nonperiodic data, the second derivative is assumed
  !! to be zero at the boundaries.
  !!
  !! @param y Data array.
  !! @param dx Data spacing.
  !! @param y2 The resulting second derivative.
  subroutine spline_natural(y, dx, y2)
    implicit none
    double precision, intent(in) :: y(:), dx
    double precision, intent(out) :: y2(size(y))

    double precision, dimension(size(y)) :: y_prime, a,b,c
    integer :: n
    n = size(y)

    y_prime(2:n-1) = (y(3:n)-y(2:n-1))/dx - (y(2:n-1)-y(1:n-2))/dx
    y_prime(1) = 0.d0
    y_prime(n) = 0.d0

    a(2:n) = dx/6.d0
    c(1:n-1) = dx/6.d0
    b(2:n-1) = 2.d0*dx/3.d0
    b(1) = b(2)
    b(n) = b(n-1)
    call tridiag(a,b,c,y2,y_prime)

  end subroutine spline_natural


  !> Computes the second derivative of the data array y.
  !! This subroutine is suited for periodic data.
  !!
  !! @param y Data array.
  !! @param dx Data spacing.
  !! @param y2 The resulting second derivative.
  subroutine spline_periodic(y,dx,y2)
    implicit none
    double precision, intent(in) :: y(:), dx
    double precision, intent(out) :: y2(size(y))

    double precision, dimension(size(y)) :: y_prime, u, z, a, b, c
    double precision :: factor, alpha, beta, gamma
    integer :: n
    n = size(y)

    y_prime(2:n-1) = (y(3:n)-y(2:n-1))/dx - (y(2:n-1)-y(1:n-2))/dx
    y_prime(1) = (y(2)-y(1))/dx - (y(1)-y(n))/dx
    y_prime(n) = (y(1)-y(n))/dx - (y(n)-y(n-1))/dx

    alpha = dx/6.d0
    beta = alpha
    gamma = -2.d0*dx/3.d0
    a(1:n) = dx/6.d0
    c(1:n) = dx/6.d0
    b(1:n) = 2.d0*dx/3.d0
    b(1) = b(1) - gamma
    b(n) = b(n) - alpha*beta/gamma
    call tridiag(a,b,c,y2,y_prime)
    u = 0.d0 ; u(1) = gamma ; u(n) = alpha
    a(1:n) = dx/6.d0
    c(1:n) = dx/6.d0
    b(1:n) = 2.d0*dx/3.d0
    b(1) = b(1) - gamma
    b(n) = b(n) - alpha*beta/gamma
    call tridiag(a,b,c,z,u)
    factor = ( y2(1)+ y2(n)*beta/gamma ) / (1.d0 + z(1) + z(n)*beta/gamma)
    y2 = y2 - factor * z

  end subroutine spline_periodic

  !> Given a data array, its step size and the second derivative data array, returns
  !! the cubic spline at point x.
  !! @param y Data points.
  !! @param dx The grid stepsize.
  !! @param y2 The second derivative of the data.
  !! @param x The point at which the evaluation is made.
  !! @return The interpolated value at x.
  function spline_2(y,dx,y2,x)
    implicit none
    double precision :: spline_2
    double precision, intent(in) :: y(:),dx,y2(:), x

    integer :: i
    double precision :: my_a, my_b, my_c, my_d

    i = floor(x/dx)+1
    my_b = x/dx-dble(i-1)
    my_a = 1.d0-my_b
    my_c = (my_a**3-my_a)*dx**2/6.d0
    my_d = (my_b**3-my_b)*dx**2/6.d0

    spline_2 = y(i)*my_a + y(i+1)*my_b + y2(i)*my_c + y2(i+1)*my_d

  end function spline_2

end module spline_module
