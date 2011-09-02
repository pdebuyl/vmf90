program runFEL
  use FEL_module
  use ParseText
  implicit none
  
  type(FEL) :: F
  
  type(PTo) :: HCF
  type(datafile_h5) :: h5fel

  integer :: Nx, Nv
  double precision :: width, bag, vmin, vmax, DT, x1, x2, x3, Mx, My, p0
  integer :: i,m,t,t_top, n_images, t_images
  integer :: n_steps, n_top
  
  double precision :: masse, entropy, loopf, fzero, temp, L2
  integer, parameter :: time_number = 12
  double precision :: vals(time_number)
  character(len=12) :: time_names(time_number), IC
  
! declaration for the routines allowing the inversion of m = m(width)
  external zbrent, ftbz
  double precision zbrent

  call PTparse(HCF,'FEL_in',7)
  Nx = PTread_i(HCF,'Nx')
  Nv = PTread_i(HCF,'Nv')
  vmax = PTread_d(HCF,'vmax')
  vmin = PTread_d(HCF,'vmin')
  n_steps = PTread_i(HCF,'n_steps')  !100
  n_top = PTread_i(HCF,'n_top')      !100
  DT = PTread_d(HCF,'DT')
  IC = PTread_s(HCF, 'IC')
  n_images = PTread_i(HCF, 'n_images')
  t_images = 1

  call newFEL(F,Nx,Nv,vmax,vmin, delta=PTread_d(HCF,'delta'))
  F%V%DT = DT
  F%phi = 0.d0
  F%I = 0.0001d0
  

  write(*,*) 'launching FEL with ', Nx, 'x', Nv,  'points'
  F%Ax = 0.0001d0;
  F%Ay = 0.d0;

  if (IC.eq.'waterbag') then
     width = PTread_d(HCF, 'width')
     bag = PTread_d(HCF, 'bag')
     call init_carre(F%V, width, bag)
     fzero = .25d0/(width*bag)
  else if (IC.eq.'wb_p0') then
     width = PTread_d(HCF, 'width')
     bag = PTread_d(HCF, 'bag')
     p0 = PTread_d(HCF, 'p0')
     call init_carre(F%V, width, bag, p0 = p0)
     fzero = .25d0/(width*bag)
  else if (IC.eq.'wb_b0') then
     width = zbrent(ftbz, 1d-8,4.d0*atan(1.d0), PTread_d(HCF,'b0'),1d-6)
     bag = sqrt(6.d0*PTread_d(HCF,'e0'))
     fzero = .25d0/(width*bag)
     call init_carre(F%V, width, bag)
  else if (IC.eq.'LB') then
     x1 = PTread_d(HCF, 'x1')
     x2 = PTread_d(HCF, 'x2')
     x3 = PTread_d(HCF, 'x3')
     fzero = PTread_d(HCF, 'f0')
     do i=1,F%V%Nx
        do m=1,F%V%Nv
           temp = x3*exp(-x2 * (.5d0*get_v(F%V,m)**2 + 2.d0*x1*sin(get_x(F%V,i)) + x1**2*get_v(F%V,m) + 0.5d0*x1**4) )
           F%V%f(i,m) = fzero*temp/(1+temp)
        end do
     end do
     F%I = x1
     if (x1.lt.1.d-8) then
        F%I=0.00001d0
     end if
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
  masse = sum(F%V%rho)*F%V%dx
  F%V%en_int = 2.d0*(F%Ax(1)*F%My + F%Ay(1)*F%Mx)
  F%V%en_kin = 0.d0
  F%V%momentum = 0.d0
  entropy=0.d0
  L2=0.d0
  do i=1,F%V%Nx
     do m=1,F%V%Nv
        loopf = F%V%f(i,m)
        F%V%en_kin = F%V%en_kin + get_v(F%V,m)**2 * loopf
        F%V%momentum = F%V%momentum + get_v(F%V,m) * loopf
        L2 = L2+loopf**2
        loopf=loopf/fzero
        if (loopf.gt.0.d0) then
           entropy = entropy + loopf*log(loopf)
        end if
        if (loopf.lt.1.d0) then
           entropy = entropy + (1.d0-loopf)*log(1.d0-loopf)
        end if
     end do
  end do
  F%V%en_kin = F%V%en_kin * 0.5d0 * F%V%dx * F%V%dv
  F%V%momentum = F%V%momentum     * F%V%dx * F%V%dv
  entropy = -entropy              * F%V%dx * F%V%dv
  L2           = L2               * F%V%dx * F%V%dv
  F%V%energie = F%V%en_int + F%V%en_kin
  call write_data_group_h5(h5fel, F%V, 0)
  vals = (/t_top*DT, masse, F%V%energie, F%V%en_int, F%V%en_kin,F%V%momentum, F%I, F%phi, entropy, F%Ax(1), F%Ay(1), L2/)
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

     masse = sum(F%V%rho)*F%V%dx
     F%V%en_int = 2.d0*(F%Ax(1)*F%My + F%Ay(1)*F%Mx)
     F%V%en_kin = 0.d0
     F%V%momentum = 0.d0
     entropy=0.d0
     L2=0.d0
     do i=1,F%V%Nx
        do m=1,F%V%Nv
           loopf = F%V%f(i,m)
           F%V%en_kin = F%V%en_kin + get_v(F%V,m)**2 * loopf
           F%V%momentum = F%V%momentum + get_v(F%V,m) * loopf
           L2=L2+loopf**2
           loopf=loopf/fzero
           if (loopf.gt.0.d0) then
              entropy = entropy + loopf*log(loopf)
           end if
           if (loopf.lt.1.d0) then
              entropy = entropy + (1.d0-loopf)*log(1.d0-loopf)
           end if
        end do
     end do
     F%V%en_kin = F%V%en_kin * 0.5d0 * F%V%dx * F%V%dv
     F%V%momentum = F%V%momentum     * F%V%dx * F%V%dv
     entropy = -entropy              * F%V%dx * F%V%dv
     L2      =  L2                   * F%V%dx * F%V%dv
     F%V%energie = F%V%en_int + F%V%en_kin
     vals = (/t_top*n_steps*DT, masse, F%V%energie, F%V%en_int, F%V%en_kin,F%V%momentum, F%I, F%phi, entropy, F%Ax(1), F%Ay(1), L2/)
     call write_time_slice_h5(h5fel, t_top, vals)

     if (t_top*n_images/n_top.ge.t_images) then
        call write_data_group_h5(h5fel, F%V, t_top)
        t_images = t_images + 1
     end if
     
  end do



  call close_h5(h5fel)

end program runFEL

function ftbz(s_theta,bunching)
  implicit none
  double precision s_theta, bunching, ftbz

  ftbz = sin(s_theta)/s_theta - bunching

end function ftbz
