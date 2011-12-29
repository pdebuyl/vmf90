program runAB
  use AB_module
  use ParseText
  use vmf90
  implicit none

  type(AB) :: H
  type(datafile_h5) :: h5AB
  type(PTo) :: HCF

  integer :: n_top, n_steps, t_top, t_step, n_images, t_images, t
  double precision :: DT, width, bag, e0

  integer, parameter :: time_number = 19
  double precision :: vals(time_number)
  character(len=32) :: time_names(time_number)
  character(len=32) :: IC

  
  call PTparse(HCF, 'AB_in', 10)

  call newAB(H, PTread_i(HCF,'NAx'), PTread_i(HCF,'NAv'), PTread_i(HCF,'NBx'), PTread_i(HCF,'NBv'), PTread_d(HCF, 'vmax'))
  DT = PTread_d(HCF,'DT')
  H%A%DT = DT
  H%B%DT = DT
  n_top = PTread_i(HCF,'n_top')
  n_steps = PTread_i(HCF,'n_steps')
  n_images = PTread_i(HCF, 'n_images')
  t_images = 1

  H%gamma = PTread_d(HCF, 'gamma')


  call create_AB_h5(h5AB, H, 'AB.h5', n_top, time_number)
  time_names = (/'time    ', 'masseA  ', 'masseB  ', 'energy  ', &
                 'intAB   ', 'kinA    ', 'kinB    ', 'MAx     ', &
                 'MAy     ', 'MBx     ', 'MBy     ', 'Mx      ', &
                 'My      ', 'dAx     ', 'dAy     ', 'dBx     ', &
                 'dBy     ', 'pFA     ', 'pFB     ' /)
  call write_info_string_array_h5(h5AB, 'time_names', time_names)

  call write_info_double_h5(h5AB, 'gamma', H%gamma)

  call vmf90_info()

  IC = PTread_s(HCF, 'ICA')
  if (IC.eq.'wb') then
     width = PTread_d(HCF,'widthA')
     if (width > PI) width = PI
     bag = PTread_d(HCF,'bagA')
     call init_carre(H%A, width, bag)
  else if (IC.eq.'fromHMF') then
     call load_data_from_h5(H%A, 'hmf.h5', 'data/t00400/f')
  else
     write(*,*) 'unknown IC', IC
     stop
  end if

  IC = PTread_s(HCF, 'ICB')
  if (IC.eq.'wb') then
     width = PTread_d(HCF,'widthB')
     if (width > PI) width = PI
     bag = PTread_d(HCF,'bagB')
     call init_carre(H%B, width, bag)
  else
     write(*,*) 'unknown IC', IC
     stop
  end if

  call write_grid_info(h5AB, H%A, 'AB')

  call PTkill(HCF)

  t_top = 0
  call compute_rho(H%A)
  call compute_rho(H%B)
  call compute_phi(H%A)
  call compute_phi(H%B)
  call compute_M(H)
  call write_data_group_h5(h5AB, H%A, t_top)
  call write_data_group_h5(h5AB, H%B, t_top, group2=.true.)
  call compute_phys(H, vals, t_top*n_steps*DT)
  call write_time_slice_h5(h5AB, t_top, vals)

  do t_top = 1,n_top

     call spline_x(H%A)
     call advection_x_demi(H%A)
     H%A%f = H%A%g
     call spline_x(H%B)
     call advection_x_demi(H%B)
     H%B%f = H%B%g
     do t=1,n_steps-1
        
        call compute_rho(H%A)
        call compute_rho(H%B)
        call compute_M(H)
        call compute_force(H)
        call spline_v(H%A)
        call advection_v(H%A)
        H%A%f = H%A%g
        call spline_v(H%B)
        call advection_v(H%B)
        H%B%f = H%B%g
     
        call spline_x(H%A)
        call advection_x(H%A)
        H%A%f = H%A%g
        call spline_x(H%B)
        call advection_x(H%B)
        H%B%f = H%B%g
     end do

     call compute_rho(H%A)
     call compute_rho(H%B)
     call compute_M(H)
     call compute_force(H)
     call spline_v(H%A)
     call advection_v(H%A)
     H%A%f = H%A%g
     call spline_v(H%B)
     call advection_v(H%B)
     H%B%f = H%B%g
     call spline_x(H%A)
     call advection_x_demi(H%A)
     H%A%f = H%A%g
     call spline_x(H%B)
     call advection_x_demi(H%B)
     H%B%f = H%B%g
     
     call compute_rho(H%A)
     call compute_rho(H%B)
     call compute_phi(H%A)
     call compute_phi(H%B)
     call compute_M(H)
     if (t_top*n_images/n_top.ge.t_images) then
        call write_data_group_h5(h5AB, H%A, t_top)
        call write_data_group_h5(h5AB, H%B, t_top, group2=.true.)
        t_images = t_images + 1
     end if
     call compute_phys(H, vals, t_top*n_steps*DT)
     call write_time_slice_h5(h5AB, t_top, vals)

  end do



  call close_AB_h5(h5AB)

end program runAB
