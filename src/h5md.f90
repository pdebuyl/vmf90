!! Copyright 2011-2013 Pierre de Buyl
!!
!! This file is part of f90h5md
!!
!! f90h5md is free software and is licensed under the modified BSD license (see
!! LICENSE file).

!> This module allows to handle H5MD files.

module h5md
  use hdf5
  implicit none

  !> Global variable to keep the error from HDF5 instructions.
  integer :: h5_error

  !> A type to hold an observable.
  !! Provides a buffering facility.
  type h5md_obs
     !> The ID of the observable in the file
     integer(HID_T) :: id
     !> The double precision buffer.
     double precision, allocatable :: d_buffer(:)
     !> The integer buffer.
     integer, allocatable :: i_buffer(:)
     !> The length of the buffer.
     integer :: buffer_len
     !> The index in the buffer.
     integer :: buffer_i
  end type h5md_obs

  !> A type to hold a reference to a dataset alongside with the step and time information.
  type h5md_t
     !> The ID of the dataset in the file.
     integer(HID_T) :: d_id
     !> The ID of the related step dataset.
     integer(HID_T) :: s_id
     !> The ID of the related time dataset.
     integer(HID_T) :: t_id
  end type h5md_t

  !> Interface for overloaded create_obs routines.
  !! This routine creates a dataset in a H5MD file that corresponds to the type and shape of
  !! the data argument, with a variable first dimension that corresponds to the time.
  interface h5md_create_obs
     module procedure h5md_create_obs_is
     module procedure h5md_create_obs_i1
     module procedure h5md_create_obs_i2
     module procedure h5md_create_obs_i3
     module procedure h5md_create_obs_i4
     module procedure h5md_create_obs_ds
     module procedure h5md_create_obs_d1
     module procedure h5md_create_obs_d2
     module procedure h5md_create_obs_d3
     module procedure h5md_create_obs_d4
  end interface

  !> Interface for the overloaded write_obs routines.
  !! Accepts integer or double precision, scalar, rank-1 or 2 arrays.
  !! This routine appends to the 1st dimension (0-th rank) of the appropriate dataset
  !! the actual values of the "data" argument.
  !! \code
  !! ! Assuming that ID has been created by h5md_create_obs with the same data, that is 
  !! ! that x has been passed as an argument to h5md_create_obs for ID.
  !! ! i_time is the integer timestep and time is the real time.
  !! call h5md_write_obs(ID, x, i_time, time)
  !! \endcode
  interface h5md_write_obs
     module procedure h5md_write_obs_is
     module procedure h5md_write_obs_i1
     module procedure h5md_write_obs_i2
     module procedure h5md_write_obs_i3
     module procedure h5md_write_obs_i4
     module procedure h5md_write_obs_ds
     module procedure h5md_write_obs_d1
     module procedure h5md_write_obs_d2
     module procedure h5md_write_obs_d3
     module procedure h5md_write_obs_d4
  end interface

  !> Interface for the overloaded read_obs routines.
  !! Accepts integer or double precision, scalar, rank-1 or 2 arrays.
  !! This routine appends looks up the requested integer timestep and returns the data
  !! corresponding to that timestep. It also returns the corresponding real-time.
  !! \code
  !! ! Assuming that ID has been opened by h5md_open_ID.
  !! ! i_time is the integer timestep and time is the real time.
  !! call h5md_read_obs(ID, x, i_time, time)
  !! \endcode
  interface h5md_read_obs
     module procedure h5md_read_obs_is
     module procedure h5md_read_obs_i1
     module procedure h5md_read_obs_i2
     module procedure h5md_read_obs_i3
     module procedure h5md_read_obs_i4
     module procedure h5md_read_obs_ds
     module procedure h5md_read_obs_d1
     module procedure h5md_read_obs_d2
     module procedure h5md_read_obs_d3
     module procedure h5md_read_obs_d4
  end interface

  !> Interface for the overloaded parameter routines.
  !! Accepts an integer, double precision or logical parameter that is a scalar or 
  !! a rank 1 or 2 array and also a single string of characters.
  !! \code
  !! ! Assuming timestep is declared as double precision
  !! ! and that file_id refers to an open H5MD file.
  !! timestep = 0.1d0
  !! call h5md_write_par(file_id, 'timestep', timestep)
  !! \endcode
  interface h5md_write_par
     module procedure h5md_write_par_is
     module procedure h5md_write_par_i1
     module procedure h5md_write_par_i2
     module procedure h5md_write_par_ds
     module procedure h5md_write_par_d1
     module procedure h5md_write_par_d2
     module procedure h5md_write_par_cs
     module procedure h5md_write_par_ls
     module procedure h5md_write_par_l1
     module procedure h5md_write_par_l2
  end interface

  !> Interface for the overloaded parameter routines.
  !! Accepts an integer, double precision or logical parameter that is a scalar or 
  !! a rank 1 or 2 array. Does not accept strings of characters.
  !! \code
  !! ! Assuming timestep is declared as double precision
  !! ! and that file_id refers to an open H5MD file.
  !! call h5md_read_par(file_id, 'timestep', timestep)
  !! write(*,*) 'timestep = ', timestep
  !! \endcode
  !! The data to be read is expected to be of the exact same shape than the data in
  !! the file.
  interface h5md_read_par
     module procedure h5md_read_par_is
     module procedure h5md_read_par_i1
     module procedure h5md_read_par_i2
     module procedure h5md_read_par_ds
     module procedure h5md_read_par_d1
     module procedure h5md_read_par_d2
     module procedure h5md_read_par_ls
     module procedure h5md_read_par_l1
     module procedure h5md_read_par_l2
  end interface
  
contains

  !> Creates a h5md file
  !!
  !! Creates the '/h5md' group and gives it the attributes 'creator' and 
  !! 'version'. currently, only supports creating a new file.
  !! also creates 'trajectory', 'observables' and 'parameters' groups.
  !! @param file_id the returned hdf5 location of the file.
  !! @param filename name of the file.
  !! @param creator name that appears in the 'creator' global attribute.
  subroutine h5md_create_file(file_id, filename, author, creator, creator_version)
    integer(HID_T), intent(out) :: file_id
    character(len=*), intent(in) :: filename, author, creator, creator_version

    integer :: h5md_version(2), creation_time, val(8)
    integer(HID_T) :: h5_g_id, g_id
    integer(HID_T) :: a_type, a_space, a_id
    integer(SIZE_T) :: a_size(1)
    integer :: months(12), i

    call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, h5_error)

    call h5gcreate_f(file_id, 'h5md', h5_g_id, h5_error)

    ! write author attribute
    call h5screate_f(H5S_SCALAR_F, a_space, h5_error)
    a_size(1) = len(author)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, a_type, h5_error)
    call h5tset_size_f(a_type, a_size(1), h5_error)
    call h5acreate_f(h5_g_id, 'author', a_type, a_space, a_id, h5_error)
    call h5awrite_f(a_id, a_type, author, a_size, h5_error)
    call h5aclose_f(a_id, h5_error)
    call h5sclose_f(a_space, h5_error)
    call h5tclose_f(a_type, h5_error)

    ! write creator attribute
    call h5screate_f(H5S_SCALAR_F, a_space, h5_error)
    a_size(1) = len(creator)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, a_type, h5_error)
    call h5tset_size_f(a_type, a_size(1), h5_error)
    call h5acreate_f(h5_g_id, 'creator', a_type, a_space, a_id, h5_error)
    call h5awrite_f(a_id, a_type, creator, a_size, h5_error)
    call h5aclose_f(a_id, h5_error)
    call h5sclose_f(a_space, h5_error)
    call h5tclose_f(a_type, h5_error)

    ! write creator_version attribute
    call h5screate_f(H5S_SCALAR_F, a_space, h5_error)
    a_size(1) = len(creator_version)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, a_type, h5_error)
    call h5tset_size_f(a_type, a_size(1), h5_error)
    call h5acreate_f(h5_g_id, 'creator_version', a_type, a_space, a_id, h5_error)
    call h5awrite_f(a_id, a_type, creator_version, a_size, h5_error)
    call h5aclose_f(a_id, h5_error)
    call h5sclose_f(a_space, h5_error)
    call h5tclose_f(a_type, h5_error)
    
    ! write h5md/version
    h5md_version = (/ 0, 1 /)
    a_size(1) = 2
    call h5screate_simple_f(1, a_size, a_space, h5_error)
    call h5acreate_f(h5_g_id, 'version', H5T_NATIVE_INTEGER, a_space, a_id, h5_error)
    call h5awrite_f(a_id, H5T_NATIVE_INTEGER, h5md_version, a_size, h5_error)
    call h5aclose_f(a_id, h5_error)
    call h5sclose_f(a_space, h5_error)

    ! write h5md/creation_time as seconds since Epoch
    months = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
    call date_and_time(values=val)
    creation_time = sum(months) * 24*60*60 * (val(1)-1970)
    if (val(2).gt.1) then
       creation_time = creation_time + sum(months(1:val(2)-1)) * 24*60*60
    end if
    do i=1970,val(1)
       if (mod(i,4).eq.0) then
          if (mod(i,100).eq.0) then
             if (mod(i,400).eq.0) creation_time = creation_time + 24*60*60
          else
             creation_time = creation_time + 24*60*60
          end if
       end if
    end do
    creation_time = creation_time + (val(3)-1)*24*60*60
    creation_time = creation_time - val(4)*60  
    creation_time = creation_time + val(5)*60*60
    creation_time = creation_time + val(6)*60
    creation_time = creation_time + val(7)

    a_size(1) = 1
    call h5screate_f(H5S_SCALAR_F, a_space, h5_error)
    call h5acreate_f(h5_g_id, 'creation_time', H5T_NATIVE_INTEGER, a_space, a_id, h5_error)
    call h5awrite_f(a_id, H5T_NATIVE_INTEGER, creation_time, a_size, h5_error)
    call h5aclose_f(a_id, h5_error)
    call h5sclose_f(a_space, h5_error)

    call h5gclose_f(h5_g_id, h5_error)

    call h5gcreate_f(file_id, 'trajectory', g_id, h5_error)
    call h5gclose_f(g_id, h5_error)

    call h5gcreate_f(file_id, 'observables', g_id, h5_error)
    call h5gclose_f(g_id, h5_error)

    call h5gcreate_f(file_id, 'parameters', g_id, h5_error)
    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_file

  !> Opens a h5md file
  !!
  !! Returns file_id as the location of the opened file.
  !! @param file_id the returned hdf5 location of the file.
  !! @param filename name of the file.
  !! @param rw flag that allows to open the file in read/write mode.
  subroutine h5md_open_file(file_id, filename, rw)
    integer(HID_T), intent(out) :: file_id
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: rw
    
    if (present(rw) .and. rw) then
       call h5fopen_f(filename, H5F_ACC_RDWR_F, file_id, h5_error)
    else
       call h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, h5_error)
    end if

  end subroutine h5md_open_file

  !> Adds a trajectory group in a h5md file
  !! @param file_id the HDF5 ID of the file.
  !! @param group_name name of a subgroup of 'trajectory'.
  subroutine h5md_create_trajectory_group(file_id, group_name)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: group_name

    integer(HID_T) :: traj_g_id
    
    call h5gcreate_f(file_id, 'trajectory/'//group_name , traj_g_id, h5_error)

    call h5gclose_f(traj_g_id, h5_error)

  end subroutine h5md_create_trajectory_group

  !> Creates the step and time datasets in the specified group
  !! @param group_id location to place the datasets.
  subroutine h5md_create_step_time(group_id)
    integer(HID_T), intent(inout) :: group_id

    integer :: rank
    integer(HSIZE_T) :: dims(1), max_dims(1), chunk_dims(1)
    integer(HID_T) :: s_id, d_id, plist

    rank = 1
    dims = (/ 0 /)
    max_dims = (/ H5S_UNLIMITED_F /)
    call h5screate_simple_f(rank, dims, s_id, h5_error, max_dims)
    chunk_dims = (/ 1024 /)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(group_id, 'step', H5T_NATIVE_INTEGER, s_id, d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5dclose_f(d_id, h5_error)
    call h5sclose_f(s_id, h5_error)
    
    rank = 1
    dims = (/ 0 /)
    max_dims = (/ H5S_UNLIMITED_F /)
    call h5screate_simple_f(rank, dims, s_id, h5_error, max_dims)
    chunk_dims = (/ 1024 /)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(group_id, 'time', H5T_NATIVE_DOUBLE, s_id, d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5dclose_f(d_id, h5_error)
    call h5sclose_f(s_id, h5_error)
    
  end subroutine h5md_create_step_time

  !> Adds a trajectory dataset to a trajectory group of name trajectory_name
  !! @param file_id the HDF5 ID of the file.
  !! @param trajectory_name can be 'position', 'velocity', 'force' or 'species'
  !! @param N the number of atoms
  !! @param D the spatial dimension
  !! @param ID the returned h5md_t value of the created dataset.
  !! @param group_name optional group name for the trajectory group
  !! @param species_react optional argument. if set to .true., 'species' will be
  !! time dependent, if set to .false., 'species' will not possess the time
  !! dimension
  !! @param link_from is the name of another trajectory from which the time can be copied
  !! @param compress Switch to toggle GZIP compression.
  !! @param force_kind Force the supplied kind, "integer" or "double", for the dataset.
  !! @param force_rank Force the dataset rank to 1, 2 or 3.
  subroutine h5md_add_trajectory_data(file_id, trajectory_name, N, D, ID, group_name, species_react, &
       link_from, compress, force_kind, force_rank)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: trajectory_name
    integer, intent(in) :: N, D
    type(h5md_t), intent(out) :: ID
    character(len=*), intent(in), optional :: group_name
    logical, intent(in), optional :: species_react
    character(len=*), intent(in), optional :: link_from
    logical, intent(in), optional :: compress
    character(len=*), intent(in), optional :: force_kind
    integer, intent(in), optional :: force_rank
    
    character(len=128) :: path
    integer :: rank
    integer(HSIZE_T) :: dims(3), max_dims(3), chunk_dims(3)
    integer(HID_T) :: traj_g_id, g_id, s_id, plist
    logical :: gzip_avail, compress_var
    integer :: filter_info, gzip_encode

    if (.not.present(force_kind)) then
    if ( (trajectory_name .ne. 'position') .and. (trajectory_name .ne. 'velocity') &
         .and. (trajectory_name .ne. 'jumps') &
         .and. (trajectory_name .ne. 'force') .and. (trajectory_name .ne. 'species') ) then
       write(*,*) 'non conforming trajectory name in h5md_add_trajectory_data'
       stop
    end if
    end if
    
    if (present(compress)) then
       compress_var = compress
    else
       compress_var = .false.
    end if

    if (present(group_name)) then
       call h5gopen_f(file_id, 'trajectory/'//group_name, traj_g_id, h5_error)
    else
       call h5gopen_f(file_id, 'trajectory', traj_g_id, h5_error)
    end if
    path = trajectory_name

    ! g_id is opened as the container of trajectory_name
    call h5gcreate_f(traj_g_id, path, g_id, h5_error)
       
    if (present(force_rank)) then
       rank = force_rank
    else
       if (trajectory_name .eq. 'species') then
          if (present(species_react) .and. (species_react) ) then
             rank=2
          else
             rank=1
          end if
       else
          rank=3
       end if
    end if
    if (rank.eq.3) then
       dims = (/ D, N, 0 /)
       max_dims = (/ D, N, H5S_UNLIMITED_F /)
       chunk_dims = (/ D, N, 1 /)
    else if (rank.eq.2) then
       dims = (/ N, 0, 0 /)
       max_dims = (/ N, H5S_UNLIMITED_F, 0 /)
       chunk_dims = (/ N, 1, 0 /)
    else if (rank.eq.1) then
       dims = (/ N, 0, 0 /)
       max_dims = (/ N, 0, 0 /)
       chunk_dims = (/ N, 0, 0 /)
    else
       write(*,*) 'the rank ', rank, ' is inappropriate in h5md_add_trajectory_data'
    end if

    call h5screate_simple_f(rank, dims, s_id, h5_error, max_dims)
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    if (compress_var) then
       call h5zfilter_avail_f(H5Z_FILTER_DEFLATE_F,gzip_avail, h5_error)
       if (.not.gzip_avail) stop 'GZIP filter not available'
       CALL h5zget_filter_info_f(H5Z_FILTER_DEFLATE_F, filter_info, h5_error)
       gzip_encode = IOR(H5Z_FILTER_ENCODE_ENABLED_F,H5Z_FILTER_DECODE_ENABLED_F)
       if (gzip_encode .ne. filter_info) stop 'GZIP filter not available for encoding and decoding'
       call h5pset_deflate_f(plist, 6, h5_error)
    end if
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    if (present(force_kind)) then
       if (force_kind .eq. 'integer') then
          call h5dcreate_f(g_id, 'value', H5T_NATIVE_INTEGER, s_id, ID% d_id, h5_error, plist)
       else if (force_kind .eq. 'double') then
          call h5dcreate_f(g_id, 'value', H5T_NATIVE_DOUBLE, s_id, ID% d_id, h5_error, plist)
       else
          write(*,*) 'non supported force_kind ', force_kind
       end if
    else
    if (trajectory_name .ne. 'species' .and. trajectory_name .ne. 'jumps') then
       call h5dcreate_f(g_id, 'value', H5T_NATIVE_DOUBLE, s_id, ID% d_id, h5_error, plist)
    else
       call h5dcreate_f(g_id, 'value', H5T_NATIVE_INTEGER, s_id, ID% d_id, h5_error, plist)
    end if
    end if
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(s_id, h5_error)

    if (present(link_from)) then
       call h5lcreate_hard_f(traj_g_id, link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(traj_g_id, link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)
    call h5gclose_f(traj_g_id, h5_error)
    
  end subroutine h5md_add_trajectory_data

  !> Sets the minimum and maximum attributes to a position trajectory dataset.
  !! @param ID h5md_t of the dataset
  !! @param xmin lower coordinates of the box
  !! @param xmax upper coordinates of the box
  subroutine h5md_set_box_size(ID, xmin, xmax)
    type(h5md_t), intent(inout) :: ID
    double precision, intent(in) :: xmin(:)
    double precision, intent(in) :: xmax(:)

    integer(HID_T) :: att_id, att_s
    integer(HSIZE_T) :: att_size(1)

    ! Checking that xmin and xmax have the same size
    if (size(xmin) /= size(xmax)) stop 'non equal dimensions for xmin and xmax in h5md_set_box_size'
    
    ! Setting the size for the attribute
    att_size(1) = size(xmin)

    ! Writing the minimum attribute
    call h5screate_simple_f(1, att_size, att_s, h5_error)
    call h5acreate_f(ID%d_id, 'minimum', H5T_NATIVE_DOUBLE, att_s, att_id, h5_error)
    call h5awrite_f(att_id, H5T_NATIVE_DOUBLE, xmin, att_size, h5_error)
    call h5aclose_f(att_id, h5_error)
    call h5sclose_f(att_s, h5_error)

    ! Writing the maximum attribute
    call h5screate_simple_f(1, att_size, att_s, h5_error)
    call h5acreate_f(ID%d_id, 'maximum', H5T_NATIVE_DOUBLE, att_s, att_id, h5_error)
    call h5awrite_f(att_id, H5T_NATIVE_DOUBLE, xmax, att_size, h5_error)
    call h5aclose_f(att_id, h5_error)
    call h5sclose_f(att_s, h5_error)

  end subroutine h5md_set_box_size

  !> Opens a trajectory dataset and its associated step and time datasets.
  !! @param file_id the HDF5 ID of the file.
  !! @param trajectory_name The name of the trajectory: position, velocity, force or species.
  !! @param ID The h5md_t variable for the dataset.
  !! @param group_name The optional group name.
  subroutine h5md_open_trajectory(file_id, trajectory_name, ID, group_name)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: trajectory_name
    type(h5md_t), intent(out) :: ID
    character(len=*), intent(in), optional :: group_name

    if (present(group_name)) then
       call h5dopen_f(file_id, 'trajectory/'//group_name//'/'//trajectory_name//'/coordinates', ID%d_id, h5_error)
       call h5dopen_f(file_id, 'trajectory/'//group_name//'/'//trajectory_name//'/step', ID%s_id, h5_error)
       call h5dopen_f(file_id, 'trajectory/'//group_name//'/'//trajectory_name//'/time', ID%t_id, h5_error)
    else
       call h5dopen_f(file_id, 'trajectory/'//trajectory_name//'/coordinates', ID%d_id, h5_error)
       call h5dopen_f(file_id, 'trajectory/'//trajectory_name//'/step', ID%s_id, h5_error)
       call h5dopen_f(file_id, 'trajectory/'//trajectory_name//'/time', ID%t_id, h5_error)
    end if

  end subroutine h5md_open_trajectory


  !> Appends step and time information.
  !! If the last step is the present step, does nothing. Else, append the data.
  !! @param s_id ID of the step dataset
  !! @param t_id ID of the time dataset
  !! @param step The present step
  !! @param time The present time
  subroutine h5md_append_step_time(s_id, t_id, step, time)
    integer(HID_T), intent(inout) :: s_id, t_id
    integer, intent(in) :: step
    double precision, intent(in) :: time

    integer(HID_T) :: file_s, mem_s
    integer(HSIZE_T) :: dims(1), max_dims(1), start(1), num(1)
    integer :: last_step(1)

    ! open step
    call h5dget_space_f(s_id, file_s, h5_error)
    call h5sget_simple_extent_dims_f(file_s, dims, max_dims, h5_error)

    if (dims(1) .le. 0) then
       last_step(1) = -1
    else
       start(1) = dims(1)-1
       num(1) = 1
       call h5screate_simple_f(1, num, mem_s, h5_error)
       call h5sselect_hyperslab_f(file_s, H5S_SELECT_SET_F, start, num, h5_error)
       call h5dread_f(s_id, H5T_NATIVE_INTEGER, last_step, num, h5_error, mem_space_id=mem_s, file_space_id=file_s)
       call h5sclose_f(mem_s, h5_error)
    end if

    call h5sclose_f(file_s, h5_error)

    ! check last
    if (last_step(1) .gt. step) then ! if last > present_step -> error
       write(*,*) 'error, last step is bigger than present step'
    else if (last_step(1) .lt. step) then ! else if last < present_step -> extend step and append present_step, same for time
       ! add step value to the end of the step dataset

       call h5dget_space_f(s_id, file_s, h5_error)
       dims(1) = 1
       call h5screate_simple_f(1, dims, mem_s, h5_error)
       call h5sget_simple_extent_dims_f(file_s, dims, max_dims, h5_error)
       call h5sclose_f(file_s, h5_error)
       start(1) = dims(1)
       num(1) = 1
       dims(1) = dims(1) + 1
       call h5dset_extent_f(s_id, dims, h5_error)
       call h5dget_space_f(s_id, file_s, h5_error)

       call h5sselect_hyperslab_f(file_s, H5S_SELECT_SET_F, start, num, h5_error)
       call h5dwrite_f(s_id, H5T_NATIVE_INTEGER, step, num, h5_error, mem_space_id=mem_s, file_space_id=file_s)
       call h5sclose_f(file_s, h5_error)
       call h5sclose_f(mem_s, h5_error)

       ! add time value to the end of the time dataset
       dims(1) = 1
       call h5screate_simple_f(1, dims, mem_s, h5_error)

       call h5dget_space_f(t_id, file_s, h5_error)
       call h5sget_simple_extent_dims_f(file_s, dims, max_dims, h5_error)
       call h5sclose_f(file_s, h5_error)
       dims(1) = dims(1) + 1
       call h5dset_extent_f(t_id, dims, h5_error)
       call h5dget_space_f(t_id, file_s, h5_error)
       call h5sget_simple_extent_dims_f(file_s, dims, max_dims, h5_error)
       start(1) = dims(1) - 1
       num(1) = 1
       call h5sselect_hyperslab_f(file_s, H5S_SELECT_SET_F, start, num, h5_error)
       call h5dwrite_f(t_id, H5T_NATIVE_DOUBLE, time, num, h5_error, mem_space_id=mem_s, file_space_id=file_s)
       call h5sclose_f(file_s, h5_error)
       call h5sclose_f(mem_s, h5_error)

    end if ! else if last = present_step, do nothing

  end subroutine h5md_append_step_time


  !> Close a h5md_obs variable.
  !! @param obs h5md_obs variable.
  subroutine h5md_close_obs(obs)
    type(h5md_obs), intent(inout) :: obs

    if (obs% buffer_i .gt. 0) then
       !call h5md_append_obs_value_d(obs, force_dump = .true.)
    end if

    call h5dclose_f(obs% id, h5_error)

  end subroutine h5md_close_obs
  
  !> Close a h5md_t variable.
  !! @param ID h5md_t variable.
  subroutine h5md_close_ID(ID)
    type(h5md_t), intent(inout) :: ID

    call h5dclose_f(ID% d_id, h5_error)
    call h5dclose_f(ID% s_id, h5_error)
    call h5dclose_f(ID% t_id, h5_error)

  end subroutine h5md_close_ID

  !> Opens a h5md_t variable.
  !! @param file_ID The ID of the HDF5 file. 
  !! @param ID h5md_t variable.
  !! @param h5md_group The name of the group, it can be either 'trajectory' or 'observables'.
  !! @param full_name The name of the data, including its group, e.g. 'solvent/position'.
  subroutine h5md_open_ID(file_ID, ID, h5md_group, full_name)
    integer(HID_T), intent(in) :: file_ID
    type(h5md_t), intent(inout) :: ID
    character(len=*), intent(in) :: h5md_group, full_name

    if ( (h5md_group .ne. 'trajectory') &
         .and. &
         (h5md_group .ne. 'observables') &
         ) then
       write(*,*) h5md_group, ' is not accepted in h5md_open_ID'
       stop
    end if
    call h5dopen_f(file_ID, h5md_group//'/'//full_name//'/value',ID% d_id, h5_error)
    call h5dopen_f(file_ID, h5md_group//'/'//full_name//'/step',ID% s_id, h5_error)
    call h5dopen_f(file_ID, h5md_group//'/'//full_name//'/time',ID% t_id, h5_error)
       
  end subroutine h5md_open_ID

  !> Provides the idx that corresponds to the integer timestep step.
  !! Also provides the corresponding real-valued time.
  !! @param ID The h5md_t variable for the data.
  !! @param step The integer timestep requested.
  !! @param idx Is set to the idx of the array to fetch.
  !! @param time Is set to the corresponding real time.
  subroutine h5md_get_step_time(ID, step, idx, time)
    implicit none
    type(h5md_t), intent(inout) :: ID
    integer, intent(in) :: step
    integer, intent(out) :: idx
    double precision, intent(out) :: time
    
    integer(HSIZE_T) :: dims(1), max_dims(1)
    integer(HID_T) :: space
    integer, allocatable :: step_data(:)
    double precision, allocatable :: time_data(:)
    integer :: rank
    logical :: ok

    call h5dget_space_f(ID% s_id, space, h5_error)
    call h5sget_simple_extent_ndims_f(space, rank, h5_error)
    if (rank .ne. 1) stop 'wrong rank for step in h5md_get_step_time'
    call h5sget_simple_extent_dims_f(space, dims, max_dims, h5_error)
    call h5sclose_f(space, h5_error)
    allocate(step_data(dims(1)))
    call h5dread_f(ID% s_id, H5T_NATIVE_INTEGER, step_data, dims, h5_error)

    call h5md_lookup(step_data, step, idx, ok)
    if (.not. ok) stop 'could not find step in h5md_get_step_time'
    deallocate(step_data)

    call h5dget_space_f(ID% t_id, space, h5_error)
    call h5sget_simple_extent_ndims_f(space, rank, h5_error)
    if (rank .ne. 1) stop 'wrong rank for step in h5md_get_step_time'
    call h5sget_simple_extent_dims_f(space, dims, max_dims, h5_error)
    call h5sclose_f(space, h5_error)
    allocate(time_data(dims(1)))
    call h5dread_f(ID% t_id, H5T_NATIVE_DOUBLE, time_data, dims, h5_error)
    time = time_data(idx)
    deallocate(time_data)

  end subroutine h5md_get_step_time

  !> Looks up an integer in an ordered non-repeating list.
  !! @param list An ordered non-repeating list of integers.
  !! @param wish The integer whose index is looked for.
  !! @param idx The resulting index.
  !! @param success Results in .true. if the operation was successful.
  subroutine h5md_lookup(list, wish, idx, success)
    implicit none
    integer, intent(in) :: list(:), wish
    integer, intent(out) :: idx
    logical, intent(out) :: success

    integer :: low, up, i
    integer, allocatable :: order_check(:)    
    
    success = .false.

    if (size(list) .eq. 1) then
       if (list(1) .eq. wish) then
          success = .true.
          idx = 1
       end if
       return
    end if

    allocate(order_check(size(list)-1))
    order_check = list(2:size(list))-list(1:size(list)-1)
    low = minval(order_check)
    deallocate(order_check)
    if (low <= 0) then 
       write(*,*) 'nonordered or repeating list in lookup'
       return
    end if

    if ( wish < list(1) .or. wish > list(size(list)) ) then
       write(*,*) wish, ' out of range in lookup'
       stop
    end if

    low = 1
    up = size(list)
    lookup_loop: do
       
       i = (low+up)/2
       if (list(i) >= wish) then
          up =i
       else
          low=i
       end if
       if ( up-low .eq. 1 ) then
          if (list(up).eq.wish) i=up
          if (list(low).eq.wish) i=low
          exit lookup_loop
       end if
    end do lookup_loop
  
    if (list(i).eq.wish) success=.true.
    idx = i
    
  end subroutine h5md_lookup

  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param data The data that will fit into the observable.
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  !! @private
  subroutine h5md_create_obs_i1(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    integer, intent(in) :: data(:)
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 2
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1:1) = shape(data)
    max_dims(1:1) = shape(data)
    chunk_dims(1:1) = shape(data)
    

    dims(2)     = 0
    max_dims(2) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(2) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'value', H5T_NATIVE_INTEGER, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_i1


  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param data The data that will fit into the observable.
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  !! @private
  subroutine h5md_create_obs_is(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    integer, intent(in) :: data
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 1
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1)     = 0
    max_dims(1) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(1) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'value', H5T_NATIVE_INTEGER, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_is


  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param data The data that will fit into the observable.
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  !! @private
  subroutine h5md_create_obs_i2(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    integer, intent(in) :: data(:,:)
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 3
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1:2) = shape(data)
    max_dims(1:2) = shape(data)
    chunk_dims(1:2) = shape(data)
    

    dims(3)     = 0
    max_dims(3) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(3) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'value', H5T_NATIVE_INTEGER, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_i2


  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param data The data that will fit into the observable.
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  !! @private
  subroutine h5md_create_obs_i3(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    integer, intent(in) :: data(:,:,:)
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 4
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1:3) = shape(data)
    max_dims(1:3) = shape(data)
    chunk_dims(1:3) = shape(data)
    

    dims(4)     = 0
    max_dims(4) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(4) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'value', H5T_NATIVE_INTEGER, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_i3


  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param data The data that will fit into the observable.
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  !! @private
  subroutine h5md_create_obs_i4(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    integer, intent(in) :: data(:,:,:,:)
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 5
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1:4) = shape(data)
    max_dims(1:4) = shape(data)
    chunk_dims(1:4) = shape(data)
    

    dims(5)     = 0
    max_dims(5) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(5) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'value', H5T_NATIVE_INTEGER, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_i4


  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param data The data that will fit into the observable.
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  !! @private
  subroutine h5md_create_obs_d1(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    double precision, intent(in) :: data(:)
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 2
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1:1) = shape(data)
    max_dims(1:1) = shape(data)
    chunk_dims(1:1) = shape(data)
    

    dims(2)     = 0
    max_dims(2) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(2) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'value', H5T_NATIVE_DOUBLE, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_d1


  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param data The data that will fit into the observable.
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  !! @private
  subroutine h5md_create_obs_ds(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    double precision, intent(in) :: data
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 1
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1)     = 0
    max_dims(1) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(1) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'value', H5T_NATIVE_DOUBLE, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_ds


  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param data The data that will fit into the observable.
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  !! @private
  subroutine h5md_create_obs_d2(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    double precision, intent(in) :: data(:,:)
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 3
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1:2) = shape(data)
    max_dims(1:2) = shape(data)
    chunk_dims(1:2) = shape(data)
    

    dims(3)     = 0
    max_dims(3) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(3) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'value', H5T_NATIVE_DOUBLE, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_d2


  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param data The data that will fit into the observable.
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  !! @private
  subroutine h5md_create_obs_d3(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    double precision, intent(in) :: data(:,:,:)
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 4
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1:3) = shape(data)
    max_dims(1:3) = shape(data)
    chunk_dims(1:3) = shape(data)
    

    dims(4)     = 0
    max_dims(4) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(4) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'value', H5T_NATIVE_DOUBLE, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_d3


  !> Sets up a h5md_t variable.
  !! @param file_id ID of the file.
  !! @param name Name of the observable
  !! @param ID Resulting h5md_t variable
  !! @param data The data that will fit into the observable.
  !! @param link_from Indicates if the step and time for this observable should be linked from another one.
  !! @private
  subroutine h5md_create_obs_d4(file_id, name, ID, data, link_from)
    integer(HID_T), intent(inout) :: file_id
    character(len=*), intent(in) :: name
    type(h5md_t), intent(out) :: ID
    double precision, intent(in) :: data(:,:,:,:)
    character(len=*), intent(in), optional :: link_from

    integer(HID_T) :: file_s, plist, g_id
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), chunk_dims(:)
    integer :: rank

    call h5gcreate_f(file_id, 'observables/'//name, g_id, h5_error)

    rank = 5
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(chunk_dims(rank))

    dims(1:4) = shape(data)
    max_dims(1:4) = shape(data)
    chunk_dims(1:4) = shape(data)
    

    dims(5)     = 0
    max_dims(5) = H5S_UNLIMITED_F
    call h5screate_simple_f(rank, dims, file_s, h5_error, max_dims)
    chunk_dims(5) = 128
    call h5pcreate_f(H5P_DATASET_CREATE_F, plist, h5_error)
    call h5pset_chunk_f(plist, rank, chunk_dims, h5_error)
    call h5dcreate_f(g_id, 'value', H5T_NATIVE_DOUBLE, file_s, ID% d_id, h5_error, plist)
    call h5pclose_f(plist, h5_error)
    call h5sclose_f(file_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(chunk_dims)

    if (present(link_from)) then
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/step', g_id, 'step', h5_error)
       call h5lcreate_hard_f(file_id, 'observables/'//link_from//'/time', g_id, 'time', h5_error)
    else
       call h5md_create_step_time(g_id)
    end if

    call h5dopen_f(g_id, 'step', ID% s_id, h5_error)
    call h5dopen_f(g_id, 'time', ID% t_id, h5_error)

    call h5gclose_f(g_id, h5_error)

  end subroutine h5md_create_obs_d4

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  !! @private
  subroutine h5md_write_obs_i1(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    integer, intent(in) :: data(:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T) :: dims(2), max_dims(2), start(2), num(2)


    dims(1:1) = shape(data)

    call h5screate_simple_f(1, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)

    start(1:1) = 0
    num(1:1) = dims(1:1)
    start(2) = dims(2)
    num(2) = 1
    dims(2) = dims(2) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_INTEGER, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_i1

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  !! @private
  subroutine h5md_write_obs_is(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    integer, intent(in) :: data
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T) :: dims(1), max_dims(1), start(1), num(1)


    dims(1) = 1

    call h5screate_simple_f(0, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)
    start(1) = dims(1)
    num(1) = 1
    dims(1) = dims(1) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_INTEGER, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_is

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  !! @private
  subroutine h5md_write_obs_i2(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    integer, intent(in) :: data(:,:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T) :: dims(3), max_dims(3), start(3), num(3)


    dims(1:2) = shape(data)

    call h5screate_simple_f(2, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)

    start(1:2) = 0
    num(1:2) = dims(1:2)
    start(3) = dims(3)
    num(3) = 1
    dims(3) = dims(3) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_INTEGER, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_i2

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  !! @private
  subroutine h5md_write_obs_i3(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    integer, intent(in) :: data(:,:,:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T) :: dims(4), max_dims(4), start(4), num(4)


    dims(1:3) = shape(data)

    call h5screate_simple_f(3, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)

    start(1:3) = 0
    num(1:3) = dims(1:3)
    start(4) = dims(4)
    num(4) = 1
    dims(4) = dims(4) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_INTEGER, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_i3

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  !! @private
  subroutine h5md_write_obs_i4(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    integer, intent(in) :: data(:,:,:,:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T) :: dims(5), max_dims(5), start(5), num(5)


    dims(1:4) = shape(data)

    call h5screate_simple_f(4, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)

    start(1:4) = 0
    num(1:4) = dims(1:4)
    start(5) = dims(5)
    num(5) = 1
    dims(5) = dims(5) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_INTEGER, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_i4

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  !! @private
  subroutine h5md_write_obs_d1(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    double precision, intent(in) :: data(:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T) :: dims(2), max_dims(2), start(2), num(2)


    dims(1:1) = shape(data)

    call h5screate_simple_f(1, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)

    start(1:1) = 0
    num(1:1) = dims(1:1)
    start(2) = dims(2)
    num(2) = 1
    dims(2) = dims(2) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_DOUBLE, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_d1

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  !! @private
  subroutine h5md_write_obs_ds(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    double precision, intent(in) :: data
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T) :: dims(1), max_dims(1), start(1), num(1)


    dims(1) = 1

    call h5screate_simple_f(0, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)
    start(1) = dims(1)
    num(1) = 1
    dims(1) = dims(1) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_DOUBLE, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_ds

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  !! @private
  subroutine h5md_write_obs_d2(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    double precision, intent(in) :: data(:,:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T) :: dims(3), max_dims(3), start(3), num(3)


    dims(1:2) = shape(data)

    call h5screate_simple_f(2, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)

    start(1:2) = 0
    num(1:2) = dims(1:2)
    start(3) = dims(3)
    num(3) = 1
    dims(3) = dims(3) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_DOUBLE, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_d2

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  !! @private
  subroutine h5md_write_obs_d3(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    double precision, intent(in) :: data(:,:,:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T) :: dims(4), max_dims(4), start(4), num(4)


    dims(1:3) = shape(data)

    call h5screate_simple_f(3, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)

    start(1:3) = 0
    num(1:3) = dims(1:3)
    start(4) = dims(4)
    num(4) = 1
    dims(4) = dims(4) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_DOUBLE, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_d3

  !> Takes a single value and appends it to the appropriate dataset.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param present_step integer time step.
  !! @param time real-valued time.
  !! @private
  subroutine h5md_write_obs_d4(ID, data, present_step, time)
    type(h5md_t), intent(inout) :: ID
    double precision, intent(in) :: data(:,:,:,:)
    integer, intent(in) :: present_step
    double precision, intent(in) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T) :: dims(5), max_dims(5), start(5), num(5)


    dims(1:4) = shape(data)

    call h5screate_simple_f(4, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    call h5sclose_f(obs_s, h5_error)

    start(1:4) = 0
    num(1:4) = dims(1:4)
    start(5) = dims(5)
    num(5) = 1
    dims(5) = dims(5) + 1
    call h5dset_extent_f(ID% d_id, dims, h5_error)
    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dwrite_f(ID% d_id, H5T_NATIVE_DOUBLE, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    call h5md_append_step_time(ID% s_id, ID% t_id, present_step, time)
       
  end subroutine h5md_write_obs_d4

  !> Reads a time frame for the given observable.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param step integer time step.
  !! @param time The time corresponding to the step.
  !! @private
  subroutine h5md_read_obs_i1(ID, data, step, time)
    implicit none
    type(h5md_t), intent(inout) :: ID
    integer, intent(out) :: data(:)
    integer, intent(in) :: step
    double precision, intent(out) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank, idx

    rank = 2
    call h5md_get_step_time(ID, step, idx, time)
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1:1) = shape(data)

    call h5screate_simple_f(1, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)

    start(1:1) = 0
    num(1:1) = dims(1:1)
    start(2) = idx-1
    num(2) = 1
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dread_f(ID% d_id, H5T_NATIVE_INTEGER, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)
       
  end subroutine h5md_read_obs_i1

  !> Reads a time frame for the given observable.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param step integer time step.
  !! @param time The time corresponding to the step.
  !! @private
  subroutine h5md_read_obs_is(ID, data, step, time)
    implicit none
    type(h5md_t), intent(inout) :: ID
    integer, intent(out) :: data
    integer, intent(in) :: step
    double precision, intent(out) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank, idx

    rank = 1
    call h5md_get_step_time(ID, step, idx, time)
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1) = 1

    call h5screate_simple_f(0, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    start(1) = idx-1
    num(1) = 1
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dread_f(ID% d_id, H5T_NATIVE_INTEGER, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)
       
  end subroutine h5md_read_obs_is

  !> Reads a time frame for the given observable.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param step integer time step.
  !! @param time The time corresponding to the step.
  !! @private
  subroutine h5md_read_obs_i2(ID, data, step, time)
    implicit none
    type(h5md_t), intent(inout) :: ID
    integer, intent(out) :: data(:,:)
    integer, intent(in) :: step
    double precision, intent(out) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank, idx

    rank = 3
    call h5md_get_step_time(ID, step, idx, time)
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1:2) = shape(data)

    call h5screate_simple_f(2, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)

    start(1:2) = 0
    num(1:2) = dims(1:2)
    start(3) = idx-1
    num(3) = 1
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dread_f(ID% d_id, H5T_NATIVE_INTEGER, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)
       
  end subroutine h5md_read_obs_i2

  !> Reads a time frame for the given observable.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param step integer time step.
  !! @param time The time corresponding to the step.
  !! @private
  subroutine h5md_read_obs_i3(ID, data, step, time)
    implicit none
    type(h5md_t), intent(inout) :: ID
    integer, intent(out) :: data(:,:,:)
    integer, intent(in) :: step
    double precision, intent(out) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank, idx

    rank = 4
    call h5md_get_step_time(ID, step, idx, time)
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1:3) = shape(data)

    call h5screate_simple_f(3, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)

    start(1:3) = 0
    num(1:3) = dims(1:3)
    start(4) = idx-1
    num(4) = 1
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dread_f(ID% d_id, H5T_NATIVE_INTEGER, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)
       
  end subroutine h5md_read_obs_i3

  !> Reads a time frame for the given observable.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param step integer time step.
  !! @param time The time corresponding to the step.
  !! @private
  subroutine h5md_read_obs_i4(ID, data, step, time)
    implicit none
    type(h5md_t), intent(inout) :: ID
    integer, intent(out) :: data(:,:,:,:)
    integer, intent(in) :: step
    double precision, intent(out) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank, idx

    rank = 5
    call h5md_get_step_time(ID, step, idx, time)
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1:4) = shape(data)

    call h5screate_simple_f(4, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)

    start(1:4) = 0
    num(1:4) = dims(1:4)
    start(5) = idx-1
    num(5) = 1
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dread_f(ID% d_id, H5T_NATIVE_INTEGER, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)
       
  end subroutine h5md_read_obs_i4

  !> Reads a time frame for the given observable.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param step integer time step.
  !! @param time The time corresponding to the step.
  !! @private
  subroutine h5md_read_obs_d1(ID, data, step, time)
    implicit none
    type(h5md_t), intent(inout) :: ID
    double precision, intent(out) :: data(:)
    integer, intent(in) :: step
    double precision, intent(out) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank, idx

    rank = 2
    call h5md_get_step_time(ID, step, idx, time)
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1:1) = shape(data)

    call h5screate_simple_f(1, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)

    start(1:1) = 0
    num(1:1) = dims(1:1)
    start(2) = idx-1
    num(2) = 1
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dread_f(ID% d_id, H5T_NATIVE_DOUBLE, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)
       
  end subroutine h5md_read_obs_d1

  !> Reads a time frame for the given observable.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param step integer time step.
  !! @param time The time corresponding to the step.
  !! @private
  subroutine h5md_read_obs_ds(ID, data, step, time)
    implicit none
    type(h5md_t), intent(inout) :: ID
    double precision, intent(out) :: data
    integer, intent(in) :: step
    double precision, intent(out) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank, idx

    rank = 1
    call h5md_get_step_time(ID, step, idx, time)
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1) = 1

    call h5screate_simple_f(0, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)
    start(1) = idx-1
    num(1) = 1
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dread_f(ID% d_id, H5T_NATIVE_DOUBLE, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)
       
  end subroutine h5md_read_obs_ds

  !> Reads a time frame for the given observable.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param step integer time step.
  !! @param time The time corresponding to the step.
  !! @private
  subroutine h5md_read_obs_d2(ID, data, step, time)
    implicit none
    type(h5md_t), intent(inout) :: ID
    double precision, intent(out) :: data(:,:)
    integer, intent(in) :: step
    double precision, intent(out) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank, idx

    rank = 3
    call h5md_get_step_time(ID, step, idx, time)
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1:2) = shape(data)

    call h5screate_simple_f(2, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)

    start(1:2) = 0
    num(1:2) = dims(1:2)
    start(3) = idx-1
    num(3) = 1
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dread_f(ID% d_id, H5T_NATIVE_DOUBLE, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)
       
  end subroutine h5md_read_obs_d2

  !> Reads a time frame for the given observable.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param step integer time step.
  !! @param time The time corresponding to the step.
  !! @private
  subroutine h5md_read_obs_d3(ID, data, step, time)
    implicit none
    type(h5md_t), intent(inout) :: ID
    double precision, intent(out) :: data(:,:,:)
    integer, intent(in) :: step
    double precision, intent(out) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank, idx

    rank = 4
    call h5md_get_step_time(ID, step, idx, time)
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1:3) = shape(data)

    call h5screate_simple_f(3, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)

    start(1:3) = 0
    num(1:3) = dims(1:3)
    start(4) = idx-1
    num(4) = 1
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dread_f(ID% d_id, H5T_NATIVE_DOUBLE, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)
       
  end subroutine h5md_read_obs_d3

  !> Reads a time frame for the given observable.
  !! @param ID h5md_t variable.
  !! @param data value of the observable.
  !! @param step integer time step.
  !! @param time The time corresponding to the step.
  !! @private
  subroutine h5md_read_obs_d4(ID, data, step, time)
    implicit none
    type(h5md_t), intent(inout) :: ID
    double precision, intent(out) :: data(:,:,:,:)
    integer, intent(in) :: step
    double precision, intent(out) :: time

    integer(HID_T) :: obs_s, mem_s
    integer(HSIZE_T), allocatable :: dims(:), max_dims(:), start(:), num(:)
    integer :: rank, idx

    rank = 5
    call h5md_get_step_time(ID, step, idx, time)
    allocate(dims(rank)) ; allocate(max_dims(rank)) ; allocate(start(rank)) ; allocate(num(rank))

    dims(1:4) = shape(data)

    call h5screate_simple_f(4, dims, mem_s, h5_error)

    call h5dget_space_f(ID% d_id , obs_s, h5_error)
    call h5sget_simple_extent_dims_f(obs_s, dims, max_dims, h5_error)

    start(1:4) = 0
    num(1:4) = dims(1:4)
    start(5) = idx-1
    num(5) = 1
    call h5sselect_hyperslab_f(obs_s, H5S_SELECT_SET_F, start, num, h5_error)
    call h5dread_f(ID% d_id, H5T_NATIVE_DOUBLE, data, num, h5_error, mem_space_id=mem_s, file_space_id=obs_s)
    call h5sclose_f(obs_s, h5_error)
    call h5sclose_f(mem_s, h5_error)

    deallocate(dims) ; deallocate(max_dims) ; deallocate(start) ; deallocate(num)
       
  end subroutine h5md_read_obs_d4



  !> Writes a parameter to the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_write_par_i1(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: data(:)

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(1)
    dims = shape(data)
    call h5screate_simple_f(1, dims, par_s, h5_error)

    call h5dcreate_f(file_id, 'parameters/'//name, H5T_NATIVE_INTEGER, par_s, par_d, h5_error)

    call h5dwrite_f(par_d, H5T_NATIVE_INTEGER, data, dims, h5_error)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_write_par_i1

  !> Writes a parameter to the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_write_par_is(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: data

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(1)
    call h5screate_f(H5S_SCALAR_F, par_s, h5_error)

    call h5dcreate_f(file_id, 'parameters/'//name, H5T_NATIVE_INTEGER, par_s, par_d, h5_error)

    call h5dwrite_f(par_d, H5T_NATIVE_INTEGER, data, dims, h5_error)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_write_par_is

  !> Writes a parameter to the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_write_par_i2(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    integer, intent(in) :: data(:,:)

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(2)
    dims = shape(data)
    call h5screate_simple_f(2, dims, par_s, h5_error)

    call h5dcreate_f(file_id, 'parameters/'//name, H5T_NATIVE_INTEGER, par_s, par_d, h5_error)

    call h5dwrite_f(par_d, H5T_NATIVE_INTEGER, data, dims, h5_error)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_write_par_i2

  !> Writes a parameter to the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_write_par_cs(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: data

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(1)
    integer(SIZE_T) :: a_size(1)
    integer(HID_T) :: a_type
    a_size(1) = len(data)
    call h5tcopy_f(H5T_NATIVE_CHARACTER, a_type, h5_error)
    call h5tset_size_f(a_type, a_size(1), h5_error)
    call h5screate_f(H5S_SCALAR_F, par_s, h5_error)

    call h5dcreate_f(file_id, 'parameters/'//name, a_type, par_s, par_d, h5_error)

    call h5dwrite_f(par_d, a_type, data, dims, h5_error)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_write_par_cs

  !> Writes a parameter to the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_write_par_d1(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    double precision, intent(in) :: data(:)

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(1)
    dims = shape(data)
    call h5screate_simple_f(1, dims, par_s, h5_error)

    call h5dcreate_f(file_id, 'parameters/'//name, H5T_NATIVE_DOUBLE, par_s, par_d, h5_error)

    call h5dwrite_f(par_d, H5T_NATIVE_DOUBLE, data, dims, h5_error)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_write_par_d1

  !> Writes a parameter to the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_write_par_ds(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    double precision, intent(in) :: data

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(1)
    call h5screate_f(H5S_SCALAR_F, par_s, h5_error)

    call h5dcreate_f(file_id, 'parameters/'//name, H5T_NATIVE_DOUBLE, par_s, par_d, h5_error)

    call h5dwrite_f(par_d, H5T_NATIVE_DOUBLE, data, dims, h5_error)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_write_par_ds

  !> Writes a parameter to the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_write_par_d2(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    double precision, intent(in) :: data(:,:)

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(2)
    dims = shape(data)
    call h5screate_simple_f(2, dims, par_s, h5_error)

    call h5dcreate_f(file_id, 'parameters/'//name, H5T_NATIVE_DOUBLE, par_s, par_d, h5_error)

    call h5dwrite_f(par_d, H5T_NATIVE_DOUBLE, data, dims, h5_error)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_write_par_d2

  !> Writes a parameter to the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_write_par_l1(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    logical, intent(in) :: data(:)

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(1)
    integer, allocatable :: data_int(:)
    allocate(data_int(size(data)))
    where (data)
        data_int = 1
    elsewhere
        data_int=0
    endwhere

    dims = shape(data)
    call h5screate_simple_f(1, dims, par_s, h5_error)

    call h5dcreate_f(file_id, 'parameters/'//name, H5T_NATIVE_INTEGER, par_s, par_d, h5_error)

    call h5dwrite_f(par_d, H5T_NATIVE_INTEGER, data_int, dims, h5_error)

    deallocate(data_int)
    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_write_par_l1

  !> Writes a parameter to the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_write_par_ls(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    logical, intent(in) :: data

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(1)
    integer :: data_int
    if (data) then
        data_int = 1
    else
        data_int = 0
    end if

    call h5screate_f(H5S_SCALAR_F, par_s, h5_error)

    call h5dcreate_f(file_id, 'parameters/'//name, H5T_NATIVE_INTEGER, par_s, par_d, h5_error)

    call h5dwrite_f(par_d, H5T_NATIVE_INTEGER, data_int, dims, h5_error)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_write_par_ls

  !> Writes a parameter to the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_write_par_l2(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    logical, intent(in) :: data(:,:)

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(2)
    integer, allocatable :: data_int(:,:)
    allocate(data_int(size(data,dim=1),size(data,dim=2)))
    where (data)
        data_int = 1
    elsewhere
        data_int=0
    endwhere

    dims = shape(data)
    call h5screate_simple_f(2, dims, par_s, h5_error)

    call h5dcreate_f(file_id, 'parameters/'//name, H5T_NATIVE_INTEGER, par_s, par_d, h5_error)

    call h5dwrite_f(par_d, H5T_NATIVE_INTEGER, data_int, dims, h5_error)

    deallocate(data_int)
    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_write_par_l2

  !> Reads a parameter from the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_read_par_i1(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    integer, intent(out) :: data(:)

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(1)
    dims = shape(data)
    call h5screate_simple_f(1, dims, par_s, h5_error)

    call h5dopen_f(file_id, 'parameters/'//name, par_d, h5_error)

    call h5dread_f(par_d, H5T_NATIVE_INTEGER, data, dims, h5_error)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_read_par_i1

  !> Reads a parameter from the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_read_par_is(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    integer, intent(out) :: data

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(1)
    call h5screate_f(H5S_SCALAR_F, par_s, h5_error)

    call h5dopen_f(file_id, 'parameters/'//name, par_d, h5_error)

    call h5dread_f(par_d, H5T_NATIVE_INTEGER, data, dims, h5_error)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_read_par_is

  !> Reads a parameter from the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_read_par_i2(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    integer, intent(out) :: data(:,:)

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(2)
    dims = shape(data)
    call h5screate_simple_f(2, dims, par_s, h5_error)

    call h5dopen_f(file_id, 'parameters/'//name, par_d, h5_error)

    call h5dread_f(par_d, H5T_NATIVE_INTEGER, data, dims, h5_error)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_read_par_i2

  !> Reads a parameter from the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_read_par_d1(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    double precision, intent(out) :: data(:)

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(1)
    dims = shape(data)
    call h5screate_simple_f(1, dims, par_s, h5_error)

    call h5dopen_f(file_id, 'parameters/'//name, par_d, h5_error)

    call h5dread_f(par_d, H5T_NATIVE_DOUBLE, data, dims, h5_error)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_read_par_d1

  !> Reads a parameter from the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_read_par_ds(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    double precision, intent(out) :: data

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(1)
    call h5screate_f(H5S_SCALAR_F, par_s, h5_error)

    call h5dopen_f(file_id, 'parameters/'//name, par_d, h5_error)

    call h5dread_f(par_d, H5T_NATIVE_DOUBLE, data, dims, h5_error)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_read_par_ds

  !> Reads a parameter from the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_read_par_d2(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    double precision, intent(out) :: data(:,:)

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(2)
    dims = shape(data)
    call h5screate_simple_f(2, dims, par_s, h5_error)

    call h5dopen_f(file_id, 'parameters/'//name, par_d, h5_error)

    call h5dread_f(par_d, H5T_NATIVE_DOUBLE, data, dims, h5_error)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_read_par_d2

  !> Reads a parameter from the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_read_par_l1(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    logical, intent(out) :: data(:)

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(1)
    integer, allocatable :: data_int(:)
    allocate(data_int(size(data)))
    where (data)
        data_int = 1
    elsewhere
        data_int=0
    endwhere

    dims = shape(data)
    call h5screate_simple_f(1, dims, par_s, h5_error)

    call h5dopen_f(file_id, 'parameters/'//name, par_d, h5_error)

    call h5dread_f(par_d, H5T_NATIVE_INTEGER, data_int, dims, h5_error)

    where (data_int .eq. 0)
        data = .false.
    elsewhere
        data = .true.
    endwhere
    deallocate(data_int)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_read_par_l1

  !> Reads a parameter from the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_read_par_ls(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    logical, intent(out) :: data

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(1)
    integer :: data_int
    if (data) then
        data_int = 1
    else
        data_int = 0
    end if

    call h5screate_f(H5S_SCALAR_F, par_s, h5_error)

    call h5dopen_f(file_id, 'parameters/'//name, par_d, h5_error)

    call h5dread_f(par_d, H5T_NATIVE_INTEGER, data_int, dims, h5_error)

    if (data_int.eq.0) then
        data = .false.
    else
        data = .true.
    end if

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_read_par_ls

  !> Reads a parameter from the parameter group.
  !! @param file_id hdf5 file ID.
  !! @param name name of the parameter.
  !! @param data value of the parameter.
  !! @private
  subroutine h5md_read_par_l2(file_id, name, data)
    integer(HID_T), intent(in) :: file_id
    character(len=*), intent(in) :: name
    logical, intent(out) :: data(:,:)

    integer(HID_T) :: par_d, par_s
    integer(HSIZE_T) :: dims(2)
    integer, allocatable :: data_int(:,:)
    allocate(data_int(size(data,dim=1),size(data,dim=2)))
    where (data)
        data_int = 1
    elsewhere
        data_int=0
    endwhere

    dims = shape(data)
    call h5screate_simple_f(2, dims, par_s, h5_error)

    call h5dopen_f(file_id, 'parameters/'//name, par_d, h5_error)

    call h5dread_f(par_d, H5T_NATIVE_INTEGER, data_int, dims, h5_error)

    where (data_int .eq. 0)
        data = .false.
    elsewhere
        data = .true.
    endwhere
    deallocate(data_int)

    call h5sclose_f(par_s, h5_error)

       
  end subroutine h5md_read_par_l2



end module h5md
