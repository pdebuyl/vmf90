!> The following fortran 90 module helps parsing small text files for use
!! as input in fortran programs.
!!
!! note : fixed length : length of a line -> elements et tempchar
!!                     : length of the PTread_s function

module ParseText

  !> A defined type that contains all the info from the input file
  type PTo
     private
     !> the filename
     character(len=128) :: filename
     !> the number of lines in the file, will be computed at parsing
     integer :: nlines
     !> the character used as the equality relation, defaults to '='
     character(len=8) :: equalsign 
     !> the character used as the comment identifier, defaults to '#'
     character(len=8) :: comment
     !> the array containing the lines of the file; the maximum length
     !! of a line is 144 (can be modified); elements will be allocated
     !! at parsing
     character(len=144), allocatable :: elements(:)
  end type PTo

contains

  !> Reads a file into a PTo variable.
  !!
  !! @param PTin A PTo variable.
  !! @param filename The file name.
  !! @param fileunit A Fortran file unit.
  !! @param equalsign Optional, allows to redefine an assignment symbol.
  !! @param comment Optional, allows to redefine the comment symbol.
  subroutine PTparse(PTin, filename, fileunit, equalsign, comment)
    type(PTo), intent(out) :: PTin
    character(len=*), intent(in) :: filename
    integer, intent(in), optional :: fileunit
    character(len=*), intent(in), optional :: equalsign, comment

    integer :: i, iostat
    character(len=144) :: tempchar

    PTin%filename = filename
    if (present(fileunit)) then
       write(*,*) ''
    else
       stop 'please give a fileunit to PTparse'
    end if

    if (present(equalsign)) then
       PTin%equalsign = equalsign
    else
       PTin%equalsign = '='
    end if
    if (present(comment)) then
       PTin%comment = comment
    else
       PTin%comment = '#'
    end if
           

    open(fileunit,file=filename)
    PTin%nlines = 0
    count_loop : do
       read(fileunit, '(A)', iostat=iostat) tempchar
       if (iostat.lt.0) exit count_loop
       if (iostat.gt.0) cycle
       
       if (index(tempchar, PTin%comment).ne.1) then
          PTin%nlines = PTin%nlines + 1
       end if
    end do count_loop
    close(fileunit)

    if (PTin%nlines.gt.0) then
       allocate(PTin%elements(PTin%nlines))
    else
       write(*,*) 'filename = ', filename
       stop 'No line to parse in PTparse'
    end if

    i = 1
    open(fileunit,file=filename)
    read_loop : do
       read(fileunit, '(A)', iostat=iostat) tempchar
       if (iostat.lt.0) exit read_loop
       if (iostat.gt.0) cycle
       
       if (index(tempchar, PTin%comment).ne.1) then
          PTin%elements(i) = adjustl(tempchar)
          i = i + 1
       end if
    end do read_loop
    close(fileunit)

  end subroutine PTparse

  !> Prints all the content of the PTo object from memory.
  !!
  !! @param PTin A PTo variable.
  subroutine PTprint(PTin)
    type(PTo), intent(in) :: PTin
    
    integer :: i

    if (allocated(PTin%elements)) then
       do i=1, PTin%nlines
          write(*,'(A)') PTin%elements(i)
       end do
    else
       write(*,*) 'nothing to print in PTprint'
    end if
       
  end subroutine PTprint

  !> Resets a PTo object.
  !!
  !! @param PTin A PTo variable.
  subroutine PTkill(PTin)
    type(PTo), intent(inout) :: PTin

    if (allocated(PTin%elements)) then
       deallocate(PTin%elements)
       PTin%nlines = 0
    else
       write(*,*) 'nothing to kill in PTkill'
    end if
  end subroutine PTkill

  !> Reads a double precision value.
  !!
  !! @param PTin A PTo variable.
  !! @param var_name The name to look up.
  double precision function PTread_d(PTin, var_name)
    type(PTo), intent(in) :: PTin
    character(len=*), intent(in) :: var_name

    integer :: i
    character(len=144) :: tempchar
    double precision :: value
    logical :: found

    found = .false.
    do i=1,PTin%nlines
       if (index(PTin%elements(i),var_name).eq.1) then
          tempchar = PTin%elements(i)(len(trim(var_name))+1:)
          tempchar = adjustl(tempchar)
          if (index(tempchar,trim(PTin%equalsign)).eq.1) then
             tempchar = adjustl(tempchar)
             tempchar = tempchar(len(trim(PTin%equalsign))+1:)
             read(tempchar,*) value
             found = .true.
             exit
          end if
       end if
    end do

    if (found) then
       PTread_d = value
    else
       write(*,*) 'variable ',var_name,' not found in the file ',PTin%filename
       stop
    end if

  end function PTread_d

  !> Reads a real value.
  !!
  !! @param PTin A PTo variable.
  !! @param var_name The name to look up.
  real function PTread_r(PTin, var_name)
    type(PTo), intent(in) :: PTin
    character(len=*), intent(in) :: var_name

    integer :: i
    character(len=144) :: tempchar
    real :: value
    logical :: found

    found = .false.
    do i=1,PTin%nlines
       if (index(PTin%elements(i),var_name).eq.1) then
          tempchar = PTin%elements(i)(len(trim(var_name))+1:)
          tempchar = adjustl(tempchar)
          if (index(tempchar,trim(PTin%equalsign)).eq.1) then
             tempchar = adjustl(tempchar)
             tempchar = tempchar(len(trim(PTin%equalsign))+1:)
             read(tempchar,*) value
             found = .true.
             exit
          end if
       end if
    end do

    if (found) then
    PTread_r = value
    else
       write(*,*) 'variable ',var_name,' not found in the file ',PTin%filename
       stop
    end if

  end function PTread_r


  !> Reads an integer value.
  !!
  !! @param PTin A PTo variable.
  !! @param var_name The name to look up.
  integer function PTread_i(PTin, var_name)
    type(PTo), intent(in) :: PTin
    character(len=*), intent(in) :: var_name

    integer :: i
    character(len=144) :: tempchar
    integer :: value
    logical :: found

    found = .false.
    do i=1,PTin%nlines
       if (index(PTin%elements(i),var_name).eq.1) then
          tempchar = PTin%elements(i)(len(trim(var_name))+1:)
          tempchar = adjustl(tempchar)
          if (index(tempchar,trim(PTin%equalsign)).eq.1) then
             tempchar = adjustl(tempchar)
             tempchar = tempchar(len(trim(PTin%equalsign))+1:)
             read(tempchar,*) value
             found = .true.
             exit
          end if
       end if
    end do

    if (found) then
    PTread_i = value
    else
       write(*,*) 'variable ',var_name,' not found in the file ',PTin%filename
       stop
    end if

  end function PTread_i

  !> Reads a string value.
  !!
  !! @param PTin A PTo variable.
  !! @param var_name The name to look up.
  character*144 function PTread_s(PTin, var_name)
    type(PTo), intent(in) :: PTin
    character(len=*), intent(in) :: var_name

    integer :: i
    character(len=144) :: tempchar
    character(len=144) :: value
    logical :: found

    found = .false.
    do i=1,PTin%nlines
       if (index(PTin%elements(i),var_name).eq.1) then
          tempchar = PTin%elements(i)(len(trim(var_name))+1:)
          tempchar = adjustl(tempchar)
          if (index(tempchar,trim(PTin%equalsign)).eq.1) then
             tempchar = adjustl(tempchar)
             tempchar = tempchar(len(trim(PTin%equalsign))+1:)
             read(tempchar,*) value
             found = .true.
             exit
          end if
       end if
    end do

    if (found) then
    PTread_s = value
    else
       write(*,*) 'variable ',var_name,' not found in the file ',PTin%filename
       stop
    end if

  end function PTread_s

  !> Reads a logical value.
  !!
  !! @param PTin A PTo variable.
  !! @param var_name The name to look up.
  logical function PTread_l(PTin, var_name)
    type(PTo), intent(in) :: PTin
    character(len=*), intent(in) :: var_name

    integer :: i
    character(len=144) :: tempchar
    logical :: value
    logical :: found

    found = .false.
    do i=1,PTin%nlines
       if (index(PTin%elements(i),var_name).eq.1) then
          tempchar = PTin%elements(i)(len(trim(var_name))+1:)
          tempchar = adjustl(tempchar)
          if (index(tempchar,trim(PTin%equalsign)).eq.1) then
             tempchar = adjustl(tempchar)
             tempchar = tempchar(len(trim(PTin%equalsign))+1:)
             read(tempchar,*) value
             found = .true.
             exit
          end if
       end if
    end do

    if (found) then
    PTread_l = value
    else
       write(*,*) 'variable ',var_name,' not found in the file ',PTin%filename
       stop
    end if

  end function PTread_l

  !> Reads a vector of integers.
  !!
  !! @param PTin A PTo variable.
  !! @param var_name The name to look up.
  !! @param n The length of the vector.
  integer function PTread_ivec(PTin, var_name, n)
    type(PTo), intent(in) :: PTin
    character(len=*), intent(in) :: var_name
    integer, intent(in) :: n
    dimension :: PTread_ivec(n)

    integer :: i
    character(len=144) :: tempchar
    integer, dimension(n) :: value
    logical :: found

    found = .false.
    do i=1,PTin%nlines
       if (index(PTin%elements(i),var_name).eq.1) then
          tempchar = PTin%elements(i)(len(trim(var_name))+1:)
          tempchar = adjustl(tempchar)
          if (index(tempchar,trim(PTin%equalsign)).eq.1) then
             tempchar = adjustl(tempchar)
             tempchar = tempchar(len(trim(PTin%equalsign))+1:)
             read(tempchar,*) value
             found = .true.
             exit
          end if
       end if
    end do

    if (found) then
    PTread_ivec = value
    else
       write(*,*) 'variable ',var_name,' not found in the file ',PTin%filename
       stop
    end if

  end function PTread_ivec

  !> Reads a vector of double precision.
  !!
  !! @param PTin A PTo variable.
  !! @param var_name The name to look up.
  !! @param n The length of the vector.
  double precision function PTread_dvec(PTin, var_name, n)
    type(PTo), intent(in) :: PTin
    character(len=*), intent(in) :: var_name
    integer, intent(in) :: n
    dimension :: PTread_dvec(n)

    integer :: i
    character(len=144) :: tempchar
    double precision, dimension(n) :: value
    logical :: found

    found = .false.
    do i=1,PTin%nlines
       if (index(PTin%elements(i),var_name).eq.1) then
          tempchar = PTin%elements(i)(len(trim(var_name))+1:)
          tempchar = adjustl(tempchar)
          if (index(tempchar,trim(PTin%equalsign)).eq.1) then
             tempchar = adjustl(tempchar)
             tempchar = tempchar(len(trim(PTin%equalsign))+1:)
             read(tempchar,*) value
             found = .true.
             exit
          end if
       end if
    end do

    if (found) then
    PTread_dvec = value
    else
       write(*,*) 'variable ',var_name,' not found in the file ',PTin%filename
       stop
    end if

  end function PTread_dvec

  !> Reads a vector of logical.
  !!
  !! @param PTin A PTo variable.
  !! @param var_name The name to look up.
  !! @param n The length of the vector.
  logical function PTread_lvec(PTin, var_name, n)
    type(PTo), intent(in) :: PTin
    character(len=*), intent(in) :: var_name
    integer, intent(in) :: n
    dimension :: PTread_lvec(n)

    integer :: i
    character(len=144) :: tempchar
    logical, dimension(n) :: value
    logical :: found

    found = .false.
    do i=1,PTin%nlines
       if (index(PTin%elements(i),var_name).eq.1) then
          tempchar = PTin%elements(i)(len(trim(var_name))+1:)
          tempchar = adjustl(tempchar)
          if (index(tempchar,trim(PTin%equalsign)).eq.1) then
             tempchar = adjustl(tempchar)
             tempchar = tempchar(len(trim(PTin%equalsign))+1:)
             read(tempchar,*) value
             found = .true.
             exit
          end if
       end if
    end do

    if (found) then
    PTread_lvec = value
    else
       write(*,*) 'variable ',var_name,' not found in the file ',PTin%filename
       stop
    end if

  end function PTread_lvec

end module ParseText



