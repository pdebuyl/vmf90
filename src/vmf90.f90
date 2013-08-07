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

!> This module provides information on the package version.
module vmf90
  implicit none

  private

  public :: vmf90_version, vmf90_info

  include 'vmf90_version.h'

contains
  
  !> Prints version information about the vmf90 package.
  subroutine vmf90_info()
    implicit none
    
    write(*,*) 'vmf90> vmf90 software'
    write(*,*) 'vmf90> (C) 2009-2013 P. de Buyl'
    write(*,*) 'vmf90> Version    : ', trim(adjustl(git_describe))
    write(*,*) 'vmf90> Date       : ', trim(adjustl(git_date))
    write(*,*) 'vmf90> git commit : ', trim(adjustl(git_sha1))
    write(*,*) 'vmf90> git status : ', trim(adjustl(git_status))

  end subroutine vmf90_info

  !> Returns version information about the vmf90 package.
  function vmf90_version()
    implicit none
    character(len=80) :: vmf90_version

    if ( trim(adjustl(git_status)).eq.'Clean' ) then
       vmf90_version = trim(adjustl(git_describe))//' SHA1 '//trim(adjustl(git_sha1))
    else
       vmf90_version = 'Unclean based on '//trim(adjustl(git_describe))
    end if

  end function vmf90_version

end module vmf90
