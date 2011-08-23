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
contains
  
  !> Prints version information about the vmf90 package.
  subroutine vmf90_info()
    implicit none
    include 'vmf90_version.h'
    
    write(*,*) 'vmf90> vmf90 software'
    write(*,*) 'vmf90> (C) 2009-2011 P. de Buyl'
    write(*,*) 'vmf90> Version/date : ', trim(adjustl(git_describe))
    write(*,*) 'vmf90> git commit   : ', trim(adjustl(git_sha1))
    write(*,*) 'vmf90> git status   : ', trim(adjustl(git_status))

  end subroutine vmf90_info

end module vmf90
