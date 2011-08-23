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
