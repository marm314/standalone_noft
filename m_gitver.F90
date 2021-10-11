!!****m* DoNOF/m_gitver
!! NAME
!!  m_gitver
!!
!! FUNCTION
!!  Module to print git version information 
!! 
!! COPYRIGHT 
!! This file is distributed under the terms of the 
!! GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!! 
!! 
!! PARENTS 
!! 
!! CHILDREN 
!! 
!! 
!! SOURCE 
  
module m_gitver 
 
 implicit none
 
! private ::
!!***
 
 public :: gitversion
!!***
 
contains 
!!*** 
 
!!****f* DoNOF/gitversion
!! NAME 
!! gitversion
!!
!! FUNCTION
!!  Get the SHA string from Github version 
!! 
!! INPUTS 
!!   sha = string that on exit contains the SHA Github Key
!! 
!! OUTPUT 
!! 
!! PARENTS 
!!   
!! CHILDREN 
!! 
!! SOURCE 
 
subroutine gitversion(sha)
!Arguments ------------------------------------ 
!scalars 
!arrays 
character(100),intent(inout)::sha
!Local variables ------------------------------ 
!scalars 
!arrays 
!************************************************************************ 
 
  write(sha,'(a)')'e5da8cf748c0a981872654f5805e946e31a98261'
 
end subroutine gitversion
!!***
 
end module m_gitver
!!***
