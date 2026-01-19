!!****m* DoNOF/m_fcidump_nof
!! NAME
!! m_fcidump_nof
!!
!! FUNCTION
!!  Parameters and arrays used when reading a FCIDUMP file
!!
!! COPYRIGHT
!!
!! NOTES
!!
!! SOURCE

module m_fcidump_nof
 use m_definitions
 implicit none

 logical::fort_fcidump=.false.
 real(dp),allocatable,dimension(:,:)::hCORE_IN_NOF
 real(dp),allocatable,dimension(:,:,:,:)::ERI_IN_NOF

 public :: transformERI_NOF
!!***

contains

!!***
!!****f* DoNOF/transformERI_NOF
!! NAME
!! transformERI_NOF
!!
!! FUNCTION
!!  Transform the ERI from SITE basis to the 'molecular basis'
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine transformERI_NOF(NBF,NO_COEF,ERImol)
 use m_definitions
 implicit none
 integer,intent(in)::NBF
 double precision,dimension(NBF,NBF,NBF,NBF),intent(inout)::ERImol
 double precision,dimension(NBF,NBF),intent(in)::NO_COEF
 double precision,dimension(:,:,:,:),allocatable::TMP_ERI
 integer::i,j,k,l,m
 allocate(TMP_ERI(NBF,NBF,NBF,NBF))
 ! L -> S
 TMP_ERI=zero
 do i=1,NBF
  do j=1,NBF
   do k=1,NBF
    do m=1,NBF
     do l=1,NBF
      TMP_ERI(i,j,k,m)=TMP_ERI(i,j,k,m)+NO_COEF(l,m)*ERImol(i,j,k,l)
     enddo
    enddo
   enddo
  enddo
 enddo
 ! K -> R
 ERImol=zero
 do i=1,NBF
  do j=1,NBF
   do m=1,NBF
    do l=1,NBF
     do k=1,NBF
      ERImol(i,j,m,l)=ERImol(i,j,m,l)+NO_COEF(k,m)*TMP_ERI(i,j,k,l)
     enddo
    enddo
   enddo
  enddo
 enddo
 ! J -> Q
 TMP_ERI=zero
 do i=1,NBF
  do m=1,NBF
   do k=1,NBF
    do l=1,NBF
     do j=1,NBF
      TMP_ERI(i,m,k,l)=TMP_ERI(i,m,k,l)+NO_COEF(j,m)*ERImol(i,j,k,l)
     enddo
    enddo
   enddo
  enddo
 enddo
 ! I -> P
 ERImol=zero
 do m=1,NBF
  do j=1,NBF
   do k=1,NBF
    do l=1,NBF
     do i=1,NBF
      ERImol(m,j,k,l)=ERImol(m,j,k,l)+NO_COEF(i,m)*TMP_ERI(i,j,k,l)
     enddo
    enddo
   enddo
  enddo
 enddo
 deallocate(TMP_ERI)

end subroutine transformERI_NOF
!!***

end module m_fcidump_nof
!!***

