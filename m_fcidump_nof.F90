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

 integer::iskip_nof=0
 logical::fort_fcidump=.false.,standard_fcidump=.true.
 real(dp)::Vnn_nof
 real(dp),allocatable,dimension(:,:)::hCORE_IN_NOF
 real(dp),allocatable,dimension(:,:,:,:)::ERI_IN_NOF

 public :: transformERI_NOF,read_fcidump_NOF
!!***

contains

!!***
!!****f* DoNOF/transformERI_NOF
!! NAME
!! transformERI_NOF
!!
!! FUNCTION
!!  Transform the ERI from MO to new MO basis
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

!!***
!!****f* DoNOF/read_fcidump_NOF
!! NAME
!! read_fcidump_NOF
!!
!! FUNCTION
!!  Read the FCIDUMP file
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

subroutine read_fcidump_NOF()
 use m_definitions
 implicit none
 logical::nbf_FCIDUMP_found
 integer::iorb,jorb,korb,lorb
 integer::nbf_unit=25
 real(dp)::Val
 character(len=100)::arg

 inquire(file='FCIDUMP',exist=nbf_FCIDUMP_found)
 if(nbf_FCIDUMP_found) then
  open(newunit=nbf_unit,file='FCIDUMP',status='old',form='formatted')
  if(iskip_nof/=0) then
   do iorb=1,iskip_nof
    read(nbf_unit,'(a)') arg
   enddo 
  endif 
  do
   if(standard_fcidump) then ! Standard FCIDUMP
    read(nbf_unit,*) Val,iorb,jorb,korb,lorb
   else ! CMZ style
    read(nbf_unit,*) iorb,jorb,korb,lorb,Val
   endif
   if(iorb==-1 .and. jorb==-1 .and. korb==-1 .and. lorb==-1) exit
   if(iorb==0 .and. jorb==0 .and. korb==0 .and. lorb==0) Vnn_nof=Val
   if(iorb/=0 .and. jorb/=0 .and. korb==0 .and. lorb==0) then

    hCORE_IN_NOF(iorb,jorb)=Val
    hCORE_IN_NOF(jorb,iorb)=Val

   endif
   if(iorb/=0 .and. jorb/=0 .and. korb/=0 .and. lorb/=0) then

    ERI_IN_NOF(iorb,korb,jorb,lorb)=Val
    ERI_IN_NOF(iorb,lorb,jorb,korb)=Val
    ERI_IN_NOF(jorb,lorb,iorb,korb)=Val
    ERI_IN_NOF(jorb,korb,iorb,lorb)=Val

    ERI_IN_NOF(korb,iorb,lorb,jorb)=Val
    ERI_IN_NOF(lorb,iorb,korb,jorb)=Val
    ERI_IN_NOF(lorb,jorb,korb,iorb)=Val
    ERI_IN_NOF(korb,jorb,lorb,iorb)=Val

   endif

  enddo
  close(nbf_unit)
 else
  write(*,*) 'FCIDUMP file not found!'
  stop
 endif

end subroutine read_fcidump_NOF
!!***

end module m_fcidump_nof
!!***
