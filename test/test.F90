!!****m* DoNOF/m_hubbard
!! NAME
!! m_hubbard
!!
!! FUNCTION
!!  Parameters used by Hubbard test
!!
!! COPYRIGHT
!!
!! NOTES
!!
!! SOURCE

module m_hubbard
 use m_vars
 implicit none

 integer::Nsites
 real(dp)::t,U

end module m_hubbard
!!***

!!****m* DoNOF/noft_hubbard
!! NAME
!!  noft_hubbard
!!
!! FUNCTION
!!  Testing the module using NOFTs (in perfect pairing) for N-sites Hubbard model.
!!  For PNOF7 we use PNOF7 not PNOF7s
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

program noft_hubbard
 use m_vars
 use m_hubbard
 use m_noft_driver
 implicit none
 integer::INOF,Ista=0,NBF_tot,NBF_occ,Nfrozen,Npairs
 integer::Ncoupled=1,Nbeta_elect,Nalpha_elect,iERItyp=1
 integer::imethocc=1,imethorb=1,itermax=10000,iprintdmn=0,iprintints=0
 integer::itolLambda=5,ndiis=5
 real(dp)::Enof,tolE=1.0d-9,Vnn=0.0d0
 real(dp),allocatable,dimension(:)::Occ,Work
 real(dp),allocatable,dimension(:,:)::NO_COEF,SITE_Overlap
 character(len=100)::ofile_name
 CHARACTER(len=32)::arg
 external::mo_ints
 integer::isite,isite1,lwork,info

 ! Read parameters
 if(iargc()/=4) then
  write(*,'(a)') 'Include the arguments:'
  write(*,'(a)') ' INOF(integer) Nsites(integer) Nalpha_electrons(integer) U(real)'
  stop
 endif
 call getarg(1,arg)
 read(arg,'(i4)') INOF
 call getarg(2,arg)
 read(arg,'(i4)') Nsites 
 call getarg(3,arg)
 read(arg,'(i4)') Nalpha_elect
 call getarg(4,arg)
 read(arg,'(f20.8)') U

 ! Hubbard parameteres
 t=1.0          ! We fix this one and tune U
 ! NOFT parameteres
 NBF_tot=Nsites
 NBF_occ=NBF_tot
 Nfrozen=0
 Nbeta_elect=Nalpha_elect
 Npairs=(Nbeta_elect+Nalpha_elect)/2
 ofile_name='hubbard.noft'
 ! Allocate arrays
 allocate(Occ(NBF_tot),Work(1))
 allocate(NO_COEF(NBF_tot,NBF_tot),SITE_Overlap(NBF_tot,NBF_tot))
 Occ=0.0d0; NO_COEF=0.0d0; SITE_Overlap=0.0d0;
 ! SITE_Overlap is the identity matrix
 do isite=1,NBF_tot
  SITE_Overlap(isite,isite)=1.0d0
 enddo
 ! Hcore (saved in NO_COEF initially before the diagonalization)
 NO_COEF=0.0d0
 do isite=1,NBF_tot
  do isite1=1,NBF_tot
   if(isite/=isite1) then
    NO_COEF(isite,isite1)=-t
   endif
  enddo
 enddo
 ! Use as initial GUESS for the Nat. orb. coefs. the Hcore matrix
 lwork=-1
 call DSYEV('V','L',NBF_tot,NO_COEF,NBF_tot,Occ,Work,lwork,info)
 lwork=nint(Work(1))
 if(info==0) then
  deallocate(Work)
  allocate(Work(lwork))
  call DSYEV('V','L',NBF_tot,NO_COEF,NBF_tot,Occ,Work,lwork,info)
 endif
 deallocate(Work)
 ! Initialize and run the optimization
 Occ=0.0d0
 call run_noft(INOF,Ista,NBF_tot,NBF_occ,Nfrozen,Npairs,Ncoupled,Nbeta_elect,Nalpha_elect,&
 &  iERItyp,imethocc,imethorb,itermax,iprintdmn,iprintints,itolLambda,ndiis,Enof,tolE,Vnn,&
 &  NO_COEF,SITE_Overlap,Occ,mo_ints,ofile_name)
 ! Print the optimal energy 
 write(*,'(a)') ' '
 write(*,'(a,f12.6)') 'NOFT SCF energy',Enof
 write(*,'(a)') ' '
 deallocate(Occ,NO_COEF,SITE_Overlap)

end program noft_hubbard
!!***

!!***
!!****f* DoNOF/mo_ints
!! NAME
!! mo_ints
!!
!! FUNCTION
!!  Subroutine provided by the user to compute ONEBODY and ERImol integrals from the sites basis to the 
!!  'molecular' basis NO_COEF
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

subroutine mo_ints(NBF_tot,NBF_occ,NBF_jkl,NO_COEF,ONEBODY,ERImol,ERImolv)
 use m_vars
 use m_hubbard
 implicit none
 integer,intent(in)::NBF_tot,NBF_occ,NBF_jkl
 real(dp),dimension(NBF_tot,NBF_tot),intent(inout)::ONEBODY
 real(dp),dimension(NBF_tot,NBF_jkl,NBF_jkl,NBF_jkl),optional,intent(inout)::ERImol
 real(dp),dimension(NBF_tot*NBF_jkl*NBF_jkl*NBF_jkl),optional,intent(inout)::ERImolv
 real(dp),dimension(NBF_tot,NBF_tot),intent(in)::NO_COEF
 real(dp),allocatable,dimension(:,:)::TMP_ONEBODY
 integer::isite,isite1

 ! Compute ONEBODY (initially SITE_ONEBODY, in the end ONEBODY)
 allocate(TMP_ONEBODY(NBF_tot,NBF_tot))
 ONEBODY=0.0d0
 do isite=1,NBF_tot
  do isite1=1,NBF_tot
   if(isite/=isite1) then
    ONEBODY(isite,isite1)=-t
   endif
  enddo
 enddo
 TMP_ONEBODY=matmul(ONEBODY,NO_COEF)
 ONEBODY=matmul(transpose(NO_COEF),TMP_ONEBODY)
 deallocate(TMP_ONEBODY)

 ! Compute ERImol (initially SITE_ERI, in the end ERImol)
 if(present(ERImol)) then
  ERImol=0.0d0
  do isite=1,NBF_tot
   ERImol(isite,isite,isite,isite)=U
  enddo
  call transformERI(NBF_tot,NO_COEF,ERImol) ! Site -> Nat. orbs.
 endif

end subroutine mo_ints
!!***

!!***
!!****f* DoNOF/transformERI
!! NAME
!! transformERI
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

subroutine transformERI(NBF,NO_COEF,ERImol)
 implicit none
 integer,intent(in)::NBF
 double precision,dimension(NBF,NBF,NBF,NBF),intent(inout)::ERImol
 double precision,dimension(NBF,NBF),intent(in)::NO_COEF
 double precision,dimension(:,:,:,:),allocatable::TMP_ERI
 integer::i,j,k,l,m
 allocate(TMP_ERI(NBF,NBF,NBF,NBF))
 do i=1,NBF
  do j=1,NBF
   do k=1,NBF
    do m=1,NBF
     TMP_ERI(i,j,k,m)=0.0d0
     do l=1,NBF
      TMP_ERI(i,j,k,m)=TMP_ERI(i,j,k,m)+NO_COEF(l,m)*ERImol(i,j,k,l)
     enddo
    enddo
   enddo
  enddo
  enddo
  ! K -> R
  do i=1,NBF
   do j=1,NBF
    do m=1,NBF
     do l=1,NBF
      ERImol(i,j,m,l)=0.0d0
      do k=1,NBF
      ERImol(i,j,m,l)=ERImol(i,j,m,l)+NO_COEF(k,m)*TMP_ERI(i,j,k,l)
     enddo
    enddo
   enddo
  enddo
 enddo
 ! J -> Q
 do i=1,NBF
  do m=1,NBF
   do k=1,NBF
    do l=1,NBF
     TMP_ERI(i,m,k,l)=0.0d0
     do j=1,NBF
      TMP_ERI(i,m,k,l)=TMP_ERI(i,m,k,l)+NO_COEF(j,m)*ERImol(i,j,k,l)
     enddo
    enddo
   enddo
  enddo
 enddo
 ! I -> P
 do m=1,NBF
  do j=1,NBF
   do k=1,NBF
    do l=1,NBF
     ERImol(m,j,k,l)=0.0d0
     do i=1,NBF
      ERImol(m,j,k,l)=ERImol(m,j,k,l)+NO_COEF(i,m)*TMP_ERI(i,j,k,l)
     enddo
    enddo
   enddo
  enddo
 enddo
 deallocate(TMP_ERI)

end subroutine transformERI
!!***
