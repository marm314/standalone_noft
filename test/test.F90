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
 use m_definitions
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
 use m_definitions
 use m_hubbard
 use m_noft_driver
 implicit none
 logical::restart=.false.
 integer::INOF,Ista=0,NBF_tot,NBF_occ,Nfrozen,Npairs
 integer::Ncoupled=1,Nbeta_elect,Nalpha_elect
 integer::imethocc=1,imethorb=1,itermax=10000,iprintdmn=0,iprintswdmn=0,iprintints=0
 integer::itolLambda=5,ndiis=5,iguess=0
 real(dp)::Enof,tolE=tol9,Vnn=zero
 real(dp),allocatable,dimension(:)::Occ,Work
 real(dp),allocatable,dimension(:,:)::NO_COEF,SITE_Overlap
 character(len=100)::ofile_name
 CHARACTER(len=32)::arg
 external::mo_ints
 integer::isite,isite1,lwork,info

 ! Read parameters
 if(iargc()/=5 .and. iargc()/=6) then
  write(*,'(a)') 'Include the arguments:'
  write(*,'(a)') ' INOF(integer) Nsites(integer) Nalpha_electrons(integer) U(real) iguess(0 Hcore, 1 site basis)'
  write(*,'(a)') ' or'
  write(*,'(a)') ' INOF(integer) Nsites(integer) Nalpha_electrons(integer) U(real) iguess(0 Hcore, 1 site basis) R(restart)'
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
 call getarg(5,arg)
 read(arg,'(i4)') iguess
 restart=.false.
 if(iargc()==6) restart=.true.

 ! Hubbard parameteres
 t=one          ! We fix this one and tune U
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
 Occ=zero; NO_COEF=zero; SITE_Overlap=zero;
 ! SITE_Overlap is the identity matrix
 do isite=1,NBF_tot
  SITE_Overlap(isite,isite)=one
 enddo
 if(iguess==0) then ! Hcore (saved in NO_COEF initially before the diagonalization)  
  NO_COEF=zero
  do isite=1,NBF_tot-1
   isite1=isite+1
   NO_COEF(isite,isite1)=-t
   NO_COEF(isite1,isite)=-t
  enddo
  NO_COEF(1,NBF_tot)=-t
  NO_COEF(NBF_tot,1)=-t
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
 else ! Site basis (identity) 
  do isite=1,NBF_tot
   NO_COEF(isite,isite)=one
  enddo
 endif
 ! Initialize and run the optimization
 Occ=zero
 if(.not.restart) then
  write(*,*) 'running NOFT module'
  call run_noft(INOF,Ista,NBF_tot,NBF_occ,Nfrozen,Npairs,Ncoupled,Nbeta_elect,Nalpha_elect,&
 &   imethocc,imethorb,itermax,iprintdmn,iprintswdmn,iprintints,itolLambda,ndiis,&
 &   Enof,tolE,Vnn,SITE_Overlap,Occ,mo_ints,ofile_name,NO_COEF=NO_COEF)
 else
  write(*,*) 'running NOFT module with restart'
  call run_noft(INOF,Ista,NBF_tot,NBF_occ,Nfrozen,Npairs,Ncoupled,Nbeta_elect,Nalpha_elect,&
  & imethocc,imethorb,itermax,0,0,0,itolLambda,ndiis,Enof,tolE,Vnn,SITE_overlap,Occ,&
  & mo_ints,ofile_name,NO_COEF=NO_COEF,restart=.true.,ireadGAMMAS=1,ireadOCC=1,ireadCOEF=1,&
  & ireadFdiag=1,iNOTupdateOCC=0,iNOTupdateORB=0)
 endif
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
!!  Subroutine provided by the user to compute hCORE and ERImol integrals from the sites basis to the 
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

subroutine mo_ints(NBF_tot,NBF_occ,NBF_jkl,Occ,NO_COEF,hCORE,ERImol)
 use m_definitions
 use m_hubbard
 implicit none
 integer,intent(in)::NBF_tot,NBF_occ,NBF_jkl
 real(dp),dimension(NBF_occ),intent(in)::Occ
 real(dp),optional,dimension(NBF_tot,NBF_tot),intent(inout)::hCORE
 real(dp),optional,dimension(NBF_tot,NBF_jkl,NBF_jkl,NBF_jkl),intent(inout)::ERImol
 real(dp),optional,dimension(NBF_tot,NBF_tot),intent(in)::NO_COEF
 real(dp),allocatable,dimension(:,:)::TMP_hCORE
 integer::isite,isite1

 !write(*,*) ' Starting transformation of hCORE and ERI integrals'
 ! Compute hCORE (initially SITE_hCORE, in the end hCORE)
 allocate(TMP_hCORE(NBF_tot,NBF_tot))
 hCORE=zero
 do isite=1,NBF_tot-1
  isite1=isite+1
  hCORE(isite,isite1)=-t
  hCORE(isite1,isite)=-t
 enddo
 hCORE(1,NBF_tot)=-t
 hCORE(NBF_tot,1)=-t
 TMP_hCORE=matmul(hCORE,NO_COEF)
 hCORE=matmul(transpose(NO_COEF),TMP_hCORE)
 deallocate(TMP_hCORE)

 ! Compute ERImol (initially SITE_ERI, in the end ERImol)
 ERImol=zero
 do isite=1,NBF_tot
  ERImol(isite,isite,isite,isite)=U
 enddo
 call transformERI(NBF_tot,NO_COEF,ERImol) ! Site -> Nat. orbs.

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

end subroutine transformERI
!!***
