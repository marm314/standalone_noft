!!****m* DoNOF/m_fcidump
!! NAME
!! m_fcidump
!!
!! FUNCTION
!!  Parameters used by Hubbard test
!!
!! COPYRIGHT
!!
!! NOTES
!!
!! SOURCE

module m_fcidump
 use m_definitions
 implicit none

 real(dp),allocatable,dimension(:,:)::hCORE_IN
 real(dp),allocatable,dimension(:,:,:,:)::ERI_IN

end module m_fcidump
!!***

!!****m* DoNOF/noft_fcidump
!! NAME
!!  noft_fcidump
!!
!! FUNCTION
!!  Testing the module using NOFTs (in perfect pairing) for N-MO Hubbard model.
!!  For PNOF7 we use PNOF7 not PNOF7s
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

program noft_fcidump
 use m_definitions
 use m_fcidump
 use m_noft_driver
 implicit none
 logical::restart=.false.
 logical::nbf_FCIDUMP_found,not_orb_opt
 integer::INOF,Ista=0,NBF_tot,NBF_occ,Nfrozen,Npairs
 integer::Ncoupled=1,Nbeta_elect,Nalpha_elect
 integer::imethocc=1,imethorb=1,itermax=10000,iprintdmn=0,iprintswdmn=0,iprintints=0
 integer::itolLambda=5,ndiis=5,iguess=0,irestart=0,iskip=0,freeze_orb=0
 real(dp)::Enof,tolE=tol9,Vnn=zero,Val
 real(dp),allocatable,dimension(:)::Occ,Work
 real(dp),allocatable,dimension(:,:)::NO_COEF,MO_Overlap,DM1
 character(len=100)::ofile_name
 character(len=100)::arg
 external::mo_ints
 integer::ibasis,ibasis1,ibasis2,ibasis3,lwork,info
 integer::nbf_unit=25

 ! Read parameters
 if(iargc()/=4 .and. iargc()/=5 .and. iargc()/=6) then
  write(*,'(a)') 'Include the arguments:'
  write(*,'(a)') ' INOF(integer) Nbasis(integer) Nalpha_electrons(integer) iguess(0 Hcore, 1 init MO basis)'
  write(*,'(a)') ' or '
  write(*,'(a)') ' INOF(integer) Nbasis(integer) Nalpha_electrons(integer) iguess(0 Hcore, 1 init MO basis) restart(integer)'
  write(*,'(a)') ' or '
  write(*,'(a)') ' INOF(integer) Nbasis(integer) Nalpha_electrons(integer) iguess(0 Hcore, 1 init MO basis) restart(integer) & 
                   skip_lines(integer)'
  stop
 endif
 call getarg(1,arg)
 read(arg,'(i4)') INOF
 call getarg(2,arg)
 read(arg,'(i4)') NBF_tot
 call getarg(3,arg)
 read(arg,'(i4)') Nalpha_elect
 call getarg(4,arg)
 read(arg,'(i4)') iguess
 restart=.false.
 if(iargc()==5 .or. iargc()==6) then
  call getarg(5,arg)
  read(arg,'(i4)') irestart
  if(irestart==1) restart=.true.
 endif
 if(iargc()==6) then
  call getarg(6,arg)
  read(arg,'(i4)') iskip
 endif

 ! NOFT parameteres
 NBF_occ=NBF_tot
 Nfrozen=0
 Nbeta_elect=Nalpha_elect
 Npairs=(Nbeta_elect+Nalpha_elect)/2
 Ncoupled=(NBF_tot/Npairs)-1
 ofile_name='fcidump.noft'
 ! Allocate arrays
 allocate(Occ(NBF_tot),Work(1))
 allocate(NO_COEF(NBF_tot,NBF_tot),MO_Overlap(NBF_tot,NBF_tot),DM1(NBF_tot,NBF_tot))
 Occ=zero; NO_COEF=zero; MO_Overlap=zero; DM1=zero;
 ! MO_Overlap is the identity matrix
 do ibasis=1,NBF_tot
  MO_Overlap(ibasis,ibasis)=one
 enddo
 ! Read and store FCIDUMP info
 allocate(hCORE_IN(NBF_tot,NBF_tot),ERI_IN(NBF_tot,NBF_tot,NBF_tot,NBF_tot))
 hCORE_IN=0d0;ERI_IN=0d0;
 inquire(file='NOT_ORB',exist=not_orb_opt)
 if(not_orb_opt) freeze_orb=1
 inquire(file='FCIDUMP',exist=nbf_FCIDUMP_found)
 if(nbf_FCIDUMP_found) then
  open(newunit=nbf_unit,file='FCIDUMP',status='old',form='formatted')
  if(iskip/=0) then
   do ibasis=1,iskip
    read(nbf_unit,'(a)') arg 
   enddo 
  endif
  do
   if(iskip==0) then
    read(nbf_unit,*) ibasis,ibasis1,ibasis2,ibasis3,Val
   else
    read(nbf_unit,*) Val,ibasis,ibasis1,ibasis2,ibasis3
   endif
   if(ibasis==-1 .and. ibasis1==-1 .and. ibasis2==-1 .and. ibasis3==-1) exit
   if(ibasis==0 .and. ibasis1==0 .and. ibasis2==0 .and. ibasis3==0) Vnn=Val
   if(ibasis/=0 .and. ibasis1/=0 .and. ibasis2==0 .and. ibasis3==0) then

    hCORE_IN(ibasis,ibasis1)=Val
    hCORE_IN(ibasis1,ibasis)=Val

   endif
   if(ibasis/=0 .and. ibasis1/=0 .and. ibasis2/=0 .and. ibasis3/=0) then

    ERI_IN(ibasis ,ibasis2,ibasis1,ibasis3)=Val
    ERI_IN(ibasis ,ibasis3,ibasis1,ibasis2)=Val
    ERI_IN(ibasis1,ibasis3,ibasis ,ibasis2)=Val
    ERI_IN(ibasis1,ibasis2,ibasis ,ibasis3)=Val

    ERI_IN(ibasis2,ibasis ,ibasis3,ibasis1)=Val
    ERI_IN(ibasis3,ibasis ,ibasis2,ibasis1)=Val
    ERI_IN(ibasis3,ibasis1,ibasis2,ibasis )=Val
    ERI_IN(ibasis2,ibasis1,ibasis3,ibasis )=Val

   endif
   
  enddo
  close(nbf_unit)
 else
  write(*,*) 'FCIDUMP file not found!'
  stop
 endif


 if(iguess==0) then ! Hcore (saved in NO_COEF initially before the diagonalization)  
     
  ! Use as initial GUESS for the MO that diagonalize the Hcore matrix
  NO_COEF=hCORE_IN
  lwork=-1
  call DSYEV('V','L',NBF_tot,NO_COEF,NBF_tot,Occ,Work,lwork,info)
  lwork=nint(Work(1))
  if(info==0) then
   deallocate(Work)
   allocate(Work(lwork))
   call DSYEV('V','L',NBF_tot,NO_COEF,NBF_tot,Occ,Work,lwork,info)
  endif
  deallocate(Work)
 else ! MO basis (identity) 
  do ibasis=1,NBF_tot
   NO_COEF(ibasis,ibasis)=one
  enddo
 endif
 ! Initialize and run the optimization
 Occ=zero
 if(.not.restart) then
   write(*,*) 'running NOFT module'
   call run_noft(INOF,Ista,NBF_tot,NBF_occ,Nfrozen,Npairs,Ncoupled,Nbeta_elect,Nalpha_elect,&
  &   imethocc,imethorb,itermax,iprintdmn,iprintswdmn,iprintints,itolLambda,ndiis,&
  &   Enof,tolE,Vnn,MO_Overlap,Occ,mo_ints,ofile_name,NO_COEF=NO_COEF,iNOTupdateORB=freeze_orb)
 else
   write(*,*) 'running NOFT module (restart)'
   call run_noft(INOF,Ista,NBF_tot,NBF_occ,Nfrozen,Npairs,Ncoupled,Nbeta_elect,Nalpha_elect,&
  &   imethocc,imethorb,itermax,iprintdmn,iprintswdmn,iprintints,itolLambda,ndiis,&
  &   Enof,tolE,Vnn,MO_Overlap,Occ,mo_ints,ofile_name,NO_COEF=NO_COEF,restart=.true.,&
  &   ireadGAMMAS=1,ireadOCC=1,ireadCOEF=1,ireadFdiag=1,iNOTupdateOCC=0,iNOTupdateORB=freeze_orb)
 endif
 ! Print the optimal energy 
 write(*,'(a)') ' '
 write(*,'(a,f20.8)') 'NOFT SCF energy',Enof
 write(*,'(a)') ' '
 write(*,'(a)') ' '
 write(*,*) 'Final occupation numbers'
 do ibasis=1,NBF_tot
  write(*,'(*(f15.5))') Occ(ibasis)
 enddo
 write(*,'(a)') ' '
 write(*,*) 'Final density matrix on the initial basis'
 write(*,'(a)') ' '
 do ibasis=1,NBF_tot
  DM1(ibasis,ibasis)=0.5d0*Occ(ibasis)
 enddo
 DM1=matmul(matmul(NO_COEF,DM1),transpose(NO_COEF))
 do ibasis=1,NBF_tot
  write(*,'(*(f15.5))') DM1(ibasis,:)
 enddo
 deallocate(Occ,NO_COEF,MO_Overlap,DM1)
 deallocate(hCORE_IN,ERI_IN)

end program noft_fcidump
!!***

!!***
!!****f* DoNOF/mo_ints
!! NAME
!! mo_ints
!!
!! FUNCTION
!!  Subroutine provided by the user to compute hCORE and ERImol integrals from the MO basis to the 
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

subroutine mo_ints(NBF_tot, NBF_occ, NBF_jkl, Occ, DM2_JK, NO_COEF, hCORE, ERImol, ERImolJsr, ERImolLsr, &
     &             NO_COEF_cmplx, hCORE_cmplx, ERImol_cmplx, ERImolJsr_cmplx, ERImolLsr_cmplx, all_ERIs, &
     &             Edft_xc, do_xc_dft)
 use m_definitions
 use m_fcidump
 implicit none
 logical, optional, intent(in)     :: all_ERIs, do_xc_dft
 integer, intent(in)              :: NBF_tot, NBF_occ, NBF_jkl
 real(dp), intent(in)             :: Occ(NBF_occ)
 real(dp), optional, intent(inout) :: Edft_xc
 real(dp), optional, intent(in)    :: DM2_JK(2, NBF_occ, NBF_occ)
 real(dp), optional, intent(in)    :: NO_COEF(NBF_tot, NBF_tot)
 real(dp), optional, intent(inout) :: hCORE(NBF_tot, NBF_tot)
 real(dp), optional, intent(inout) :: ERImol(NBF_tot, NBF_jkl, NBF_jkl, NBF_jkl)
 real(dp), optional, intent(inout) :: ERImolJsr(NBF_tot, NBF_jkl, NBF_jkl)
 real(dp), optional, intent(inout) :: ERImolLsr(NBF_tot, NBF_jkl, NBF_jkl)
 complex(dp), optional, intent(in)    :: NO_COEF_cmplx(NBF_tot, NBF_tot)
 complex(dp), optional, intent(inout) :: hCORE_cmplx(NBF_tot, NBF_tot)
 complex(dp), optional, intent(inout) :: ERImol_cmplx(NBF_tot, NBF_jkl, NBF_jkl, NBF_jkl)
 complex(dp), optional, intent(inout) :: ERImolJsr_cmplx(NBF_tot, NBF_jkl, NBF_jkl)
 complex(dp), optional, intent(inout) :: ERImolLsr_cmplx(NBF_tot, NBF_jkl, NBF_jkl)
 real(dp),allocatable,dimension(:,:)::TMP_hCORE

 !write(*,*) ' Starting transformation of hCORE and ERI integrals'
 ! Compute hCORE (initially SITE_hCORE, in the end hCORE)
 allocate(TMP_hCORE(NBF_tot,NBF_tot))
 TMP_hCORE=matmul(hCORE_IN,NO_COEF)
 hCORE=matmul(transpose(NO_COEF),TMP_hCORE)
 deallocate(TMP_hCORE)

 ! Compute ERImol (initially SITE_ERI, in the end ERImol)
 ERImol=ERI_IN
 call transformERI(NBF_tot,NO_COEF,ERImol) ! MO -> New NOs.

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
