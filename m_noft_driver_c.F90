!!****m* DoNOF/m_fcidump_c
!! NAME
!! m_fcidump_c
!!
!! FUNCTION
!!  FCIDUMP info
!!
!! COPYRIGHT
!!
!! NOTES
!!
!! SOURCE

module m_fcidump_c
 use m_definitions
 implicit none

 real(dp),allocatable,dimension(:,:)::hCORE_in
 real(dp),allocatable,dimension(:,:,:,:)::ERI_in

end module m_fcidump_c
!!***

!!****m* DoNOF/m_noft_driver_c
!! NAME
!!  m_noft_driver_c
!!
!! FUNCTION
!!  Call the m_noft_driver and serve as an interface for C
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
module m_noft_driver_c

 use iso_c_binding
 use m_fcidump_c
 use m_noft_driver

 implicit none

!!***

 public :: run_noft_c
!!***

contains

!!***
!!****f* DoNOF/run_noft_c
!! NAME
!! run_noft_c
!!
!! FUNCTION
!!  Interface to call the run_noft procedure
!!
!! INPUTS
!! Simplified arguments to be used in run_noft
!! 
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine run_noft_c(INOF,Ista,NBF_tot,NBF_occ,Nfrozen,Npairs,Ncoupled,Nbeta_elect,Nalpha_elect, &
&  imethocc,imethorb,itermax,iprintdmn,iprintswdmn,iprintints,itolLambda,ndiis,Enof,tolE,Vnn,Occ, &
&  Overlap_in,NO_COEF_in,restart,ireadGAMMAS,ireadocc,ireadCOEF,ireadFdiag,iNOTupdateocc,iNOTupdateORB) bind(C,name="run_noft_c")
 use m_definitions
 implicit none
!Arguments ------------------------------------
!scalars
 integer(c_int),intent(in)::restart
 integer(c_int),intent(in)::ireadGAMMAS,ireadocc,ireadCOEF,ireadFdiag,iNOTupdateocc,iNOTupdateORB
 integer(c_int),intent(in)::INOF,Ista,imethocc,imethorb,itermax,iprintdmn,iprintints,iprintswdmn
 integer(c_int),intent(in)::NBF_tot,NBF_occ,Nfrozen,Npairs,Ncoupled,itolLambda,ndiis  
 integer(c_int),intent(in)::Nbeta_elect,Nalpha_elect
 real(c_double),intent(inout)::Vnn,tolE
 real(c_double),intent(inout)::Enof
!arrays
 real(c_double),dimension(NBF_tot),intent(inout)::Occ
 real(c_double),dimension(NBF_tot*NBF_tot),intent(inout)::Overlap_in
 real(c_double),dimension(NBF_tot*NBF_tot),intent(inout)::NO_COEF_in
!Local variables ------------------------------
!scalars
 logical::FCIDUMP_found
 integer::iorb,jorb,korb,lorb
 integer::nbf_unit=2026
 real(c_double)::Val
 external::mo_ints_c
!arrays
 real(dp),allocatable::Overlap(:,:)
 real(dp),allocatable::NO_COEF(:,:)
 character(len=100)::ofile_name
!************************************************************************
 
 ofile_name='tmp.noft'

write(*,*) INOF,Npairs
do iorb=1,NBF_tot
 do jorb=1,NBF_tot
  write(*,*) iorb,jorb,NO_COEF_in(jorb+NBF_tot*(iorb-1))
 enddo
enddo
write(*,*) 'Hello fortran'
goto 111
 ! Allocate and initialize arrays
 allocate(Overlap(NBF_tot,NBF_tot),NO_COEF(NBF_tot,NBF_tot))
 allocate(hCORE_in(NBF_tot,NBF_tot),ERI_in(NBF_tot,NBF_tot,NBF_tot,NBF_tot))
 Overlap=zero;NO_COEF=zero;hCORE_in=zero;ERI_in=zero;
 do iorb=1,NBF_tot
  NO_COEF(iorb,iorb)=one
  do jorb=1,NBF_tot
   Overlap(jorb,iorb)=Overlap_in(jorb+NBF_tot*(iorb-1))
  enddo
 enddo

 ! Read the FCIDUMP
 inquire(file='imp_model.fcidump',exist=FCIDUMP_found)
 if(FCIDUMP_found) then
  call system('echo "-1 -1 -1 -1  0.0" >> imp_model.fcidump')
  open(newunit=nbf_unit,file='imp_model.fcidump',status='old',form='formatted')
  do
   read(nbf_unit,*) iorb,jorb,korb,lorb,Val
   if(iorb==-1 .and. jorb==-1 .and. korb==-1 .and. lorb==-1) exit
   if(iorb==0 .and. jorb==0 .and. korb==0 .and. lorb==0) Vnn=Val
   if(iorb/=0 .and. jorb/=0 .and. korb==0 .and. lorb==0) then
    hCORE_in(iorb,jorb)=Val
    hCORE_in(jorb,iorb)=Val
   endif
   if(iorb/=0 .and. jorb/=0 .and. korb/=0 .and. lorb/=0) then
    ERI_in(iorb,korb,jorb,lorb)=Val
    ERI_in(iorb,lorb,jorb,korb)=Val
    ERI_in(jorb,lorb,iorb,korb)=Val
    ERI_in(jorb,korb,iorb,lorb)=Val
    ERI_in(korb,iorb,lorb,jorb)=Val
    ERI_in(lorb,iorb,korb,jorb)=Val
    ERI_in(lorb,jorb,korb,iorb)=Val
    ERI_in(korb,jorb,lorb,iorb)=Val
   endif
  enddo
  close(nbf_unit)
 else
  write(*,*) 'imp_model.fcidump file not found!'
  stop
 endif

 ! Call the run_noft module
 if(restart/=1) then
   write(*,*) 'running NOFT module'
   call run_noft(INOF,Ista,NBF_tot,NBF_occ,Nfrozen,Npairs,Ncoupled,Nbeta_elect,Nalpha_elect, &
  &   imethocc,imethorb,itermax,iprintdmn,iprintswdmn,iprintints,itolLambda,ndiis,           &
  &   Enof,tolE,Vnn,Overlap,Occ,mo_ints_c,ofile_name,NO_COEF=NO_COEF,iNOTupdateORB=iNOTupdateORB)
 else
   write(*,*) 'running NOFT module (restart)'
   call run_noft(INOF,Ista,NBF_tot,NBF_occ,Nfrozen,Npairs,Ncoupled,Nbeta_elect,Nalpha_elect, &
  &   imethocc,imethorb,itermax,iprintdmn,iprintswdmn,iprintints,itolLambda,ndiis,           &
  &   Enof,tolE,Vnn,Overlap,Occ,mo_ints_c,ofile_name,NO_COEF=NO_COEF,restart=.true.,         &
  &   ireadGAMMAS=ireadGAMMAS,ireadOCC=ireadocc,ireadCOEF=ireadCOEF,ireadFdiag=ireadFdiag,   &
  &   iNOTupdateOCC=iNOTupdateocc,iNOTupdateORB=iNOTupdateORB)
 endif

 ! Transfer optimized coefs
 korb=1
 do iorb=1,NBF_tot
  do jorb=1,NBF_tot
   NO_COEF_in(korb)=NO_COEF(jorb,iorb)
   korb=korb+1
  enddo
 enddo

 ! Deallocate arrays
 deallocate(Overlap,NO_COEF)
 deallocate(hCORE_in,ERI_in)

111 continue
end subroutine run_noft_c
!!***

end module m_noft_driver_c
!!***

!!***
!!****f* DoNOF/mo_ints_c
!! NAME
!! mo_ints_c
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

subroutine mo_ints_c(NBF_tot, NBF_occ, NBF_jkl, Occ, DM2_JK, NO_COEF, hCORE, ERImol, ERImolJsr, ERImolLsr, &
     &              NO_COEF_cmplx, hCORE_cmplx, ERImol_cmplx, ERImolJsr_cmplx, ERImolLsr_cmplx, all_ERIs,  &
     &              Edft_xc, do_xc_dft)
 use m_definitions
 use iso_c_binding
 use m_fcidump_c
 implicit none
!Arguments ------------------------------------
!scalars
 logical, optional, intent(in) :: all_ERIs, do_xc_dft
 integer, intent(in)           :: NBF_tot, NBF_occ, NBF_jkl
!arrays
 real(dp), intent(in)              :: Occ(NBF_occ)
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
!Local variables ------------------------------
!scalars
 integer(c_int)::iorb,jorb,korb,lorb,NBF
 real(c_double)::Val
 interface
  subroutine coef_2_hcore(NO_COEF_v,NBF) bind(C,name="coef_2_hcore")
  use iso_c_binding
  real(c_double),dimension(*)::NO_COEF_v
  integer(c_int), value :: NBF
  end subroutine
  subroutine coef_2_ERI(NO_COEF_v,NBF) bind(C,name="coef_2_ERI")
  use iso_c_binding
  real(c_double),dimension(*)::NO_COEF_v
  integer(c_int), value :: NBF
  end subroutine
  subroutine hcore_ij(iorb,jorb,Val) bind(C,name="hcore_ij")
  use iso_c_binding
  integer(c_int), value :: iorb
  integer(c_int), value :: jorb
  real(c_double), value :: Val
  end subroutine
  subroutine ERI_ijkl(iorb,jorb,korb,lorb,Val) bind(C,name="ERI_ijkl")
  use iso_c_binding
  integer(c_int), value :: iorb
  integer(c_int), value :: jorb
  integer(c_int), value :: korb
  integer(c_int), value :: lorb
  real(c_double), value :: Val
  end subroutine
 end interface
!arrays
 real(c_double),allocatable::NO_COEF_vec(:)
!************************************************************************

 Val=zero
 NBF=NBF_tot
 ! Allocate arrays
 allocate(NO_COEF_vec(NBF_tot*NBF_tot))
 korb=1
 do iorb=1,NBF_tot
  do jorb=1,NBF_tot
   NO_COEF_vec(korb)=NO_COEF(jorb,iorb)
   korb=korb+1
  enddo
 enddo

 ! 'Use' the unsed variables to avoid warnings!
 if(Occ(1)==one) then
 endif
 if(present(all_ERIs)) then
 endif
 if(present(DM2_JK)) then
 endif
 if(present(do_xc_dft)) then
 endif
 if(present(Edft_xc)) then
 endif
 if(present(ERImolJsr)) then
 endif
 if(present(ERImolLsr)) then
 endif
 if(present(ERImolJsr_cmplx)) then
 endif
 if(present(ERImolLsr_cmplx)) then
 endif
 if(present(hCORE_cmplx)) then
 endif
 if(present(ERImol_cmplx)) then
 endif
 if(present(NO_COEF_cmplx)) then
 endif

 ! Transformation of hCORE
 call coef_2_hcore(NO_COEF_vec,NBF)
 do iorb=1,NBF_tot
  do jorb=1,NBF_tot
   call hcore_ij(iorb,jorb,Val)
   hCORE(iorb,jorb)=Val 
  enddo
 enddo
 !hCORE=matmul(transpose(NO_COEF),matmul(hCORE_in,NO_COEF))

 ! Transformation of ERI
 call coef_2_ERI(NO_COEF_vec,NBF)
 do iorb=1,NBF_tot
  do jorb=1,NBF_jkl
   do korb=1,NBF_jkl
    do lorb=1,NBF_jkl
     call ERI_ijkl(iorb,jorb,korb,lorb,Val)
     ERImol(iorb,jorb,korb,lorb)=Val
    enddo
   enddo
  enddo
 enddo

 ! Deallocate arrays
 deallocate(NO_COEF_vec)

end subroutine mo_ints_c
!!***
