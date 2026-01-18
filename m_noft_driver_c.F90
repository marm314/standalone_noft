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
 integer::iorb,jorb,korb
 external::mo_ints_c
!arrays
 real(dp),allocatable::Overlap(:,:)
 real(dp),allocatable::NO_COEF(:,:)
 character(len=100)::ofile_name
!************************************************************************
 
 ofile_name='tmp.noft'

 ! Allocate and initialize arrays
 allocate(Overlap(NBF_tot,NBF_tot),NO_COEF(NBF_tot,NBF_tot))
 Overlap=zero;NO_COEF=zero;
 do iorb=1,NBF_tot
  do jorb=1,NBF_tot
   Overlap(jorb,iorb)=Overlap_in(jorb+NBF_tot*(iorb-1))
   NO_COEF(jorb,iorb)=NO_COEF_in(jorb+NBF_tot*(iorb-1))
  enddo
 enddo

 ! Call the run_noft module
 if(restart==0) then
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
 integer(c_int)::iiorb,jjorb,kkorb,llorb
 real(c_double)::Val
 interface
  subroutine coef_2_hcore(NO_COEF_v,NBF) bind(C,name="coef_2_hcore")
  use iso_c_binding
  real(c_double),dimension(*)::NO_COEF_v
  integer(c_int), intent(inout) :: NBF
  end subroutine
  subroutine coef_2_ERI(NO_COEF_v,NBF) bind(C,name="coef_2_ERI")
  use iso_c_binding
  real(c_double),dimension(*)::NO_COEF_v
  integer(c_int), intent(inout) :: NBF
  end subroutine
  subroutine hcore_ij(iorb,jorb,NBF,Val) bind(C,name="hcore_ij")
  use iso_c_binding
  integer(c_int), intent(inout) :: iorb
  integer(c_int), intent(inout) :: jorb
  integer(c_int), intent(inout) :: NBF
  real(c_double), intent(inout) :: Val
  end subroutine
  subroutine ERI_ijkl(iorb,jorb,korb,lorb,NBF,Val) bind(C,name="ERI_ijkl")
  use iso_c_binding
  integer(c_int), intent(inout) :: iorb
  integer(c_int), intent(inout) :: jorb
  integer(c_int), intent(inout) :: korb
  integer(c_int), intent(inout) :: lorb
  integer(c_int), intent(inout) :: NBF
  real(c_double), intent(inout) :: Val
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
  iiorb=iorb
  do jorb=1,NBF_tot
   jjorb=jorb
   call hcore_ij(iiorb,jjorb,NBF,Val)
   hCORE(iorb,jorb)=Val
  enddo
 enddo

 ! Transformation of ERI
 call coef_2_ERI(NO_COEF_vec,NBF)
 do iorb=1,NBF_tot
  iiorb=iorb
  do jorb=1,NBF_jkl
   jjorb=jorb
   do korb=1,NBF_jkl
    kkorb=korb
    do lorb=1,NBF_jkl
     llorb=lorb
     call ERI_ijkl(iiorb,jjorb,kkorb,llorb,NBF,Val)
     ERImol(iorb,jorb,korb,lorb)=Val
    enddo
   enddo
  enddo
 enddo

 ! Deallocate arrays
 deallocate(NO_COEF_vec)

end subroutine mo_ints_c
!!***
