!!****m* DoNOF/m_tz_pCCD_amplitudes
!! NAME
!!  m_tz_pCCD_amplitudes
!!
!! FUNCTION
!!  Module prepared to compute the pCCD t- and z-amplitudes for fixed orbitals and the 2-RDM
!!
!! PARENTS
!!  m_noft_driver
!!
!! CHILDREN
!!  m_elag
!!  m_diagf
!!  m_e_grad_occ
!!
!! SOURCE
module m_tz_pCCD_amplitudes

 use m_nofoutput
 use m_gammatodm2 
 use m_rdmd
 use m_integd
 use m_elag
 use m_lbfgs_intern
 use m_diagf
 use m_e_grad_occ
 use m_e_grad_occ_cpx

 implicit none

 real(dp),parameter,private::pccd_tol=1.0d-10 ! Amplitudes converged when the residue 2-norm is below pccd_tol
 integer,parameter,private::pccd_max_iter=200,pccd_krylov_max_iter=200,pccd_diis_max_iter=1000,pccd_diis_size=8
 logical,parameter,private::pccd_lm_fallback=.false. ! First t solve: LM from the starting amplitudes if Newton fails

 private :: calc_t_residues,calc_E_sd,Grad_r_wrt_t,amplitudes_from_occs,mat_to_vec,vec_to_mat,is_finite
 private :: pccd_res,pccd_jac,solve_amplitudes,solve_newton,solve_lm,solve_diis,solve_krylov,solve_lbfgs
!!***

 public :: calc_tz_pCCD_amplitudes
!!***

contains

!!***
!!****f* DoNOF/calc_tz_pCCD_amplitudes
!! NAME
!!  calc_tz_pCCD_amplitudes
!!
!! FUNCTION
!!  Solve the pCCD t- and z-amplitudes for fixed orbitals and compute the energy with the resulting 1-RDM and 2-RDM.
!!  The amplitude solvers follow the ReRDMFT implementation (Occ_opt/pCCD.cpp):
!!   - t: R_ia(t) = 0 (quadratic in t) solved to |R|_2 < 1e-10.
!!   - z: R_z(z) = G_ov + J(t)^T z = 0 (linear in z), where J=dR/dt is the Jacobian at the converged t.
!!   - imethod = 0 (default) or 2 DIIS (Newton if it fails), 1 Newton-Raphson, 3 Jacobian-free Newton-Krylov (GMRES +
!!               Armijo), 4 Levenberg-Marquardt (Newton for z), otherwise (e.g. 5) L-BFGS on 0.5 |R|^2.
!!  If RDMd%pccd_fixed_occs (READ_OCCUPANCIES) the amplitudes are built from the input occs. and are not solved.
!!
!! INPUTS
!!  Vnn=Nuclear-nuclear rep. energy
!!
!! OUTPUT
!!  Energy=Sum of nuclear and electronic energy
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine calc_tz_pCCD_amplitudes(ELAGd,RDMd,INTEGd,Vnn,Energy,Phases,iter_global,imethod,keep_occs,only_phases) 
!Arguments ------------------------------------
!scalars
 logical,optional,intent(in)::only_phases
 logical,intent(in)::keep_occs
 integer,intent(in)::imethod
 integer,intent(inout)::iter_global
 real(dp),intent(in)::Vnn
 real(dp),intent(inout)::Energy
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
!arrays
 real(dp),dimension(RDMd%NBF_occ,RDMd%NBF_occ)::Phases
!Local variables ------------------------------
!scalars
 logical::converged_t,converged_z,first_solve,only_phases_,solve_amp
 integer::iter_t,iter_z,iorb,iorb1,iorb2,iorb3,nOcc,nVir
 real(dp)::rnorm_t,rnorm_z
 real(dp)::Ecorr_new,Ediff,Esingle_det,Energy_dm
!arrays
 character(len=200)::msg
 real(dp),allocatable,dimension(:)::x_vec,G_vec
 real(dp),allocatable,dimension(:,:)::Jacobian_t
!************************************************************************

 Ecorr_new=zero; Ediff=Energy; rnorm_t=zero; rnorm_z=zero;
 iter_t=0; iter_z=0; converged_t=.true.; converged_z=.true.;
 nOcc=RDMd%Npairs; nVir=RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs)
 only_phases_=.false.; first_solve=(iter_global==-1)
 if(present(only_phases)) only_phases_=only_phases
 solve_amp=(.not.keep_occs .and. .not.only_phases_ .and. .not.RDMd%pccd_fixed_occs)

 ! Build diag elements of the Lambda matrix (with HF 1-RDM and 2-RDM) and compute SD energy
 call ELAGd%build_sd_diag(RDMd,INTEGd)
 call calc_E_sd(RDMd,INTEGd,Esingle_det)

 ! READ_OCCUPANCIES: t and z are built from the occs. and kept fixed along the run
 if(RDMd%pccd_fixed_occs .and. first_solve) call amplitudes_from_occs(RDMd)

 ! The starting amplitudes are the checkpoint ones (if read) or zero. Newton, L-BFGS, Krylov and LM start from them, while
 ! DIIS builds the MP2-like start x=-R(0)/D if they are all zero.
 if(solve_amp) then

  allocate(x_vec(nOcc*nVir),G_vec(nOcc*nVir),Jacobian_t(nOcc*nVir,nOcc*nVir))

  ! G_ia = (ia|ai)
  do iorb=1,nOcc
   iorb1=iorb+RDMd%Nfrozen
   do iorb2=1,nVir
    iorb3=iorb2+RDMd%Nfrozen+nOcc
    if(INTEGd%complex_ints) then
     G_vec(iorb2+(iorb-1)*nVir)=real(INTEGd%ERImol_cmplx(iorb1,iorb3,iorb3,iorb1))
    else
     G_vec(iorb2+(iorb-1)*nVir)=INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)
    endif
   enddo
  enddo

  ! t-amplitudes
  call mat_to_vec(RDMd%t_pccd_old,x_vec,nOcc,nVir)
  call solve_amplitudes(1,ELAGd,RDMd,INTEGd,imethod,nOcc*nVir,x_vec,Jacobian_t,G_vec,first_solve,converged_t,iter_t,rnorm_t)
  call vec_to_mat(x_vec,RDMd%t_pccd_old,nOcc,nVir)
  RDMd%t_pccd=RDMd%t_pccd_old

  ! E_corr = sum_ia G_ia t_ia (exactly linear in t)
  Ecorr_new=dot_product(G_vec,x_vec)

  ! z-amplitudes (with fixed t-amplitudes): J(t)^T z = -G, the Jacobian is evaluated once at the converged t
  call Grad_r_wrt_t(ELAGd,RDMd,INTEGd,Jacobian_t,nOcc,nVir)
  call mat_to_vec(RDMd%z_pccd_old,x_vec,nOcc,nVir)
  call solve_amplitudes(2,ELAGd,RDMd,INTEGd,imethod,nOcc*nVir,x_vec,Jacobian_t,G_vec,first_solve,converged_z,iter_z,rnorm_z)
  call vec_to_mat(x_vec,RDMd%z_pccd_old,nOcc,nVir)
  RDMd%z_pccd=RDMd%z_pccd_old

  deallocate(x_vec,G_vec,Jacobian_t)

 endif

 ! Calc. the final Energy using new RDMs
 if(only_phases_) then
  RDMd%INOF=-1 ! Call pCCD RDM constructors
  if(INTEGd%complex_ints) then
   call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_dm,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
   & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx,only_phases=only_phases_)
  else
   call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy_dm,Phases,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
   & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr,only_phases=only_phases_)
  endif
  RDMd%INOF=70 ! Recover call to PNOF7-PHASES
 else
  iter_global=iter_global+1
  if(INTEGd%complex_ints) then
   call calc_E_occ_cmplx(RDMd,RDMd%GAMMAs_old,Energy_dm,Phases,INTEGd%hCORE_cmplx,INTEGd%ERI_J_cmplx,INTEGd%ERI_K_cmplx, &
   & INTEGd%ERI_L_cmplx,INTEGd%ERI_Jsr_cmplx,INTEGd%ERI_Lsr_cmplx)
  else
   call calc_E_occ(RDMd,RDMd%GAMMAs_old,Energy_dm,Phases,INTEGd%hCORE,INTEGd%ERI_J,INTEGd%ERI_K, &
   & INTEGd%ERI_L,INTEGd%ERI_Jsr,INTEGd%ERI_Lsr)
  endif
  if(solve_amp) then
   write(msg,'(a,f15.6)') 'Single-Det. Energy |0>        =',Esingle_det
   call write_output(msg)
   write(msg,'(a,f15.6)') 'Correlation Energy w.r.t. |0> =',Ecorr_new
   call write_output(msg)
   write(msg,'(a,f15.6)') 'Energy from 1-RDM and 2-RDM   =',Energy_dm+Vnn
   call write_output(msg)
   Energy=Esingle_det+Ecorr_new
  else
   write(msg,'(a,f15.6)') 'Energy from 1-RDM and 2-RDM   =',Energy_dm+Vnn
   call write_output(msg)
   Energy=Energy_dm
  endif
  write(msg,'(a,f15.6,a,i6,a,i6,a)') 'T-,Z-amp. opt. energy= ',Energy+Vnn,' after ',iter_t,' t-iter. and',&
   & iter_z,' z-iter.'
  call write_output(msg)
  if(iter_global>0) then
   Ediff=Energy-Ediff
   write(msg,'(a,es15.6)') 'Error t-residues        =      ',rnorm_t
   call write_output(msg)
   write(msg,'(a,es15.6)') 'Error z-residues        =      ',rnorm_z
   call write_output(msg)
   write(msg,'(a,f19.10)') 'Energy difference amp. opt.=',Ediff
   call write_output(msg)
   write(msg,'(a)') 'Current occ. numbers '
   call write_output(msg)
   do iorb=1,(RDMd%NBF_occ/10)*10,10
    write(msg,'(f12.6,9f11.6)') two*RDMd%occ(iorb:iorb+9)
    call write_output(msg)
   enddo
   iorb=(RDMd%NBF_occ/10)*10+1
   write(msg,'(a,*(f11.6))') ' ',two*RDMd%occ(iorb:)
   call write_output(msg)
  endif
  write(msg,'(a,i6)') 'Number of global iter. ',iter_global
  call write_output(msg)
  write(msg,'(a)') ' '
  call write_output(msg)
 endif

end subroutine calc_tz_pCCD_amplitudes
!!***

!!***
!!****f* DoNOF/amplitudes_from_occs
!! NAME
!!  amplitudes_from_occs
!!
!! FUNCTION
!!  Build the pCCD amplitudes that reproduce the given window occupations (READ_OCCUPANCIES). Pair k couples the occupied
!!  HOMO-k with the virtual LUMO+k: t = -sqrt(n_a/n_i), z = -sqrt(n_i n_a), so that n_i = 1 - t z and n_a = t z.
!!  Every pair must have n_i + n_a = 1, the unpaired occupied n = 1 and the unpaired virtuals n = 0.
!!
!! INPUTS
!!  RDMd%pccd_window_occ=Occs. in [0:1] of orbitals Nfrozen+1:NBF_occ (occupied lowest..HOMO, then virtuals)
!!
!! OUTPUT
!!  RDMd%t_pccd,RDMd%t_pccd_old,RDMd%z_pccd,RDMd%z_pccd_old
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine amplitudes_from_occs(RDMd)
!Arguments ------------------------------------
!scalars
 type(rdm_t),intent(inout)::RDMd
!Local variables ------------------------------
!scalars
 logical::has_i,has_a
 integer::k,nOcc,nVir
 real(dp)::n_i,n_a
 real(dp),parameter::kTol=1.0d-8
!arrays
 character(len=200)::msg
!************************************************************************

 nOcc=RDMd%Npairs; nVir=RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs)
 RDMd%t_pccd=zero; RDMd%z_pccd=zero;
 do k=0,max(nOcc,nVir)-1
  has_i=(k<nOcc); has_a=(k<nVir)
  n_i=zero; n_a=one
  if(has_i) n_i=RDMd%pccd_window_occ(nOcc-k)
  if(has_a) n_a=RDMd%pccd_window_occ(nOcc+k+1)
  if(has_i .and. has_a) then
   if(abs(n_i+n_a-one)>kTol .or. n_i<-kTol .or. n_a<-kTol .or. n_i<kTol) then
    write(msg,'(a,i5,a,f12.6,a,f12.6,a)') 'READ_OCCUPANCIES (pCCD): pair ',k,' (occupied HOMO-k with virtual LUMO+k) has occs. ', &
    & two*n_i,' and ',two*n_a,'; pCCD needs n_occupied + n_virtual = 1 (and n_occupied > 0)'
    call write_output(msg)
    stop
   endif
   RDMd%t_pccd(nOcc-k,k+1)=-dsqrt(max(zero,n_a)/n_i)
   RDMd%z_pccd(nOcc-k,k+1)=-dsqrt(max(zero,n_a)*n_i)
  elseif(has_i .and. abs(n_i-one)>kTol) then
   write(msg,'(a,i5,a,f12.6)') 'READ_OCCUPANCIES (pCCD): occupied HOMO-',k,' has no virtual partner and needs occ. 2, got ',two*n_i
   call write_output(msg)
   stop
  elseif(has_a .and. .not.has_i .and. abs(n_a)>kTol) then
   write(msg,'(a,i5,a,f12.6)') 'READ_OCCUPANCIES (pCCD): virtual LUMO+',k,' has no occupied partner and needs occ. 0, got ',two*n_a
   call write_output(msg)
   stop
  endif
 enddo
 RDMd%t_pccd_old=RDMd%t_pccd; RDMd%z_pccd_old=RDMd%z_pccd;

end subroutine amplitudes_from_occs
!!***

!!***
!!****f* DoNOF/mat_to_vec
!! NAME
!!  mat_to_vec / vec_to_mat
!!
!! FUNCTION
!!  Flatten (nOcc,nVir) amplitudes row-major, i.e. index a+(i-1)*nVir like the Jacobian Grad_r_wrt_t, and back.
!!
!! SOURCE

subroutine mat_to_vec(mat,vec,nOcc,nVir)
!Arguments ------------------------------------
 integer,intent(in)::nOcc,nVir
 real(dp),dimension(nOcc,nVir),intent(in)::mat
 real(dp),dimension(nOcc*nVir),intent(out)::vec
!Local variables ------------------------------
 integer::iorb,iorb1
!************************************************************************
 do iorb=1,nOcc
  do iorb1=1,nVir
   vec(iorb1+(iorb-1)*nVir)=mat(iorb,iorb1)
  enddo
 enddo
end subroutine mat_to_vec

subroutine vec_to_mat(vec,mat,nOcc,nVir)
!Arguments ------------------------------------
 integer,intent(in)::nOcc,nVir
 real(dp),dimension(nOcc*nVir),intent(in)::vec
 real(dp),dimension(nOcc,nVir),intent(out)::mat
!Local variables ------------------------------
 integer::iorb,iorb1
!************************************************************************
 do iorb=1,nOcc
  do iorb1=1,nVir
   mat(iorb,iorb1)=vec(iorb1+(iorb-1)*nVir)
  enddo
 enddo
end subroutine vec_to_mat
!!***

!!***
!!****f* DoNOF/pccd_res
!! NAME
!!  pccd_res
!!
!! FUNCTION
!!  Residual vector r(x). mode=1: t residue R_ia(t). mode=2: z residue G + J(t)^T z (Jt=dR/dt at the converged t).
!!
!! SOURCE

subroutine pccd_res(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,r)
!Arguments ------------------------------------
 integer,intent(in)::mode,nAmp
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
 real(dp),dimension(nAmp),intent(in)::x,Gvec
 real(dp),dimension(nAmp,nAmp),intent(in)::Jt
 real(dp),dimension(nAmp),intent(out)::r
!Local variables ------------------------------
 real(dp),allocatable,dimension(:,:)::y_ij
!************************************************************************
 if(mode==1) then
  allocate(y_ij(RDMd%Npairs,RDMd%Npairs))
  call vec_to_mat(x,RDMd%t_pccd_old,RDMd%Npairs,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs))
  call calc_t_residues(ELAGd,RDMd,INTEGd,y_ij)
  call mat_to_vec(RDMd%tz_residue,r,RDMd%Npairs,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs))
  deallocate(y_ij)
 else
  r=Gvec
  call DGEMV('T',nAmp,nAmp,one,Jt,nAmp,x,1,one,r,1)
 endif
end subroutine pccd_res
!!***

!!***
!!****f* DoNOF/pccd_jac
!! NAME
!!  pccd_jac
!!
!! FUNCTION
!!  Jacobian dr/dx. mode=1: dR_ia/dt_kd. mode=2: dR_z/dz = J(t)^T (independent of z).
!!
!! SOURCE

subroutine pccd_jac(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Jac)
!Arguments ------------------------------------
 integer,intent(in)::mode,nAmp
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
 real(dp),dimension(nAmp),intent(in)::x
 real(dp),dimension(nAmp,nAmp),intent(in)::Jt
 real(dp),dimension(nAmp,nAmp),intent(out)::Jac
!************************************************************************
 if(mode==1) then
  call vec_to_mat(x,RDMd%t_pccd_old,RDMd%Npairs,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs))
  call Grad_r_wrt_t(ELAGd,RDMd,INTEGd,Jac,RDMd%Npairs,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs))
 else
  Jac=transpose(Jt)
 endif
end subroutine pccd_jac
!!***

!!***
!!****f* DoNOF/solve_amplitudes
!! NAME
!!  solve_amplitudes
!!
!! FUNCTION
!!  Drive r(x)=0 to |r|_2 < pccd_tol for the t (mode=1) or z (mode=2) amplitudes, with the strategy selected by imethod
!!  (see calc_tz_pCCD_amplitudes). DIIS falls back to Newton from the starting amplitudes. In the first t solve, a
!!  non-converged Newton can fall back to Levenberg-Marquardt (pccd_lm_fallback). A non-converged Krylov stops the run.
!!
!! SOURCE

subroutine solve_amplitudes(mode,ELAGd,RDMd,INTEGd,imethod,nAmp,x,Jt,Gvec,first_solve,converged,iter,rnorm)
!Arguments ------------------------------------
 logical,intent(in)::first_solve
 logical,intent(out)::converged
 integer,intent(in)::mode,imethod,nAmp
 integer,intent(out)::iter
 real(dp),intent(out)::rnorm
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
 real(dp),dimension(nAmp),intent(inout)::x
 real(dp),dimension(nAmp),intent(in)::Gvec
 real(dp),dimension(nAmp,nAmp),intent(in)::Jt
!Local variables ------------------------------
 logical::diis_failed,converged_lm
 integer::method,iter_lm
 real(dp)::rnorm_lm
 character(len=200)::msg
 real(dp),allocatable,dimension(:)::x0,x_lm
!************************************************************************

 method=imethod
 if(method==0) method=2 ! Default (imethocc=0): DIIS
 if(mode==2 .and. method==4) method=1 ! R_z is linear: LM damping buys nothing, one exact Newton step solves it
 diis_failed=.false.; converged=.false.; iter=0; rnorm=zero;
 allocate(x0(nAmp),x_lm(nAmp))
 x0=x

 if(method==4) then
  write(msg,'(a)') 'Using Levenberg-Marquardt to optimize amplitudes'
  call write_output(msg)
  call solve_lm(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,converged,iter,rnorm)
  if(.not.converged) then
   write(msg,'(a,i5,a,es12.4,a)') 'WARNING: pCCD amplitude solve: LM did not converge (',iter,' iterations, |R| = ',rnorm, &
   & '); the unconverged amplitudes are used'
   call write_output(msg)
  endif
 elseif(method==3) then
  write(msg,'(a)') 'Using Newton-Krylov to optimize amplitudes'
  call write_output(msg)
  call solve_krylov(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,converged,iter,rnorm)
  if(.not.converged) then
   write(msg,'(a,i5,a,es12.4,a)') 'Error: pCCD amplitude solve failed: KRYLOV did not converge (',iter,' iterations, |R| = ', &
   & rnorm,')'
   call write_output(msg)
   stop
  endif
 else
  if(method==2) then
   write(msg,'(a)') 'Using DIIS to optimize amplitudes'
   call write_output(msg)
   call solve_diis(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,converged,iter,rnorm)
   if(.not.converged) then
    diis_failed=.true.
    if(is_finite(rnorm)) then
     write(msg,'(a,i5,a,es12.4,a)') 'pCCD amplitude solve: DIIS did not converge (',iter,' iterations, |R| = ',rnorm, &
     & '); NEWTON reactivated'
    else
     write(msg,'(a,i5,a)') 'pCCD amplitude solve: DIIS diverged (|R| not finite after ',iter,' iterations); NEWTON reactivated'
    endif
    call write_output(msg)
    x=x0
   endif
  endif
  if(method==1 .or. diis_failed) then
   if(.not.diis_failed) then
    write(msg,'(a)') 'Using Newton-Raphson to optimize amplitudes'
    call write_output(msg)
   endif
   call solve_newton(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,converged,iter,rnorm)
   if(.not.converged .and. pccd_lm_fallback .and. first_solve .and. mode==1) then
    write(msg,'(a,i5,a,es12.4,a)') 'pCCD amplitude solve: NEWTON did not converge (',iter,' iterations, |R| = ',rnorm, &
    & '); LM reactivated'
    call write_output(msg)
    x_lm=x0
    call solve_lm(mode,ELAGd,RDMd,INTEGd,nAmp,x_lm,Jt,Gvec,converged_lm,iter_lm,rnorm_lm)
    if(converged_lm .or. rnorm_lm<rnorm) then
     x=x_lm; converged=converged_lm; iter=iter_lm; rnorm=rnorm_lm;
    endif
   endif
   if(.not.converged .and. (diis_failed .or. (pccd_lm_fallback .and. first_solve .and. mode==1))) then
    write(msg,'(a,i5,a,es12.4,a)') 'WARNING: pCCD amplitude solve: NEWTON did not converge (',iter,' iterations, |R| = ',rnorm, &
    & '); the unconverged amplitudes are used'
    call write_output(msg)
   endif
  elseif(method/=2) then
   write(msg,'(a)') 'Using L-BFGS to optimize amplitudes'
   call write_output(msg)
   call solve_lbfgs(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,converged,iter,rnorm)
  endif
 endif

 deallocate(x0,x_lm)

end subroutine solve_amplitudes
!!***

!!***
!!****f* DoNOF/is_finite
!! NAME
!!  is_finite
!!
!! SOURCE

function is_finite(val) result(finite)
 logical::finite
 real(dp),intent(in)::val
 finite=(val==val .and. abs(val)<huge(val))
end function is_finite
!!***

!!***
!!****f* DoNOF/solve_newton
!! NAME
!!  solve_newton
!!
!! FUNCTION
!!  Newton-Raphson: J dx = r with the analytic Jacobian rebuilt at every iteration, x <- x - dx. It converges in a single
!!  iteration for the linear z equation.
!!
!! SOURCE

subroutine solve_newton(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,converged,iter,rnorm)
!Arguments ------------------------------------
 logical,intent(out)::converged
 integer,intent(in)::mode,nAmp
 integer,intent(out)::iter
 real(dp),intent(out)::rnorm
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
 real(dp),dimension(nAmp),intent(inout)::x
 real(dp),dimension(nAmp),intent(in)::Gvec
 real(dp),dimension(nAmp,nAmp),intent(in)::Jt
!Local variables ------------------------------
 integer::info
 integer,allocatable,dimension(:)::IPIV
 real(dp),allocatable,dimension(:)::r,step
 real(dp),allocatable,dimension(:,:)::Jac
 character(len=200)::msg
!************************************************************************

 allocate(r(nAmp),step(nAmp),IPIV(nAmp),Jac(nAmp,nAmp))
 iter=0
 call pccd_res(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,r)
 rnorm=dsqrt(sum(r*r))
 do while(iter<pccd_max_iter .and. rnorm>=pccd_tol)
  call pccd_jac(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Jac)
  step=r
  call DGETRF(nAmp,nAmp,Jac,nAmp,IPIV,info)
  if(info==0) call DGETRS('N',nAmp,1,Jac,nAmp,IPIV,step,nAmp,info)
  if(info/=0) then
   write(msg,'(a,i5)') 'Error solving the Newton equations (singular Jacobian) at iter ',iter
   call write_output(msg)
   exit
  endif
  x=x-step
  call pccd_res(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,r)
  rnorm=dsqrt(sum(r*r))
  iter=iter+1
 enddo
 converged=(rnorm<pccd_tol)
 deallocate(r,step,IPIV,Jac)

end subroutine solve_newton
!!***

!!***
!!****f* DoNOF/solve_lm
!! NAME
!!  solve_lm
!!
!! FUNCTION
!!  Levenberg-Marquardt on 0.5|r|^2 with the analytic Jacobian: (J^T J + lambda diag(J^T J)) dx = -J^T r, lambda adapted
!!  by the gain ratio (Nielsen); the steps that do not lower |r| are rejected. Meant for a (nearly) singular Jacobian.
!!
!! SOURCE

subroutine solve_lm(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,converged,iter,rnorm)
!Arguments ------------------------------------
 logical,intent(out)::converged
 integer,intent(in)::mode,nAmp
 integer,intent(out)::iter
 real(dp),intent(out)::rnorm
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
 real(dp),dimension(nAmp),intent(inout)::x
 real(dp),dimension(nAmp),intent(in)::Gvec
 real(dp),dimension(nAmp,nAmp),intent(in)::Jt
!Local variables ------------------------------
 logical::accepted
 integer::k,attempt,info
 integer,allocatable,dimension(:)::IPIV
 real(dp)::lambda,nu,new_norm,predicted,actual,rho
 real(dp),allocatable,dimension(:)::f,fnew,xn,dx,g
 real(dp),allocatable,dimension(:,:)::Jac,jtj,amat
!************************************************************************

 allocate(f(nAmp),fnew(nAmp),xn(nAmp),dx(nAmp),g(nAmp),IPIV(nAmp),Jac(nAmp,nAmp),jtj(nAmp,nAmp),amat(nAmp,nAmp))
 call pccd_res(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,f)
 rnorm=dsqrt(sum(f*f))
 lambda=1.0d-3; nu=two; iter=0;
 do while(iter<pccd_max_iter .and. is_finite(rnorm) .and. rnorm>=pccd_tol)
  call pccd_jac(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Jac)
  call DGEMM('T','N',nAmp,nAmp,nAmp,one,Jac,nAmp,Jac,nAmp,zero,jtj,nAmp) ! J^T J
  call DGEMV('T',nAmp,nAmp,one,Jac,nAmp,f,1,zero,g,1)                    ! J^T r
  accepted=.false.
  do attempt=1,40
   if(accepted) exit
   amat=jtj
   do k=1,nAmp
    amat(k,k)=amat(k,k)+lambda*max(jtj(k,k),1.0d-12)
   enddo
   dx=g
   call DGESV(nAmp,1,amat,nAmp,IPIV,dx,nAmp,info)
   if(info/=0) then ! Singular damped system: increase the damping
    lambda=lambda*nu; nu=nu*two;
    cycle
   endif
   dx=-dx
   xn=x+dx
   call pccd_res(mode,ELAGd,RDMd,INTEGd,nAmp,xn,Jt,Gvec,fnew)
   new_norm=dsqrt(sum(fnew*fnew))
   ! predicted reduction of 0.5|r|^2 for the damped step: 0.5 dx^T (lambda D dx - J^T r)
   predicted=zero
   do k=1,nAmp
    predicted=predicted+dx(k)*(lambda*max(jtj(k,k),1.0d-12)*dx(k)-g(k))
   enddo
   predicted=half*predicted
   actual=half*(rnorm*rnorm-new_norm*new_norm)
   if(is_finite(new_norm) .and. predicted>zero .and. actual>zero) then
    rho=actual/predicted
    x=xn; f=fnew; rnorm=new_norm;
    lambda=lambda*max(one/3.0d0,one-(two*rho-one)**3)
    nu=two
    accepted=.true.
   else
    lambda=lambda*nu; nu=nu*two;
   endif
  enddo
  if(.not.accepted) exit ! No descent step found
  iter=iter+1
 enddo
 converged=(is_finite(rnorm) .and. rnorm<pccd_tol)
 deallocate(f,fnew,xn,dx,g,IPIV,Jac,jtj,amat)

end subroutine solve_lm
!!***

!!***
!!****f* DoNOF/solve_diis
!! NAME
!!  solve_diis
!!
!! FUNCTION
!!  Jacobi update x <- x - r(x)/D with FIXED denominators D (diagonal of the Jacobian at x=0), the MP2-like start
!!  x=-r(0)/D when the starting amplitudes are all zero, and DIIS (pccd_diis_size vectors, the update x_new-x as error
!!  vector) on the iterates.
!!
!! SOURCE

subroutine solve_diis(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,converged,iter,rnorm)
!Arguments ------------------------------------
 logical,intent(out)::converged
 integer,intent(in)::mode,nAmp
 integer,intent(out)::iter
 real(dp),intent(out)::rnorm
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
 real(dp),dimension(nAmp),intent(inout)::x
 real(dp),dimension(nAmp),intent(in)::Gvec
 real(dp),dimension(nAmp,nAmp),intent(in)::Jt
!Local variables ------------------------------
 integer::k,i,j,nhist,info
 integer,allocatable,dimension(:)::IPIV
 real(dp),allocatable,dimension(:)::denom,r,xnew,err,xzero,coef
 real(dp),allocatable,dimension(:,:)::Jac,Xh,Eh,Bmat
!************************************************************************

 allocate(denom(nAmp),r(nAmp),xnew(nAmp),err(nAmp),xzero(nAmp),Jac(nAmp,nAmp))
 xzero=zero
 call pccd_jac(mode,ELAGd,RDMd,INTEGd,nAmp,xzero,Jt,Jac)
 do k=1,nAmp
  denom(k)=Jac(k,k)
  if(abs(denom(k))<1.0d-12) denom(k)=sign(1.0d-12,denom(k))
 enddo
 if(all(x==zero)) then ! MP2-like start x = -R(0) / D
  call pccd_res(mode,ELAGd,RDMd,INTEGd,nAmp,xzero,Jt,Gvec,r)
  x=-r/denom
 endif
 if(pccd_diis_size>=2) then
  allocate(Xh(nAmp,pccd_diis_size),Eh(nAmp,pccd_diis_size))
 endif
 nhist=0; iter=0;
 call pccd_res(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,r)
 rnorm=dsqrt(sum(r*r))
 do while(iter<pccd_diis_max_iter .and. is_finite(rnorm) .and. rnorm>=pccd_tol)
  xnew=x-r/denom
  err=xnew-x
  if(pccd_diis_size>=2) then
   ! Store (xnew,err), discarding the oldest pair if the history is full
   if(nhist==pccd_diis_size) then
    do i=1,nhist-1
     Xh(:,i)=Xh(:,i+1); Eh(:,i)=Eh(:,i+1);
    enddo
    nhist=nhist-1
   endif
   nhist=nhist+1
   Xh(:,nhist)=xnew; Eh(:,nhist)=err;
   x=xnew
   if(nhist>=2) then
    ! Pulay: min |sum_i c_i e_i| subject to sum_i c_i = 1
    allocate(Bmat(nhist+1,nhist+1),coef(nhist+1),IPIV(nhist+1))
    Bmat=zero
    do i=1,nhist
     do j=1,nhist
      Bmat(i,j)=dot_product(Eh(:,i),Eh(:,j))
     enddo
     Bmat(i,nhist+1)=-one; Bmat(nhist+1,i)=-one;
    enddo
    coef=zero; coef(nhist+1)=-one;
    call DGESV(nhist+1,1,Bmat,nhist+1,IPIV,coef,nhist+1,info)
    if(info==0) then
     x=zero
     do i=1,nhist
      x=x+coef(i)*Xh(:,i)
     enddo
    endif
    deallocate(Bmat,coef,IPIV)
   endif
  else
   x=xnew
  endif
  call pccd_res(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,r)
  rnorm=dsqrt(sum(r*r))
  iter=iter+1
 enddo
 converged=(is_finite(rnorm) .and. rnorm<pccd_tol)
 deallocate(denom,r,xnew,err,xzero,Jac)
 if(allocated(Xh)) deallocate(Xh,Eh)

end subroutine solve_diis
!!***

!!***
!!****f* DoNOF/solve_krylov
!! NAME
!!  solve_krylov
!!
!! FUNCTION
!!  Jacobian-free Newton-Krylov: each Newton step solves J dx = -r by restarted GMRES (Arnoldi with modified Gram-Schmidt
!!  and Givens rotations) using finite-difference J v products and the forcing term eta=min(0.1,|r|), followed by an
!!  Armijo backtracking line search on 0.5|r|^2.
!!
!! SOURCE

subroutine solve_krylov(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,converged,iter,rnorm)
!Arguments ------------------------------------
 logical,intent(out)::converged
 integer,intent(in)::mode,nAmp
 integer,intent(out)::iter
 real(dp),intent(out)::rnorm
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
 real(dp),dimension(nAmp),intent(inout)::x
 real(dp),dimension(nAmp),intent(in)::Gvec
 real(dp),dimension(nAmp,nAmp),intent(in)::Jt
!Local variables ------------------------------
 integer::m,j,i,ii,jj,restart,ls,used
 real(dp)::rdiff,xnorm,eta,bnorm,beta,dotp,wn,t1,den,alpha,fn_norm
 real(dp),allocatable,dimension(:)::f,dx,rr,w,xn,fn,cs,sn,g,y
 real(dp),allocatable,dimension(:,:)::Q,H
!************************************************************************

 allocate(f(nAmp),dx(nAmp),rr(nAmp),w(nAmp),xn(nAmp),fn(nAmp))
 m=min(nAmp,100)
 allocate(Q(nAmp,m+1),H(m+1,m),cs(m),sn(m),g(m+1),y(m))
 rdiff=dsqrt(epsilon(one))
 iter=0
 call pccd_res(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,f)
 rnorm=dsqrt(sum(f*f))
 do while(iter<pccd_krylov_max_iter .and. is_finite(rnorm) .and. rnorm>=pccd_tol)
  xnorm=dsqrt(sum(x*x))
  eta=min(0.1d0,rnorm)
  dx=zero; rr=-f; bnorm=dsqrt(sum(rr*rr));
  do restart=1,10
   beta=dsqrt(sum(rr*rr))
   if(beta<=eta*bnorm) exit
   Q(:,1)=rr/beta
   H=zero; g=zero; g(1)=beta; used=0;
   do j=1,m
    call matvec(Q(:,j),w)
    do i=1,j
     dotp=dot_product(w,Q(:,i))
     H(i,j)=dotp
     w=w-dotp*Q(:,i)
    enddo
    wn=dsqrt(sum(w*w))
    H(j+1,j)=wn
    do i=1,j-1 ! Apply previous rotations
     t1=cs(i)*H(i,j)+sn(i)*H(i+1,j)
     H(i+1,j)=-sn(i)*H(i,j)+cs(i)*H(i+1,j)
     H(i,j)=t1
    enddo
    den=dsqrt(H(j,j)**2+H(j+1,j)**2)
    if(den>zero) then
     cs(j)=H(j,j)/den; sn(j)=H(j+1,j)/den;
    else
     cs(j)=one; sn(j)=zero;
    endif
    H(j,j)=den
    g(j+1)=-sn(j)*g(j)
    g(j)=cs(j)*g(j)
    used=j
    if(abs(g(j+1))<=eta*bnorm .or. H(j,j)==zero .or. wn<1.0d-14) exit
    Q(:,j+1)=w/wn
   enddo
   y=zero ! Back substitution on the rotated Hessenberg
   do ii=used,1,-1
    dotp=g(ii)
    do jj=ii+1,used
     dotp=dotp-H(ii,jj)*y(jj)
    enddo
    if(H(ii,ii)/=zero) y(ii)=dotp/H(ii,ii)
   enddo
   do jj=1,used
    dx=dx+y(jj)*Q(:,jj)
   enddo
   call matvec(dx,w)
   rr=-f-w
  enddo
  ! Armijo backtracking on 0.5|r|^2
  alpha=one
  do ls=1,30
   xn=x+alpha*dx
   call pccd_res(mode,ELAGd,RDMd,INTEGd,nAmp,xn,Jt,Gvec,fn)
   fn_norm=dsqrt(sum(fn*fn))
   if(is_finite(fn_norm) .and. fn_norm*fn_norm<=(one-2.0d-4*alpha)*rnorm*rnorm) exit
   alpha=half*alpha
  enddo
  x=xn; f=fn; rnorm=fn_norm;
  iter=iter+1
 enddo
 converged=(is_finite(rnorm) .and. rnorm<pccd_tol)
 deallocate(f,dx,rr,w,xn,fn,Q,H,cs,sn,g,y)

contains

 subroutine matvec(v,jv)
 ! Finite-difference J v = (r(x+h v)-r(x))/h
 real(dp),dimension(nAmp),intent(in)::v
 real(dp),dimension(nAmp),intent(out)::jv
 real(dp)::vnorm,h
 real(dp),dimension(nAmp)::xp,fp
  vnorm=dsqrt(sum(v*v))
  if(vnorm==zero) then
   jv=zero
  else
   h=rdiff*max(one,xnorm)/vnorm
   xp=x+h*v
   call pccd_res(mode,ELAGd,RDMd,INTEGd,nAmp,xp,Jt,Gvec,fp)
   jv=(fp-f)/h
  endif
 end subroutine matvec

end subroutine solve_krylov
!!***

!!***
!!****f* DoNOF/solve_lbfgs
!! NAME
!!  solve_lbfgs
!!
!! FUNCTION
!!  L-BFGS on f(x)=0.5|r(x)|^2 with the exact chain-rule gradient J(x)^T r(x), never solving the Jacobian system.
!!
!! SOURCE

subroutine solve_lbfgs(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,converged,iter,rnorm)
!Arguments ------------------------------------
 logical,intent(out)::converged
 integer,intent(in)::mode,nAmp
 integer,intent(out)::iter
 real(dp),intent(out)::rnorm
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
 real(dp),dimension(nAmp),intent(inout)::x
 real(dp),dimension(nAmp),intent(in)::Gvec
 real(dp),dimension(nAmp,nAmp),intent(in)::Jt
!Local variables ------------------------------
 logical::diagco
 integer,parameter::msave=7
 integer::Mtosave,Nwork,iflag
 integer,dimension(2)::info_print
 real(dp)::fval
 real(dp),allocatable,dimension(:)::Work,diag,r,grad
 real(dp),allocatable,dimension(:,:)::Jac
!************************************************************************

 Nwork=nAmp*(2*msave+1)+2*msave
 Mtosave=5; info_print(1)= -1; info_print(2)= 0; diagco= .false.;
 iter=0; iflag=0;
 allocate(Work(Nwork),diag(nAmp),r(nAmp),grad(nAmp),Jac(nAmp,nAmp))
 do
  call pccd_res(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,r)
  fval=half*sum(r*r)
  call pccd_jac(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Jac)
  call DGEMV('T',nAmp,nAmp,one,Jac,nAmp,r,1,zero,grad,1)
  call LBFGS_INTERN(nAmp,Mtosave,x,fval,grad,diagco,diag,info_print,pccd_tol,tol16,Work,iflag)
  if(iflag<=0) exit
  iter=iter+1
  !  We allow at most 2000 evaluations of function and gradient
  if(iter==2000) exit
 enddo
 call pccd_res(mode,ELAGd,RDMd,INTEGd,nAmp,x,Jt,Gvec,r)
 rnorm=dsqrt(sum(r*r))
 converged=(rnorm<pccd_tol)
 deallocate(Work,diag,r,grad,Jac)

end subroutine solve_lbfgs
!!***

!!***
!!****f* DoNOF/Grad_r_wrt_t
!! NAME
!!  Grad_r_wrt_t
!!
!! FUNCTION
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

subroutine Grad_r_wrt_t(ELAGd,RDMd,INTEGd,Grad_res_wrt_tia,nOcc,nVir)
!Arguments ------------------------------------
!scalars
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
 integer,intent(in)::nOcc,nVir
!arrays
 real(dp),dimension(nOcc*nVir,nOcc*nVir),intent(out)::Grad_res_wrt_tia
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorb3,iorb4,iorb5,iorbi,iorba,iorbk,iorbd
 integer::index_ia,index_jb
 real(dp)::sum_tmp
!arrays
!************************************************************************

 Grad_res_wrt_tia=zero

 ! d r_ia / d t_kd
 ! Residue r_ia = r_iorb,iorb1
 do iorb=1,nOcc
  iorbi=iorb+RDMd%Nfrozen
  do iorb1=1,nVir
   iorba=iorb1+nOcc+RDMd%Nfrozen
   index_ia=iorb1+(iorb-1)*nVir
   ! d r_iorb,iorb1 / d t_iorb2,iorb3 with t_kd = t_iorb2,iorb3
   do iorb2=1,nOcc
    iorbk=iorb2+RDMd%Nfrozen
    do iorb3=1,nVir
     iorbd=iorb3+nOcc+RDMd%Nfrozen
     index_jb=iorb3+(iorb2-1)*nVir
     ! i=k and a=d
     if(iorb==iorb2 .and. iorb1==iorb3) then
      if(INTEGd%complex_ints) then
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb)   &
       &   -two*(two*real(INTEGd%ERImol_cmplx(iorbi,iorba,iorbi,iorba))          &
       &   -real(INTEGd%ERImol_cmplx(iorba,iorbi,iorbi,iorba))                   &
       &   -real(INTEGd%ERImol_cmplx(iorba,iorbi,iorbi,iorba))*RDMd%t_pccd_old(iorb,iorb1))
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb)   &
       &   +two*real(INTEGd%ERImol_cmplx(iorba,iorbi,iorbi,iorba))*RDMd%t_pccd_old(iorb,iorb1)
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb)   &
       &   +two*(ELAGd%Lambdas_pp(iorba)-ELAGd%Lambdas_pp(iorbi))
       sum_tmp=zero
       do iorb4=1,nOcc
        iorb5=iorb4+RDMd%Nfrozen
        sum_tmp=sum_tmp+real(INTEGd%ERImol_cmplx(iorba,iorb5,iorb5,iorba))*RDMd%t_pccd_old(iorb4,iorb1)
       enddo
       do iorb4=1,nVir
        iorb5=iorb4+nOcc+RDMd%Nfrozen
        sum_tmp=sum_tmp+real(INTEGd%ERImol_cmplx(iorb5,iorbi,iorbi,iorb5))*RDMd%t_pccd_old(iorb,iorb4)
       enddo
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb)-two*sum_tmp
      else
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb)   &
       &   -two*(two*INTEGd%ERImol(iorbi,iorba,iorbi,iorba)                      &
       &   -INTEGd%ERImol(iorba,iorbi,iorbi,iorba)                               &
       &   -INTEGd%ERImol(iorba,iorbi,iorbi,iorba)*RDMd%t_pccd_old(iorb,iorb1))
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb)   &
       &   +two*INTEGd%ERImol(iorba,iorbi,iorbi,iorba)*RDMd%t_pccd_old(iorb,iorb1)
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb)   &
       &   +two*(ELAGd%Lambdas_pp(iorba)-ELAGd%Lambdas_pp(iorbi))
       sum_tmp=zero
       do iorb4=1,nOcc
        iorb5=iorb4+RDMd%Nfrozen
        sum_tmp=sum_tmp+INTEGd%ERImol(iorba,iorb5,iorb5,iorba)*RDMd%t_pccd_old(iorb4,iorb1)
       enddo
       do iorb4=1,nVir
        iorb5=iorb4+nOcc+RDMd%Nfrozen
        sum_tmp=sum_tmp+INTEGd%ERImol(iorb5,iorbi,iorbi,iorb5)*RDMd%t_pccd_old(iorb,iorb4)
       enddo
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb)-two*sum_tmp
      endif
     endif
     ! a=d
     if(iorb1==iorb3) then
      if(INTEGd%complex_ints) then
       do iorb4=1,nVir
        iorb5=iorb4+nOcc+RDMd%Nfrozen
        Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb) &
        &   +real(INTEGd%ERImol_cmplx(iorb5,iorbk,iorbk,iorb5))*RDMd%t_pccd_old(iorb,iorb4)
       enddo
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb) &
       &   +real(INTEGd%ERImol_cmplx(iorbi,iorbk,iorbk,iorbi))
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb) &
       &   -two*real(INTEGd%ERImol_cmplx(iorba,iorbk,iorbk,iorba))*RDMd%t_pccd_old(iorb,iorb1)
      else
       do iorb4=1,nVir
        iorb5=iorb4+nOcc+RDMd%Nfrozen
        Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb) &
        &   +INTEGd%ERImol(iorb5,iorbk,iorbk,iorb5)*RDMd%t_pccd_old(iorb,iorb4)
       enddo
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb) &
       &   +INTEGd%ERImol(iorbi,iorbk,iorbk,iorbi)
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb) &
       &   -two*INTEGd%ERImol(iorba,iorbk,iorbk,iorba)*RDMd%t_pccd_old(iorb,iorb1)
      endif
     endif
     ! i=k
     if(iorb==iorb2) then
      if(INTEGd%complex_ints) then
       do iorb4=1,nOcc
        iorb5=iorb4+RDMd%Nfrozen
        Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb) &
        &   +real(INTEGd%ERImol_cmplx(iorbd,iorb5,iorb5,iorbd))*RDMd%t_pccd_old(iorb4,iorb1)
       enddo
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb) &
       &   +real(INTEGd%ERImol_cmplx(iorbd,iorba,iorba,iorbd))
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb) &
       &   -two*real(INTEGd%ERImol_cmplx(iorbd,iorbi,iorbi,iorbd)*RDMd%t_pccd_old(iorb,iorb1))
      else
       do iorb4=1,nOcc
        iorb5=iorb4+RDMd%Nfrozen
        Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb) &
        &   +INTEGd%ERImol(iorbd,iorb5,iorb5,iorbd)*RDMd%t_pccd_old(iorb4,iorb1)
       enddo
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb) &
       &   +INTEGd%ERImol(iorbd,iorba,iorba,iorbd)
       Grad_res_wrt_tia(index_ia,index_jb)=Grad_res_wrt_tia(index_ia,index_jb) &
       &   -two*INTEGd%ERImol(iorbd,iorbi,iorbi,iorbd)*RDMd%t_pccd_old(iorb,iorb1)
      endif
     endif
    enddo
   enddo
  enddo
 enddo

end subroutine Grad_r_wrt_t
!!***

!!***
!!****f* DoNOF/calc_t_residues
!! NAME
!!  calc_t_residues
!!
!! FUNCTION
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

subroutine calc_t_residues(ELAGd,RDMd,INTEGd,y_ij) 
!Arguments ------------------------------------
!scalars
 type(elag_t),intent(inout)::ELAGd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
!arrays
 real(dp),dimension(RDMd%Npairs,RDMd%Npairs)::y_ij
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1,iorb2,iorb3,iorb4,iorb5
 real(dp)::sum_tmp
!arrays
!************************************************************************

 ! Build intermediate y_ij = sum_b v_bb^jj t_i^b 
 y_ij=zero 
 do iorb=1,RDMd%Npairs ! Occ
  iorb1=iorb+RDMd%Nfrozen
  do iorb2=1,RDMd%Npairs ! Occ
   iorb3=iorb2+RDMd%Nfrozen
   do iorb4=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
    iorb5=iorb4+RDMd%Nfrozen+RDMd%Npairs
    if(INTEGd%complex_ints) then
     y_ij(iorb,iorb2)=y_ij(iorb,iorb2)+RDMd%t_pccd_old(iorb,iorb4)*real(INTEGd%ERImol_cmplx(iorb3,iorb5,iorb5,iorb3))
    else
     y_ij(iorb,iorb2)=y_ij(iorb,iorb2)+RDMd%t_pccd_old(iorb,iorb4)*INTEGd%ERImol(iorb3,iorb5,iorb5,iorb3)
    endif 
   enddo
  enddo
 enddo

 ! Build the t residue
 RDMd%tz_residue=zero
 do iorb=1,RDMd%Npairs ! Occ
  iorb1=iorb+RDMd%Nfrozen
  do iorb2=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
   iorb3=iorb2+RDMd%Nfrozen+RDMd%Npairs
   if(INTEGd%complex_ints) then
    RDMd%tz_residue(iorb,iorb2)=real(INTEGd%ERImol_cmplx(iorb1,iorb3,iorb3,iorb1))
   else
    RDMd%tz_residue(iorb,iorb2)=INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)
   endif
   sum_tmp=ELAGd%Lambdas_pp(iorb3)-ELAGd%Lambdas_pp(iorb1)
   do iorb4=1,RDMd%Npairs ! Occ
    iorb5=iorb4+RDMd%Nfrozen
    if(INTEGd%complex_ints) then
     RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+RDMd%t_pccd_old(iorb4,iorb2)* &
     &        real(INTEGd%ERImol_cmplx(iorb1,iorb5,iorb5,iorb1))
     sum_tmp=sum_tmp-RDMd%t_pccd_old(iorb4,iorb2)*real(INTEGd%ERImol_cmplx(iorb3,iorb5,iorb5,iorb3))
    else       
     RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+RDMd%t_pccd_old(iorb4,iorb2)* &
     &        INTEGd%ERImol(iorb1,iorb5,iorb5,iorb1)
     sum_tmp=sum_tmp-RDMd%t_pccd_old(iorb4,iorb2)*INTEGd%ERImol(iorb3,iorb5,iorb5,iorb3)
    endif 
    RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+y_ij(iorb,iorb4)*RDMd%t_pccd_old(iorb4,iorb2)
   enddo
   do iorb4=1,RDMd%NBF_occ-(RDMd%Nfrozen+RDMd%Npairs) ! Virt
    iorb5=iorb4+RDMd%Nfrozen+RDMd%Npairs
    if(INTEGd%complex_ints) then
     RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+RDMd%t_pccd_old(iorb,iorb4)* &
     &        real(INTEGd%ERImol_cmplx(iorb3,iorb5,iorb5,iorb3))
     sum_tmp=sum_tmp-RDMd%t_pccd_old(iorb,iorb4)*real(INTEGd%ERImol_cmplx(iorb1,iorb5,iorb5,iorb1))
    else
     RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+RDMd%t_pccd_old(iorb,iorb4)* &
     &        INTEGd%ERImol(iorb3,iorb5,iorb5,iorb3)
     sum_tmp=sum_tmp-RDMd%t_pccd_old(iorb,iorb4)*INTEGd%ERImol(iorb1,iorb5,iorb5,iorb1)
    endif
   enddo
   if(INTEGd%complex_ints) then
    sum_tmp=sum_tmp-real(two*INTEGd%ERImol_cmplx(iorb1,iorb3,iorb1,iorb3)-INTEGd%ERImol_cmplx(iorb1,iorb3,iorb3,iorb1) &
   &         -INTEGd%ERImol_cmplx(iorb1,iorb3,iorb3,iorb1)*RDMd%t_pccd_old(iorb,iorb2))
   else
    sum_tmp=sum_tmp-(two*INTEGd%ERImol(iorb1,iorb3,iorb1,iorb3)-INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1) &
   &         -INTEGd%ERImol(iorb1,iorb3,iorb3,iorb1)*RDMd%t_pccd_old(iorb,iorb2))
   endif
   sum_tmp=two*sum_tmp*RDMd%t_pccd_old(iorb,iorb2)
   RDMd%tz_residue(iorb,iorb2)=RDMd%tz_residue(iorb,iorb2)+sum_tmp
  enddo
 enddo

end subroutine calc_t_residues
!!***

!!***
!!****f* DoNOF/calc_E_sd
!! NAME
!!  calc_E_sd
!!
!! FUNCTION
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

subroutine calc_E_sd(RDMd,INTEGd,Energy_sd)
!Arguments ------------------------------------
!scalars
 real(dp),intent(inout)::Energy_sd
 type(rdm_t),intent(inout)::RDMd
 type(integ_t),intent(inout)::INTEGd
!arrays
!Local variables ------------------------------
!scalars
 integer::iorb,iorb1
!arrays
!************************************************************************

 Energy_sd=zero
 if(INTEGd%complex_ints) then
  do iorb=1,RDMd%Nfrozen+RDMd%Npairs
   Energy_sd=Energy_sd+two*real(INTEGd%hCORE_cmplx(iorb,iorb))
   do iorb1=1,RDMd%Nfrozen+RDMd%Npairs
    Energy_sd=Energy_sd+real(two*INTEGd%ERImol_cmplx(iorb,iorb1,iorb,iorb1)-INTEGd%ERImol_cmplx(iorb,iorb1,iorb1,iorb))
   enddo
  enddo
 else
  do iorb=1,RDMd%Nfrozen+RDMd%Npairs
   Energy_sd=Energy_sd+two*INTEGd%hCORE(iorb,iorb)
   do iorb1=1,RDMd%Nfrozen+RDMd%Npairs
    Energy_sd=Energy_sd+two*INTEGd%ERImol(iorb,iorb1,iorb,iorb1)-INTEGd%ERImol(iorb,iorb1,iorb1,iorb)
   enddo
  enddo
 endif

end subroutine calc_E_sd
!!***

end module m_tz_pCCD_amplitudes
!!***
