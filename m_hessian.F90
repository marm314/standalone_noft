!!****m* DoNOF/m_hessian
!! NAME
!!  m_hessian
!!
!! FUNCTION
!!  Module to build the Hessian 
!!
!! COPYRIGHT
!! This file is distributed under the terms of the
!! GNU General Public License, see http://www.gnu.org/copyleft/gpl.txt .
!!
!!
!! PARENTS
!!  m_optorb
!!  m_noft_driver
!!
!! CHILDREN
!!  m_rdmd
!!  m_elag
!!
!! SOURCE

module m_hessian

 use m_nofoutput
 use m_definitions
 use m_rdmd
 use m_elag
 use m_integd

 implicit none


!!***
!!****t* m_hessian/hessian_t
!! NAME
!! hessian_t
!!
!! FUNCTION
!! Datatype building the (orbitals) Hessian matrix and gradients
!!
!! SOURCE

 type,public :: hessian_t

  logical::cpx_hessian=.false.  ! True for complex Hessian (i.e. complex orbitals)
  integer::NDIM_hess            ! Size of the HESSIAN
! arrays 
  real(dp),allocatable,dimension(:)::Gradient_vec             ! F_pq - F_qp (Gradient matrix)
  complex(dp),allocatable,dimension(:)::Gradient_vec_cmplx    ! F_pq - F_qp* (Gradient matrix, complex)
  real(dp),allocatable,dimension(:,:)::Hessian_mat            ! Hessian matrix 
  complex(dp),allocatable,dimension(:,:)::Hessian_mat_cmplx   ! Hessian matrix (complex)


 contains 
   procedure :: free => hessian_free
   ! Destructor.

   !procedure :: build => build_hessian
   ! Use integrals and the 1,2-RDM to build Hessian.

 end type hessian_t

 public :: hessian_init 
!!***

CONTAINS  !==============================================================================

!!***
!!****f* DoNOF/hessian_init
!! NAME
!! hessian_init
!!
!! FUNCTION
!!  Initialize the data type hessian_t 
!!
!! INPUTS
!! NBF_tot=Number of total orbitals
!!
!! OUTPUT
!!
!! PARENTS
!!  
!! CHILDREN
!!
!! SOURCE

subroutine hessian_init(HESSIANd,NBF_tot,cpx_mos)
!Arguments ------------------------------------
!scalars
 logical,intent(in)::cpx_mos
 integer,intent(in)::NBF_tot
 type(hessian_t),intent(inout)::HESSIANd
!Local variables ------------------------------
!scalars
 real(dp)::totMEM
!arrays
 character(len=200)::msg
!************************************************************************

 HESSIANd%cpx_hessian=cpx_mos
 HESSIANd%NDIM_hess=NBF_tot*NBF_tot
 ! Calculate memory needed
 if(cpx_mos) then
  totMEM=2*(NBF_tot*NBF_tot*NBF_tot*NBF_tot+NBF_tot*NBF_tot)
 else
  totMEM=NBF_tot*NBF_tot*NBF_tot*NBF_tot+NBF_tot*NBF_tot
 endif
 totMEM=8*totMEM       ! Bytes
 totMEM=totMEM*tol6    ! Bytes to Mb  
 if(totMEM>thousand) then  ! Mb to Gb
  write(msg,'(a,f10.3,a)') 'Mem. required for storing HESSIANd object ',totMEM*tol3,' Gb'
 elseif(totMEM<one) then   ! Mb to Kb
  write(msg,'(a,f10.3,a)') 'Mem. required for storing HESSIANd object ',totMEM*thousand,' Kb'
 else                      ! Mb
  write(msg,'(a,f10.3,a)') 'Mem. required for storing HESSIANd object ',totMEM,' Mb'
 endif
 call write_output(msg)
 ! Allocate arrays
 if(cpx_mos) then
  allocate(HESSIANd%Gradient_vec_cmplx(HESSIANd%NDIM_hess))
  allocate(HESSIANd%Hessian_mat_cmplx(HESSIANd%NDIM_hess,HESSIANd%NDIM_hess)) 
 else 
  allocate(HESSIANd%Gradient_vec(HESSIANd%NDIM_hess))
  allocate(HESSIANd%Hessian_mat(HESSIANd%NDIM_hess,HESSIANd%NDIM_hess)) 
 endif
 
end subroutine hessian_init
!!***

!!***
!!****f* DoNOF/hessian_free
!! NAME
!! hessian_free
!!
!! FUNCTION
!!  Free allocated arrays of the data type hessian_t 
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

subroutine hessian_free(HESSIANd)
!Arguments ------------------------------------
!scalars
 class(hessian_t),intent(inout)::HESSIANd
!Local variables ------------------------------
!scalars
!arrays
!************************************************************************

 if(HESSIANd%cpx_hessian) then
  deallocate(HESSIANd%Gradient_vec_cmplx) 
  deallocate(HESSIANd%Hessian_mat_cmplx) 
 else
  deallocate(HESSIANd%Gradient_vec) 
  deallocate(HESSIANd%Hessian_mat) 
 endif

end subroutine hessian_free
!!***

end module m_hessian
!!***
