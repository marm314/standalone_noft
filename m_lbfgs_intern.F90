!======================================================================!
!  m_lbfgs_intern.f90                                                  !
!                                                                       !
!  Limited-memory BFGS optimizer (Jorge Nocedal, July 1990) with the   !
!  More'-Thuente line search (MCSRCH / MCSTEP), modernized from fixed- !
!  form Fortran 77 to free-form Fortran 90+.                           !
!                                                                       !
!  WHAT CHANGED (dialect only -- the algorithm is untouched):          !
!    - Free source form, IMPLICIT NONE everywhere, explicit typing     !
!      of every variable (no more "IMPLICIT DOUBLE PRECISION(A-H,O-Z)")!
!    - The COMMON /LB3/ MP,LP,GTOL,STPMIN,STPMAX block is replaced by  !
!      public module variables of the same name. Any external driver  !
!      that used to poke that COMMON block must now do:                !
!          use m_lbfgs_intern, only: gtol, stpmin, stpmax, lp, mp      !
!          gtol = 0.1d0   ! etc.                                       !
!    - Old-style numbered DO/CONTINUE loops -> DO / END DO, and a      !
!      couple of manual dot-product loops -> DOT_PRODUCT / array       !
!      section syntax.                                                 !
!    - The computed GOTO used for reverse-communication dispatch is    !
!      replaced by SELECT CASE.                                        !
!    - Plain error exits (bad N/M, non-positive DIAG element, failed   !
!      line search) are now structured IF blocks with a direct RETURN  !
!      instead of GOTO to a shared error label.                        !
!    - D-prefixed intrinsics (DSQRT, DABS, DMAX1, ...) -> generic ones !
!      (SQRT, ABS, MAX, ...).                                          !
!                                                                       !
!  WHAT DELIBERATELY DID NOT CHANGE:                                   !
!    LBFGS_INTERN and MCSRCH both use *reverse communication*: they    !
!    return to the caller asking it to compute F/G (or DIAG) and are   !
!    then called again to resume exactly where they left off. That     !
!    pattern needs a handful of re-entry points into the middle of a   !
!    subroutine, which in Fortran still has to be expressed with a     !
!    small number of GOTOs/labels (this is exactly what every other    !
!    modern reverse-communication solver, e.g. L-BFGS-B, still does).  !
!    Those GOTOs are kept, minimal, and heavily commented below; all   !
!    the *other* GOTOs in the original code (plain error exits, loops) !
!    have been removed as described above.                             !
!                                                                       !
!  NOTE ON PRECISION: the working precision "wp" is defined locally    !
!  below as double precision (kind(1.0d0)), matching the original.     !
!  If the host project's m_definitions module already exports a        !
!  working-precision kind parameter, prefer that one instead and       !
!  delete the local "wp" definition here.                              !
!======================================================================!
module m_lbfgs_intern

  use m_definitions

  implicit none

  private
  public :: lbfgs_intern
  public :: mp, lp, gtol, stpmin, stpmax

  integer, parameter :: wp = kind(1.0d0)

  ! Replaces the old COMMON /LB3/ MP,LP,GTOL,STPMIN,STPMAX.
  !
  !   MP      unit number for monitoring output controlled by IPRINT
  !           (default 6... but see the note on LP below).
  !   LP      unit number for error messages; <=0 suppresses printing.
  !   GTOL    controls the accuracy of the line search MCSRCH; must be
  !           > 1.0D-04 (0.1 is a typical "loose" value).
  !   STPMIN, STPMAX
  !           lower/upper bounds for the line-search step length.
  integer,  save :: mp     = 11
  integer,  save :: lp     = 11
  real(wp), save :: gtol   = 9.0d-01
  real(wp), save :: stpmin = 1.0d-20
  real(wp), save :: stpmax = 1.0d+20

contains

!======================================================================!
!                                                                      !
!                           LBFGS SUBROUTINES                          !
!                                                                      !
!======================================================================!

  subroutine lbfgs_intern(n, m, x, f, g, diagco, diag, iprint, eps, xtol, w, iflag)
    !--------------------------------------------------------------------
    !     LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
    !                       JORGE NOCEDAL, *** July 1990 ***
    !
    !     This subroutine solves the unconstrained minimization problem
    !
    !                      min F(x),    x = (x1,x2,...,xN),
    !
    !     using the limited memory BFGS method. The routine is especially
    !     effective on problems involving a large number of variables. In
    !     a typical iteration of this method an approximation Hk to the
    !     inverse of the Hessian is obtained by applying M BFGS updates to
    !     a diagonal matrix Hk0, using information from the previous M
    !     steps. The user specifies the number M, which determines the
    !     amount of storage required by the routine. The user may also
    !     provide the diagonal matrices Hk0 if not satisfied with the
    !     default choice. The algorithm is described in "On the limited
    !     memory BFGS method for large scale optimization", by D. Liu and
    !     J. Nocedal, Mathematical Programming B 45 (1989) 503-528.
    !
    !     The user is required to calculate the function value F and its
    !     gradient G. In order to allow the user complete control over
    !     these computations, reverse communication is used. The routine
    !     must be called repeatedly under the control of the parameter
    !     IFLAG.
    !
    !     The steplength is determined at each iteration by means of the
    !     line search routine MCSRCH, which is a slight modification of
    !     the routine CSRCH written by More' and Thuente.
    !
    !     ARGUMENTS
    !     ---------
    !     N       number of variables (N > 0), unchanged by the routine.
    !     M       number of BFGS corrections to keep (M > 0); 3 <= M <= 7
    !             is typically recommended.
    !     X       length-N array. On initial entry, the starting point.
    !             On exit with IFLAG=0, the best point found.
    !     F       on initial entry and on re-entry with IFLAG=1, must
    !             hold F(X).
    !     G       on initial entry and on re-entry with IFLAG=1, must
    !             hold the gradient of F at X.
    !     DIAGCO  .TRUE. if the caller wants to supply the diagonal Hk0
    !             at every iteration (returned via IFLAG=2); otherwise
    !             .FALSE. and LBFGS uses a default.
    !     DIAG    length-N array; if DIAGCO, must hold Hk0 on initial
    !             entry and on re-entry with IFLAG=2. All elements must
    !             be positive.
    !     IPRINT  IPRINT(1) < 0: no output; = 0: first/last iteration
    !             only; > 0: every IPRINT(1) iterations.
    !             IPRINT(2) selects how much is printed (0..3), see LB1.
    !     EPS     terminate when ||G|| < EPS * max(1,||X||).
    !     XTOL    estimate of machine precision; the line search
    !             terminates if the interval of uncertainty is narrower
    !             (relatively) than XTOL.
    !     W       workspace of length N*(2M+1)+2M; must not be altered
    !             by the caller.
    !     IFLAG   set to 0 on the initial call. On return:
    !               0  converged, no errors
    !               1  the caller must evaluate F and G at X and call
    !                  again with everything else unchanged
    !               2  (DIAGCO only) the caller must evaluate DIAG and
    !                  call again
    !              -1  the line search MCSRCH failed (see the write-up
    !                  printed on unit LP for the INFO code)
    !              -2  a diagonal element of DIAG was not positive
    !              -3  N or M was not positive
    !
    !     Other routines called directly: DAXPY, DDOT, LB1, MCSRCH
    !--------------------------------------------------------------------
    integer,  intent(in)    :: n, m
    integer,  intent(in)    :: iprint(2)
    integer,  intent(inout) :: iflag
    real(wp), intent(inout) :: x(n), diag(n)
    real(wp), intent(inout) :: w(n*(2*m+1)+2*m)
    real(wp), intent(in)    :: f, g(n), eps, xtol
    logical,  intent(in)    :: diagco

    real(wp), parameter :: one = 1.0_wp, zero = 0.0_wp

    integer,  save :: iter, nfun, point, ispt, iypt, maxfev, info
    integer,  save :: bound, npt, cp, i, nfev, inmc, iycn, iscn
    real(wp), save :: gnorm, stp1, ftol, stp, ys, yy, sq, yr, beta, xnorm
    logical,  save :: finish

    real(wp), external :: ddot

    ! NOTE: the original code forced LP = 6 on *every* call, silently
    ! overriding whatever the driver had set through the old COMMON
    ! block (whose DATA statement actually defaulted LP to 11). That
    ! quirk is preserved here for strict behavioural compatibility.
    lp = 6

    if (iflag == 0) then
       goto 10
    else
       select case (iflag)
       case (1)
          goto 172
       case (2)
          goto 100
       end select
    end if

    !------------------------------------------------------------------
    ! INITIALIZE
    !------------------------------------------------------------------
10  continue
    iter = 0
    if (n <= 0 .or. m <= 0) then
       iflag = -3
       if (lp > 0) write(lp, 240)
       return
    end if
    if (gtol <= 1.0d-04) then
       if (lp > 0) write(lp, 245)
       gtol = 9.0d-01
    end if
    nfun   = 1
    point  = 0
    finish = .false.
    if (diagco) then
       do i = 1, n
          if (diag(i) <= zero) then
             iflag = -2
             if (lp > 0) write(lp, 235) i
             return
          end if
       end do
    else
       diag(1:n) = 1.0_wp
    end if

    ! THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
    ! ---------------------------------------
    ! THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND OTHER
    !     TEMPORARY INFORMATION.
    ! LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
    ! LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED IN THE
    !     FORMULA THAT COMPUTES H*G.
    ! LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH STEPS.
    ! LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M GRADIENT
    !     DIFFERENCES.
    ! THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
    ! CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
    ispt = n + 2*m
    iypt = ispt + n*m
    w(ispt+1:ispt+n) = -g(1:n) * diag(1:n)
    gnorm = sqrt(ddot(n, g, 1, g, 1))
    stp1  = one / gnorm

    ! PARAMETERS FOR LINE SEARCH ROUTINE
    ftol   = 1.0d-4
    maxfev = 20

    if (iprint(1) >= 0) call lb1(iprint, iter, nfun, gnorm, n, m, x, f, g, stp, finish)

    !------------------------------------------------------------------
    !     MAIN ITERATION LOOP
    !------------------------------------------------------------------
80  continue
    iter  = iter + 1
    info  = 0
    bound = iter - 1
    if (iter == 1) goto 165
    if (iter > m) bound = m

    ys = ddot(n, w(iypt+npt+1), 1, w(ispt+npt+1), 1)
    if (.not. diagco) then
       yy = ddot(n, w(iypt+npt+1), 1, w(iypt+npt+1), 1)
       diag(1:n) = ys / yy
    else
       iflag = 2
       return
    end if

100 continue
    if (diagco) then
       do i = 1, n
          if (diag(i) <= zero) then
             iflag = -2
             if (lp > 0) write(lp, 235) i
             return
          end if
       end do
    end if

    ! COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
    ! "Updating quasi-Newton matrices with limited storage",
    ! Mathematics of Computation, Vol.24, No.151, pp. 773-782.
    ! ---------------------------------------------------------
    cp = point
    if (point == 0) cp = m
    w(n+cp) = one / ys
    w(1:n)  = -g(1:n)
    cp = point
    do i = 1, bound
       cp = cp - 1
       if (cp == -1) cp = m - 1
       sq   = ddot(n, w(ispt+cp*n+1), 1, w, 1)
       inmc = n + m + cp + 1
       iycn = iypt + cp*n
       w(inmc) = w(n+cp+1) * sq
       call daxpy(n, -w(inmc), w(iycn+1), 1, w, 1)
    end do

    w(1:n) = diag(1:n) * w(1:n)

    do i = 1, bound
       yr   = ddot(n, w(iypt+cp*n+1), 1, w, 1)
       beta = w(n+cp+1) * yr
       inmc = n + m + cp + 1
       beta = w(inmc) - beta
       iscn = ispt + cp*n
       call daxpy(n, beta, w(iscn+1), 1, w, 1)
       cp = cp + 1
       if (cp == m) cp = 0
    end do

    ! STORE THE NEW SEARCH DIRECTION
    ! ------------------------------
    w(ispt+point*n+1:ispt+point*n+n) = w(1:n)

    ! OBTAIN THE ONE-DIMENSIONAL MINIMIZER OF THE FUNCTION BY USING
    ! THE LINE SEARCH ROUTINE MCSRCH
    ! ----------------------------------------------------
165 continue
    nfev = 0
    stp  = one
    if (iter == 1) stp = stp1
    w(1:n) = g(1:n)

172 continue
    call mcsrch(n, x, f, g, w(ispt+point*n+1), stp, ftol, xtol, maxfev, info, nfev, diag)
    if (info == -1) then
       iflag = 1
       return
    end if
    if (info /= 1) then
       iflag = -1
       if (lp > 0) write(lp, 200) info
       return
    end if
    nfun = nfun + nfev

    ! COMPUTE THE NEW STEP AND GRADIENT CHANGE
    ! -----------------------------------------
    npt = point*n
    do i = 1, n
       w(ispt+npt+i) = stp * w(ispt+npt+i)
       w(iypt+npt+i) = g(i) - w(i)
    end do
    point = point + 1
    if (point == m) point = 0

    ! TERMINATION TEST
    ! ----------------
    gnorm = sqrt(ddot(n, g, 1, g, 1))
    xnorm = sqrt(ddot(n, x, 1, x, 1))
    xnorm = max(1.0_wp, xnorm)
    if (gnorm/xnorm <= eps) finish = .true.

    if (iprint(1) >= 0) call lb1(iprint, iter, nfun, gnorm, n, m, x, f, g, stp, finish)
    if (finish) then
       iflag = 0
       return
    end if
    goto 80

    ! FORMATS
    ! -------
200 format(/' IFLAG= -1 ',/' LINE SEARCH FAILED. SEE', &
            ' DOCUMENTATION OF ROUTINE MCSRCH',/' ERROR RETURN', &
            ' OF LINE SEARCH: INFO= ',i2,/, &
            ' POSSIBLE CAUSES: FUNCTION OR GRADIENT ARE INCORRECT',/, &
            ' OR INCORRECT TOLERANCES')
235 format(/' IFLAG= -2',/' THE',i5,'-TH DIAGONAL ELEMENT OF THE',/, &
            ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
240 format(/' IFLAG= -3',/' IMPROPER INPUT PARAMETERS (N OR M', &
            ' ARE NOT POSITIVE)')
245 format(/'  GTOL IS LESS THAN OR EQUAL TO 1.D-04', &
            / ' IT HAS BEEN RESET TO 9.D-01')

  end subroutine lbfgs_intern

  !--------------------------------------------------------------------
  subroutine lb1(iprint, iter, nfun, gnorm, n, m, x, f, g, stp, finish)
    !--------------------------------------------------------------------
    ! Prints monitoring information; the frequency and amount of output
    ! are controlled by IPRINT (see LBFGS_INTERN's header).
    !--------------------------------------------------------------------
    integer,  intent(in) :: iprint(2), iter, nfun, n, m
    real(wp), intent(in) :: x(n), g(n), f, gnorm, stp
    logical,  intent(in) :: finish

    integer :: i

    if (iter == 0) then
       write(mp, 10)
       write(mp, *)
       write(mp, 20) n, m
       write(mp, 30) f, gnorm
       if (iprint(2) >= 1) then
          write(mp, 40)
          write(mp, 50) (x(i), i = 1, n)
          write(mp, 60)
          write(mp, 50) (g(i), i = 1, n)
       end if
    else
       if ((iprint(1) == 0) .and. (iter /= 1 .and. .not. finish)) return
       if (iprint(1) /= 0) then
          if (mod(iter-1, iprint(1)) == 0 .or. finish) then
             if (iprint(2) > 1 .and. iter > 1) write(mp, 70)
             write(mp, 80) iter, nfun, f, gnorm, stp
          else
             return
          end if
       else
          if (iprint(2) > 1 .and. finish) write(mp, 70)
          write(mp, 80) iter, nfun, f, gnorm, stp
       end if
       if (iprint(2) == 2 .or. iprint(2) == 3) then
          if (finish) then
             write(mp, 90) f
          else
             write(mp, 40)
          end if
          write(mp, 50) (x(i), i = 1, n)
          if (iprint(2) == 3) then
             write(mp, 60)
             write(mp, 50) (g(i), i = 1, n)
          end if
       end if
       if (finish) then
          write(mp, 100)
          write(mp, *)
          write(mp, 90) f
       end if
    end if

10  format('  Start of LBFGS optimization details')
20  format('  N=', i5, '   NUMBER OF CORRECTIONS=', i2, &
            /, '       INITIAL VALUES')
30  format(' F= ', 1pd10.3, '   GNORM= ', 1pd10.3)
40  format(' VECTOR X= ')
50  format(6(2x, 1pd10.3))
60  format(' GRADIENT VECTOR G= ')
70  format(/'   I   NFN', 4x, 'FUNC', 8x, 'GNORM', 7x, 'STEPLENGTH'/)
80  format(2(i4, 1x), 3x, 3(1pd10.3, 2x))
90  format(' Final objective value = ', 1pd10.3)
100 format(/' THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS.', &
            /' IFLAG = 0')

  end subroutine lb1

  !--------------------------------------------------------------------
  subroutine mcsrch(n, x, f, g, s, stp, ftol, xtol, maxfev, info, nfev, wa)
    !--------------------------------------------------------------------
    !                     LINE SEARCH SUBROUTINE MCSRCH
    !
    !     A slight modification of the subroutine CSRCH of More' and
    !     Thuente. The changes allow reverse communication and do not
    !     affect the performance of the routine.
    !
    !     Finds a step which satisfies a sufficient decrease condition
    !     and a curvature condition:
    !
    !         F(X+STP*S) <= F(X) + FTOL*STP*(GRADF(X)'S)                 (sufficient decrease)
    !         |GRADF(X+STP*S)'S| <= GTOL*|GRADF(X)'S|                    (curvature)
    !
    !     using safeguarded cubic/quadratic interpolation (see MCSTEP).
    !
    !     Reverse communication: this routine RETURNs with INFO = -1 to
    !     ask the caller to evaluate F and G at the new X, then must be
    !     called again (with F, G updated and everything else
    !     unchanged) to resume -- exactly like the outer LBFGS_INTERN
    !     driver that calls it.
    !
    !     INFO on return:
    !       0  improper input parameters
    !      -1  a return was made to compute F and G (call again)
    !       1  sufficient decrease and curvature conditions both hold
    !       2  relative width of the interval of uncertainty <= XTOL
    !       3  number of calls to the function has reached MAXFEV
    !       4  the step is at the lower bound STPMIN
    !       5  the step is at the upper bound STPMAX
    !       6  rounding errors prevent further progress
    !
    !     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
    !     JORGE J. MORE', DAVID J. THUENTE
    !--------------------------------------------------------------------
    integer,  intent(in)    :: n, maxfev
    integer,  intent(inout) :: info, nfev
    real(wp), intent(inout) :: x(n), stp
    real(wp), intent(in)    :: f, g(n), s(n), ftol, xtol
    real(wp), intent(inout) :: wa(n)

    real(wp), parameter :: p5 = 0.5_wp, p66 = 0.66_wp, xtrapf = 4.0_wp, zero = 0.0_wp

    integer,  save :: infoc
    logical,  save :: brackt, stage1
    real(wp), save :: dg, dgm, dginit, dgtest, dgx, dgxm, dgy, dgym
    real(wp), save :: finit, ftest1, fm, fx, fxm, fy, fym
    real(wp), save :: stx, sty, stmin, stmax, width, width1

    if (info /= -1) then
       infoc = 1

       ! CHECK THE INPUT PARAMETERS FOR ERRORS.
       if (n <= 0 .or. stp <= zero .or. ftol < zero .or. &
           gtol < zero .or. xtol < zero .or. stpmin < zero .or. &
           stpmax < stpmin .or. maxfev <= 0) return

       ! COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION AND CHECK
       ! THAT S IS A DESCENT DIRECTION.
       dginit = dot_product(g, s)
       if (dginit >= zero) then
          if (lp > 0) write(lp, 15)
15        format(/'  THE SEARCH DIRECTION IS NOT A DESCENT DIRECTION')
          return
       end if

       brackt = .false.
       stage1 = .true.
       nfev   = 0
       finit  = f
       dgtest = ftol * dginit
       width  = stpmax - stpmin
       width1 = width / p5
       wa(1:n) = x(1:n)

       ! THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
       ! FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
       ! THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
       ! FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF THE
       ! INTERVAL OF UNCERTAINTY. THE VARIABLES STP, F, DG CONTAIN THE
       ! VALUES OF THE STEP, FUNCTION, AND DERIVATIVE AT THE CURRENT
       ! STEP.
       stx = zero; fx = finit; dgx = dginit
       sty = zero; fy = finit; dgy = dginit
    else
       goto 45
    end if

    ! START OF ITERATION.
30  continue

    ! SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND TO THE PRESENT
    ! INTERVAL OF UNCERTAINTY.
    if (brackt) then
       stmin = min(stx, sty)
       stmax = max(stx, sty)
    else
       stmin = stx
       stmax = stp + xtrapf * (stp - stx)
    end if

    ! FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
    stp = max(stp, stpmin)
    stp = min(stp, stpmax)

    ! IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET STP BE THE LOWEST
    ! POINT OBTAINED SO FAR.
    if ((brackt .and. (stp <= stmin .or. stp >= stmax)) &
         .or. nfev >= maxfev-1 .or. infoc == 0 &
         .or. (brackt .and. stmax-stmin <= xtol*stmax)) stp = stx

    ! EVALUATE THE FUNCTION AND GRADIENT AT STP AND COMPUTE THE
    ! DIRECTIONAL DERIVATIVE. We return to the caller to obtain F and G.
    x(1:n) = wa(1:n) + stp * s(1:n)
    info = -1
    return

45  continue
    info = 0
    nfev = nfev + 1
    dg = dot_product(g, s)
    ftest1 = finit + stp * dgtest

    ! TEST FOR CONVERGENCE.
    if ((brackt .and. (stp <= stmin .or. stp >= stmax)) &
         .or. infoc == 0) info = 6
    if (stp == stpmax .and. f <= ftest1 .and. dg <= dgtest) info = 5
    if (stp == stpmin .and. (f > ftest1 .or. dg >= dgtest)) info = 4
    if (nfev >= maxfev) info = 3
    if (brackt .and. stmax-stmin <= xtol*stmax) info = 2
    if (f <= ftest1 .and. abs(dg) <= gtol*(-dginit)) info = 1

    ! CHECK FOR TERMINATION.
    if (info /= 0) return

    ! IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED FUNCTION
    ! HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
    if (stage1 .and. f <= ftest1 .and. dg >= min(ftol, gtol)*dginit) stage1 = .false.

    ! A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF WE HAVE
    ! NOT OBTAINED A STEP FOR WHICH THE MODIFIED FUNCTION HAS A
    ! NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE, AND IF A
    ! LOWER FUNCTION VALUE HAS BEEN OBTAINED BUT THE DECREASE IS NOT
    ! SUFFICIENT.
    if (stage1 .and. f <= fx .and. f > ftest1) then
       fm   = f  - stp*dgtest
       fxm  = fx - stx*dgtest
       fym  = fy - sty*dgtest
       dgm  = dg - dgtest
       dgxm = dgx - dgtest
       dgym = dgy - dgtest

       ! CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY AND TO
       ! COMPUTE THE NEW STEP.
       call mcstep(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, &
                    brackt, stmin, stmax, infoc)

       ! RESET THE FUNCTION AND GRADIENT VALUES FOR F.
       fx  = fxm + stx*dgtest
       fy  = fym + sty*dgtest
       dgx = dgxm + dgtest
       dgy = dgym + dgtest
    else
       call mcstep(stx, fx, dgx, sty, fy, dgy, stp, f, dg, &
                    brackt, stmin, stmax, infoc)
    end if

    ! FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE INTERVAL OF
    ! UNCERTAINTY.
    if (brackt) then
       if (abs(sty-stx) >= p66*width1) stp = stx + p5*(sty - stx)
       width1 = width
       width  = abs(sty-stx)
    end if

    ! END OF ITERATION.
    goto 30

  end subroutine mcsrch

  !--------------------------------------------------------------------
  subroutine mcstep(stx, fx, dx, sty, fy, dy, stp, fp, dp, brackt, stpmin, stpmax, info)
    !--------------------------------------------------------------------
    !                          SUBROUTINE MCSTEP
    !
    !     Computes a safeguarded step for a line search and updates an
    !     interval of uncertainty for a minimizer of the function.
    !
    !     STX, FX, DX specify the step, function, and derivative at the
    !     best step obtained so far; DX must be negative in the
    !     direction of the step. STY, FY, DY specify the same at the
    !     other endpoint of the interval of uncertainty. STP, FP, DP
    !     specify the current step. BRACKT is true once a minimizer has
    !     been bracketed between STX and STY. STPMIN/STPMAX bound the
    !     step. INFO is set to 1..4 according to which of the four cases
    !     below was used, or left 0 for improper input.
    !
    !     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
    !     JORGE J. MORE', DAVID J. THUENTE
    !--------------------------------------------------------------------
    real(wp), intent(inout) :: stx, fx, dx, sty, fy, dy, stp
    real(wp), intent(in)    :: fp, dp, stpmin, stpmax
    logical,  intent(inout) :: brackt
    integer,  intent(out)   :: info

    logical  :: bound
    real(wp) :: gamma, p, q, r, s, sgnd, stpc, stpf, stpq, theta

    info = 0

    ! CHECK THE INPUT PARAMETERS FOR ERRORS.
    if ((brackt .and. (stp <= min(stx,sty) .or. stp >= max(stx,sty))) &
         .or. dx*(stp-stx) >= 0.0_wp .or. stpmax < stpmin) return

    ! DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
    sgnd = dp * (dx/abs(dx))

    if (fp > fx) then
       ! FIRST CASE. A HIGHER FUNCTION VALUE. THE MINIMUM IS BRACKETED.
       ! IF THE CUBIC STEP IS CLOSER TO STX THAN THE QUADRATIC STEP,
       ! THE CUBIC STEP IS TAKEN, ELSE THE AVERAGE OF THE CUBIC AND
       ! QUADRATIC STEPS IS TAKEN.
       info  = 1
       bound = .true.
       theta = 3*(fx - fp)/(stp - stx) + dx + dp
       s     = max(abs(theta), abs(dx), abs(dp))
       gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
       if (stp < stx) gamma = -gamma
       p = (gamma - dx) + theta
       q = ((gamma - dx) + gamma) + dp
       r = p/q
       stpc = stx + r*(stp - stx)
       stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp - stx)
       if (abs(stpc-stx) < abs(stpq-stx)) then
          stpf = stpc
       else
          stpf = stpc + (stpq - stpc)/2
       end if
       brackt = .true.

    else if (sgnd < 0.0_wp) then
       ! SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
       ! OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS
       ! CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP, THE CUBIC
       ! STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
       info  = 2
       bound = .false.
       theta = 3*(fx - fp)/(stp - stx) + dx + dp
       s     = max(abs(theta), abs(dx), abs(dp))
       gamma = s*sqrt((theta/s)**2 - (dx/s)*(dp/s))
       if (stp > stx) gamma = -gamma
       p = (gamma - dp) + theta
       q = ((gamma - dp) + gamma) + dx
       r = p/q
       stpc = stp + r*(stx - stp)
       stpq = stp + (dp/(dp-dx))*(stx - stp)
       if (abs(stpc-stp) > abs(stpq-stp)) then
          stpf = stpc
       else
          stpf = stpq
       end if
       brackt = .true.

    else if (abs(dp) < abs(dx)) then
       ! THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE SAME
       ! SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES. THE
       ! CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY IN THE
       ! DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC IS
       ! BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE EITHER
       ! STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
       ! COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE STEP
       ! CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS
       ! TAKEN.
       info  = 3
       bound = .true.
       theta = 3*(fx - fp)/(stp - stx) + dx + dp
       s     = max(abs(theta), abs(dx), abs(dp))

       ! THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND TO
       ! INFINITY IN THE DIRECTION OF THE STEP.
       gamma = s*sqrt(max(0.0_wp, (theta/s)**2 - (dx/s)*(dp/s)))
       if (stp > stx) gamma = -gamma
       p = (gamma - dp) + theta
       q = (gamma + (dx - dp)) + gamma
       r = p/q
       if (r < 0.0_wp .and. gamma /= 0.0_wp) then
          stpc = stp + r*(stx - stp)
       else if (stp > stx) then
          stpc = stpmax
       else
          stpc = stpmin
       end if
       stpq = stp + (dp/(dp-dx))*(stx - stp)
       if (brackt) then
          if (abs(stp-stpc) < abs(stp-stpq)) then
             stpf = stpc
          else
             stpf = stpq
          end if
       else
          if (abs(stp-stpc) > abs(stp-stpq)) then
             stpf = stpc
          else
             stpf = stpq
          end if
       end if

    else
       ! FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE SAME
       ! SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES NOT DECREASE.
       ! IF THE MINIMUM IS NOT BRACKETED, THE STEP IS EITHER STPMIN OR
       ! STPMAX, ELSE THE CUBIC STEP IS TAKEN.
       info  = 4
       bound = .false.
       if (brackt) then
          theta = 3*(fp - fy)/(sty - stp) + dy + dp
          s     = max(abs(theta), abs(dy), abs(dp))
          gamma = s*sqrt((theta/s)**2 - (dy/s)*(dp/s))
          if (stp > sty) gamma = -gamma
          p = (gamma - dp) + theta
          q = ((gamma - dp) + gamma) + dy
          r = p/q
          stpc = stp + r*(sty - stp)
          stpf = stpc
       else if (stp > stx) then
          stpf = stpmax
       else
          stpf = stpmin
       end if
    end if

    ! UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT DEPEND
    ! ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
    if (fp > fx) then
       sty = stp
       fy  = fp
       dy  = dp
    else
       if (sgnd < 0.0_wp) then
          sty = stx
          fy  = fx
          dy  = dx
       end if
       stx = stp
       fx  = fp
       dx  = dp
    end if

    ! COMPUTE THE NEW STEP AND SAFEGUARD IT.
    stpf = min(stpmax, stpf)
    stpf = max(stpmin, stpf)
    stp  = stpf
    if (brackt .and. bound) then
       if (sty > stx) then
          stp = min(stx+0.66_wp*(sty-stx), stp)
       else
          stp = max(stx+0.66_wp*(sty-stx), stp)
       end if
    end if

  end subroutine mcstep

end module m_lbfgs_intern
