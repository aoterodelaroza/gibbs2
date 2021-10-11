! Copyright (c) 2011 Alberto Otero de la Roza <aoterodelaroza@gmail.com>,
! Víctor Luaña <victor@carbono.quimica.uniovi.es> and David
! Abbasi <david@carbono.quimica.uniovi.es>. Universidad de Oviedo.
!
! gibbs2 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or (at
! your option) any later version.
!
! gibbs2 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module fit
  implicit none
  private

  public :: mpar, mparmin, mmpar, ndel
  public :: fit_init, fit_pshift
  public :: fitt_polygibbs, fit_ev
  public :: polfit

  public :: fitinfo

  integer, parameter :: mmpar = 20
  integer :: mpar, mparmin, ndel
  integer, parameter :: mdata = 15000

  type fitinfo
     integer :: nfit = 0
     integer :: npar(mmpar) = 0
     integer :: mode(mmpar) = 0
     real*8 :: wei(mmpar) = 0d0
     real*8 :: apar(0:mmpar,mmpar) = 0d0
     real*8 :: veq(mmpar) = 0d0
     real*8 :: beq(mmpar) = 0d0
  end type fitinfo

contains

  ! Initialize the variables in this routine
  subroutine fit_init()
    use evfunc, only: pfit_mode, pfit_slatec, pweigh_mode, pweigh_gibbs2
    pfit_mode = pfit_slatec
    pweigh_mode = pweigh_gibbs2
    !
    mpar = 12
    mparmin = 2
    ndel = 3

  end subroutine fit_init

  ! Given a volume grid (v, au), an energy (npol, cpol) with some
  ! fitting mode (mode), and a pressure (p_in, GPa), calculate the
  ! volume that corresponds to that pressure (vx, au), the bulk
  ! modulus (bx, GPa), and the value of the energy (ex, Ha) and
  ! enthalpy (hx, Ha) at that point. Return ierr /= 0 if the fit
  ! fails.
  subroutine fit_pshift(mode,v,p_in,npol,cpol,vx,bx,ex,hx,ierr)
    use param, only: au2gpa, faterr
    use tools, only: error
    use evfunc, only: fv0, fv1, fv2
    integer, intent(in) :: mode
    real*8, intent(in) :: v(:)
    real*8, intent(in) :: p_in
    integer, intent(in) :: npol
    real*8, intent(in) :: cpol(0:npol)
    real*8, intent(out) :: vx, bx, ex, hx
    integer, intent(out) :: ierr

    real*8, parameter :: tol = 1d-10
    real*8, parameter :: tolv = 1d-7
    integer, parameter :: mstep = 100

    integer :: i, istep, ndat
    real*8 :: p1, p2, p, px, ppx
    logical :: found
    real*8 :: v1, v2, vxx
    integer :: nstep

    ierr = 0
    p = p_in / au2gpa

    ! bracket the p range
    ndat = size(v)
    found = .false.
    a: do istep = 1, ndat/2
       p2 = -fv1(mode,v(1),npol,cpol)
       do i = istep+1, ndat, istep
          p1 = p2
          p2 = -fv1(mode,v(i),npol,cpol)
          found = (p <= p1) .and. (p2 <= p)
          if (found) then
             v1 = v(i-istep)
             v2 = v(i)
             p1 = p1 - p
             p2 = p2 - p
             exit a
          end if
       end do
    end do a
    if (.not. found) then
       ierr = 2
       if (fv0(mode,v(1),npol,cpol) < fv0(mode,v(ndat),npol,cpol)) then
          vx = v(1)
       else
          vx = v(ndat)
       end if
       bx = vx * fv2(mode,vx,npol,cpol) * au2gpa
       ex = fv0(mode,vx,npol,cpol)
       hx = (ex + p * vx)
       return
    end if

    ! find p(V) = p
    vx = 0.5d0*(v1+v2)
    px = 1d30
    ierr = 1
    do nstep = 1, mstep
       ! funcall
       px = -fv1(mode,vx,npol,cpol)-p
       if (abs(px) < tol .or. abs(v1-v2) < tolv) then
          ierr = 0
          exit
       end if

       ! update interval
       if (px * p1 > 0) then
          v1 = vx
       else
          v2 = vx
       end if

       ! step
       ppx = -fv2(mode,vx,npol,cpol)
       vxx = vx - px / ppx
       if (vxx >= v1 .and. vxx <= v2) then
          vx = vxx
       else
          vx = 0.5d0*(v1+v2)
       end if
    end do

    bx = vx * fv2(mode,vx,npol,cpol) * au2gpa
    ex = fv0(mode,vx,npol,cpol)
    hx = (ex + p * vx)

    if (ierr == 1) then
       call error('fit_pshift','Error finding minimum',faterr)
    end if

  end subroutine fit_pshift

  ! Wrapper routine for different E(V) fit equations.
  subroutine fit_ev(fmode, rmode, var, func, nparpro, aparpro, ierrout, ispv,&
     nfix, idfix, obelix, pfit)
    use param, only: au2gpa, warning
    use tools, only: error
    use evfunc, only: fit_polygibbs
    integer, intent(in) :: fmode, rmode
    real*8, intent(in) :: var(:), func(:)
    integer, intent(out) :: nparpro
    real*8, intent(out) :: aparpro(0:mmpar)
    integer, intent(out) :: ierrout
    logical, intent(in) :: ispv
    integer, intent(in), optional :: nfix, idfix(0:mmpar)
    real*8, intent(in), optional :: obelix(0:mmpar)
    type(fitinfo), intent(out), optional :: pfit

    integer :: i, i1, idx, idxx(1)
    real*8 :: mslope, rfunc(size(func))
    logical :: dopshift

    idxx = minloc(func)
    idx = idxx(1)
    dopshift = (idx <= 2 .or. idx >= size(func)-1)

    ! initialize
    nparpro = 0
    aparpro = 0d0

    ! obtain a +pV scaling that shifts the minimum to the center of the grid
    if (dopshift) then
       i1 = size(var) / 2
       mslope = - (func(i1+1) - func(i1)) / (var(i1+1) - var(i1))
    else
       mslope = 0d0
    end if

    rfunc = func + mslope * var

    if (fmode == fit_polygibbs .or. fmode > 100) then
       call fitt_polygibbs(fmode, var, rfunc, nparpro, aparpro, ierrout, ispv, pfit)

       if (ispv .and. ierrout == 0 .and. fmode == fit_polygibbs) then
          if (abs(aparpro(0)) * au2gpa > 2d0) then
             call error('fit_ev','fitted p_0 > 2GPa',warning)
          end if

          ! integrate
          aparpro(nparpro+1) = aparpro(nparpro)
          do i = nparpro, 1, -1
             aparpro(i) = aparpro(i-1) / real(i,8) * aparpro(nparpro+1)
          end do
          aparpro(0) = 0d0
          aparpro(0:nparpro) = -aparpro(0:nparpro)
          nparpro = nparpro + 1
       end if
    else
       call fitt_eos(fmode, rmode, var, rfunc, nparpro, aparpro, ierrout, ispv,&
          nfix, idfix, obelix)
    end if

    nparpro = nparpro + 1
    aparpro(nparpro) = mslope
    if (present(pfit)) then
       if (pfit%nfit > 0) then
          pfit%npar(1:pfit%nfit) = pfit%npar(1:pfit%nfit) + 1
          do i = 1, pfit%nfit
             pfit%apar(pfit%npar(i),i) = mslope
          end do
       end if
    end if

  end subroutine fit_ev

  ! Fits polynomials to (var,function) data and averages them,
  ! weighted by its chi-square test probabilities. It returns the
  ! averaged polynomial coefficients in aparpro, and the coefficients
  ! of its square in a2parpro.
  subroutine fitt_polygibbs (mode, var, func, nparpro, aparpro, ierrout, ispv, pfit)
    use param, only: warning
    use tools, only: error
    use evfunc, only: fit_strain, pweigh_gibbs1, pweigh_gibbs2, pweigh_mode, pweigh_slatec,&
       v2str, fv1
    integer, intent(in) :: mode
    real*8, intent(in) :: var(:), func(:)
    integer, intent(out) :: nparpro
    real*8, intent(out) :: aparpro(0:mmpar)
    integer, intent(out) :: ierrout
    logical, intent(in) :: ispv
    type(fitinfo), intent(out), optional :: pfit

    integer :: imin, iamin(1), ndata
    real*8  :: apar(0:mmpar), rms
    real*8, allocatable :: w(:)
    integer :: npar, nfit
    integer :: mmfit
    integer, allocatable :: nparfit(:), ndatafit(:)
    integer :: npar2(mpar+1)
    real*8 :: rms2(mpar+1), pwei(mpar+1), rms2min, apar2(0:mmpar,mpar+1)
    integer :: nparmax, i, ifit
    real*8  :: wnorm, wtmp, rmsmin, eps
    integer, parameter :: limit = 4
    integer :: nparmin, ndatamin, ndataact, iinf, isup
    real*8 :: tmp
    integer :: ierr
    real*8 :: r(size(var)), a(3*(size(var)+mmpar)+3)
    integer :: idxx(1), idx
    real*8 :: x(size(var)), vscal
    integer :: ndel_now
    integer :: strain, order

    ! initialize
    ndata = size(var)
    ierrout = 0
    nparpro = 0
    aparpro = 0d0
    allocate(w(mdata))

    idxx = minloc(func)
    idx = idxx(1)
    vscal = var(idx)
    if (mode > 100) then
       strain = (mode - fit_strain * 10000) / 100
       order =  (mode - fit_strain * 10000 - 100 * strain)
       x = v2str(strain,var,vscal)
    else
       ! scale the volumes
       x = var / vscal
       strain = 0
    end if

    if (strain > 0 .and. order > 0) then
       w = 1d0
       npar = order
       call polfit(ndata,1,ndata,x,func,w,rms,npar,apar)
       nparpro = npar
       aparpro = apar
    else
       if (pweigh_mode == pweigh_gibbs2) then
          ! save the minimum of the energy
          iamin = minloc(func)
          imin = iamin(1)

          !.initialize the weights for the polynomial fit
          w(1:ndata) = 1d0

          !.fitting loops through different polynomials:
          nfit = 0
          ! nparmax = min(ndata-5,mpar)
          nparmax = max(min(ndata/2-1,mpar),3)
          mmfit = mpar + 1
          rms2min = 1d30
          apar2 = 0d0
          do npar = mparmin, nparmax
             ! fit this polynomial
             nfit = nfit + 1
             call polfit (ndata,1,ndata,x,func,w,rms,npar,apar2(:,nfit))

             !.discard fitts that don't give a minimum in the input bracket
             if (nfit.gt.mmfit) then
                nfit = mmfit
                call error('fitt_polygibbs','max number of fitts exceeded',warning)
             else
                ! save order, rms
                npar2(nfit) = npar
                rms2(nfit) = rms * rms
                rms2min = min(rms2min,rms2(nfit))
             endif
          end do

          ! check at least one polynomial has been fitted
          if (nfit == 0) then
             nparpro = 0
             aparpro = 0d0
             ierrout = 1
             return
          end if

          ! polynomial weights
          pwei(1:nfit) = exp(-rms2(1:nfit) / (rms2min+1d-15) * real(npar2(1:nfit),8) / real(ndata,8))
          pwei = pwei / sum(pwei(1:nfit))


          ! average polynomial
          nparpro = nparmax
          aparpro(1:nparmax) = 0d0
          do i = 1, nfit
             aparpro(0:npar2(i)) = aparpro(0:npar2(i)) + pwei(i) * apar2(0:npar2(i),i)
          end do

          ! save partial fit information
          if (present(pfit)) then
             pfit%nfit = nfit
             pfit%npar(1:nfit) = npar2(1:nfit)
             pfit%mode(1:nfit) = fit_strain * 10000 + strain * 100 + pfit%npar(1:nfit)
             pfit%apar(:,1:nfit) = apar2(:,1:nfit)
             pfit%wei(1:nfit) = pwei(1:nfit)
          end if

       elseif (pweigh_mode == pweigh_gibbs1) then
          iamin = minloc(func)
          imin = iamin(1)
          if (ispv) then
             ndel_now = 0
          else
             ndel_now = ndel
          end if

          !.initialize the weights
          w(1:ndata) = 1d0

          !.fitting loops through different polynomials:
          nfit = 0
          nparpro = 0
          nparmin = 3
          iinf = 1
          isup = ndata
          ndatamin = max (ndata-2*ndel_now, nparmin+1)
          ndataact = ndata
          rmsmin = 1d33
          mmfit = ((ndataact-ndatamin+1)/2 + 1) * (mpar - nparmin + 1) + 1
          allocate(nparfit(mmfit),ndatafit(mmfit))

          !.loop eliminating pairs of outermost elements ndel times:
          do while (ndataact.ge.ndatamin .and. (iinf.lt.imin .and. isup.gt.imin.or.ispv))
             !.loop increasing the number of parameters until limit is reached:
             nparmax = min (ndataact-limit-1, mpar)
             if (ispv) nparmax = nparmax - 1
             do npar = nparmin, nparmax
                nfit = nfit + 1
                call polfit (ndata,iinf,isup,x,func,w,rms,npar,apar)
                apar(npar+1) = 1d0
                apar(npar+2) = 0d0

                !.discard fitts that don't give a minimum in the input bracket
                tmp = fv1(mode,x(imin-1),npar+2,apar) * fv1(mode,x(imin+1),npar+2,apar)
                if (tmp.ge.0d0 .and..not.ispv) then
                   nfit = nfit - 1

                   !.discard fitts over the maximum declared dimensions
                else if (nfit.gt.mmfit) then
                   nfit = mmfit
                   call error('fitt_polygibbs','max number of fitts exceeded',warning)
                   !.save fitt parameters
                else
                   nparpro = max (nparpro,npar)
                   nparfit(nfit) = npar
                   ndatafit(nfit) = ndataact
                   rmsmin = min(rmsmin,rms*npar/ndataact)
                endif
             enddo
             iinf = iinf + 1
             isup = isup - 1
             ndataact = ndataact - 2
          enddo

          !.number of fitts control
          if (nfit.eq.0) then
             nparpro = 0
             aparpro = 0d0
             ierrout = 1
             return
          endif

          !.average the polynomial coefficients (repeating the fitts)
          wnorm = 0d0
          aparpro = 0d0
          do ifit = 1, nfit
             ndataact = ndatafit(ifit)
             npar = nparfit(ifit)
             iinf = 1 + (ndata-ndataact) / 2
             isup = ndata - iinf + 1
             call polfit (ndata,iinf,isup,x,func,w,rms,npar,apar)
             wtmp = rms*npar/(rmsmin*ndataact)
             wtmp = exp(-wtmp*wtmp)
             do i = 0, npar
                aparpro(i) = aparpro(i) + wtmp * apar(i)
             enddo
             wnorm = wnorm + wtmp
          enddo

          !.put the proper weight into the polynomials
          aparpro = aparpro / wnorm

          deallocate(nparfit,ndatafit)

       else if (pweigh_mode == pweigh_slatec) then
          eps = -1d0
          w = 1d0
          nparmax = min(ndata-limit-1, mpar)
          if (ispv) nparmax = nparmax - 1
          call dpolft(ndata,x,func,w(1:ndata),nparmax,nparpro,eps,r,ierr,a)
          call dpcoef(nparpro,0d0,aparpro,a)
       end if
    end if

    nparpro = nparpro + 1
    aparpro(nparpro) = vscal
    if (present(pfit)) then
       if (pfit%nfit > 0) then
          pfit%npar(1:pfit%nfit) = pfit%npar(1:pfit%nfit) + 1
          do i = 1, pfit%nfit
             pfit%apar(pfit%npar(i),i) = vscal
          end do
       end if
    end if

  end subroutine fitt_polygibbs

  ! Fit equation of state of type fmode and regression mode rmode to data var (x)
  ! and func (x). Returns the fitted parameters aparpro(0:nparpro). ierrout = 0
  ! if success, non-zero if failure. ispv = .true. if these are pv data. nfix,
  ! idfix, obelix - fix some EOS parameters.
  subroutine fitt_eos(fmode,rmode,var,func,nparpro,aparpro,ierrout,ispv,nfix,idfix,obelix)
    use evfunc, only: evfunc_xm, evfunc_ym, evfunc_mask, evfunc_fixval, evfunc_domask,&
       evfunc_minpack_mode, evfunc_npar, evfunc_reg_mode, fcn_minpack, fcn_minpack1,&
       fit_antons, fit_order
    use tools, only: error
    use param, only: faterr
    integer, intent(in) :: fmode, rmode
    real*8, intent(in) :: var(:), func(:)
    integer, intent(out) :: nparpro
    real*8, intent(out) :: aparpro(0:mmpar)
    integer, intent(out) :: ierrout
    logical, intent(in) :: ispv
    integer, intent(in), optional :: nfix, idfix(0:mmpar)
    real*8, intent(in), optional :: obelix(0:mmpar)

    real*8, parameter :: minpack_tol = 1d-8

    integer :: i, count
    integer :: idxx(1), idx, m, iwa(7)
    real*8 :: ffit(size(func))
    real*8 :: wa(size(var)*(7+1) + 7*5), atemp(0:mmpar)
    integer :: lwa

    ! output order
    nparpro = fit_order(fmode)
    aparpro = 0d0

    ! fill in data
    m = size(var)
    evfunc_minpack_mode = fmode
    evfunc_reg_mode = rmode
    allocate(evfunc_xm(m),evfunc_ym(m))
    evfunc_xm = var
    evfunc_ym = func

    ! initial parameters
    lwa = size(var)*(nparpro+1) + 5*(nparpro+1)
    idxx = minloc(func)
    idx = idxx(1)
    if (ispv) then
       aparpro(1) = 0d0
       aparpro(2) = var(m)
       aparpro(3) = (func(m-2)+func(m)-2*func(m-1)) / &
          (var(m-2)-var(m-1)) * (var(m-1)-var(m))
    else
       if (fmode == fit_antons) then
          aparpro(1) = func(size(func))
       else
          aparpro(1) = minval(func)
       end if
       aparpro(2) = var(idx)
       if (idx /= 1 .and. idx /= m) then
          aparpro(3) = var(idx) * (func(idx+1)+func(idx-1)-2*func(idx)) / &
             (var(idx+1)-var(idx)) * (var(idx)-var(idx-1))
       else
          aparpro = 0d0
          ierrout = 1
          deallocate(evfunc_xm,evfunc_ym)
          return
       end if
    end if
    if (nparpro > 3) aparpro(4) = 4d0
    if (nparpro > 4) aparpro(5) = -5d-2
    if (nparpro > 5) aparpro(6) = 5d-4

    ! optional fixing of some parameters
    evfunc_npar = nparpro
    if (present(nfix) .and. present(idfix) .and. present(obelix)) then
       evfunc_domask = .true.
       evfunc_mask = 0
       evfunc_fixval = 0d0
       do i = 1, nfix
          if (idfix(i) < 2 .or. idfix(i) > nparpro+i-1) &
             call error('fitt_eos','Fixed parameter not allowed (is it V0?)',faterr)
          evfunc_mask(idfix(i)) = nparpro
          evfunc_fixval(idfix(i)) = obelix(i)
          nparpro = nparpro - 1
       end do
       count = 0
       do i = 1, nparpro + nfix
          if (evfunc_mask(i) == 0) then
             count = count + 1
             evfunc_mask(i) = count
          end if
       end do
    else
       evfunc_domask = .false.
    end if

    ! minpack (approximate jacobian)
    if (ispv) then
       call lmdif1(fcn_minpack1,m,nparpro,aparpro(1:nparpro),ffit,minpack_tol,&
          ierrout,iwa(1:nparpro),wa(1:lwa),lwa)
    else
       call lmdif1(fcn_minpack,m,nparpro,aparpro(1:nparpro),ffit,minpack_tol,&
          ierrout,iwa(1:nparpro),wa(1:lwa),lwa)
    end if

    if (ierrout >= 1 .and. ierrout <= 3 .or. ierrout >= 6 .and. ierrout <= 7) then
       ierrout = 0
    else
       ierrout = 1
    end if

    ! undo the masking
    if (evfunc_domask) then
       atemp = aparpro
       do i = 1, nfix
          aparpro(idfix(i)) = obelix(i)
       end do
       nparpro = evfunc_npar
       do i = 1, nparpro
          if (evfunc_mask(i) > nparpro-nfix) cycle
          aparpro(i) = atemp(evfunc_mask(i))
       end do
    end if

    deallocate(evfunc_xm,evfunc_ym)

  end subroutine fitt_eos

  ! Fit the (x,y) ndata pairs to a polynomial:
  !   y(x) = SUM(i=0,npar) a(i) * x**i
  ! Only the points from iinf to isup are used.
  !
  ! a(i) are linear parameters obtained by a least squares (absolute
  ! deviation) technique, and w(i) are the weights of each point (if
  ! the weights are the inverse of the variance at each data point,
  ! then rms would be the square root of the "chi square" estimator
  ! divided by the sum of errors)
  !
  subroutine polfit (ndata, iinf, isup, x, y, w, rms, npar, apar)
    use param, only: warning
    use tools, only: error, gauss
    use evfunc, only: fit_polygibbs, pfit_gauss, pfit_mode, fv0
    integer, intent(in) :: ndata, iinf, isup
    real*8, intent(in) :: x(ndata), y(ndata), w(ndata)
    integer, intent(in) :: npar
    real*8, intent(out) :: apar(0:npar), rms

    real*8 :: c(0:mmpar,0:mmpar+1), wnorm, apar2(0:npar+1)
    integer :: k, j, i, ij, ierr
    real*8 :: s2
    integer :: nfit, nparout
    real*8 :: eps, a(3*(ndata+npar)+3), r(ndata)

    if (pfit_mode == pfit_gauss) then
       !.calculate the norm of the weights
       wnorm = 0d0
       do k = iinf, isup
          wnorm = wnorm + w(k)
       enddo
       wnorm = 1d0 / wnorm

       !.construct the linear system matrix c:
       do j = 0, npar
          c(j,npar+1) = 0d0
          if (j.gt.0) then
             do k = iinf, isup
                c(j,npar+1) = c(j,npar+1) + w(k) * y(k) * x(k)**j
             enddo
          else
             do k = iinf, isup
                c(j,npar+1) = c(j,npar+1) + w(k) * y(k)
             enddo
          endif
          c(j,npar+1) = wnorm * c(j,npar+1)
          do i = j, npar
             c(i,j) = 0d0
             ij = i + j
             if (ij .gt. 0) then
                do k = iinf, isup
                   c(i,j) = c(i,j) + w(k) * x(k)**ij
                enddo
                c(i,j) = wnorm * c(i,j)
             else
                c(i,j) = 1d0
             endif
             c(j,i) = c(i,j)
          enddo
       enddo

       !.Solve the linear system for the best A()'s:
       call gauss (c, npar, mpar, apar, ierr)

       !.compute the rms deviation:
       apar2(0:npar) = apar
       apar2(npar+1) = 1d0
       apar2(npar+2) = 0d0
       s2 = 0d0
       do k = iinf, isup
          s2 = s2 + w(k) * (y(k)-fv0(fit_polygibbs,x(k),npar+2,apar2))**2
       enddo
       rms = sqrt(s2*wnorm)
    else
       nfit = isup - iinf + 1
       nparout = npar
       eps = 0d0

       call dpolft(nfit,x(iinf:isup),y(iinf:isup),w(iinf:isup),npar,nparout,&
          eps,r,ierr,a)
       call dpcoef(nparout,0d0,apar,a)
       rms = abs(eps)

       if (ierr /= 1) then
          call error('polfit_slatec','polynomial fitting error',warning)
       end if
    end if

!    if (rmode == reg_lad) then
!       lwa = ndata*(npar+1) + 5*(npar+1)
!       evfunc_minpack_mode = fit_polygibbs
!       evfunc_reg_mode = reg_lad
!       evfunc_xm(1:isup-iinf+1) = x(iinf:isup)
!       evfunc_ym(1:isup-iinf+1) = y(iinf:isup)
!       ! minpack (approximate jacobian)
!       call lmdif1(fcn_minpack,isup-iinf+1,npar+1,apar,ffit,1d-6,&
!          ierr,iwa(1:npar+1),wa(1:lwa),lwa)
!       if (ierr >= 1 .and. ierr <= 3) ierr = 0
!       !.compute the rms deviation:
!       apar2(0:npar) = apar
!       apar2(npar+1) = 1d0
!       s2 = 0d0
!       do k = iinf, isup
!          s2 = s2 + (y(k)-fv0(fit_polygibbs,x(k),npar+1,apar2))**2
!       enddo
!       rms = sqrt(s2)
!    end if

  end subroutine polfit

end module fit
