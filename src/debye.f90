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

module debye
  implicit none
  private

  ! public
  public :: fill_thetad
  public :: thermal, debeins, thermal_qha, thermal_debye_extended

  ! if temperature is lower, then gamma(t) = gamma(tlim)
  real*8, parameter, public :: tlim_gamma = 1d0 + 1d-6
  real*8, parameter, public :: cvlim = 1d-30

contains

  ! Fill Debye-related information for phase p: %td(:).  Optionally, verbose output.
  subroutine fill_thetad(p,verbose)
    use evfunc, only: fv0, fv2, fv3
    use fit, only: fitt_polygibbs
    use tools, only: leng, realloc, error
    use varbas, only: phase, tm_debye_input, tm_debye_poisson_input, mm, vfree, &
       phase_realloc_volume, tm_debyegrun
    use param, only: mline, mline_fmt, faterr, ifmt_v, ifmt_t, pckbau, pi, third, uout, warning,&
       format_string, format_string_header, au2gpa
    type(phase), intent(inout) :: p
    logical, intent(in) :: verbose

    integer :: j
    character*(mline) :: msg
    integer :: idx, idxx(1)
    integer :: ierr
    character*(mline_fmt) :: fm
    real*8 :: gamma, td, td0, f2s, f3s, b, v
    real*8 :: fx, gx, hx, poi, pofunc

    if (mm < 0d0 .or. vfree < 0d0) return

    ! remove non-convex region
    idxx = minloc(p%e)
    idx = idxx(1)
    do j = max(idx-1,1), p%nv
       f2s = fv2(p%fit_mode,p%v(j),p%npol,p%cpol)
       if (f2s < 0) then
          write (msg,&
             '("Removed ",I3," high-V points where Epp<0 , phase ",A)')&
             p%nv-j+1, trim(adjustl(p%name(1:leng(p%name))))
          call error('gibbs2',msg,warning)

          call phase_realloc_volume(p,j-1)
          p%nv = j-1
          exit
       end if
    end do

    if (.not.allocated(p%td)) allocate(p%td(p%nv))

    ! If the model is debye with poisson, calculate thetad(V) here
    if (p%tmodel == tm_debye_poisson_input) then
       do j = 1, p%nv
          poi = p%poissonv(j)
          fx=2*(1+poi)/3d0/(1-2*poi)
          gx=(1+poi)/3d0/(1-poi)
          hx=2d0*sqrt(fx**3)+sqrt(gx**3)
          pofunc=exp(-log(hx/3)/3)
          f2s = fv2(p%fit_mode,p%v(j),p%npol,p%cpol)
          p%td(j) = (6*pi*pi*vfree*p%v(j)*p%v(j))**third / pckbau * pofunc * sqrt(f2s/mm)
       end do
    end if

    ! td at the equilibrium volume
    p%td0 = (6*pi*pi*vfree*p%veq_static*p%veq_static)**third / pckbau * p%pofunc * sqrt(f2s/mm)

    ! header
    if (verbose) then
       write (uout,'("+ Phase ",A)') trim(adjustl(p%name(1:leng(p%name))))
       write (uout,'("# ThetaD at static eq. volume: ",F10.2)') p%td0
       fm = format_string_header((/1,ifmt_v,ifmt_t,ifmt_t/),(/1,9,9,17/))
       write (uout,fm) "#","V(bohr^3)","Tdebye(K)","Tdebye_slater(K)"
       fm = format_string((/ifmt_v,ifmt_t,ifmt_t/),1)
    end if

    ! get thetad and output
    do j = 1, p%nv
       v = p%v(j)
       f2s = fv2(p%fit_mode,v,p%npol,p%cpol)
       if (f2s < 0d0) then
          p%dyn_active(j) = .false.
          write (uout,'(" ",F10.4," deactivated because E'''' < 0")') p%v(j)
          cycle
       end if
       f3s = fv3(p%fit_mode,v,p%npol,p%cpol)

       if (p%tmodel == tm_debyegrun) then
          b = v * f2s * au2gpa
          p%td(j) = p%td0 * (b / p%beq_static)**p%b_grun / (v / p%veq_static)**p%a_grun
       elseif (p%tmodel /= tm_debye_input .and. p%tmodel /= tm_debye_poisson_input) then
          p%td(j) = (6*pi*pi*vfree*v*v)**third / pckbau * p%pofunc * sqrt(f2s/mm)
       end if

       if (verbose) then
          td0 = (6*pi*pi*vfree*p%v(j)*p%v(j))**third / pckbau * p%pofunc * sqrt(f2s/mm)
          write (uout,fm) p%v(j), p%td(j), td0
       end if
    end do

  end subroutine fill_thetad

  ! Compute Debye model vibrational properties.
  subroutine thermal(ThetaD,T,debye,xabs,en,cv,he,ent)
    use param, only: faterr, pckbau, pi, zero
    use tools, only: error, gauleg
    use varbas, only: vfree
    !
    ! This routine obtains the molar vibrational properties of a given
    ! crystal by means of the Debye model: internal energy (U), heat
    ! capacity at constant volume (Cv), Helmholtz's free energy (F),
    ! and vibrational entropy (S).
    !
    ! To evaluate this properties, the following integral is needed:
    !
    !                                  |    x^3     |
    ! Debye (y) = 3*y^(-3) * INT (0,y) | ---------- | dx
    !                                  | exp(x) - 1 |
    !
    ! where y=ThetaD/T, being ThetaD Debye's temperature (K) and T the
    ! absolute (thermodynamic) temperature. The integral is evaluated
    ! using a Gauss-Legendre quadrature.
    !
    !-INPUT-------------------------------------------------------------
    !   ThetaD : Debye's temperature (K).
    !        T : Absolute temperature (K).
    !-OUTPUT-------------------------------------------------------------
    !       en : Vibrational internal energy, U (hartree/molecule).
    !       cv : Constant V heat capacity, Cv (hartree/K molecule).
    !       he : Helmholtz's free energy (hartree/molecule).
    !      ent : Entropy (hartree/K molecule).
    !    Debye : Debye's integral.
    !     xabs : Maximum error in Debye's integral evaluation.
    !------------------------------------------------------------------------

    real*8, intent(in) :: ThetaD, T
    real*8, intent(out) :: debye, xabs, en, cv, he, ent

    real*8, parameter :: eps=1D-12
    integer, parameter :: maxnl=100
    real*8 :: x(maxnl), w(maxnl), y, z, sum, debye0, debyeout
    integer :: i, nl

    !.error condition controls
    debyeout = 0d0
    xabs = 0d0
    if (t<1d-5) then
       en = vfree*9d0*pckbau*ThetaD/8d0
       cv = 0d0
       he = vfree*9d0*pckbau*ThetaD/8d0
       ent = 0d0
       debye = 0d0
       return
    endif
    if (ThetaD.eq.0d0) then
       en = vfree*3d0*pckbau*T
       cv = vfree*3d0*pckbau
       ent = 1d100
       he = en - T*ent
       debye = 0d0
       return
    endif
    if (thetad.lt.0d0 .or. t.lt.0d0) then
       call error('thermal','Negative ThetaD or T',faterr)
    endif
    y=ThetaD/T
    debyeout=3d0*pi*pi*pi*pi/y/y/y/15d0
    if (y.gt.250d0) goto 22

    !.Loop with increasing number of Legendre points.
    debye0=1d30
    do nl=5,maxnl,5
       call gauleg (0D0,y,x,w,nl)
       sum=0d0
       do i=1,nl
          sum=sum+w(i)*fdebye(x(i))
       end do
       debyeout=sum*3d0/y/y/y
       xabs = abs(debye-debye0)
       if (xabs.lt.eps) then
          exit
       else
          debye0 = debyeout
       endif
    end do

    !.thermodynamic vibrational properties
22  en  = vfree * 3d0 * pckbau * (ThetaD*3d0/8d0 + T*debyeout)
    if (y > 100d0) then
       cv  = vfree * 12d0 * pckbau * debyeout
    else
       cv  = vfree * 3d0 * pckbau * (4d0*debyeout - 3d0*y/(exp(y)-1d0))
    endif
    ent = vfree * 3d0 * pckbau * (debyeout*4d0/3d0 - log(1d0-exp(-y)))
    he  = en - T * ent
    debye = debyeout

  contains
    function fdebye(z)
      real*8, intent(in) :: z
      real*8 :: fdebye

      fdebye = z*z*z/(exp(z)-1)

    end function fdebye

  end subroutine thermal

  ! Compute Debye-Einstein model vibrational properties.
  subroutine debeins(p,ThetaDinp,T,vol,debye,xabs,en,cv,he,ent,cv_ac,cv_op,ifit)
    use evfunc, only: fv1, fv2
    use varbas, only: phase, vfree, vbracket
    use tools, only: error, gauleg
    use param, only: au2gpa, faterr, half, pckbau, pi, third, twothird, zero
    !--------------------------------------------------------------------------
    ! Calcularemos las propiedades tratando la rama de frecuencias
    ! acusticas por el modelo cuasiarmonico de Debye y la rama de
    ! frecuencias opticas por el modelo de Einstein.
    ! vol : volumen de entrada para el q se calcularan las propiedades
    ! Bstat: es el Bstat.
    ! Kstat= (dBstat/dPstat): Es la Kstat.
    ! Pstat: (dEstat/dV)
    !--------------------------------------------------------------------------

    ! NOTA: Estoy dividiendo el volumen entre Vref que es el minimo de la curva
    ! ajustada estatica=> considerar normalizarlo respecto a otro volumen.
    type(phase), intent(in) :: p
    real*8, intent(in) :: ThetaDinp, T, vol
    real*8, intent(out) :: debye, xabs, en, cv, he, ent, cv_ac, cv_op
    integer, intent(in), optional :: ifit

    real*8, parameter :: eps=1D-12
    integer, parameter :: maxnl=100

    real*8 :: nmolec
    integer :: i, nl
    real*8 :: x(maxnl), w(maxnl), y, z, fdebye
    real*8 :: ThetaD, pstat, bstat
    real*8 :: tmpVol, tmpB, tmpPB, tmpTot
    real*8 :: sum_F_op, summ, debye0
    real*8, allocatable :: freq_i(:), x_i(:), ex_i(:)
    real*8 :: sum_S_op, sum_Cv_op
    real*8 :: en_ac, ent_ac, he_ac
    real*8 :: ent_op, he_op, veq_static, beq_static
    real*8 :: prefac

    real*8, parameter :: explimit = 250d0

    fdebye(z)=z*z*z/(exp(z)-1)

    ! number of molecules per cell and normalization factor (3*vfree-3)/(3*vfree*Z-3)
    nmolec = (p%nfreq + 3)/ vfree / 3
    if (abs(nmolec-p%z) > 1d-6) then
       call error('debeins','(nfreq+3)/3 != vfree*z ; check Z in PHASE',faterr)
    end if
    prefac = real(3*vfree-3,8) / real(3*vfree*p%z-3,8)

    ! Usamos la Theta calculada por Gibbs, pero la normalizamos a las 3 frecuencias
    ! acusticas por lo que hay que dividirla por el factor (n·Z)**(1/3)
    ThetaD = ThetaDinp/vfree**third

    !.error condition controls
    debye=zero
    xabs=zero

    ! presion y b estaticos
    if (present(ifit)) then
       pstat = -fv1(p%pfit%mode(ifit), vol, p%pfit%npar(ifit), p%pfit%apar(:,ifit)) * au2gpa
       bstat = vol * fv2(p%pfit%mode(ifit), vol, p%pfit%npar(ifit), p%pfit%apar(:,ifit)) * au2gpa
       veq_static = p%pfit%veq(ifit)
       beq_static = p%pfit%beq(ifit)
    else
       pstat = -fv1(p%fit_mode, vol, p%npol, p%cpol) * au2gpa
       bstat = vol * fv2(p%fit_mode, vol, p%npol, p%cpol) * au2gpa
       veq_static = p%veq_static
       beq_static = p%beq_static
    end if

    ! allocate frequency arrays
    allocate(freq_i(p%nfreq), x_i(p%nfreq), ex_i(p%nfreq))

    ! frecuencias en gamma a volume vol
    if (size(p%freqg,2) == 1) then
       tmpVol=(vol/p%veq_static)**(1d0/6d0)
       tmpB = (bstat/p%beq_static)**(half)
       tmpPB = (1d0 - twothird*pstat/bstat)**(half)
       tmpTot = tmpVol*tmpB*tmpPB
       freq_i = p%freqg(:,1) * tmpTot
    else
       freq_i = freqg_interpolate(p,vol)
    end if

    if (abs(t) < 1d-5) then
       sum_F_op = sum(freq_i / pckbau / 2d0) * prefac

       en = 9d0*pckbau*ThetaD/8d0 + pckbau*sum_F_op
       cv = zero
       he = 9d0*pckbau*ThetaD/8d0 + pckbau*sum_F_op
       ent = zero
       return
    else
       x_i = freq_i / pckbau / T
       if (any(x_i > (log(huge(x_i(1))-2) ))) then
          ex_i = exp(-x_i)
          sum_Cv_op = sum(x_i**2*ex_i)
       else
          ! calcular cv
          ex_i = exp(x_i)
          sum_Cv_op = sum(x_i**2/(ex_i-1)/(1-1d0/ex_i))
       end if
       cv_op = sum_Cv_op*pckbau*prefac
    end if

    if (thetad.lt.0d0 .or. t.lt.0d0) then
       call error('debeins','ThetaD or T < 0',faterr)
    endif
    y=ThetaD/T
    debye=3d0*pi*pi*pi*pi/y/y/y/15d0

    ! Numerical Debye integral
    if (y <= explimit) then
       !.Loop with increasing number of Legendre points.
       debye0=1d30
       do nl=5,maxnl,5
          call gauleg (0D0,y,x,w,nl)
          summ = 0d0
          do i=1,nl
             summ=summ+w(i)*fdebye(x(i))
          end do
          debye=summ*3d0/y/y/y
          xabs=abs(debye-debye0)
          if (xabs.lt.eps) then
             exit
          else
             debye0=debye
          endif
       end do
       cv_ac = 3d0*pckbau*(4d0*debye - 3d0*y/(exp(y)-1d0))
    else
       cv_ac = 12d0*pckbau*debye
    end if

    !.thermodynamic acoustic vibrational properties
    en_ac = 3d0*pckbau*(ThetaD*3d0/8d0 + T*debye)
    ent_ac= 3d0*pckbau*(debye*4d0/3d0-log(1d0-exp(-y)))
    he_ac = en_ac - T * ent_ac

    ! calculamos el sumatorio para la F,S y Cv optica con las
    ! frecuencias  opticas
    sum_S_op = 0d0
    do i = 1, p%nfreq
       if (x_i(i) <= explimit) then
          sum_S_op = sum_S_op + x_i(i)/(exp(x_i(i))-1)-log(1-exp(-x_i(i)))
       end if
    end do
    sum_S_op = sum_S_op * prefac
    sum_F_op = sum(x_i/2d0 + log(1-exp(-x_i))) * prefac
    he_op = sum_F_op * pckbau * T
    ent_op= sum_S_op * pckbau

    ! sumamos ambas contribuciones para obtener las totales:
    !.thermodynamic vibrational properties
    en  = he_ac  + he_op + T * (ent_ac + ent_op)
    cv  = cv_ac  + cv_op
    ent = ent_ac + ent_op
    he  = he_ac  + he_op

  end subroutine debeins

  ! For phase p, calculate the Helmholtz free energy, entropy, and
  ! CV at temperature T and the volume given by index iv using
  ! the full quasiharmonic approximation with phonon DOS.
  subroutine thermal_qha(p,T,iv,Fvib,S,CV)
    use varbas, only: phase
    use tools, only: quad1
    use param, only: pckbau
    type(phase), intent(in) :: p
    real*8, intent(in) :: T
    integer, intent(in) :: iv
    real*8, intent(out) :: Fvib, S, CV

    integer :: nf
    real*8 :: step, kt
    real*8, allocatable :: emfkt(:), aux(:), tmin(:)
    real*8, allocatable :: l1emfkt(:)
    integer :: i

    real*8, parameter :: logtiny = log(tiny(1d0))

    ! initialize
    kt = pckbau * T

    ! interpolate the phonon DOS at volume v.
    ! f = frequencies (Hartree), d = phDOS, on the same grid.
    nf = size(p%phdos_f,1)
    step = p%phstep

    ! calculate the minimum temperature
    allocate(tmin(nf))
    tmin = -p%phdos_f / (pckbau * logtiny)

    ! calculate emfkt = exp(-omega / (k*T))
    allocate(emfkt(nf),l1emfkt(nf))
    where (T < tmin)
       emfkt = 0d0
    elsewhere
       emfkt = exp(-p%phdos_f / kt)
    end where
    l1emfkt = log(1d0 - emfkt)

    ! note: frequency = 0 is gone from f(:)
    allocate(aux(nf))

    ! calculate fvib
    aux = p%phdos_d(:,1,iv)*(0.5d0 * p%phdos_f + kt * l1emfkt)
    Fvib = quad1(p%phdos_f,aux,step)

    ! calculate entropy (T < tmin implies aux = 0)
    where (T < tmin)
       aux = 0d0
    elsewhere
       aux = p%phdos_d(:,1,iv)*(-pckbau * l1emfkt + p%phdos_f/T * emfkt / (1d0 - emfkt))
    end where
    S = quad1(p%phdos_f,aux,step)

    ! calculate constant-volume heat capacity (T < tmin implies aux = 0)
    where (T < tmin)
       aux = 0d0
    elsewhere
       aux = p%phdos_d(:,1,iv) * (pckbau * (p%phdos_f / kt)**2 * emfkt / (1d0 - emfkt)**2)
    end where
    CV = quad1(p%phdos_f,aux,step)

  end subroutine thermal_qha

  ! Compute Debye extended model vibrational properties (Fvib, S, CV).
  subroutine thermal_debye_extended(p,T,iv,Fvib,S,CV)
    use varbas, only: phase, vfree
    use param, only: pckbau

    type(phase), intent(in) :: p
    real*8, intent(in) :: T
    integer, intent(in) :: iv
    real*8, intent(out) :: Fvib, S, CV

    real*8 :: V, TD, termf, terms, termcv, Fein, Sein, CVein
    real*8 :: UVib, D3, xabs, x, sumc, ex, emx, l1emx
    integer :: i

    real*8, parameter :: vsmall = 1d-80
    real*8, parameter :: lhuge = log(huge(1d0))
    real*8, parameter :: ltiny = log(tiny(1d0))

    ! input properties
    V = p%v(iv)
    TD = p%tde(iv)

    ! calculate the debye contribution
    call thermal(TD,T,D3,xabs,Uvib,CV,Fvib,S)

    ! remove the Debye zero point energy
    Fvib = Fvib - 9d0/8d0 * TD * pckbau * vfree

    ! add the anharmonic contribution
    if (p%tde_nanh > 0) then
       x = T / TD
       termf = 0d0
       terms = 0d0
       termcv = 0d0
       do i = p%tde_nanh, 1, -1
          termf = x * termf - p%tde_anh(i,iv) / real(i+1,8)
          terms = x * terms + p%tde_anh(i,iv)
          termcv = x * termcv + p%tde_anh(i,iv) * real(i,8)
       end do
       termf = x * termf
       terms = x * terms
       termcv = x * termcv
       Fvib = Fvib + vfree * pckbau * T * termf
       S = S + vfree * pckbau * terms
       CV = CV + vfree * pckbau * termcv
    end if

    ! add the einstein contribution
    if (p%tde_nein > 0) then
       Fein = 0d0
       Sein = 0d0
       CVein = 0d0
       sumc = 0d0
       do i = 1, p%tde_nein
          ! skip if the temperature is exactly zero because no contributions
          if (T < vsmall) cycle

          ! calculate TD/T
          x = p%tde_tein(i,iv) / T

          ! F and S in the low temperature limit both go to zero
          if (x < -ltiny) then
             emx = exp(-x)
             l1emx = log(1 - emx)
             Fein = Fein + p%tde_cein(i,iv) * vfree * pckbau * T * l1emx
             Sein = Sein - p%tde_cein(i,iv) * vfree * pckbau * (l1emx - x * emx / (1 - emx))
          end if

          ! CV, in the low temperature limit: CV -> 0
          if (x <= 0.5 * lhuge) then
             ex = exp(x)
             CVein = CVein + p%tde_cein(i,iv) * vfree * pckbau * x**2 * ex / (ex - 1d0)**2
          end if
          sumc = sumc + p%tde_cein(i,iv)
       end do

       Fvib = (1-sumc) * Fvib + Fein
       S = (1-sumc) * S + Sein
       CV = (1-sumc) * CV + CVein
    end if

    ! add the zero-point contribution
    Fvib = Fvib + p%f0(iv)

  end subroutine thermal_debye_extended

  ! Interpolate the frequencies at Gamma to volume v. Returns the
  ! interpolated frequencies in d.
  function freqg_interpolate(p,v) result(ff)
    use varbas, only: phase, vbracket
    use tools, only: leng, error
    use param, only: faterr, uout
    type(phase), intent(in) :: p
    real*8, intent(in) :: v
    real*8 :: ff(p%nfreq)

    integer :: id
    real*8 :: fac

    call vbracket(p,v,id,.true.)
    if (id == 0) then
       write (uout,'("Phase = ",A)') trim(adjustl(p%name(1:leng(p%name))))
       write (uout,'("Volume = ",F17.7)') v
       call error('freqg_interpolate','Requested volume out of grid bounds',faterr)
    else if (id < 0) then
       ff = p%freqg(:,-id)
       return
    end if

    ! linear interpolation of frequencies to this volume
    fac = (v-p%v(id)) / (p%v(id+1) - p%v(id))
    ff = (1d0-fac) * p%freqg(:,id) + fac * p%freqg(:,id+1)

  end function freqg_interpolate

  !> Calculate Debye's integral,
  !
  !                                  |    x^3     |
  ! Debye (y) = 3*y^(-3) * INT (0,y) | ---------- | dx
  !                                  | exp(x) - 1 |
  ! with y = TD/T. The integral is evaluated using a Gauss-Legendre
  ! quadrature.
  function debye_d3(T,TD)
    use tools, only: gauleg
    use param, only: pi
    real*8, intent(in) :: T,TD
    real*8 :: debye_d3

    real*8, parameter :: epslo = 1d-5
    real*8, parameter :: epshi = 250d0
    real*8, parameter :: epsconv = 1d-12
    integer, parameter :: maxnl=100

    real*8 :: y, debye0, summ, xabs
    integer :: nl
    real*8 :: x(maxnl), w(maxnl)

    debye_d3 = 0d0
    if (T < epslo) return

    y = TD/T
    if (y > epslo) then
       debye_d3 = 3d0 * pi**4 / y**3 / 15d0
    else
       debye_d3 = 3d0 * pi**4 / 15d0 * (T/TD)**3
    end if
    if (y > epshi) return

    ! loop with increasing number of Legendre points.
    debye0 = 1d30
    do nl = 5, maxnl, 5
       call gauleg(0d0, y, x, w, nl)
       summ = sum(w(1:nl) * x(1:nl)**3 / (exp(x(1:nl))-1))
       if (y > epslo) then
          debye_d3 = summ * 3d0 / y**3
       else
          debye_d3 = summ * 3d0 * (T/TD)**3
       end if
       xabs = abs(debye_d3 - debye0)
       if (xabs < epsconv) then
          return
       else
          debye0 = debye_d3
       endif
    end do

  end function debye_d3

end module debye
