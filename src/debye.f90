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
  use param
  use tools
  use varbas
  use eos
  implicit none
  private

  ! public
  public :: fill_thetad, get_thetad
  public :: thermal, debeins, thermalphon, thermalomega

  ! if temperature is lower, then gamma(t) = gamma(tlim) 
  real*8, parameter, public :: tlim_gamma = 50d0 + 1d-6
  real*8, parameter, public :: cvlim = 1d-30

contains

  ! Fill Debye-related information: %td(:), %td0 and ntpol/tpol fit.
  ! Optionally, verbose output.
  subroutine fill_thetad(p,verbose)
    use evfunc

    type(phase), intent(inout) :: p
    logical, intent(in) :: verbose

    integer :: j
    character*(mline) :: msg
    integer :: idx, idxx(1)
    integer :: ierr
    character*(mline_fmt) :: fm
    real*8 :: gamma, td, td0, f2s, f3s
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

    ! If the model is debye with poisson, p%td contains the poisson
    ! coefficient; convert them to thetad(V).
    if (p%tmodel == tm_debye_poisson_input) then
       do j = 1, p%nv
          poi = p%td(j)
          fx=2*(1+poi)/3d0/(1-2*poi)
          gx=(1+poi)/3d0/(1-poi)
          hx=2d0*sqrt(fx**3)+sqrt(gx**3)
          pofunc=exp(-log(hx/3)/3)
          f2s = fv2(p%fit_mode,p%v(j),p%npol,p%cpol)
          p%td(j) = (6*pi*pi*vfree*p%v(j)*p%v(j))**third / pckbau * pofunc * sqrt(f2s/mm)
       end do
    end if

    ! fit thetad(V) if debye_input or debye_poisson_input
    if (p%tmodel == tm_debye_input .or. p%tmodel == tm_debye_poisson_input) then
       call fitt_polygibbs(p%tdfit_mode,log(p%v),log(p%td),p%ntpol,p%tpol,ierr,.false.)
       p%ntpol = p%ntpol + 1
       p%tpol(p%ntpol) = 0d0
       if (ierr > 0) call error('get_thetad','Can not fit logTd vs. logV',faterr)
       p%td0 = exp(fv0(p%tdfit_mode,log(p%veq_static),p%ntpol,p%tpol))
    else
       ! fill td0
       f2s = fv2(p%fit_mode,p%veq_static,p%npol,p%cpol)
       p%td0 = (6*pi*pi*vfree*p%veq_static*p%veq_static)**third / pckbau * p%pofunc * sqrt(f2s/mm)
    end if

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
       f2s = fv2(p%fit_mode,p%v(j),p%npol,p%cpol)
       f3s = fv3(p%fit_mode,p%v(j),p%npol,p%cpol)
       call get_thetad(p,p%v(j),f2s,f3s,td,gamma)
       if (p%tmodel /= tm_debye_input .and. p%tmodel /= tm_debye_poisson_input) then
          p%td(j) = td
       end if

       if (verbose) then
          td0 = (6*pi*pi*vfree*p%v(j)*p%v(j))**third / pckbau * p%pofunc * sqrt(f2s/mm)
          write (uout,fm) p%v(j), p%td(j), td0
       end if
    end do

  end subroutine fill_thetad

  subroutine get_thetad(p,v,f2o,f3,td,gamma)
    
    type(phase), intent(in) :: p
    real*8, intent(in) :: v, f2o, f3
    real*8, intent(out) :: td, gamma

    real*8 :: b, f2

    td = 0d0
    gamma = 0d0

    if (f2o < 0d0) then
       call error('fill_thetad','Epp < 0 in get_thetad',warning)
       f2 = 0d0
    else
       f2 = f2o
    end if

    select case(p%tmodel)
    case(tm_debye_input, tm_debye_poisson_input)
       td = exp(fv0(p%tdfit_mode,log(v),p%ntpol,p%tpol))
       gamma = -fv1(p%tdfit_mode,log(v),p%ntpol,p%tpol)
       
    case(tm_debye,tm_debye_einstein)
       td = (6*pi*pi*vfree*v*v)**third / pckbau * p%pofunc * sqrt(f2/mm)
       gamma = -1d0/6d0 - 0.5d0 * (1+v*f3/f2)

    case(tm_debyegrun)
       b = v * f2 * au2gpa
       td = p%td0 * (b / p%beq_static)**p%b_grun / (v / p%veq_static)**p%a_grun
       gamma = p%a_grun - p%b_grun * (1+v*f3/f2)

    case default
       return
    end select

  end subroutine get_thetad

  subroutine thermal (ThetaD,T,debye,xabs,en,cv,he,ent)
    !-------------------------------------------------------------------
    !
    !.thermal - compute Debye model vibrational properties.
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
    real*8 :: x(maxnl), w(maxnl), y, z, sum, debye0, fdebye
    integer :: i, nl
    fdebye(z)=z*z*z/(exp(z)-1)

    !.error condition controls
    !
    debye=zero
    xabs=zero
    if (t<1d-5) then
       en = vfree*9d0*pckbau*ThetaD/8d0 
       cv = zero
       he = vfree*9d0*pckbau*ThetaD/8d0 
       ent = zero
       return
    endif
    if (ThetaD.eq.0d0) then
       en = vfree*3d0*pckbau*T
       cv = vfree*3d0*pckbau
       ent = 1d100
       he = en - T*ent
       return
    endif
    if (thetad.lt.0d0 .or. t.lt.0d0) then
       call error('thermal','Negative ThetaD or T',faterr)
    endif
    y=ThetaD/T
    debye=3d0*pi*pi*pi*pi/y/y/y/15d0
    if (y.gt.250d0) goto 22

    !.Loop with increasing number of Legendre points.
    debye0=1d30
    do nl=5,maxnl,5
       call gauleg (0D0,y,x,w,nl)
       sum=0d0
       do i=1,nl
          sum=sum+w(i)*fdebye(x(i))
       end do
       debye=sum*3d0/y/y/y
       xabs=abs(debye-debye0)
       if (xabs.lt.eps) then
          go to 22
       else
          debye0=debye
       endif
    end do

    !.thermodynamic vibrational properties
22  en  = vfree * 3d0 * pckbau * (ThetaD*3d0/8d0 + T*debye)
    cv  = vfree * 3d0 * pckbau * (4d0*debye - 3d0*y/(exp(y)-1d0))
    ent = vfree * 3d0 * pckbau * (debye*4d0/3d0 - log(1d0-exp(-y)))
    he  = en - T * ent

  end subroutine thermal

  subroutine debeins(p,ThetaDinp,T,vol,debye,xabs,en,cv,he,ent,cv_ac,cv_op,ifit)
    use evfunc
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
    real*8 :: sum_F_op, freq_i(p%nfreq), x_i(p%nfreq), ex_i(p%nfreq), summ, debye0
    real*8 :: sum_S_op, sum_Cv_op
    real*8 :: en_ac, ent_ac, he_ac
    real*8 :: ent_op, he_op, veq_static, beq_static
    real*8 :: prefac

    fdebye(z)=z*z*z/(exp(z)-1)

    ! number of molecules per cell and normalization factor (3*vfree-3)/(3*vfree*Z-3)
    nmolec = (p%nfreq + 3)/ vfree / 3
    if (abs(nmolec-p%z) > 1d-6) then
       call error('debeins','(nfreq+3)/3 != vfree*z ; check Z in PHASE',faterr)
    end if
    prefac = real(3*vfree-3,8) / real(3*vfree*p%z-3,8)

    ! Usamos la Theta calculada por Gibbs, pero la normalizamos a las 3 frecuencias
    ! acusticas por lo que hay que dividirla por el factor (n·Z)**(1/3)
    ! ThetaD = ThetaDinp/((vfree*nmolec)**(third))
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

    ! frecuencias en gamma a volume vol
    tmpVol=(vol/p%veq_static)**(1d0/6d0)
    tmpB = (bstat/p%beq_static)**(half)
    tmpPB = (1d0 - twothird*pstat/bstat)**(half)
    tmpTot = tmpVol*tmpB*tmpPB

    if (size(freq_i) /= size(p%freqg0)) then
       write (uout,'("size(freq_i) = ",I8," size(p%freqg0) = ",I8)') size(freq_i), size(p%freqg0)
       call error('debeins','size(freq_i) /= size(p%freqg0)',faterr)
    end if

    ! calcular x_i a partir de las frecuencias
    freq_i = p%freqg0 * tmpTot 
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
    if (y <= 250d0) then
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
    end if

    !.thermodynamic acoustic vibrational properties
    en_ac = 3d0*pckbau*(ThetaD*3d0/8d0 + T*debye)      
    cv_ac = 3d0*pckbau*(4d0*debye - 3d0*y/(exp(y)-1d0))
    ent_ac= 3d0*pckbau*(debye*4d0/3d0-log(1d0-exp(-y)))
    he_ac = en_ac - T * ent_ac

    ! calculamos el sumatorio para la F,S y Cv optica con las 
    ! frecuencias  opticas 
    sum_F_op = sum(x_i/2d0 + log(1-exp(-x_i))) * prefac
    sum_S_op = sum(x_i/(exp(x_i)-1)-log(1-exp(-x_i))) * prefac
    he_op = sum_F_op * pckbau * T 
    ent_op= sum_S_op * pckbau     

    ! sumamos ambas contribuciones para obtener las totales:
    !.thermodynamic vibrational properties
    en  = he_ac  + he_op + T * (ent_ac + ent_op)
    cv  = cv_ac  + cv_op
    ent = ent_ac + ent_op
    he  = he_ac  + he_op

  end subroutine debeins

  ! For phase p, calculate thermodynamic non-gamma-dependent
  ! thermodynamic properties at V and T using phonon DOS QHA.
  ! Returns: uvib = vibrational contribution to the internal energy,
  ! cv = constant-volume heat capacity, fvib = vibrational
  ! contribution to the free energy, ent = entropy, cv_lowt = 
  ! low-temeprature cv (for the calculation of gamma).
  subroutine thermalphon (p,T,v,uvib,cv,fvib,ent,cv_lowt)
    type(phase), intent(in) :: p
    real*8, intent(in) :: T, v
    real*8, intent(out) :: uvib, cv, fvib, ent, cv_lowt

    integer :: i, nf
    real*8 :: step, kt
    real*8, allocatable :: emfkt(:), aux(:), tmin(:), f(:), d(:)
    real*8, allocatable :: l1emfkt(:)

    real*8, parameter :: hvol = 1d-7
    real*8, parameter :: logtiny = log(tiny(1d0))

    ! xxxx !
    uvib = 0d0
    cv = 0d0
    fvib = 0d0
    ent = 0d0
    cv_lowt = 1d0

    ! initialize
    kt = pckbau * T

    ! interpolate the phonon DOS at volume v. 
    ! f = frequencies (Hartree), d = phDOS, on the same grid.
    nf = size(p%phdos_f,1)
    allocate(f(nf),d(nf))
    step = p%phstep
    f = p%phdos_f
    d = phdos_interpolate(p,v)

    ! calculate the minimum temperature
    allocate(tmin(nf))
    tmin = -f / (pckbau * logtiny)
       
    ! calculate emfkt = exp(-omega / (k*T))
    allocate(emfkt(nf),l1emfkt(nf))
    where (T < tmin)
       emfkt = 0d0
    elsewhere
       emfkt = exp(-f / kt)
    end where
    l1emfkt = log(1d0 - emfkt)

    ! note: frequency = 0 is gone from f(:)
    allocate(aux(nf))

    ! calculate fvib
    do i = 1, nf
       aux = d*(0.5d0 * f + kt * l1emfkt)
    end do
    fvib = quad1(f,aux,step)

    ! calculate entropy (T < tmin implies aux = 0)
    where (T < tmin)
       aux = 0d0
    elsewhere
       aux = d*(-pckbau * l1emfkt + f/T * emfkt / (1d0 - emfkt))
    end where
    ent = quad1(f,aux,step)

    ! calculate constant-volume heat capacity (T < tmin implies aux = 0)
    where (T < tmin)
       aux = 0d0
    elsewhere
       aux = d * (pckbau * (f / kt)**2 * emfkt / (1d0 - emfkt)**2)
    end where
    cv = quad1(f,aux,step)

    ! internal energy
    uvib = fvib + T * ent

    ! Cv is needed in the calculation of gamma -> 0/0 at low temp.
    cv_lowt = cv
    if (T < tlim_gamma) then
       where (tlim_gamma < tmin)
          emfkt = 0d0
       elsewhere
          emfkt = exp(-f / (pckbau * tlim_gamma))
       end where
       where (f > log(huge(hvol))/2*pckbau*tlim_gamma)
          aux = emfkt
       elsewhere
          aux = 1d0 / (1d0 - 1d0/emfkt)
       end where
       cv_lowt = quad1(f,d*(pckbau*(f/pckbau/tlim_gamma)**2 / (emfkt-1) * aux),step)
    end if
    deallocate(f,d,tmin,emfkt,l1emfkt,aux)

  end subroutine thermalphon

  function phdos_interpolate(p,v) result(d)

    type(phase), intent(in) :: p
    real*8, intent(in) :: v
    real*8 :: d(size(p%phdos_f))

    real*8, parameter :: deps = 1d-6

    integer :: i, id
    real*8 :: fac, ppvalu

    call vbracket(p,v,id,.true.)
    if (id == 0) then
       write (uout,'("Phase = ",A)') trim(adjustl(p%name(1:leng(p%name))))
       write (uout,'("Volume = ",F17.7)') v
       call error('phdos_interpolate','Requested volume out of grid bounds',faterr)
    else if (id < 0) then
       d = p%phdos_d(:,1,-id)
       return
    end if

    ! interpolate phonon DOS to this volume
    if (phonsplin) then
       ! not-a-knot cubic spline 
       d = 0d0
       do i = 1, size(p%phdos_f)
          d(i) = ppvalu(p%v,p%phdos_d(i,:,:),p%nv-1,4,v,0)
          if ((p%phdos_d(i,1,id) < deps .or. p%phdos_d(i,1,id+1) < deps) .and. d(i) < deps) then
             d(i) = 0d0
             exit
          end if
       end do
    else
       ! linear 
       fac = (v-p%v(id)) / (p%v(id+1) - p%v(id))
       d = (1d0-fac) * p%phdos_d(:,1,id) + fac * p%phdos_d(:,1,id+1)
    end if

  end function phdos_interpolate
  
  subroutine thermalomega (p,T,v,uvib,cv,he,ent,cv_lowt,gamma)

    type(phase), intent(in) :: p
    real*8, intent(in) :: T, v
    real*8, intent(out) :: uvib, cv, he, ent, cv_lowt
    real*8, intent(out), optional :: gamma

    integer :: nq
    real*8, dimension(size(p%omega,1)) :: omega, gammaj, aux, aux2
    character*(mline) :: msg

    nq = size(p%omega,1)

    call omega_interpolate(p,v,omega,gammaj)

    if (abs(T) < 1d-5) then
       he = 0.5d0 * dot_product(p%wei,omega)
       cv = 0d0
       ent = 0d0
    else
       ! careful with Infs and NaNs
       !   4400 cm^-1 ---> tmin = 17.84 K
       aux = exp(-omega / pckbau / T)
       where (omega > log(huge(cv))/2*pckbau*T)
          aux2 = aux
       elsewhere
          aux2 = 1d0 / (1d0-1d0/aux)
       end where
       he  = dot_product(p%wei,0.5d0*omega+pckbau*T*log(1-aux))
       ent = dot_product(p%wei,-pckbau*log(1-aux)-omega/T*aux2)
       aux2 = pckbau*(omega/pckbau/T)**2 / (aux-1) * aux2
       cv  = dot_product(p%wei,aux2)
    end if
    uvib = he + T * ent

    ! use linear extrapolation for 0 K gamma.
    if (present(gamma)) then
       if (T > tlim_gamma) then
          gamma = dot_product(p%wei,gammaj*aux2) / cv
       else
          write (msg,'(" T = ",F10.2,", gamma extrapolation (low temperature).")') T
          call error('thermalomega',msg,warning)
          aux = exp(-omega / pckbau / tlim_gamma)
          where (omega > log(huge(cv))/2*pckbau*tlim_gamma)
             aux2 = aux
          elsewhere
             aux2 = 1d0 / (1d0-1d0/aux)
          end where
          aux2 = pckbau*(omega/pckbau/tlim_gamma)**2 / (aux-1) * aux2
          gamma = dot_product(p%wei,gammaj*aux2) / dot_product(p%wei,aux2)
       end if
    end if

    ! low temperature cv
    if (T < tlim_gamma) then
       aux = exp(-omega / pckbau / tlim_gamma)
       where (omega > log(huge(cv))/2*pckbau*tlim_gamma)
          aux2 = aux
       elsewhere
          aux2 = 1d0 / (1d0-1d0/aux)
       end where
       aux2 = pckbau*(omega/pckbau/tlim_gamma)**2 / (aux-1) * aux2
       cv_lowt  = dot_product(p%wei,aux2)
    end if

  end subroutine thermalomega

  subroutine omega_interpolate(p,v,omega,gamma)

    type(phase), intent(in) :: p
    real*8, intent(in) :: v
    real*8, dimension(size(p%omega,1)) :: omega, gamma

    integer :: i, j, id

    call vbracket(p,v,id,.true.)
    if (id == 0) then
       write (uout,'("Phase = ",A)') trim(adjustl(p%name(1:leng(p%name))))
       write (uout,'("Volume = ",F17.7)') v
       call error('omega_interpolate','Requested volume out of grid bounds',faterr)
    elseif (id < 0) then
       omega = p%omega(:,-id)
       do j = 1, size(p%omega,1)
          gamma(j) = -p%v(-id)/omega(j)*fv1(omega_fitmode,p%v(-id),p%omega_npol,p%omega_cpol(:,j))
       end do
    end if

    ! polynomial interpolation
    do i = 1, size(p%omega,1)
       omega(i) = fv0(omega_fitmode,v,p%omega_npol,p%omega_cpol(:,i))
       gamma(i) = -v/omega(i)*fv1(omega_fitmode,v,p%omega_npol,p%omega_cpol(:,i))
    end do

  end subroutine omega_interpolate

end module debye
