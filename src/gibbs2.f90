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

!.gibbs2 - (P,T) thermodynamics of crystals.
program gibbs2
  use evfunc, only: nelectrons, evfunc_init, punch_params
  use topcalc, only: fint, iint, interp_input, mxint, nxint, printfreqs, write_energy,&
     dyn_transp, stablevbg, deltag, dyneos, topcalc_init, popinput, eshift_vexp,&
     popenergyfit, plotdh, interpolate, static_transp
  use debye, only: fill_thetad, thermal
  use fit, only: fit_init, fit_pshift
  use eos, only: eosfit_ev_fitt
  use varbas, only: nps, plist, nvs, vlist, nts, tlist, nph, ph, mm, pdefault, &
     phase_max, pstep, tstep, vstep, tdefault, vdefault, temp_pmax, temp_tmax, &
     vfree, temp_vmax, tm_debye, tm_debye_einstein_v, &
     tm_debye_input, tm_debye_einstein, tm_debye_poisson_input, tm_debyegrun, tm_static, &
     phase_init, help_me_and_exit, inplace_sort, varbas_init, process_argv, setup_phases, &
     props_staticeq, phase_punch_pfit, phase_popinfo, quiet
  use tools, only: lower, equal, fopen, fclose, isinteger, ffdate, getword, isreal, fgetline, &
     leng, error, realloc, ioinit, getargs, stdargs
  use param, only: amu2au, faterr, ioerror, ncomms, nwarns, null, pct0, uin, uout, &
     header, mline, marg, title, fileroot, param_init
  implicit none

  interface
     subroutine setvariables (line, lp)
       character*(*), intent(in) :: line
       integer, intent(inout) :: lp
     end subroutine setvariables
  end interface

  integer :: argc
  character*(mline) :: argv(marg), optv, sdate, line, line2, word
  character*(mline) :: fileout
  integer :: ipid, lp, lp2, onps, onts, onvs
  integer :: i, j, iph
  logical :: ok, ok2, ok3, verbose
  logical :: callpf, calleout, dothermal
  logical, allocatable :: dodebye(:)
  integer :: nhouse
  real*8 :: pini, pend, tini, tend, tmaxmin, vini, vend
  real*8 :: vout_ini, vout_end, vout_step, e0
  integer :: ierr, isv, onxint, nid
  real*8 :: adebye, axabs, aen, acv, ahe, aent

  ! initialize
  call param_init()
  call ioinit()
  call getargs(argc,argv)
  call stdargs(argc,argv,optv,uin,uout)

  ! get date
  sdate = ffdate()

  ! initialization of default values for variables
  call evfunc_init()
  call fit_init()
  call varbas_init()
  call topcalc_init()
  fileout = "stdout"//null
  pini = 0d0
  pend = 0d0
  tini = 0d0
  tend = 0d0
  vini = 0d0
  vend = 0d0
  callpf = .false.
  calleout = .false.

  ! Interpret the input options -> changes default options
  call process_argv(argc,argv)

  ! input file root
  fileroot = argv(1)
  fileroot = fileroot(1:index(fileroot,'.')-1)
  if (len(trim(fileroot)) == 0) fileroot = "stdout"
  if (uin.eq.ioerror .or. uout.eq.ioerror) then
     call help_me_and_exit()
  endif

  ! header
  call header()
  if (.not.quiet) then
     write (uout,'("tictac -- ",A)') trim(sdate)
  end if

  ! Start reading
  write (uout,'(/"* Reading main input file")')
  do while (fgetline(uin,line))
     lp=1
     word = getword(word,line,lp)
     word = lower(word)

     ! title title.s
     if (equal(word,'title'//null)) then
        title = line(lp:)

     ! {nat|vfree} nat.i
     elseif (equal(word,'nat'//null).or.equal(word,'vfree'//null)) then
        ok = isinteger(vfree,line,lp)
        if (.not.ok) call error('gibbs2','Error input, VFREE/NAT keyword',faterr)

     ! mm mm.r
     elseif (equal(word,'mm'//null)) then
        ok = isreal(mm,line,lp)
        mm = mm*amu2au
        if (.not.ok) call error('gibbs2','Error input, MM keyword',faterr)

     ! nelectrons nelec.i
     elseif (equal(word,'nelectrons'//null).or.equal(word,'nelec'//null)) then
        ok = isinteger(nelectrons,line,lp)
        if (.not.ok) call error('gibbs2','Error input, NELECTRONS/NELEC keyword',faterr)

     ! pressure 0
     ! pressure npres.i
     ! pressure invpstep.r
     ! pressure pini.r pstep.r pend.r
     ! pressure
     !  p1 p2 p3 p4 ...
     !  p5 p6 ...
     ! endpressure
     elseif (equal(word,'pressure'//null)) then
        pdefault = .false.
        ok = isinteger(nps,line,lp)
        ok2 = isreal(pstep,line,lp)
        ok3 = isreal(pend,line,lp)
        if (ok2 .and. ok3) then
           ! pressure pini.r pstep.r pend.r
           if (ok) then
              pini = real(nps,8)
           else
              pini = pstep
              pstep = pend
              ok2 = isreal(pend,line,lp)
           end if
           if (.not. ok2) call error('gibbs2','Error input, PRESSURE keyword, exp. pini.r pstep.r pend.r',faterr)
           nps = -1

           nps = ceiling((pend - pini) / pstep) + 1
           if (nps <= 0) call error('gibbs2','Error input, PRESSURE keyword: wrong range',faterr)
           if (allocated(plist)) deallocate(plist)
           allocate(plist(nps))
           nps = 0
           do while(pini < pend+1d-12)
              nps = nps + 1
              plist(nps) = pini
              pini = pini + pstep
           end do
        else if (ok2 .and..not.ok3) then
           ! pressure invpstep.r
           ! pstep contains the inv. pstep
           ! allocate when min_ph{pmax} is known
           nps = -1
        else
           if (ok) then
              if (nps == 0) then
                 allocate(plist(1))
                 nps = 1
                 plist(1) = 0d0
              else
                 ! pressure npres.i
                 ! allocate when min_ph{pmax} is known
              end if
           else
              ! pressure
              !  p1 p2 p3 p4 ...
              !  p5 p6 ...
              ! endpressure
              allocate(plist(temp_pmax))
              nps = 1
              do while (.true.)
                 ok = fgetline(uin,line)
                 lp2 = 1
                 if (.not.ok) &
                    call error('gibbs2','Error input, PRESSURE..ENDPRESSURE: missing ENDPRESSURE',faterr)

                 onps = nps
                 do while(isreal(plist(nps),line,lp2))
                    nps = nps + 1
                    if (nps == size(plist)) then
                       call realloc(plist,2*nps)
                    end if
                 end do

                 if (nps == onps) then
                    word = getword(word,line,lp2)
                    word = lower(word)
                    if (equal(word,'endpressure'//null)) then
                       nps = nps - 1
                       call realloc(plist,nps)
                       exit
                    else if (word(1:1) == '#') then
                       cycle
                    else
                       call error('gibbs2','Error input, PRESSURE..ENDPRESSURE: missing ENDPRESSURE',faterr)
                    end if
                 end if
              end do
           end if
        end if

        ! volume input
        ! volume nvol.i
        ! volume invvstep.r
        ! volume vini.r vstep.r vend.r
        ! volume
        !  v1 v2 v3 v4 ...
        !  v5 v6 ...
        ! endvolume
     elseif (equal(word,'volume'//null)) then
        vdefault = .false.
        ok = isinteger(nvs,line,lp)
        ok2 = isreal(vstep,line,lp)
        ok3 = isreal(vend,line,lp)
        if (ok2 .and. ok3) then
           ! volume vini.r vstep.r vend.r
           if (ok) then
              vini = real(nvs,8)
           else
              vini = vstep
              vstep = vend
              ok2 = isreal(vend,line,lp)
           end if
           if (.not. ok2) call error('gibbs2','Error input, VOLUME keyword, expected vini.r vstep.r vend.r',faterr)
           nvs = -1

           nvs = ceiling((vend - vini) / vstep) + 1
           if (nvs <= 0) call error('gibbs2','Error input, VOLUME keyword: wrong range',faterr)
           allocate(vlist(nvs))
           nvs = 0
           do while(vini < vend+1d-12)
              nvs = nvs + 1
              vlist(nvs) = vini
              vini = vini + vstep
           end do
        else if (ok2 .and..not.ok3) then
           ! volume invvstep.r
           ! vstep contains the inv. vstep
           ! set volumes when max and min are known
           nvs = -1
        else
           if (ok) then
              ! volume nvol.i
              ! allocate when min_ph{vmin,vmax} is known
           else
              word = getword(word,line,lp)
              word = lower(word)
              if (equal(word,'input'//null)) then
                 ! volume input
                 nvs = -2
              else
                 ! volume
                 !  v1 v2 v3 v4 ...
                 !  v5 v6 ...
                 ! endvolume
                 allocate(vlist(temp_vmax))
                 nvs = 1
                 do while (.true.)
                    ok = fgetline(uin,line)
                    lp2 = 1
                    if (.not.ok) &
                       call error('gibbs2','Error input, VOLUME..ENDVOLUME: missing ENDVOLUME',faterr)

                    onvs = nvs
                    do while(isreal(vlist(nvs),line,lp2))
                       nvs = nvs + 1
                       if (nvs == size(vlist)) then
                          call realloc(vlist,2*nvs)
                       end if
                    end do

                    if (nvs == onvs) then
                       word = getword(word,line,lp2)
                       word = lower(word)
                       if (equal(word,'endvolume'//null)) then
                          nvs = nvs - 1
                          call realloc(vlist,nvs)
                          exit
                       else if (word(1:1) == '#') then
                          cycle
                       else
                          call error('gibbs2','Error input, VOLUME..ENDVOLUME: missing ENDVOLUME',faterr)
                       end if
                    end if
                 end do
              end if
           end if
        end if

        ! interpolate input
        ! interpolate [def: p]
        !  V
        !  v1 v2
        !  v3 ...
        !  PT
        !  p1 t1 p2 t2
        !  p3 t3 ...
        !  P
        !  p1 p2
        !  p3 ...
        ! endinterpolate
     elseif (equal(word,'interpolate'//null)) then
        word = getword(word,line,lp)
        word = lower(word)
        if (equal(word,'input'//null)) then
           word = getword(word,line,lp)
           word = lower(word)
           if (equal(word,'static'//null)) then
              interp_input = 1
           else
              interp_input = 2
           end if
           cycle
        end if

        allocate(fint(mxint),iint(mxint))
        nxint = 1
        isv = 2
        do while (.true.)
           ok = fgetline(uin,line)
           lp2 = 1
           if (.not.ok) &
              call error('gibbs2','Error input, INTERPOLATE..ENDINTERPOLATE: missing ENDINTERPOLATE',faterr)

           onxint = nxint
           do while(isreal(fint(nxint),line,lp2))
              iint(nxint) = isv
              nxint = nxint + 1
              if (nxint == size(fint)) then
                 call realloc(fint,2*nxint)
                 call realloc(iint,2*nxint)
              end if
           end do

           if (nxint == onxint) then
              word = getword(word,line,lp2)
              word = lower(word)
              if (equal(word,'endinterpolate'//null)) then
                 nxint = nxint - 1
                 call realloc(fint,nxint)
                 call realloc(iint,nxint)
                 exit
              else if (word(1:1) == '#') then
                 cycle
              else if (word(1:1) == 'v') then
                 isv = 1
                 cycle
              else if (word(1:2) == 'pt') then
                 isv = 3
                 cycle
              else if (word(1:1) == 'p') then
                 isv = 2
                 cycle
              else
                 call error('gibbs2','Error input, INTERPOLATE..ENDINTERPOLATE: missing ENDINTERPOLATE',faterr)
              end if
           end if
        end do

     ! temperature 0
     ! temperature ntemp.i
     ! temperature invtstep.r
     ! temperature tini.r tstep.r tend.r
     ! temperature
     !  t1 t2 t3 t4 ...
     !  t5 t6 ...
     ! endtemperature
     elseif (equal(word,'temperature'//null)) then
        tdefault = .false.
        ok = isinteger(nts,line,lp)
        ok2 = isreal(tstep,line,lp)
        ok3 = isreal(tend,line,lp)
        if (ok2 .and. ok3) then
           ! temperature tini.r tstep.r tend.r
           if (ok) then
              tini = real(nts,8)
           else
              tini = tstep
              tstep = tend
              ok2 = isreal(tend,line,lp)
           end if
           if (.not. ok2) call error('gibbs2','Error input, TEMPERATURE keyword, expected tini.r tstep.r tend.r',faterr)
           nts = -1

           nts = ceiling((tend - tini) / tstep) + 1
           if (nts <= 0) call error('gibbs2','Error input, TEMPERATURE keyword: wrong range',faterr)
           allocate(tlist(nts))
           nts = 0
           do while(tini < tend+1d-12)
              nts = nts + 1
              tlist(nts) = tini
              tini = tini + tstep
           end do
        else if (ok2 .and..not.ok3) then
           ! temperature invtstep.r
           ! tstep contains the inv. tstep
           ! allocate when min_ph{pmax} is known
           nts = -1
        else
           if (ok) then
              if (nts == 0) then
                 allocate(tlist(1))
                 nts = 1
                 tlist(1) = 0d0
              else if (nts == -1) then
                 allocate(tlist(1))
                 nts = 1
                 tlist(1) = pct0
              else
                 ! temperature ntemp.i
                 ! allocate when min_ph{pmax} is known
              end if
           else
              ! temperature
              !  t1 t2 t3 t4 ...
              !  t5 t6 ...
              ! endtemperature
              allocate(tlist(temp_tmax))
              nts = 1
              do while (.true.)
                 ok = fgetline(uin,line)
                 lp2 = 1
                 if (.not.ok) &
                    call error('gibbs2','Error input, TEMPERATURE..ENDTEMPERATURE: missing ENDTEMPERATURE',faterr)

                 onts = nts
                 do while(isreal(tlist(nts),line,lp2))
                    nts = nts + 1
                    if (nts == size(tlist)) then
                       call realloc(tlist,2*nts)
                    end if
                 end do

                 if (nts == onts) then
                    word = getword(word,line,lp2)
                    word = lower(word)
                    if (equal(word,'endtemperature'//null)) then
                       nts = nts - 1
                       call realloc(tlist,nts)
                       exit
                    else if (word(1:1) == '#') then
                       cycle
                    else
                       call error('gibbs2','Error input, TEMPERATURE..TEMPERATURE: missing ENDTEMPERATURE',faterr)
                    end if
                 end if
              end do
           end if
        end if

     ! phase name.s [Z z.r] [poisson sigma.r] [file file.s]
     ! tmodel {static|debye_input|debye_static|debye_sc|debye_staticbv|
     !        debye_einstein|debye_poisson_input|einstein} [sigma.r] ...
     !   v1 e1 [td1]
     !   ...
     ! endphase
     elseif (equal(word,'phase'//null)) then
        nph = nph + 1
        if (nph > phase_max) then
           call error('gibbs2','Too many phases, increase phase_max',faterr)
        end if

        line2 = line(lp:)
        call phase_init(ph(nph),line2)

     ! activate {name.s|num.i} {n1 n2 ...|ALL}
     elseif (equal(word,'activate'//null)) then

        iph = 0
        ok = isinteger(iph,line,lp)
        if (.not.ok) then
           lp2 = lp
           word = getword(word,line,lp)
           word = lower(word)
           do i = 1, nph
              if (equal(word,lower(ph(i)%name))) then
                 iph = i
                 exit
              endif
           end do
        end if

        if (iph == 0) then
           lp = lp2
           word = getword(word,line,lp)
           word = lower(word)
           if (equal(word,'all'//null)) then
              do i = 1, nph
                 ph(i)%dyn_active = .true.
              end do
           else
              call error('gibbs2','Unknown phase in ACTIVATE keyword',faterr)
           end if
        else
           nid = 0
           ok = isinteger(nid,line,lp)
           if (ok) then
              do while (ok)
                 if (nid < 1 .or. nid > ph(iph)%nv) &
                    call error('gibbs2','Volume id out of range in ACTIVATE keyword',faterr)
                 ph(iph)%dyn_active(nid) = .true.
                 ok = isinteger(nid,line,lp)
              end do
           end if
        end if

     ! printfreq/printfreqs (topcalc)
     elseif (equal(word,'printfreqs'//null) .or. equal(word,'printfreq'//null)) then
        callpf = .true.

     ! printfreq/printfreqs (topcalc)
     elseif (equal(word,'eoutput'//null)) then
        calleout = .true.
        ok = isreal(vout_ini,line,lp)
        if (ok) then
           ok = ok .and. isreal(vout_step,line,lp)
           ok = ok .and. isreal(vout_end,line,lp)
           if (.not.ok) call error('gibbs2','Error input, EOUTPUT keyword, expected vini.r vstep.r vend.r',faterr)
        else
           vout_ini = -1d0
           vout_end = -1d0
           vout_step = -1d0
        end if

     elseif (equal(word,'root'//null)) then
        word = getword(word,line,lp)
        fileroot = word(1:leng(word))

     ! set ...
     elseif (equal(word,'set'//null)) then
        call setvariables(line,lp)

     ! end
     elseif (equal(word,'end'//null)) then
        exit

     ! unknown -> let setvariables handle it
     else
        lp = 1
        ! setvariables
        write (uout,'("Unknown keyword : ",A)') word(1:leng(word))
        call error('gibbs2','Unknown keyword in line: '//line(1:leng(line)-1),faterr)
     endif
  enddo
  write (uout,*)

  ! check we have at least one phase
  if (nph <= 0) &
     call error('gibbs2','at least one phase is required',faterr)

  ! some initial info about the run
  call popinput(fileout)

  ! check phase data, set pressure range (varbas)
  call setup_phases()

  ! shift to fit experimental volume (topcalc)
  call eshift_vexp()

  !!! BEGIN the STATIC run !!!

  write (uout,'("** STATIC RUN **"/)')

  ! obtain static properties (varbas)
  call props_staticeq()

  ! Phase information
  write (uout,'("* Phase information after initial setup")')
  do i = 1, nph
     call phase_popinfo(ph(i),i)
     write (uout,'("  Fit to static E(V) data:")')
     if (.not.ph(i)%staticmin) then
        write (uout,'("# Static volume (V0) corresponds to a extreme of the grid, not the real eq. volume")')
     end if
     call eosfit_ev_fitt(ph(i))
     call punch_params(uout,ph(i)%fit_mode,ph(i)%npol,ph(i)%cpol)
     call phase_punch_pfit(ph(i))
  end do

  ! print the (new) static energy
  if (calleout) then
     call write_energy(vout_ini,vout_step,vout_end)
  end if

  ! Output fitted static energies (topcalc), enthalpy plot
  call popenergyfit()
  call plotdh()

  ! Calculate or fit debye temperatures at input volumes (verbose if any phase is debye)
  allocate(dodebye(nph))
  do i = 1, nph
     dodebye(i) = (ph(i)%tmodel == tm_debye_input .or. ph(i)%tmodel == tm_debye .or.&
        ph(i)%tmodel == tm_debyegrun .or. ph(i)%tmodel == tm_debye_einstein .or.&
        ph(i)%tmodel == tm_debye_einstein_v .or. ph(i)%tmodel == tm_debye_poisson_input)
  end do
  if (any(dodebye)) write (uout,'("* Computed Debye temperatures from static data")')
  do i = 1, nph
     if (dodebye(i)) &
        call fill_thetad(ph(i),.true.)
  end do
  if (any(dodebye)) write (uout,*)

  ! transition pressures, static (topcalc)
  call static_transp()

  ! Print Debye-Eisntein frequencies (topcalc)
  if (callpf) call printfreqs()

  !!! END of the STATIC run !!!

  ! determine if we want the thermal run
  dothermal = .false.
  do i = 1, nph
     dothermal = dothermal .or. (ph(i)%tmodel /= tm_static)
  end do

  if (dothermal) then
     ! prepare for the THERMAL run
     ! check that mm and vfree were given in input
     if (mm < 0d0) &
        call error('gibbs2','thermal run requested but no molecular mass (MM) not found',faterr)
     if (vfree <= 0) &
        call error('gibbs2','thermal run requested but no VFREE found',faterr)

     ! Set temperature range if not given in input
     write (uout,'("* Temperature range examined")')
     if (any(dodebye)) then
        tmaxmin = 1d30
        do i = 1, nph
           if (dodebye(i)) then
              do j = 1, ph(i)%nv
                 tmaxmin = min(tmaxmin,ph(i)%td(j))
              end do
           end if
        end do
        tmaxmin = tmaxmin * 1.5d0
        tmaxmin = max(tmaxmin,pct0)
        write (uout,'("  Min_{DebyeT} (K): ",F12.3)') tmaxmin / 1.5d0
     else
        tmaxmin = pct0
     end if
     if (tdefault) then
        nts = 100
        allocate(tlist(nts))
        tlist(1) = 0d0
        tstep = tmaxmin / real(nts-1,8)
        do i = 2, nts
           tlist(i) = tlist(i-1) + tstep
        end do
     else if (.not.allocated(tlist)) then
        if (nts > 0) then
           allocate(tlist(nts))
           tlist(1) = 0d0
           tstep = tmaxmin / real(nts-1,8)
           do i = 2, nts
              tlist(i) = tlist(i-1) + tstep
           end do
        else
           nts = floor(tmaxmin / tstep) + 1
           allocate(tlist(nts))
           tlist(1) = 0d0
           do i = 2, nts
              tlist(i) = tlist(i-1) + tstep
           end do
        end if
     end if
     if (allocated(dodebye)) deallocate(dodebye)

     ! sort the temperatures - externalfvib is already sorted
     call inplace_sort(tlist(1:nts))
     write (uout,'("  Temperature range (K): ",F12.3," -> ",F12.3)') &
        tlist(1), tlist(nts)
     write (uout,'("  Number of T points: ",I6)') nts
     write (uout,*)

     ! Properties at input temperatures (topcalc)
     call dyneos()

     ! delta_G(T,p) (topcalc)
     call deltag()

     ! G(T,p), V(T,p) and B(T,p) for the stable phase (topcalc)
     call stablevbg()

     ! transitions pressures (T) (topcalc).
     call dyn_transp()
  end if

  ! user-requested interpolations
  call interpolate()

  ! end
  write (uout,'("GIBBS2 ended successfully (",I3," WARNINGS, ",I3," COMMENTS)"/)')&
     nwarns, ncomms
  if (.not.quiet) then
     sdate = ffdate()
     write (uout,'("tictac -- ",A/)') trim(sdate)
  end if

  ! deallocate
  if (allocated(vlist)) deallocate(vlist)
  if (allocated(plist)) deallocate(plist)
  if (allocated(tlist)) deallocate(tlist)
  if (allocated(fint)) deallocate(fint)
  if (allocated(iint)) deallocate(iint)
  if (allocated(ph)) deallocate(ph)

end program gibbs2
