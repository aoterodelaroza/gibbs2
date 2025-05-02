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

module varbas
  use fit, only: fitinfo, mmpar
  use evfunc, only: fit_strain_x1, fit_strain_v, fit_strain
  use param, only: mline, half, &
     ifmt_p, ifmt_v, ifmt_x, ifmt_g, ifmt_b, ifmt_bp, ifmt_bpp, ifmt_e, ifmt_bm, ifmt_t,&
     ifmt_cv, ifmt_s, ifmt_dpdt, ifmt_alpha, ifmt_interp, ifmt_eprec, ifmt_aic, ifmt_order,&
     ifmt_ef, ifmt_nef
  implicit none
  public

  integer :: vfree   !< number of atoms per primitive cell
  real*8 :: mm       !< molecular mass (atomic units)

  ! phases
  integer, parameter :: phase_max = 50
  integer, parameter :: phase_vmax = 200
  integer, parameter :: phase_freqgmax = 100
  integer, parameter :: phdos_fmax = 40000
  integer, parameter :: phdos_qmax = 10000
  integer, parameter :: edos_emax = 10000
  integer, parameter :: omega_fitorder = 6
  integer, parameter :: omega_fitmode = fit_strain * 10000 + fit_strain_x1 * 100 + omega_fitorder
  integer, parameter :: ftsel_fitmode = fit_strain * 10000 + fit_strain_v * 100 + 4
  type phase
     ! basic info
     character*(mline) :: name = ""
     logical :: pvdata = .false.
     integer :: nv = 0
     real*8 :: pmin, pmax
     real*8 :: veq_static, eeq_static, beq_static
     real*8 :: z
     integer :: tmodel
     real*8, allocatable :: v(:), e(:)
     real*8, allocatable :: gtp(:,:), vtp(:,:), btp(:,:)
     logical, allocatable :: didtp(:,:)

     ! input units
     integer :: units_e, units_v, units_p, units_f

     ! energy fits
     integer :: fit_mode, reg_mode, sfit_mode, cvfit_mode, tdfit_mode
     integer :: npol
     real*8 :: cpol(0:mmpar)
     type(fitinfo) :: pfit
     real*8 :: rms, maxdev, r2, aic, bic
     integer :: nfix = 0
     integer :: idfix(0:mmpar)
     real*8 :: obelix(0:mmpar)
     logical :: staticmin

     ! scaling
     integer :: scaltype
     real*8 :: vscal, bscal
     real*8 :: scale_a1, scale_a2
     integer :: iscal
     real*8 :: eec_t, eec_p

     ! interpolations
     integer :: ninterp = 0
     real*8, allocatable :: interp(:,:)

     ! extrapolations
     logical :: extend = .false.

     ! dyn (all models)
     logical, allocatable :: dyn_active(:)

     ! debye
     real*8 :: poisson = 0.25d0
     real*8 :: pofunc
     real*8 :: a_grun = -1d0/6d0
     real*8 :: b_grun = half
     real*8 :: td0
     real*8, allocatable :: td(:)
     integer :: ntpol
     real*8 :: tpol(0:mmpar)

     ! debye-einstein
     integer :: nfreq = 0
     real*8, allocatable :: freqg(:,:)

     ! qha (using phDOS)
     real*8 :: phstep
     real*8, allocatable :: phdos_f(:), phdos_d(:,:,:)

     ! external fvib
     integer :: nfvib_t
     real*8, allocatable :: fvib_t(:)
     real*8, allocatable :: fvib_f(:,:)
     real*8, allocatable :: fvib_s(:,:)
     real*8, allocatable :: fvib_cv(:,:)

     ! debye extended
     real*8, allocatable :: tde(:)
     real*8, allocatable :: f0(:)
     integer :: tde_nanh
     integer :: tde_nein
     real*8, allocatable :: tde_anh(:,:)
     real*8, allocatable :: tde_cein(:,:)
     real*8, allocatable :: tde_tein(:,:)

     ! electronic contribution to the free energy
     integer :: nelec = 0
     integer :: emodel
     integer :: eunits_e
     logical :: efree = .false.
     real*8, allocatable :: nefermi(:)
     real*8, allocatable :: fel_cpol(:,:), tsel_cpol(:,:)
  end type phase
  integer :: nph
  type(phase), allocatable :: ph(:)

  ! pressure
  integer, parameter :: temp_pmax = 500
  integer :: nps
  real*8, allocatable :: plist(:)    !< pressure list
  real*8, parameter :: warn_pmax = 500d0
  logical :: pdefault
  real*8 :: pstep

  ! temperature
  integer, parameter :: temp_tmax = 500
  integer :: nts
  real*8, allocatable :: tlist(:)    !< temperature list
  logical :: tdefault
  real*8 :: tstep

  ! volume
  integer, parameter :: temp_vmax = 500
  integer :: nvs
  real*8, allocatable :: vlist(:)    !< volume list
  logical :: vdefault
  real*8 :: vstep

  ! frequencies
  real*8, parameter :: fnegcrit = -1d-2 !< negative freq. criterion (cm_1)
  real*8 :: ignore_neg_cutoff = 20d0 !< negative freqs below this trigger deactivation
  real*8 :: fgrid_step = -1d0
  logical :: renormalize = .true. !< renormalize the phDOS in QHA?

  ! temperature models
  integer, parameter :: tm_static = 1
  integer, parameter :: tm_debye_input = 2
  integer, parameter :: tm_debye = 3
  integer, parameter :: tm_debye_einstein = 4
  integer, parameter :: tm_debye_einstein_v = 5
  integer, parameter :: tm_qhafull = 6
  integer, parameter :: tm_debyegrun = 7
  integer, parameter :: tm_debye_poisson_input = 8
  integer, parameter :: tm_debye_extended = 9
  integer, parameter :: tm_externalfvib = 10

  ! scaling modes
  integer, parameter :: scal_noscal = 1
  integer, parameter :: scal_pshift = 2
  integer, parameter :: scal_bpscal = 3
  integer, parameter :: scal_apbaf = 4
  integer, parameter :: scal_use = 5

  ! electronic models
  integer, parameter :: em_no = 1
  integer, parameter :: em_sommerfeld = 2
  integer, parameter :: em_pol4 = 4

  ! p(v) input data fit
  integer :: newpts
  real*8 :: facexpand

  ! dyn calc. properties in output -> coordinated with dyneos_calc (topcalc)
  integer, parameter :: mpropout = 29
  character*14, parameter :: propname(mpropout) = (/&
     "p(GPa)        ",&
     "T(K)          ",&
     "V(bohr^3)     ",&
     "Estatic(Ha)   ",&
     "G(kJ/mol)     ",&
     "Gerr(kJ/mol)  ",&
     "psta(GPa)     ",&
     "pth(GPa)      ",&
     "B(GPa)        ",&
     "U-Esta(kJ/mol)",&
     "Cv(J/molK)    ",&
     "F-Esta(kJ/mol)",&
     "S(J/molK)     ",&
     "ThetaD(K)     ",&
     "gamma         ",&
     "alpha(10^-5/K)",&
     "dp/dT(GPa/K)  ",&
     "Bs(GPa)       ",&
     "Cp(J/molK)    ",&
     "BTp           ",&
     "BTpp(GPa-1)   ",&
     "Fvib(kJ/mol)  ",&
     "Fel(kJ/mol)   ",&
     "Uvib(kJ/mol)  ",&
     "Uel(kJ/mol)   ",&
     "Svib(J/molK)  ",&
     "Sel(J/molK)   ",&
     "Cv,vib(J/molK)",&
     "Cv,el(J/molK) "&
     /)
  integer, parameter :: propfmt(0:mpropout) = (/&
     1,&
     ifmt_p,&
     ifmt_t,&
     ifmt_v,&
     ifmt_e,&
     ifmt_g,&
     ifmt_g,&
     ifmt_p,&
     ifmt_p,&
     ifmt_b,&
     ifmt_e,&
     ifmt_cv,&
     ifmt_e,&
     ifmt_s,&
     ifmt_t,&
     ifmt_x,&
     ifmt_alpha,&
     ifmt_dpdt,&
     ifmt_b,&
     ifmt_cv,&
     ifmt_bp,&
     ifmt_bpp,&
     ifmt_e,&
     ifmt_e,&
     ifmt_e,&
     ifmt_e,&
     ifmt_s,&
     ifmt_s,&
     ifmt_cv,&
     ifmt_cv&
     /)

  ! operation flags
  logical :: dotrans      !< calculate transition pressures?
  logical :: doerrorbar   !< calculate error bars?
  logical :: doefit       !< write input and fitted energy?
  logical :: doplotdh     !< write static dH plot?
  logical :: phonsplin    !< use splines to fit phonon DOS?
  integer :: writelevel   !< 0 - nothing ; 1 - eos only ; 2 - all
  logical :: quiet        !< if true, do not print timestamps

  ! private functions
  private :: phase_sort, phase_checkconvex

  ! private parameters
  integer, parameter :: sfit_order = 10 !< maximum fit order for entropy vs. volume fits
  integer, parameter :: tdfit_order = 10 !< maximum fit order for ThetaD vs. volume fits

contains

  ! Initialize the variables in this module
  subroutine varbas_init()

    mm = -1d0
    vfree = -1
    !
    nps = -1
    nts = -1
    nvs = -1
    nph = 0
    pstep = 0.1d0
    tstep = 1d0
    !
    pdefault = .true.
    tdefault = .true.
    vdefault = .true.
    !
    dotrans = .true.
    doerrorbar = .false.
    doefit = .true.
    doplotdh = .true.
    phonsplin = .false.
    writelevel = 2
    quiet = .false.
    !
    newpts = 20
    facexpand = 0.40d0

    ! allocate phases
    allocate(ph(phase_max))

  end subroutine varbas_init

  ! Process the command-line parameters and set the approrpriate flags.
  subroutine process_argv(argc,argv)
    use param, only: mline, marg, null
    use tools, only: equal
    integer, intent(inout) :: argc
    character*(mline), intent(inout) :: argv(marg)

    integer :: i, nargc

    nargc = 0
    do i = 1, argc
       if (equal(argv(i)//null,"-q"//null) .or. equal(argv(i)//null,"--quiet"//null)) then
          quiet = .true.

       elseif (equal(argv(i)//null,"-n"//null) .or. equal(argv(i)//null,"--noplot"//null)) then
          writelevel = 0

       else if (equal(argv(i)//null,"-e"//null) .or. equal(argv(i)//null,"--eos"//null)) then
          writelevel = 1

       else if (equal(argv(i)//null,"-b"//null) .or. equal(argv(i)//null,"--errorbar"//null)) then
          doerrorbar = .true.

       else if (equal(argv(i)//null,"-t"//null) .or. equal(argv(i)//null,"--notrans"//null)) then
          dotrans = .false.

       else if (equal(argv(i)//null,"-d"//null) .or. equal(argv(i)//null,"--noplotdh"//null)) then
          doplotdh = .false.

       else if (equal(argv(i)//null,"-f"//null) .or. equal(argv(i)//null,"--noefit"//null)) then
          doefit = .false.

       else if (equal(argv(i)//null,"-h"//null) .or. equal(argv(i)//null,"-?"//null) .or. &
          equal(argv(i)//null,"--help"//null)) then
          call help_me_and_exit()
       else
          nargc = nargc + 1
          argv(nargc) = argv(i)
       end if
    end do
    argc = nargc

  end subroutine process_argv

  ! Write help to the standard output and exit.
  subroutine help_me_and_exit()
    use param, only: uout, header

    call header()
    write (uout,'("")')
    write (uout,'("Usage: gibbs2 [-nebtdfh] [input [output]]")')
    write (uout,'("")')
    write (uout,'(" -n  --noplot")')
    write (uout,'("         Inhibits all the auxiliary files and plots written by gibbs2.")')
    write (uout,'("         The only output written goes to stdout and stderr.")')
    write (uout,'("         (SET WRITELEVEL 0) ")')
    write (uout,'("")')
    write (uout,'(" -q  --quiet")')
    write (uout,'("         Do not print timestamps to the output.")')
    write (uout,'("")')
    write (uout,'(" -e --eos")')
    write (uout,'("         Same as -n, but the .eos and .eos_static files are written.")')
    write (uout,'("         (SET WRITELEVEL 1)")')
    write (uout,'("")')
    write (uout,'(" -b --errorbar")')
    write (uout,'("         Calculate and output the error bars for each thermodynamic quantity. ")')
    write (uout,'("         The error values are marked by an ''e'' at the beginning of the line in ")')
    write (uout,'("         the .eos file.")')
    write (uout,'("         (SET ERRORBAR)")')
    write (uout,'(" ")')
    write (uout,'(" -t --notrans")')
    write (uout,'("         Do not compute transition pressures.")')
    write (uout,'("         (SET NOTRANS)")')
    write (uout,'("")')
    write (uout,'(" -d --noplotdh")')
    write (uout,'("         Do not produce plots of enthalpy differences to the first phase in input.")')
    write (uout,'("         (SET NOPLOTDH)")')
    write (uout,'("")')
    write (uout,'(" -f --noefit")')
    write (uout,'("         Do not produce plots of the input and fitted static energy.")')
    write (uout,'("         (SET NOEFIT)")')
    write (uout,'("")')
    write (uout,'(" -h --help -?")')
    write (uout,'("         This help.")')
    write (uout,'("")')
    stop

  end subroutine help_me_and_exit

  ! Returns the number of phases that are not pressure-volume.
  function n_not_pv()
    integer :: n_not_pv

    integer :: i

    n_not_pv = 0
    do i = 1, nph
       if (.not.ph(i)%pvdata) n_not_pv = n_not_pv + 1
    end do

  end function n_not_pv

  ! Returns the number of phases that are E(V) and have a minimum.
  function n_not_pv_min()

    integer :: n_not_pv_min

    integer :: i

    n_not_pv_min = 0
    do i = 1, nph
       if (ph(i)%staticmin.and..not.ph(i)%pvdata) n_not_pv_min = n_not_pv_min + 1
    end do

  end function n_not_pv_min

  !> Read colums i1 and i2 from the phDOS file file. Interpret it
  !> using units units and return the frequencies f(:) and phDOS
  !> d(:). nread = number of read frequencies. fstep = frequency
  !> step. didinterp = whether interpolation had to be used to satisfy
  !> fstep.
  subroutine read_phdos(file,i1,i2,units,f,d,fstep,nread,didinterp)
    use param, only: mline, uout, faterr, ha2cm_1, ioread, thz2cm_1, units_f_ha, units_f_thz
    use tools, only: fopen, leng, isreal, leng, getword, fgetline, error, fclose
    character*(mline), intent(in) :: file
    integer, intent(in) :: i1, i2
    integer, intent(in) :: units
    real*8, intent(out) :: f(:), d(:)
    real*8, intent(inout) :: fstep
    integer, intent(out) :: nread
    logical, intent(out) :: didinterp

    real*8, parameter :: eps_nointerp = 1d-2 ! cm-1

    real*8 :: ff(size(f)), dd(size(d))
    integer :: lu, lp, nn, idum, i, nn2
    character*(mline) :: line, word
    logical :: ok, doneg
    integer :: msize
    real*8 :: fdum
    real*8 :: ffnegcrit, fignore, eps
    integer :: nfmax

    ! cutoffs with units
    ffnegcrit = fnegcrit
    fignore = abs(ignore_neg_cutoff)
    if (units == units_f_ha) then
       ffnegcrit = ffnegcrit / ha2cm_1
       fignore = fignore / ha2cm_1
    elseif (units == units_f_thz) then
       ffnegcrit = ffnegcrit / thz2cm_1
       fignore = fignore / thz2cm_1
    end if

    msize = size(f)
    lu = fopen(lu,file,ioread)
    nn = 0
    dd = 0d0
    ff = 0d0
    doneg = .false.
    do while (fgetline(lu,line))
       lp = 1
       nn = nn + 1
       if (nn > msize) then
          write (uout,'("File: ",A)') file(1:leng(file))
          call error('read_phdos','size of array exceeded (input grid) -> phdos_fmax',faterr)
       end if
       idum = 0
       do while(lp < leng(line))
          idum = idum + 1
          if (idum == i1) then
             ok = isreal(fdum,line,lp)
             if (fdum < ffnegcrit) then
                if (fdum < -fignore) &
                   doneg = .true.
                nn = nn - 1
                exit
             end if
             fdum = max(fdum,0d0)
             if (.not. ok) then
                write (uout,'("File: ",A)') file(1:leng(file))
                call error('read_phdos','could not read freq field',faterr)
             end if
             ff(nn) = fdum
          else if (idum == i2) then
             ok = isreal(dd(nn),line,lp)
             dd(nn) = max(dd(nn),0d0)
             if (.not. ok) then
                write (uout,'("File: ",A)') file(1:leng(file))
                call error('read_phdos','could not read DOS field',faterr)
             end if
          else if (idum > i1 .and. idum > i2) then
             exit
          else
             word = getword(word,line,lp)
          end if
       end do
    end do

    call fclose(lu)

    ! set the frequency step and number of freqs needed in the interpolated f
    if (fstep < 0d0) then
       fstep = ff(nn)-ff(nn-1)
    end if
    nfmax = ceiling(maxval(ff(1:nn)) / fstep) + 1
    if (nfmax > phdos_fmax) then
       call error('read_phdos','size of array exceeded (fine grid) -> phdos_fmax',faterr)
    end if
    do i = 1, nfmax
       f(i) = (i-1) * fstep
    end do

    ! linear interpolation to the fixed grid
    d(1) = 0d0
    nn2 = 1
    didinterp = .false.
    eps = eps_nointerp
    if (units == units_f_ha) then
       eps = eps / ha2cm_1
    elseif (units == units_f_thz) then
       eps = eps / thz2cm_1
    end if
    do i = 1, nfmax-1
       do while(ff(nn2) < f(i) - eps)
          nn2 = nn2 + 1
          if (nn2 > nn) exit
       end do

       if (nn2 > nn .or. nn2 <= 1) then
          d(i) = 0d0
       else
          if (abs(ff(nn2)-f(i)) < eps) then
             d(i) = dd(nn2)
          else
             didinterp = .true.
             d(i) = dd(nn2-1) + (dd(nn2)-dd(nn2-1)) * (f(i)-ff(nn2-1)) / (ff(nn2)-ff(nn2-1))
          end if
       end if
    end do
    d(nfmax) = 0d0
    nread = nfmax
    if (doneg) nread = -nread

  end subroutine read_phdos

  ! Read nfreq frequencies at gamma from file file. Return them in
  ! array f.
  subroutine read_freqg(file,f)
    use param, only: ioread
    use tools, only: fopen, fclose, fgetline, isreal, realloc
    character*(mline), intent(in) :: file
    real*8, allocatable, intent(inout) :: f(:)

    integer :: lu, nn, lp
    logical :: ok
    character*(mline) :: line
    real*8 :: aux

    if (allocated(f)) deallocate(f)
    allocate(f(300))

    lu = fopen(lu,file,ioread)

    nn = 0
    do while (.true.)
       ok = fgetline(lu,line)
       if (.not.ok) exit

       lp = 1
       do while(isreal(aux,line,lp))
          nn = nn + 1
          if (nn > size(f,1)) &
             call realloc(f,2*nn)
          f(nn) = aux
       end do
    end do
    call realloc(f,nn)
    call fclose(lu)

  end subroutine read_freqg

  ! Read the Fvib, Svib, and Cv as a function of volume and
  ! temperature from an external set of files. Used in the
  ! externalfvib temperature model. The file "file" contains two
  ! columns with the temperature in K and another file. The file in
  ! each line, prefixed by "prefix", has the G(V) data at the
  ! corresponding temperature. It must have two columns, volume and
  ! G(V,T), and the volumes must be equal in number and value to the
  ! volumes in the static curve.  Requires the static volumes (nv, v)
  ! and the static energies (e).  Returns the lits of temperatures and
  ! the Fvib(V,T), Svib(V,T), and Cv(V,T).
  subroutine read_externalfvib(file,prefix,nv,v,nfvib_t,fvib_t,fvib_f,fvib_s,fvib_cv)
    use tools, only: fopen, fclose, fgetline, error, isreal, getword,&
       leng, cat, realloc
    use param, only: ioread, faterr, null
    character*(*), intent(in) :: file
    character*(*), intent(in) :: prefix
    integer, intent(in) :: nv
    real*8, intent(in) :: v(nv)
    integer, intent(out) :: nfvib_t
    real*8, allocatable, intent(inout) :: fvib_t(:)
    real*8, allocatable, intent(inout) :: fvib_f(:,:)
    real*8, allocatable, intent(inout) :: fvib_s(:,:)
    real*8, allocatable, intent(inout) :: fvib_cv(:,:)

    integer :: lu, lu2, lp, lp2
    integer :: iv
    character*(mline) :: line, line2, word
    real*8 :: tdum, vdum, rdum
    logical :: ok
    character*(len(prefix)+1) :: prefix_

    real*8, parameter :: veps = 1d-5

    ! whether the prefix is absolute
    if (prefix(leng(prefix):leng(prefix)) /= '/') then
       prefix_ = cat(prefix,'/'//null)
    else
       prefix_ = prefix
    end if

    ! initialize
    nfvib_t = 0
    if (allocated(fvib_t)) deallocate(fvib_t)
    if (allocated(fvib_f)) deallocate(fvib_f)
    if (allocated(fvib_s)) deallocate(fvib_s)
    if (allocated(fvib_cv)) deallocate(fvib_cv)
    allocate(fvib_t(10))
    allocate(fvib_f(nv,10))
    allocate(fvib_s(nv,10))
    allocate(fvib_cv(nv,10))

    ! loop over the lines in the file
    lu = fopen(lu,file,ioread)
    do while (fgetline(lu,line))
       ! get the temperature
       lp = 1
       ok = isreal(tdum,line,lp)
       if (.not.ok) &
          call error('read_externalfvib','error reading temperature from external Fvib file',faterr)

       ! realloc if necessary and record the temperature
       nfvib_t = nfvib_t + 1
       if (nfvib_t > size(fvib_t,1)) then
          call realloc(fvib_t,2*nfvib_t)
          call realloc(fvib_f,nv,2*nfvib_t)
          call realloc(fvib_s,nv,2*nfvib_t)
          call realloc(fvib_cv,nv,2*nfvib_t)
       end if
       fvib_t(nfvib_t) = tdum

       ! build the path
       word = getword(word,line,lp)
       word = cat(trim(prefix_),word)

       ! open and read the file
       iv = 0
       lu2 = fopen(lu2,word,ioread)
       do while (fgetline(lu2,line2))
          ! read the two values
          lp2 = 1
          ok = isreal(vdum,line2,lp2)
          ok = ok.and.isreal(rdum,line2,lp2)
          if (.not.ok) &
             call error('read_externalfvib','error reading external file: '//&
             word(1:leng(word)),faterr)

          ! check for consistency
          iv = iv + 1
          if (iv > nv) &
             call error('read_externalfvib','extra V/G lines in file: '//&
             word(1:leng(word)),faterr)
          if (abs(v(iv)-vdum) > veps) &
             call error('read_externalfvib','volumes do not match static volumes in file: '//&
             word(1:leng(word)),faterr)

          ! Write down the Fvib. Assume same units as E
          fvib_f(iv,nfvib_t) = rdum

          ! try to read the entropy
          ok = isreal(rdum,line2,lp2)
          if (ok) then
             fvib_s(iv,nfvib_t) = rdum

             ! try to read the heat capacity
             ok = isreal(rdum,line2,lp2)
             if (ok) then
                fvib_cv(iv,nfvib_t) = rdum
             else
                fvib_cv(iv,nfvib_t) = 0d0
             end if
          else
             fvib_s(iv,nfvib_t) = 0d0
          end if
       end do
       call fclose(lu2)
    end do
    call fclose(lu)
    if (iv /= nv) &
       call error('read_externalfvib','inconsistent number of volumes in externalfvib',faterr)

    ! final realloc
    call realloc(fvib_t,nfvib_t)
    call realloc(fvib_f,nv,nfvib_t)
    call realloc(fvib_s,nv,nfvib_t)
    call realloc(fvib_cv,nv,nfvib_t)

  end subroutine read_externalfvib

  !> Initialize a phase structure by parsing input line (line_in).
  !> Returns the phase.
  subroutine phase_init(p,line_in)
    use param, only: uout, mline, faterr, ioread, noerr, null, pct0, uin, units_e_ev, units_e_ha, &
       units_e_ry, units_f_cm1, units_f_ha, units_f_thz, units_p_au, units_p_gpa, units_v_ang3, &
       units_v_bohr3
    use tools, only: getword, lower, isreal, isinteger, leng, equal, fopen, fgetline, cat, error,&
       fclose, realloc
    use evfunc, only: fit_antons, fit_ap2, fit_bm2, fit_bm3, fit_bm4, fit_murn, fit_polygibbs, &
       fit_pt2, fit_pt3, fit_pt4, fit_pt5, fit_strain_bm, fit_strain_inf, fit_strain_lagr,&
       fit_strain_pt, fit_strain_x3, fit_strain_xinv3, fit_vinet, reg_lad, reg_lsq
    type(phase), intent(out) :: p
    character*(mline), intent(in) :: line_in

    integer, parameter :: minterp = 20

    integer :: lp, lp2
    integer :: i, j
    integer :: icol_v, icol_e, icol_td, icol_nef, icol_ph, icol_pol4(8)
    integer :: icol_int(minterp), icol_f0, icol_tde
    integer, allocatable :: icol_anh(:), icol_cein(:), icol_tein(:)
    integer :: idum, ifound, iphdos_1, iphdos_2, isep, isep2, nn, nq, numax, uuin
    character*(mline) :: line, word, prefix, file, linedum, msg, extfvibfile
    real*8 :: fx, gx, hx, zz, eshift
    real*8, allocatable :: ffreq(:)
    logical :: ok, d0, didinterp, havev, havee, haveph

    ! Set defaults for this phase. also, check type definition
    ! basic info
    p%tmodel = tm_debye
    p%emodel = em_no
    p%z = 1d0
    ! input units
    p%units_v = units_v_bohr3
    p%units_e = units_e_ha
    p%units_p = units_p_gpa
    p%units_f = units_f_ha
    p%eunits_e = units_e_ha
    ! energy fits
    p%fit_mode = fit_strain * 10000 + fit_strain_bm * 100 + 0
    p%sfit_mode = p%fit_mode
    p%cvfit_mode = p%fit_mode
    p%tdfit_mode = p%fit_mode
    p%reg_mode = reg_lsq
    ! scaling
    p%scaltype = scal_noscal
    p%vscal = 0d0
    p%eec_p = 0d0
    p%eec_t = pct0
    ! frequency step
    p%phstep = -1d0

    ! Start parsing the input line
    line = line_in
    lp = 1

    ! read phase identifier
    p%name = getword(p%name,line,lp)

    ! input file field management
    icol_v = 1
    icol_e = 2
    icol_td = -1
    icol_nef = -1
    icol_ph = -1
    icol_pol4 = -1
    icol_f0 = -1
    icol_tde = -1

    ! local initializations
    file = ""
    zz = 1d0
    eshift = 0d0
    prefix = "./" // null
    extfvibfile = ""
    allocate(icol_anh(1),icol_cein(1),icol_tein(1))
    icol_anh = -1
    icol_cein = -1
    icol_tein = -1
    p%tde_nanh = 0
    p%tde_nein = 0
    p%ninterp = 0

    ! start reading
    ok = .true.
    do while(.true.)
       word = getword(word,line,lp)
       word = lower(word)
       if (equal(word,'z'//null)) then
          ok = isreal(zz,line,lp)
          p%z = zz
          if (.not.ok) call error('phase_init','Wrong Z in PHASE',faterr)
       elseif (equal(word,'poisson'//null)) then
          ok = isreal(p%poisson,line,lp)
          if (.not.ok) call error('phase_init','Wrong POISSON in PHASE',faterr)
       elseif (equal(word,'fstep'//null)) then
          ok = isreal(p%phstep,line,lp)
          if (.not.ok) call error('phase_init','Wrong FSTEP in PHASE',faterr)
       elseif (equal(word,'file'//null)) then
          file = getword(file,line,lp)
       elseif (equal(word,'prefix'//null)) then
          prefix = getword(prefix,line,lp)
       elseif (equal(word,'fit'//null)) then
          word = getword(word,line,lp)
          word = lower(word)
          if (equal(word,'polygibbs'//null)) then
             p%fit_mode = fit_polygibbs
          elseif (equal(word,'bm2'//null)) then
             p%fit_mode = fit_bm2
          elseif (equal(word,'bm3'//null)) then
             p%fit_mode = fit_bm3
          elseif (equal(word,'bm4'//null)) then
             p%fit_mode = fit_bm4
          elseif (equal(word,'pt2'//null)) then
             p%fit_mode = fit_pt2
          elseif (equal(word,'pt3'//null)) then
             p%fit_mode = fit_pt3
          elseif (equal(word,'pt4'//null)) then
             p%fit_mode = fit_pt4
          elseif (equal(word,'pt5'//null)) then
             p%fit_mode = fit_pt5
          elseif (equal(word,'murn'//null)) then
             p%fit_mode = fit_murn
          elseif (equal(word,'antons'//null)) then
             p%fit_mode = fit_antons
          elseif (equal(word,'vinet'//null)) then
             p%fit_mode = fit_vinet
          elseif (equal(word,'ap2'//null)) then
             p%fit_mode = fit_ap2
          elseif (equal(word,'strain'//null)) then
             word = getword(word,line,lp)
             word = lower(word)
             p%fit_mode = fit_strain * 10000
             if (equal(word,'eulerian'//null).or.equal(word,'bm'//null)) then
                p%fit_mode = p%fit_mode + fit_strain_bm * 100
             elseif (equal(word,'natural'//null).or.equal(word,'pt'//null)) then
                p%fit_mode = p%fit_mode + fit_strain_pt * 100
             elseif (equal(word,'lagrangian'//null).or.equal(word,'lagr'//null)) then
                p%fit_mode = p%fit_mode + fit_strain_lagr * 100
             elseif (equal(word,'infinitesimal'//null).or.equal(word,'inf'//null)) then
                p%fit_mode = p%fit_mode + fit_strain_inf * 100
             elseif (equal(word,'quotient'//null).or.equal(word,'x1'//null)) then
                p%fit_mode = p%fit_mode + fit_strain_x1 * 100
             elseif (equal(word,'x3'//null)) then
                p%fit_mode = p%fit_mode + fit_strain_x3 * 100
             elseif (equal(word,'xinv3'//null) .or. equal(word,'x3inv'//null)) then
                p%fit_mode = p%fit_mode + fit_strain_xinv3 * 100
             elseif (equal(word,'v'//null)) then
                p%fit_mode = p%fit_mode + fit_strain_v * 100
             else
                call error('phase_init','unknown STRAIN in PHASE/FIT keyword',faterr)
             end if
             ok = isinteger(idum,line,lp)
             if (.not.ok) call error('phase_init','Wrong syntax in STRAIN',faterr)
             p%fit_mode = p%fit_mode + idum
          else
             call error('phase_init','unknown FIT in PHASE keyword',faterr)
          end if
       elseif (equal(word,'tmodel'//null)) then
          word = getword(word,line,lp)
          word = lower(word)
          if (equal(word,'static'//null)) then
             p%tmodel = tm_static
          elseif (equal(word,'debye_input'//null)) then
             p%tmodel = tm_debye_input
             icol_td = 0
          elseif (equal(word,'debye'//null)) then
             p%tmodel = tm_debye
          elseif (equal(word,'debye_einstein'//null)) then
             lp2 = lp
             word = getword(word,line,lp)
             word = lower(word)
             if (equal(word,'freqg0'//null)) then
                p%tmodel = tm_debye_einstein
                word = getword(word,line,lp)
                call read_freqg(word,ffreq)
                if (allocated(p%freqg)) deallocate(p%freqg)
                allocate(p%freqg(size(ffreq,1),1))
                p%freqg(:,1) = ffreq
                deallocate(ffreq)
             else
                p%tmodel = tm_debye_einstein_v
                icol_ph = 0
                lp = lp2
             end if
          elseif (equal(word,'debye_poisson_input'//null)) then
             p%tmodel = tm_debye_poisson_input
             icol_td = 0
          elseif (equal(word,'debye_gruneisen'//null)) then
             p%tmodel = tm_debyegrun
             lp2 = lp
             word = getword(word,line,lp)
             word = lower(word)
             if (equal(word,'slater'//null)) then
                p%a_grun = -1d0/6d0
                p%b_grun = half
             elseif (equal(word,'dm'//null)) then
                p%a_grun = -half
                p%b_grun = half
             elseif (equal(word,'vz'//null)) then
                p%a_grun = -5d0/6d0
                p%b_grun = half
             elseif (equal(word,'mfv'//null)) then
                p%a_grun = -0.95d0
                p%b_grun = half
             else
                lp = lp2
                ok = isreal(p%a_grun,line,lp)
                ok = ok .and. isreal(p%b_grun,line,lp)
                if (.not.ok) call error('phase_init','wrong DEBYE_GRUNEISEN input',faterr)
             end if
          elseif (equal(word,'qha'//null) .or. equal(word,'qhafull'//null)) then
             p%tmodel = tm_qhafull
             icol_ph = 0
             iphdos_1 = 1
             iphdos_2 = 2
             do while (.true.)
                lp2 = lp
                word = getword(word,line,lp)
                word = lower(word)
                if (equal(word,'phfield'//null)) then
                   ok = isinteger(icol_ph,line,lp)
                   if (.not.ok) call error('phase_init','wrong PHFIELD number',faterr)
                else if (equal(word,'dosfield'//null)) then
                   ok = isinteger(iphdos_1,line,lp)
                   ok = ok .and. isinteger(iphdos_2,line,lp)
                   if (.not.ok) call error('phase_init','wrong DOSFIELD numbers',faterr)
                else
                   lp = lp2
                   exit
                end if
             end do
          elseif (equal(word,'externalfvib'//null)) then
             p%tmodel = tm_externalfvib
             extfvibfile = getword(extfvibfile,line,lp)
          elseif (equal(word,'debye_extended'//null)) then
             p%tmodel = tm_debye_extended
             ok = isinteger(p%tde_nanh,line,lp)
             ok = ok .and. isinteger(p%tde_nein,line,lp)
             if (.not.ok) call error('phase_init','wrong DEBYE_EXTENDED input',faterr)

             icol_f0 = 0
             icol_tde = 0
             if (allocated(icol_anh)) deallocate(icol_anh)
             if (allocated(icol_cein)) deallocate(icol_cein)
             if (allocated(icol_tein)) deallocate(icol_tein)
             allocate(icol_anh(p%tde_nanh),icol_cein(p%tde_nein),icol_tein(p%tde_nein))
             icol_anh = 0
             icol_cein = 0
             icol_tein = 0
          else
             call error('phase_init','unknown TMODEL in PHASE keyword',faterr)
          end if

       elseif (equal(word,'reg'//null)) then
          word = getword(word,line,lp)
          word = lower(word)
          if (equal(word,'lsq'//null)) then
             p%reg_mode = reg_lsq
          elseif (equal(word,'lad'//null)) then
             p%reg_mode = reg_lad
          else
             call error('phase_init','Wrong REG in PHASE',faterr)
          end if

       elseif (equal(word,'u'//null) .or. equal(word,'using'//null)) then
          isep = lp
          do isep = lp, leng(line)
             if (line(isep:isep) == ":") exit
          end do
          ok = .not.(isep >= leng(line))
          if (.not.ok) call error('phase_init','Wrong USING in PHASE',faterr)

          do isep2 = isep+1, leng(line)
             if (line(isep2:isep2) == ":") exit
          end do
          if (isep2 >= leng(line)) isep2 = -1

          idum = 1
          linedum = adjustl(line(lp:isep-1)//null)
          ok = isinteger(icol_v,linedum,idum)
          if (isep2 > 0) then
             idum = 1
             linedum = adjustl(line(isep+1:isep2-1)//null)
             ok = ok .and. isinteger(icol_e,linedum,idum)
             lp = isep2+1
             ok = ok .and. isinteger(icol_td,line,lp)
             if (.not.ok) call error('phase_init','Wrong USING in PHASE',faterr)
          else
             lp = isep+1
             ok = ok .and. isinteger(icol_e,line,lp)
             if (.not.ok) call error('phase_init','Wrong USING in PHASE',faterr)
          end if

       elseif (equal(word,'interpolate'//null)) then
          p%ninterp = 0
          do while(isinteger(icol_int(p%ninterp+1),line,lp))
             p%ninterp = p%ninterp + 1
             if (p%ninterp >= minterp) &
                call error('phase_init','too many interpolations (increase minterp)',faterr)
          end do
       elseif (equal(word,'fix'//null)) then
          do while(isinteger(idum,line,lp))
             p%nfix = p%nfix + 1
             p%idfix(p%nfix) = idum
             ok = isreal(p%obelix(p%nfix),line,lp)
             if (.not.ok) call error('phase_init','wrong FIX in PHASE keyword',faterr)
          end do
       elseif (equal(word,'eec'//null)) then
          word = getword(word,line,lp)
          word = lower(word)
          if (equal(word,'noscal'//null)) then
             p%scaltype = scal_noscal
             p%iscal = 0
          elseif (equal(word,'pshift'//null)) then
             p%scaltype = scal_pshift
             ok = isreal(p%vscal,line,lp)
             p%iscal = 0
             if (.not.ok) call error('phase_init','wrong PSHIFT in PHASE keyword',faterr)
          elseif (equal(word,'bpscal'//null)) then
             p%scaltype = scal_bpscal
             ok = isreal(p%vscal,line,lp)
             ok = ok .and. isreal(p%bscal,line,lp)
             p%iscal = 0
             if (.not.ok) call error('phase_init','wrong BPSCAL in PHASE keyword',faterr)
          elseif (equal(word,'apbaf'//null)) then
             p%scaltype = scal_apbaf
             ok = isreal(p%vscal,line,lp)
             p%iscal = 0
             if (.not.ok) call error('phase_init','wrong APBAF in PHASE keyword',faterr)
          elseif (equal(word,'use'//null)) then
             p%scaltype = scal_use
             ok = isinteger(p%iscal,line,lp)
             if (.not.ok) call error('phase_init','wrong USE in PHASE keyword',faterr)
          else
             call error('phase_init','wrong EEC in PHASE keyword',faterr)
          end if
       elseif (equal(word,'eec_t'//null)) then
          ok = isreal(p%eec_t,line,lp)
          if (.not.ok) call error('phase_init','wrong EEC_T in PHASE keyword',faterr)
       elseif (equal(word,'eec_p'//null)) then
          ok = isreal(p%eec_p,line,lp)
          if (.not.ok) call error('phase_init','wrong EEC_P in PHASE keyword',faterr)
       elseif (equal(word,'eshift'//null)) then
          ok = isreal(eshift,line,lp)
          if (.not.ok) call error('phase_init','wrong ESHIFT in PHASE keyword',faterr)
       elseif (equal(word,'pvdata'//null)) then
          p%pvdata = .true.
       elseif (equal(word,'units'//null)) then
          do while(.true.)
             lp2 = lp
             word = getword(word,line,lp)
             word = lower(word)
             if (equal(word,'volume'//null)) then
                word = getword(word,line,lp)
                word = lower(word)
                if (equal(word,'bohr3'//null).or.equal(word,'bohr^3'//null).or.&
                   equal(word,'bohr'//null)) then
                   p%units_v = units_v_bohr3
                elseif (equal(word,'ang3'//null).or.equal(word,'ang^3'//null).or.&
                   equal(word,'ang'//null)) then
                   p%units_v = units_v_ang3
                else
                   call error('phase_init','unknown VOLUME unit in PHASE',faterr)
                end if
             else if (equal(word,'pressure'//null)) then
                word = getword(word,line,lp)
                word = lower(word)
                if (equal(word,'au'//null).or.equal(word,'a.u.'//null)) then
                   p%units_p = units_p_au
                elseif (equal(word,'gpa'//null)) then
                   p%units_p = units_p_gpa
                else
                   call error('phase_init','unknown PRESSURE unit in PHASE',faterr)
                end if
             else if (equal(word,'energy'//null)) then
                word = getword(word,line,lp)
                word = lower(word)
                if (equal(word,'hy'//null).or.equal(word,'ha'//null).or.&
                   equal(word,'hartree'//null)) then
                   p%units_e = units_e_ha
                elseif (equal(word,'ry'//null).or.equal(word,'rydberg'//null)) then
                   p%units_e = units_e_ry
                elseif (equal(word,'ev'//null).or.equal(word,'evolt'//null).or.&
                   equal(word,'electronvolt'//null)) then
                   p%units_e = units_e_ev
                else
                   call error('phase_init','unknown ENERGY unit in PHASE',faterr)
                end if
             else if (equal(word,'freq'//null).or.equal(word,'frequency'//null)) then
                word = getword(word,line,lp)
                word = lower(word)
                if (equal(word,'hy'//null).or.equal(word,'ha'//null).or.&
                   equal(word,'hartree'//null)) then
                   p%units_f = units_f_ha
                elseif (equal(word,'cm-1'//null).or.equal(word,'cm^-1'//null).or.&
                   equal(word,'cm_1'//null)) then
                   p%units_f = units_f_cm1
                elseif (equal(word,'thz'//null)) then
                   p%units_f = units_f_thz
                else
                   call error('phase_init','unknown FREQUENCY unit in PHASE',faterr)
                end if
             else if (equal(word,'edos'//null)) then
                word = getword(word,line,lp)
                word = lower(word)
                if (equal(word,'hy'//null).or.equal(word,'ha'//null).or.&
                   equal(word,'hartree'//null)) then
                   p%eunits_e = units_e_ha
                elseif (equal(word,'ry'//null).or.equal(word,'rydberg'//null)) then
                   p%eunits_e = units_e_ry
                elseif (equal(word,'ev'//null).or.equal(word,'evolt'//null).or.&
                   equal(word,'electronvolt'//null)) then
                   p%eunits_e = units_e_ev
                else
                   call error('phase_init','unknown ENERGY unit in PHASE',faterr)
                end if
             else
                lp = lp2
                exit
             end if
          end do
       elseif (equal(word,'extend'//null)) then
          p%extend = .true.
       elseif (equal(word,'nelec'//null)) then
          ok = isinteger(p%nelec,line,lp)
          if (.not.ok) call error('phase_init','wrong NELEC in PHASE keyword',faterr)
       elseif (equal(word,'elec'//null)) then
          word = getword(word,line,lp)
          word = lower(word)
          if (equal(word,'sommerfeld'//null)) then
             ! Sommerfeld model
             p%emodel = em_sommerfeld
             lp2 = lp
             word = getword(word,line,lp)
             word = lower(word)
             if (equal(word,'free'//null)) then
                ! Use free electrons N(ef) value
                p%efree = .true.
             else
                ! Read N(ef) from input file
                p%efree = .false.
                lp = lp2
                ok = isinteger(icol_nef,line,lp)
                if (.not.ok) icol_nef = 0
             end if
          else if (equal(word,'pol4'//null)) then
             ! Smear electrons over calculated density of states
             p%emodel = em_pol4
             ok = isinteger(icol_pol4(1),line,lp)
             if (.not.ok) then
                icol_pol4 = 0
             else
                icol_pol4(2:8) = (/(icol_pol4(1)+j,j=1,7)/)
             end if
          else
             call error('phase_init','unknown ELEC in PHASE',faterr)
          end if
       else
          if (leng(word) == 0) then
             exit
          else
             write (uout,'("Unknown keyword : ",A)') word(1:leng(word))
             write (uout,'("Rest of line : ",A)') line(lp:leng(line))
             call error('phase_init','Error input, PHASE keyword',faterr)
          end if
       end if
    end do

    ! calculate the f(poisson)
    fx=2*(1+p%poisson)/3d0/(1-2*p%poisson)
    gx=(1+p%poisson)/3d0/(1-p%poisson)
    hx=2d0*sqrt(fx**3)+sqrt(gx**3)
    p%pofunc=exp(-log(hx/3)/3)

    ! set column identifiers
    i = 0
    do while(.true.)
       i = i + 1
       if (icol_v == i .or. icol_e == i .or. icol_td == i .or. icol_nef == i .or. &
           icol_ph == i .or. any(icol_int(1:p%ninterp) == i) .or. &
           any(icol_pol4 == i)) cycle
       if (icol_v == 0) then
          icol_v = i
          cycle
       else if (icol_e == 0) then
          icol_e = i
          cycle
       else if (icol_td == 0) then
          icol_td = i
          cycle
       else if (icol_nef == 0) then
          icol_nef = i
          cycle
       else if (any(icol_pol4 == 0)) then
          do j = 1, 8
             if (icol_pol4(j) == 0) then
                icol_pol4(j) = i
                exit
             end if
          end do
          cycle
       else if (icol_ph == 0) then
          icol_ph = i
          cycle
       else if (icol_f0 == 0) then
          !! debye_extended
          icol_f0 = i
          cycle
       else if (icol_tde == 0) then
          icol_tde = i
          cycle
       else if (any(icol_anh(1:p%tde_nanh) == 0)) then
          do j = 1, p%tde_nanh
             if (icol_anh(j) == 0) then
                icol_anh(j) = i
                exit
             end if
          end do
          cycle
       else if (any(icol_cein(1:p%tde_nein) == 0)) then
          do j = 1, p%tde_nein
             if (icol_cein(j) == 0) then
                icol_cein(j) = i
                exit
             end if
          end do
          cycle
       else if (any(icol_tein(1:p%tde_nein) == 0)) then
          do j = 1, p%tde_nein
             if (icol_tein(j) == 0) then
                icol_tein(j) = i
                exit
             end if
          end do
          cycle
       end if
       exit
    end do

    ! open input file
    if (file /= "") then
       uuin = fopen(uuin,file,ioread)
    else
       uuin = uin
    end if

    ! volume-related memory allocation
    allocate(p%v(phase_vmax),p%e(phase_vmax))
    allocate(p%dyn_active(phase_vmax))
    p%dyn_active = .true.

    if (p%ninterp > 0) allocate(p%interp(phase_vmax,p%ninterp))
    if (icol_td > 0) allocate(p%td(phase_vmax))
    if (icol_ph > 0) then
       if (p%tmodel == tm_qhafull) then
          allocate(p%phdos_f(phdos_fmax))
          allocate(p%phdos_d(phdos_fmax,4,phase_vmax))
          p%phdos_f = -1d0
          p%phdos_d = 0d0
       elseif (p%tmodel == tm_debye_einstein_v) then
          if (allocated(p%freqg)) deallocate(p%freqg)
       end if
       numax = -1
    end if
    if (icol_nef > 0) then
       allocate(p%nefermi(phase_vmax))
    end if
    if (any(icol_pol4 > 0)) then
       allocate(p%fel_cpol(4,phase_vmax))
       allocate(p%tsel_cpol(4,phase_vmax))
       p%fel_cpol = 0d0
       p%tsel_cpol = 0d0
    end if
    !! debye_extended
    if (icol_tde > 0) then
       allocate(p%f0(phase_vmax))
       allocate(p%tde(phase_vmax))
       allocate(p%tde_anh(p%tde_nanh,phase_vmax))
       allocate(p%tde_cein(p%tde_nein,phase_vmax))
       allocate(p%tde_tein(p%tde_nein,phase_vmax))
       p%tde_anh = 0d0
       p%tde_cein = 0d0
       p%tde_tein = 0d0
    end if

    ! call phase_init(ph(nph),line2)
    write (uout,'("+ Initializing phase: ",A)') trim(adjustl(p%name(1:leng(p%name))))
    if (len_trim(file) > 0) then
       write (uout,'("  Reading data from file: ",A)') trim(adjustl(file(1:leng(file))))
    else
       write (uout,'("  Reading data from main input")')
    end if
    write (uout,'("  Interpretation of columns in the data file:")')
    write (uout,'("    Volume: ",I3)') icol_v
    write (uout,'("    Energy: ",I3)') icol_e
    if (icol_td > 0) &
       write (uout,'("    Debye temperature: ",I3)') icol_td
    if (icol_nef > 0) &
       write (uout,'("    N(Ef): ",I3)') icol_nef
    if (icol_ph > 0) &
       write (uout,'("    Phonon DOS/frequency file: ",I3)') icol_ph
    if (any(icol_int(1:p%ninterp) > 0)) &
       write (uout,'("    Interpolation fields: ",999(I3,X))') (icol_int(i),i=1,p%ninterp)
    if (icol_pol4(1) > 0) &
       write (uout,'("    Electronic contribution polynomial: ",8(I3,X))') (icol_pol4(i),i=1,8)
    if (icol_f0 > 0) &
       write (uout,'("    Zero-point free energy: ",I3)') icol_f0
    if (icol_tde > 0) &
       write (uout,'("    Debye temperature (extended): ",I3)') icol_tde
    if (any(icol_anh(1:p%tde_nanh) > 0)) &
       write (uout,'("    Anharmonicity polynomial coefs.: ",99(I3,X))') (icol_anh(i),i=1,p%tde_nanh)
    if (any(icol_cein(1:p%tde_nein) > 0)) &
       write (uout,'("    Einstein polynomial coefs.: ",99(I3,X))') (icol_cein(i),i=1,p%tde_nein)
    if (any(icol_tein(1:p%tde_nein) > 0)) &
       write (uout,'("    Einstein polynomial temps.: ",99(I3,X))') (icol_tein(i),i=1,p%tde_nein)

    ! run over the input file
    nn = 0
    didinterp = .false.
    do while (.true.)
       ok = fgetline(uuin,line)
       if (.not.ok) exit

       ! skip comments and blank lines and read endphase
       lp = 1
       word = getword(word,line,lp)
       word = lower(word)
       if (equal(word,'endphase'//null)) then
          exit
       else if (word(1:1) == '#' .or. lp >= leng(line)) then
          cycle
       end if

       ! read the line
       nn = nn + 1
       if (nn > size(p%v)) call phase_realloc_volume(p,2*nn)

       ! parse fields
       lp = 1
       idum = 0
       havev = .false.
       havee = .false.
       haveph = .false.
       ok = .true.
       do while(ok .and. lp < leng(line))
          idum = idum + 1

          if (idum == icol_v) then
             ! volume field
             ok = isreal(p%v(nn),line,lp)
             havev = .true.
          else if (idum == icol_e) then
             ! energy field
             ok = isreal(p%e(nn),line,lp)
             havee = .true.
          else if (idum == icol_td) then
             ! thetad field
             ok = isreal(p%td(nn),line,lp)
          else if (idum == icol_ph) then
             ! phonon DOS or frequencies file name
             ! get name
             word = getword(word,line,lp)
             if (prefix(leng(prefix):leng(prefix)) == '/') then
                word = cat(prefix,word)
             else
                word = cat(cat(prefix,'/'//null),word)
             end if

             ! read frequency information
             nq = numax
             if (p%tmodel == tm_qhafull) then
                call read_phdos(word,iphdos_1,iphdos_2,p%units_f,&
                   p%phdos_f(:),p%phdos_d(:,1,nn),p%phstep,nq,d0)
                didinterp = didinterp .or. d0
             elseif (p%tmodel == tm_debye_einstein_v) then
                call read_freqg(word,ffreq)
                if (.not.allocated(p%freqg)) then
                   allocate(p%freqg(size(ffreq,1),phase_vmax))
                else
                   if (size(p%freqg,1) /= size(ffreq,1)) then
                      write (uout,'("In file: ",A)') trim(file(1:leng(file)))
                      write (uout,'("Line: ",A)') trim(line(1:leng(line)))
                      write (uout,'("Line number:",I5)') nn
                      call error('phase_init','Number of frequencies not consistent',faterr)
                   end if
                   if (nn > size(p%freqg,2)) &
                      call realloc(p%freqg,size(ffreq,1),2*nn)
                end if
                p%freqg(:,nn) = ffreq
                if (any(ffreq < 0d0)) then
                   nq = -count(ffreq < 0d0)
                else
                   nq = 0
                end if
             end if

             if (nq < 0) then
                ! negative freqs. -> deactivate for thermal
                p%dyn_active(nn) = .false.
                nq = abs(nq)
             end if
             numax = max(nq,numax)
             haveph = .true.
          else if (idum == icol_nef) then
             ! nefermi field
             ok = isreal(p%nefermi(nn),line,lp)
          else if (any(idum == icol_pol4)) then
             do j = 1, 4
                if (idum == icol_pol4(j)) then
                   ok = isreal(p%fel_cpol(j,nn),line,lp)
                   exit
                end if
             end do
             do j = 5, 8
                if (idum == icol_pol4(j)) then
                   ok = isreal(p%tsel_cpol(j-4,nn),line,lp)
                   exit
                end if
             end do
          else if (idum == icol_f0) then
             ! debye_extended: free energy
             ok = isreal(p%f0(nn),line,lp)
          else if (idum == icol_tde) then
             ! debye_extended: debye temperature
             ok = isreal(p%tde(nn),line,lp)
          else
             ! debye-extended: anharmonic coefficients
             ifound = 0
             do i = 1, p%tde_nanh
                if (idum == icol_anh(i)) ifound = i
             end do
             if (ifound > 0) then
                ok = isreal(p%tde_anh(ifound,nn),line,lp)
                if (.not.ok) exit
                cycle
             end if

             ! debye-extended: einstein coefficients
             ifound = 0
             do i = 1, p%tde_nein
                if (idum == icol_cein(i)) ifound = i
             end do
             if (ifound > 0) then
                ok = isreal(p%tde_cein(ifound,nn),line,lp)
                if (.not.ok) exit
                cycle
             end if

             ! debye-extended: einstein temperature
             ifound = 0
             do i = 1, p%tde_nein
                if (idum == icol_tein(i)) ifound = i
             end do
             if (ifound > 0) then
                ok = isreal(p%tde_tein(ifound,nn),line,lp)
                if (.not.ok) exit
                cycle
             end if

             ! interpolation field
             ifound = 0
             do i = 1, p%ninterp
                if (idum == icol_int(i)) ifound = i
             end do
             if (ifound > 0) then
                ok = isreal(p%interp(nn,ifound),line,lp)
                if (.not.ok) exit
                cycle
             end if

             ! none of the above -> skip this field
             word = getword(word,line,lp)
          end if
       end do

       ! some error was found during the read
       if (.not.ok) then
          if (uuin /= uin) then
             write (uout,'("In file: ",A)') trim(file(1:leng(file)))
             write (uout,'("Line: ",A)') trim(line(1:leng(line)))
             write (uout,'("Line number:",I5)') nn
          end if
          call error('phase_init','Error input, PHASE..ENDPHASE environment',faterr)
       end if

       if (.not.havev.or..not.havee) then
          if (uuin /= uin) then
             write (uout,'("In file: ",A)') trim(file(1:leng(file)))
          end if
          call error('phase_init','Error reading file: volume or energy missing.',faterr)
       end if

       if ((p%tmodel == tm_qhafull.or.p%tmodel == tm_debye_einstein_v) .and..not.haveph) then
          if (uuin /= uin) then
             write (uout,'("In file: ",A)') trim(file(1:leng(file)))
          end if
          call error('phase_init','Error reading file: phonon DOS or frequencies file missing.',faterr)
       end if
    end do

    if (didinterp) then
       write (msg,'("There were some interpolated phDOS, with step (input units) = ",F20.12)') p%phstep
       call error('phase_init',msg,noerr)
    end if

    if (uuin /= uin) then
       call fclose(uuin)
    endif

    ! reallocate
    p%nv = nn
    call phase_realloc_volume(p,nn,numax)

    ! if the temperature model is external fvib, read the external Fvib files
    if (p%tmodel == tm_externalfvib) &
       call read_externalfvib(extfvibfile,prefix,p%nv,p%v,p%nfvib_t,p%fvib_t,p%fvib_f,p%fvib_s,p%fvib_cv)

    ! apply energy shift
    p%e(1:p%nv) = p%e(1:p%nv) + eshift

    ! input units
    call phase_inputdata(p,zz)

  end subroutine phase_init

  !> Prepare all phases for calculation.
  subroutine setup_phases()
    use param, only: uout, warning, mline_fmt, au2gpa, faterr, ha2cm_1, ha2thz, units_f_cm1,&
       units_f_thz, format_string, format_string_header
    use tools, only: leng, error, realloc
    use evfunc, only: fit_polygibbs, fv0, fv1, fv2, fv3, fv4, punch_params
    use fit, only: fit_ev
    integer :: i, j
    character*(mline) :: msg
    character*(mline_fmt) :: fm
    real*8 :: pmaxmin, vmin_setv, vmax_setv
    logical :: anyactive, tim1, tsetexternal
    integer :: ierr, idx(1)
    real*8 :: vmin, vmax, v
    real*8 :: f1, f2, f3, f4, pt, bk, b1, b2
    integer :: strain, order

    real*8, parameter :: teps = 1d-5

    pmaxmin = 1d30
    vmin_setv = -1d30
    vmax_setv = 1d30
    anyactive = .false.
    tim1 = .true.
    tsetexternal = .false.
    do i = 1, nph
       ! sort
       call phase_sort(ph(i))

       ! check repeated points
       do j = 2, ph(i)%nv
          if (abs(ph(i)%v(j)-ph(i)%v(j-1)) < 1d-5) then
             write (msg,'("Phase ",A," -- repeated volumes at ",F12.6," bohr^3")') &
                trim(adjustl(ph(i)%name(1:leng(ph(i)%name)))), ph(i)%v(j)
             call error('setup_phases',msg,faterr)
          end if
       end do

       ! transform p(V) to E(V)
       if (ph(i)%pvdata) then
          if (tim1) then
             write (uout,'("* Results of fit to input p(V) data")')
             tim1 = .false.
          end if
          call fit_ev(fit_polygibbs, ph(i)%reg_mode, ph(i)%v, ph(i)%e,&
             ph(i)%npol, ph(i)%cpol, ierr, .true., ph(i)%nfix, ph(i)%idfix,&
             ph(i)%obelix)

          write (uout,'("+ Phase ",I2," (",A,")")') i, &
             trim(adjustl(ph(i)%name(1:leng(ph(i)%name))))
          fm = format_string_header( &
             (/1,ifmt_p,ifmt_v,ifmt_x,ifmt_p,ifmt_b,ifmt_bp,ifmt_bpp/),&
             (/1,6,9,4,10,6,2,10/))
          write (uout,fm) "#","p(GPa)","V(bohr^3)","V/V0","p_fit(GPa)","B(GPa)","Bp","Bpp(GPa-1)"
          fm = format_string((/ifmt_p,ifmt_v,ifmt_x,ifmt_p,ifmt_b,ifmt_bp,ifmt_bpp/),1)

          do j = ph(i)%nv, 1, -1
             v = ph(i)%v(j)
             f1 = fv1(ph(i)%fit_mode,v,ph(i)%npol,ph(i)%cpol)
             f2 = fv2(ph(i)%fit_mode,v,ph(i)%npol,ph(i)%cpol)
             f3 = fv3(ph(i)%fit_mode,v,ph(i)%npol,ph(i)%cpol)
             f4 = fv4(ph(i)%fit_mode,v,ph(i)%npol,ph(i)%cpol)
             pt = -f1 * au2gpa
             bk = v * f2 * au2gpa
             b1 = -(1+v*f3/f2)
             b2 = ((f3+v*f4)/f2**2 - v*f3**2/f2**3) / au2gpa
             write (uout,fm) ph(i)%e(j)*au2gpa,v,v/ph(i)%v(ph(i)%nv),pt,bk,b1,b2
          end do
          call punch_params(uout,ph(i)%fit_mode,ph(i)%npol,ph(i)%cpol)
          call phase_checkfiterr(ph(i),.true.)
          write (uout,'("#  Error RMS (GPa): ",1p,E15.7)') ph(i)%rms * au2gpa
          write (uout,'("#  max|error| (GPa): ",1p,E15.7)') ph(i)%maxdev * au2gpa
          write (uout,'("#  r2 of the fit: ",1p,E20.12)') ph(i)%r2
          write (uout,'("#  Akaike information criterion: ",1p,E20.12)') ph(i)%aic
          write (uout,'("#  Bayesian information (Schwarz) criterion: ",1p,E20.12)') ph(i)%bic

          vmin = ph(i)%v(1)
          vmax = ph(i)%v(ph(i)%nv) + facexpand * (ph(i)%v(ph(i)%nv)-ph(i)%v(1))
          write (uout,'("! Generating E(V) data using ",I3," points in the V-range: ",F12.4,X,F12.4/)') &
             newpts, vmin, vmax
          ph(i)%nv = newpts
          call phase_realloc_volume(ph(i),ph(i)%nv)

          do j = 0, newpts-1
             v = vmin + real(j,8)/real(newpts-1,8) * (vmax-vmin)
             ph(i)%v(j+1) = v
             ph(i)%e(j+1) = fv0(ph(i)%fit_mode,v,ph(i)%npol,ph(i)%cpol)
          end do
       end if

       ! check convexity (warnings)
       call phase_checkconvex(ph(i),.false.,ierr)

       ! check that minimum exists
       idx = minloc(ph(i)%e)
       if (idx(1) == 1 .or. idx(1) == ph(i)%nv) then
          write (msg,'("Phase ",A,": E(V) data: minimum not found")')&
             trim(adjustl(ph(i)%name(1:leng(ph(i)%name))))
          call error('setup_phases',msg,warning)
          ph(i)%staticmin = .false.
       else
          ph(i)%staticmin = .true.
       end if

       ! static fit
       call fit_ev(ph(i)%fit_mode, ph(i)%reg_mode, ph(i)%v, ph(i)%e, ph(i)%npol,&
          ph(i)%cpol, ierr, .false., ph(i)%nfix, ph(i)%idfix, ph(i)%obelix, ph(i)%pfit)
       if (ierr > 0) then
          write (uout,'("Phase ",A)') trim(adjustl(ph(i)%name(1:leng(ph(i)%name))))
          call error('setup_phases','E(V) fit: minimum not found',faterr)
       end if
       call phase_checkfiterr(ph(i),.false.)

       ! pmin and pmax
       ph(i)%pmax = -fv1(ph(i)%fit_mode,ph(i)%v(2),ph(i)%npol,ph(i)%cpol) * au2gpa
       ph(i)%pmin = -fv1(ph(i)%fit_mode,ph(i)%v(ph(i)%nv),ph(i)%npol,ph(i)%cpol) * au2gpa
       pmaxmin = min(pmaxmin,ph(i)%pmax)
       vmin_setv = max(vmin_setv,ph(i)%v(1))
       vmax_setv = min(vmax_setv,ph(i)%v(ph(i)%nv))

       ! tmodel-dependent setup
       if (ph(i)%tmodel == tm_debye_einstein .or. ph(i)%tmodel == tm_debye_einstein_v) then
          ! check freqg is allocated
          if (.not.allocated(ph(i)%freqg)) &
             call error('setup_phases','Debye-Einstein requires frequencies at Gamma',faterr)

          ! transform units
          if (ph(i)%units_f == units_f_cm1) then
             ph(i)%freqg = ph(i)%freqg / ha2cm_1
          else if (ph(i)%units_f == units_f_thz) then
             ph(i)%freqg = ph(i)%freqg / ha2thz
          end if

          ! check nfreq vs. vfree * Z
          ph(i)%nfreq = size(ph(i)%freqg,1)
          if (ph(i)%nfreq /= 3*(vfree*ph(i)%z)-3) then
             write (uout,'("* No. of read frequencies: ",I6)') ph(i)%nfreq
             write (uout,'("* No. of expected frequencies 3*(vfree*Z)-3: ",I6)') 3*(vfree*nint(ph(i)%z))-3
             call error('gibbs2','nfreq /= 3*(vfree*z)-3, check input',faterr)
          end if
       else if (ph(i)%tmodel == tm_qhafull) then
          call phase_phdos(ph(i))
       else if (ph(i)%tmodel == tm_externalfvib) then
          ! tsetexternal, override temperature list and check consistency
          if (tdefault .or. .not.tdefault.and..not.tsetexternal) then
             if (.not.tdefault) &
                call error('setup_phases','externaltvib overrides user-defined temperature list',warning)
             tdefault = .false.
             tsetexternal = .true.
             nts = ph(i)%nfvib_t
             if (allocated(tlist)) deallocate(tlist)
             allocate(tlist(nts))
             tlist = ph(i)%fvib_t
          else if (tsetexternal) then
             if (nts /= ph(i)%nfvib_t) &
                call error('setup_phases','inconsistent num. of temperatures in two externalfvib phases',faterr)
             if (any(abs(tlist - ph(i)%fvib_t) > teps)) &
                call error('setup_phases','inconsistent temperatures in two externalfvib phases',faterr)
          end if
       end if

       ! entropy and ThetaD ratio fit modes
       if (ph(i)%fit_mode > 10000) then
          strain = (ph(i)%fit_mode - fit_strain * 10000) / 100
          order =  (ph(i)%fit_mode - fit_strain * 10000 - 100 * strain)
          if (order == 0) then
             ph(i)%sfit_mode = ph(i)%fit_mode
             ph(i)%tdfit_mode = ph(i)%fit_mode
          else
             ph(i)%sfit_mode = fit_strain * 10000 + fit_strain_x1 * 100 + min(ph(i)%nv - 5, sfit_order)
             ph(i)%tdfit_mode = fit_strain * 10000 + fit_strain_x1 * 100 + min(ph(i)%nv - 5, tdfit_order)
          end if
       else
          ph(i)%sfit_mode = fit_strain * 10000 + fit_strain_x1 * 100 + min(ph(i)%nv - 5, sfit_order)
          ph(i)%tdfit_mode = fit_strain * 10000 + fit_strain_x1 * 100 + min(ph(i)%nv - 5, tdfit_order)
       end if
    end do

    ! set pressure range if not given in input
    if (pdefault) then
       nps = 100
       allocate(plist(nps))
       plist(1) = 0d0
       pstep = max(pmaxmin,500d0) / real(nps-1,8)
       do i = 2, nps
          plist(i) = plist(i-1) + pstep
       end do
    else
       if (.not.allocated(plist)) then
          if (nps > 0) then
             allocate(plist(nps))
             plist(1) = 0d0
             pstep = pmaxmin / real(nps-1,8)
             do i = 2, nps
                plist(i) = plist(i-1) + pstep
             end do
          else
             nps = floor(pmaxmin / pstep) + 1
             allocate(plist(nps))
             plist(1) = 0d0
             do i = 2, nps
                plist(i) = plist(i-1) + pstep
             end do
          end if
       end if
    end if
    call inplace_sort(plist(1:nps))

    ! Make sure that the first pressure is p=0
    if (abs(plist(1)) > 1d-5) then
       call realloc(plist,nps+1)
       plist(2:nps+1) = plist(1:nps)
       plist(1) = 0d0
       nps = nps + 1
       call error('setup_phases','The first element of the pressure list must be zero.',warning)
    end if

    write (uout,'("* Pressure range examined")')
    write (uout,'("  Minimum p_max across all phases (GPa): ",F12.3)') pmaxmin
    write (uout,'("  Pressure range (GPa): ",F12.3," -> ",F12.3)') &
       plist(1), plist(nps)
    write (uout,'("  Number of pressure points: ",I6)') nps
    write (uout,*)

    ! set volume range if not given in input
    if (.not.vdefault) then
       if (allocated(vlist)) deallocate(vlist)
       if (nvs > 0) then
          allocate(vlist(nvs))
          vlist(1) = vmin_setv
          vstep = (vmax_setv-vmin_setv) / real(nvs-1,8)
          do i = 2, nvs
             vlist(i) = vlist(i-1) + vstep
          end do
       else if (nvs == -1) then
          nvs = floor((vmax_setv - vmin_setv) / vstep) + 1
          allocate(vlist(nvs))
          vlist(1) = vmin_setv
          do i = 2, nvs
             vlist(i) = vlist(i-1) + vstep
          end do
       else
          ! let each phase handle it
          nvs = 0
       end if

       if (allocated(vlist)) then
          call inplace_sort(vlist(1:nvs))

          write (uout,'("* Volume range examined")')
          write (uout,'("  Volume range (bohr^3): ",F14.5," -> ",F14.5)') &
             vlist(1), vlist(nvs)
          write (uout,'("  Number of V points: ",I6)') nvs
          write (uout,*)
       else
          write (uout,'("* Input volumes used for the EOS")')
       end if
    end if

  end subroutine setup_phases

  !> Calculate the static properties (V0, B0, and E0) for each phase.
  subroutine props_staticeq()
    use fit, only: fit_pshift
    use tools, only: leng, error, realloc
    use param, only: uout, warning, mline_fmt, format_string, ifmt_p, ifmt_v, ifmt_x,&
       ifmt_integer5
    integer :: i, j, ipcut
    character*(mline) :: msg
    real*8 :: vk, bk, ek, gk
    integer :: ierr
    logical :: changedplist
    character*(mline_fmt) :: fm

    write (uout,'("* Calculating static properties on pressure grid")')
    changedplist = .false.
    fm = format_string((/ifmt_integer5,ifmt_p,ifmt_v,ifmt_x/),1)

    ! check that it is possible to minimize G(static) = E(static) + pV vs. V
    ! get static equilibrium properties.
    do i = 1, nph
       write (uout,'("+ Phase ",I2," (",A,")")') i, trim(adjustl(ph(i)%name(1:leng(ph(i)%name))))
       write (uout,'("# Volumes per formula unit.")')
       write (uout,'("# xV indicates position of V on the volume grid (0 = beginning, 1 = end)")')
       write (uout,'("# Num    p(GPa)    V(bohr^3)    xV")')

       ! check that it is possible to minimize G(static) = E(static) + pV vs. V
       ipcut = nps+1
       do j = 1, nps
          call fit_pshift(ph(i)%fit_mode, ph(i)%v,plist(j),ph(i)%npol,ph(i)%cpol,vk,bk,ek,gk,ierr)
          if (ierr > 0) then
             write (uout,'(X,I5,X,F10.4,X," -- minimum outside volume grid --")') j, plist(j)
             if (plist(j) > ph(i)%pmax.and..not.ph(i)%extend) then
                ipcut = min(j-1,ipcut)
                changedplist = .true.
             end if
          else
             write (uout,fm) j, plist(j), vk, (vk - ph(i)%v(1)) / (ph(i)%v(ph(i)%nv) - ph(i)%v(1))
          end if

          if (j == 1) then
             ph(i)%veq_static = vk
             ph(i)%eeq_static = gk
             ph(i)%beq_static = bk
          end if
       end do
    end do

    if (changedplist) then
       nps = ipcut
       call realloc(plist,nps)
       call error('props_staticeq',"Pressure range reduced because no minimum in static curve",warning)
       write (uout,'("  New pressure range (GPa): ",F12.3," -> ",F12.3)') &
          plist(1), plist(nps)
       write (uout,'("  New number of p points: ",I6)') nps
       write (uout,*)
    end if

  end subroutine props_staticeq

  !> In-place sort of array a using qcksort.
  subroutine inplace_sort(a)
    use tools, only: qcksort
    real*8, intent(inout) :: a(:)

    integer :: i, idx(size(a))

    do i = 1, size(a)
       idx(i) = i
    end do
    call qcksort(a,idx,1,size(a))
    a = a(idx)

  end subroutine inplace_sort

  !!!!!!!!!!!!!!!!!!!!!!!!
  !xx! type(phase) methods
  !!!!!!!!!!!!!!!!!!!!!!!!

  ! Sort the volume and energy arrays and re-arrange all the other
  ! arrays.
  subroutine phase_sort(p)
    use tools, only: qcksort
    type(phase), intent(inout) :: p

    integer, allocatable :: idx(:)
    integer :: i

    if (.not.allocated(p%v).or..not.allocated(p%e)) return

    ! volume sort
    allocate(idx(p%nv))
    do i = 1, p%nv
       idx(i) = i
    end do
    call qcksort(p%v,idx,1,p%nv)
    p%v = p%v(idx)
    p%e = p%e(idx)
    if (allocated(p%td)) p%td = p%td(idx)
    if (allocated(p%dyn_active)) p%dyn_active = p%dyn_active(idx)
    if (allocated(p%interp)) p%interp = p%interp(idx,:)
    if (allocated(p%phdos_d)) p%phdos_d = p%phdos_d(:,:,idx)
    if (allocated(p%nefermi)) p%nefermi = p%nefermi(idx)
    if (allocated(p%fel_cpol)) p%fel_cpol = p%fel_cpol(:,idx)
    if (allocated(p%tsel_cpol)) p%tsel_cpol = p%tsel_cpol(:,idx)
    if (allocated(p%fvib_f)) p%fvib_f = p%fvib_f(idx,:)
    if (allocated(p%fvib_s)) p%fvib_s = p%fvib_s(idx,:)
    if (allocated(p%fvib_cv)) p%fvib_cv = p%fvib_cv(idx,:)
    if (allocated(p%f0)) p%f0 = p%f0(idx)
    if (allocated(p%tde)) p%tde = p%tde(idx)
    if (allocated(p%tde_anh)) p%tde_anh = p%tde_anh(:,idx)
    if (allocated(p%tde_cein)) p%tde_cein = p%tde_cein(:,idx)
    if (allocated(p%tde_tein)) p%tde_tein = p%tde_tein(:,idx)

    ! temperature sort in external fvib
    if (p%tmodel == tm_externalfvib .and. allocated(p%fvib_t) .and. allocated(p%fvib_f)&
       .and. allocated(p%fvib_s) .and. allocated(p%fvib_cv)) then
       deallocate(idx)
       allocate(idx(p%nfvib_t))
       do i = 1, p%nfvib_t
          idx(i) = i
       end do
       call qcksort(p%fvib_t,idx,1,p%nfvib_t)
       p%fvib_t = p%fvib_t(idx)
       p%fvib_f = p%fvib_f(:,idx)
       p%fvib_s = p%fvib_s(:,idx)
       p%fvib_cv = p%fvib_cv(:,idx)
    end if

  end subroutine phase_sort

  ! Realloc volume arrays in phase. If doshift, move the data to the
  ! right and shrink. If numax, increase the number of frequencies in
  ! the QHA arrays.
  subroutine phase_realloc_volume(p,n,numax,lshift)
    use tools, only: realloc
    type(phase), intent(inout) :: p
    integer, intent(in) :: n
    integer, intent(in), optional :: numax
    logical, intent(in), optional :: lshift

    integer :: n2, nv
    logical :: doshift

    doshift = .false.
    if (present(lshift)) doshift = lshift

    nv = size(p%v)

    if (doshift) p%v(1:n) = p%v(nv-n+1:nv)
    call realloc(p%v,n)

    if (doshift) p%e(1:n) = p%e(nv-n+1:nv)
    call realloc(p%e,n)

    if (allocated(p%td)) then
       if (doshift) p%td(1:n) = p%td(nv-n+1:nv)
       call realloc(p%td,n)
    end if

    if (allocated(p%dyn_active)) then
       if (doshift) p%dyn_active(1:n) = p%dyn_active(nv-n+1:nv)
       call realloc(p%dyn_active,n)
    end if

    if (allocated(p%interp)) then
       if (doshift) p%interp(1:n,:) = p%interp(nv-n+1:nv,:)
       call realloc(p%interp,n,size(p%interp,2))
    end if

    if (allocated(p%phdos_d)) then
       n2 = size(p%phdos_d,1)
    else
       n2 = 0
    end if
    if (present(numax)) then
       if (numax > 0) then
          n2 = numax
          if (allocated(p%phdos_f)) call realloc(p%phdos_f,numax)
       end if
    end if

    if (allocated(p%phdos_d)) then
       if (doshift) then
          p%phdos_d(:,:,1:n) = p%phdos_d(:,:,nv-n+1:nv)
       end if
       call realloc(p%phdos_d,n2,size(p%phdos_d,2),n)
    end if

    if (allocated(p%nefermi)) then
       if (doshift) p%nefermi(1:n) = p%nefermi(nv-n+1:nv)
       call realloc(p%nefermi,n)
    end if

    if (allocated(p%fel_cpol)) then
       if (doshift) p%fel_cpol(:,1:n) = p%fel_cpol(:,nv-n+1:nv)
       call realloc(p%fel_cpol,size(p%fel_cpol,1),n)
    end if

    if (allocated(p%tsel_cpol)) then
       if (doshift) p%tsel_cpol(:,1:n) = p%tsel_cpol(:,nv-n+1:nv)
       call realloc(p%tsel_cpol,size(p%tsel_cpol,1),n)
    end if

    if (allocated(p%fvib_f)) then
       if (doshift) p%fvib_f(1:n,:) = p%fvib_f(nv-n+1:nv,:)
       call realloc(p%fvib_f,n,size(p%fvib_f,2))
    end if

    if (allocated(p%fvib_s)) then
       if (doshift) p%fvib_s(1:n,:) = p%fvib_s(nv-n+1:nv,:)
       call realloc(p%fvib_s,n,size(p%fvib_s,2))
    end if

    if (allocated(p%fvib_cv)) then
       if (doshift) p%fvib_cv(1:n,:) = p%fvib_cv(nv-n+1:nv,:)
       call realloc(p%fvib_cv,n,size(p%fvib_cv,2))
    end if

    if (allocated(p%tde)) then
       if (doshift) p%tde(1:n) = p%tde(nv-n+1:nv)
       call realloc(p%tde,n)
    end if

    if (allocated(p%f0)) then
       if (doshift) p%f0(1:n) = p%f0(nv-n+1:nv)
       call realloc(p%f0,n)
    end if

    if (allocated(p%tde_anh)) then
       if (doshift) p%tde_anh(:,1:n) = p%tde_anh(:,nv-n+1:nv)
       call realloc(p%tde_anh,size(p%tde_anh,1),n)
    end if

    if (allocated(p%tde_cein)) then
       if (doshift) p%tde_cein(:,1:n) = p%tde_cein(:,nv-n+1:nv)
       call realloc(p%tde_cein,size(p%tde_cein,1),n)
    end if

    if (allocated(p%tde_tein)) then
       if (doshift) p%tde_tein(:,1:n) = p%tde_tein(:,nv-n+1:nv)
       call realloc(p%tde_tein,size(p%tde_tein,1),n)
    end if

  end subroutine phase_realloc_volume

  !> Write information about phase i to the standard output.
  subroutine phase_popinfo(p,i)
    use param, only: uout, units_v_bohr3, units_v_ang3, units_e_ha, units_e_ev, units_e_ry,&
       units_p_au, units_p_gpa, units_f_ha, units_f_cm1, units_f_thz, units_e_ha, &
       units_e_ev, units_e_ry, ha2kjmol, warning, faterr
    use tools, only: leng, error
    type(phase), intent(in) :: p
    integer, intent(in) :: i

    integer :: k, ic
    logical :: dowarn_neg, dowarn_min
    character*(mline) :: msg

    dowarn_neg = .false.
    dowarn_min = .false.

    write (uout,'("+ Phase ",I2," (",A,")")') i, &
       trim(adjustl(p%name(1:leng(p%name))))
    write (uout,'("  Number of volume points: ",I4)') p%nv
    write (uout,'("  Number of vfree units (Z): ",F12.3)') p%z
    write (uout,'("  p(V) input data? ",L1)') p%pvdata
    write (uout,'("  Pressure range (GPa): ",F12.3," -> ",F12.3)') &
       p%pmin, p%pmax
    write (uout,'("  Number of interpolated fields : ",I4)') p%ninterp

    ! input units
    write (uout,'("  Input units: ")')
    select case(p%units_v)
    case(units_v_bohr3)
       write (uout,'("    Volume : bohr^3")')
    case(units_v_ang3)
       write (uout,'("    Volume : ang^3")')
    end select
    select case(p%units_e)
    case(units_e_ha)
       write (uout,'("    Energy : Hartree")')
    case(units_e_ev)
       write (uout,'("    Energy : eV")')
    case(units_e_ry)
       write (uout,'("    Energy : Ry")')
    end select
    select case(p%units_p)
    case(units_p_au)
       write (uout,'("    Pressure : Atomic units")')
    case(units_p_gpa)
       write (uout,'("    Pressure : GPa")')
    end select
    select case(p%units_f)
    case(units_f_ha)
       write (uout,'("    Frequency : Hartree")')
    case(units_f_cm1)
       write (uout,'("    Frequency : cm^(-1)")')
    case(units_f_thz)
       write (uout,'("    Frequency : Thz")')
    end select
    select case(p%eunits_e)
    case(units_e_ha)
       write (uout,'("    DOS energy : Hartree")')
    case(units_e_ev)
       write (uout,'("    DOS energy : eV")')
    case(units_e_ry)
       write (uout,'("    DOS energy : Ry")')
    end select
    write (uout,'("  Output units are atomic units (Ha), except where noted.")')

    ! static energy fit
    write (uout,'("  First/last volume (bohr^3): ",1p,2(E20.12,2X))') &
       p%v(1), p%v(p%nv)
    write (uout,'("  First/last energy (Ha): ",1p,2(E20.12,2X))') &
       p%e(1), p%e(p%nv)
    if (allocated(p%td)) then
       write (uout,'("  First/last Debye temp. (K): ",1p,2(E20.12,2X))') &
          p%td(1), p%td(p%nv)
    end if
    write (uout,'("  Poisson ratio from input (sigma): ",F14.6)') p%poisson
    write (uout,'("  Poisson function, f(sigma): ",F14.6)') p%pofunc
    if (allocated(p%freqg)) then
       write (uout,'("  Number of freq. at G (p=0) : ",I4)') p%nfreq
    end if
    if (.not.p%staticmin) then
       write (uout,'("  Beware!! Static properties are EXTRAPOLATED ")')
       dowarn_min = .true.
    end if
    select case(p%scaltype)
    case(scal_noscal)
       write (uout,'("  Correction of static energy: no correction")')
    case(scal_pshift)
       write (uout,'("  Correction of static energy: pshift (Vexp0=",1p,E12.4,")")') &
          p%vscal
    case(scal_bpscal)
       write (uout,'("  Correction of static energy: bpscal (Vexp0=",1p,E12.4," Bexp0=",E12.4")")') &
          p%vscal, p%bscal
    case(scal_apbaf)
       write (uout,'("  Correction of static energy: apbaf (Vexp0=",1p,E12.4,")")') &
          p%vscal
    case(scal_use)
       write (uout,'("  Correction of static energy: from phase ",I3)') p%iscal
    end select
    write (uout,'("  Energy fit mode: ",I8)') p%fit_mode
    write (uout,'("  S(V) fit mode: ",I8)') p%sfit_mode
    write (uout,'("  ThetaD(V) fit mode: ",I8)') p%tdfit_mode
    write (uout,'("  omega(V) fit mode: ",I8)') omega_fitmode
    write (uout,'("  Number of fixed fit parameters: ",I3)') p%nfix
    if (p%nfix > 0) then
       write (uout,'("  Id. of fixed parameters: ",99(I2,X))') p%idfix(1:p%nfix)
       write (uout,'("  Val. of fixed parameters: ",1p,99(E12.4,X))') p%obelix(1:p%nfix)
    end if
    write (uout,'("  Static equilibrium volume (bohr^3): ",F20.10)') p%veq_static
    write (uout,'("  Static equilibrium energy (Ha): ",F20.10)') p%eeq_static
    write (uout,'("  Static equilibrium energy (kJ/mol): ",F20.10)') p%eeq_static * ha2kjmol
    write (uout,'("  Static bulk modulus (GPa): ",F15.6)') p%beq_static
    write (uout,'("  Static EOS fit, error RMS (Ha): ",1p,E15.7)') p%rms
    write (uout,'("  Static EOS fit, max|error| (Ha): ",1p,E15.7)') p%maxdev
    write (uout,'("  r2 of the fit: ",1p,E20.12)') p%r2
    write (uout,'("  Akaike information criterion: ",1p,E20.12)') p%aic
    write (uout,'("  Bayesian information (Schwarz) criterion: ",1p,E20.12)') p%bic

    ! dynamic info
    select case(p%tmodel)
    case(tm_static)
       write (uout,'("  Temperature model: static.")')
    case(tm_debye_input)
       write (uout,'("  Temperature model: Debye, Td read from input.")')
    case(tm_debye)
       write (uout,'("  Temperature model: Debye, Td from static B(V).")')
    case(tm_debye_einstein)
       write (uout,'("  Temperature model: Debye-Einstein (one set of frequencies at Gamma).")')
    case(tm_debye_einstein_v)
       write (uout,'("  Temperature model: Debye-Einstein (frequencies at Gamma at every volume).")')
    case(tm_debye_poisson_input)
       write (uout,'("  Temperature model: Debye with Poisson ratio in input.")')
    case(tm_qhafull)
       write (uout,'("  Temperature model: QHA (phonon DOS).")')
       write (uout,'("  Frequency step: ",1p,E12.4)') p%phstep
    case(tm_debyegrun)
       write (uout,'("  Temperature model: Debye-Gruneisen with gamma = a + b*B'' .")')
       write (uout,'("  Gamma a coefficient: ",F12.3)') p%a_grun
       write (uout,'("  Gamma b coefficient: ",F12.3)') p%b_grun
    end select
    if (all(p%dyn_active)) then
       write (uout,'("  All data points are ACTIVE for dynamic calculations")')
    else
       write (uout,'("  Negative frequencies, volumes INACTIVE in dynamic calc. [id (volume)]:")')
       msg = ""
       ic = 0
       do k = 1, p%nv
          if (.not.p%dyn_active(k)) then
             ic = ic + 1
             write (msg,'(A,X,I4," (",1p,E12.4,0p,")")') trim(msg), k, p%v(k)
          end if
          if (ic == 3) then
             ic = 0
             write (uout,'(A)') trim(msg)
             msg = ""
          end if
       end do
       if (ic /= 0) write (uout,'(A)') trim(msg)
       dowarn_neg = .true.

       if (count(p%dyn_active) < 5 .and. p%tmodel /= tm_static) then
          write (uout,'(/"The number of active volumes for dynamic calc. is: ", I1 )') count(p%dyn_active)
          write (uout,'("Please, check the temperature model and source data.")')
          call error('phase_popinfo','Not enough active volumes for dynamic calc. fits (< 5).',faterr)
       end if
    end if
    if (p%emodel /= em_no) then
       select case(p%emodel)
       case(em_sommerfeld)
          if (p%efree) then
             write (uout,'("  Electronic contribution model: Sommerfeld with free electron N(Ef).")')
          else
             write (uout,'("  Electronic contribution model: Sommerfeld with input N(Ef).")')
          end if
          write (uout,'("  Number of conduction electrons: ",I3)') p%nelec
       case(em_pol4)
          write (uout,'("  Electronic contribution model: F_el(T) and -TS_el(T) fitted with 4th order polynomials.")')
       end select
    end if

    ! Warnings
    if (dowarn_neg) then
       call error('phase_popinfo','Negative frequencies found.',warning)
    end if
    if (dowarn_min) then
       call error('phase_popinfo','Static equilibrium volume out of grid bounds.',warning)
    end if
    if (p%pmax > warn_pmax) then
       write (msg,'("Max. pressure (",F12.3,") exceeds ",F12.3)') ph(i)%pmax, warn_pmax
       call error('phase_popinfo',msg,warning)
    end if

    write (uout,*)

  end subroutine phase_popinfo

  !> Remove points with negative E(V) curvature from the static E(V)
  !> curve.
  subroutine phase_checkconvex(p,usecpol,ierr)
    use param, only: uout, noerr, warning
    use tools, only: leng, error
    use evfunc, only: fv2
    type(phase), intent(inout) :: p
    logical, intent(in) :: usecpol
    integer, intent(out) :: ierr

    integer :: i, iremove
    logical :: conv, firstconv, oconv
    integer :: inflidx(p%nv), ninf
    character*(mline) :: msg
    real*8 :: f2

    ierr = 0
    firstconv = .true.
    oconv = .false.
    iremove = 0
    ninf = 0
    if (.not.usecpol) then
       do i = 2, p%nv-1
          conv = (p%e(i)-p%e(i-1))/(p%v(i)-p%v(i-1)) < &
             (p%e(i+1)-p%e(i))/(p%v(i+1)-p%v(i))
          if (firstconv) then
             if (conv) then
                firstconv = .false.
                oconv = conv
             else
                iremove = iremove + 1
             end if
          else
             if (conv .neqv. oconv) then
                ninf = ninf + 1
                inflidx(ninf) = i
                oconv = conv
             end if
          end if
       end do
    else
       do i = p%nv, 1, -1
          f2 = fv2(p%fit_mode, p%v(i), p%npol, p%cpol)
          if (f2 < 0) then
             iremove = iremove + 1
          else
             exit
          end if
       end do
    end if

    if (iremove > 0) then
       write (msg,&
          '("Removed ",I3," points at the beginning of E(V), phase ",A)')&
          iremove, trim(adjustl(p%name(1:leng(p%name))))
       call error('phase_checkconvex',msg,warning)

       call phase_realloc_volume(p,p%nv-iremove,-1,.true.)
       p%nv = p%nv - iremove
       ierr = 1
    end if
    if (ninf > 0) then
       write (msg,&
          '("Found ",I3," inflection points in (E,V) of phase ",A)')&
          ninf, trim(adjustl(p%name(1:leng(p%name))))
       call error('phase_checkconvex',msg,noerr)
       do i = 1, ninf
          write (uout,'("  Near (E,V) point number ",I2," : ",F15.7)') &
             inflidx(i), p%v(inflidx(i))
       end do
       write (uout,*)
    end if

  end subroutine phase_checkconvex

  !> Convert input units and apply the Z keyword (zz).
  subroutine phase_inputdata(p,zz)
    use param, only: au2gpa, bohr2angstrom, ha2ev, units_e_ev, units_e_ry, units_p_gpa,&
       units_v_ang3
    type(phase), intent(inout) :: p
    real*8, intent(in) :: zz

    ! volume
    p%v = p%v / zz
    p%vscal = p%vscal / zz
    if (p%units_v == units_v_ang3) then
       p%v = p%v / bohr2angstrom**3
       p%vscal = p%vscal / bohr2angstrom**3
    end if

    ! pressure is stored in p%e until fit.
    if (.not.p%pvdata) then
       p%e = p%e / zz
       if (p%units_e == units_e_ev) then
          p%e = p%e / ha2ev
       else if (p%units_e == units_e_ry) then
          p%e = p%e / 2d0
       end if
    else
       if (p%units_p == units_p_gpa) then
          p%e = p%e / au2gpa
       end if
    end if

    ! energy
    if (allocated(p%fel_cpol)) p%fel_cpol = p%fel_cpol / zz
    if (allocated(p%tsel_cpol)) p%tsel_cpol = p%tsel_cpol / zz
    if (allocated(p%fvib_f)) p%fvib_f = p%fvib_f / zz
    if (allocated(p%fvib_s)) p%fvib_s = p%fvib_s / zz
    if (allocated(p%fvib_cv)) p%fvib_cv = p%fvib_cv / zz
    if (allocated(p%f0)) p%f0 = p%f0 / zz

    if (p%units_e == units_e_ry) then
       if (allocated(p%fel_cpol)) p%fel_cpol = p%fel_cpol / 2d0
       if (allocated(p%tsel_cpol)) p%tsel_cpol = p%tsel_cpol / 2d0
       if (allocated(p%fvib_f)) p%fvib_f = p%fvib_f / 2d0
       if (allocated(p%fvib_s)) p%fvib_s = p%fvib_s / 2d0
       if (allocated(p%fvib_cv)) p%fvib_cv = p%fvib_cv / 2d0
       if (allocated(p%f0)) p%f0 = p%f0 / 2d0
    else if (p%units_e == units_e_ev) then
       if (allocated(p%fel_cpol)) p%fel_cpol = p%fel_cpol / ha2ev
       if (allocated(p%tsel_cpol)) p%tsel_cpol = p%tsel_cpol / ha2ev
       if (allocated(p%fvib_f)) p%fvib_f = p%fvib_f / ha2ev
       if (allocated(p%fvib_s)) p%fvib_s = p%fvib_s / ha2ev
       if (allocated(p%fvib_cv)) p%fvib_cv = p%fvib_cv / ha2ev
       if (allocated(p%f0)) p%f0 = p%f0 / ha2ev
    end if

    ! edos input units
    if (p%eunits_e == units_e_ry) then
       if (allocated(p%nefermi)) p%nefermi = p%nefermi * 2d0
    else if (p%eunits_e == units_e_ev) then
       if (allocated(p%nefermi)) p%nefermi = p%nefermi * ha2ev
    end if

  end subroutine phase_inputdata

  !> Calculate static E(V) fit errors. If ispv, assume the data is p(V).
  subroutine phase_checkfiterr(p,ispv)
    use evfunc, only: fv0, fv1
    type(phase), intent(inout) :: p
    logical, intent(in) :: ispv

    integer :: i, ntot
    real*8 :: d, ymean, y2mean, yvar, sse

    ymean = sum(p%e) / p%nv
    y2mean = sum(p%e**2) / p%nv
    yvar = y2mean - ymean**2

    sse = 0d0
    p%maxdev = -1d30
    do i = 1, p%nv
       if (ispv) then
          d = p%e(i)+fv1(p%fit_mode,p%v(i),p%npol,p%cpol)
       else
          d = p%e(i)-fv0(p%fit_mode,p%v(i),p%npol,p%cpol)
       end if
       sse = sse + d**2
       if (d > p%maxdev) then
          p%maxdev = d
       end if
    end do
    p%rms = sqrt(sse / p%nv)

    if (abs(yvar) < 1d-12) then
       p%r2 = 1
    else
       p%r2 = 1 - sse / yvar
    end if
    ntot = p%npol - p%nfix
    p%aic = 2 * ntot + p%nv * log(max(sse,1d-40))
    p%bic = log(real(p%nv,8)) * ntot + p%nv * log(max(sse,1d-40))

  end subroutine phase_checkfiterr

  !> Process the phDOS data for phase p. Convert units, renormalize,
  !> fit spline, etc.
  subroutine phase_phdos(p)
    use param, only: ha2cm_1, ha2thz, units_f_cm1, units_f_thz, warning
    use tools, only: quad1, error, realloc
    type(phase), intent(inout) :: p

    integer :: iv
    integer :: nfreq, i, ini, numax
    real*8 :: sumn
    character*(mline) :: msg

    real*8, parameter :: deps = 1d-6
    real*8, parameter :: fsmallcrit = 0.1d0 / ha2cm_1

    if (p%tmodel /= tm_qhafull) return

    ! convert input units
    if (p%units_f == units_f_cm1) then
       p%phdos_f = p%phdos_f / ha2cm_1
       p%phdos_d = p%phdos_d * ha2cm_1
       p%phstep = p%phstep / ha2cm_1
    elseif (p%units_f == units_f_thz) then
       p%phdos_f = p%phdos_f / ha2thz
       p%phdos_d = p%phdos_d * ha2thz
       p%phstep = p%phstep / ha2thz
    end if

    ! number of frequencies
    nfreq = size(p%phdos_f)

    ! clean zero-frequencies at the beginning
    do i = 1, nfreq
       if (p%phdos_f(i) > fsmallcrit) then
          ini = i
          exit
       end if
    end do

    ! clean zero-DOS at the end
    do i = nfreq,1,-1
       if (any(abs(p%phdos_d(i,1,:)) > deps)) then
          numax = i
          exit
       end if
    end do

    p%phdos_f(1:numax-ini+1) = p%phdos_f(ini:numax)
    p%phdos_d(1:numax-ini+1,:,:) = p%phdos_d(ini:numax,:,:)

    ! reallocate
    nfreq = numax-ini+1
    call realloc(p%phdos_f,nfreq)
    call realloc(p%phdos_d,nfreq,4,size(p%phdos_d,3))

    ! ! check if frequency step is constant
    ! step = p%phdos_f(2) - p%phdos_f(1)
    ! do i = 3, nfreq
    !    if (p%phdos_f(i) < 0d0) exit
    !    step1 = p%phdos_f(i) - p%phdos_f(i-1)
    !    if (abs(step-step1)/abs(step+step1) > stepcrit) then
    !       write (uout,'("Phase ",I2," (",A,")")') i, &
    !          trim(adjustl(p%name(1:leng(p%name))))
    !       write (uout,'("step=",E12.4,X,"step1=",E12.4,X,"stepcrit=",E12.4)') step, step1, stepcrit
    !       call error('phase_phdos','Frequency step must be the same for all volumes',faterr)
    !       exit
    !    end if
    ! end do
    ! p%phstep = step

    ! normalization
    do iv = 1, p%nv
       if (.not.p%dyn_active(iv)) cycle
       sumn = quad1(p%phdos_f,p%phdos_d(:,1,iv),p%phstep)
       if (renormalize) then
          if (abs(sumn - real(3*vfree*p%z,8)) > 1d-2) then
             write (msg,'(" Volume num. ",I3,1p," [",E12.4,"] phDOS renormalized from ",E20.10," to ",E20.10)') &
                iv, p%v(iv), sumn, 3*vfree*p%z
             call error('phase_phdos',msg,warning)
          end if
          p%phdos_d(:,1,iv) = p%phdos_d(:,1,iv) / sumn * (3d0 * vfree)
       else
          if (abs(sumn - real(3*vfree*p%z,8)) > 1d-2) then
             write (msg,'(" Volume num. ",I3,1p," [",E12.4,"] phDOS norm. = ",E20.10," instead of ",E20.10)') &
                iv, p%v(iv), sumn, 3*vfree*p%z
          end if
       end if
    end do

    ! spline fit to phdos(V)
    p%phdos_d(:,2:4,:) = 0d0
    do iv = 1, nfreq
       call cubspl(p%v,p%phdos_d(iv,:,:),p%nv,0,0)
    end do

  end subroutine phase_phdos

  !> Write information to the standard output about the average
  !> polynomial fit to the static data.
  subroutine phase_punch_pfit(p)
    use fit, only: fit_pshift
    use evfunc, only: fv1, fv2, fv3, fv4
    use tools, only: error
    use param, only: mline_fmt, au2gpa, uout, warning, format_string, format_string_header
    type(phase), intent(inout) :: p

    integer :: i
    character*(mline_fmt) :: fm
    real*8 :: vk, bk, ek, gk, f1, f2, f3, f4, b1, b2
    integer :: ierr
    real*8 :: prop(5), prop2(5)

    real*8, parameter :: bfrac_warn = 0.05d0

    if (p%pfit%nfit <= 0) return

    write (uout,'("# Composition of the average polynomial and equilibrium static properties: ")')
    fm = format_string_header( &
       (/1,ifmt_order,ifmt_bp,ifmt_v,ifmt_e,ifmt_b,ifmt_bp,ifmt_bpp/),&
       (/1,1,6,9,5,6,2,10/))
    write (uout,fm) "#","n","weight","V(bohr^3)","E(Ha)","B(GPa)","Bp","Bpp(GPa-1)"
    fm = format_string((/ifmt_order,ifmt_bp,ifmt_v,ifmt_e,ifmt_b,ifmt_bp,ifmt_bpp/),1)

    prop = 0d0
    prop2 = 0d0
    do i = 1, p%pfit%nfit
       call fit_pshift(p%pfit%mode(i),p%v,0d0,p%pfit%npar(i),p%pfit%apar(:,i),vk,bk,ek,gk,ierr)
       f1 = fv1(p%pfit%mode(i),vk,p%pfit%npar(i),p%pfit%apar(:,i))
       f2 = fv2(p%pfit%mode(i),vk,p%pfit%npar(i),p%pfit%apar(:,i))
       f3 = fv3(p%pfit%mode(i),vk,p%pfit%npar(i),p%pfit%apar(:,i))
       f4 = fv4(p%pfit%mode(i),vk,p%pfit%npar(i),p%pfit%apar(:,i))
       b1 = -(1+vk*f3/f2)
       b2 = ((f3+vk*f4)/f2**2 - vk*f3**2/f2**3) / au2gpa
       write (uout,fm) p%pfit%npar(i)-2, p%pfit%wei(i), vk, gk, bk, b1, b2

       p%pfit%veq(i) = vk
       p%pfit%beq(i) = bk

       prop = prop + (/vk, gk, bk, b1, b2/) * p%pfit%wei(i)
       prop2 = prop2 + (/vk, gk, bk, b1, b2/)**2 * p%pfit%wei(i)
    end do
    prop2 = sqrt(max(prop2 - prop*prop,0d0))
    fm = format_string((/15,ifmt_v,ifmt_e,ifmt_b,ifmt_bp,ifmt_bpp/),1)

    call fit_pshift(p%fit_mode,p%v,0d0,p%npol,p%cpol,vk,bk,ek,gk,ierr)
    f1 = fv1(p%fit_mode,vk,p%npol,p%cpol)
    f2 = fv2(p%fit_mode,vk,p%npol,p%cpol)
    f3 = fv3(p%fit_mode,vk,p%npol,p%cpol)
    f4 = fv4(p%fit_mode,vk,p%npol,p%cpol)
    b1 = -(1+vk*f3/f2)
    b2 = ((f3+vk*f4)/f2**2 - vk*f3**2/f2**3) / au2gpa
    write (uout,fm) "-average pol.--", vk, gk, bk, b1, b2
    write (uout,fm) "--dir. average-", prop
    write (uout,fm) "---std. dev.---", prop2
    if (prop2(3) > bfrac_warn * prop(3)) then
       call error('phase_punch_pfit','Errorbar in B is too high. Check E(V) data and fit.',warning)
    end if
    write (uout,*)

  end subroutine phase_punch_pfit

  ! Find bracketing volumes for volume v in phase p. If dyn, only the
  ! active volumes are valid. Returns id such that v is between v(id)
  ! and v(id+1), -id if v is within 1d-6 of v(-id), or 0 if not found.
  subroutine vbracket(p,v,id,dyn)
    type(phase), intent(in) :: p
    real*8, intent(in) :: v
    integer, intent(out) :: id
    logical, intent(in) :: dyn

    logical :: found
    integer :: i

    real*8, parameter :: veps = 1d-6

    ! find bracketing volumes
    id = 0
    found = .false.
    do i = 1, p%nv
       if (abs(p%v(i)-v) < veps) then
          if (p%dyn_active(i)) then
             id = -i
          else
             id = 0
          end if
          return
       end if
       if (i == 1) cycle
       if (p%v(i-1)-veps <= v .and. v <= p%v(i)+veps) then
          id = i
          found = .true.
          exit
       end if
    end do
    if (dyn) then
       found = found .and. p%dyn_active(i-1) .and. p%dyn_active(i)
    end if
    if (.not.found) then
       id = 0
    else
       id = i-1
    end if

  end subroutine vbracket

end module varbas
