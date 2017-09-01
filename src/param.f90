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

module param
  implicit none
  public

  ! Math constants
  real*8, parameter :: pi       = 3.14159265358979323846d0
  real*8, parameter :: sqrpi    = 1.77245385090551599275d0
  real*8, parameter :: rad      = pi / 180d0               
  real*8, parameter :: twopi    = 2d0 * pi                 
  real*8, parameter :: halfpi   = pi / 2d0                 
  real*8, parameter :: tosqrpi  = 2d0 / sqrpi              
  real*8, parameter :: pisquare = pi * pi                  
  real*8, parameter :: cte      = 2.71828182845904523536d0
  real*8, parameter :: ctsq2    = 1.41421356237309504880d0
  real*8, parameter :: ctsq3    = 1.73205080756887729352d0
  real*8, parameter :: cteuler  = 0.57721566490153286061d0
  real*8, parameter :: ctgold   = 1.61803398874989484820d0
  real*8, parameter :: zero   = 0d0
  real*8, parameter :: one    = 1d0
  real*8, parameter :: two    = 2d0
  real*8, parameter :: three  = 3d0
  real*8, parameter :: four   = 4d0
  real*8, parameter :: half   = 0.5d0
  real*8, parameter :: third  = one/three
  real*8, parameter :: twothird  = two/three
  real*8, parameter :: fourth = 0.25d0
  real*8, parameter :: undef = -9.72d21

  ! Physical constants and conversion factors.
  !
  ! Most of them involve atomic units:  hartree (Ha), bohr, ...
  ! The physical constant names begin with pc.
  ! The conversion factors are named as: <unit1>2<unit2>
  ! meaning that the factor converts <unit1> to <unit2>.
  !
  real*8, parameter :: pckbau = 3.166815d-6    !Boltzmann ct. [Ha/K] (nist2006,wikipedia 10/10/10)
  real*8, parameter :: pcamu = 1.660538782d-24 !atomic mass unit [g] (nist2006)
  real*8, parameter :: pcme = 9.10938215d-28   !electron mass [g] (nist2006)
  real*8, parameter :: pcna = 6.02214179d23    !Avogadro ct. [1/mol] (nist2006)
  real*8, parameter :: pct0 = 298.15d0         !ambient temperature (K)
  real*8, parameter :: pch = 6.62606896d-34    !Planck ct. [J.s] (nist2006)
  real*8, parameter :: pcc = 2.99792458d10     !light speed [cm/s] (nist2006)

  real*8, parameter :: bohr2cm = 0.52917720859d-8      !bohr -> cm (nist2006)
  real*8, parameter :: bohr2angstrom = 0.52917720859d0 !bohr -> angstrom (nist2006)
  real*8, parameter :: bohr2pm = 0.52917720859d2       !bohr -> pm (nist2006)
  real*8, parameter :: ha2k = 3.1577465d5              !hartree -> K (nist2006)
  real*8, parameter :: ha2ev = 27.21138386d0           !hartree -> eV (nist2006)
  real*8, parameter :: ha2cm_1 = 2.194746313705d5      !hartree -> cm**(-1) (nist2006)
  real*8, parameter :: thz2cm_1 = 33.35641d0           !THz -> cm**(-1) (nist2006)
  real*8, parameter :: ha2thz = ha2cm_1 / thz2cm_1     !hartree -> THz
  real*8, parameter :: ha2kjmol = 2625.4996d0          !hartree -> kJ/mol (nist2006)
  real*8, parameter :: au2gpa = 29421.0108037190       !at.u.(pres) --> GPa (nist2006)
  real*8, parameter :: amu2au = pcamu/pcme           !amu --> at. units

  ! dimension constants
  integer, parameter :: mline = 2048
  integer, parameter :: mline_fmt = 2048
  integer, parameter :: marg = 10

  ! logical units
  integer, parameter :: stderr = 0 !< standard error lu
  integer, parameter :: stdin = 5 !< standard input lu
  integer, parameter :: stdout = 6 !< standard output lu
  integer :: uin = stdin   !< input lu
  integer :: uout = stdout !< output lu

  ! file access modes
  integer, parameter :: ioread = -2 !< open file for reading
  integer, parameter :: iowrite = -3 !< open file for writing
  integer, parameter :: ioappend = -4 !< open file for appending
  integer, parameter :: iofortran = -5 !< open file, fortran carriage control
  integer, parameter :: ioerror = -1 !< error opening file

  ! input information
  character*(mline) :: fileroot
  character*(mline) :: title

  ! error types
  integer, parameter :: faterr = -1 !< fatal error flag
  integer, parameter :: warning = 1 !< warning flag
  integer, parameter :: noerr = 0   !< info flag
  integer :: nwarns = 0
  integer :: ncomms = 0
  
  ! constants that require initialization
  character*(1) :: null !< null character
  character*(1) :: tab !< tab character
  character*(1) :: newline !< newline character
  character*(2) :: eol !< eol char.
  character*(1) :: eof !< eof char.
  character*(1) :: blank !< blank
  character*(1) :: dquote !< "
  character*(1) :: quote !< '
  character*(1) :: backsl !< \

  ! input units
  integer, parameter :: units_v_bohr3 = 1
  integer, parameter :: units_v_ang3 = 2
  integer, parameter :: units_e_ha = 1
  integer, parameter :: units_e_ev = 2
  integer, parameter :: units_e_ry = 3
  integer, parameter :: units_p_au = 1
  integer, parameter :: units_p_gpa = 2
  integer, parameter :: units_f_ha = 1
  integer, parameter :: units_f_cm1 = 2
  integer, parameter :: units_f_thz = 3

  ! colors
  integer, parameter :: mcols = 13
  character*7 :: gplt_rgb(mcols)
  integer :: gplt_sym(mcols)

  ! formats
  integer, parameter :: afmts = 20, ifmts = 20
  integer, parameter :: ifmt_p = afmts+1
  integer, parameter :: ifmt_v = afmts+2
  integer, parameter :: ifmt_x = afmts+3
  integer, parameter :: ifmt_g = afmts+4
  integer, parameter :: ifmt_b = afmts+5
  integer, parameter :: ifmt_bp = afmts+6
  integer, parameter :: ifmt_bpp = afmts+7
  integer, parameter :: ifmt_e = afmts+8
  integer, parameter :: ifmt_bm = afmts+9
  integer, parameter :: ifmt_t = afmts+10
  integer, parameter :: ifmt_cv = afmts+11
  integer, parameter :: ifmt_s = afmts+12
  integer, parameter :: ifmt_dpdt = afmts+13
  integer, parameter :: ifmt_alpha = afmts+14
  integer, parameter :: ifmt_interp = afmts+15
  integer, parameter :: ifmt_eprec = afmts+16
  integer, parameter :: ifmt_aic = afmts+17
  integer, parameter :: ifmt_order = afmts+18
  integer, parameter :: ifmt_ef = afmts+19
  integer, parameter :: ifmt_nef = afmts+20
  character*(12) :: fmt(afmts+ifmts)
  integer :: ifmtlen(afmts+ifmts)

contains

  subroutine param_init()

    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    ! named constants
    null = char(0)
    tab  = char(9)
    eof     = char(04)
    newline = char(10)
    blank   = char(32)
    eol  = newline//null
    dquote = char(34)
    quote  = char(39)
    backsl = char(92)

    ! title
    title = ""

    ! colors and point types
    gplt_rgb = (/&
       "#000000",&
       "#0000FF",&
       "#008B00",&
       "#FF0000",&
       "#FF00FF",&
       "#643700",&
       "#787878",&
       "#00FFFF",&
       "#FFB45A",&
       "#7D26CD",&
       "#CD9B9B",&
       "#CD6D0C",&
       "#00B98B"/)
    gplt_sym = (/&
       4,&
       6,&
       8,&
       10,&
       12,&
       14,&
       3,&
       5,&
       7,&
       9,&
       11,&
       13,&
       1/)

    ! formats
    do i = 1,afmts
       write (fmt(i),'("A",I2.2)') i
       ifmtlen(i) = i
    end do
    fmt(ifmt_p) = "F10.4"
    fmt(ifmt_v) = "F10.4"
    fmt(ifmt_x) = "F10.7"
    fmt(ifmt_g) = "1p,E18.10,0p"
    fmt(ifmt_b) = "F10.4"
    fmt(ifmt_bp) = "F12.7"
    fmt(ifmt_bpp) = "1p,E12.4,0p"
    fmt(ifmt_e) = "1p,E14.6,0p"
    fmt(ifmt_bm) = "1p,E12.4,0p"
    fmt(ifmt_t) = "F10.2"
    fmt(ifmt_cv) = "1p,E14.6,0p"
    fmt(ifmt_s) = "1p,E14.6,0p"
    fmt(ifmt_dpdt) = "1p,E14.6,0p"
    fmt(ifmt_alpha) = "1p,E14.6,0p"
    fmt(ifmt_interp) = "1p,E17.9,0p"
    fmt(ifmt_eprec) = "1p,E17.9,0p"
    fmt(ifmt_aic) = "1p,E11.3,0p"
    fmt(ifmt_order) = "I2"
    fmt(ifmt_ef) = "F12.7"
    fmt(ifmt_nef) = "F14.10"
    ifmtlen(ifmt_p) = 10
    ifmtlen(ifmt_v) = 10
    ifmtlen(ifmt_x) = 10
    ifmtlen(ifmt_g) = 18
    ifmtlen(ifmt_b) = 10
    ifmtlen(ifmt_bp) = 12
    ifmtlen(ifmt_bpp) = 12
    ifmtlen(ifmt_e) = 14
    ifmtlen(ifmt_bm) = 12
    ifmtlen(ifmt_t) = 10
    ifmtlen(ifmt_cv) = 14
    ifmtlen(ifmt_s) = 14
    ifmtlen(ifmt_dpdt) = 14
    ifmtlen(ifmt_alpha) = 14
    ifmtlen(ifmt_interp) = 17
    ifmtlen(ifmt_eprec) = 17
    ifmtlen(ifmt_aic) = 11
    ifmtlen(ifmt_order) = 2
    ifmtlen(ifmt_ef) = 12
    ifmtlen(ifmt_nef) = 14

    ! random seed
    call random_seed(size = n)
    allocate(seed(n))
    call system_clock(count=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)
    deallocate(seed)
    
  end subroutine param_init

  subroutine header()

    write (uout,'("                                                        ")')
    write (uout,'("                 _ _     _         ____                 ")')
    write (uout,'("            __ _(_) |__ | |__  ___|___ \                ")')
    write (uout,'("           / _` | | ''_ \| ''_ \/ __| __) |               ")')
    write (uout,'("          | (_| | | |_) | |_) \__ \/ __/                ")')
    write (uout,'("           \__, |_|_.__/|_.__/|___/_____|               ")')
    write (uout,'("           |___/                                        ")')
    write (uout,'("      GIBBS2: (p,t) thermodynamics of solids.           ")') 
    write (uout,'("                                                        ")')
    write (uout,'("by A. Otero-de-la-Roza, V. Lua~na and D. Abbasi. "/)') 
    write (uout,'("If you find this software useful, please cite:")')
    write (uout,'("  * Comput. Phys. Commun. 182 (2011) 1708--1720 ")')
    write (uout,'("These articles describe the strain polynomial fits:")')
    write (uout,'("  * Comput. Phys. Commun. 182 (2011) 2232--2248.")')
    write (uout,'("  * Comput. Theor. Chem. doi:10.1016/j.comptc.2011.03.050 ")')
    write (uout,'("Empirical energy corrections: ")')
    write (uout,'("  * Phys. Rev. B 84 (2011) 024109 "/)')
    write (uout,'("Dedicated to the memory of Miguel Alvarez Blanco (1969--2010)."/)')

  end subroutine header

  function format_string(ifmt,pad)

    character*(mline_fmt) format_string
    integer, intent(in) :: ifmt(:)
    integer, intent(in) :: pad
    
    integer :: i

    if (pad > 0) then
       write (format_string,'("(",I2.2,"X,")') pad
    else
       format_string = "("
    end if
    do i = 1, size(ifmt)-1
       format_string = trim(format_string) // trim(fmt(ifmt(i))) // ",1x,"
    end do
    format_string = trim(format_string) // fmt(ifmt(size(ifmt))) // ")"

  end function format_string

  function format_string_header(ifmt,iout)

    character*(mline_fmt) format_string_header
    integer, intent(in) :: ifmt(:)
    integer, intent(in) :: iout(:)

    ! space -> n/d left and (d-n)/d right
    ! n/d must be irreducible
    integer, parameter :: n = 2
    integer, parameter :: d = 3

    integer :: i
    integer :: dif, ipad1, ipad2, ishift
    character*(mline_fmt) :: aux

    if (size(ifmt) /= size(iout)) then
       write (uout,'("ifmt /= iout")')
       stop 1
    end if
    format_string_header = "("
    ipad1 = 0
    ipad2 = 0
    dif = 0
    do i = 1, size(ifmt)
       dif = min(dif,0) + ifmtlen(ifmt(i)) - iout(i) 
       ipad1 = ipad2 + max(n*dif/d,0)
       ipad2 = max((d-n)*dif/d + max(min(mod(dif,d),1),0) ,0) + 1
       if (i /= size(ifmt)) then
          if (ipad1 == 0) then
             write (aux,'("A",I2.2,",")') iout(i)
          else
             write (aux,'(I2.2,"X,A",I2.2,",")') ipad1, iout(i)
          end if
       else
          write (aux,'(I2.2,"X,A",I2.2,")")') ipad1, iout(i)
       end if
       format_string_header = trim(format_string_header) // trim(adjustl(aux))
    end do

  endfunction format_string_header

end module param
