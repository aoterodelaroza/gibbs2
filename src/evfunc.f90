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

module evfunc
  use tools
  use param
  implicit none
  private

  public :: evfunc_init
  public :: punch_params
  public :: fcn_minpack, fcn_minpack1
  public :: v2str, str2v
  public :: fv0, fv1, fv2, fv3, fv4
  interface v2str
     module procedure v2strs
     module procedure v2strv
  end interface
  interface str2v
     module procedure str2vs
     module procedure str2vv
  end interface
  interface derivstrain
     module procedure derivstrains
     module procedure derivstrainv
  end interface
  interface fv0
     module procedure fv0s
     module procedure fv0v
  end interface
  interface fv1
     module procedure fv1s
     module procedure fv1v
  end interface
  interface fv2
     module procedure fv2s
     module procedure fv2v
  end interface
  interface fv3
     module procedure fv3s
     module procedure fv3v
  end interface
  interface fv4
     module procedure fv4s
     module procedure fv4v
  end interface

  ! info for minpack
  real*8, allocatable, public :: evfunc_xm(:), evfunc_ym(:)
  integer, public :: evfunc_minpack_mode, evfunc_reg_mode
  logical, public :: evfunc_domask
  integer, public :: evfunc_mask(10), evfunc_npar
  real*8, public :: evfunc_fixval(10)

  ! fit modes
  ! Meaning of the EOS parameters in apar and npar
  ! * polygibbs: npar represents a npar-1 degree polynomial, where:
  !      apar(0) -- apar(n-1), coeficients a_0 .. a_{n-1}
  !      apar(n) is a scaling constant for x
  !      E(V) = apar(0) + apar(1) * (V/apar(n)) + ...
  !             ... + apar(n-1) * (V/apar(n))^(n-1)
  !
  ! * bm2 and pt2: npar = 2
  !   apar(1) = E_0
  !   apar(2) = V_0
  !   apar(3) = B_0
  !
  ! * bm3, pt3, murn, vinet and AP2: the above (npar = 3) plus
  !   apar(4) = B_0'
  ! 
  ! * bm4, pt4: the above (npar = 4) plus
  !   apar(5) = B_0''
  ! 
  ! * pt5: the above (npar = 5) plus
  !   apar(6) = B_0'''
  !
  ! * antons: npar = 4
  !   apar(1) = E_infty
  !   apar(2) = V_0
  !   apar(3) = B_0
  !   apar(4) = B_0'
  !

  integer, parameter, public :: mfit = 12
  integer, parameter, public :: fit_polygibbs = 1
  integer, parameter, public :: fit_bm2 = 2
  integer, parameter, public :: fit_bm3 = 3
  integer, parameter, public :: fit_bm4 = 4
  integer, parameter, public :: fit_pt2 = 5
  integer, parameter, public :: fit_pt3 = 6
  integer, parameter, public :: fit_pt4 = 7
  integer, parameter, public :: fit_pt5 = 8
  integer, parameter, public :: fit_murn = 9
  integer, parameter, public :: fit_antons = 10
  integer, parameter, public :: fit_vinet = 11
  integer, parameter, public :: fit_ap2 = 12
  integer, parameter, public :: fit_strain = 13
  integer, parameter, public :: fit_strain_bm = 1
  integer, parameter, public :: fit_strain_pt = 2
  integer, parameter, public :: fit_strain_lagr = 3
  integer, parameter, public :: fit_strain_inf = 4
  integer, parameter, public :: fit_strain_x1 = 5
  integer, parameter, public :: fit_strain_x3 = 6
  integer, parameter, public :: fit_strain_xinv3 = 7
  integer, parameter, public :: fit_strain_v = 8
  character*8, parameter, public :: fit_name(mfit) = (/&
     "Wei.pol.",&
     "BirchM2 ",&
     "BirchM3 ",&
     "BirchM4 ",&
     "Poirier2",&
     "Poirier3",&
     "Poirier4",&
     "Poirier5",&
     "Murnaghn",&
     "Anton S.",&
     "Vinet   ",&
     "AP2     "/)
  integer, parameter, public :: fit_order(mfit) = &
     (/0,3,4,5,3,4,5,6,4,4,4,4/)

  ! regression models
  integer, parameter, public :: mreg = 2
  integer, parameter, public :: reg_lsq = 1
  integer, parameter, public :: reg_lad = 2
  character*8, parameter, public :: reg_name(mreg) = (/&
     " L. sqr.",&
     "   LAD  "/)

  ! polynomial fit modes
  integer, public :: pfit_mode
  integer, parameter, public :: pfit_gauss = 1
  integer, parameter, public :: pfit_slatec = 2

  ! polynomial weight modes
  integer, public :: pweigh_mode
  integer, parameter, public :: pweigh_gibbs1 = 1
  integer, parameter, public :: pweigh_gibbs2 = 2
  integer, parameter, public :: pweigh_slatec = 3

  ! number of electrons for holzapfel's ap2
  integer, public :: nelectrons

contains

  subroutine evfunc_init()
    
    nelectrons = 0

  end subroutine evfunc_init

  subroutine punch_params(lu, fmode, npar, apar)
    
    integer, intent(in) :: lu, fmode, npar
    real*8, intent(in) :: apar(0:npar)

    integer :: i

    select case(fmode)
    case(fit_polygibbs)
       write (lu,'("# Polynomial EOS, result of averaged fit: ")') 
       write (lu,'("# Degree : ",I2)') npar-2
       write (lu,'("#      p(x) = a_0 + a_1 * (V/V_scal) + ... + a_n * (V/V_scal)^n")') 
       do i = 0, npar-2
          write (lu,'("# a_",I2.2," = ",1p,E20.12,0p)') i, apar(i)
       end do
       write (lu,'("# V_scal (bohr^3) = ",1p,E20.12,0p)') apar(npar-1)
       write (lu,'("# p_scal (GPa) = ",F20.12)') apar(npar) * au2gpa
    case(fit_bm2)
       write (lu,'("# Birch-Murnaghan, 2nd order EOS parameters: ")') 
       write (lu,'("# E_0 (Ha) = ",1p,E17.9,0p)') apar(1)
       write (lu,'("# V_0 (bohr^3) = ",F14.8)') apar(2)
       write (lu,'("# B_0 (GPa) = ",F12.6)') apar(3) * au2gpa
       write (lu,'("# p_scal (GPa) = ",F20.12)') apar(4) * au2gpa
    case(fit_bm3)
       write (lu,'("# Birch-Murnaghan, 3rd order EOS parameters: ")') 
       write (lu,'("# E_0 (Ha) = ",1p,E17.9,0p)') apar(1)
       write (lu,'("# V_0 (bohr^3) = ",F14.8)') apar(2)
       write (lu,'("# B_0 (GPa) = ",F12.6)') apar(3) * au2gpa
       write (lu,'("# B_0'' = ",F12.7)') apar(4)
       write (lu,'("# p_scal (GPa) = ",F20.12)') apar(5) * au2gpa
    case(fit_bm4)
       write (lu,'("# Birch-Murnaghan, 4th order EOS parameters: ")') 
       write (lu,'("# E_0 (Ha) = ",1p,E17.9,0p)') apar(1)
       write (lu,'("# V_0 (bohr^3) = ",F14.8)') apar(2)
       write (lu,'("# B_0 (GPa) = ",F12.6)') apar(3) * au2gpa
       write (lu,'("# B_0'' = ",F12.7)') apar(4)
       write (lu,'("# B_0'''' = ",1p,E15.7,0p)') apar(5) / au2gpa
       write (lu,'("# p_scal (GPa) = ",F20.12)') apar(6) * au2gpa
    case(fit_pt2)
       write (lu,'("# Poirier-Tarantola, 2nd order EOS parameters: ")')
       write (lu,'("# E_0 (Ha) = ",1p,E17.9,0p)') apar(1)
       write (lu,'("# V_0 (bohr^3) = ",F14.8)') apar(2)
       write (lu,'("# B_0 (GPa) = ",F12.6)') apar(3) * au2gpa
       write (lu,'("# p_scal (GPa) = ",F20.12)') apar(4) * au2gpa
    case(fit_pt3)
       write (lu,'("# Poirier-Tarantola, 3rd order EOS parameters: ")')
       write (lu,'("# E_0 (Ha) = ",1p,E17.9,0p)') apar(1)
       write (lu,'("# V_0 (bohr^3) = ",F14.8)') apar(2)
       write (lu,'("# B_0 (GPa) = ",F12.6)') apar(3) * au2gpa
       write (lu,'("# B_0'' = ",F12.7)') apar(4)
       write (lu,'("# p_scal (GPa) = ",F20.12)') apar(5) * au2gpa
    case(fit_pt4)
       write (lu,'("# Poirier-Tarantola, 4th order EOS parameters: ")') 
       write (lu,'("# E_0 (Ha) = ",1p,E17.9,0p)') apar(1)
       write (lu,'("# V_0 (bohr^3) = ",F14.8)') apar(2)
       write (lu,'("# B_0 (GPa) = ",F12.6)') apar(3) * au2gpa
       write (lu,'("# B_0'' = ",F12.7)') apar(4)
       write (lu,'("# B_0'''' = ",1p,E15.7,0p)') apar(5) / au2gpa
       write (lu,'("# p_scal (GPa) = ",F20.12)') apar(6) * au2gpa
    case(fit_pt5)
       write (lu,'("# Poirier-Tarantola, 5th order EOS parameters: ")') 
       write (lu,'("# E_0 (Ha) = ",1p,E17.9,0p)') apar(1)
       write (lu,'("# V_0 (bohr^3) = ",F14.8)') apar(2)
       write (lu,'("# B_0 (GPa) = ",F12.6)') apar(3) * au2gpa
       write (lu,'("# B_0'' = ",F12.7)') apar(4)
       write (lu,'("# B_0'''' = ",1p,E15.7,0p)') apar(5) / au2gpa
       write (lu,'("# B_0'''''' = ",1p,E15.7,0p)') apar(6) / au2gpa**2
       write (lu,'("# p_scal (GPa) = ",F20.12)') apar(7) * au2gpa
    case(fit_murn)
       write (lu,'("# Murnaghan, 3th order EOS parameters: ")') 
       write (lu,'("# E_0 (Ha) = ",1p,E17.9,0p)') apar(1)
       write (lu,'("# V_0 (bohr^3) = ",F14.8)') apar(2)
       write (lu,'("# B_0 (GPa) = ",F12.6)') apar(3) * au2gpa
       write (lu,'("# B_0'' = ",F12.7)') apar(4)
       write (lu,'("# p_scal (GPa) = ",F20.12)') apar(5) * au2gpa
    case(fit_antons)
       write (lu,'("# Anton-Schmidt, 3th order EOS parameters: ")') 
       write (lu,'("# E_infty (Ha) = ",1p,E17.9,0p)') apar(1)
       write (lu,'("# V_0 (bohr^3) = ",F14.8)') apar(2)
       write (lu,'("# B_0 (GPa) = ",F12.6)') apar(3) * au2gpa
       write (lu,'("# B_0'' = ",F12.7)') apar(4)
       write (lu,'("# p_scal (GPa) = ",F20.12)') apar(5) * au2gpa
    case(fit_vinet)
       write (lu,'("# Vinet, 3th order EOS parameters: ")') 
       write (lu,'("# E_0 (Ha) = ",1p,E17.9,0p)') apar(1)
       write (lu,'("# V_0 (bohr^3) = ",F14.8)') apar(2)
       write (lu,'("# B_0 (GPa) = ",F12.6)') apar(3) * au2gpa
       write (lu,'("# B_0'' = ",F12.7)') apar(4)
       write (lu,'("# p_scal (GPa) = ",F20.12)') apar(5) * au2gpa
    case(fit_ap2)
       write (lu,'("# Holzapfel AP2, 3th order EOS parameters: ")') 
       write (lu,'("# E_0 (Ha) = ",1p,E17.9,0p)') apar(1)
       write (lu,'("# V_0 (bohr^3) = ",F14.8)') apar(2)
       write (lu,'("# B_0 (GPa) = ",F12.6)') apar(3) * au2gpa
       write (lu,'("# B_0'' = ",F12.7)') apar(4)
       write (lu,'("# p_scal (GPa) = ",F20.12)') apar(5) * au2gpa
    case default
       write (lu,'("# Polynomial fit to strain: ")') 
       write (lu,'("# Degree : ",I2)') npar-2
       write (lu,'("#      p(x) = a_0 + a_1 * f(V) + ... + a_n * f(V)^n")') 
       do i = 0, npar-2
          write (lu,'("# a_",I2.2," = ",1p,E20.12,0p)') i, apar(i)
       end do
       write (lu,'("# V_scal (bohr^3) = ",1p,E20.12,0p)') apar(npar-1)
       write (lu,'("# p_scal (GPa) = ",F20.12)') apar(npar) * au2gpa
    end select
    write (lu,*)

  end subroutine punch_params

  function v2strs(strain,v,v0)

    integer, intent(in) :: strain
    real*8, intent(in) :: v, v0
    real*8 :: v2strs

    select case(strain)
    case (fit_strain_bm)
       v2strs = 0.5d0*((v/v0)**(-twothird)-1)
    case (fit_strain_pt)
       v2strs = log(v/v0)/3d0
    case (fit_strain_lagr)
       v2strs = 0.5d0*((v/v0)**(twothird)-1)
    case (fit_strain_inf)
       v2strs = -(v/v0)**(-third)+1
    case (fit_strain_x1)
       v2strs = v/v0
    case (fit_strain_x3)
       v2strs = (v/v0)**third
    case (fit_strain_xinv3)
       v2strs = (v/v0)**(-third)
    case (fit_strain_v)
       v2strs = v
    end select

  end function v2strs

  function v2strv(strain,v,v0)

    integer, intent(in) :: strain
    real*8, intent(in) :: v(:), v0
    real*8 :: v2strv(size(v))

    select case(strain)
    case (fit_strain_bm)
       v2strv = 0.5d0*((v/v0)**(-twothird)-1)
    case (fit_strain_pt)
       v2strv = log(v/v0)/3d0
    case (fit_strain_lagr)
       v2strv = 0.5d0*((v/v0)**(twothird)-1)
    case (fit_strain_inf)
       v2strv = -(v/v0)**(-third)+1
    case (fit_strain_x1)
       v2strv = v/v0
    case (fit_strain_x3)
       v2strv = (v/v0)**third
    case (fit_strain_xinv3)
       v2strv = (v/v0)**(-third)
    case (fit_strain_v)
       v2strv = v
    end select

  end function v2strv

  function str2vs(strain,f,v0)

    integer, intent(in) :: strain
    real*8, intent(in) :: f, v0
    real*8 :: str2vs

    select case(strain)
    case (fit_strain_bm)
       str2vs = (f*2+1)**(-3d0/2d0) * v0
    case (fit_strain_pt)
       str2vs = exp(f*3) * v0
    case (fit_strain_lagr)
       str2vs = (f*2+1)**(3d0/2d0) * v0
    case (fit_strain_inf)
       str2vs = (-f+1)**(-3) * v0
    case (fit_strain_x1)
       str2vs = f * v0
    case (fit_strain_x3)
       str2vs = f**3 * v0
    case (fit_strain_xinv3)
       str2vs = f**(-3) * v0
    case (fit_strain_v)
       str2vs = f
    end select

  end function str2vs

  function str2vv(strain,f,v0)

    integer, intent(in) :: strain
    real*8, intent(in) :: f(:), v0
    real*8 :: str2vv(size(f))

    select case(strain)
    case (fit_strain_bm)
       str2vv = (f*2+1)**(-3d0/2d0) * v0
    case (fit_strain_pt)
       str2vv = exp(f*3) * v0
    case (fit_strain_lagr)
       str2vv = (f*2+1)**(3d0/2d0) * v0
    case (fit_strain_inf)
       str2vv = (-f+1)**(-3) * v0
    case (fit_strain_x1)
       str2vv = f * v0
    case (fit_strain_x3)
       str2vv = f**3 * v0
    case (fit_strain_xinv3)
       str2vv = f**(-3) * v0
    case (fit_strain_v)
       str2vv = f
    end select

  end function str2vv

  subroutine derivstrains(strain,f,v,v0,deg,f1v,f2v,f3v,f4v)
    
    integer, intent(in) :: strain
    real*8, intent(in) :: f, v, v0
    integer, intent(in) :: deg
    real*8, intent(out) :: f1v, f2v, f3v, f4v

    real*8 :: f2, ss

    f1v = 0d0
    f2v = 0d0
    f3v = 0d0
    f4v = 0d0
    select case(strain)
    case (fit_strain_bm)
       f2 = 2 * f + 1
       ss = -f2**(3d0/2d0)/(3d0*v0)
       f1v = -f2**(5d0/2d0) / (3d0 * v0)
       if (deg == 1) return
       f2v = 5d0 * f1v * ss
       if (deg == 2) return
       f3v = 8d0 * f2v * ss
       if (deg == 3) return
       f4v = 11d0 * f3v * ss
    case (fit_strain_pt)
       ss = -1/v
       f1v = 1/(3*v)
       if (deg == 1) return
       f2v = f1v * ss
       if (deg == 2) return
       f3v = f2v * ss * 2
       if (deg == 3) return
       f4v = f3v * ss * 3
    case (fit_strain_lagr)
       ss = -1/(3*v)
       f1v = (f + f + 1)**(-half) / (3*v0)
       if (deg == 1) return
       f2v = ss * f1v
       if (deg == 2) return
       f3v = ss * f2v * 4
       if (deg == 3) return
       f4v = ss * f3v * 7
    case (fit_strain_inf)
       ss = -1/(3*v)
       f1v = (-f + 1)**4 / (3*v0)
       if (deg == 1) return
       f2v = f1v * ss * 4
       if (deg == 2) return
       f3v = f2v * ss * 7
       if (deg == 3) return
       f4v = f3v * ss * 10
    case (fit_strain_x1)
       f1v = 1/v0
       if (deg == 1) return
       f2v = 0d0
       if (deg == 2) return
       f3v = 0d0
       if (deg == 3) return
       f4v = 0d0
    case (fit_strain_x3)
       ss = -1/(3*v)
       f1v = v**(-twothird)/(3d0*v0**third)
       if (deg == 1) return
       f2v = f1v * ss * 2
       if (deg == 2) return
       f3v = f2v * ss * 5
       if (deg == 3) return
       f4v = f3v * ss * 8
    case (fit_strain_xinv3)
       ss = -1/(3*v)
       f1v = -v**(-4d0/3d0) * (v0**third/3)
       if (deg == 1) return
       f2v = f1v * ss * 4
       if (deg == 2) return
       f3v = f2v * ss * 7
       if (deg == 3) return
       f4v = f3v * ss * 10
    case (fit_strain_v)
       f1v = 1
       if (deg == 1) return
       f2v = 0
       if (deg == 2) return
       f3v = 0
       if (deg == 3) return
       f4v = 0
    end select

  end subroutine derivstrains

  subroutine derivstrainv(strain,f,v,v0,deg,f1v,f2v,f3v,f4v)
    
    integer, intent(in) :: strain
    real*8, intent(in) :: f(:), v(:), v0
    integer, intent(in) :: deg
    real*8, dimension(size(f)), intent(out) :: f1v, f2v, f3v, f4v

    real*8 :: f2(size(f)), ss(size(f))

    f1v = 0d0
    f2v = 0d0
    f3v = 0d0
    f4v = 0d0
    select case(strain)
    case (fit_strain_bm)
       f2 = 2 * f + 1
       ss = -f2**(3d0/2d0) / (3d0*v0)
       f1v = -f2**(5d0/2d0) / (3d0 * v0)
       if (deg == 1) return
       f2v = 5d0 * f1v * ss
       if (deg == 2) return
       f3v = 8d0 * f2v * ss
       if (deg == 3) return
       f4v = 11d0 * f3v * ss
    case (fit_strain_pt)
       ss = -1/v
       f1v = 1/(3*v)
       if (deg == 1) return
       f2v = f1v * ss
       if (deg == 2) return
       f3v = f2v * ss * 2
       if (deg == 3) return
       f4v = f3v * ss * 3
    case (fit_strain_lagr)
       ss = -1/(3*v)
       f1v = (f + f + 1)**(-half) / (3*v0)
       if (deg == 1) return
       f2v = ss * f1v
       if (deg == 2) return
       f3v = ss * f2v * 4
       if (deg == 3) return
       f4v = ss * f3v * 7
    case (fit_strain_inf)
       ss = -1/(3*v)
       f1v = (-f + 1)**4 / (3*v0)
       if (deg == 1) return
       f2v = f1v * ss * 4
       if (deg == 2) return
       f3v = f2v * ss * 7
       if (deg == 3) return
       f4v = f3v * ss * 10
    case (fit_strain_x1)
       f1v = 1/v0
       if (deg == 1) return
       f2v = 0d0
       if (deg == 2) return
       f3v = 0d0
       if (deg == 3) return
       f4v = 0d0
    case (fit_strain_x3)
       ss = -1/(3*v)
       f1v = v**(-twothird)/(3d0*v0**third)
       if (deg == 1) return
       f2v = f1v * ss * 2
       if (deg == 2) return
       f3v = f2v * ss * 5
       if (deg == 3) return
       f4v = f3v * ss * 8
    case (fit_strain_xinv3)
       ss = -1/(3*v)
       f1v = -v**(-4d0/3d0) * (v0**third/3)
       if (deg == 1) return
       f2v = f1v * ss * 4
       if (deg == 2) return
       f3v = f2v * ss * 7
       if (deg == 3) return
       f4v = f3v * ss * 10
    case (fit_strain_v)
       f1v = 1
       if (deg == 1) return
       f2v = 0
       if (deg == 2) return
       f3v = 0
       if (deg == 3) return
       f4v = 0
    end select

  end subroutine derivstrainv

  function fv0s (mode, x, npar, apar) result(y)
    
    integer, intent(in) :: mode
    real*8 :: y
    integer, intent(in) :: npar
    real*8, intent(in) :: x, apar(0:npar)

    integer :: i

    select case(mode)
    case(fit_polygibbs,100:)
       y = polin0s(mode,x,npar-1,apar)
    case(fit_bm2)
       y = bm2s(x,apar(1:3),0)
    case(fit_bm3)
       y = bm3s(x,apar(1:4),0)
    case(fit_bm4)
       y = bm4s(x,apar(1:5),0)
    case(fit_pt2)
       y = pt2s(x,apar(1:3),0)
    case(fit_pt3)
       y = pt3s(x,apar(1:4),0)
    case(fit_pt4)
       y = pt4s(x,apar(1:5),0)
    case(fit_pt5)
       y = pt5s(x,apar(1:6),0)
    case(fit_murn)
       y = murns(x,apar(1:4),0)
    case(fit_antons)
       y = antonss(x,apar(1:4),0)
    case(fit_vinet)
       y = vinets(x,apar(1:4),0)
    case(fit_ap2)
       y = ap2s(x,apar(1:4),0)
    case default
       call error('gibbs2','fv0s: unknown E(V) equation',faterr)
    end select

    y = y - apar(npar) * x

  end function fv0s

  function fv0v (mode, x, npar, apar) result(y)
    
    integer, intent(in) :: mode
    integer, intent(in) :: npar
    real*8, intent(in) :: x(:), apar(0:npar)

    real*8 :: y(size(x))
    integer :: strain

    select case(mode)
    case(fit_polygibbs,100:)
       y = polin0v(mode,x,npar-1,apar)
    case(fit_bm2)
       y = bm2v(x,apar(1:3),0)
    case(fit_bm3)
       y = bm3v(x,apar(1:4),0)
    case(fit_bm4)
       y = bm4v(x,apar(1:5),0)
    case(fit_pt2)
       y = pt2v(x,apar(1:3),0)
    case(fit_pt3)
       y = pt3v(x,apar(1:4),0)
    case(fit_pt4)
       y = pt4v(x,apar(1:5),0)
    case(fit_pt5)
       y = pt5v(x,apar(1:6),0)
    case(fit_murn)
       y = murnv(x,apar(1:4),0)
    case(fit_antons)
       y = antonsv(x,apar(1:4),0)
    case(fit_vinet)
       y = vinetv(x,apar(1:4),0)
    case(fit_ap2)
       y = ap2v(x,apar(1:4),0)
    case default
       call error('gibbs2','fv0s: unknown E(V) equation',faterr)
    end select

    y = y - apar(npar) * x

  end function fv0v

  function fv1s (mode, x, npar, apar) result(y)
    
    integer, intent(in) :: mode
    real*8 :: y
    integer, intent(in) :: npar
    real*8, intent(in) :: x, apar(0:npar)

    integer :: i

    select case(mode)
    case(fit_polygibbs,100:)
       y = polin1s(mode,x,npar-1,apar)
    case(fit_bm2)
       y = bm2s(x,apar(1:3),1)
    case(fit_bm3)
       y = bm3s(x,apar(1:4),1)
    case(fit_bm4)
       y = bm4s(x,apar(1:5),1)
    case(fit_pt2)
       y = pt2s(x,apar(1:3),1)
    case(fit_pt3)
       y = pt3s(x,apar(1:4),1)
    case(fit_pt4)
       y = pt4s(x,apar(1:5),1)
    case(fit_pt5)
       y = pt5s(x,apar(1:6),1)
    case(fit_murn)
       y = murns(x,apar(1:4),1)
    case(fit_antons)
       y = antonss(x,apar(1:4),1)
    case(fit_vinet)
       y = vinets(x,apar(1:4),1)
    case(fit_ap2)
       y = ap2s(x,apar(1:4),1)
    case default
       call error('gibbs2','fv1s: unknown E(V) equation',faterr)
    end select

    y = y - apar(npar)

  end function fv1s

  function fv1v (mode, x, npar, apar) result(y)
    
    integer, intent(in) :: mode
    integer, intent(in) :: npar
    real*8, intent(in) :: x(:), apar(0:npar)
    real*8 :: y(size(x))

    integer :: strain

    select case(mode)
    case(fit_polygibbs,100:)
       y = polin1v(mode,x,npar-1,apar)
    case(fit_bm2)
       y = bm2v(x,apar(1:3),1)
    case(fit_bm3)
       y = bm3v(x,apar(1:4),1)
    case(fit_bm4)
       y = bm4v(x,apar(1:5),1)
    case(fit_pt2)
       y = pt2v(x,apar(1:3),1)
    case(fit_pt3)
       y = pt3v(x,apar(1:4),1)
    case(fit_pt4)
       y = pt4v(x,apar(1:5),1)
    case(fit_pt5)
       y = pt5v(x,apar(1:6),1)
    case(fit_murn)
       y = murnv(x,apar(1:4),1)
    case(fit_antons)
       y = antonsv(x,apar(1:4),1)
    case(fit_vinet)
       y = vinetv(x,apar(1:4),1)
    case(fit_ap2)
       y = ap2v(x,apar(1:4),1)
    case default
       call error('gibbs2','fv1v: unknown E(V) equation',faterr)
    end select

    y = y - apar(npar)

  end function fv1v

  function fv2s (mode, x, npar, apar) result(y)
    
    integer, intent(in) :: mode
    real*8 :: y
    integer, intent(in) :: npar
    real*8, intent(in) :: x, apar(0:npar)

    select case(mode)
    case(fit_polygibbs,100:)
       y = polin2s(mode,x,npar-1,apar)
    case(fit_bm2)
       y = bm2s(x,apar(1:3),2)
    case(fit_bm3)
       y = bm3s(x,apar(1:4),2)
    case(fit_bm4)
       y = bm4s(x,apar(1:5),2)
    case(fit_pt2)
       y = pt2s(x,apar(1:3),2)
    case(fit_pt3)
       y = pt3s(x,apar(1:4),2)
    case(fit_pt4)
       y = pt4s(x,apar(1:5),2)
    case(fit_pt5)
       y = pt5s(x,apar(1:6),2)
    case(fit_murn)
       y = murns(x,apar(1:4),2)
    case(fit_antons)
       y = antonss(x,apar(1:4),2)
    case(fit_vinet)
       y = vinets(x,apar(1:4),2)
    case(fit_ap2)
       y = ap2s(x,apar(1:4),2)
    case default
       call error('gibbs2','fv2s: unknown E(V) equation',faterr)
    end select
  end function fv2s

  function fv2v (mode, x, npar, apar) result(y)
    
    integer, intent(in) :: mode
    integer, intent(in) :: npar
    real*8, intent(in) :: x(:), apar(0:npar)
    real*8 :: y(size(x))

    select case(mode)
    case(fit_polygibbs,100:)
       y = polin2v(mode,x,npar-1,apar)
    case(fit_bm2)
       y = bm2v(x,apar(1:3),2)
    case(fit_bm3)
       y = bm3v(x,apar(1:4),2)
    case(fit_bm4)
       y = bm4v(x,apar(1:5),2)
    case(fit_pt2)
       y = pt2v(x,apar(1:3),2)
    case(fit_pt3)
       y = pt3v(x,apar(1:4),2)
    case(fit_pt4)
       y = pt4v(x,apar(1:5),2)
    case(fit_pt5)
       y = pt5v(x,apar(1:6),2)
    case(fit_murn)
       y = murnv(x,apar(1:4),2)
    case(fit_antons)
       y = antonsv(x,apar(1:4),2)
    case(fit_vinet)
       y = vinetv(x,apar(1:4),2)
    case(fit_ap2)
       y = ap2v(x,apar(1:4),2)
    case default
       call error('gibbs2','fv2v: unknown E(V) equation',faterr)
    end select
  end function fv2v

  function fv3s (mode, x, npar, apar) result(y)
    
    integer, intent(in) :: mode
    real*8 :: y
    integer, intent(in) :: npar
    real*8, intent(in) :: x, apar(0:npar)

    select case(mode)
    case(fit_polygibbs,100:)
       y = polin3s(mode,x,npar-1,apar)
    case(fit_bm2)
       y = bm2s(x,apar(1:3),3)
    case(fit_bm3)
       y = bm3s(x,apar(1:4),3)
    case(fit_bm4)
       y = bm4s(x,apar(1:5),3)
    case(fit_pt2)
       y = pt2s(x,apar(1:3),3)
    case(fit_pt3)
       y = pt3s(x,apar(1:4),3)
    case(fit_pt4)
       y = pt4s(x,apar(1:5),3)
    case(fit_pt5)
       y = pt5s(x,apar(1:6),3)
    case(fit_murn)
       y = murns(x,apar(1:4),3)
    case(fit_antons)
       y = antonss(x,apar(1:4),3)
    case(fit_vinet)
       y = vinets(x,apar(1:4),3)
    case(fit_ap2)
       y = ap2s(x,apar(1:4),3)
    case default
       call error('gibbs2','fv3s: unknown E(V) equation',faterr)
    end select
  end function fv3s

  function fv3v (mode, x, npar, apar) result(y)
    
    integer, intent(in) :: mode
    integer, intent(in) :: npar
    real*8, intent(in) :: x(:), apar(0:npar)
    real*8 :: y(size(x))

    select case(mode)
    case(fit_polygibbs,100:)
       y = polin3v(mode,x,npar-1,apar)
    case(fit_bm2)
       y = bm2v(x,apar(1:3),3)
    case(fit_bm3)
       y = bm3v(x,apar(1:4),3)
    case(fit_bm4)
       y = bm4v(x,apar(1:5),3)
    case(fit_pt2)
       y = pt2v(x,apar(1:3),3)
    case(fit_pt3)
       y = pt3v(x,apar(1:4),3)
    case(fit_pt4)
       y = pt4v(x,apar(1:5),3)
    case(fit_pt5)
       y = pt5v(x,apar(1:6),3)
    case(fit_murn)
       y = murnv(x,apar(1:4),3)
    case(fit_antons)
       y = antonsv(x,apar(1:4),3)
    case(fit_vinet)
       y = vinetv(x,apar(1:4),3)
    case(fit_ap2)
       y = ap2v(x,apar(1:4),3)
    case default
       call error('gibbs2','fv3v: unknown E(V) equation',faterr)
    end select
  end function fv3v

  function fv4s (mode, x, npar, apar) result(y)
    
    integer, intent(in) :: mode
    real*8 :: y
    integer, intent(in) :: npar
    real*8, intent(in) :: x, apar(0:npar)

    select case(mode)
    case(fit_polygibbs,100:)
       y = polin4s(mode,x,npar-1,apar)
    case(fit_bm2)
       y = bm2s(x,apar(1:3),4)
    case(fit_bm3)
       y = bm3s(x,apar(1:4),4)
    case(fit_bm4)
       y = bm4s(x,apar(1:5),4)
    case(fit_pt2)
       y = pt2s(x,apar(1:3),4)
    case(fit_pt3)
       y = pt3s(x,apar(1:4),4)
    case(fit_pt4)
       y = pt4s(x,apar(1:5),4)
    case(fit_pt5)
       y = pt5s(x,apar(1:6),4)
    case(fit_murn)
       y = murns(x,apar(1:4),4)
    case(fit_antons)
       y = antonss(x,apar(1:4),4)
    case(fit_vinet)
       y = vinets(x,apar(1:4),4)
    case(fit_ap2)
       y = ap2s(x,apar(1:4),4)
    case default
       call error('gibbs2','fv4s: unknown E(V) equation',faterr)
    end select
  end function fv4s

  function fv4v (mode, x, npar, apar) result(y)
    
    integer, intent(in) :: mode
    integer, intent(in) :: npar
    real*8, intent(in) :: x(:), apar(0:npar)
    real*8 :: y(size(x))

    select case(mode)
    case(fit_polygibbs,100:)
       y = polin4v(mode,x,npar-1,apar)
    case(fit_bm2)
       y = bm2v(x,apar(1:3),4)
    case(fit_bm3)
       y = bm3v(x,apar(1:4),4)
    case(fit_bm4)
       y = bm4v(x,apar(1:5),4)
    case(fit_pt2)
       y = pt2v(x,apar(1:3),4)
    case(fit_pt3)
       y = pt3v(x,apar(1:4),4)
    case(fit_pt4)
       y = pt4v(x,apar(1:5),4)
    case(fit_pt5)
       y = pt5v(x,apar(1:6),4)
    case(fit_murn)
       y = murnv(x,apar(1:4),4)
    case(fit_antons)
       y = antonsv(x,apar(1:4),4)
    case(fit_vinet)
       y = vinetv(x,apar(1:4),4)
    case(fit_ap2)
       y = ap2v(x,apar(1:4),4)
    case default
       call error('gibbs2','fv4v: unknown E(V) equation',faterr)
    end select
  end function fv4v

  function polin0s (mode, x, npar, apar) result(y)
    !.polin0 - Horner's evaluation of a polynomial.
    ! The polinomial is given by:
    !   y(x) = SUM(i=0,npar) apar(i) * x**i
    ! It is assumed that npar>=0.

    integer, intent(in) :: mode
    real*8 :: y
    integer, intent(in) :: npar
    real*8, intent(in) :: x, apar(0:npar)

    real*8 :: xp
    integer :: i, strain

    if (mode == fit_polygibbs) then
       strain = fit_strain_x1
    else
       strain = (mode - fit_strain * 10000) / 100
    end if
    xp = v2strs(strain,x,apar(npar))

    y = apar(npar-1)
    do i = npar-1, 1, -1
       y = y * xp + apar(i-1)
    enddo

  end function polin0s

  function polin0v (mode, x, npar, apar) result(y)
    !.polin0 - Horner's evaluation of a polynomial.
    ! The polinomial is given by:
    !   y(x) = SUM(i=0,npar) apar(i) * x**i
    ! It is assumed that npar>=0.

    integer, intent(in) :: mode, npar
    real*8, intent(in) :: x(:), apar(0:npar)
    real*8 :: y(size(x))

    real*8 :: xp(size(x))
    integer :: i, strain

    if (mode == fit_polygibbs) then
       strain = fit_strain_x1
    else
       strain = (mode - fit_strain * 10000) / 100
    end if
    xp = v2strv(strain,x,apar(npar))

    y = apar(npar-1)
    do i = npar-1, 1, -1
       y = y * xp + apar(i-1)
    enddo

  end function polin0v

  function polin1s (mode, x, npar, apar) result(y)
    !.polin1 - Horner's evaluation of the first derivative of a
    ! polynomial.
    ! It is assumed that npar>=1.

    integer, intent(in) :: mode, npar
    real*8, intent(in) :: x, apar(0:npar)
    real*8 :: y

    integer :: i, strain
    real*8 :: xp
    real*8 :: f1v, f2v, f3v, f4v

    if (mode == fit_polygibbs) then
       strain = fit_strain_x1
    else
       strain = (mode - fit_strain * 10000) / 100
    end if
    xp = v2strs(strain,x,apar(npar))

    y = apar(npar-1) * (npar-1)
    do i = npar-1, 2, -1
       y = y * xp + apar(i-1) * (i-1)
    enddo
    call derivstrains(strain,xp,x,apar(npar),1,f1v,f2v,f3v,f4v)
    y = y * f1v

  end function polin1s

  function polin1v (mode, x, npar, apar) result(y)
    !.polin1 - Horner's evaluation of the first derivative of a
    ! polynomial.
    ! It is assumed that npar>=1.

    integer, intent(in) :: mode, npar
    real*8, intent(in) :: x(:), apar(0:npar)
    real*8 :: y(size(x))

    integer :: i, strain
    real*8 :: xp(size(x))
    real*8, dimension(size(x)) :: f1v, f2v, f3v, f4v

    if (mode == fit_polygibbs) then
       strain = fit_strain_x1
    else
       strain = (mode - fit_strain * 10000) / 100
    end if
    xp = v2strv(strain,x,apar(npar))

    y = apar(npar-1) * (npar-1)
    do i = npar-1, 2, -1
       y = y * xp + apar(i-1) * (i-1)
    enddo
    call derivstrainv(strain,xp,x,apar(npar),1,f1v,f2v,f3v,f4v)
    y = y * f1v

  end function polin1v

  function polin2s (mode, x, npar, apar) result(y)
    !.polin2 - Horner's evaluation of the second derivative of a
    ! polynomial.
    ! It is assumed that npar>=2.
    
    integer, intent(in) :: npar, mode
    real*8, intent(in) :: x, apar(0:npar)
    real*8 :: y

    integer :: i, strain
    real*8 :: xp, y1
    real*8 :: f1v, f2v, f3v, f4v

    if (mode == fit_polygibbs) then
       strain = fit_strain_x1
    else
       strain = (mode - fit_strain * 10000) / 100
    end if
    xp = v2strs(strain,x,apar(npar))

    y1 = apar(npar-1) * (npar-1)
    y = apar(npar-1) * (npar-1) * (npar-2)
    do i = npar-1, 2, -1
       y1 = y1 * xp + apar(i-1) * (i-1)
       if (i < 3) cycle
       y = y * xp + apar(i-1) * (i-1) * (i-2)
    enddo
    call derivstrains(strain,xp,x,apar(npar),2,f1v,f2v,f3v,f4v)
    y = y*f1v**2 + y1*f2v

  end function polin2s

  function polin2v (mode, x, npar, apar) result(y)
    !.polin2 - Horner's evaluation of the second derivative of a
    ! polynomial.
    ! It is assumed that npar>=2.
    
    integer, intent(in) :: npar, mode
    real*8, intent(in) :: x(:), apar(0:npar)
    real*8 :: y(size(x))

    integer :: i, strain
    real*8, dimension(size(x)) :: xp, y1, f1v, f2v, f3v, f4v

    if (mode == fit_polygibbs) then
       strain = fit_strain_x1
    else
       strain = (mode - fit_strain * 10000) / 100
    end if
    xp = v2strv(strain,x,apar(npar))

    y1 = apar(npar-1) * (npar-1)
    y = apar(npar-1) * (npar-1) * (npar-2)
    do i = npar-1, 2, -1
       y = y * xp + apar(i-1) * (i-1) * (i-2)
       if (i < 3) cycle
       y = y * xp + apar(i-1) * (i-1) * (i-2)
    enddo
    call derivstrainv(strain,xp,x,apar(npar),2,f1v,f2v,f3v,f4v)
    y = y*f1v**2 + y1*f2v

  end function polin2v

  function polin3s (mode, x, npar, apar) result(y)
    !.polin3 - Horner's evaluation of the third derivative of a
    ! polynomial.
    ! It is assumed that npar>=3.

    integer :: mode, npar
    real*8, intent(in) :: x, apar(0:npar)
    real*8 :: y

    integer :: i, strain
    real*8 :: xp, y1, y2
    real*8 :: f1v, f2v, f3v, f4v

    if (mode == fit_polygibbs) then
       strain = fit_strain_x1
    else
       strain = (mode - fit_strain * 10000) / 100
    end if
    xp = v2strs(strain,x,apar(npar))

    y1 = apar(npar-1) * (npar-1)
    y2 = apar(npar-1) * (npar-1) * (npar-2)
    y = apar(npar-1) * (npar-1) * (npar-2) * (npar-3)
    do i = npar-1, 2, -1
       y1 = y1 * xp + apar(i-1) * (i-1)
       if (i < 3) cycle
       y2 = y2 * xp + apar(i-1) * (i-1) * (i-2)
       if (i < 4) cycle
       y = y * xp + apar(i-1) * (i-1) * (i-2) * (i-3)
    enddo
    call derivstrains(strain,xp,x,apar(npar),3,f1v,f2v,f3v,f4v)
    y = y*f1v**3 + y2*f1v*f2v*3 + y1*f3v

  end function polin3s

  function polin3v (mode, x, npar, apar) result(y)
    !.polin3 - Horner's evaluation of the third derivative of a
    ! polynomial.
    ! It is assumed that npar>=3.

    integer :: mode, npar
    real*8, intent(in) :: x(:), apar(0:npar)
    real*8 :: y(size(x))

    integer :: i, strain
    real*8, dimension(size(x)) :: xp, y1, y2, f1v, f2v, f3v, f4v

    if (mode == fit_polygibbs) then
       strain = fit_strain_x1
    else
       strain = (mode - fit_strain * 10000) / 100
    end if
    xp = v2strv(strain,x,apar(npar))

    y1 = apar(npar-1) * (npar-1)
    y2 = apar(npar-1) * (npar-1) * (npar-2)
    y = apar(npar-1) * (npar-1) * (npar-2) * (npar-3)
    do i = npar-1, 2, -1
       y1 = y1 * xp + apar(i-1) * (i-1)
       if (i < 3) cycle
       y2 = y2 * xp + apar(i-1) * (i-1) * (i-2)
       if (i < 4) cycle
       y = y * xp + apar(i-1) * (i-1) * (i-2) * (i-3)
    enddo
    call derivstrainv(strain,xp,x,apar(npar),3,f1v,f2v,f3v,f4v)
    y = y*f1v**3 + y2*f1v*f2v*3 + y1*f3v

  end function polin3v

  function polin4s (mode, x, npar, apar) result(y)
    !.polin4 - Horner's evaluation of the fourth derivative of a
    ! polynomial.
    ! It is assumed that npar>=4.
    integer, intent(in) :: mode, npar
    real*8, intent(in) :: x, apar(0:npar)
    real*8 :: y

    integer :: i, strain
    real*8 :: xp, y1, y2, y3
    real*8 :: f1v, f2v, f3v, f4v

    if (mode == fit_polygibbs) then
       strain = fit_strain_x1
    else
       strain = (mode - fit_strain * 10000) / 100
    end if
    xp = v2strs(strain,x,apar(npar))

    y1 = apar(npar-1) * (npar-1)
    y2 = apar(npar-1) * (npar-1) * (npar-2)
    y3 = apar(npar-1) * (npar-1) * (npar-2) * (npar-3)
    y = apar(npar-1) * (npar-1) * (npar-2) * (npar-3) * (npar-4)
    do i = npar-1, 2, -1 
       y1 = y1 * xp + apar(i-1) * (i-1)
       if (i < 3) cycle
       y2 = y2 * xp + apar(i-1) * (i-1) * (i-2)
       if (i < 4) cycle
       y3 = y3 * xp + apar(i-1) * (i-1) * (i-2) * (i-3)
       if (i < 5) cycle
       y = y * xp + apar(i-1) * (i-1) * (i-2) * (i-3) * (i-4)
    enddo
    call derivstrains(strain,xp,x,apar(npar),4,f1v,f2v,f3v,f4v)
    y = y*f1v**4 + y3*f1v**2*f2v*6 + y2*(f1v*f3v*4+f2v**2*3) + y1*f4v

  end function polin4s

  function polin4v (mode, x, npar, apar) result(y)
    !.polin4 - Horner's evaluation of the fourth derivative of a
    ! polynomial.
    ! It is assumed that npar>=4.
    integer, intent(in) :: mode, npar
    real*8, intent(in) :: x(:), apar(0:npar)
    real*8 :: y(size(x))

    integer :: i, strain
    real*8, dimension(size(x)) :: xp, y1, y2, y3, f1v, f2v, f3v, f4v
    
    if (mode == fit_polygibbs) then
       strain = fit_strain_x1
    else
       strain = (mode - fit_strain * 10000) / 100
    end if
    xp = v2strv(strain,x,apar(npar))

    y1 = apar(npar-1) * (npar-1)
    y2 = apar(npar-1) * (npar-1) * (npar-2)
    y3 = apar(npar-1) * (npar-1) * (npar-2) * (npar-3)
    y = apar(npar-1) * (npar-1) * (npar-2) * (npar-3) * (npar-4)
    do i = npar-1, 2, -1
       y1 = y1 * xp + apar(i-1) * (i-1)
       if (i < 3) cycle
       y2 = y2 * xp + apar(i-1) * (i-1) * (i-2)
       if (i < 4) cycle
       y3 = y3 * xp + apar(i-1) * (i-1) * (i-2) * (i-3)
       if (i < 5) cycle
       y = y * xp + apar(i-1) * (i-1) * (i-2) * (i-3) * (i-4)
    enddo
    call derivstrainv(strain,xp,x,apar(npar),4,f1v,f2v,f3v,f4v)
    y = y*f1v**4 + y3*f1v**2*f2v*6 + y2*(f1v*f3v*4+f2v**2*3) + y1*f4v

  end function polin4v

  subroutine fcn_minpack(m,n,x0,fvec,iflag)
    ! wrapper for minpack fit

    integer, intent(in) :: m,n,iflag
    real*8, intent(in) :: x0(n)
    real*8, intent(out) :: fvec(m)

    real*8 :: x(evfunc_npar)
    integer :: i

    if (evfunc_domask) then
       do i = 1, evfunc_npar
          if (evfunc_mask(i) <= n) then
             x(i) = x0(evfunc_mask(i))
          else
             x(i) = evfunc_fixval(i)
          end if
       end do
    else
       x = x0
    end if

    select case(evfunc_minpack_mode)
    case(fit_bm2)
       fvec = bm2v(evfunc_xm(1:m),x,0) - evfunc_ym(1:m)
    case(fit_bm3)
       fvec = bm3v(evfunc_xm(1:m),x,0) - evfunc_ym(1:m)
    case(fit_bm4)
       fvec = bm4v(evfunc_xm(1:m),x,0) - evfunc_ym(1:m)
    case(fit_pt2)
       fvec = pt2v(evfunc_xm(1:m),x,0) - evfunc_ym(1:m)
    case(fit_pt3)
       fvec = pt3v(evfunc_xm(1:m),x,0) - evfunc_ym(1:m)
    case(fit_pt4)
       fvec = pt4v(evfunc_xm(1:m),x,0) - evfunc_ym(1:m)
    case(fit_pt5)
       fvec = pt5v(evfunc_xm(1:m),x,0) - evfunc_ym(1:m)
    case(fit_murn)
       fvec = murnv(evfunc_xm(1:m),x,0) - evfunc_ym(1:m)
    case(fit_antons)
       fvec = antonsv(evfunc_xm(1:m),x,0) - evfunc_ym(1:m)
    case(fit_vinet)
       fvec = vinetv(evfunc_xm(1:m),x,0) - evfunc_ym(1:m)
    case(fit_ap2)
       fvec = ap2v(evfunc_xm(1:m),x,0) - evfunc_ym(1:m)
    case default
       call error('fcn_minpack','unknown E(V) expression',faterr)
    end select

    if (evfunc_reg_mode == reg_lad) fvec = sqrt(abs(fvec))

  end subroutine fcn_minpack

  subroutine fcn_minpack1(m,n,x,fvec,iflag)
    ! wrapper for minpack fit

    ! x(1) E_0 ; x(2) V_0 ; x(3) B_0 ; x(4) B_0'
    integer, intent(in) :: m,n,iflag
    real*8, intent(in) :: x(n)
    real*8, intent(out) :: fvec(m)

    select case(evfunc_minpack_mode)
    case(fit_bm2)
       fvec = -bm2v(evfunc_xm(1:m),x,1) - evfunc_ym(1:m)
    case(fit_bm3)
       fvec = -bm3v(evfunc_xm(1:m),x,1) - evfunc_ym(1:m)
    case(fit_bm4)
       fvec = -bm4v(evfunc_xm(1:m),x,1) - evfunc_ym(1:m)
    case(fit_pt2)
       fvec = -pt2v(evfunc_xm(1:m),x,1) - evfunc_ym(1:m)
    case(fit_pt3)
       fvec = -pt3v(evfunc_xm(1:m),x,1) - evfunc_ym(1:m)
    case(fit_pt4)
       fvec = -pt4v(evfunc_xm(1:m),x,1) - evfunc_ym(1:m)
    case(fit_pt5)
       fvec = -pt5v(evfunc_xm(1:m),x,1) - evfunc_ym(1:m)
    case(fit_murn)
       fvec = -murnv(evfunc_xm(1:m),x,1) - evfunc_ym(1:m)
    case(fit_antons)
       fvec = -antonsv(evfunc_xm(1:m),x,1) - evfunc_ym(1:m)
    case(fit_vinet)
       fvec = -vinetv(evfunc_xm(1:m),x,1) - evfunc_ym(1:m)
    case(fit_ap2)
       fvec = -ap2v(evfunc_xm(1:m),x,1) - evfunc_ym(1:m)
    case default
       call error('fcn_minpack','unknown E(V) expression',faterr)
    end select

    if (evfunc_reg_mode == reg_lad) fvec = sqrt(abs(fvec))

  end subroutine fcn_minpack1

  function bm2s(V,x,ider)
    ! birch-murnaghan 2, scalar, 0-4th derivative

    real*8, intent(in) :: V, x(3)
    integer, intent(in) :: ider
    real*8 :: bm2s
    
    real*8 :: c, d, f, ff

    f = 0.5d0*((x(2) / V)**(twothird) - 1d0)
    c = 4.5d0*x(3)*x(2)
    if (ider == 0) then
       bm2s = f * f * c + x(1)
    else if (ider == 1) then
       ff = (2*f+1)**(5d0/2d0)
       bm2s = -f*ff* 2d0*c/(3d0*x(2))
    else if (ider == 2) then
       ff = (2*f+1)**4
       bm2s = ff * (7*f+1) * (x(3)/x(2))
    else if (ider == 3) then
       ff = (2*f+1)**(11d0/2d0)
       bm2s = -ff * (14*f+3) * (2.5d0*x(3)/x(2)**2)
    else if (ider == 4) then
       ff = (2*f+1)**7
       bm2s = ff * (182*f+47) * (5*x(3)/9d0/x(2)**3)
    end if

  end function bm2s

  function bm2v(V,x,ider)
    ! birch-murnaghan 2, vector, 0-4th derivative

    real*8, intent(in) :: V(:), x(3)
    integer, intent(in) :: ider
    real*8 :: bm2v(size(V))
    
    real*8 :: c, f(size(V)), ff(size(V))

    f = 0.5d0*((x(2) / V)**(twothird) - 1d0)
    c = 4.5d0*x(3)*x(2)
    if (ider == 0) then
       bm2v = f * f * c + x(1)
    else if (ider == 1) then
       ff = (2*f+1)**(5d0/2d0)
       bm2v = -f*ff* 2d0*c/(3d0*x(2))
    else if (ider == 2) then
       ff = (2*f+1)**4
       bm2v = ff * (7*f+1) * (x(3)/x(2))
    else if (ider == 3) then
       ff = (2*f+1)**(11d0/2d0)
       bm2v = -ff * (14*f+3) * (2.5d0*x(3)/x(2)**2)
    else if (ider == 4) then
       ff = (2*f+1)**7
       bm2v = ff * (182*f+47) * (5*x(3)/9d0/x(2)**3)
    end if

  end function bm2v

  function bm3s(V,x,ider)
    ! birch-murnaghan 3, scalar, 0-4th derivative

    real*8, intent(in) :: V, x(4)
    integer, intent(in) :: ider
    real*8 :: bm3s
    
    real*8 :: c, d, f, ff

    f = 0.5d0*((x(2) / V)**(twothird) - 1d0)
    c = 4.5d0*x(3)*x(2)
    d = c * (x(4)-4d0)
    if (ider == 0) then
       bm3s = f * f * (d * f + c) + x(1)
    else if (ider == 1) then
       ff = (2*f+1)**(5d0/2d0)
       bm3s = -f*((3*d)*f+2*c)*ff / (3*x(2)) 
    else if (ider == 2) then
       ff = (2*f+1)**4
       bm3s = (2*c + f*((14*c+6*d) + 27*f))*ff / (9*x(2)**2)
    else if (ider == 3) then
       ff = (2*f+1)**(11d0/2d0)
       bm3s = -((15*c+3*d)+f*((70*c+57*d)+f*(162*d)))*ff*(2/(27*x(2)**3))
    else if (ider == 4) then
       ff = (2*f+1)**7
       bm3s = ((47*c+18*d)+f*((182*c+213*d)+f*(486*d)))*ff*(10/(81*x(2)**4))
    end if

  end function bm3s

  function bm3v(V,x,ider)
    ! birch-murnaghan 3, vector, 0-4th derivative

    real*8, intent(in) :: V(:), x(4)
    integer, intent(in) :: ider
    real*8 :: bm3v(size(V))
    
    real*8 :: c, d, f(size(V)), ff(size(V))

    f = 0.5d0*((x(2) / V)**(twothird) - 1d0)
    c = 4.5d0*x(3)*x(2)
    d = c * (x(4)-4d0)
    if (ider == 0) then
       bm3v = f * f * (d * f + c) + x(1)
    else if (ider == 1) then
       ff = (2*f+1)**(5d0/2d0)
       bm3v = -f*((3*d)*f+2*c)*ff / (3*x(2))
    else if (ider == 2) then
       ff = (2*f+1)**4
       bm3v = (2*c + f*((14*c+6*d) + 27*f))*ff / (9*x(2)**2)
    else if (ider == 3) then
       ff = (2*f+1)**(11d0/2d0)
       bm3v = -((15*c+3*d)+f*((70*c+57*d)+f*(162*d)))*ff*(2/(27*x(2)**3))
    else if (ider == 4) then
       ff = (2*f+1)**7
       bm3v = ((47*c+18*d)+f*((182*c+213*d)+f*(486*d)))*ff*(10/(81*x(2)**4))
    end if

  end function bm3v

  function bm4s(V,x,ider)
    ! birch-murnaghan 4, scalar, 0-4th derivative

    real*8, intent(in) :: V, x(5)
    integer, intent(in) :: ider
    real*8 :: bm4s
    
    real*8 :: c, d, ee, f, ff

    f = 0.5d0*((x(2) / V)**(twothird) - 1d0)
    c = 4.5d0*x(3)*x(2)
    d = c * (x(4)-4d0)
    ee = (3d0/8d0)*x(3)*x(2)*(9*(x(5)*x(3)+x(4)**2)-63*x(4)+143)

    if (ider == 0) then
       bm4s = f*f*(f*(f*ee+d)+c)+x(1)
    else if (ider == 1) then
       ff = (2*f+1)**(5d0/2d0)
       bm4s = -ff*f*(f*(f*(4*ee)+(3*d))+(2*c))/(3*x(2))
    else if (ider == 2) then
       ff = (2*f+1)**4
       bm4s = ff*(f*(f*(f*(44*ee)+(27*d+12*ee))+(14*c+6*d))+(2*c))/(9*x(2)**2)
    else if (ider == 3) then
       ff = (2*f+1)**(11d0/2d0)
       bm4s = -ff*(f*(f*(f*(308*ee)+(162*d+138*ee))+(70*c+57*d+12*ee))+(15*c+3*d))*(2/(3*x(2))**3)
    else if (ider == 4) then
       ff = (2*f+1)**7
       bm4s = ff*(f*(f*(f*(5236*ee)+(2430*d+2994*ee))+(910*c+165*d+432*ee))+(235*c+90*d+12*ee))*(2/(3*x(2))**4)
    end if

  end function bm4s
  
  function bm4v(V,x,ider)
    ! birch-murnaghan 4, vector, 0-4th derivative

    real*8, intent(in) :: V(:), x(5)
    integer, intent(in) :: ider
    real*8 :: bm4v(size(V))

    real*8 :: c, d, ee, f(size(V)), ff(size(V))

    f = 0.5d0*((x(2) / V)**(twothird) - 1d0)
    c = 4.5d0*x(3)*x(2)
    d = c * (x(4)-4d0)
    ee = (3d0/8d0)*x(3)*x(2)*(9*(x(5)*x(3)+x(4)**2)-63*x(4)+143)

    if (ider == 0) then
       bm4v = f*f*(f*(f*ee+d)+c)+x(1)
    else if (ider == 1) then
       ff = (2*f+1)**(5d0/2d0)
       bm4v = -ff*f*(f*(f*(4*ee)+(3*d))+(2*c))/(3*x(2))
    else if (ider == 2) then
       ff = (2*f+1)**4
       bm4v = ff*(f*(f*(f*(44*ee)+(27*d+12*ee))+(14*c+6*d))+(2*c))/(9*x(2)**2)
    else if (ider == 3) then
       ff = (2*f+1)**(11d0/2d0)
       bm4v = -ff*(f*(f*(f*(308*ee)+(162*d+138*ee))+(70*c+57*d+12*ee))+(15*c+3*d))*(2/(3*x(2))**3)
    else if (ider == 4) then
       ff = (2*f+1)**7
       bm4v = ff*(f*(f*(f*(5236*ee)+(2430*d+2994*ee))+(910*c+165*d+432*ee))+(235*c+90*d+12*ee))*(2/(3*x(2))**4)
    end if

  end function bm4v
  
  function pt2s(V,x,ider)
    ! poirier-tarantola 2, scalar, 0-4th derivative

    real*8, intent(in) :: V, x(3)
    integer, intent(in) :: ider
    real*8 :: pt2s

    real*8 :: f, fpv, c

    f = log(V/x(2))/3d0
    fpv = 1/(3*V)
    c = 4.5d0*x(3)*x(2)

    if (ider == 0) then
       pt2s = f**2*c + x(1)
    else if (ider == 1) then
       pt2s = f*fpv*(2*c)
    else if (ider == 2) then
       pt2s = -fpv**2*(f*3-1)*(2*c)
    else if (ider == 3) then
       pt2s = +fpv**3*(f*2-1)*(18*c)
    else if (ider == 4) then
       pt2s = -fpv**4*(f*18-11)*(18*c)
    end if

  end function pt2s

  function pt2v(V,x,ider)
    ! poirier-tarantola 2, vector, 0-4th derivative

    real*8, intent(in) :: V(:), x(3)
    integer, intent(in) :: ider
    real*8 :: pt2v(size(V))

    real*8 :: f(size(V)), fpv(size(V)), c

    f = log(V/x(2))/3d0
    fpv = 1/(3*V)
    c = 4.5d0*x(3)*x(2)

    if (ider == 0) then
       pt2v = f**2*c + x(1)
    else if (ider == 1) then
       pt2v = f*fpv*(2*c)
    else if (ider == 2) then
       pt2v = -fpv**2*(f*3-1)*(2*c)
    else if (ider == 3) then
       pt2v = +fpv**3*(f*2-1)*(18*c)
    else if (ider == 4) then
       pt2v = -fpv**4*(f*18-11)*(18*c)
    end if
  end function pt2v

  function pt3s(V,x,ider)
    ! poirier-tarantola 3, scalar, 0-4th derivative

    real*8, intent(in) :: V, x(4)
    integer, intent(in) :: ider
    real*8 :: pt3s

    real*8 :: f, fpv, c, d
    real*8 :: c0, c1, c2

    f = log(V/x(2))/3d0
    fpv = 1/(3*V)
    c = 4.5d0*x(3)*x(2)
    d=c*(x(4)+2)

    if (ider == 0) then
       pt3s = f**2*(f*d+c) + x(1)
    else if (ider == 1) then
       pt3s = f*fpv*(f*(3*d)+(2*c))
    else if (ider == 2) then
       c2=9*d
       c1=-6*d+6*c
       c0=-2*c
       pt3s = -fpv**2*(f*(f*c2+c1)+c0)
    else if (ider == 3) then
       c2=9*d
       c1=-9*d+6*c
       c0=d-3*c
       pt3s = +fpv**3*(f*(f*c2+c1)+c0)*18
    else if (ider == 4) then
       c2=27*d
       c1=-33*d+18*c
       c0=6*d-11*c
       pt3s = -fpv**4*(f*18-11)*(18*c)
    end if

  end function pt3s

  function pt3v(V,x,ider)
    ! poirier-tarantola 3, vector, 0-4th derivative

    real*8, intent(in) :: V(:), x(4)
    integer, intent(in) :: ider
    real*8 :: pt3v(size(V))
    real*8 :: f(size(V)), fpv(size(V)), c, d
    real*8 :: c0, c1, c2

    f = log(V/x(2))/3d0
    fpv = 1/(3*V)
    c = 4.5d0*x(3)*x(2)
    d=c*(x(4)+2)
    
    if (ider == 0) then
       pt3v = f**2*(f*d+c) + x(1)
    else if (ider == 1) then
       pt3v = f*fpv*(f*(3*d)+(2*c))
    else if (ider == 2) then
       c2=9*d
       c1=-6*d+6*c
       c0=-2*c
       pt3v = -fpv**2*(f*(f*c2+c1)+c0)
    else if (ider == 3) then
       c2=9*d
       c1=-9*d+6*c
       c0=d-3*c
       pt3v = +fpv**3*(f*(f*c2+c1)+c0)*18
    else if (ider == 4) then
       c2=27*d
       c1=-33*d+18*c
       c0=6*d-11*c
       pt3v = -fpv**4*(f*18-11)*(18*c)
    end if

  end function pt3v

  function pt4s(V,x,ider)
    ! poirier-tarantola 4, scalar, 0-4th derivative

    real*8, intent(in) :: V, x(5)
    integer, intent(in) :: ider
    real*8 :: pt4s

    real*8 :: f, fpv, c, d, ee
    real*8 :: c0, c1, c2, c3

    f = log(V/x(2))/3d0
    fpv = 1/(3*V)
    c = 4.5d0*x(3)*x(2)
    d=c*(x(4)+2)
    ee=(9*(d*d-c*d+c*c)*x(2)+2*c**3*x(5))/(12*c*x(2))

    if (ider == 0) then
       pt4s = f**2 * (f*(f*ee+d)+c) + x(1)
    else if (ider == 1) then
       pt4s = f * fpv * (f*(f*(4*ee)+(3*d))+(2*c))
    else if (ider == 2) then
       c3=12*ee
       c2=-12*ee+9*d
       c1=-6*d+6*c
       c0=-2*c
       pt4s = -fpv**2 * (f*(f*(f*(c3)+(c2))+(c1))+(c0))
    else if (ider == 3) then
       c3=12*ee
       c2=-18*ee+9*d
       c1=4*ee-9*d+6*c
       c0=d-3*c
       pt4s = +fpv**3 * (f*(f*(f*(c3)+(c2))+(c1))+(c0)) * 18
    else if (ider == 4) then
       c3=108*ee
       c2=-198*ee+81*d
       c1=72*ee-99*d+54*c
       c0=-4*ee+18*d-33*c
       pt4s = -fpv**4 * (f*(f*(f*(c3)+(c2))+(c1))+(c0)) * 6
    end if

  end function pt4s

  function pt4v(V,x,ider)
    ! poirier-tarantola 4, vector, 0-4th derivative

    real*8, intent(in) :: V(:), x(5)
    integer, intent(in) :: ider
    real*8 :: pt4v(size(V))

    real*8 :: f(size(V)), fpv(size(V)), c, d, ee
    real*8 :: c0, c1, c2, c3

    f = log(V/x(2))/3d0
    fpv = 1/(3*V)
    c = 4.5d0*x(3)*x(2)
    d=c*(x(4)+2)
    ee=(9*(d*d-c*d+c*c)*x(2)+2*c**3*x(5))/(12*c*x(2))

    if (ider == 0) then
       pt4v = f**2 * (f*(f*ee+d)+c) + x(1)
    else if (ider == 1) then
       pt4v = f * fpv * (f*(f*(4*ee)+(3*d))+(2*c))
    else if (ider == 2) then
       c3=12*ee
       c2=-12*ee+9*d
       c1=-6*d+6*c
       c0=-2*c
       pt4v = -fpv**2 * (f*(f*(f*(c3)+(c2))+(c1))+(c0))
    else if (ider == 3) then
       c3=12*ee
       c2=-18*ee+9*d
       c1=4*ee-9*d+6*c
       c0=d-3*c
       pt4v = +fpv**3 * (f*(f*(f*(c3)+(c2))+(c1))+(c0)) * 18
    else if (ider == 4) then
       c3=108*ee
       c2=-198*ee+81*d
       c1=72*ee-99*d+54*c
       c0=-4*ee+18*d-33*c
       pt4v = -fpv**4 * (f*(f*(f*(c3)+(c2))+(c1))+(c0)) * 6
    end if

  end function pt4v

  function pt5s(V,x,ider)
    ! poirier-tarantola 5, scalar, 0-4th derivative

    real*8, intent(in) :: V, x(6)
    integer, intent(in) :: ider
    real*8 :: pt5s

    real*8 :: f, fpv, c, d, ee, g
    real*8 :: c0, c1, c2, c3, c4

    f = log(V/x(2))/3d0
    fpv = 1/(3*V)
    c = 4.5d0*x(3)*x(2)
    d=c*(2-x(4))
    ee=(9*(d*d-c*d+c*c)*x(2)+2*c**3*x(5))/(12*c*x(2))
    g=((432*(d-c)*c*ee-243*d**3+486*(d-c)*c*d+324*c**3)*x(2)*x(2)-4*c**5*x(6))/(180*(c*x(2))**2)

    if (ider == 0) then
       pt5s = f**2 * (f*(f*(f*g+ee)+d)+c) + x(1)
    else if (ider == 1) then
       pt5s = f * fpv * (f*(f*(f*(5*g)+(4*ee))+(3*d))+(2*c))
    else if (ider == 2) then
       c4=15*g
       c3=-20*g+12*ee
       c2=-12*ee+9*d
       c1=-6*d+6*c
       c0=-2*c
       pt5s = -fpv**2 * (f*(f*(f*(f*(c4)+(c3))+(c2))+(c1))+(c0))
    else if (ider == 3) then
       c4=15*g
       c3=-30*g+12*ee
       c2=10*g-18*ee+9*d
       c1=4*ee-9*d+6*c
       c0=d-3*c
       pt5s = +fpv**3 * (f*(f*(f*(f*(c4)+(c3))+(c2))+(c1))+(c0)) * 6
    else if (ider == 4) then
       c4=135*g
       c3=-330*g+108*ee
       c2=180*g-198*ee+81*d
       c1=-20*g+72*ee-99*d+54*c
       c0=-4*ee+18*d-33*c
       pt5s = -fpv**4 * (f*(f*(f*(f*(c4)+(c3))+(c2))+(c1))+(c0)) * 6
    end if

  end function pt5s

  function pt5v(V,x,ider)
    ! poirier-tarantola 5, vector, 0-4th derivative

    real*8, intent(in) :: V(:), x(6)
    integer, intent(in) :: ider
    real*8 :: pt5v(size(V))

    real*8 :: f(size(V)), fpv(size(V)), c, d, ee, g
    real*8 :: c0, c1, c2, c3, c4

    f = log(V/x(2))/3d0
    fpv = 1/(3*V)
    c = 4.5d0*x(3)*x(2)
    d=c*(x(4)+2)
    ee=(9*(d*d-c*d+c*c)*x(2)+2*c**3*x(5))/(12*c*x(2))
    g=((432*(d-c)*c*ee-243*d**3+486*(d-c)*c*d+324*c**3)*x(2)*x(2)+4*c**5*x(6))/(180*(c*x(2))**2)

    if (ider == 0) then
       pt5v = f**2 * (f*(f*(f*g+ee)+d)+c) + x(1)
    else if (ider == 1) then
       pt5v = f * fpv * (f*(f*(f*(5*g)+(4*ee))+(3*d))+(2*c))
    else if (ider == 2) then
       c4=15*g
       c3=-20*g+12*ee
       c2=-12*ee+9*d
       c1=-6*d+6*c
       c0=-2*c
       pt5v = -fpv**2 * (f*(f*(f*(f*(c4)+(c3))+(c2))+(c1))+(c0))
    else if (ider == 3) then
       c4=15*g
       c3=-30*g+12*ee
       c2=10*g-18*ee+9*d
       c1=4*ee-9*d+6*c
       c0=d-3*c
       pt5v = +fpv**3 * (f*(f*(f*(f*(c4)+(c3))+(c2))+(c1))+(c0)) * 6
    else if (ider == 4) then
       c4=135*g
       c3=-330*g+108*ee
       c2=180*g-198*ee+81*d
       c1=-20*g+72*ee-99*d+54*c
       c0=-4*ee+18*d-33*c
       pt5v = -fpv**4 * (f*(f*(f*(f*(c4)+(c3))+(c2))+(c1))+(c0)) * 6
    end if

  end function pt5v

  function murns(V,x,ider)
    ! murnaghan (3), scalar, 0-4th derivative

    real*8, intent(in) :: V, x(4)
    integer, intent(in) :: ider
    real*8 :: murns

    if (ider == 0) then
       murns = V*((x(2)/V)**x(4)/(x(4)-1)+1)*(x(3)/x(4)) + x(1) - x(3)*x(2)/(x(4)-1)
    else if (ider == 1) then
       murns = -((x(2)/V)**x(4)-1)*(x(3)/x(4))
    else if (ider == 2) then
       murns = V**(-(x(4)+1)) * (x(3)*x(2)**x(4))
    else if (ider == 3) then
       murns = -V**(-(x(4)+2)) * (x(3)*x(2)**x(4)*(x(4)+1))
    else if (ider == 4) then
       murns = V**(-(x(4)+3)) * (x(3)*x(2)**x(4)*(x(4)+1)*(x(4)+2))
    end if

  end function murns

  function murnv(V,x,ider)
    ! murnaghan (3), vector, 0-4th derivative

    real*8, intent(in) :: V(:), x(4)
    integer, intent(in) :: ider
    real*8 :: murnv(size(V))

    if (ider == 0) then
       murnv = V*((x(2)/V)**x(4)/(x(4)-1)+1)*(x(3)/x(4)) + x(1) - x(3)*x(2)/(x(4)-1)
    else if (ider == 1) then
       murnv = -((x(2)/V)**x(4)-1)*(x(3)/x(4))
    else if (ider == 2) then
       murnv = V**(-(x(4)+1)) * (x(3)*x(2)**x(4))
    else if (ider == 3) then
       murnv = -V**(-(x(4)+2)) * (x(3)*x(2)**x(4)*(x(4)+1))
    else if (ider == 4) then
       murnv = V**(-(x(4)+3)) * (x(3)*x(2)**x(4)*(x(4)+1)*(x(4)+2))
    end if

  end function murnv

  function antonss(V,x,ider)
    ! Anton-Schmidt (3), scalar, 0-4th derivative

    real*8, intent(in) :: V, x(4)
    integer, intent(in) :: ider
    real*8 :: antonss

    real*8 :: xr, lxr, n

    xr = (V / x(2))
    lxr = log(xr)
    n = -x(4)/2d0
    if (ider == 0) then
       antonss = xr**(n+1) *(lxr-(1d0/(n+1)))*(x(3)*x(2)/(n+1)) + x(1)
    else if (ider == 1) then
       antonss = xr**n * lxr *x(3)
    else if (ider == 2) then
       antonss = xr**(n-1) *(lxr*n+1) *(x(3)/x(2))
    else if (ider == 3) then
       antonss = xr**(n-2) *(lxr*(n**2-n)+(2*n-1)) *(x(3)/x(2)**2)
    else if (ider == 4) then
       antonss = xr**(n-3) *(lxr*(n**3-3*n*n+2*n)+(3*n**2-6*n+2)) *(x(3)/x(2)**3)
    end if

  end function antonss

  function antonsv(V,x,ider)
    ! murnaghan (3), vector, 0-4th derivative

    real*8, intent(in) :: V(:), x(4)
    integer, intent(in) :: ider
    real*8 :: antonsv(size(V))

    real*8 :: xr(size(V)), lxr(size(V)), n

    xr = (V / x(2))
    lxr = log(xr)
    n = -x(4)/2d0
    if (ider == 0) then
       antonsv = xr**(n+1) *(lxr-(1d0/(n+1)))*(x(3)*x(2)/(n+1)) + x(1)
    else if (ider == 1) then
       antonsv = xr**n * lxr *x(3)
    else if (ider == 2) then
       antonsv = xr**(n-1) *(lxr*n+1) *(x(3)/x(2))
    else if (ider == 3) then
       antonsv = xr**(n-2) *(lxr*(n**2-n)+(2*n-1)) *(x(3)/x(2)**2)
    else if (ider == 4) then
       antonsv = xr**(n-3) *(lxr*(n**3-3*n*n+2*n)+(3*n**2-6*n+2)) *(x(3)/x(2)**3)
    end if

  end function antonsv

  function vinets(V,x,ider)
    ! Vinet (3), scalar, 0-4th derivative

    real*8, intent(in) :: V, x(4)
    integer, intent(in) :: ider
    real*8 :: vinets

    real*8 :: xr, xr1, xr2, expo, b0p1

    xr = (V/x(2))**third
    xr1 = xr - 1
    xr2 = xr1 - 1
    b0p1 = x(4) - 1
    expo = exp(-xr1 * (1.5d0*b0p1))
    if (ider == 0) then
       vinets = (-expo * (xr1*(3*b0p1)+2) + 2) * (2*x(3)*x(2)/b0p1**2) + x(1)
    else if (ider == 1) then
       vinets = (expo * xr1 / xr**2) * (3*x(3))
    else if (ider == 2) then
       vinets = -expo*(xr*xr1*(3*b0p1) + xr2*2)/xr**5*(x(3)/(2*x(2)))
    else if (ider == 3) then
       vinets = expo*(xr*(xr*(xr*(9*b0p1**2)-9d0*b0p1**2+24d0*b0p1)+52d0-36d0*x(4))-40d0)/&
          xr**8*(x(3)/(12*x(2)**2))
    else if (ider == 4) then
       vinets = -expo*(xr*(xr*(xr*(xr*27d0*b0p1**3&
          - x(4)*(x(4)*(x(4)*27d0-243d0)+405d0)+ 189d0)-x(4)*(x(4)*216d0-768d0)-552d0)&
          - 624d0*x(4)+848d0)-640)/xr**11*(x(3)/(72*x(2)**3))
    end if

  end function vinets

  function vinetv(V,x,ider)
    ! Vinet (3), vector, 0-4th derivative

    real*8, intent(in) :: V(:), x(4)
    integer, intent(in) :: ider
    real*8 :: vinetv(size(V))

    real*8 :: xr(size(V)), xr1(size(V)), xr2(size(V)), expo(size(V)), b0p1

    xr = (V/x(2))**third
    xr1 = xr - 1
    xr2 = xr1 - 1
    b0p1 = x(4) - 1
    expo = exp(-xr1 * (1.5d0*b0p1))
    if (ider == 0) then
       vinetv = (-expo * (xr1*(3*b0p1)+2) + 2) * (2*x(3)*x(2)/b0p1**2) + x(1)
    else if (ider == 1) then
       vinetv = (expo * xr1 / xr**2) * (3*x(3))
    else if (ider == 2) then
       vinetv = -expo*(xr*xr1*(3*b0p1) + xr2*2)/xr**5*(x(3)/(2*x(2)))
    else if (ider == 3) then
       vinetv = expo*(xr*(xr*(xr*(9*b0p1**2)-9d0*b0p1**2+24d0*b0p1)+52d0-36d0*x(4))-40d0)/&
          xr**8*(x(3)/(12*x(2)**2))
    else if (ider == 4) then
       vinetv = -expo*(xr*(xr*(xr*(xr*27d0*b0p1**3&
          - x(4)*(x(4)*(x(4)*27d0-243d0)+405d0)+ 189d0)-x(4)*(x(4)*216d0-768d0)-552d0)&
          - 624d0*x(4)+848d0)-640)/xr**11*(x(3)/(72*x(2)**3))
    end if

  end function vinetv

  function ap2s(V,x,ider)
    ! AP2 (3), scalar, 0-4th derivative

    real*8, intent(in) :: V, x(4)
    integer, intent(in) :: ider
    real*8 :: ap2s

    real*8 :: eta, pfg, c0, c2, z, ec0, g0, g1, g2, g3, g4, g5, g6
    real*8 :: gg, gexp, gfac, dum1, dum2
    real*8 :: gfz, gfc
    integer :: ierr, ierr2

    interface
       subroutine dgam(A, X, ACC, G, GSTAR, IFLG, IFLGST)
         double precision a, x, g, gstar
         real acc
         integer iflg, iflgst
       end subroutine dgam
    end interface

    if (nelectrons == 0) &
       call error('ap2s','NELECTRONS is necessary for AP2 fits',faterr)

    eta = (V/x(2))**third
    pfg = (3*pisquare)**twothird/5d0 * (nelectrons/x(2))**(5d0/3d0)
    c0 = -log(3*x(3)/pfg)
    c2 = 1.5d0 * (x(4)-3) - c0

    if (ider == 0) then
       z = eta*c0
       ec0 = exp(c0)
       gfz = gammai(-2d0,z)
       gfc = gammai(-2d0,c0)
       g2 = (gfz-gfc) * (c0**2 * ec0)
       gfz = gammai(-1d0,z)
       gfc = gammai(-1d0,c0)
       g1 = (gfz-gfc) * (c0 * (c2-1) * ec0)
       gfz = gammai(0d0,z)
       gfc = gammai(0d0,c0)
       g0 = (gfz-gfc) * (-2 * c2 * ec0)
       gg = (exp(-z+c0)-1)*(c2/c0)
       ap2s = (g2+g1+g0+gg) * (9*x(2)*x(3)) + x(1)
    else if (ider == 1) then
       gexp = -exp((-eta+1)*c0)*(3*x(3)) / eta**5
       ap2s = gexp * (-eta+1) * ((-eta+1)*eta*c2+1)
    else if (ider == 2) then
       gfac = -1d0/(eta**3 * 3 * x(2))
       gexp = -exp((-eta+1)*c0)*(3*x(3)) / eta**5 * gfac
       g1 = 4*c2+c0-4
       g2 = (c0-6)*c2-c0
       g3 = 2*c2*(1-c0)
       g4 = c0*c2
       ap2s = gexp * (eta*(eta*(eta*(eta*g4+g3)+g2)+g1)+5)
    else if (ider == 3) then
       gfac = -1d0/(eta**3 * 3 * x(2))
       gexp = -exp((-eta+1)*c0)*(3*x(3)) / eta**5 * gfac**2
       g1 = 28*(c2-1)+12*c0
       g2 = 10*c0*(c2-1)-36*c2+c0**2
       g3 = c0**2*(c2-1)-16*c0*c2+10*c2
       g4 = 2*c0*c2*(3-c0)
       g5 = c0**2*c2
       ap2s = gexp * (eta*(eta*(eta*(eta*(eta*g5+g4)+g3)+g2)+g1)+40)
    else if (ider == 4) then
       gfac = -1d0/(eta**3 * 3 * x(2))
       gexp = -exp((-eta+1)*c0)*(3*x(3)) / eta**5 * gfac**3
       g1 = 280*(c2-1)+160*c0
       g2 = 118*c0*(c2-1)-324*c2+21*c0**2
       g3 = 18*c0**2*(c2-1)+c0*(c0**2-164*c2)+80*c2
       g4 = c0*(52*c2+c0*(c0*(c2-1)-30*c2))
       g5 = 2*c0**2*c2*(6-c0)
       g6 = c0**3*c2
       ap2s = gexp * (eta*(eta*(eta*(eta*(eta*(eta*g6+g5)+g4)+g3)+g2)+g1)+440)
    end if

  end function ap2s

  function ap2v(v,x,ider)
    ! ap2 (3), vector, 0-4th derivative

    real*8, intent(in) :: v(:), x(4)
    integer, intent(in) :: ider
    real*8 :: ap2v(size(v))

    real*8, dimension(size(v)) :: eta, z, gexp, gfac
    real*8, dimension(size(v)) :: g0, g1, g2, g3, g4, g5, g6, gg, gfz
    real*8 :: pfg, c0, c2, ec0, gfc, dum1, dum2
    integer :: i, ierr, ierr2

    interface
       subroutine dgam(A, X, ACC, G, GSTAR, IFLG, IFLGST)
         double precision a, x, g, gstar
         real acc
         integer iflg, iflgst
       end subroutine dgam
    end interface

    if (nelectrons == 0) &
       call error('ap2v','nelectrons is necessary for ap2 fits',faterr)

    eta = (v/x(2))**third
    pfg = (3*pisquare)**twothird/5d0 * (nelectrons/x(2))**(5d0/3d0)
    c0 = -log(3*x(3)/pfg)
    c2 = 1.5d0 * (x(4)-3) - c0

    if (ider == 0) then
       z = eta*c0
       ec0 = exp(c0)
       gfz = gammai(-2d0,z)
       gfc = gammai(-2d0,c0)
       g2 = (gfz-gfc) * (c0**2 * ec0)
       gfz = gammai(-1d0,z)
       gfc = gammai(-1d0,c0)
       g1 = (gfz-gfc) * (c0 * (c2-1) * ec0)
       gfz = gammai(0d0,z)
       gfc = gammai(0d0,c0)
       g0 = (gfz-gfc) * (-2 * c2 * ec0)
       gg = (exp(-z+c0)-1)*(c2/c0)
       ap2v = (g2+g1+g0+gg) * (9*x(2)*x(3)) + x(1)
    else if (ider == 1) then
       gexp = -exp((-eta+1)*c0)*(3*x(3)) / eta**5
       ap2v = gexp * (-eta+1) * ((-eta+1)*eta*c2+1)
    else if (ider == 2) then
       gfac = -1d0/(eta**3 * 3 * x(2))
       gexp = -exp((-eta+1)*c0)*(3*x(3)) / eta**5 * gfac
       g1 = 4*c2+c0-4
       g2 = (c0-6)*c2-c0
       g3 = 2*c2*(1-c0)
       g4 = c0*c2
       ap2v = gexp * (eta*(eta*(eta*(eta*g4+g3)+g2)+g1)+5)
    else if (ider == 3) then
       gfac = -1d0/(eta**3 * 3 * x(2))
       gexp = -exp((-eta+1)*c0)*(3*x(3)) / eta**5 * gfac**2
       g1 = 28*(c2-1)+12*c0
       g2 = 10*c0*(c2-1)-36*c2+c0**2
       g3 = c0**2*(c2-1)-16*c0*c2+10*c2
       g4 = 2*c0*c2*(3-c0)
       g5 = c0**2*c2
       ap2v = gexp * (eta*(eta*(eta*(eta*(eta*g5+g4)+g3)+g2)+g1)+40)
    else if (ider == 4) then
       gfac = -1d0/(eta**3 * 3 * x(2))
       gexp = -exp((-eta+1)*c0)*(3*x(3)) / eta**5 * gfac**3
       g1 = 280*(c2-1)+160*c0
       g2 = 118*c0*(c2-1)-324*c2+21*c0**2
       g3 = 18*c0**2*(c2-1)+c0*(c0**2-164*c2)+80*c2
       g4 = c0*(52*c2+c0*(c0*(c2-1)-30*c2))
       g5 = 2*c0**2*c2*(6-c0)
       g6 = c0**3*c2
       ap2v = gexp * (eta*(eta*(eta*(eta*(eta*(eta*g6+g5)+g4)+g3)+g2)+g1)+440)
    end if

  end function ap2v

end module evfunc
