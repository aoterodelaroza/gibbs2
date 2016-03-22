! Copyright (c) 2011 Alberto Otero de la Roza <alberto@carbono.quimica.uniovi.es>,
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

module eos
  use param
  use tools
  use varbas
  implicit none
  private

  ! public
  public :: eosfit_ev_fitt

  ! info for optm, spinodal EOS
  logical :: spinodal_opt
  real*8 :: bcnt_vstatic, bcnt_gbao
  real*8, allocatable :: bcnt_p(:), bcnt_b(:)
  real*8 :: xkopt

  ! info for minpack
  real*8, allocatable :: xm(:), ym(:)

contains

  !> Fit by weighed polynomials to static data. Output results.
  subroutine eosfit_ev_fitt(p)
    use evfunc

    type(phase), intent(inout) :: p

    integer :: i, j, ierr
    real*8 :: vk, gk, ek, bk
    real*8 :: f1, f2, f3, f4, pt, tmp, tmp2, b1, b2
    real*8 :: a(p%nv)
    character*(mline_fmt) :: fm, fme
    integer :: luw
    logical, save :: luw_open = .false.
    real*8, dimension(8) :: prop, prop2

    if (writelevel > 0) then
       if (luw_open) then
          luw = fopen(luw,trim(fileroot)//".eos_static"//null,ioappend)
       else
          luw = fopen(luw,trim(fileroot)//".eos_static"//null,iowrite)
          luw_open = .true.
       end if
    end if
    write (uout,'("# Copy in file : ",A)') trim(fileroot)//".eos_static"

    if (p%pfit%nfit > 0) then
       write (uout,'("# Lines beginning with ''e'' contain fit error estimation.")')
       if (writelevel > 0) then
          write (luw,'("# Lines beginning with ''e'' contain fit error estimation.")')
       end if
    end if
    fm = format_string_header( &
       (/2,ifmt_p,ifmt_eprec,ifmt_eprec,ifmt_v,ifmt_x,ifmt_p,ifmt_b,ifmt_bp,ifmt_bpp/),&
       (/2,6,6,6,9,4,10,6,2,10/))
    write (uout,fm) "#","p(GPa)","E(Ha)","H(Ha)","V(bohr^3)","V/V0","p_fit(GPa)","B(GPa)",&
       "Bp","Bpp(GPa-1)"

    if (writelevel > 0) then
       write (luw,'("# Phase ",A)') trim(adjustl(p%name(1:leng(p%name))))
       write (luw,fm) "#","p(GPa)","E(Ha)","H(Ha)","V(bohr^3)","V/V0","p_fit(GPa)","B(GPa)",&
          "Bp","Bpp(GPa-1)"
    end if

    fm  = format_string((/ifmt_p,ifmt_eprec,ifmt_eprec,ifmt_v,ifmt_x,&
       ifmt_p,ifmt_b,ifmt_bp,ifmt_bpp/),2)
    fme  = format_string((/1,ifmt_p,ifmt_eprec,ifmt_eprec,ifmt_v,ifmt_x,&
       ifmt_p,ifmt_b,ifmt_bp,ifmt_bpp/),0)
    do i = 1, nps
       call fit_pshift(p%fit_mode,p%v,plist(i),p%npol,p%cpol,vk,bk,ek,gk,ierr)
       if (ierr > 0) cycle
       f1 = fv1(p%fit_mode,vk,p%npol,p%cpol)
       f2 = fv2(p%fit_mode,vk,p%npol,p%cpol)
       f3 = fv3(p%fit_mode,vk,p%npol,p%cpol)
       f4 = fv4(p%fit_mode,vk,p%npol,p%cpol)
       pt = -f1 * au2gpa
       b1 = -(1+vk*f3/f2)
       b2 = ((f3+vk*f4)/f2**2 - vk*f3**2/f2**3) / au2gpa

       write (uout,fm) plist(i),ek,gk,vk,vk/p%veq_static,pt,bk,b1,b2
       if (writelevel > 0) then
          write (luw,fm) plist(i),ek,gk,vk,vk/p%veq_static,pt,bk,b1,b2
       end if

       if (p%pfit%nfit > 0) then
          prop = 0d0
          prop2 = 0d0
          do j = 1, p%pfit%nfit
             call fit_pshift(p%pfit%mode(j),p%v,plist(i),p%pfit%npar(j),p%pfit%apar(:,j),vk,bk,ek,gk,ierr)
             f1 = fv1(p%pfit%mode(j),vk,p%pfit%npar(j),p%pfit%apar(:,j))
             f2 = fv2(p%pfit%mode(j),vk,p%pfit%npar(j),p%pfit%apar(:,j))
             f3 = fv3(p%pfit%mode(j),vk,p%pfit%npar(j),p%pfit%apar(:,j))
             f4 = fv4(p%pfit%mode(j),vk,p%pfit%npar(j),p%pfit%apar(:,j))
             pt = -f1 * au2gpa
             b1 = -(1+vk*f3/f2)
             b2 = ((f3+vk*f4)/f2**2 - vk*f3**2/f2**3) / au2gpa
             prop = prop + (/ek,gk,vk,vk/p%veq_static,pt,bk,b1,b2/) * p%pfit%wei(j)
             prop2 = prop2 + (/ek,gk,vk,vk/p%veq_static,pt,bk,b1,b2/)**2 * p%pfit%wei(j)
          end do
          prop = sqrt(max(prop2 - prop*prop,0d0))
          write (uout,fme) "e",plist(i),prop
          if (writelevel > 0) then
             write (luw,fme) "e",plist(i),prop
          end if
       end if
    end do
    
    if (allocated(vlist)) then
       write (uout,'("# --- Volume data --- ")')
       if (writelevel > 0) then
          write (luw,'("# --- Volume data --- ")')
       end if
       do i = 1, nvs
          vk = vlist(i)
          f1 = fv1(p%fit_mode,vk,p%npol,p%cpol)
          f2 = fv2(p%fit_mode,vk,p%npol,p%cpol)
          f3 = fv3(p%fit_mode,vk,p%npol,p%cpol)
          f4 = fv4(p%fit_mode,vk,p%npol,p%cpol)

          ek = fv0(p%fit_mode,vk,p%npol,p%cpol)
          gk = ek - f1 * vk
          pt = -f1 * au2gpa
          bk = vk * f2 * au2gpa
          b1 = -(1+vk*f3/f2)
          b2 = ((f3+vk*f4)/f2**2 - vk*f3**2/f2**3) / au2gpa

          write (uout,fm) pt,ek,gk,vk,vk/p%veq_static,pt,bk,b1,b2
          if (writelevel > 0) then
             write (luw,fm) pt,ek,gk,vk,vk/p%veq_static,pt,bk,b1,b2
          end if
       end do
    end if
    if (writelevel > 0) then
       write (luw,'(/)') 
       call fclose(luw)
    end if

  end subroutine eosfit_ev_fitt

end module eos
