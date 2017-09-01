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

subroutine setvariables (line, lp)
  use topcalc
  use varbas
  use evfunc
  use fit
  use tools
  use param
  implicit none
  character*(*), intent(inout) :: line
  integer, intent(inout) :: lp

  integer           :: ipid
  character*(mline) :: word
  logical           :: ok
  
  word = getword(word,line,lp)
  word = lower(word)

  if (equal (word,'root'//null)) then
     fileroot = getword(fileroot,line,lp)
     fileroot = trim(adjustl(fileroot(1:leng(fileroot))))
  elseif (equal (word,'pfit_mode'//null)) then
     word = getword(word,line,lp)
     word = lower(word)
     if (equal(word,'gauss'//null)) then
        pfit_mode = pfit_gauss
     elseif (equal(word,'slatec'//null)) then
        pfit_mode = pfit_slatec
     else
        call error ('setvariables','Unknown SET PFIT_MODE option',faterr)
     end if
  elseif (equal (word,'pweigh_mode'//null)) then
     word = getword(word,line,lp)
     word = lower(word)
     if (equal(word,'gibbs1'//null)) then
        pweigh_mode = pweigh_gibbs1
     elseif (equal(word,'gibbs2'//null)) then
        pweigh_mode = pweigh_gibbs2
     elseif (equal(word,'slatec'//null)) then
        pweigh_mode = pweigh_slatec
     else
        call error ('setvariables','Unknown SET PWEIGH_MODE option',faterr)
     end if
  elseif (equal (word,'noefit'//null)) then
     doefit = .false.
  elseif (equal (word,'noplotdh'//null)) then
     doplotdh = .false.
  elseif (equal (word,'mpar'//null)) then
     ok = isinteger(mpar,line,lp)
     if (.not.ok) &
        call error('setvariables','Wrong MPAR value',faterr)
     if (mpar+2 > mmpar) then
        write (uout,'(" MPAR = ",I4)') mpar
        call error('setvariables','MPAR+2 can not exceed mmpar',faterr)
     end if
  elseif (equal (word,'mparmin'//null)) then
     ok = isinteger(mparmin,line,lp)
     if (.not.ok) &
        call error('setvariables','Wrong MPARMIN value',faterr)
  elseif (equal (word,'ndel'//null)) then
     ok = isinteger(ndel,line,lp)
     if (.not.ok) &
        call error('setvariables','Wrong NDEL value',faterr)
  elseif (equal (word,'newpts'//null)) then
     ok = isinteger(newpts,line,lp)
     if (.not.ok) &
        call error('setvariables','Wrong NEWPTS value',faterr)
  elseif (equal (word,'facexpand'//null)) then
     ok = isreal(facexpand,line,lp)
     if (.not.ok) &
        call error('setvariables','Wrong FACEXPAND value',faterr)
  elseif (equal (word,'notrans'//null)) then
     dotrans = .false.
  elseif (equal(word,'errorbar'//null) .or. equal(word,'errorbars'//null) .or. &
          equal(word,'error_bar'//null) .or. equal(word,'error_bars'//null)) then
     doerrorbar = .true.
  elseif (equal (word,'phonfit'//null)) then
     word = getword(word,line,lp)
     word = lower(word)
     phonsplin = .not.(equal(word,'linear'//null))
  elseif (equal (word,'writelevel'//null)) then
     ok = isinteger(writelevel,line,lp)
     if (.not.ok) &
        call error('setvariables','Wrong WRITELEVEL value',faterr)
  elseif (equal(word,'ignore_neg_cutoff'//null)) then
     ok = isreal(ignore_neg_cutoff,line,lp)
     if (.not.ok) &
        call error('setvariables','Wrong IGNORE_NEG_CUTOFF',faterr)
  else
     call error ('setvariables','Unknown set option or keyword',faterr)
  endif
  
end subroutine setvariables
