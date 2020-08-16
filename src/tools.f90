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

module tools
  use param, only: mline
  implicit none
  private

  ! public
  !xx! input/output
  public :: ioinit
  public :: getargs, stdargs
  public :: fopen, fclose
  public :: falloc, fdealloc
  public :: fgetline, getword
  public :: islogical, isreal, isinteger
  !xx! timer routines
  public :: timer
  public :: ffdate
  !xx! errors, warnings, comments
  public :: error
  !xx! string tools
  public :: equal, cat, leng, lower
  !xx! memory management
  public :: realloc
  !xx! sorting
  public :: qcksort
  !xx! math
  public :: gauleg
  public :: gauss
  public :: quad1, trapezoidal, simpson
  public :: gammai
  !xx! symmetry
  public :: laue
  !xx! system
  public :: spawn

  ! overloaded functions
  interface realloc
     module procedure realloc1r
     module procedure realloc1i
     module procedure realloc1l
     module procedure realloc2r
     module procedure realloc3r
  end interface
  interface gammai
     module procedure gammais
     module procedure gammaiv
  end interface gammai

  ! private
  !xx! input/output
  private :: fflush
  private :: fputstr
  private :: isdigit, atoi, atof
  !xx! timer routines
  private :: atimer
  !xx! errors, warnings, comments
  !xx! string tools

  ! private
  integer, parameter :: mopen = 100  !< maximum number of open lus
  logical, private :: alloc (0:mopen) !< allocation flag array
  integer, private :: access(0:mopen) !< access type array
  character*(mline), private :: buffer(0:mopen) !< i/o line buffer array
  integer, private :: bp(0:mopen) !< buffer pointer array
  integer, private :: length(0:mopen) !< lengths of line buffers
  
contains

  !xx! input/output

  !> Initialize file system and connect standard units
  subroutine ioinit ()
    use param, only: ioread, iowrite, stderr, stdin, stdout
    !.initialize i/o buffers:
    bp = 0
    length = 0
    alloc = .false.
    access = 0

    !.connect standard units:
    alloc(stderr) = .true.
    alloc(stdin) = .true.
    alloc(stdout) = .true.
    access(stderr) = iowrite
    access(stdin) = ioread
    access(stdout) = iowrite

  end subroutine ioinit
  
  !> Get command line arguments and argument count
  !> on return, argc will contain the number of arguments and
  !> argv will be a character array containing the arguments.
  subroutine getargs (argc, argv)
    use param, only: null
    integer, intent(out) :: argc
    character*(*), intent(out) :: argv(*)
    
    character*(mline) :: line  !local string for parsing
    integer :: length         !total length of command line

    argc=0
10  continue
    call getarg (argc+1,line)
    length=leng(line)
    if (length.gt.0) then
       argc=argc+1
       argv(argc)=line(1:length)//null
    endif
    if (length.gt.0) goto 10
  end subroutine getargs

  !> Connect input files to units with standard defaults
  !> local_uin and local_uout will be set to the input and output units.
  !>   command line is interpreted as follows:
  !>   1. if two file names are provided as argv(1) and argv(2), the
  !>      first is taken as the input and the second as the output file.
  !>   2. if only one file is provided as argv(1) it will be opened
  !>      for input, and output will be set to stdout.
  !>   3. if no files are provided, local_uin will be set to stdin, and
  !>      local_uout to stdout.
  !>   4. optv provides user options
  !>   5. More than two file names will be ignored.
  subroutine stdargs (argc, argv, optv, local_uin, local_uout)
    use param, only: ioread, iowrite, stdin, stdout
    character*(*), intent(in) :: argv(*)
    character*(*), intent(out) :: optv
    integer, intent(in) :: argc
    integer, intent(out) :: local_uin, local_uout

    integer :: u(2)
    integer :: a, nopti, n 

    nopti=0
    optv=''
    a=0
    n=0
    if (argc.gt.0) then
       do while (n.lt.argc) 
          n=n+1
          if (argv(n)(1:1).ne.'-') then
             a = a + 1
             if (a.le.2) u(a)=n
          else
             optv=optv(1:leng(optv))//argv(n)(2:leng(argv(n)))
             nopti=nopti+1
          endif
       enddo
    endif
    if (a.eq.0) then
       local_uin = stdin
       local_uout = stdout
    else if (a.eq.1) then
       local_uin = fopen (local_uin, argv(u(1)), ioread)
       local_uout = stdout
    else if (a.ge.2) then
       local_uin = fopen (local_uin, argv(u(1)), ioread)
       local_uout = fopen (local_uout, argv(u(2)), iowrite)
    endif
    return
  end subroutine stdargs

  !> Open file and assign unique unit number.
  !> assigned number returned thru both, fopen and u.
  !> files can be opened in one of several modes:
  !>  ioread ......... read only access
  !>  iowrite ........ write access with ascii carriage control
  !>  ioappend ....... write access with append to existing file
  !>  iofortran ...... write access with fortran carriage control
  integer function fopen (u, filename, mode)
    use param, only: eol, faterr, ioappend, ioerror, iofortran, ioread, iowrite, null, stderr, &
       uout
    character*(*), intent(in) :: filename
    integer, intent(in) :: mode
    integer, intent(out) :: u
    
    integer :: flength,ios

    flength = index(filename,null) - 1
    u = falloc(u)
    fopen = u
    !
    if (u.eq.ioerror) then
       continue
       !
    else if (mode.eq.ioread) then
       open (unit=u,file=filename(1:flength),status='old',iostat=ios,err=10)
       access(u)=ioread
       bp(u)=0
       length(u)=0
       !
    else if (mode.eq.iowrite) then
       open (unit=u,file=filename(1:flength),status='unknown',iostat=ios,err=10)
       access(u)=iowrite
       bp(u)=0
       length(u)=0
       !
    else if (mode.eq.ioappend) then
       open (unit=u,file=filename(1:flength),status='old',iostat=ios,access='append',err=10)
       access(u)=iowrite
       bp(u)=0
       length(u)=0
       !
    else if (mode.eq.iofortran) then
       open (unit=u,file=filename(1:flength),status='unknown',iostat=ios,err=10)
       access(u)=iofortran
       bp(u)=0
       length(u)=0
       !
    else
       call fdealloc (u)
       call fputstr (stderr,'    E: fopen: illegal access mode'//eol)
       fopen=ioerror
    endif

    return

10  call fdealloc (u)
    fopen=ioerror
    call fputstr (stderr,'    E: fopen: '//null)
    if (ios > 0) then
       write (uout,'("File access error : ",A)') filename(1:flength)
       call error("fopen","file access error",faterr)
    end if
    return
  end function fopen

  !> Close and deallocate a file unit.
  subroutine fclose (u)
    use param, only: eol, ioread, null, stderr
    integer, intent(in) :: u
    integer :: ios

    if (u.le.0 .or. u.gt.mopen) then
       call fputstr (stderr,'fclose: illegal unit number'//eol)
    else
       if (access(u) .ne. ioread) then
          call fflush (u)
       endif
       close (unit=u, iostat=ios, err=10)
       call fdealloc (u)
    endif
    return
    !
10  call fputstr (stderr,'fclose: '//null)
    return
  end subroutine fclose

  !> Allocate a unique unit number for i/o unit number is returned
  !> trhu falloc and u simultaneously.
  integer function falloc (u)
    use param, only: eol, ioerror, stderr
    integer, intent(out) :: u

    u = min(stderr,1)
10  if (u.le.mopen .and. alloc(u)) then
       u = u + 1
       goto 10
    endif
    !
    if (u.gt.mopen) then
       u = ioerror
       call fputstr (stderr,'falloc: exceded open file limit'//eol)
    else
       alloc(u) = .true.
    endif
    falloc = u
    return
  end function falloc

  !> Deallocate a unit
  subroutine fdealloc (u)
    use param, only: eol, stderr
    integer, intent(in) :: u

    if (u.le.min(0,stderr-1) .or. u.gt.mopen) then
       call fputstr (stderr,'fdealloc: illegal unit number'//eol)
    else
       alloc(u)=.false.
    endif
    return
  end subroutine fdealloc

  !> Dump buffer associated with unit u (private).
  subroutine fflush (u)

    integer, intent(in) :: u

    if (bp(u).gt.0) then
       write (u,'(a)') buffer(u)(1:bp(u))
    endif
    return
  end subroutine fflush

  !> Get line from file with unit u. Returns true if read was successful.
  logical function fgetline (u,oline)
    use param, only: backsl, eof, newline, null, tab
    character*(mline), intent(out) :: oline
    integer, intent(in) :: u

    integer :: i,lenu
    character*(mline) :: line
    logical :: notfirst

    oline = ""
    notfirst = .false.
    do while (.true.)
       ! read the line
       read (u,'(a)',end=10) line

       ! remove tabs
       do i = 1, len(line)
          if (line(i:i) == tab) line(i:i) = " "
       end do
       ! remove blanks
       line = trim(adjustl(line))
       lenu = len(trim(adjustl(line)))

       ! comments
       if (line(1:1) == "#") then
          if (notfirst) then
             exit
          else
             cycle
          end if
       end if
       ! blank lines
       if (lenu == 0) then
          if (notfirst) then
             exit
          else
             cycle
          end if
       end if

       ! continuation
       lenu = leng(trim(adjustl(line)))
       if (line(lenu:lenu) /= backsl) then
          line(lenu+1:lenu+1) = newline
          line(lenu+2:lenu+2) = null
          oline = trim(adjustl(oline)) // " " // trim(adjustl(line))
          exit
       end if
       oline = trim(adjustl(oline)) // " " // line(1:lenu-1)
       notfirst = .true.
    end do

    fgetline = .true.
    return
    ! end-of-u-file

10  fgetline = .false.
    line(1:1) = eof
    line(2:2) = null
    oline = line
    return
  end function fgetline

  !> Get next word from line at lp and increase lp
  !>   - a word is defined as any sequence on nonblanks
  !>   - word and getword will return the same
  function getword (word, line, lp)
    use param, only: blank, newline, null, tab
    character*(*), intent(out) :: word
    character*(*), intent(in) :: line
    integer, intent(inout) :: lp
    character*(len(word)) :: getword

    integer           wp, l, i

    l = len(word)
    do while (line(lp:lp) .eq. blank .or. line(lp:lp).eq.tab)   !skip blanks
       lp = lp + 1
    enddo
    if (line(lp:lp).eq.newline) then
       !        word = newline//null
       word = null
    else
       wp = 1
       do while ( wp.lt.l .and. line(lp:lp).ne.blank .and. line(lp:lp).ne.tab .and.&
          line(lp:lp).ne.null .and. line(lp:lp).ne.newline )
          word(wp:wp) = line(lp:lp)
          wp = wp + 1
          lp = lp + 1
       enddo
       word(wp:wp) = null

       !.clean rest of the word:
       !
       do i = wp+1, l
          word(i:i) = blank
       enddo
    endif
    getword = word

  end function getword

  !> Set lval to the logical value of the next word in
  !> line, if any. Only fortran style logical constants are recognized,
  !> i.e. .TRUE. or .FALSE.
  logical function islogical (lval, line, lp)
    use param, only: blank, null
    logical, intent(out) :: lval
    integer, intent(inout) :: lp
    character*(*), intent(in) :: line

    character*(len(line)) :: tempstr
    integer           i, j
    logical           found

    i = lp
    do while (line(i:i) .eq. blank)
       i = i + 1
    enddo
    j = index (line, null) - 1
    found = .false.

    if (j-i .ge. 6) then
       tempstr = line(i:i+6) // null
       if (lower (tempstr) .eq. ('.false.' // null)) then
          found = .true.
          lval = .false.
          lp = i + 7
       else 
          tempstr = line(i:i+5) // null
          if (lower (tempstr) .eq. ('.true.' // null)) then
             found = .true.
             lval = .true.
             lp = i + 6
          end if
       endif
    else if (j-i .ge. 5) then
       tempstr = line(i:i+5) // null
       if (lower (tempstr) .eq. ('.true.' // null)) then
          found = .true.
          lval = .true.
          lp = i + 6
       endif
    endif
    islogical = found

    return
  end function islogical

  !> Get integer value from input text. If a valid integer is not
  !> found, then return .false.
  logical function isinteger (ival,line,lp)
    use param, only: blank, newline, null
    integer, intent(out) :: ival
    character*(mline), intent(in) :: line
    integer, intent(inout) :: lp

    integer :: i

    do while (line(lp:lp) .eq. blank)
       lp=lp+1
    enddo
    i=lp
    if (line(i:i) .eq. '+' .or. line(i:i) .eq. '-') i=i+1
    if (isdigit(line(i:i))) then 
       do while (isdigit(line(i:i)))
          i=i+1
       enddo
       if (line(i:i) .eq. blank .or. line(i:i) .eq. null .or. line(i:i).eq.newline ) then
          ival=atoi(line(lp:i-1))
          lp = i
          isinteger=.true.
       else
          ival=0
          isinteger=.false.
       endif
    else
       ival=0
       isinteger=.false.
    endif

    return
  end function isinteger

  !> Get a real number from line and sets rval to it.
  !> If a valid real number is not found, isreal returns .false.
  logical function isreal (rval, line, lp)
    use param, only: blank, newline, null
    real*8, intent(out) :: rval
    character*(*), intent(in) :: line
    integer, intent(inout) :: lp

    character*(mline) dumchar
    integer           tp, i
    character*(1)     ch
    logical           matched, isdigit
    isdigit(ch) = ch.ge.'0' .and. ch.le.'9'


    do while (line(lp:lp) .eq. blank)
       lp = lp + 1
    end do

    i = lp
    if (line(i:i) .eq. '+' .or. line(i:i) .eq. '-') i = i + 1 
    if (isdigit(line(i:i))) then 
       do while (isdigit(line(i:i)))
          i = i + 1
       enddo
       if (line(i:i) .eq. '.') then 
          i = i + 1
          do while (isdigit(line(i:i)))
             i = i + 1
          enddo
       endif
       matched = .true.
    else if (line(i:i) .eq. '.') then 
       i = i + 1
       if (isdigit(line(i:i))) then 
          do while (isdigit(line(i:i)))
             i = i + 1
          enddo
          matched = .true.
       else
          matched = .false.
       endif
    else
       matched = .false.
    endif

    !.....get optional exponent
    tp = i - 1
    if (matched) then 
       if (line(i:i)=='e' .or. line(i:i)=='E' .or. line(i:i)=='d' .or. line(i:i)=='D'.or.&
           line(i:i)=='-' .or. line(i:i)=='+') then 
          i = i + 1
          if (line(i:i) .eq. '+' .or. line(i:i) .eq. '-') i = i + 1 
          if (isdigit (line(i:i))) then 
             do while (isdigit(line(i:i)))
                i = i + 1
             enddo
             if (index(blank//','//null//newline, line(i:i)) .gt. 0)          then
                dumchar=line(lp:i-1)//null
                rval = atof (dumchar)
                lp = i
             else
                matched = .false.
                rval = 0d0
             endif
          else 
             matched = .false.
          endif
       else
          if (index(blank//','//null//newline, line(i:i)) .gt. 0)          then
             dumchar=line(lp:tp)//null
             rval = atof(dumchar)
             lp = i
          else
             matched = .false.
             rval = 0d0
          endif
       endif
    else
       rval = 0d0
    endif
    !
    isreal = matched
    return
  end function isreal

  !> Return true if c is a digit (private).
  logical function isdigit (c)

    character*(1), intent(in) :: c

    isdigit = c.ge.'0' .and. c.le.'9'

  end function isdigit

  !> Convert ascii string to integer (private).
  integer function atoi (string)
    use param, only: blank
    character*(*), intent(in) :: string

    integer           i, sign

    atoi=0
    i=1
    do while (string(i:i).eq.blank)
       i=i+1
    end do
    sign=1
    if(string(i:i) .eq. '+' .or. string(i:i) .eq. '-') then 
       if (string(i:i) .eq. '-') sign=-1
       i=i+1
    endif
    do while (isdigit(string(i:i)))
       atoi=10*atoi+ichar(string(i:i))-ichar('0')
       i=i+1
       if (i > len(string)) then
          exit
       end if
    end do
    atoi=atoi*sign
    return

  end function atoi

  !> Convert ascii string to real (private).
  real*8 function atof(str)
    use param, only: blank
    character*(*),intent(in) :: str

    real*8, parameter :: ten=10d0

    real*8 val,power
    integer exponent,sign,esign,i

    sign=1
    val=0d0
    power=1d0
    exponent=0
    esign=1
    i=1
    do while (str(i:i) .eq. blank)
       i=i+1 
    end do
    if (str (i:i) .eq. '+' .or. str(i:i) .eq.'-') then 
       if (str(i:i) .eq. '-') sign=-1
       i=i+1
    endif
    do while (isdigit(str(i:i)))
       val=ten*val+ichar(str(i:i))-ichar('0')
       i=i+1
    enddo
    if (str(i:i) .eq. '.') then 
       i=i+1
       do while (isdigit(str(i:i)))
          val=ten*val+ichar(str(i:i))-ichar('0')
          i=i+1
          power=power*ten
       enddo
    endif
    if (str(i:i).eq.'e' .or. str(i:i).eq.'E' .or. str(i:i).eq.'d' .or. str(i:i).eq.'D') then
       i=i+1
       if (str (i:i) .eq. '+' .or. str(i:i) .eq.'-') then 
          if (str(i:i) .eq. '-') esign=-1
          i=i+1
       endif
       do while (isdigit(str(i:i)))
          exponent=ten*exponent+ichar(str(i:i))-ichar('0')
          i=i+1
       end do
    elseif (str(i:i).eq.'-' .or. str(i:i).eq.'+') then
       esign = 1
       if (str(i:i) .eq. '-') esign=-1
       i = i + 1 
       do while (isdigit(str(i:i)))
          exponent=ten*exponent+ichar(str(i:i))-ichar('0')
          i=i+1
       end do
    endif

    atof=(sign*val/power)*ten**(esign*exponent)
    return
  end function atof

  !> Output string to the u file unit (private).
  subroutine fputstr (u, string)
    use param, only: newline, null
    character*(*), intent(in) :: string
    integer, intent(in) :: u
    integer i

    i=1
10  if (string(i:i).ne.null) then
       bp(u)=bp(u)+1
       if (bp(u).eq.mline) then
          write (u,'(a)') buffer(u)(1:mline)
          bp(u)=0
       else if (string(i:i).eq.newline) then
          write (u,'(a)') buffer(u)(1:bp(u)-1)
          i=i+1
          bp(u)=0
       else
          buffer(u)(bp(u):bp(u))=string(i:i)
          i=i+1
       endif
       goto 10
    endif
    return

  end subroutine fputstr

  !xx! timer routines

  !> Accumulates and prints out the elapsed times of a series of up to
  !> MPROCESS processes.
  !> Input parameters are:
  !>   key ........... 0 = reset all time tables.
  !>                   1 = reset pid entry and begin counting for it.
  !>                   2 = continue the cont for pid process. Do not
  !>                       reset previous entry times.
  !>                   3 = end of pid process. Free pid entry.
  !>                   4 = end of pid process. Do not free entry.
  !>                   5 = end of run. Print out all time tables.
  !>                   6 = close all processes and print out tables.
  !>   pid ........... process identification number (1..MPROCESS).
  !>   name .......... process name (used only in the print out).
  !>   lw ............ printer logical unit. Negative if print out is not
  !>                   desired.
  !> key controls what data are used in the run, as the following table
  !> resumes:
  !>   key value    pid      name       lw
  !>   ---------  -------   -------   -------
  !>      0       ignored   ignored   ignored
  !>     1,2      output    input     ignored
  !>     3,4      input     ignored   input
  !>     5,6      ignored   ignored   input
  subroutine timer (key,pid,name,lw)
    use param, only: warning
    integer, intent(in) :: key, lw
    integer, intent(inout) :: pid
    character*(*), intent(in) :: name

    character*(80)   message
    integer          mprocess
    parameter        (mprocess=200)
    character*(20)    pname(mprocess)
    real*4           time(mprocess),cumtime(mprocess),timedum
    logical          popen(mprocess),pused(mprocess)
    logical          firsttime
    integer*4        pcalls(mprocess)
    save             time,cumtime,pname,popen,pused,pcalls
    save             firsttime
    integer          i

    data firsttime /.true./

    if (key.eq.0 .or. firsttime) then

       ! initiallize all entries:
       firsttime=.false.
       do  i=1,mprocess
          time(i)=0.0
          cumtime(i)=0.0
          popen(i)=.false.
          pused(i)=.false.
          pname(i)='        '
          pcalls(i)=0
       end do
    else if (key.eq.1 .or. key.eq.2) then

       ! begin pid count:
       call atimer (timedum)
       i=1
15     if (pused(i)) then
          if (pname(i).ne.name) then
             i=i+1
             if (i.gt.mprocess) then
                call error ('timer','too many processes',warning)
                return
             endif
             goto 15
          endif
       endif
       pid=i
       if (key.eq.1) then
          cumtime(pid)=0.0
          pcalls(pid)=1
       else
          pcalls(pid)=pcalls(pid)+1
       endif
       time(pid)=timedum
       popen(pid)=.true.
       pused(pid)=.true.
       pname(pid)=name
    else if (key.eq.3 .or. key.eq.4) then

       ! end pid accounting:
       if (pid.le.0 .or. pid.gt.mprocess) then
          write (message,130) pid, name
          call error ('timer',message,warning)
       else if (.not.popen(pid)) then
          call error ('timer','pid unused or closed',warning)
       else
          call atimer (timedum)
          time(pid)=timedum-time(pid)
          cumtime(pid)=cumtime(pid)+time(pid)
          popen(pid)=.false.
          if (lw.gt.0) write (lw,100) pname(pid),time(pid)
          if (key.eq.3) then
             pused(pid)=.false.
             pcalls(pid)=0
             cumtime(pid)=0.0
             time(pid)=0.0
          endif
       endif
    else if (key.eq.5 .or. key.eq.6) then

       ! print out the time tables:
       write (lw,105)
       call atimer (timedum)
       do i=1,mprocess
          if (pused(i)) then
             if (popen(i)) then
                time(i)=timedum-time(i)
                cumtime(i)=cumtime(i)+time(i)
                if (key.eq.6) popen(i)=.false.
             endif
             write (lw,110) i,pname(i),cumtime(i),pcalls(i),popen(i)
          endif
       end do
       write (lw,115)
    else
       call error ('timer','key value not recognized',warning)
    endif

    return

100 format (/1x,'*** timer:'/1x,'*** process name:',a10,5x,'elapsed time (sec):',f10.3/)
105 format (/1x,'***'/1x,'*** timer:'/1x,'***'/1x,'*** -pid----------name-----------cumtime----pcalls----popen-')
110 format (1x,'***',i4,2x,a20,f11.3,1x,i9,6x,l2)
115 format (1x,'***'/)
130 format ('pid (', i3, ') out of bonds for: ', a)
  end subroutine timer

  !> Interface to the timer routines of different computers. Returns
  !> machine seconds at the calling time. (private)
  subroutine atimer (tiempo)

    real*4          tiempo
    real*4          etime
    real*4          tarray(2)

    tiempo = etime (tarray)

  end subroutine atimer

  !> Outputs date and time to a string.
  function ffdate()

    character*(mline) :: ffdate

    integer :: values(8)

    call date_and_time(values=values)

    write (ffdate,'(I4.4,".",I2.2,".",I2.2,", ",&
       &I2.2,":",I2.2,":",I2.2,".",I3.3)') values(1:3),values(5:8)

  end function ffdate

  !xx! errors, warnings, comments

  !> Send an error message 'message' to stdout, coming from routine
  !> 'routine'. errortype is the error code (see mod_param.f90).
  subroutine error (routine,message,errortype)
    use param, only: faterr, ncomms, nwarns, noerr, stderr, uout, warning
    character*(*), intent(in) :: routine !< routine calling the error
    character*(*), intent(in) :: message !< the message
    integer, intent(in) :: errortype !< fatal, warning or info

    character*(20)      chtype

    ! message styles shamelessly copied from abinit.
    if (errortype.eq.faterr) then
       chtype='ERROR'
    else if (errortype.eq.warning) then
       chtype='WARNING'
    else if (errortype.eq.noerr) then
       chtype='COMMENT'
    else
       chtype='UNKNOWN'
    endif
    write (uout,100) trim(adjustl(chtype)),trim(adjustl(routine)),&
       trim(adjustl(message))
    if (errortype.eq.faterr) then
       write (stderr,100) trim(adjustl(chtype)),trim(adjustl(routine)),&
          trim(adjustl(message))
       stop 1
    else if(errortype.eq.warning) then
       nwarns = nwarns + 1
    else if(errortype.eq.noerr) then
       ncomms = ncomms + 1
    endif

100 format (A,"(",A,"): ",A)

  end subroutine error

  !xx! string tools

  !> Compare two null-terminated strings for equality.
  logical function equal (s,t)
    use param, only: null
    character*(*), intent(in) :: s, t

    integer           i

    i = 1
    do while (s(i:i) .ne. null)
       if (s(i:i).ne.t(i:i)) then
          equal = .false.
          return
       endif
       i = i + 1
    enddo
    !
    if (t(i:i) .eq. null) then
       equal = .true.
    else
       equal = .false.
    endif
    return

  end function equal

  !> Concatenates null-terminated strings
  function cat (str1, str2)
    use param, only: null, eol, stderr
    character*(*), intent(in) :: str1 !< First string
    character*(*), intent(in) :: str2 !< Second string
    integer           i1, i2
    character*(index(str1,null)+index(str2,null)-1) :: cat
    
    i1 = index (str1, null) - 1
    i2 = index (str2, null) - 1
    if (i1+i2 .ge. len(cat)) then
       call fputstr (stderr,'cat: WARNING: cat size less than glued strings'//eol)
    endif
    cat = str1(1:i1)//str2(1:i2)//null

  end function cat

  !> Obtains the leng of the string, assuming that the blanks at
  !> the end are dummy.
  integer function leng (string)
    use param, only: blank, null
    character*(*), intent(in) :: string

    integer :: i

    do i=len(string),1,-1
       if (string(i:i).ne.blank .and. string(i:i).ne.null) then
          leng=i
          return
       endif
    end do
    leng=0
    return
  end function leng

  !> Convert string to lowercase except where quoted
  !> string and lower will return the same
  function lower (string)
    use param, only: null, quote
    character*(*), intent(inout) :: string
    character*(len(string)) :: lower

    integer           i, iadd
    logical           inquote

    iadd = ichar('A') - ichar('a')
    i = 1
    inquote = .false.
    do while (string(i:i) .ne. null)
       if (string(i:i) .eq. quote) then
          inquote = .not. inquote
       else if (.not. inquote) then
          if (string(i:i).ge.'A' .and. string(i:i).le.'Z') then
             string(i:i) = char (ichar(string(i:i)) - iadd)
          endif
       endif
       i = i + 1
    enddo
    lower = string
    return
  end function lower

  !> Adapt the size of an allocatable 1D real*8 array
  subroutine realloc1r(a,nnew)
    use param, only: faterr
    real*8, intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew
    
    real*8, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call error('realloc1r','array not allocated',faterr)
    nold = size(a)
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1r

  !> Adapt the size of an allocatable 1D integer array
  subroutine realloc1i(a,nnew)
    use param, only: faterr
    integer, intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew
    
    integer, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call error('realloc1i','array not allocated',faterr)
    nold = size(a)
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1i

  !> Adapt the size of an allocatable 1D integer array
  subroutine realloc1l(a,nnew)
    use param, only: faterr
    logical, intent(inout), allocatable :: a(:)
    integer, intent(in) :: nnew
    
    logical, allocatable :: temp(:)
    integer :: nold
    
    if (.not.allocated(a)) &
       call error('realloc1l','array not allocated',faterr)
    nold = size(a)
    allocate(temp(nnew))
    
    temp(1:min(nnew,nold)) = a(1:min(nnew,nold))
    call move_alloc(temp,a)

  end subroutine realloc1l

  !> Adapt the size of an allocatable 2D real*8 array
  subroutine realloc2r(a,n1,n2)
    use param, only: faterr
    real*8, intent(inout), allocatable :: a(:,:)
    integer, intent(in) :: n1,n2
    
    real*8, allocatable :: temp(:,:)
    integer :: nold(2)
    
    if (.not.allocated(a)) &
       call error('realloc2r','array not allocated',faterr)
    nold = size(a)
    allocate(temp(n1,n2))
    
    temp(1:min(n1,nold(1)),1:min(n2,nold(2))) = a(1:min(n1,nold(1)),1:min(n2,nold(2)))
    call move_alloc(temp,a)

  end subroutine realloc2r

  !> Adapt the size of an allocatable 3D real*8 array
  subroutine realloc3r(a,n1,n2,n3)
    use param, only: faterr
    real*8, intent(inout), allocatable :: a(:,:,:)
    integer, intent(in) :: n1,n2,n3
    
    real*8, allocatable :: temp(:,:,:)
    integer :: nold(3)
    
    if (.not.allocated(a)) &
       call error('realloc3r','array not allocated',faterr)
    nold = size(a)
    allocate(temp(n1,n2,n3))
    
    temp(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3))) = &
       a(1:min(n1,nold(1)),1:min(n2,nold(2)),1:min(n3,nold(3)))
    call move_alloc(temp,a)

  end subroutine realloc3r

  !> Sort the elements of the array arr in ascending order using the
  !> quicksort algorithm. iord is the initial order of data in arr and 
  !> first and last the intervals for elements to be analyzed. In the output,
  !> iord contains the final order in the array arr.
  subroutine qcksort (arr, iord, first, last)
    use param, only: faterr, zero
    !.....Maximum number of elements to be sorted depends on nstack value:
    !     nstack...... ~2log2(last-first+1)
    real*8, dimension(:), intent(in) :: arr !< array to be sorted
    integer, dimension(:), intent(inout) :: iord !< permutation array
    integer, intent(in) :: first !< First index
    integer, intent(in) :: last  !< Last index

    integer           m, nstack
    parameter         (m=7, nstack=50)
    real*8            fm, fa, fc, fmi
    parameter         (fm=7875d0, fa=211d0, fc=1663d0, fmi=1.2698413d-4)
    integer           istack(nstack), jstack, l, ir, j, na, i, iq
    real*8            fx, a

    jstack=0
    l=first
    ir=last
    fx=zero
10  if(ir-l.lt.m)then
       do j=l+1,ir
          na=iord(j)
          a=arr(na)
          do i=j-1,first,-1
             if(arr(iord(i)).le.a) go to 12
             iord(i+1)=iord(i)
          enddo
          i=first-1
12        iord(i+1)=na
       enddo
       if(jstack.eq.0) then
          return
       end if
       ir=istack(jstack)
       l=istack(jstack-1)
       jstack=jstack-2
    else
       i=l
       j=ir
       fx=mod(fx*fa+fc,fm)
       iq=l+(ir-l+1)*(fx*fmi)
       na=iord(iq)
       a=arr(na)
       iord(iq)=iord(l)
20     continue
21     if (j.ge.first) then 
          if (a.lt.arr(iord(j))) then
             j=j-1
             goto 21
          endif
       endif
       if(j.le.i)then
          iord(i)=na
          go to 30
       endif
       iord(i)=iord(j)
       i=i+1
22     if (i.le.last) then
          if (a.gt.arr(iord(i))) then
             i=i+1
             goto 22
          endif
       endif
       if(j.le.i)then
          iord(j)=na
          i=j
          go to 30
       endif
       iord(j)=iord(i)
       j=j-1
       go to 20
30     jstack=jstack+2
       if(jstack.gt.nstack) call error ('qcksort','Increase nstack',faterr) 
       if(ir-i.ge.i-l)then
          istack(jstack)=ir
          istack(jstack-1)=i+1
          ir=i-1
       else
          istack(jstack)=i-1
          istack(jstack-1)=l
          l=i+1
       endif
    endif
    go to 10
  end subroutine qcksort

  !> Find the Gauss-Legendre nodes and weights for an interval.
  subroutine gauleg (x1,x2,x,w,n)
    use param, only: pi
    real*8, intent(in) :: x1 !< Left limit of the interval
    real*8, intent(in) :: x2 !< Right limit of the interval
    real*8, dimension(n), intent(out) :: x !< Position of the nodes
    real*8, dimension(n), intent(out) :: w !< Weights
    integer, intent(in) :: n !< Number of nodes

    real*8, parameter :: eps = 3.0d-16

    integer :: m
    real*8 :: xm, xl
    real*8 :: z, p1, p2, p3, pp, z1
    integer :: i, j

    !.The roots are symmetric in the interval, so we only have to find
    ! half of them.

    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)

    ! Loop over the desired roots.
    do i=1,m
       z=cos(pi*(i-.25d0)/(n+.5d0))
1      continue
       p1=1.d0
       p2=0.d0

       ! Loop up the recurrence relation to get the Legendre polynomial
       ! evaluated at z.
       do j=1,n
          p3=p2
          p2=p1
          p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
       end do

       !.p1 is now the desired Legendre polynomial. We next compute pp,
       ! derivative , by a standard relation involving p2, the polyn-
       ! omial of one lower order.
       pp=n*(z*p1-p2)/(z*z-1.d0)
       z1=z

       !.Newton's method.
       z=z1-p1/pp
       if(abs(z-z1).gt.eps)go to 1

       !.Scale the root to the desired interval.
       x(i)=xm-xl*z

       !.and put in its symmetric counterpart.
       x(n+1-i)=xm+xl*z

       !.compute the weight.
       w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)

       !.and its symmetric counterpart.
       w(n+1-i)=w(i)
    end do

  end subroutine gauleg

  subroutine gauss (c, n, mn, x, ierr)
    !.gauss - solves a system of N linear equations with N variables.
    !
    ! B = A X     ;   B & A are stored together in C = [A ...] [b]
    !                                                  [. ...] [.]
    !                                                  [. ...] [.]
    !
    ! where B is a column vector of N components (independent terms),
    ! A is the N*N matrix of coefficients, and X the N-column vector
    ! of variables to be obtained. The routine uses the Gauss elimina-
    ! tion method, including pivotal substitution.
    !
    !.INPUT parameters:
    ! c ....... c(0:mn,0:mn+1) matrix containing A and b. b must be in
    !           the LAST column of c.
    ! n ....... number of independent equations and variables.
    ! mn ...... physical dimensions on c and x.
    !
    !.OUTPUT parameters:
    ! c ....... triangularized c matrix.
    ! x ....... solutions.
    ! ierr .... Error code: 0=no error, -1=singular matrix, -2=MN>mnt,
    !           -3=less than 1 variable and equation.

    integer, intent(in) :: n, mn
    real*8, intent(inout) :: c(0:mn,0:mn+1)
    real*8, intent(out) :: x(0:mn)
    integer, intent(out) :: ierr

    real*8 :: pmax, cdum
    integer :: i, j, k, ip, ii, ic

    ierr = 0
    if (n.lt.0) then
       ierr = -3
    else if (n.eq.0) then
       x(0) = c(0,1) / c(0,0)
    else
       do i = 0, n-1
          !.search for the pivot:
          pmax = abs(c(i,i))
          ip = i
          do ii = i+1, n
             if (abs(c(ii,i)).gt.pmax) then
                pmax = abs(c(ii,i))
                ip = ii
             endif
          enddo
          if (pmax.eq.0d0) then
             ierr = -1
             return
          endif

          !.interchange the column pivot with the i-th one:
          if (ip.ne.i) then
             do ic = i, n+1
                cdum = c(i,ic)
                c(i,ic) = c(ip,ic)
                c(ip,ic) = cdum
             enddo
          endif

          !.Gauss elimination procedure:
          do j = i+1, n
             pmax = c(j,i) / c(i,i)
             do k = i, n+1
                c(j,k) = c(j,k) - pmax * c(i,k)
             enddo
          enddo
       enddo

       !.inverse substitution for obtain the x's:
       x(n) = c(n,n+1) / c(n,n)
       do i = n-1, 0, -1
          x(i) = c(i,n+1)
          do j = i+1, n
             x(i) = x(i) - c(i,j) * x(j)
          enddo
          x(i) = x(i) / c(i,i)
       enddo
    endif
    return
  end subroutine gauss

  ! One-dimensional quadrature. If step > 0, assume x is equally
  ! spaced and use Simpson's rule. If step < 0, use x and y with the 
  ! trapezoidal rule
  function quad1(x,y,step)
    real*8, intent(in) :: x(:), y(:), step
    real*8 :: quad1

    if (step > 0d0) then
       quad1 = simpson(y,step)
    else
       quad1 = trapezoidal(x,y)
    end if

  endfunction quad1

  ! Integrate a function given by y(1:n) at non-equally-spaced points
  ! x(1:n) using the trapezoidal rule.
  function trapezoidal(x,y)

    real*8, intent(in) :: x(:), y(:)
    real*8 :: trapezoidal

    integer :: n

    n = size(x,1)
    if (n < 2) then
       trapezoidal = 0d0
       return
    end if

    trapezoidal = 0.5d0 * sum((x(2:n)-x(1:n-1)) * (y(2:n)+y(1:n-1)))

  end function trapezoidal

  ! Integrate a function given by y(1:n) at equally-spaced points,
  ! separated by step, using Simpson's rule (1424...241)/3. If the
  ! number of points is even, use trapezoidal rule for the last one.
  function simpson(y,step)

    real*8, intent(in) :: y(:), step
    real*8 :: simpson

    integer :: i, n

    n = size(y,1)
    if (n < 2) then
       simpson = 0d0
       return
    elseif (n == 2) then
       simpson = 0.5d0 * step * (y(1) + y(n))
       return
    end if

    simpson = 2 * sum(y(2:n-1)) + y(1) + y(n)
    do i = 2, n-1, 2
       simpson = simpson + 2 * y(i)
    end do
    simpson = simpson * step / 3d0
    if (modulo(n,2) == 0) then
       ! correct the last interval, if necessary
       simpson = simpson + (y(n)+y(n-1)) * step / 6d0
    end if

  end function simpson

  function gammais(a,x)
    use param, only: faterr
    real*8, intent(in) :: a, x
    real*8 :: gammais

    real*8 :: dum1
    integer :: ierr, ierr2

    call dgam(a,x,16.,gammais,dum1,ierr,ierr2)
    if (ierr > 0) call error('gammais','ierr > 0 in dgam',faterr)

    if (a <= 0d0) then
       gammais = gammais / (exp(x)*x**(-a))
    else
       call error('gammais','gammai for a<0 not implemented',faterr)
    end if

  end function gammais

  function gammaiv(a,x)
    use param, only: faterr
    real*8, intent(in) :: a, x(:)
    real*8 :: gammaiv(size(x))

    real*8 :: dum1
    integer :: i, ierr, ierr2

    do i = 1, size(x)
       call dgam(a,x(i),16.,gammaiv(i),dum1,ierr,ierr2)
       if (ierr > 0) call error('gammaiv','ierr > 0 in dgam',faterr)
    end do

    if (a <= 0d0) then
       gammaiv = gammaiv / (exp(x)*x**(-a))
    else
       call error('gammaiv','gammai for a<0 not implemented',faterr)
    end if

  end function gammaiv

  !> Call the operating system to run the command contained in
  !> 'line' in a separate shell. 
  subroutine spawn(line)

    character*(*)     line

    call system (line(1:leng(line)))
    !mdc*else if VMS
    !      integer           status, lib$spawn, flags
    !c
    !      flags = CLI$M_NOKEYPAD
    !      if (leng(line) .gt. 0) then
    !         status = lib$spawn(line(1:leng(line)),,,flags)
    !         if (status .ne. 0) call lib$signal(%val(status))
    !      endif
    !      return
    !mdc*endif
  end subroutine spawn

  subroutine laue(strin,nrot,rot)
    use param, only: faterr, null, uout

    character*(*), intent(in) :: strin
    integer, intent(out) :: nrot
    real*8, intent(out) :: rot(3,3,48)

    character*(len(strin)) :: straux, str
    integer :: i, ic, ilenn

    straux = strin
    straux = lower(straux)
    straux = adjustl(straux)
    ic = 0
    do i = 1, len(trim(straux))
       if (straux(i:i) /= " ") then
          ic = ic + 1
          str(ic:ic) = straux(i:i)
          ilenn = ic
       end if
    end do
    str = str(1:ilenn)

    if (equal(trim(str)//null,"-1"//null) .or. equal(trim(str)//null,"ci"//null)) then
       nrot = 2
       ! 1 : x,y,z           => 1                   
       rot(1,:,1) = (/ 1d0,  0d0,  0d0 /)
       rot(2,:,1) = (/ 0d0,  1d0,  0d0 /)
       rot(3,:,1) = (/ 0d0,  0d0,  1d0 /)
       ! 2 : -x,-y,-z        => -1                  
       rot(1,:,2) = (/-1d0,  0d0,  0d0 /)
       rot(2,:,2) = (/ 0d0, -1d0,  0d0 /)
       rot(3,:,2) = (/ 0d0,  0d0, -1d0 /)
    else if (equal(trim(str)//null,"2/m"//null) .or. equal(trim(str)//null,"c2h"//null)) then
       nrot = 4
       ! 1 : x,y,z           => 1                   
       rot(1,:,1) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,1) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,1) = (/ 0d0,  0d0,  1d0/)
       ! 2 : -x,y,-z         => 2 [ 0 1 0 ]         
       rot(1,:,2) = (/-1d0,  0d0,  0d0/)
       rot(2,:,2) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,2) = (/ 0d0,  0d0, -1d0/)
       ! 3 : -x,-y,-z        => -1                  
       rot(1,:,3) = (/-1d0,  0d0,  0d0/)
       rot(2,:,3) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,3) = (/ 0d0,  0d0, -1d0/)
       ! 4 : x,-y,z          => m [ 0 1 0 ]         
       rot(1,:,4) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,4) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,4) = (/ 0d0,  0d0,  1d0/)
    else if (equal(trim(str)//null,"mmm"//null).or.equal(trim(str)//null,"d2h"//null)) then
       nrot = 8
       ! 1 : x,y,z           => 1                   
       rot(1,:,1) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,1) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,1) = (/ 0d0,  0d0,  1d0/)
       ! 2 : -x,-y,z         => 2 [ 0 0 1 ]         
       rot(1,:,2) = (/-1d0,  0d0,  0d0/)
       rot(2,:,2) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,2) = (/ 0d0,  0d0,  1d0/)
       ! 3 : -x,y,-z         => 2 [ 0 1 0 ]         
       rot(1,:,3) = (/-1d0,  0d0,  0d0/)
       rot(2,:,3) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,3) = (/ 0d0,  0d0, -1d0/)
       ! 4 : x,-y,-z         => 2 [ 1 0 0 ]         
       rot(1,:,4) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,4) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,4) = (/ 0d0,  0d0, -1d0/)
       ! 5 : -x,-y,-z        => -1                  
       rot(1,:,5) = (/-1d0,  0d0,  0d0/)
       rot(2,:,5) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,5) = (/ 0d0,  0d0, -1d0/)
       ! 6 : x,y,-z          => m [ 0 0 1 ]         
       rot(1,:,6) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,6) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,6) = (/ 0d0,  0d0, -1d0/)
       ! 7 : x,-y,z          => m [ 0 1 0 ]         
       rot(1,:,7) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,7) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,7) = (/ 0d0,  0d0,  1d0/)
       ! 8 : -x,y,z          => m [ 1 0 0 ]         
       rot(1,:,8) = (/-1d0,  0d0,  0d0/)
       rot(2,:,8) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,8) = (/ 0d0,  0d0,  1d0/)
    else if (equal(trim(str)//null,"-3"//null).or.equal(trim(str)//null,"c3i"//null)) then
       nrot = 6
       ! 1 : x,y,z           => 1                   
       rot(1,:,1) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,1) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,1) = (/ 0d0,  0d0,  1d0/)
       ! 2 : -y,x-y,z        => 3+ [ 0 0 1 ]        
       rot(1,:,2) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,2) = (/ 1d0, -1d0,  0d0/)
       rot(3,:,2) = (/ 0d0,  0d0,  1d0/)
       ! 3 : -x+y,-x,z       => 3- [ 0 0 1 ]        
       rot(1,:,3) = (/-1d0,  1d0,  0d0/)
       rot(2,:,3) = (/-1d0,  0d0,  0d0/)
       rot(3,:,3) = (/ 0d0,  0d0,  1d0/)
       ! 4 : -x,-y,-z        => -1                  
       rot(1,:,4) = (/-1d0,  0d0,  0d0/)
       rot(2,:,4) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,4) = (/ 0d0,  0d0, -1d0/)
       ! 5 : y,-x+y,-z       => -3+ [ 0 0 1 ]       
       rot(1,:,5) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,5) = (/-1d0,  1d0,  0d0/)
       rot(3,:,5) = (/ 0d0,  0d0, -1d0/)
       ! 6 : x-y,x,-z        => -3- [ 0 0 1 ]       
       rot(1,:,6) = (/ 1d0, -1d0,  0d0/)
       rot(2,:,6) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,6) = (/ 0d0,  0d0, -1d0/)
    else if (equal(trim(str)//null,"4/m"//null).or.equal(trim(str)//null,"c4h"//null)) then
       nrot = 8
       ! 1 : x,y,z           => 1                   
       rot(1,:,1) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,1) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,1) = (/ 0d0,  0d0,  1d0/)
       ! 2 : -x,-y,z         => 2 [ 0 0 1 ]         
       rot(1,:,2) = (/-1d0,  0d0,  0d0/)
       rot(2,:,2) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,2) = (/ 0d0,  0d0,  1d0/)
       ! 3 : -y,x,z          => 4+ [ 0 0 1 ]        
       rot(1,:,3) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,3) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,3) = (/ 0d0,  0d0,  1d0/)
       ! 4 : y,-x,z          => 4- [ 0 0 1 ]        
       rot(1,:,4) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,4) = (/-1d0,  0d0,  0d0/)
       rot(3,:,4) = (/ 0d0,  0d0,  1d0/)
       ! 5 : -x,-y,-z        => -1                  
       rot(1,:,5) = (/-1d0,  0d0,  0d0/)
       rot(2,:,5) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,5) = (/ 0d0,  0d0, -1d0/)
       ! 6 : x,y,-z          => m [ 0 0 1 ]         
       rot(1,:,6) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,6) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,6) = (/ 0d0,  0d0, -1d0/)
       ! 7 : y,-x,-z         => -4+ [ 0 0 1 ]       
       rot(1,:,7) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,7) = (/-1d0,  0d0,  0d0/)
       rot(3,:,7) = (/ 0d0,  0d0, -1d0/)
       ! 8 : -y,x,-z         => -4- [ 0 0 1 ]       
       rot(1,:,8) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,8) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,8) = (/ 0d0,  0d0, -1d0/)
    else if (equal(trim(str)//null,"4/mmm"//null).or.equal(trim(str)//null,"d4h"//null)) then
       nrot = 16
       ! 1 : x,y,z           => 1                   
       rot(1,:,1) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,1) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,1) = (/ 0d0,  0d0,  1d0/)
       ! 2 : -x,-y,z         => 2 [ 0 0 1 ]         
       rot(1,:,2) = (/-1d0,  0d0,  0d0/)
       rot(2,:,2) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,2) = (/ 0d0,  0d0,  1d0/)
       ! 3 : -y,x,z          => 4+ [ 0 0 1 ]        
       rot(1,:,3) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,3) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,3) = (/ 0d0,  0d0,  1d0/)
       ! 4 : y,-x,z          => 4- [ 0 0 1 ]        
       rot(1,:,4) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,4) = (/-1d0,  0d0,  0d0/)
       rot(3,:,4) = (/ 0d0,  0d0,  1d0/)
       ! 5 : -x,y,-z         => 2 [ 0 1 0 ]         
       rot(1,:,5) = (/-1d0,  0d0,  0d0/)
       rot(2,:,5) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,5) = (/ 0d0,  0d0, -1d0/)
       ! 6 : x,-y,-z         => 2 [ 1 0 0 ]         
       rot(1,:,6) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,6) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,6) = (/ 0d0,  0d0, -1d0/)
       ! 7 : y,x,-z          => 2 [ 1 1 0 ]         
       rot(1,:,7) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,7) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,7) = (/ 0d0,  0d0, -1d0/)
       ! 8 : -y,-x,-z        => 2 [ 1 -1 0 ]        
       rot(1,:,8) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,8) = (/-1d0,  0d0,  0d0/)
       rot(3,:,8) = (/ 0d0,  0d0, -1d0/)
       ! 9 : -x,-y,-z        => -1                  
       rot(1,:,9) = (/-1d0,  0d0,  0d0/)
       rot(2,:,9) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,9) = (/ 0d0,  0d0, -1d0/)
       ! 10: x,y,-z          => m [ 0 0 1 ]         
       rot(1,:,10) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,10) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,10) = (/ 0d0,  0d0, -1d0/)
       ! 11: y,-x,-z         => -4+ [ 0 0 1 ]       
       rot(1,:,11) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,11) = (/-1d0,  0d0,  0d0/)
       rot(3,:,11) = (/ 0d0,  0d0, -1d0/)
       ! 12: -y,x,-z         => -4- [ 0 0 1 ]       
       rot(1,:,12) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,12) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,12) = (/ 0d0,  0d0, -1d0/)
       ! 13: x,-y,z          => m [ 0 1 0 ]         
       rot(1,:,13) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,13) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,13) = (/ 0d0,  0d0,  1d0/)
       ! 14: -x,y,z          => m [ 1 0 0 ]         
       rot(1,:,14) = (/-1d0,  0d0,  0d0/)
       rot(2,:,14) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,14) = (/ 0d0,  0d0,  1d0/)
       ! 15: -y,-x,z         => m [ 1 1 0 ]         
       rot(1,:,15) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,15) = (/-1d0,  0d0,  0d0/)
       rot(3,:,15) = (/ 0d0,  0d0,  1d0/)
       ! 16: y,x,z           => m [ 1 -1 0 ]        
       rot(1,:,16) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,16) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,16) = (/ 0d0,  0d0,  1d0/)
    else if (equal(trim(str)//null,"6/m"//null).or.equal(trim(str)//null,"c6h"//null)) then
       nrot = 12
       ! 1 : x,y,z           => 1                   
       rot(1,:,1) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,1) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,1) = (/ 0d0,  0d0,  1d0/)
       ! 2 : -y,x-y,z        => 3+ [ 0 0 1 ]        
       rot(1,:,2) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,2) = (/ 1d0, -1d0,  0d0/)
       rot(3,:,2) = (/ 0d0,  0d0,  1d0/)
       ! 3 : -x+y,-x,z       => 3- [ 0 0 1 ]        
       rot(1,:,3) = (/-1d0,  1d0,  0d0/)
       rot(2,:,3) = (/-1d0,  0d0,  0d0/)
       rot(3,:,3) = (/ 0d0,  0d0,  1d0/)
       ! 4 : -x,-y,z         => 2 [ 0 0 1 ]         
       rot(1,:,4) = (/-1d0,  0d0,  0d0/)
       rot(2,:,4) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,4) = (/ 0d0,  0d0,  1d0/)
       ! 5 : y,-x+y,z        => 6- [ 0 0 1 ]        
       rot(1,:,5) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,5) = (/-1d0,  1d0,  0d0/)
       rot(3,:,5) = (/ 0d0,  0d0,  1d0/)
       ! 6 : x-y,x,z         => 6+ [ 0 0 1 ]        
       rot(1,:,6) = (/ 1d0, -1d0,  0d0/)
       rot(2,:,6) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,6) = (/ 0d0,  0d0,  1d0/)
       ! 7 : -x,-y,-z        => -1                  
       rot(1,:,7) = (/-1d0,  0d0,  0d0/)
       rot(2,:,7) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,7) = (/ 0d0,  0d0, -1d0/)
       ! 8 : y,-x+y,-z       => -3+ [ 0 0 1 ]       
       rot(1,:,8) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,8) = (/-1d0,  1d0,  0d0/)
       rot(3,:,8) = (/ 0d0,  0d0, -1d0/)
       ! 9 : x-y,x,-z        => -3- [ 0 0 1 ]       
       rot(1,:,9) = (/ 1d0, -1d0,  0d0/)
       rot(2,:,9) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,9) = (/ 0d0,  0d0, -1d0/)
       ! 10: x,y,-z          => m [ 0 0 1 ]         
       rot(1,:,10) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,10) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,10) = (/ 0d0,  0d0, -1d0/)
       ! 11: -y,x-y,-z       => -6- [ 0 0 1 ]       
       rot(1,:,11) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,11) = (/ 1d0, -1d0,  0d0/)
       rot(3,:,11) = (/ 0d0,  0d0, -1d0/)
       ! 12: -x+y,-x,-z      => -6+ [ 0 0 1 ]       
       rot(1,:,12) = (/-1d0,  1d0,  0d0/)
       rot(2,:,12) = (/-1d0,  0d0,  0d0/)
       rot(3,:,12) = (/ 0d0,  0d0, -1d0/)
    else if (equal(trim(str)//null,"6/mmm"//null).or.equal(trim(str)//null,"d6h"//null)) then
       nrot = 24
       ! 1 : x,y,z           => 1                   
       rot(1,:,1) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,1) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,1) = (/ 0d0,  0d0,  1d0/)
       ! 2 : -y,x-y,z        => 3+ [ 0 0 1 ]        
       rot(1,:,2) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,2) = (/ 1d0, -1d0,  0d0/)
       rot(3,:,2) = (/ 0d0,  0d0,  1d0/)
       ! 3 : -x+y,-x,z       => 3- [ 0 0 1 ]        
       rot(1,:,3) = (/-1d0,  1d0,  0d0/)
       rot(2,:,3) = (/-1d0,  0d0,  0d0/)
       rot(3,:,3) = (/ 0d0,  0d0,  1d0/)
       ! 4 : -x,-y,z         => 2 [ 0 0 1 ]         
       rot(1,:,4) = (/-1d0,  0d0,  0d0/)
       rot(2,:,4) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,4) = (/ 0d0,  0d0,  1d0/)
       ! 5 : y,-x+y,z        => 6- [ 0 0 1 ]        
       rot(1,:,5) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,5) = (/-1d0,  1d0,  0d0/)
       rot(3,:,5) = (/ 0d0,  0d0,  1d0/)
       ! 6 : x-y,x,z         => 6+ [ 0 0 1 ]        
       rot(1,:,6) = (/ 1d0, -1d0,  0d0/)
       rot(2,:,6) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,6) = (/ 0d0,  0d0,  1d0/)
       ! 7 : y,x,-z          => 2 [ 1 1 0 ]         
       rot(1,:,7) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,7) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,7) = (/ 0d0,  0d0, -1d0/)
       ! 8 : x-y,-y,-z       => 2 [ 1 0 0 ]         
       rot(1,:,8) = (/ 1d0, -1d0,  0d0/)
       rot(2,:,8) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,8) = (/ 0d0,  0d0, -1d0/)
       ! 9 : -x,-x+y,-z      => 2 [ 0 -1 0 ]        
       rot(1,:,9) = (/-1d0,  0d0,  0d0/)
       rot(2,:,9) = (/-1d0,  1d0,  0d0/)
       rot(3,:,9) = (/ 0d0,  0d0, -1d0/)
       ! 10: -y,-x,-z        => 2 [ 1 -1 0 ]        
       rot(1,:,10) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,10) = (/-1d0,  0d0,  0d0/)
       rot(3,:,10) = (/ 0d0,  0d0, -1d0/)
       ! 11: -x+y,y,-z       => 2 [ 0d05 1 0 ]       
       rot(1,:,11) = (/-1d0,  1d0,  0d0/)
       rot(2,:,11) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,11) = (/ 0d0,  0d0, -1d0/)
       ! 12: x,x-y,-z        => 2 [ 1 0d05 0 ]       
       rot(1,:,12) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,12) = (/ 1d0, -1d0,  0d0/)
       rot(3,:,12) = (/ 0d0,  0d0, -1d0/)
       ! 13: -x,-y,-z        => -1                  
       rot(1,:,13) = (/-1d0,  0d0,  0d0/)
       rot(2,:,13) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,13) = (/ 0d0,  0d0, -1d0/)
       ! 14: y,-x+y,-z       => -3+ [ 0 0 1 ]       
       rot(1,:,14) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,14) = (/-1d0,  1d0,  0d0/)
       rot(3,:,14) = (/ 0d0,  0d0, -1d0/)
       ! 15: x-y,x,-z        => -3- [ 0 0 1 ]       
       rot(1,:,15) = (/ 1d0, -1d0,  0d0/)
       rot(2,:,15) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,15) = (/ 0d0,  0d0, -1d0/)
       ! 16: x,y,-z          => m [ 0 0 1 ]         
       rot(1,:,16) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,16) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,16) = (/ 0d0,  0d0, -1d0/)
       ! 17: -y,x-y,-z       => -6- [ 0 0 1 ]       
       rot(1,:,17) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,17) = (/ 1d0, -1d0,  0d0/)
       rot(3,:,17) = (/ 0d0,  0d0, -1d0/)
       ! 18: -x+y,-x,-z      => -6+ [ 0 0 1 ]       
       rot(1,:,18) = (/-1d0,  1d0,  0d0/)
       rot(2,:,18) = (/-1d0,  0d0,  0d0/)
       rot(3,:,18) = (/ 0d0,  0d0, -1d0/)
       ! 19: -y,-x,z         => m [ 1 1 0 ]         
       rot(1,:,19) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,19) = (/-1d0,  0d0,  0d0/)
       rot(3,:,19) = (/ 0d0,  0d0,  1d0/)
       ! 20: -x+y,y,z        => m [ 1 0 0 ]         
       rot(1,:,20) = (/-1d0,  1d0,  0d0/)
       rot(2,:,20) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,20) = (/ 0d0,  0d0,  1d0/)
       ! 21: x,x-y,z         => m [ 0 -1 0 ]        
       rot(1,:,21) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,21) = (/ 1d0, -1d0,  0d0/)
       rot(3,:,21) = (/ 0d0,  0d0,  1d0/)
       ! 22: y,x,z           => m [ 1 -1 0 ]        
       rot(1,:,22) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,22) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,22) = (/ 0d0,  0d0,  1d0/)
       ! 23: x-y,-y,z        => m [ 0d05 1 0 ]       
       rot(1,:,23) = (/ 1d0, -1d0,  0d0/)
       rot(2,:,23) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,23) = (/ 0d0,  0d0,  1d0/)
       ! 24: -x,-x+y,z       => m [ 1 0d05 0 ]       
       rot(1,:,24) = (/-1d0,  0d0,  0d0/)
       rot(2,:,24) = (/-1d0,  1d0,  0d0/)
       rot(3,:,24) = (/ 0d0,  0d0,  1d0/)
    else if (equal(trim(str)//null,"m-3"//null).or.equal(trim(str)//null,"th"//null)) then
       nrot = 24
       ! 1 : x,y,z           => 1                   
       rot(1,:,1) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,1) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,1) = (/ 0d0,  0d0,  1d0/)
       ! 2 : -x,-y,z         => 2 [ 0 0 1 ]         
       rot(1,:,2) = (/-1d0,  0d0,  0d0/)
       rot(2,:,2) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,2) = (/ 0d0,  0d0,  1d0/)
       ! 3 : -x,y,-z         => 2 [ 0 1 0 ]         
       rot(1,:,3) = (/-1d0,  0d0,  0d0/)
       rot(2,:,3) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,3) = (/ 0d0,  0d0, -1d0/)
       ! 4 : x,-y,-z         => 2 [ 1 0 0 ]         
       rot(1,:,4) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,4) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,4) = (/ 0d0,  0d0, -1d0/)
       ! 5 : z,x,y           => 3+ [ 1 1 1 ]        
       rot(1,:,5) = (/ 0d0,  0d0,  1d0/)
       rot(2,:,5) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,5) = (/ 0d0,  1d0,  0d0/)
       ! 6 : z,-x,-y         => 3- [ 1 -1 1 ]       
       rot(1,:,6) = (/ 0d0,  0d0,  1d0/)
       rot(2,:,6) = (/-1d0,  0d0,  0d0/)
       rot(3,:,6) = (/ 0d0, -1d0,  0d0/)
       ! 7 : -z,-x,y         => 3+ [ 1 -1 -1 ]      
       rot(1,:,7) = (/ 0d0,  0d0, -1d0/)
       rot(2,:,7) = (/-1d0,  0d0,  0d0/)
       rot(3,:,7) = (/ 0d0,  1d0,  0d0/)
       ! 8 : -z,x,-y         => 3- [ 1 1 -1 ]       
       rot(1,:,8) = (/ 0d0,  0d0, -1d0/)
       rot(2,:,8) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,8) = (/ 0d0, -1d0,  0d0/)
       ! 9 : y,z,x           => 3- [ 1 1 1 ]        
       rot(1,:,9) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,9) = (/ 0d0,  0d0,  1d0/)
       rot(3,:,9) = (/ 1d0,  0d0,  0d0/)
       ! 10: -y,z,-x         => 3- [ 1 -1 -1 ]      
       rot(1,:,10) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,10) = (/ 0d0,  0d0,  1d0/)
       rot(3,:,10) = (/-1d0,  0d0,  0d0/)
       ! 11: y,-z,-x         => 3+ [ 1 1 -1 ]       
       rot(1,:,11) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,11) = (/ 0d0,  0d0, -1d0/)
       rot(3,:,11) = (/-1d0,  0d0,  0d0/)
       ! 12: -y,-z,x         => 3+ [ 1 -1 1 ]       
       rot(1,:,12) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,12) = (/ 0d0,  0d0, -1d0/)
       rot(3,:,12) = (/ 1d0,  0d0,  0d0/)
       ! 13: -x,-y,-z        => -1                  
       rot(1,:,13) = (/-1d0,  0d0,  0d0/)
       rot(2,:,13) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,13) = (/ 0d0,  0d0, -1d0/)
       ! 14: x,y,-z          => m [ 0 0 1 ]         
       rot(1,:,14) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,14) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,14) = (/ 0d0,  0d0, -1d0/)
       ! 15: x,-y,z          => m [ 0 1 0 ]         
       rot(1,:,15) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,15) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,15) = (/ 0d0,  0d0,  1d0/)
       ! 16: -x,y,z          => m [ 1 0 0 ]         
       rot(1,:,16) = (/-1d0,  0d0,  0d0/)
       rot(2,:,16) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,16) = (/ 0d0,  0d0,  1d0/)
       ! 17: -z,-x,-y        => -3+ [ 1 1 1 ]       
       rot(1,:,17) = (/ 0d0,  0d0, -1d0/)
       rot(2,:,17) = (/-1d0,  0d0,  0d0/)
       rot(3,:,17) = (/ 0d0, -1d0,  0d0/)
       ! 18: -z,x,y          => -3- [ 1 -1 1 ]      
       rot(1,:,18) = (/ 0d0,  0d0, -1d0/)
       rot(2,:,18) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,18) = (/ 0d0,  1d0,  0d0/)
       ! 19: z,x,-y          => -3+ [ 1 -1 -1 ]     
       rot(1,:,19) = (/ 0d0,  0d0,  1d0/)
       rot(2,:,19) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,19) = (/ 0d0, -1d0,  0d0/)
       ! 20: z,-x,y          => -3- [ 1 1 -1 ]      
       rot(1,:,20) = (/ 0d0,  0d0,  1d0/)
       rot(2,:,20) = (/-1d0,  0d0,  0d0/)
       rot(3,:,20) = (/ 0d0,  1d0,  0d0/)
       ! 21: -y,-z,-x        => -3- [ 1 1 1 ]       
       rot(1,:,21) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,21) = (/ 0d0,  0d0, -1d0/)
       rot(3,:,21) = (/-1d0,  0d0,  0d0/)
       ! 22: y,-z,x          => -3- [ 1 -1 -1 ]     
       rot(1,:,22) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,22) = (/ 0d0,  0d0, -1d0/)
       rot(3,:,22) = (/ 1d0,  0d0,  0d0/)
       ! 23: -y,z,x          => -3+ [ 1 1 -1 ]      
       rot(1,:,23) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,23) = (/ 0d0,  0d0,  1d0/)
       rot(3,:,23) = (/ 1d0,  0d0,  0d0/)
       ! 24: y,z,-x          => -3+ [ 1 -1 1 ]      
       rot(1,:,24) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,24) = (/ 0d0,  0d0,  1d0/)
       rot(3,:,24) = (/-1d0,  0d0,  0d0/)
    else if (equal(trim(str)//null,"m-3m"//null).or.equal(trim(str)//null,"oh"//null)) then
       nrot = 48
       ! 1 : x,y,z           => 1                   
       rot(1,:,1) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,1) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,1) = (/ 0d0,  0d0,  1d0/)
       ! 2 : -x,-y,z         => 2 [ 0 0 1 ]         
       rot(1,:,2) = (/-1d0,  0d0,  0d0/)
       rot(2,:,2) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,2) = (/ 0d0,  0d0,  1d0/)
       ! 3 : x,-y,-z         => 2 [ 1 0 0 ]         
       rot(1,:,3) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,3) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,3) = (/ 0d0,  0d0, -1d0/)
       ! 4 : -x,y,-z         => 2 [ 0 1 0 ]         
       rot(1,:,4) = (/-1d0,  0d0,  0d0/)
       rot(2,:,4) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,4) = (/ 0d0,  0d0, -1d0/)
       ! 5 : z,x,y           => 3+ [ 1 1 1 ]        
       rot(1,:,5) = (/ 0d0,  0d0,  1d0/)
       rot(2,:,5) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,5) = (/ 0d0,  1d0,  0d0/)
       ! 6 : z,-x,-y         => 3- [ 1 -1 1 ]       
       rot(1,:,6) = (/ 0d0,  0d0,  1d0/)
       rot(2,:,6) = (/-1d0,  0d0,  0d0/)
       rot(3,:,6) = (/ 0d0, -1d0,  0d0/)
       ! 7 : -z,x,-y         => 3- [ 1 1 -1 ]       
       rot(1,:,7) = (/ 0d0,  0d0, -1d0/)
       rot(2,:,7) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,7) = (/ 0d0, -1d0,  0d0/)
       ! 8 : -z,-x,y         => 3+ [ 1 -1 -1 ]      
       rot(1,:,8) = (/ 0d0,  0d0, -1d0/)
       rot(2,:,8) = (/-1d0,  0d0,  0d0/)
       rot(3,:,8) = (/ 0d0,  1d0,  0d0/)
       ! 9 : y,z,x           => 3- [ 1 1 1 ]        
       rot(1,:,9) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,9) = (/ 0d0,  0d0,  1d0/)
       rot(3,:,9) = (/ 1d0,  0d0,  0d0/)
       ! 10: -y,z,-x         => 3- [ 1 -1 -1 ]      
       rot(1,:,10) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,10) = (/ 0d0,  0d0,  1d0/)
       rot(3,:,10) = (/-1d0,  0d0,  0d0/)
       ! 11: -y,-z,x         => 3+ [ 1 -1 1 ]       
       rot(1,:,11) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,11) = (/ 0d0,  0d0, -1d0/)
       rot(3,:,11) = (/ 1d0,  0d0,  0d0/)
       ! 12: y,-z,-x         => 3+ [ 1 1 -1 ]       
       rot(1,:,12) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,12) = (/ 0d0,  0d0, -1d0/)
       rot(3,:,12) = (/-1d0,  0d0,  0d0/)
       ! 13: y,x,-z          => 2 [ 1 1 0 ]         
       rot(1,:,13) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,13) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,13) = (/ 0d0,  0d0, -1d0/)
       ! 14: -y,-x,-z        => 2 [ 1 -1 0 ]        
       rot(1,:,14) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,14) = (/-1d0,  0d0,  0d0/)
       rot(3,:,14) = (/ 0d0,  0d0, -1d0/)
       ! 15: -y,x,z          => 4+ [ 0 0 1 ]        
       rot(1,:,15) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,15) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,15) = (/ 0d0,  0d0,  1d0/)
       ! 16: y,-x,z          => 4- [ 0 0 1 ]        
       rot(1,:,16) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,16) = (/-1d0,  0d0,  0d0/)
       rot(3,:,16) = (/ 0d0,  0d0,  1d0/)
       ! 17: x,z,-y          => 4- [ 1 0 0 ]        
       rot(1,:,17) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,17) = (/ 0d0,  0d0,  1d0/)
       rot(3,:,17) = (/ 0d0, -1d0,  0d0/)
       ! 18: -x,z,y          => 2 [ 0 1 1 ]         
       rot(1,:,18) = (/-1d0,  0d0,  0d0/)
       rot(2,:,18) = (/ 0d0,  0d0,  1d0/)
       rot(3,:,18) = (/ 0d0,  1d0,  0d0/)
       ! 19: x,-z,y          => 4+ [ 1 0 0 ]        
       rot(1,:,19) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,19) = (/ 0d0,  0d0, -1d0/)
       rot(3,:,19) = (/ 0d0,  1d0,  0d0/)
       ! 20: -x,-z,-y        => 2 [ 0 1 -1 ]        
       rot(1,:,20) = (/-1d0,  0d0,  0d0/)
       rot(2,:,20) = (/ 0d0,  0d0, -1d0/)
       rot(3,:,20) = (/ 0d0, -1d0,  0d0/)
       ! 21: z,y,-x          => 4+ [ 0 1 0 ]        
       rot(1,:,21) = (/ 0d0,  0d0,  1d0/)
       rot(2,:,21) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,21) = (/-1d0,  0d0,  0d0/)
       ! 22: z,-y,x          => 2 [ 1 0 1 ]         
       rot(1,:,22) = (/ 0d0,  0d0,  1d0/)
       rot(2,:,22) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,22) = (/ 1d0,  0d0,  0d0/)
       ! 23: -z,-y,-x        => 2 [ 1 0 -1 ]        
       rot(1,:,23) = (/ 0d0,  0d0, -1d0/)
       rot(2,:,23) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,23) = (/-1d0,  0d0,  0d0/)
       ! 24: -z,y,x          => 4- [ 0 1 0 ]        
       rot(1,:,24) = (/ 0d0,  0d0, -1d0/)
       rot(2,:,24) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,24) = (/ 1d0,  0d0,  0d0/)
       ! 25: -x,-y,-z        => -1                  
       rot(1,:,25) = (/-1d0,  0d0,  0d0/)
       rot(2,:,25) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,25) = (/ 0d0,  0d0, -1d0/)
       ! 26: x,y,-z          => m [ 0 0 1 ]         
       rot(1,:,26) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,26) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,26) = (/ 0d0,  0d0, -1d0/)
       ! 27: -x,y,z          => m [ 1 0 0 ]         
       rot(1,:,27) = (/-1d0,  0d0,  0d0/)
       rot(2,:,27) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,27) = (/ 0d0,  0d0,  1d0/)
       ! 28: x,-y,z          => m [ 0 1 0 ]         
       rot(1,:,28) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,28) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,28) = (/ 0d0,  0d0,  1d0/)
       ! 29: -z,-x,-y        => -3+ [ 1 1 1 ]       
       rot(1,:,29) = (/ 0d0,  0d0, -1d0/)
       rot(2,:,29) = (/-1d0,  0d0,  0d0/)
       rot(3,:,29) = (/ 0d0, -1d0,  0d0/)
       ! 30: -z,x,y          => -3- [ 1 -1 1 ]      
       rot(1,:,30) = (/ 0d0,  0d0, -1d0/)
       rot(2,:,30) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,30) = (/ 0d0,  1d0,  0d0/)
       ! 31: z,-x,y          => -3- [ 1 1 -1 ]      
       rot(1,:,31) = (/ 0d0,  0d0,  1d0/)
       rot(2,:,31) = (/-1d0,  0d0,  0d0/)
       rot(3,:,31) = (/ 0d0,  1d0,  0d0/)
       ! 32: z,x,-y          => -3+ [ 1 -1 -1 ]     
       rot(1,:,32) = (/ 0d0,  0d0,  1d0/)
       rot(2,:,32) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,32) = (/ 0d0, -1d0,  0d0/)
       ! 33: -y,-z,-x        => -3- [ 1 1 1 ]       
       rot(1,:,33) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,33) = (/ 0d0,  0d0, -1d0/)
       rot(3,:,33) = (/-1d0,  0d0,  0d0/)
       ! 34: y,-z,x          => -3- [ 1 -1 -1 ]     
       rot(1,:,34) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,34) = (/ 0d0,  0d0, -1d0/)
       rot(3,:,34) = (/ 1d0,  0d0,  0d0/)
       ! 35: y,z,-x          => -3+ [ 1 -1 1 ]      
       rot(1,:,35) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,35) = (/ 0d0,  0d0,  1d0/)
       rot(3,:,35) = (/-1d0,  0d0,  0d0/)
       ! 36: -y,z,x          => -3+ [ 1 1 -1 ]      
       rot(1,:,36) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,36) = (/ 0d0,  0d0,  1d0/)
       rot(3,:,36) = (/ 1d0,  0d0,  0d0/)
       ! 37: -y,-x,z         => m [ 1 1 0 ]         
       rot(1,:,37) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,37) = (/-1d0,  0d0,  0d0/)
       rot(3,:,37) = (/ 0d0,  0d0,  1d0/)
       ! 38: y,x,z           => m [ 1 -1 0 ]        
       rot(1,:,38) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,38) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,38) = (/ 0d0,  0d0,  1d0/)
       ! 39: y,-x,-z         => -4+ [ 0 0 1 ]       
       rot(1,:,39) = (/ 0d0,  1d0,  0d0/)
       rot(2,:,39) = (/-1d0,  0d0,  0d0/)
       rot(3,:,39) = (/ 0d0,  0d0, -1d0/)
       ! 40: -y,x,-z         => -4- [ 0 0 1 ]       
       rot(1,:,40) = (/ 0d0, -1d0,  0d0/)
       rot(2,:,40) = (/ 1d0,  0d0,  0d0/)
       rot(3,:,40) = (/ 0d0,  0d0, -1d0/)
       ! 41: -x,-z,y         => -4- [ 1 0 0 ]       
       rot(1,:,41) = (/-1d0,  0d0,  0d0/)
       rot(2,:,41) = (/ 0d0,  0d0, -1d0/)
       rot(3,:,41) = (/ 0d0,  1d0,  0d0/)
       ! 42: x,-z,-y         => m [ 0 1 1 ]         
       rot(1,:,42) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,42) = (/ 0d0,  0d0, -1d0/)
       rot(3,:,42) = (/ 0d0, -1d0,  0d0/)
       ! 43: -x,z,-y         => -4+ [ 1 0 0 ]       
       rot(1,:,43) = (/-1d0,  0d0,  0d0/)
       rot(2,:,43) = (/ 0d0,  0d0,  1d0/)
       rot(3,:,43) = (/ 0d0, -1d0,  0d0/)
       ! 44: x,z,y           => m [ 0 1 -1 ]        
       rot(1,:,44) = (/ 1d0,  0d0,  0d0/)
       rot(2,:,44) = (/ 0d0,  0d0,  1d0/)
       rot(3,:,44) = (/ 0d0,  1d0,  0d0/)
       ! 45: -z,-y,x         => -4+ [ 0 1 0 ]       
       rot(1,:,45) = (/ 0d0,  0d0, -1d0/)
       rot(2,:,45) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,45) = (/ 1d0,  0d0,  0d0/)
       ! 46: -z,y,-x         => m [ 1 0 1 ]         
       rot(1,:,46) = (/ 0d0,  0d0, -1d0/)
       rot(2,:,46) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,46) = (/-1d0,  0d0,  0d0/)
       ! 47: z,y,x           => m [ 1 0 -1 ]        
       rot(1,:,47) = (/ 0d0,  0d0,  1d0/)
       rot(2,:,47) = (/ 0d0,  1d0,  0d0/)
       rot(3,:,47) = (/ 1d0,  0d0,  0d0/)
       ! 48: z,-y,-x         => -4- [ 0 1 0 ]       
       rot(1,:,48) = (/ 0d0,  0d0,  1d0/)
       rot(2,:,48) = (/ 0d0, -1d0,  0d0/)
       rot(3,:,48) = (/-1d0,  0d0,  0d0/)
    else
       write (uout,'("Input group : ",A)') trim(str)
       call error('laue','unknown group',faterr)
    end if
    
  end subroutine laue

end module tools

