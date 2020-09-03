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

module gnuplot_templates
  implicit none

contains

  ! Open the gnuplot script root.gnu. If doout, set the output file
  ! in the gnuplot script. Returns the logical unit for the open
  ! file.
  function opengnu(root,noout)
    use tools, only: fopen
    use param, only: mline, iowrite, mcols, null, gplt_rgb, gplt_sym
    integer :: opengnu
    character*(*), intent(in) :: root
    logical, intent(in), optional :: noout

    integer :: lu
    character*(mline) :: file
    integer :: i
    logical :: doout

    doout = .true.
    if (present(noout)) doout = .not.noout

    file = trim(adjustl(root)) // ".gnu" // null
    lu = fopen(lu,file,iowrite)
    write (lu,'("set terminal postscript eps color enhanced ""Helvetica"" 14")')
    if (doout) then
       write (lu,'("set output """,A,".eps""")') root
    end if
    write (lu,*)
    do i = 1, mcols
       write (lu,'("set style line ",I2," lt 1 lc rgb """,A7,""" pt ",I2," ps 0.75")')&
          i, gplt_rgb(i), gplt_sym(i)
    end do
    write (lu,'("set style increment user ")')
    write (lu,*)
    opengnu = lu

  end function opengnu

  ! Close the gnuplot file root.gnu in logical unit lu. If noout,
  ! generate the pdf file from the eps file at the end.
  subroutine closegnu(root,lu,noout)
    use tools, only: fclose
    character*(*), intent(in) :: root
    integer, intent(in) :: lu
    logical, intent(in), optional :: noout

    logical :: doout

    doout = .true.
    if (present(noout)) doout = .not.noout

    if (doout) then
       write (lu,'("!epstopdf ",A,".eps")') root
       write (lu,'("!pdfcrop ",A,".pdf")') root
       write (lu,'("!mv ",A,"-crop.pdf ",A,".pdf")') root, root
       write (lu,'("!rm ",A,".eps")') root
    end if

    call fclose(lu)

  end subroutine closegnu

  ! Write the gnuplot file with all thermodynamic properties as a
  ! function of temperature in file root.gnu.
  subroutine gen_allgnu_t(root)
    use varbas, only: ph, nph, mpropout, nts, tlist, tm_static, writelevel, propname
    use tools, only: leng
    use param, only: mline, fileroot, uout
    character*(*), intent(in) :: root

    integer :: lu
    integer :: i, ip
    character*(mline) :: oroot
    integer :: ndo, imask(nph)

    integer, parameter :: tlistcrit = 4

    if (nts < tlistcrit .or. writelevel < 2) return

    ndo = 0
    imask = 0
    do i = 1, nph
       if (ph(i)%tmodel == tm_static) cycle
       ndo = ndo + 1
       imask(ndo) = i
    end do
    if (ndo == 0) return
          
    lu = opengnu(root,.true.)

    write (uout,'("  Writing file: ",A/)') trim(root)//".gnu"

    write (lu,'("!sed ''/^ *$/d'' ",A," | awk ''$1+0==0'' | awk ''/# Phase/{print """"; print """"} {print}''> temp.dat")') &
       trim(fileroot)//".eos"
    write (lu,*)
    write (lu,'("set xrange [0:",F20.1,"]")') tlist(nts)
    write (lu,'("set xlabel ""T (K)""")')
    write (lu,*)

    do ip = 1, mpropout
       if (ip == 1 .or. ip == 2) cycle ! skip temperature and pressure
       write (oroot,'(A,"_t_",I2.2)') trim(fileroot), ip
       write (lu,'("set output ''",A,"''")') trim(oroot) // ".eps"
       write (lu,'("set ylabel """,A,"""")') trim(adjustl(propname(ip)))
       write (lu,'("plot \")')
       do i = 1, ndo-1
          write (lu,'("     ''temp.dat'' u 2:",I2," index ",I3," w lines title ''",A10,"''  ,\")') &
             ip, i-1, ph(imask(i))%name(1:leng(ph(imask(i))%name))
       end do
       write (lu,'("     ''temp.dat'' u 2:",I2," index ",I3," w lines title ''",A10,"''")') &
          ip, ndo-1, ph(imask(ndo))%name(1:leng(ph(imask(ndo))%name))
       write (lu,'("!epstopdf ",A)') trim(oroot) // ".eps"
       write (lu,'("!pdfcrop ",A)') trim(oroot) // ".pdf"
       write (lu,'("!mv ",A)') trim(oroot)//"-crop.pdf "//&
          adjustl(trim(oroot))//".pdf"
       write (lu,'("!rm ",A)') trim(oroot) // ".eps"
       write (lu,*)
    end do

    write (lu,'("!rm temp.dat")')

    call closegnu(root,lu,.true.)

  end subroutine gen_allgnu_t

  ! Write the gnuplot file with all thermodynamic properties as a
  ! function of pressure in file root.gnu.
  subroutine gen_allgnu_p(root)
    use param, only: mline, fileroot, uout
    use tools, only: leng
    use varbas, only: ph, nph, mpropout, nps, plist, nts, tm_static, writelevel, propname, tlist
    character*(*), intent(in) :: root

    integer, parameter :: ntemp = 5 

    integer :: lu
    integer :: i, j, ip
    character*(mline) :: oroot, titleplot
    integer :: outidx(ntemp)
    real*8 :: step
    integer :: idotemp
    integer :: ndo, imask(nph)

    integer, parameter :: plistcrit = 4

    if (nps < plistcrit .or. writelevel < 2) return

    ndo = 0
    imask = 0
    do i = 1, nph
       if (ph(i)%tmodel == tm_static) cycle
       ndo = ndo + 1
       imask(ndo) = i
    end do
    if (ndo == 0) return
          
    lu = opengnu(root,.true.)

    write (uout,'("  Writing file: ",A/)') trim(root)//".gnu"

    step = max(real(nts,8)/real(ntemp-1,8),1d0)
    do i = 1, ntemp
       outidx(i) = max(min(nint((i-1) * step + 1),nts),1) - 1
    end do
    idotemp = min(ntemp,nts)

    write (lu,'("!cp ",A," temp.dat")') trim(fileroot)//".eos"
    write (lu,*)
    write (lu,'("set xrange [0:",F20.1,"]")') plist(nps)
    write (lu,'("set xlabel ""p (GPa)""")')
    write (lu,*)

    do ip = 1, mpropout
       if (ip == 1 .or. ip == 2) cycle ! skip temperature and pressure
       write (oroot,'(A,"_p_",I2.2)') trim(fileroot), ip
       write (lu,'("set output ''",A,"''")') trim(oroot) // ".eps"
       write (lu,'("set ylabel """,A,"""")') trim(adjustl(propname(ip)))
       write (lu,'("plot \")')
       a: do i = 1, ndo
          do j = 1, idotemp
             write (titleplot,'(A10,",",F8.2,"K")') ph(imask(i))%name(1:leng(ph(imask(i))%name)), tlist(outidx(j)+1)
             if (i == ndo .and. j == idotemp) exit a
             write (lu,'(" ''temp.dat'' u 1:",I2," index ",I3," w lines ls ",I2," title ''",A20,"''  ,\")') &
                ip, (i-1)*nts+outidx(j), j, trim(titleplot)
          end do
       end do a
       write (lu,'(" ''temp.dat'' u 1:",I2," index ",I3," w lines ls ",I2," title ''",A20,"''")') &
          ip, (ndo-1)*nts+outidx(idotemp), idotemp, trim(titleplot)
       write (lu,'("!epstopdf ",A)') trim(oroot) // ".eps"
       write (lu,'("!pdfcrop ",A)') trim(oroot) // ".pdf"
       write (lu,'("!mv ",A)') trim(oroot)//"-crop.pdf "//&
          adjustl(trim(oroot))//".pdf"
       write (lu,'("!rm ",A)') trim(oroot) // ".eps"
       write (lu,*)
    end do

    write (lu,'("!rm temp.dat")')

    call closegnu(root,lu,.true.)
  end subroutine gen_allgnu_p

end module gnuplot_templates
