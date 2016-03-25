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

module topcalc
  use varbas
  use tools
  use param
  implicit none

  public
  integer, parameter, private :: mtrans = 50

  ! interpolation
  integer, parameter :: mxint = 30
  real*8, allocatable :: fint(:)
  integer, allocatable :: iint(:)
  integer :: nxint
  integer :: interp_input

  ! zero-temperature cutoff
  real*8, parameter :: tsmall_el = 0.1d0

  ! pack of volume grid fits
  type fitpack
     integer :: nepol = 0    ! static energy
     integer :: napol = 0    ! total helmholtz free energy
     integer :: nspol = 0    ! total entropy
     ! fit modes
     integer :: emode = 0    
     integer :: amode = 0    
     integer :: smode = 0 
     ! coefficients
     real*8 :: epol(0:mmpar) = 0d0
     real*8 :: apol(0:mmpar) = 0d0
     real*8 :: spol(0:mmpar) = 0d0
  end type fitpack
  type(fitpack), save :: ft_null

contains

  subroutine topcalc_init()

    nxint = 0
    interp_input = 0

  end subroutine topcalc_init

  subroutine popinput(sdate,fileout)
    use debye
    use eos
    use evfunc

    character*(mline), intent(in) :: sdate,fileout

    integer :: i

    write (uout,'("* Input ")') 
    write (uout,'("  Title: ",A)') trim(adjustl(title(1:leng(title)-1)))
    write (uout,'("  Output file (lu=",I2,"): ",A)') uout, &
       fileout(1:leng(fileout))
    write (uout,'("  Units: output is in atomic units, except where noted.")')
    write (uout,'("  Number of atoms per primitive cell: ",I3)') vfree
    write (uout,'("  Molecular mass (amu): ",F17.8)') mm/amu2au
    write (uout,'("  Infinite V energy (hy): ",1p,E20.12)') einf
    write (uout,'("  Number of phases: ",I3)') nph
    do i = 1, nph
       write (uout,'("  Phase ",I3,": ",A)') i, ph(i)%name(1:leng(ph(i)%name))
    end do
    select case(pfit_mode)
    case(pfit_gauss)
       write (uout,'("  Polynomial fit mode: gibbs1 (gauss)")')
    case(pfit_slatec)
       write (uout,'("  Polynomial fit mode: slatec")')
    end select
    select case(pweigh_mode)
    case(pweigh_gibbs2)
       write (uout,'("  Polynomial weight mode: gibbs2")')
       write (uout,'("  Max. polynomial degree: ",I2)') mpar
    case(pweigh_gibbs1)
       write (uout,'("  Polynomial weight mode: gibbs1")')
       write (uout,'("  Max. polynomial degree: ",I2)') mpar
    case(pweigh_slatec)
       write (uout,'("  Polynomial weight mode: slatec")')
    end select
    write (uout,*)

  end subroutine popinput

  subroutine popenergyfit()
    use gnuplot_templates
    use evfunc

    integer, parameter :: igrid = 1001

    integer :: lu
    integer :: i, j
    real*8 :: xdum, f0, f1, fac
    character*5 :: starts
    character*3 :: ends
    integer :: nact, iact
    character*(mline_fmt) :: fm
    real*8 :: vmin, vmax, emin, emax

    if (.not.doefit .or. writelevel < 2) return

    ! input and fitted energies
    write (uout,'("* Input and fitted static energy")')
    write (uout,'("  Writing file: ",A/)') trim(fileroot)//".efit"
    lu = fopen(lu,trim(fileroot)//".efit"//null,iowrite)

    vmin = 1d30
    vmax = -1d30
    emin = 1d30
    emax = -1d30
    do i = 1, nph
       write (lu,'("# Phase ",A)') trim(adjustl(ph(i)%name(1:leng(ph(i)%name))))
       fm = format_string_header((/1,ifmt_v,ifmt_eprec,ifmt_eprec,ifmt_p/),(/1,9,9,9,10/))
       write (lu,fm) "#", "V(bohr^3)", "E_inp(Ha)", "E_fit(Ha)","p_fit(GPa)"
       fm = format_string((/ifmt_v,ifmt_eprec,ifmt_eprec,ifmt_p/),1)
       do j = 1, ph(i)%nv
          f0 = fv0(ph(i)%fit_mode, ph(i)%v(j), ph(i)%npol, ph(i)%cpol)
          f1 = fv1(ph(i)%fit_mode, ph(i)%v(j), ph(i)%npol, ph(i)%cpol)
          write (lu,fm) ph(i)%v(j),ph(i)%e(j),f0,-f1*au2gpa
       end do
       write (lu,'(/)')
       vmin = min(vmin,ph(i)%v(1))
       vmax = max(vmax,ph(i)%v(ph(i)%nv))
       emin = min(emin,minval(ph(i)%e))
       emax = max(emax,maxval(ph(i)%e))
    end do
    call fclose(lu)

    ! auxiliary file
    write (uout,'("  Writing file : ",A/)') trim(fileroot)//"_efit.aux"
    lu = fopen(lu,trim(fileroot)//"_efit.aux"//null,iowrite)
    write (lu,'("# E(V) plot, auxiliary file. ")') 
    write (lu,'("# Contains fitted E(V) calculated on a finer grid. ")') 
    do i = 1, nph
       write (lu,'("# Phase ",A)') trim(adjustl(ph(i)%name(1:leng(ph(i)%name))))
       fm = format_string_header((/1,ifmt_v,ifmt_eprec/),(/1,9,9/))
       write (lu,fm) "#", "V(bohr^3)", "E_fit(Ha)"
       fm = format_string((/ifmt_v,ifmt_eprec/),1)
       do j = igrid, 0, -1
          fac = real(j,8) / real(igrid,8)
          xdum = ph(i)%v(1) * fac + ph(i)%v(ph(i)%nv) * (1d0 - fac)
          f0 = fv0 (ph(i)%fit_mode, xdum, ph(i)%npol, ph(i)%cpol)
          write (lu,fm) xdum, f0
       end do
       write (lu,'(/)')
    end do
    call fclose(lu)

    ! write the gnu file
    write (uout,'("  Writing file : ",A/)') trim(fileroot)//"_efit.gnu"
    lu = opengnu(trim(fileroot)//"_efit")
    write (lu,'("set xrange ["F20.4":"F20.4"]")') vmin, vmax
    write (lu,'("set yrange ["F20.4":"F20.4"]")') emin, emax
    write (lu,'("set ylabel ""Energy (Ha)""")')
    write (lu,'("set xlabel ""Volume (bohr^3)""")')
    ends = ",\ "
    iact = 0
    do i = 1, nph
       iact = iact + 1
       if (iact == nph) then
          ends = "   "
       else if (iact == 1) then
          ends = ",\ "
       end if
       if (iact == 1) then
          starts = "plot "
       else if (iact == nph) then
          starts = "     "
       else
          starts = "     "
       end if
       write (lu,'(A5,"''",A,"'' u 1:2 index ",I2," w points ls ",I2," notitle,\")') &
          starts, trim(fileroot)//".efit", i-1, i
       write (lu,'("     ''",A,"'' u 1:2 index ",I2," w lines ls ",I2," title ''",A10,"'' ",A2)') &
          trim(fileroot)//"_efit.aux",i-1, i, ph(i)%name(1:leng(ph(i)%name)), ends
    end do
    call closegnu(trim(fileroot)//"_efit",lu)

  end subroutine popenergyfit

  subroutine plotdh(ga)
    use gnuplot_templates

    real*8, intent(in) :: ga(:,:)

    integer :: i, j, k
    character*5 :: starts
    character*3 :: ends
    integer :: lu
    integer :: imasknph, imaskph(nph)

    if (.not.doplotdh .or. writelevel < 2) return
    if (n_not_pv_min() < 2) then
       if (nph >= 2) then
          call error('plotdh','Not enough phases with p = 0 minimum to plot enthalpies.',warning)
       end if
       return
    end if

    imasknph = 0
    do i = 1, nph
       imasknph = imasknph + 1
       imaskph(imasknph) = i
    end do

    ! write the auxiliary file
    write (uout,'("* Plotting static DeltaH"/)')
    write (uout,'("  Writing file : ",A/)') trim(fileroot) // "_dH.aux"
    lu = fopen(lu,trim(fileroot) // "_dH.aux"//null,iowrite)
    write (lu,'("# H(p) plot, auxiliary file. ")') 
    write (lu,'("# Contains dH(p) (Ha) calculated on input pressures. ")') 
    write (lu,'("# p(GPa)",999(2X,A13,4X))') &
       (trim(adjustl(ph(imaskph(k))%name(1:leng(ph(imaskph(k))%name)))),k=2,imasknph)
    do j = 1, nps
       write (lu,'(F12.4,1X,1p,999(E20.10,2X))') &
          plist(j), (ga(j,imaskph(k))-ga(j,imaskph(1)),k=2,imasknph)
    end do
    write (lu,*)
    call fclose(lu)
    
    ! write the gnu file
    write (uout,'("  Writing file : ",A/)') trim(fileroot) // "_dH.gnu"
    lu = opengnu(trim(fileroot) // "_dH")
    write (lu,'("set ylabel ""{/Symbol D}H (Ha)""")')
    write (lu,'("set xlabel ""p (GPa)""")')
    write (lu,'("set xzeroaxis")')
    starts = "plot "
    ends = ",\ "
    do i = 1, imasknph
       if (i == 1) cycle
       if (i == imasknph) ends = "   "
       if (i /= 2) starts = "     "

       if (i < 9) then
          write (lu,'(A5,"''",A,"'' u 1:",I1," w points ls ",I2," title ''",A10,"'' ",A2)') &
             starts, trim(fileroot)//"_dH.aux", i, imaskph(i), ph(i)%name(1:leng(ph(i)%name)), ends
       else
          write (lu,'(A5,"''",A,"'' u 1:",I2," w points ls ",I2," title ''",A10,"'' ",A2)') &
             starts, trim(fileroot)//"_dH.aux", i, imaskph(i), ph(i)%name(1:leng(ph(i)%name)), ends
       end if
    end do
    call closegnu(trim(fileroot) // "_dH",lu)

  end subroutine plotdh

  subroutine static_transp(ga)

    real*8, intent(in) :: ga(nps,nph)

    integer :: i, j, iact, idmin, iold
    real*8 :: ptrans(mtrans), gmin, pnew, pold
    integer :: i1trans(mtrans), i2trans(mtrans)
    integer :: ntrans

    if (.not.dotrans) return
    if (n_not_pv_min() < 2) then
       if (nph >= 2) then
          call error('static_transp','Not enough phases with p = 0 minimum to find transition pressures.',warning)
       end if
       return
    end if

    i1trans = 0
    i2trans = 0
    ptrans = 0d0
    ntrans = 0
    write (uout,'("* Static transition pressures (linear interpolation)")')
    write (uout,'("#        Pressure range (GPa)      Stable phase")')
    do j = 1, nps
       idmin = -1
       gmin = 1d30
       do i = 1, nph
          if (ga(j,i) < gmin) then
             gmin = ga(j,i)
             idmin = i
          end if
       end do

       if (j == 1) then
          iold = idmin
          pold = 0d0
       else
          if (idmin /= iold) then
             pnew = plist(j-1) - (ga(j-1,idmin)-ga(j-1,iold)) * (plist(j) - plist(j-1)) /&
                (ga(j,idmin)-ga(j,iold) - (ga(j-1,idmin)-ga(j-1,iold)))
             !
             ntrans = ntrans + 1
             if (ntrans > mtrans) then
                call error('gibbs2','maximum static transitions exceeded, increase mtrans',faterr)
             end if
             ptrans(ntrans) = pnew
             i1trans(ntrans) = iold
             i2trans(ntrans) = idmin
             !
             write (uout,'(2X,F12.4," --> ",F12.4,2X,A10)')&
                pold, pnew, trim(adjustl(ph(iold)%name(1:leng(ph(iold)%name))))
             iold = idmin
             pold = pnew
          end if
       end if
    end do
    write (uout,'(2X,F12.4," --> ",F12.4,2X,A10)')&
       pold, plist(nps), trim(adjustl(ph(iold)%name(1:leng(ph(iold)%name))))
    write (uout,*)

  end subroutine static_transp

  subroutine dyn_transp()
    use gnuplot_templates

    real*8 :: gmin
    real*8 :: ptrans(mtrans)
    integer :: i1trans(mtrans), i2trans(mtrans)
    integer :: blknum(mtrans)
    integer :: ntrans
    logical :: trchange, ntranszero
    integer :: i, j, k, idmin, iold, jj
    integer :: iblk
    integer :: lu, lu2, count, n
    character*5 :: starts
    character*3 :: ends
    real*8 :: pold, pnew
    integer :: ndo, imask(nph)

    call mask_trans(ndo,imask)
    if (ndo < 2) return

    write (uout,'("  Writing file : ",A/)') trim(fileroot) // "_ptrans.gnu"
    lu = opengnu(trim(fileroot)//"_ptrans")
    write(lu,'("set ylabel ""p (GPa)""")')
    write(lu,'("set xlabel ""T (K)""")')
    write(lu,'("set yrange [0:",F20.4,"]")') plist(nps)
    write(lu,'("set xrange [0:",F20.4,"]")') tlist(nts)
    write(lu,'("unset key")')
    
    ! transition pressures (T)
    write (uout,'("* Dynamic transition pressures (linear interpolation)")')
    lu2 = fopen(lu2,trim(fileroot)//".ptrans"//null,iowrite)
    write (uout,'("  Writing file : ",A/)') trim(fileroot)//".ptrans"

    ! header
    iblk = 0
    ntrans = 0
    ntranszero = .true.
    do i = 1, nts
       ntrans = 0
       trchange = .false.
       iold = -1
       do j = 1, nps
          idmin = -1
          gmin = 1d30
          do k = 1, ndo
             if (ph(imask(k))%didtp(i,j) .and. ph(imask(k))%gtp(i,j) < gmin) then
                gmin = ph(imask(k))%gtp(i,j)
                idmin = imask(k)
             end if
          end do
          if (idmin == -1) cycle

          if (iold < 0) then
             iold = idmin
             pold = 0d0
          else if (idmin /= iold) then
             pnew = plist(j-1) - (ph(idmin)%gtp(i,j-1)-ph(iold)%gtp(i,j-1)) * (plist(j) - plist(j-1)) /&
                (ph(idmin)%gtp(i,j)-ph(iold)%gtp(i,j) - (ph(idmin)%gtp(i,j-1)-ph(iold)%gtp(i,j-1)))
             !
             ntrans = ntrans + 1
             if (ntrans > mtrans) then
                call error('gibbs2','maximum static transitions exceeded, increase mtrans',faterr)
             end if
             trchange = .not.(i1trans(ntrans) == iold .and. i2trans(ntrans) == idmin)
             ptrans(ntrans) = pnew
             i1trans(ntrans) = iold
             i2trans(ntrans) = idmin

             iold = idmin
             pold = pnew
          end if
       end do
       if (trchange .or. (ntrans == 0 .and. .not.ntranszero)) then
          ntranszero = (ntrans == 0)
          iblk = iblk + 1
          blknum(iblk) = ntrans
          if (i /= 1) write (lu2,'(/)')
          write (lu2,'("#      T(K)     ",4X,999(A10,6X))') &
             (trim(adjustl(ph(i1trans(j))%name(1:leng(ph(i1trans(j))%name)))),j=1,ntrans), &
             trim(adjustl(ph(i2trans(ntrans))%name(1:leng(ph(i2trans(ntrans))%name))))
          write (lu2,'(1X,F12.4,2X,F12.4,999(A4,F12.4))') tlist(i), 0d0, &
             (" -> ",ptrans(j),j=1,ntrans), " -> ", plist(nps)
       else
          write (lu2,'(1X,F12.4,2X,F12.4,999(A4,F12.4))') tlist(i), 0d0, &
             (" -> ",ptrans(j),j=1,ntrans), " -> ", plist(nps)
       end if
    end do
    call fclose(lu2)

    ! complete gnu file
    starts = "plot "
    ends = ",\ "
    count = sum(blknum(1:iblk))
    n = 0
    do i = 1, iblk
       do j = 1, blknum(i)
          n = n + 1
          if (n == count) ends = "   "
          if (i /= 1 .or. j/= 1) starts = "     "
          
          jj = 2 + 2*j
          if (jj > 9) then
             write (lu,'(A5,"''",A,"'' u 1:",I2.2," index ",I2," w linespoints ls 1 ",A2)') &
                starts, trim(fileroot)//".ptrans", 2+2*j, i-1, ends
          else
             write (lu,'(A5,"''",A,"'' u 1:",I1," index ",I1," w linespoints ls 1 ",A2)') &
                starts, trim(fileroot)//".ptrans", 2+2*j, i-1, ends
          end if
       end do
    end do

    call closegnu(trim(fileroot)//"_ptrans",lu)

  end subroutine dyn_transp

  subroutine dyneos()
    use gnuplot_templates
    use debye
    use eos
    use fit
    use evfunc

    integer :: lu
    integer :: i, j, k, ierr, id
    real*8 :: v, b, e, g
    character*(mline) :: msg, msg2
    character*(mline_fmt) :: fm, fme
    type(fitinfo) :: pfit
    real*8 :: proplist(mpropout), errlist(mpropout)
    type(fitpack) :: ft

    if (writelevel < 1) return

    write (uout,'("* Calculated temperature effects ")')
    write (uout,'("  Writing file : ",A/)') trim(fileroot)//".eos"
    lu = fopen(lu,trim(fileroot)//".eos"//null,iowrite)

    ! Loop over phases
    do i = 1, nph
       if (ph(i)%tmodel == tm_static) cycle

       ! header with property names
       write (uout,'("+ Phase ",A)') trim(adjustl(ph(i)%name(1:leng(ph(i)%name))))
       write (lu,'("# Phase ",A)') trim(adjustl(ph(i)%name(1:leng(ph(i)%name))))
       msg = "# "
       do j = 1, mpropout
          write (msg,'(A," ",I2.2,":",A)') trim(adjustl(msg)), j, trim(adjustl(propname(j)))
          if (mod(j,5) == 0) then
             write (lu,'(A)') trim(adjustl(msg))
             msg = "# "
          end if
       end do
       if (mod(mpropout,5) /= 0) then
          write (lu,'(A)') trim(adjustl(msg))
       end if

       ! prepare output formats
       fm = format_string_header(propfmt,(/2,(len(trim(adjustl(propname(j)))),j=1,mpropout) /) )
       write (lu,fm) "# ",(trim(adjustl(propname(j))),j=1,mpropout)
       fm = format_string(propfmt(1:mpropout),2)
       fme = format_string(propfmt(0:mpropout),0)

       ! allocate G(T,p), V(T,p) and B(T,p) arrays
       allocate(ph(i)%gtp(nts,nps), ph(i)%vtp(nts,nps), ph(i)%btp(nts,nps))
       allocate(ph(i)%didtp(nts,nps))

       ! fill in the static energy fit
       ft = ft_null
       ft%emode = ph(i)%fit_mode
       ft%nepol = ph(i)%npol
       ft%epol = ph(i)%cpol

       ! Loop over temperatures
       do j = 1, nts

          ! fit the total helmholtz free energy on the volume grid
          call fit_agrid_t(ph(i),tlist(j),ft%napol,ft%apol,ft%amode,ierr,.true.,.true.,.true.,pfit)
          if (ierr > 0) then
             write (uout,'(" T = ",F12.4," P = ",F12.4,"")') tlist(j), 0d0
             call error('dyneos','fit for A(x) not found',faterr)
          end if

          ! fit the vibrational entropy on the volume grid
          call fit_sgrid_t(ph(i),tlist(j),ft%nspol,ft%spol,ft%smode,ierr,.true.,.true.)
          if (ierr > 0) then
             write (uout,'(" T = ",F12.4," P = ",F12.4,"")') tlist(j), 0d0
             call error('dyneos','fit for -TS(x) not found',faterr)
          end if

          ! sum the pV term corresponding to each pressure and calculate
          ! the rest of properties, from min(G): etc(p,T).
          do k = 1, nps
             ! equilibrium volume (p,T)
             call fit_pshift(ft%amode,ph(i)%v,plist(k),ft%napol,ft%apol,v,b,e,g,ierr)
             if (ierr > 0 .or..not.is_dyn_v_in(ph(i),v)) then
                write (msg,'(" T = ",F10.2," P = ",F10.2,", G minimum not found, skipping.")') &
                   tlist(j), plist(k)
                call error('dyneos',msg,warning)
                ph(i)%didtp(j,k) = .false.
                ph(i)%gtp(j,k) = 0d0
                ph(i)%vtp(j,k) = 0d0
                ph(i)%btp(j,k) = 0d0
                cycle
             end if

             call dyneos_outprop(ph(i),v,tlist(j),plist(k),ft,pfit,proplist,errlist)
             write (lu,fm) proplist
             if (pfit%nfit > 0 .and. doerrorbar) then
                write (lu,fme) "e", errlist
             end if
             
             ! coordinated with preamble of varbas.f90
             ph(i)%didtp(j,k) = .true.
             ph(i)%gtp(j,k) = proplist(5) ! G(T,p)
             ph(i)%vtp(j,k) = proplist(3) ! V(T,p)
             ph(i)%btp(j,k) = proplist(9) ! BT(T,p)
          end do
          write (lu,'(/)')

          if (allocated(vlist)) then
             ! (V,T) properties
             do k = 1, nvs
                ! check volume is inside the known region
                call vbracket(ph(i),vlist(k),id,.true.)
                if (id == 0) then
                   write (msg,'(" T = ",F10.2," V = ",F12.4,", volume out of known region.")') &
                      tlist(j), vlist(k)
                   call error('dyneos',msg,warning)
                   cycle
                end if

                call dyneos_outprop(ph(i),vlist(k),tlist(j),undef,ft,pfit,proplist,errlist)
                write (lu,fm) proplist
                if (pfit%nfit > 0 .and. doerrorbar) then
                   write (lu,fme) "e", errlist
                end if
             end do
             write (lu,'(/)')
          end if

       end do
       write (lu,'(/)')

    end do
    call fclose(lu)

    write (uout,*)

    ! generate gnu files for all properties
    call gen_allgnu_t(trim(fileroot)//"_all_t")
    call gen_allgnu_p(trim(fileroot)//"_all_p")

  end subroutine dyneos

  ! Helper routine for dyneos. Calculates the output properties at a given
  ! volume/temperature/pressure
  subroutine dyneos_outprop(p,v,t,pres,ft,pfit,proplist,errlist)
    use fit

    type(phase), intent(in) :: p
    type(fitpack), intent(in) :: ft
    real*8, intent(in) :: v, t, pres
    type(fitinfo), intent(in) :: pfit
    real*8, dimension(mpropout), intent(out) :: proplist, errlist
    
    real*8, dimension(mpropout) :: prop, prop2, aux
    real*8 :: vol, b0, e, g
    integer :: i, id, ierr
    type(fitpack) :: ft_aux

    ! obtain properties
    call dyneos_calc(p,v,t,ft,proplist)

    errlist = 0d0
    ! error bars
    if (pfit%nfit > 0 .and. doerrorbar) then
       prop = 0d0
       prop2 = 0d0
       do i = 1, pfit%nfit
          ! is this coming from a fixed pressure or fixed volume calc.?
          if (pres /= undef) then
             ! obtain equilibrium volume for this fit
             call fit_pshift(pfit%mode(i),p%v,pres,pfit%npar(i),pfit%apar(:,i),vol,b0,e,g,ierr)
             if (ierr > 0) cycle

             ! check the volume is not out of the region active for T-calc.
             call vbracket(p,vol,id,.true.)
             if (id == 0) then
                write (uout,'(" Temperature = ",F12.4)') t
                write (uout,'(" Pressure = ",F12.4)') pres
                write (uout,'(" Volume = ",F17.7)') vol
                call error('dyneos_outprop','No error because of out-of-bounds v.',warning)
                return
             end if
          else
             vol = v
          end if

          ! insert the i-th polynomial into the auxiliary fitpack
          ft_aux = ft
          ft_aux%emode = p%pfit%mode(i)
          ft_aux%nepol = p%pfit%npar(i)
          ft_aux%epol  = p%pfit%apar(:,i)
          ft_aux%amode = pfit%mode(i)
          ft_aux%napol = pfit%npar(i)
          ft_aux%apol  = pfit%apar(:,i)

          ! calculate properties for this fit
          call dyneos_calc(p,vol,t,ft_aux,aux)
          
          ! add to averages
          prop = prop + aux * pfit%wei(i)
          prop2 = prop2 + aux * aux * pfit%wei(i)
       end do
       errlist = sqrt(max(prop2 - prop * prop,0d0))
    end if

  end subroutine dyneos_outprop

  ! Calculate peroperties at a given volume and temperature.
  subroutine dyneos_calc(p,v,t,ft,proplist)
    use debye
    use eos
    use evfunc

    type(phase), intent(in) :: p
    type(fitpack), intent(in) :: ft
    real*8, intent(in) :: v, t
    real*8, intent(out) :: proplist(mpropout)

    integer :: i, j
    integer :: ntpol, ierr
    real*8 :: f0, f1, f2, f3, f4, b0, b1, b2, tmp
    real*8 :: f0s, f1s, f2s, f3s, b0s, b1s, g
    real*8 :: nef, ef, theta
    real*8 :: tpol(0:mmpar)    
    real*8 :: pbeta, alpha, cp, bs, D, Derr
    real*8 :: cv_ac, cv_op, cv_lowt
    real*8 :: fvib, svib, uvib, cv_vib
    real*8 :: fel, sel, uel, cv_el, rdum
    real*8 :: fsum, ssum, usum, cv_sum
    real*8 :: gam_ac, gam_op, gamma
    real*8 :: pext, psta, pth, dg
    integer :: ini, end, id
    real*8 :: auxcpol(0:mmpar)
    logical :: gamma_from_s

    real*8, parameter :: feps = 1d-13
    
    ! calculate gamma from entropy fit?
    gamma_from_s = (p%tmodel == tm_qhafull) .or. (p%emodel /= em_no)

    ! static energy and helmholtz free energy volume derivatives
    if (ft%nepol == 0) call error('dyneos_calc','nepol = 0',faterr)
    if (ft%napol == 0) call error('dyneos_calc','napol = 0',faterr)
    f0s = fv0(ft%emode,v,ft%nepol,ft%epol)
    f1s = fv1(ft%emode,v,ft%nepol,ft%epol)
    f2s = fv2(ft%emode,v,ft%nepol,ft%epol)
    f3s = fv3(ft%emode,v,ft%nepol,ft%epol)
    f0  = fv0(ft%amode,v,ft%napol,ft%apol)
    f1  = fv1(ft%amode,v,ft%napol,ft%apol)
    f2  = fv2(ft%amode,v,ft%napol,ft%apol)
    f3  = fv3(ft%amode,v,ft%napol,ft%apol)
    f4  = fv4(ft%amode,v,ft%napol,ft%apol)

    ! pressure derivatives of the isothermal bulk modulus
    b0 = v * f2 
    b1 = -(1+v*f3/f2)
    b2 = ((f3+v*f4)/f2**2 - v*f3**2/f2**3) 

    ! thermal pressure 
    pext = -f1
    psta = -f1s
    pth = pext - psta

    ! vibrational contribution
    !  sets Uvib, Cv_vib, Svib, Fvib, theta and gamma
    call get_thetad(p,v,f2s,f3s,theta,gamma)
    select case(p%tmodel)
    case(tm_debye_input, tm_debye, tm_debyegrun, tm_debye_poisson_input)
       if (gamma_from_s .and. t < tlim_gamma) then
          call thermal (theta,t,D,Derr,uvib,cv_lowt,fvib,svib)
       end if
       call thermal (theta,t,D,Derr,uvib,cv_vib,fvib,svib)

    case(tm_debye_einstein) 
       ! Einstein or Debye-Einstein
       if (gamma_from_s .and. t < tlim_gamma) then
          call debeins (p,theta,tlim_gamma,v,D,Derr,uvib,cv_lowt,fvib,svib,cv_ac,cv_op)
       end if
       call debeins (p,theta,t,v,D,Derr,uvib,cv_vib,fvib,svib,cv_ac,cv_op)

       ! calculate gamma, handle low-temperature case
       if (cv_op > cvlim) then
          gam_ac = gamma
          b0s = v*f2s
          b1s = -(1+v*f3s/f2s)
          gam_op = (9d0*b0s*(b1s-1)+2*psta)/(6*(3*b0s-2*psta))
          gamma = (cv_ac * gam_ac + cv_op * gam_op) / cv_vib
       end if
    case(tm_qhafull) 
       ! use full phonon spectra at given V by linear interpolation
       ! then, gamma = V / Cv * (ds/dV)_T
       ! the region below tlim_gamma is avoided
       call thermalphon(p,t,v,uvib,cv_vib,fvib,svib,cv_lowt)

    case(tm_qhafull_espresso) 
       ! use full set of phonons at given V by linear interpolation
       call thermalomega(p,t,v,uvib,cv_vib,fvib,svib,cv_lowt,gamma)

    end select

    ! electronic contribution
    !  sets Uel, Cv_el, Sel, Fel
    if (t < tsmall_el) then
       uel = 0d0
       fel = 0d0
       sel = 0d0
       cv_el = 0d0
    else
       select case(p%emodel)
       case(em_sommerfeld)
          if (p%efree) then
             ef = (3d0 * pisquare * p%nelec / v)**twothird / 2d0
          else
             call vbracket(p,v,id,.false.)
             if (id < 0) then
                nef = p%nefermi(-id)
             else
                nef = p%nefermi(id) + (v-p%v(id)) / (p%v(id+1)-p%v(id)) * (p%nefermi(id+1)-p%nefermi(id))
             end if
             ef = 3d0 / 2d0 * p%nelec / nef
          end if
          sel = pi**2 / 2d0 * pckbau**2 * T * p%nelec / ef
          uel = T * sel / 2d0
          fel = -uel
          cv_el = sel
       case(em_pol4)
          call vbracket(p,v,id,.false.)
          auxcpol = 0d0
          if (id < 0) then
             auxcpol(1:4) = p%fel_cpol(1:4,-id) 
             fel = fv0(ftsel_fitmode,T,6,auxcpol)
             auxcpol(1:4) = p%tsel_cpol(1:4,-id) 
             sel = fv0(ftsel_fitmode,T,6,auxcpol)
             auxcpol(1:4) = p%fel_cpol(1:4,-id) - p%tsel_cpol(1:4,-id) 
             cv_el = fv1(ftsel_fitmode,T,6,auxcpol)
          else if (id > 0) then
             auxcpol(1:4) = p%fel_cpol(1:4,id) 
             fel  = fv0(ftsel_fitmode,T,6,auxcpol)
             auxcpol(1:4) = p%fel_cpol(1:4,id+1) 
             rdum = fv0(ftsel_fitmode,T,6,auxcpol)
             fel = fel + (v-p%v(id)) / (p%v(id+1)-p%v(id)) * (rdum-fel)
             !
             auxcpol(1:4) = p%tsel_cpol(1:4,id) 
             sel  = fv0(ftsel_fitmode,T,6,auxcpol)
             auxcpol(1:4) = p%tsel_cpol(1:4,id+1) 
             rdum = fv0(ftsel_fitmode,T,6,auxcpol)
             sel = sel + (v-p%v(id)) / (p%v(id+1)-p%v(id)) * (rdum-sel)
             !
             auxcpol(1:4) = p%fel_cpol(1:4,id) - p%tsel_cpol(1:4,id) 
             cv_el  = fv1(ftsel_fitmode,T,6,auxcpol)
             auxcpol(1:4) = p%fel_cpol(1:4,id+1) - p%tsel_cpol(1:4,id+1) 
             rdum = fv1(ftsel_fitmode,T,6,auxcpol)
             cv_el = cv_el + (v-p%v(id)) / (p%v(id+1)-p%v(id)) * (rdum-cv_el)
          else
             call error('dyneos_topcalc','volume out of bounds in em_pol4',faterr)
          end if
          uel = fel - sel
          sel = sel / (-T)
       case default
          fel = 0d0
          uel = 0d0
          sel = 0d0
          cv_el = 0d0
       end select
    end if

    ! sum up
    fsum = fvib + fel
    usum = uvib + uel
    ssum = svib + sel
    cv_sum = cv_vib + cv_el

    ! calculate gamma from entropy fit
    if (gamma_from_s) then
       if (ft%nspol == 0) call error('dyneos_calc','nspol = 0',faterr)
       if (t < tlim_gamma) then
          gamma = - v / cv_lowt * fv1(ft%smode,v,ft%nspol,ft%spol) / tlim_gamma
       else
          gamma = - v / cv_sum * fv1(ft%smode,v,ft%nspol,ft%spol) / t
       end if
    end if

    ! Rest of thermodynamic properties -- note g can be calculate in 2 ways:
    !  with Fsum from two fits or from quasiharmonic formula+.... dg is a measure 
    !  this inaccuracy.
    ! pbeta = (dp/dT)_V
    b0 = b0 * au2gpa
    b2 = b2 / au2gpa
    Pbeta = cv_sum * gamma / v * au2gpa
    alpha = Pbeta / b0
    tmp = 1d0 + gamma * alpha * t
    Cp = cv_sum * tmp 
    Bs = b0 * tmp
    g = f0 + pext * v
    dg = f0 - f0s - fsum

    ! conversion to output units
    g = g * hy2kjmol 
    dg = dg * hy2kjmol
    pext = pext * au2gpa
    pth = pth * au2gpa
    psta = psta * au2gpa
    cp = cp * hy2kjmol * 1000
    alpha = alpha * 1d5
    fvib = fvib * hy2kjmol
    uvib = uvib * hy2kjmol
    svib = svib * hy2kjmol * 1000
    cv_vib = cv_vib * hy2kjmol * 1000
    fel = fel * hy2kjmol
    uel = uel * hy2kjmol
    sel = sel * hy2kjmol * 1000
    cv_el = cv_el * hy2kjmol * 1000
    fsum = fsum * hy2kjmol
    usum = usum * hy2kjmol
    ssum = ssum * hy2kjmol * 1000
    cv_sum = cv_sum * hy2kjmol * 1000

    ! output properties list -> coordinated with preamble of varbas.f90 
    proplist( 1) = pext
    proplist( 2) = t
    proplist( 3) = v
    proplist( 4) = f0s
    proplist( 5) = g
    proplist( 6) = dg
    proplist( 7) = psta
    proplist( 8) = pth
    proplist( 9) = b0
    proplist(10) = usum
    proplist(11) = cv_sum
    proplist(12) = fsum
    proplist(13) = ssum
    proplist(14) = theta
    proplist(15) = gamma
    proplist(16) = alpha
    proplist(17) = pbeta
    proplist(18) = bs
    proplist(19) = cp
    proplist(20) = b1
    proplist(21) = b2
    proplist(22) = fvib
    proplist(23) = fel
    proplist(24) = uvib
    proplist(25) = uel
    proplist(26) = svib
    proplist(27) = sel
    proplist(28) = cv_vib
    proplist(29) = cv_el
    
  end subroutine dyneos_calc

  subroutine deltag()

    integer :: lu, ll
    integer :: i, j, k
    character*(1000) :: blin
    character*(mline) :: alin
    integer :: ndo, imask(nph)

    call mask_trans(ndo,imask)
    if (ndo < 2) return

    ! delta_G(T,p)
    write (uout,'("* DeltaG(T,p) = G_k(T,p) - G_1(T,p) data ")')
    write (uout,'("* Reference: ",I3)') imask(1)

    write (uout,'("  Writing file : ",A/)') trim(fileroot)//".dgtp"
    lu = fopen(lu,trim(fileroot)//".dgtp"//null,iowrite)

    write (lu,'("# DeltaG of phases 2...n referred to phase num.:",I3)') imask(1)
    write (lu,'("# DeltaG(T,p) = G_k(T,p) - G_1(T,p) (kJ/mol)")')
    write (lu,'("#    T(K)      p(GPa)",999(4X,A13,6X))') &
       (trim(adjustl(ph(imask(k))%name(1:leng(ph(imask(k))%name)))),k=2,ndo)
    do i = 1, nts
       do j = 1, nps
          write (blin,'(2(F10.4,1X))') tlist(i), plist(j)
          ll = 22
          do k = 2, ndo
             if (.not.(ph(imask(k))%didtp(i,j)).or..not.ph(imask(1))%didtp(i,j)) cycle
             if (ph(imask(k))%didtp(i,j)) then
                write(alin,'(1p,E17.9,2X)') ph(imask(k))%gtp(i,j)-ph(imask(1))%gtp(i,j)
             else
                write(alin,'(4X,A6,9X)') "n/a"
             end if
             blin = blin(1:ll) // "   " // alin(1:20)
             ll = ll + 23
          end do
          write (lu,'(A)') trim(blin)
       end do
       write (lu,*)
    end do
    call fclose(lu)

  end subroutine deltag

  subroutine stablevbg()

    real*8 :: gmin
    integer :: lu, idmin
    integer :: i, j, k
    integer :: ndo, imask(nph)

    call mask_trans(ndo,imask)
    if (ndo < 2) return

    write (uout,'("* G(T,p), V(T,p) and B(T,p) of the stable phase")')
    write (uout,'("  Writing file : ",A/)') trim(fileroot)//".tpstab"
    lu = fopen(lu,trim(fileroot)//".tpstab"//null,iowrite)

    do i = 1, nts
       do j = 1, nps
          gmin = 1d30
          idmin = -1
          do k = 1, ndo
             if (ph(imask(k))%didtp(i,j) .and. ph(imask(k))%gtp(i,j) < gmin) then
                gmin = ph(imask(k))%gtp(i,j)
                idmin = imask(k)
             end if
          end do
          if (idmin == -1) cycle
          write (lu,'(2(F10.4,1X),A10,3(F17.7,2X))') tlist(i), plist(j), &
             trim(adjustl(ph(idmin)%name(1:leng(ph(idmin)%name)))),&
             ph(idmin)%gtp(i,j), ph(idmin)%vtp(i,j), ph(idmin)%btp(i,j)
       end do
       write (lu,*)
    end do
    call fclose(lu)

  end subroutine stablevbg

  subroutine interpolate()
    use debye
    use evfunc

    integer :: i, j, k
    integer :: mint, n, lu
    real*8 :: v, p, t, e, x, b, g, fac, f1
    real*8, allocatable :: fi(:)
    integer, allocatable :: ifm(:), ilens(:)
    character*(mline_fmt) :: fm1, fms
    integer :: napol, ierr, imode
    real*8 :: apol(0:mmpar)

    if (writelevel < 2) return

    if (interp_input > 0) then
       if (.not.allocated(fint)) allocate(fint(mxint))
       if (.not.allocated(iint)) allocate(iint(mxint))
       do i = 1, nps
          if (interp_input > 1 .and. ph(i)%tmodel /= tm_static) then
             do j = 1, nts
                nxint = nxint + 1
                if (nxint > size(fint)) then
                   call realloc(fint,2*nxint)
                   call realloc(iint,2*nxint)
                end if
                fint(nxint) = plist(i)
                iint(nxint) = 3
                nxint = nxint + 1
                if (nxint > size(fint)) then
                   call realloc(fint,2*nxint)
                   call realloc(iint,2*nxint)
                end if
                fint(nxint) = tlist(j)
                iint(nxint) = 3
             end do
          else
             nxint = nxint + 1
             if (nxint > size(fint)) then
                call realloc(fint,2*nxint)
                call realloc(iint,2*nxint)
             end if
             fint(nxint) = plist(i)
             iint(nxint) = 2
          end if
       end do
    end if
    if (nxint <= 0) return

    mint = 0
    do i = 1, nph
       mint = max(mint,ph(i)%ninterp)
    end do
    if (mint == 0) return
    allocate(fi(mint),ifm(mint+4))

    write (uout,'("* Linear interpolation of phase satellite data")')
    lu = fopen(lu,trim(fileroot)//".interp"//null,iowrite)
    write (uout,'("  Writing file : ",A/)') trim(fileroot)//".interp"

    do i = 1, nph
       if (ph(i)%ninterp <= 0) cycle
       write (lu,'("# Phase ",A)') trim(adjustl(ph(i)%name(1:leng(ph(i)%name))))
       
       n = ph(i)%ninterp
       ifm(1:4) = (/1,ifmt_v,ifmt_p,ifmt_t/)
       ifm(5:n+4) = ifmt_interp
       fm1 = format_string_header(ifm(1:4),(/1,9,6,4/))
       write (lu,fm1) "#","V(bohr^3)","p(GPa)","T(K)"
       fm1 = format_string(ifm(2:n+4),1)
       ifm(4) = 10
       fms = format_string(ifm(2:n+4),1)

       j = 0
       do while (j < nxint)
          j = j + 1
          if (iint(j) == 1) then
             v = fint(j)
             t = 0d0
             f1 = fv1(ph(i)%fit_mode,v,ph(i)%npol,ph(i)%cpol)
             p = - f1 * au2gpa
          else if (iint(j) == 2) then
             p = fint(j)
             t = 0d0
             call fit_pshift(ph(i)%fit_mode,ph(i)%v,p,ph(i)%npol,ph(i)%cpol,v,b,e,g,ierr)
             if (ierr > 0) then
                write (uout,'(" Pressure = ",F12.4)') p
                call error('dyneos','minimum of static H not found',warning)
                cycle
             end if
          else if (ph(i)%tmodel /= tm_static) then
             p = fint(j)
             t = fint(j+1)
             j = j + 1

             call fit_agrid_t(ph(i),t,napol,apol,imode,ierr,.true.,.true.,.true.)
             call fit_pshift(imode,ph(i)%v,p,napol,apol,v,b,e,g,ierr)
             if (ierr > 0) then
                write (uout,'(" Temperature = ",F12.4)') t
                write (uout,'(" Pressure = ",F12.4)') p
                call error('dyneos','minimum of G(x) not found',warning)
                cycle
             end if
          else
             cycle
          end if
          
          if (v <= ph(i)%v(1) .or. v >= ph(i)%v(ph(i)%nv)) then
             write (uout,'("Volume = ",F12.4)') v
             call error('interpolate','Volume out of V-range for this phase',warning)
             cycle
          else
             do k = 1, ph(i)%nv-1
                if (ph(i)%v(k) <= v .and. v <= ph(i)%v(k+1)) exit
             end do
             fac = (v-ph(i)%v(k)) / (ph(i)%v(k+1)-ph(i)%v(k))
             fi(1:n) = (1d0-fac) * ph(i)%interp(k,1:n) + fac * ph(i)%interp(k+1,1:n)
          end if
          
          if (iint(j) == 1 .or. iint(j) == 2) then
             write (lu,fms) v, p, "  static  ", fi(1:n)
          else
             write (lu,fm1) v, p, t, fi(1:n)
          end if

       end do
       
       write (lu,'(/)')
    end do

    deallocate(fi,ifm)
    call fclose(lu)

  end subroutine interpolate

  subroutine eshift_vexp()
    use debye
    use evfunc
    
    real*8, parameter :: facprec = 1d-10

    integer :: i, j, ierr
    integer :: napol
    real*8 :: apol(0:mmpar), psum, f1
    character*(mline_fmt) :: fm
    character*(mline) :: msg
    
    integer :: imode, niter
    real*8 :: vexpt, bexpt, fac
    real*8 :: v0, b0, e0, g0, e0new, g0new, v0t, b0t, e0t, g0t
    real*8 :: vexp, bexp
    real*8 :: q1, q2
    real*8 :: fa, fb, qfa, qfb, qfx
    real*8 :: psta_vexpt, pth_vexpt, bsta_vexpt
    real*8 :: bt_vexpt, bpobj
    real*8 :: vold

    real*8, parameter :: vcycle = 1d-7
    integer, parameter :: miter = 20

    ! scaling of E(V) to experimental data
    do i = 1, nph
       if (ph(i)%scaltype == scal_noscal) cycle

       ! fill static properties to uncorrected static energy
       call fit_pshift(ph(i)%fit_mode,ph(i)%v,0d0,ph(i)%npol,ph(i)%cpol,v0,b0,e0,g0,ierr)
       ph(i)%veq_static = v0
       ph(i)%eeq_static = g0 ! because the fit_pshift had p = 0.
       ph(i)%beq_static = b0
       b0 = b0 / au2gpa

       ! header
       write (uout,'("* Scaling of phase ",I2," (",A,")")') i, &
          trim(adjustl(ph(i)%name(1:leng(ph(i)%name))))
       write (uout,'("  Experimental T (K): ",F10.2)') ph(i)%eec_t
       write (uout,'("  Experimental p (GPa): ",F10.2)') ph(i)%eec_p
       write (uout,'("  Vexp at T and p (bohr^3): ",F12.3)') ph(i)%vscal
       select case(ph(i)%scaltype)
       case(scal_pshift)
          write (uout,'("+ Type: pshift")')
       case(scal_apbaf)
          write (uout,'("+ Type: apbaf")')
       case(scal_bpscal)
          write (uout,'("+ Type: bpscal")')
          write (uout,'("  Bexp at T and p (bohr^3): ",F12.3)') ph(i)%bscal
       end select

       niter = 0
1      continue
       niter = niter + 1
       if (niter > miter) call error('eshift_vexp','too many iterations in EEC',faterr)
       vold = v0
       write (uout,'("  - Iteration ",I2)') niter
       write (uout,'("      V0 (static,bohr^3) = ",F12.3)') vold
       write (uout,'("      B0 (static,bohr^3) = ",F14.2)') b0 * au2gpa

       ! fill the debye temperature, for debye models
       call fill_thetad(ph(i),.false.)

       ! find v0, b0 and g0 at temperature t0
       call fit_agrid_t(ph(i),ph(i)%eec_t,napol,apol,imode,ierr,.true.,.true.,.true.)
       call fit_pshift(imode,ph(i)%v,ph(i)%eec_p,napol,apol,v0t,b0t,e0t,g0t,ierr)
       if (ierr > 0) then
          write (msg,'("Phase ",A,": Can not find minimum of F(V,T). Using NOSCAL.")')&
             trim(adjustl(ph(i)%name(1:leng(ph(i)%name))))
          call error('eshift_vexp',msg,warning)
          cycle
       end if

       ! experimental input
       vexpt = ph(i)%vscal
       bexpt = ph(i)%bscal / au2gpa

       if (ph(i)%scaltype == scal_use) then
          ph(i)%scaltype = ph(ph(i)%iscal)%scaltype
          ph(i)%scale_a1 = ph(ph(i)%iscal)%scale_a1
          ph(i)%scale_a2 = ph(ph(i)%iscal)%scale_a2
       else
          if (vexpt > ph(i)%v(ph(i)%nv)) call error('eshift_vexp','Experimental volume outside of input volume grid.',faterr)
       end if

       ! Calculate scaling parameters if not known
       if (ph(i)%iscal == 0) then
          select case(ph(i)%scaltype)
          case(scal_pshift)
             ! shift the energy, e' = e + p*v
             psta_vexpt = -fv1(ph(i)%fit_mode,vexpt,ph(i)%npol,ph(i)%cpol)
             psum = -fv1(imode,vexpt,napol,apol) - ph(i)%eec_p / au2gpa
             pth_vexpt = psum - psta_vexpt
             ph(i)%scale_a1 = psum 
          case(scal_apbaf)
             ! shift the energy, e' = e + p*v
             psum = -fv1(imode,vexpt,napol,apol) - ph(i)%eec_p / au2gpa
             ph(i)%scale_a1 = psum 
             ph(i)%scale_a2 = vexpt
          case(scal_bpscal)
             ! propeties at vexpt
             psta_vexpt = -fv1(ph(i)%fit_mode,vexpt,ph(i)%npol,ph(i)%cpol)
             pth_vexpt = -fv1(imode,vexpt,napol,apol) - psta_vexpt - ph(i)%eec_p / au2gpa
             bsta_vexpt = vexpt * fv2(ph(i)%fit_mode,vexpt,ph(i)%npol,ph(i)%cpol)
             bt_vexpt = vexpt * fv2(imode,vexpt,napol,apol) 
             bpobj = -(bexpt - bt_vexpt + bsta_vexpt) / pth_vexpt

             ! bracket
             if (bpobj < 0d0) then
                fa = 1d0 + 1d-5
                fb = min(1.3d0,ph(i)%v(ph(i)%nv)/v0)
             else
                fa = 1d0 - 1d-5
                fb = max(0.7d0,ph(i)%v(1)/v0)
             end if
             qfa = qfac(fa)
             qfb = qfac(fb)
             do while(qfa*qfb > 0)
                if (bpobj < 0d0) then
                   fb = 0.99d0 * fb
                else
                   fb = 1.01d0 * fb
                end if
                qfb = qfac(fb)
             end do

             ! bisection
             do while(abs(fa-fb) > facprec)
                fac = 0.5d0 * (fa+fb)
                qfx = qfac(fac)
                if (qfx * qfa < 0) then
                   fb = fac
                   qfb = qfx
                else
                   fa = fac
                   qfa = qfx
                end if
             end do
             fac = 0.5d0 * (fa+fb)
             qfx = qfac(fac)

             ! final static vexp and bexp scaling parameters
             vexp = vexpt / fac
             bexp = b0 * pth_vexpt / fv1(ph(i)%fit_mode,fac * v0,ph(i)%npol,ph(i)%cpol) 
             ph(i)%scale_a1 = vexp * bexp / v0 / b0
             ph(i)%scale_a2 = v0 / vexp
             
          end select
       end if

       ! do the scaling using the calculated coefficients
       if (ph(i)%scaltype == scal_bpscal) then
          ph(i)%e = g0 + ph(i)%scale_a1 * &
             (fv0(ph(i)%fit_mode,ph(i)%v*ph(i)%scale_a2,ph(i)%npol,ph(i)%cpol) - g0)
       else if (ph(i)%scaltype == scal_apbaf) then
          ph(i)%e = ph(i)%e - ph(i)%scale_a1 * ph(i)%scale_a2**2 / ph(i)%v
       else if (ph(i)%scaltype == scal_pshift) then
          ph(i)%e = ph(i)%e + ph(i)%scale_a1 * ph(i)%v
       end if

       ! recalculate the static energy fit and the energy minimum
       ! static fit
       call fit_ev(ph(i)%fit_mode, ph(i)%reg_mode, ph(i)%v, ph(i)%e, ph(i)%npol,&
          ph(i)%cpol, ierr, .false., ph(i)%nfix, ph(i)%idfix, ph(i)%obelix)
       if (ierr > 0) call error('eshift_vexp','Minimum E(x) not found',faterr)
       call phase_checkfiterr(ph(i),.false.)
       call fit_pshift(ph(i)%fit_mode,ph(i)%v,0d0,ph(i)%npol,ph(i)%cpol,v0,b0,e0new,g0new,ierr)
       ph(i)%e = ph(i)%e + g0 - g0new
       write (uout,'("  Energy step: ",1p,E20.12)') g0new-g0

       ! recalculate again to set the energy zero and write the ph(i) fields
       call fit_ev(ph(i)%fit_mode, ph(i)%reg_mode, ph(i)%v, ph(i)%e, ph(i)%npol,&
          ph(i)%cpol, ierr, .false., ph(i)%nfix, ph(i)%idfix, ph(i)%obelix)
       if (ierr > 0) then
          call error('eshift_vexp','Minimum E(x) not found',faterr)
       end if
       call phase_checkfiterr(ph(i),.false.)
       call fit_pshift(ph(i)%fit_mode,ph(i)%v,0d0,ph(i)%npol,ph(i)%cpol,v0,b0,e0new,g0new,ierr)
       ph(i)%veq_static = v0
       ph(i)%eeq_static = g0new
       ph(i)%beq_static = b0
       write (uout,'("  Energy step after E0 correction: ",1p,E20.12)') g0new-g0

       ! pmin and pmax
       ph(i)%pmax = -fv1(ph(i)%fit_mode,ph(i)%v(2),ph(i)%npol,ph(i)%cpol) * au2gpa
       ph(i)%pmin = -fv1(ph(i)%fit_mode,ph(i)%v(ph(i)%nv),ph(i)%npol,ph(i)%cpol) * au2gpa

       if (abs(vold-v0)/(vold+v0+1d-12) > vcycle .and. (ph(i)%tmodel == tm_debye .or.&
          ph(i)%tmodel == tm_debye_einstein .or. ph(i)%tmodel == tm_debyegrun)) goto 1

       ! ! output new energy and static p
       ! fm = format_string_header((/1,ifmt_v,ifmt_e,ifmt_e,ifmt_p/),(/1,9,9,9,10/))
       ! write (uout,fm) "#", "V(bohr^3)", "E_inp(Ha)", "E_new(Ha)","p_sta(GPa)"
       ! fm = format_string((/ifmt_v,ifmt_e,ifmt_e,ifmt_p/),1)
       ! do j = 1, ph(i)%nv
       !    f1 = fv1(ph(i)%fit_mode, ph(i)%v(j), ph(i)%npol, ph(i)%cpol)
       !    write (uout,fm) ph(i)%v(j), ph(i)%aux(j), ph(i)%e(j),-f1*au2gpa
       ! end do
       write (uout,*)
    end do

  contains
    function qfac(ff)

      real*8, intent(in) :: ff
      real*8 :: qfac
      
      real*8 :: facv0, psta_fv0, bsta_fv0

      facv0 = ff * v0
      psta_fv0 = -fv1(ph(i)%fit_mode,facv0,ph(i)%npol,ph(i)%cpol) 
      bsta_fv0 = facv0 * fv2(ph(i)%fit_mode,facv0,ph(i)%npol,ph(i)%cpol) 
      
      qfac = bsta_fv0 / psta_fv0 - bpobj

    end function qfac
  end subroutine eshift_vexp

  !< Returns fitting parameters of F(V,T).
  subroutine fit_agrid_t(p,T,napol,apol,mode,ierr,dostatic,dovib,doel,pfit)
    use debye

    type(phase), intent(in) :: p 
    real*8, intent(in) :: T
    integer, intent(out) :: napol
    real*8, intent(out) :: apol(0:mmpar)
    integer, intent(out) :: mode
    integer, intent(out) :: ierr
    logical, intent(in) :: dostatic, dovib, doel
    type(fitinfo), intent(out), optional :: pfit

    integer :: i, j, rnv, ini, end, nend
    real*8 :: uvib, cv, cv2, ent, D, Derr, cv_ac, cv_op, uel, sel
    real*8 :: realv(p%nv), reala(p%nv), aux(p%nv), dum, ef(p%nv)
    real*8 :: auxcpol(0:mmpar)

    real*8, parameter :: feps = 1d-13

    rnv = count(p%dyn_active)

    aux = 0d0

    ! find Fvib(V,T) at the grid volumes and store in aux
    if (dovib) then
       do i = 1, p%nv
          if (.not.p%dyn_active(i)) cycle
          if (p%tmodel == tm_qhafull) then
             call thermalphon(p,T,p%v(i),uvib,cv,aux(i),dum,cv2)
          else if (p%tmodel == tm_debye_einstein) then
             call debeins (p,p%td(i),T,p%v(i),D,Derr,uvib,cv,aux(i),dum,cv_ac,cv_op)
          else if (p%tmodel == tm_qhafull_espresso) then
             call thermalomega(p,T,p%v(i),uvib,cv,aux(i),dum,cv2)
          else 
             call thermal (p%td(i),T,D,Derr,uvib,cv,aux(i),dum)
          end if
       end do
    end if

    ! build Helmholtz free energy
    if (dostatic) then
       aux = p%e + aux 
    end if
    
    if (doel .and. T > tsmall_el) then
       if (p%emodel == em_sommerfeld) then
          if (p%efree) then
             ef = (3d0 * pisquare * p%nelec / p%v)**twothird / 2d0
          else
             ef = 3d0 / 2d0 * p%nelec / p%nefermi
          end if
          aux = aux - pi**2 / 2d0 * pckbau**2 * T * p%nelec / ef
       else if (p%emodel == em_pol4) then
          do i = 1, p%nv
             auxcpol = 0d0
             auxcpol(1:4) = p%fel_cpol(1:4,i) 
             aux(i) = aux(i) + fv0(ftsel_fitmode,T,6,auxcpol)
          end do
       end if
    end if

    ! apply mask to remove points with negative frequencies
    realv(1:rnv) = pack(p%v,p%dyn_active)
    reala(1:rnv) = pack(aux,p%dyn_active)

    ! Numerical fit -> napol and apol
    call fit_ev(p%fit_mode, p%reg_mode, realv(1:rnv), reala(1:rnv), napol, apol,&
       ierr, .false., pfit=pfit)
    mode = p%fit_mode

  end subroutine fit_agrid_t

  !< Returns fitting parameters of F(V,T).
  subroutine fit_sgrid_t(p,T,nspol,spol,mode,ierrs,dovib,doel)
    use debye

    type(phase), intent(in) :: p
    real*8, intent(in) :: T
    integer, intent(out) :: nspol
    real*8, intent(out) :: spol(0:mmpar)
    integer, intent(out) :: mode
    integer, intent(out) :: ierrs
    logical, intent(in) :: dovib, doel

    integer :: i, rnv
    real*8 :: uvib, cv, cv2, ent, D, Derr, cv_ac, cv_op
    real*8 :: realv(p%nv), reals(p%nv), aux(p%nv), dum
    real*8 :: t0, ef(p%nv)
    real*8 :: auxcpol(0:mmpar)

    ! only active
    rnv = count(p%dyn_active)

    ! 0 K calculations of gamma are tricky -> use ~30 K instead
    t0 = max(T,tlim_gamma)
    
    ! find Svib(V,T) at the grid volumes and store in aux
    aux = 0d0

    if (dovib) then
       ! vibrational contribution
       do i = 1, p%nv
          if (.not.p%dyn_active(i)) cycle
          if (p%tmodel == tm_qhafull) then
             call thermalphon(p,t0,p%v(i),uvib,cv,dum,aux(i),cv2)
          else if (p%tmodel == tm_debye_einstein) then
             call debeins (p,p%td(i),t0,p%v(i),D,Derr,uvib,cv,dum,aux(i),cv_ac,cv_op)
          else if (p%tmodel == tm_qhafull_espresso) then
             call thermalomega(p,t0,p%v(i),uvib,cv,dum,cv2,aux(i),cv2)
          else 
             call thermal (p%td(i),t0,D,Derr,uvib,cv,dum,aux(i))
          end if
       end do
       aux = aux * (-t0)
    end if

    if (doel) then
       ! electronic contribution
       if (p%emodel == em_sommerfeld) then
          if (p%efree) then
             ef = (3d0 * pisquare * p%nelec / p%v)**twothird / 2d0
          else
             ef = 3d0 / 2d0 * p%nelec / p%nefermi
          end if
          aux = aux + pi**2 / 2d0 * pckbau**2 * T * p%nelec / ef
       else if (p%emodel == em_pol4) then
          do i = 1, p%nv
             auxcpol = 0d0
             auxcpol(1:4) = p%tsel_cpol(1:4,i) 
             aux(i) = aux(i) + fv0(ftsel_fitmode,t0,6,auxcpol)
          end do
       end if
    end if

    ! apply mask to remove points with negative frequencies
    realv(1:rnv) = pack(p%v,p%dyn_active)
    reals(1:rnv) = pack(aux,p%dyn_active)

    ! -T*S fit for the calculation 
    call fitt_polygibbs(p%sfit_mode,realv(1:rnv),reals(1:rnv),nspol,spol,ierrs,.false.)
    mode = p%sfit_mode

    ! no centering slope in entropy fits
    nspol = nspol + 1
    spol(nspol) = 0d0

  end subroutine fit_sgrid_t

  !< Returns fitting parameters of Cv_el(V,T).
  subroutine fit_cvelgrid_t(p,T,ncvpol,cvpol,mode,ierr)
    use debye

    type(phase), intent(in) :: p 
    real*8, intent(in) :: T
    integer, intent(out) :: ncvpol
    real*8, intent(out) :: cvpol(0:mmpar)
    integer, intent(out) :: mode
    integer, intent(out) :: ierr

    integer :: i
    integer :: rnv
    real*8 :: aux(p%nv), auxcpol(0:mmpar)
    real*8 :: realv(p%nv), realcv(p%nv)

    rnv = count(p%dyn_active)

    aux = 0d0

    if (p%emodel == em_pol4) then
       do i = 1, p%nv
          auxcpol = 0d0
          auxcpol(1:4) = p%fel_cpol(1:4,i) - p%tsel_cpol(1:4,i)
          aux(i) = fv1(ftsel_fitmode,T,6,auxcpol)
       end do
    end if

    ! apply mask to remove points with negative frequencies
    realv(1:rnv) = pack(p%v,p%dyn_active)
    realcv(1:rnv) = pack(aux,p%dyn_active)

    ! Numerical fit -> napol and apol
    call fit_ev(p%cvfit_mode, p%reg_mode, realv(1:rnv), realcv(1:rnv), ncvpol, cvpol,&
       ierr, .false.)
    mode = p%cvfit_mode

  end subroutine fit_cvelgrid_t

  function is_dyn_v_in(p,v)

    type(phase), intent(in) :: p
    real*8, intent(in) :: v
    logical :: is_dyn_v_in

    real*8, parameter :: veps2 = 1d-6

    integer :: i

    is_dyn_v_in = .false.
    if (v < p%v(1)-veps2 .or. v > p%v(p%nv)+veps2) return
    do i = 1, p%nv-1
       is_dyn_v_in = v >= p%v(i)-veps2    .and. v <=p%v(i+1)+veps2 &
                    .and. p%dyn_active(i) .and. p%dyn_active(i+1)
       if (is_dyn_v_in) return
    end do

  end function is_dyn_v_in

  subroutine drhouse(nsamples)

    integer, intent(in) :: nsamples

    type polyfit
       real*8 :: c(0:mmpar)
       real*8 :: rms2
       integer :: deg, npar
       integer :: ndat
       real*8 :: vmin, emin, b0, b1, b2
       real*8 :: w
    end type polyfit

    integer :: i, j, k, n, mv, minv
    real*8, allocatable :: v(:), x(:), y(:), w(:)
    logical, allocatable :: lx(:)
    integer :: nx(nsamples), nv
    integer :: nparmin, nparmax
    integer :: mode, deg, idx, idxx(1), ierr
    real*8 :: vscal, rms, rms2min, fdum1, fdum2, fdum3
    real*8 :: v0, f0, f1, f2, f3, f4
    character*(mline) :: fm, fme, msg
    real*8 :: prop(5), prop2(5), quo(5), amean, astd, yleft, yright
    type(polyfit), allocatable :: pol(:)
    real*8, allocatable :: emean(:), emean2(:), eeval(:,:)
    type(polyfit) :: pavg
    integer :: npol, nrej, nsrej, nsteps
    real*8, allocatable, dimension(:) :: xstep, ystep

    real*8 :: outlier_fac = 10d0
    real*8 :: step_fac = 100d0

    ! count max. number of volumes
    mv = 0
    minv = 999999999
    do i = 1, nph
       mv = max(mv,ph(i)%nv)
       minv = min(minv,ph(i)%nv)
    end do
    nparmin = mparmin
    nparmax = min(minv-5,mpar)

    ! prepare memory and weights
    allocate(x(mv),y(mv),v(mv),lx(mv),w(mv))
    allocate(pol(nsamples*(nparmax-nparmin+1)),emean(mv),emean2(mv))
    allocate(eeval(mv,nsamples*(nparmax-nparmin+1)))
    allocate(xstep(mv),ystep(mv))
    w = 1d0

    write (uout,'("* Diagnosis of problems in static data")')
    write (uout,'("  Average of polynomials, with degree range: ",I3,I3)') nparmin, nparmax
    write (uout,'("  Number of samples: ",I9)') nsamples
    write (uout,*)
    do i = 1, nph
       write (uout,'("+ Phase ",I3,": ",A)') i, ph(i)%name(1:leng(ph(i)%name))

       ! initialization
       nv = ph(i)%nv

       ! scaling volume
       idxx = minloc(ph(i)%e)
       idx = idxx(1)
       vscal = ph(i)%v(idx)

       ! statistics of polynomial fits
       lx = .false.
       x = 0d0
       rms2min = 1d30
       npol = 0
       nrej = 0
       nsrej = 0
       n = 0
       do while (n < nsamples)
          ! generate the dataset for this sample
          n = n + 1
          call random_number(x(1:nv))
          lx(1:nv) = (x(1:nv) > 0.5d0)
          if (count(lx(1:nv)) < nparmax+10) then
             nsrej = nsrej + 1
             n = n - 1
             if (nsrej > 2 * nsamples) then
                call error('drhouse','Too many datasets rejected.',faterr)
             end if
             cycle
          end if

          nx(n) = count(lx(1:nv))
          v(1:nx(n)) = pack(ph(i)%v,lx(1:nv))
          x(1:nx(n)) = v2str(fit_strain_bm,v(1:nx(n)),vscal)
          y(1:nx(n)) = pack(ph(i)%e,lx(1:nv))

          do deg = nparmin, nparmax
             ! fit the new polynomial
             npol = npol + 1
             pol(npol)%deg = deg
             pol(npol)%ndat = nx(n)

             mode = fit_strain * 10000 + fit_strain_bm * 100 + deg
             pol(npol)%c = 0d0
             call polfit(nx(n),1,nx(n),x(1:nx(n)),y(1:nx(n)),w,rms,deg,pol(npol)%c)
             pol(npol)%rms2 = rms*rms
             pol(npol)%npar = deg + 2
             pol(npol)%c(deg+1) = vscal
             pol(npol)%c(deg+2) = 0d0
             rms2min = min(rms2min,pol(npol)%rms2)

             call fit_pshift(mode,v(1:nx(n)),0d0,pol(npol)%npar,pol(npol)%c,&
                v0,fdum1,fdum2,fdum3,ierr)
             if (ierr > 0) then
                npol = npol - 1
                nrej = nrej + 1
                cycle
             end if

             pol(npol)%vmin = v0
             f0 = fv0(mode,v0,pol(npol)%npar,pol(npol)%c)
             f1 = fv1(mode,v0,pol(npol)%npar,pol(npol)%c)
             f2 = fv2(mode,v0,pol(npol)%npar,pol(npol)%c)
             f3 = fv3(mode,v0,pol(npol)%npar,pol(npol)%c)
             f4 = fv4(mode,v0,pol(npol)%npar,pol(npol)%c)
             pol(npol)%emin = f0
             pol(npol)%b0 = v0 * f2 * au2gpa
             pol(npol)%b1 = -(1+v0*f3/f2)
             pol(npol)%b2 = ((f3+v0*f4)/f2**2 - v0*f3**2/f2**3) / au2gpa
          end do
       end do 

       ! weights
       pol(1:npol)%w = exp(-pol(1:npol)%rms2 / rms2min * &
          real(pol(1:npol)%deg,8) / real(pol(1:npol)%ndat,8))
       pol(1:npol)%w = pol(1:npol)%w / sum(pol(1:npol)%w)

       ! average polynomial
       pavg%deg = nparmax
       pavg%npar = pavg%deg + 2
       pavg%c = 0d0
       prop = 0d0
       prop2 = 0d0
       eeval = 0d0
       do j = 1, npol
          pavg%c(0:pol(j)%deg) = pavg%c(0:pol(j)%deg) + pol(j)%w * pol(j)%c(0:pol(j)%deg)
          prop = prop + pol(j)%w * (/pol(j)%vmin, pol(j)%emin, pol(j)%b0, pol(j)%b1, pol(j)%b2/)
          prop2 = prop2 + pol(j)%w * (/pol(j)%vmin, pol(j)%emin, pol(j)%b0, pol(j)%b1, pol(j)%b2/)**2
          eeval(1:nv,j) = fv0(mode,ph(i)%v,pol(j)%npar,pol(j)%c)
       end do
       prop2 = sqrt(max(prop2 - prop*prop,0d0))
       pavg%c(nparmax+1) = vscal
       pavg%c(nparmax+2) = 0d0

       emean2 = 0d0
       do j = 1, npol
          emean = 0d0
          do k = 1, npol
             emean(1:nv) = emean(1:nv) + pol(k)%w * (eeval(1:nv,j)-eeval(1:nv,k))
          end do
          emean2(1:nv) = emean2(1:nv) + pol(j)%w * emean*emean
       end do
       emean2 = sqrt(max(emean2,0d0))
       
       pavg%ndat = ph(i)%nv
       pavg%w = 1d0
       call fit_pshift(mode,ph(i)%v,0d0,pavg%npar,pavg%c,v0,fdum1,fdum2,fdum3,ierr)
       if (ierr > 0) then
          call error('drhouse','No minimum found for average polynomial.',faterr)
       end if
       pavg%vmin = v0
       f0 = fv0(mode,v0,pavg%npar,pavg%c)
       f1 = fv1(mode,v0,pavg%npar,pavg%c)
       f2 = fv2(mode,v0,pavg%npar,pavg%c)
       f3 = fv3(mode,v0,pavg%npar,pavg%c)
       f4 = fv4(mode,v0,pavg%npar,pavg%c)
       pavg%emin = f0
       pavg%b0 = v0 * f2 * au2gpa
       pavg%b1 = -(1+v0*f3/f2)
       pavg%b2 = ((f3+v0*f4)/f2**2 - v0*f3**2/f2**3) / au2gpa
       pavg%rms2 = sum((ph(i)%e - fv0(mode,ph(i)%v,pavg%npar,pavg%c))**2)

       write (uout,'("  Bootstrap process.")') 
       write (uout,'("  Polynomials accepted: ",I9)') npol
       write (uout,'("  Polynomials rejected: ",I9)') nrej
       fm = format_string_header( &
          (/1,ifmt_order,ifmt_bp,ifmt_v,ifmt_e,ifmt_b,ifmt_bp,ifmt_bpp/),&
          (/1,1,6,9,5,6,2,10/))
       write (uout,fm) "#","n","weight","V(bohr^3)","E(Ha)","B(GPa)","Bp","Bpp(GPa-1)"
       fm = format_string((/ifmt_order,ifmt_bp,ifmt_v,ifmt_e,ifmt_b,ifmt_bp,ifmt_bpp/),1)
       do j = 1, npol
          write (uout,fm) pol(j)%deg, pol(j)%w, pol(j)%vmin, &
             pol(j)%emin, pol(j)%b0, pol(j)%b1, pol(j)%b2
       end do
       
       quo = abs(prop2 / prop)
       fm = format_string((/15,ifmt_v,ifmt_e,ifmt_b,ifmt_bp,ifmt_bpp/),1)
       write (uout,fm) "-average pol.--", pavg%vmin, pavg%emin, pavg%b0, pavg%b1, pavg%b2
       write (uout,fm) "--dir. average-", prop
       write (uout,fm) "---std. dev.---", prop2
       write (uout,fm) "--|std./avg.|--", quo
       write (uout,*)
       write (uout,'("* OK (|STD./AVG.| < 1)?    ",5(L8,X))') (quo<=1d0)
       write (uout,'("* OK (|STD./AVG.| < 0.01)? ",5(L8,X))') (quo<=0.01d0)

       if (any(quo > 1d0)) call error('drhouse','For one of {V,E,B,Bp,Bpp}, stdev/mean > 1',warning)
       if (any(quo > 0.01d0)) call error('drhouse','For one of {V,E,B,Bp,Bpp}, stdev/mean > 0.01',warning)
       write (uout,*)
    end do

    deallocate(x,y,v,lx,w,pol,emean,emean2,eeval,xstep,ystep)

  end subroutine drhouse

  subroutine mask_trans(n,imask)
    
    integer, intent(out) :: n
    integer, intent(out) :: imask(nph)

    integer :: i

    n = 0
    imask = 0
    if (writelevel < 2 .or..not.dotrans) return

    do i = 1, nph
       if (ph(i)%tmodel == tm_static) cycle
       if (.not.allocated(ph(i)%didtp)) cycle
       n = n + 1
       imask(n) = i
    end do

  end subroutine mask_trans

  subroutine printfreqs()

    integer :: i, j
    real*8 :: vol, pstat, bstat, veq, beq
    real*8 :: tmpvol, tmpb, tmppb, tmptot, f0
    real*8, allocatable :: freq(:)
    integer :: mfreq, n
    integer :: lu

    mfreq = 0
    do i = 1, nph
       if (ph(i)%tmodel == tm_debye_einstein) then
          mfreq = max(mfreq,ph(i)%nfreq)
       end if
    end do
    if (mfreq == 0) return
    allocate(freq(mfreq))

    write (uout,'("* Gamma frequencies in the Debye-Einstein model/")')
    write (uout,'("  Writing file: ",A/)') trim(fileroot)//".gammafreq"

    lu = fopen(lu,trim(fileroot)//".gammafreq"//null,iowrite)
    do i = 1, nph
       if (ph(i)%tmodel == tm_debye_einstein) then
          write (lu,'("# Phase ",I2," (",A,")")') i, &
             trim(adjustl(ph(i)%name(1:leng(ph(i)%name))))
          write (lu,'("# V (bohr^3) freq1 freq2 ... freq(3n-3) (cm^-1)")')
          n = ph(i)%nfreq

          do j = 1, ph(i)%nv
             vol = ph(i)%v(j)
             pstat = -fv1(ph(i)%fit_mode, vol, ph(i)%npol, ph(i)%cpol) * au2gpa
             bstat = vol * fv2(ph(i)%fit_mode, vol, ph(i)%npol, ph(i)%cpol) * au2gpa
             veq = ph(i)%veq_static
             beq = ph(i)%beq_static

             tmpVol=(vol/veq)**(1d0/6d0)
             tmpB = (bstat/beq)**(half)
             tmpPB = (1d0 - twothird*pstat/bstat)**(half)
             tmpTot = tmpVol*tmpB*tmpPB

             freq(1:n) = ph(i)%freqg0 * tmpTot * hy2cm_1
             
             write (lu,'(F10.4,9999(F16.8))') ph(i)%v(j), freq(1:n) 
          end do
          write (lu,'(/)')
       end if
    end do

    call fclose(lu)
    deallocate(freq)

  end subroutine printfreqs

  subroutine write_energy(ini,step,end)

    real*8, intent(in) :: ini, step, end

    integer :: i, j, lu, nstep
    character*(mline) :: name
    real*8 :: vfac, efac, v, vini, vstep, vend, f0

    vini = ini
    vend = end
    vstep = step
    write (uout,'("* Output of static energy to external files"/)')
    if (step < 0d0 .and. ini < 0d0 .and. end < 0d0) then
       do i = 1, nph
          ! open external file
          write (name,'(A,"_",I2.2,".edat")') trim(fileroot), i
          write (uout,'("  Writing file : ",A/)') trim(name)
          lu = fopen(lu,trim(name)//null,iowrite)

          ! choose the approrpriate units, header
          write (lu,'("# Phase ",I3,": ",A)') i, ph(i)%name(1:leng(ph(i)%name))
          if (ph(i)%units_v == units_v_ang3) then
             vfac = bohr2angstrom**3
             write (lu,'("# volume ang^3")')
          else
             write (lu,'("# volume bohr^3")')
             vfac = 1d0
          end if
          if (ph(i)%units_e == units_e_ev) then
             write (lu,'("# energy ev")')
             efac = hy2ev
          else if (ph(i)%units_e == units_e_ry) then
             write (lu,'("# energy ry")')
             efac = 2d0
          else
             write (lu,'("# energy hartree")')
             efac = 1d0
          end if

          ! write to external file
          do j = 1, ph(i)%nv
             write (lu,'(1p,2(E22.12,2X))') ph(i)%v(j)*vfac, ph(i)%e(j)*efac
          end do
          call fclose(lu)
       end do
    else
       do i = 1, nph
          ! open external file
          write (name,'(A,"_",I2.2,".edat")') trim(fileroot), i
          write (uout,'("  Writing file : ",A/)') trim(name)
          lu = fopen(lu,trim(name)//null,iowrite)

          ! choose the approrpriate units, header
          write (lu,'("# Phase ",I3,": ",A)') i, ph(i)%name(1:leng(ph(i)%name))
          if (ph(i)%units_v == units_v_ang3) then
             vfac = bohr2angstrom**3
             write (lu,'("# volume ang^3")')
          else
             write (lu,'("# volume bohr^3")')
             vfac = 1d0
          end if
          vini = vini / vfac
          vstep = vstep / vfac
          vend = vend / vfac
          if (ph(i)%units_e == units_e_ev) then
             write (lu,'("# energy ev")')
             efac = hy2ev
          else if (ph(i)%units_e == units_e_ry) then
             write (lu,'("# energy ry")')
             efac = 2d0
          else
             write (lu,'("# energy hartree")')
             efac = 1d0
          end if

          ! write to external file
          nstep = floor((vend+1d-5-vini)/vstep)
          do j = 0, nstep
             v = vini + real(j,8) * vstep
             f0 = fv0(ph(i)%fit_mode, v, ph(i)%npol, ph(i)%cpol)
             write (lu,'(1p,2(E22.12,2X))') v*vfac, f0*efac
          end do
          call fclose(lu)
       end do
    end if

  end subroutine write_energy

end module topcalc
