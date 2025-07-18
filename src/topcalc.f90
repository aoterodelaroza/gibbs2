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

module topcalc
  use fit, only: mmpar
  implicit none

  public

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
     integer :: nspol = 0    ! -T*Svib
     integer :: ncvpol = 0   ! Cv (for externalfvib)
     ! fit modes
     integer :: emode = 0
     integer :: amode = 0
     integer :: smode = 0
     integer :: cvmode = 0
     ! coefficients
     real*8 :: epol(0:mmpar) = 0d0
     real*8 :: apol(0:mmpar) = 0d0
     real*8 :: spol(0:mmpar) = 0d0
     real*8 :: cvpol(0:mmpar) = 0d0
  end type fitpack
  type(fitpack), save :: ft_null

contains

  ! Initialize variables for this module.
  subroutine topcalc_init()

    nxint = 0
    interp_input = 0

  end subroutine topcalc_init

  ! Write basic information about the run to the standard output.
  subroutine popinput(fileout)
    use evfunc, only: pfit_gauss, pfit_slatec, pfit_mode, pweigh_mode, pweigh_slatec, &
       pweigh_gibbs1, pweigh_gibbs2
    use varbas, only: nph, ph, mm, vfree
    use fit, only: mpar
    use tools, only: leng
    use param, only: mline, uout, title, amu2au
    character*(mline), intent(in) :: fileout

    integer :: i

    write (uout,'("* Input ")')
    if (leng(title) > 0) &
       write (uout,'("  Title: ",A)') trim(adjustl(title(1:leng(title)-1)))
    write (uout,'("  Output file (lu=",I2,"): ",A)') uout, &
       fileout(1:leng(fileout))
    write (uout,'("  Units: output is in atomic units, except where noted.")')
    write (uout,'("  Number of atoms per primitive cell: ",I3)') vfree
    write (uout,'("  Molecular mass (amu): ",F17.8)') mm/amu2au
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

  ! Write information about the static fit to the standard output.
  subroutine popenergyfit()
    use evfunc, only: fv0, fv1
    use gnuplot_templates, only: opengnu, closegnu
    use varbas, only: nph, ph, writelevel, doefit
    use tools, only: leng, fopen, fclose
    use param, only: mline_fmt, uout, format_string, format_string_header, fileroot, iowrite,&
       null, au2gpa, ifmt_eprec, ifmt_v, ifmt_p
    integer :: lu
    integer :: i, j
    real*8 :: xdum, f0, f1, fac
    character*5 :: starts
    character*3 :: ends
    integer :: iact
    character*(mline_fmt) :: fm
    real*8 :: vmin, vmax, emin, emax

    integer, parameter :: igrid = 1001

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
    write (lu,'("set xrange ["F20.4":"F20.4"]")') vmin - 0.01d0 * (vmax-vmin), &
       vmax + 0.01d0 * (vmax-vmin)
    write (lu,'("set yrange ["F20.4":"F20.4"]")') emin - 0.01d0 * (emax-emin), &
       emax + 0.01d0 * (emax-emin)
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

  ! Make the delta-H vs. p plot using the information in the ga(:) array
  subroutine plotdh()
    use gnuplot_templates, only: opengnu, closegnu
    use varbas, only: nph, ph, nps, plist, writelevel, doplotdh, n_not_pv_min
    use tools, only: error, leng, fopen, fclose
    use param, only: uout, fileroot, iowrite, warning, null

    integer :: i, j, k
    character*5 :: starts
    character*3 :: ends
    integer :: lu

    if (.not.doplotdh .or. writelevel < 2) return
    if (n_not_pv_min() < 2) then
       if (nph >= 2) then
          call error('plotdh','Not enough phases with p = 0 minimum to plot enthalpies.',warning)
       end if
       return
    end if

    ! write the auxiliary file
    write (uout,'("* Plotting static DeltaH"/)')
    write (uout,'("  Writing file : ",A/)') trim(fileroot) // "_dH.aux"
    lu = fopen(lu,trim(fileroot) // "_dH.aux"//null,iowrite)
    write (lu,'("# Relative H(p) plot, auxiliary file. ")')
    write (lu,'("# Phase 1 taken as reference: ",A)') trim(adjustl(ph(1)%name(1:leng(ph(1)%name))))
    write (lu,'("# Contains delta-H(p) (in Hartree per formula unit) calculated on input pressures. ")')
    write (lu,'("# p(GPa)",999(2X,A13,4X))') &
       (trim(adjustl(ph(k)%name(1:leng(ph(k)%name)))),k=2,nph)
    do j = 1, nps
       if (ph(1)%static_v(j) < 0d0) cycle
       write (lu,'(F12.4,1X)',advance='no') plist(j)
       do k = 2, nph
          if (ph(k)%static_v(j) < 0d0) then
             write (lu,'(1p,A20,2X)',advance='no') "n/a"
          else
             write (lu,'(1p,E20.10,2X)',advance='no') ph(k)%static_g(j) - ph(1)%static_g(j)
          end if
       end do
       write (lu,*)
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
    do i = 1, nph
       if (i == 1) cycle
       if (i == nph) ends = "   "
       if (i /= 2) starts = "     "

       if (i < 9) then
          write (lu,'(A5,"''",A,"'' u 1:",I1," w points ls ",I2," title ''",A10,"'' ",A2)') &
             starts, trim(fileroot)//"_dH.aux", i, i, ph(i)%name(1:leng(ph(i)%name)), ends
       else
          write (lu,'(A5,"''",A,"'' u 1:",I2," w points ls ",I2," title ''",A10,"'' ",A2)') &
             starts, trim(fileroot)//"_dH.aux", i, i, ph(i)%name(1:leng(ph(i)%name)), ends
       end if
    end do
    call closegnu(trim(fileroot) // "_dH",lu)

  end subroutine plotdh

  ! Calculate the static transition pressures and write to output.
  subroutine static_transp()
    use varbas, only: nph, ph, nps, nph, plist, dotrans, n_not_pv_min
    use tools, only: error, leng, realloc, leng_null
    use param, only: uout, warning, mline

    integer :: i, j, nphase, idmin
    real*8, allocatable :: ptrans(:)
    integer, allocatable :: iphase(:), idtrans(:)
    character(mline) :: phname
    real*8 :: pold
    real*8 :: gajp, gaj1p, gajp1, gaj1p1

    integer, parameter :: idtrans_not_set = 0
    integer, parameter :: idtrans_invalid_to_valid = 1
    integer, parameter :: idtrans_valid_to_invalid = 2
    integer, parameter :: idtrans_phase_appeared = 3
    integer, parameter :: idtrans_phase_disappeared = 4
    integer, parameter :: idtrans_normal = 5
    integer, parameter :: idtrans_end_of_range = 6

    character*45, parameter :: reason(0:6) = (/&
       "End of pressure range                        ",&
       "Pressure Range of invalid phase finished     ",&
       "Minimum of G not found                       ",&
       "Pressure range of a more stable phase started",&
       "Pressure range of the stable phase ended     ",&
       "                                             ",&
       "End of pressure range                        "/)

    ! Do we need to run the calculation of the transition pressures?
    if (.not.dotrans) return
    if (n_not_pv_min() < 2) then
       if (nph >= 2) &
          call error('static_transp','Not enough phases with p = 0 minimum to find transition pressures.',warning)
       return
    end if

    ! write the header
    write (uout,'("* Static transition pressures (linear interpolation)")')
    write (uout,'("#        Pressure range (GPa)      Stable phase     Reason for the transition ")')

    ! allocate arrays to save the results
    allocate(ptrans(10),idtrans(10),iphase(10))
    ptrans = 0d0
    idtrans = idtrans_not_set
    iphase = 0
    nphase = 0

    ! do the first pressure by hand
    call find_min_phase(1,idmin)
    nphase = 1
    iphase(nphase) = idmin

    ! run the rest of the pressures
    do j = 2, nps
       ! calculate the new phase
       call find_min_phase(j,idmin)

       ! is this a new phase?
       if (idmin /= iphase(nphase)) then
          nphase = nphase + 1
          if (nphase > size(ptrans,1)) then
             call realloc(ptrans,2*nphase)
             call realloc(idtrans,2*nphase)
             call realloc(iphase,2*nphase)
          end if
          iphase(nphase) = idmin

          if (iphase(nphase-1) <= 0) then
             ! from unknown/invalid to valid
             ptrans(nphase-1) = plist(j)
             idtrans(nphase-1) = idtrans_invalid_to_valid
          elseif (iphase(nphase) <= 0) then
             ! from valid to unknown/invalid
             ptrans(nphase-1) = plist(j-1)
             idtrans(nphase-1) = idtrans_valid_to_invalid
          elseif (plist(j-1) < ph(iphase(nphase))%pmin) then
             ! a more stable phase just appeared
             ptrans(nphase-1) = plist(j)
             idtrans(nphase-1) = idtrans_phase_appeared
          elseif (plist(j) > ph(iphase(nphase-1))%pmax) then
             ! the stable phase just disappeared
             ptrans(nphase-1) = plist(j)
             idtrans(nphase-1) = idtrans_phase_disappeared
          else
             ! a usual phase to phase transition
             gajp   = ph(iphase(nphase))%static_g(j)
             gaj1p  = ph(iphase(nphase))%static_g(j-1)
             gajp1  = ph(iphase(nphase-1))%static_g(j)
             gaj1p1 = ph(iphase(nphase-1))%static_g(j-1)
             ptrans(nphase-1) = plist(j-1) - (gaj1p - gaj1p1) * (plist(j) - plist(j-1)) / (gajp - gajp1 - (gaj1p-gaj1p1))
             idtrans(nphase-1) = idtrans_normal
          end if
       end if
    end do
    ptrans(nphase) = plist(nps)
    idtrans(nphase) = idtrans_end_of_range

    ! write the results, first phase
    pold = plist(1)
    do i = 1, nphase
       if (iphase(i) <= 0) then
          phname = "invalid"
       else
          phname = trim(ph(iphase(i))%name)
       end if
       write (uout,'(2X,F12.4," --> ",F12.4,2X,A13,7X,A)') pold, ptrans(i), &
          trim(adjustl(phname(1:leng_null(phname)))),&
          trim(reason(idtrans(i)))
       pold = ptrans(i)
    end do
    write (uout,*)

  contains
    subroutine find_min_phase(jj,idmin_)
      integer, intent(in) :: jj
      integer, intent(out) :: idmin_

      real*8 :: gmin
      integer :: ii

      idmin_ = 0
      gmin = huge(1d0)
      do ii = 1, nph
         if (ph(ii)%static_g(jj) < gmin .and. ph(ii)%static_v(jj) > 0d0) then
            gmin = ph(ii)%static_g(jj)
            idmin_ = ii
         end if
      end do

    end subroutine find_min_phase

  end subroutine static_transp

  ! Calculate the dynamic transition pressures and write to the output.
  subroutine dyn_transp()
    use gnuplot_templates, only: opengnu, closegnu
    use varbas, only: nph, ph, nps, plist, nts, tlist
    use tools, only: error, leng, fopen, fclose
    use param, only: uout, fileroot, iowrite, faterr, null

    integer, parameter :: mtrans = 50

    real*8 :: gmin
    real*8 :: ptrans(mtrans)
    integer :: i1trans(0:mtrans), i2trans(0:mtrans)
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

  ! Calculate and write the EOS file at the requested temperature and
  ! pressure list.
  subroutine dyneos()
    use evfunc, only: fv0
    use debye, only: thermal_debye_extended, thermal_qha, thermal
    use gnuplot_templates, only: gen_allgnu_t, gen_allgnu_p
    use fit, only: fit_pshift, fitinfo, mmpar, fit_ev
    use varbas, only: nph, ph, mpropout, propfmt, tm_static, tm_externalfvib,&
       nps, plist, nts, tlist, tm_debye_extended, tm_qhafull, tm_debye,&
       writelevel, nvs, vlist, doerrorbar, propname, vbracket, vdefault,&
       tm_debye_input, tm_debyegrun, tm_debye_poisson_input
    use tools, only: error, leng, fopen, fclose
    use param, only: mline, mline_fmt, uout, format_string, format_string_header, fileroot,&
       iowrite, warning, faterr, null, undef
    integer :: lu
    integer :: i, j, k, ierr, id
    real*8 :: v, b, e, g, dum
    character*(mline) :: msg
    character*(mline_fmt) :: fm, fme
    type(fitinfo) :: pfit
    real*8 :: proplist(mpropout), errlist(mpropout), rms, may
    type(fitpack) :: ft
    logical :: isalloc
    real*8, allocatable :: xfit(:), yfit(:), ynew(:), f0(:)

    if (writelevel < 1) return

    write (uout,'("* Thermodynamic properties at the chosen (p,T) values ")')
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

       ! calculate fvib, S, CV on the volume and temperature grid
       if (allocated(ph(i)%dynamic_fvib)) deallocate(ph(i)%dynamic_fvib)
       if (allocated(ph(i)%dynamic_s)) deallocate(ph(i)%dynamic_s)
       if (allocated(ph(i)%dynamic_cv)) deallocate(ph(i)%dynamic_cv)
       allocate(ph(i)%dynamic_fvib(ph(i)%nv,nts))
       allocate(ph(i)%dynamic_s(ph(i)%nv,nts))
       allocate(ph(i)%dynamic_cv(ph(i)%nv,nts))

       if (allocated(ph(i)%ffit_npol)) deallocate(ph(i)%ffit_npol)
       if (allocated(ph(i)%ffit_apol)) deallocate(ph(i)%ffit_apol)
       if (allocated(ph(i)%tsfit_npol)) deallocate(ph(i)%tsfit_npol)
       if (allocated(ph(i)%tsfit_apol)) deallocate(ph(i)%tsfit_apol)
       if (allocated(ph(i)%cvfit_npol)) deallocate(ph(i)%cvfit_npol)
       if (allocated(ph(i)%cvfit_apol)) deallocate(ph(i)%cvfit_apol)
       allocate(ph(i)%ffit_npol(nts))
       allocate(ph(i)%ffit_apol(0:mmpar,nts))
       allocate(ph(i)%tsfit_npol(nts))
       allocate(ph(i)%tsfit_apol(0:mmpar,nts))
       allocate(ph(i)%cvfit_npol(nts))
       allocate(ph(i)%cvfit_apol(0:mmpar,nts))

       ! calculate zero-point energies
       allocate(f0(ph(i)%nv))
       do j = 1, ph(i)%nv
          if (ph(i)%tmodel == tm_debye_extended) then
             call thermal_debye_extended(ph(i),0d0,j,f0(j),dum,dum)
          elseif (ph(i)%tmodel == tm_qhafull) then
             call thermal_qha(ph(i),0d0,j,f0(j),dum,dum)
          elseif (ph(i)%tmodel == tm_debye.or.ph(i)%tmodel == tm_debye_input.or.&
             ph(i)%tmodel == tm_debyegrun.or.ph(i)%tmodel == tm_debye_poisson_input) then
             call thermal(ph(i)%td(j),0d0,dum,dum,dum,dum,f0(j),dum)
          else
             write (*,*) "fixme!"
             stop 1
          end if
       end do

       allocate(xfit(ph(i)%nv),yfit(ph(i)%nv),ynew(ph(i)%nv))
       do k = 1, nts
          ! calculate thermodynamic properties on the (V,T) grid
          do j = 1, ph(i)%nv
             if (ph(i)%tmodel == tm_debye_extended) then
                call thermal_debye_extended(ph(i),tlist(k),j,ph(i)%dynamic_fvib(j,k),&
                   ph(i)%dynamic_s(j,k),ph(i)%dynamic_cv(j,k))
             elseif (ph(i)%tmodel == tm_qhafull) then
                call thermal_qha(ph(i),tlist(k),j,ph(i)%dynamic_fvib(j,k),&
                   ph(i)%dynamic_s(j,k),ph(i)%dynamic_cv(j,k))
             elseif (ph(i)%tmodel == tm_debye.or.ph(i)%tmodel == tm_debye_input.or.&
                ph(i)%tmodel == tm_debyegrun.or.ph(i)%tmodel == tm_debye_poisson_input) then
                call thermal(ph(i)%td(j),tlist(k),dum,dum,dum,&
                   ph(i)%dynamic_cv(j,k),ph(i)%dynamic_fvib(j,k),ph(i)%dynamic_s(j,k))
             else
                write (*,*) "fixme!"
                stop 1
             end if
          end do

          ! print out the results
          write (uout,'("# Temperature = ",F12.4)') tlist(k)
          write (uout,'("# All properties per unit formula.")')
          write (uout,'("# V(bohr^3)      F(Ha)              Fvib(Ha)         Fvib-F0(Ha)        S(Ha/K)        CV(Ha/K)")')
          do j = 1, ph(i)%nv
             write (uout,'(F10.4,1X,1p,3(E18.10,1X),E14.6,1X,E14.6)') ph(i)%v(j), &
                ph(i)%e(j) + ph(i)%dynamic_fvib(j,k), &
                ph(i)%dynamic_fvib(j,k), ph(i)%dynamic_fvib(j,k)-f0(j),&
                ph(i)%dynamic_s(j,k), ph(i)%dynamic_cv(j,k)
          end do

          ! fit the F(V) as a function of temperature
          xfit = ph(i)%v(1:ph(i)%nv)
          yfit = ph(i)%e + ph(i)%dynamic_fvib(:,k)
          call fit_ev(ph(i)%fit_mode,ph(i)%reg_mode,xfit,yfit,ph(i)%ffit_npol(k),&
             ph(i)%ffit_apol(:,k),ierr,.false.)
          if (ierr > 0) then
             write (uout,'(" T = ",F12.4," V = ",F12.4,"")') tlist(k), ph(i)%v(j)
             call error('dyneos','fit for F(V;T) not found',faterr)
          end if
          ynew = fv0(ph(i)%fit_mode,ph(i)%v,ph(i)%ffit_npol(k),ph(i)%ffit_apol(:,k))
          rms = sqrt(sum((ynew - yfit)**2) / real(ph(i)%nv,8))
          may = sum(abs(yfit)) / real(ph(i)%nv,8)
          write (uout,'("# F(V) fit: RMS(Ha) = ",1p,E18.10,"   avg-abs(Ha) = ",E18.10)') rms, may

          ! fit -T*S(V) as a function of temperature
          xfit = ph(i)%v(1:ph(i)%nv)
          yfit = -tlist(k) * ph(i)%dynamic_s(:,k)
          call fit_ev(ph(i)%sfit_mode,ph(i)%reg_mode,xfit,yfit,ph(i)%tsfit_npol(k),&
             ph(i)%tsfit_apol(:,k),ierr,.false.)
          if (ierr > 0) then
             write (uout,'(" T = ",F12.4," V = ",F12.4,"")') tlist(k), ph(i)%v(j)
             call error('dyneos','fit for -T*S(V;T) not found',faterr)
          end if
          ynew = fv0(ph(i)%sfit_mode,ph(i)%v,ph(i)%tsfit_npol(k),ph(i)%tsfit_apol(:,k))
          rms = sqrt(sum((ynew - yfit)**2) / real(ph(i)%nv,8))
          may = sum(abs(yfit)) / real(ph(i)%nv,8)
          write (uout,'("# -T*S(V) fit: RMS(Ha) = ",1p,E18.10,"   avg-abs(Ha) = ",E18.10)') rms, may

          ! fit CV(V) as a function of temperature
          xfit = ph(i)%v(1:ph(i)%nv)
          yfit = ph(i)%dynamic_cv(:,k)
          call fit_ev(ph(i)%cvfit_mode,ph(i)%reg_mode,xfit,yfit,ph(i)%cvfit_npol(k),&
             ph(i)%cvfit_apol(:,k),ierr,.false.)
          if (ierr > 0) then
             write (uout,'(" T = ",F12.4," V = ",F12.4,"")') tlist(k), ph(i)%v(j)
             call error('dyneos','fit for CV(V;T) not found',faterr)
          end if
          ynew = fv0(ph(i)%cvfit_mode,ph(i)%v,ph(i)%cvfit_npol(k),ph(i)%cvfit_apol(:,k))
          rms = sqrt(sum((ynew - yfit)**2) / real(ph(i)%nv,8))
          may = sum(abs(yfit)) / real(ph(i)%nv,8)
          write (uout,'("# CV(V) fit: RMS(Ha/K) = ",1p,E18.10,"   avg-abs(Ha/K) = ",E18.10)') rms, may
       end do
       deallocate(xfit,yfit,ynew,f0)

       ! allocate G(T,p), V(T,p) and B(T,p) arrays
       if (allocated(ph(i)%gtp)) deallocate(ph(i)%gtp)
       if (allocated(ph(i)%vtp)) deallocate(ph(i)%vtp)
       if (allocated(ph(i)%btp)) deallocate(ph(i)%btp)
       allocate(ph(i)%gtp(nts,nps),ph(i)%vtp(nts,nps),ph(i)%btp(nts,nps),ph(i)%didtp(nts,nps))
       ph(i)%vtp = -1d0
       ph(i)%vtp = -1d0
       ph(i)%gtp = 0d0
       ph(i)%didtp = .false.

       ! calculate the T,p properties
       do j = 1, nts
          do k = 1, nps
             ! find the equilibrium volume at (p,T)
             call fit_pshift(ph(i)%fit_mode,ph(i)%v,plist(k),ph(i)%ffit_npol(j),&
                ph(i)%ffit_apol(:,j),v,b,e,g,ierr)
             if (ierr > 0) cycle

             ! calculate the rest of the properties
             call dyneos_calc(ph(i),v,j,proplist)
             write (lu,fm) proplist

             ph(i)%didtp(j,k) = .true.
             ph(i)%gtp(j,k) = proplist(5) ! G(T,p)
             ph(i)%vtp(j,k) = proplist(3) ! V(T,p)
             ph(i)%btp(j,k) = proplist(9) ! BT(T,p)
          end do
          write (lu,'(/)')
       end do
    end do
    call fclose(lu)

    write (uout,*)

    ! generate gnu files for all properties
    call gen_allgnu_t(trim(fileroot)//"_all_t")
    call gen_allgnu_p(trim(fileroot)//"_all_p")

  end subroutine dyneos

  ! Calculate properties at a given volume and temperature.
  subroutine dyneos_calc(p,v,it,proplist)
    use evfunc, only: fv0, fv1, fv2, fv3, fv4
    use debye, only: tlim_gamma, debeins, thermalphon, thermal
    use varbas, only: mpropout, phase, vbracket, tlist, mm, vfree
    use tools, only: error
    use param, only: au2gpa, ha2kjmol, pckbau, pi, third
    type(phase), intent(in) :: p
    real*8, intent(in) :: v
    integer, intent(in) :: it
    real*8, intent(out) :: proplist(mpropout)

    real*8 :: t
    real*8 :: f0, f1, f2, f3, f4, b0, b1, b2, tmp
    real*8 :: f0s, f1s, f2s, f3s, g
    real*8 :: theta
    real*8 :: pbeta, alpha, cp, bs
    real*8 :: cv_lowt
    real*8 :: fvib, svib, uvib, cv_vib
    real*8 :: fel, sel, uel, cv_el
    real*8 :: fsum, ssum, usum, cv_sum
    real*8 :: gamma
    real*8 :: pext, psta, pth, dg, mtsvib

    ! temperature
    t = tlist(it)

    ! static energy and helmholtz free energy volume derivatives
    f0s = fv0(p%fit_mode,v,p%npol,p%cpol)
    f1s = fv1(p%fit_mode,v,p%npol,p%cpol)
    f2s = fv2(p%fit_mode,v,p%npol,p%cpol)
    f3s = fv3(p%fit_mode,v,p%npol,p%cpol)
    f0  = fv0(p%fit_mode,v,p%ffit_npol(it),p%ffit_apol(:,it))
    f1  = fv1(p%fit_mode,v,p%ffit_npol(it),p%ffit_apol(:,it))
    f2  = fv2(p%fit_mode,v,p%ffit_npol(it),p%ffit_apol(:,it))
    f3  = fv3(p%fit_mode,v,p%ffit_npol(it),p%ffit_apol(:,it))
    f4  = fv4(p%fit_mode,v,p%ffit_npol(it),p%ffit_apol(:,it))

    ! pressure derivatives of the isothermal bulk modulus
    b0 = v * f2
    b1 = -(1+v*f3/f2)
    b2 = ((f3+v*f4)/f2**2 - v*f3**2/f2**3)

    ! thermal pressure
    pext = -f1
    psta = -f1s
    pth = pext - psta

    ! Debye temperature
    theta = (6*pi*pi*vfree*v*v)**third / pckbau * p%pofunc * sqrt(f2s/mm)

    ! rest of properties
    fvib = f0 - f0s
    mtsvib = fv0(p%sfit_mode,v,p%tsfit_npol(it),p%tsfit_apol(:,it))
    if (t < tlim_gamma) then
       svib = 0d0
    else
       svib = mtsvib / (-t)
    end if
    uvib = fvib - mtsvib
    cv_vib = fv0(p%cvfit_mode,v,p%cvfit_npol(it),p%cvfit_apol(:,it))
    cv_lowt = cv_vib

    ! electronic contribution
    uel = 0d0
    fel = 0d0
    sel = 0d0
    cv_el = 0d0

    ! sum up
    fsum = fvib + fel
    usum = uvib + uel
    ssum = svib + sel
    cv_sum = cv_vib + cv_el

    ! bulk modulus
    b0 = b0 * au2gpa
    b2 = b2 / au2gpa

    ! gamma and related quantities
    if (t >= tlim_gamma .and. cv_sum > 1d-20) then
       ! non-zero temperature, usual calculation, no 0/0 here
       gamma = - v / cv_sum * fv1(p%sfit_mode,v,p%tsfit_npol(it),p%tsfit_apol(:,it)) / t
    else
       ! a meaningless gamma value - all properties other than
       ! gamma will be correct.
       gamma = 1d0
    end if

    ! gamma is finite at T->0, so no problem here, either
    Pbeta = cv_sum * gamma / v * au2gpa
    alpha = Pbeta / b0
    tmp = 1d0 + gamma * alpha * t
    Cp = cv_sum * tmp
    Bs = b0 * tmp

    ! Rest of thermodynamic properties -- note g can be calculated in
    ! 2 ways: with Fsum from two fits or from quasiharmonic
    ! formula.... dg is a measure this inaccuracy.
    ! pbeta = (dp/dT)_V
    g = f0 + pext * v
    dg = f0 - f0s - fsum

    ! conversion to output units
    g = g * ha2kjmol
    dg = dg * ha2kjmol
    pext = pext * au2gpa
    pth = pth * au2gpa
    psta = psta * au2gpa
    cp = cp * ha2kjmol * 1000
    alpha = alpha * 1d5
    fvib = fvib * ha2kjmol
    uvib = uvib * ha2kjmol
    svib = svib * ha2kjmol * 1000
    cv_vib = cv_vib * ha2kjmol * 1000
    fel = fel * ha2kjmol
    uel = uel * ha2kjmol
    sel = sel * ha2kjmol * 1000
    cv_el = cv_el * ha2kjmol * 1000
    fsum = fsum * ha2kjmol
    usum = usum * ha2kjmol
    ssum = ssum * ha2kjmol * 1000
    cv_sum = cv_sum * ha2kjmol * 1000

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

  ! Write the dgtp file.
  subroutine deltag()
    use varbas, only: nph, ph, nps, plist, nts, tlist
    use tools, only: leng, fopen, fclose
    use param, only: mline, uout, fileroot, iowrite, null
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

  ! Write the tpstab file.
  subroutine stablevbg()
    use varbas, only: nph, ph, nps, plist, nts, tlist
    use tools, only: leng, fopen, fclose
    use param, only: uout, fileroot, iowrite, null
    real*8 :: gmin
    integer :: lu, idmin
    integer :: i, j, k
    integer :: ndo, imask(nph)

    call mask_trans(ndo,imask)
    if (ndo < 2) return

    write (uout,'("* G(T,p), V(T,p) and B(T,p) of the stable phase")')
    write (uout,'("  Writing file : ",A/)') trim(fileroot)//".tpstab"
    lu = fopen(lu,trim(fileroot)//".tpstab"//null,iowrite)

    write (lu,'("## 1:T(K) 2:p(GPa) 3:stable phase 4:G(kJ/mol) 5:V(bohr^3) 6:B_T(GPa)")')
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

  ! Interpolate the input data as requested by the user.
  subroutine interpolate()
    use evfunc, only: fv1
    use fit, only: fit_pshift
    use varbas, only: nph, ph, nps, plist, nts, tlist, writelevel, tm_static,&
       tm_debye_extended, tm_qhafull
    use tools, only: error, leng, fopen, fclose, realloc
    use param, only: mline_fmt, uout, format_string, format_string_header, fileroot, iowrite, &
       warning, faterr, null, au2gpa, ifmt_v, ifmt_p, ifmt_t, ifmt_interp
    integer :: i, j, k
    integer :: mint, n, lu, it
    real*8 :: v, p, t, e, b, g, fac, f1
    real*8, allocatable :: fi(:)
    integer, allocatable :: ifm(:)
    character*(mline_fmt) :: fm1, fms
    integer :: napol, ierr, imode
    real*8 :: apol(0:mmpar)
    real*8, allocatable :: xfit(:), yfit(:)
    integer :: ffit_npol
    real*8, allocatable :: ffit_apol(:)

    if (writelevel < 2) return

    if (interp_input > 0) then
       if (.not.allocated(fint)) allocate(fint(mxint))
       if (.not.allocated(iint)) allocate(iint(mxint))
       do i = 1, nps
          if (interp_input > 1) then
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

             ! check that this is a known temperature
             it = 0
             do k = 1, nts
                if (abs(tlist(k) - t) < 1d-3) then
                   it = k
                   exit
                end if
             end do
             if (it == 0) &
                call error('interpolate','Interpolation temperature must be included in the temperature list',faterr)

             ! get the volume
             call fit_pshift(ph(i)%fit_mode,ph(i)%v,p,ph(i)%ffit_npol(it),&
                ph(i)%ffit_apol(:,it),v,b,e,g,ierr)
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

  ! Apply the empirical energy correction.
  subroutine eshift_vexp()
    use evfunc, only: fv0, fv1, fv2
    use debye, only: fill_thetad, thermal_debye_extended, thermal_qha, thermal
    use fit, only: fit_ev, fit_pshift
    use varbas, only: nph, ph, scal_bpscal, scal_apbaf, scal_pshift, scal_use, scal_noscal, &
       tm_debyegrun, tm_debye, tm_debye_einstein, tm_debye_einstein_v, phase_checkfiterr,&
       nts, tlist, tm_debye_extended, tm_qhafull, tm_debye_input, tm_debye_poisson_input
    use tools, only: error, leng
    use param, only: mline, uout, warning, faterr, au2gpa
    real*8, parameter :: facprec = 1d-10

    integer :: i, j, ierr
    integer :: napol
    real*8 :: apol(0:mmpar), psum
    character*(mline) :: msg
    integer :: imode, niter, it, k
    real*8 :: vexpt, bexpt, fac
    real*8 :: v0, b0, e0, g0, e0new, g0new, v0t, b0t, e0t, g0t
    real*8 :: vexp, bexp
    real*8 :: fa, fb, qfa, qfb, qfx
    real*8 :: psta_vexpt, pth_vexpt, bsta_vexpt
    real*8 :: bt_vexpt, bpobj
    real*8 :: vold, dum
    real*8, allocatable :: xfit(:), yfit(:)

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

       ! reserve space to run the fit for the EEC temperature
       if (allocated(xfit)) deallocate(xfit)
       if (allocated(yfit)) deallocate(yfit)
       allocate(xfit(ph(i)%nv),yfit(ph(i)%nv))
       xfit = ph(i)%v(1:ph(i)%nv)

       ! calculate thermodynamic properties on the (V,T) grid
       do j = 1, ph(i)%nv
          if (ph(i)%tmodel == tm_debye_extended) then
             call thermal_debye_extended(ph(i),ph(i)%eec_t,j,yfit(j),dum,dum)
          elseif (ph(i)%tmodel == tm_qhafull) then
             call thermal_qha(ph(i),ph(i)%eec_t,j,yfit(j),dum,dum)
          elseif (ph(i)%tmodel == tm_debye.or.ph(i)%tmodel == tm_debye_input.or.&
             ph(i)%tmodel == tm_debyegrun.or.ph(i)%tmodel == tm_debye_poisson_input) then
             call thermal(ph(i)%td(j),ph(i)%eec_t,dum,dum,dum,dum,yfit(j),dum)
          else
             write (*,*) "fixme!"
             stop 1
          end if
       end do
       yfit = ph(i)%e + yfit

       ! fit the F(V) as a function of temperature
       call fit_ev(ph(i)%fit_mode,ph(i)%reg_mode,xfit,yfit,napol,apol,ierr,.false.)
       if (ierr > 0) then
          write (uout,'(" T = ",F12.4," V = ",F12.4,"")') tlist(k), ph(i)%v(j)
          call error('dyneos','fit for F(V;T) not found',faterr)
       end if
       call fit_pshift(ph(i)%fit_mode,ph(i)%v,ph(i)%eec_p,napol,apol(0:napol),v0t,b0t,e0t,g0t,ierr)
       imode = ph(i)%fit_mode
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
          ph(i)%tmodel == tm_debye_einstein .or. ph(i)%tmodel == tm_debye_einstein_v .or.&
          ph(i)%tmodel == tm_debyegrun)) goto 1

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

  !< Returns fitting parameters of Cv_el(V,T).
  subroutine fit_cvelgrid_t(p,T,ncvpol,cvpol,mode,ierr)
    use evfunc, only: fv1
    use varbas, only: phase, em_pol4, ftsel_fitmode
    use fit, only: fit_ev
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

  ! Check wether phase p is calculable at point v.
  function is_dyn_v_in(p,v)
    use varbas, only: phase
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

  ! mask to elimiate phases for which dynamical properties cannot be calculated.
  subroutine mask_trans(n,imask)
    use varbas, only: nph, ph, tm_static, dotrans, writelevel
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

  ! Print Debye-Einstein model frequencies to the gammafreq file.
  subroutine printfreqs()
    use evfunc, only: fv1, fv2
    use varbas, only: nph, ph, tm_debye_einstein
    use tools, only: leng, fopen, fclose
    use param, only: uout, fileroot, iowrite, null, au2gpa, ha2cm_1, half, twothird
    integer :: i, j
    real*8 :: vol, pstat, bstat, veq, beq
    real*8 :: tmpvol, tmpb, tmppb, tmptot
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

             freq(1:n) = ph(i)%freqg(:,1) * tmpTot * ha2cm_1

             write (lu,'(F10.4,9999(F16.8))') ph(i)%v(j), freq(1:n)
          end do
          write (lu,'(/)')
       end if
    end do

    call fclose(lu)
    deallocate(freq)

  end subroutine printfreqs

  ! Write the edat plot files.
  subroutine write_energy(ini,step,end)
    use evfunc, only: fv0
    use varbas, only: nph, ph
    use tools, only: leng, fopen, fclose
    use param, only: mline, uout, bohr2angstrom, fileroot, ha2ev, iowrite, null, &
       units_e_ev, units_e_ry, units_v_ang3
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
             efac = ha2ev
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
             efac = ha2ev
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
