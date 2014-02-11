!  *******************************************************************  !
!  I-Disorder Fortran Code 2007-2014                                    !
!                                                                       !
!  Written by Alexandre Reily Rocha (reilya@ift.unesp.br),              !
!             Pedro Brandimarte (brandimarte@gmail.com) and             !
!             Alberto Torres (alberto.trj@gmail.com).                   !
!                                                                       !
!  Copyright (c), All Rights Reserved                                   !
!                                                                       !
!  This program is free software. You can redistribute it and/or        !
!  modify it under the terms of the GNU General Public License          !
!  (version 3 or later) as published by the Free Software Foundation    !
!  <http://fsf.org/>.                                                   !
!                                                                       !
!  This program is distributed in the hope that it will be useful, but  !
!  WITHOUT ANY WARRANTY, without even the implied warranty of           !
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU     !
!  General Public License for more details (file 'LICENSE_GPL'          !
!  distributed along with this program or at                            !
!  <http://www.gnu.org/licenses/gpl.html>).                             !
!  *******************************************************************  !
!                           MODULE idsrdr_out                           !
!  *******************************************************************  !
!  Description: contains subroutines for writing calculated values,     !
!  which are distributed over the nodes, at output files.               !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *******************************************************************  !

MODULE idsrdr_out

!
!   Modules
!
  use parallel,        only: 
  use idsrdr_options,  only: 
  use idsrdr_engrid,   only: 
  use idsrdr_units,    only: 
  use idsrdr_spectral, only: 
  use idsrdr_current,  only: 
  use idsrdr_io,       only: 
  use idsrdr_string,   only: 
  use idsrdr_power,    only: 
  use idsrdr_conduct,  only: 
 
  implicit none

  PUBLIC  :: output
  PRIVATE :: writespectral, writecurrent, writepower,                   &
             writedIdV, writed2IdV2


CONTAINS


!  *******************************************************************  !
!                                output                                 !
!  *******************************************************************  !
!  Description: interface subroutine for calling output writing         !
!  subroutines.                                                         !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode                 : True if it is the I/O node          !
!  integer NIVP                   : Number of bias potential points     !
!  *******************************************************************  !
  subroutine output

!
!   Modules
!
    use parallel,        only: IOnode
    use idsrdr_options,  only: NIVP

    if (IOnode) write (6,'(/,28("*"),a,29("*"))')                       &
            ' Writing output files '

!   Write spectral function and DOS to files.
    call writespectral

!   Write calculated currents to files.
    call writecurrent

!   Write calculated dissipated powers to files.
    call writepower

    if (NIVP < 3) return

!   Write calculated differential conductance ('dI/dV') to files.
    call writedIdV

!   Write derivative of differential conductance ('d2I/dV2') to files.
    call writed2IdV2


  end subroutine output


!  *******************************************************************  !
!                             writespectral                             !
!  *******************************************************************  !
!  Description: writes calculated spectral function and density of      !
!  states to output files ('slabel_unit.SPCTR' and 'slabel_unit.DOS').  !
!  The files indexed with 'units+1' contains the summation from all     !
!  units with e-ph.                                                     !
!                                                                       !
!  Written by Pedro Brandimarte, Nov 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    November 2013                                   !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Node                : Actual node (MPI_Comm_rank)            !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  integer nspin               : Number of spin components              !
!  integer label_length        : Length of system label                 !
!  character(label_length) slabel : System Label (for output files)     !
!  character(60) directory     : Working directory                      !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  real*8 Ei(NTenerg_div)      : Energy grid points                     !
!  integer nunitseph           : Number of units with eph               !
!  real*8 spctrl(NTenerg_div,nspin,nunitseph+1) : Spectral function     !
!  real*8 dos(NTenerg_div,nspin,nunitseph+1)    : Density of states     !
!  *******************************************************************  !
  subroutine writespectral

!
!   Modules
!
#ifdef MASTER_SLAVE
    use parallel,        only: IOnode, Node, Nodes, MPI_Comm_MyWorld
#else
    use parallel,        only: IOnode, Node, Nodes
#endif
    use idsrdr_options,  only: nspin, label_length, slabel, directory
    use idsrdr_engrid,   only: NTenerg_div, Ei
    use idsrdr_units,    only: nunitseph
    use idsrdr_spectral, only: spctrl, dos
    use idsrdr_io,       only: IOassign, IOclose
    use idsrdr_string,   only: STRconcat, STRpaste

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: J, e, s, iuSpc, iuDos
    real(8), parameter :: pi = 3.14159265358979323846264338327950241D0
    real(8), dimension(:), allocatable :: buffEn
    real(8), dimension(:,:), allocatable :: buffSpc, buffDos
    character(len=10) :: suffix
    character(len=label_length+70) :: fnSpc, fnDos
#ifdef MPI
#ifndef MASTER_SLAVE
    integer :: n
    integer, dimension(MPI_Status_Size) :: MPIstatus
#endif
    integer :: MPIerror
#endif

    do J = 1,nunitseph+1 ! over units with e-ph intereaction

       if (IOnode) then

!         Set file's names and open it.
          write (suffix,'(i3)') J
          call STRconcat (suffix, '.SPCTR', suffix)
          call STRpaste ('_', suffix, suffix)
          call STRpaste (slabel, suffix, fnSpc)
          call STRpaste (directory, fnSpc, fnSpc)
          write (suffix,'(i3)') J
          call STRconcat (suffix, '.DOS', suffix)
          call STRpaste ('_', suffix, suffix)
          call STRpaste (slabel, suffix, fnDos)
          call STRpaste (directory, fnDos, fnDos)
          call IOassign (iuSpc)
          open (iuSpc, FILE=fnSpc, FORM='FORMATTED', STATUS='REPLACE')
          call IOassign (iuDos)
          open (iuDos, FILE=fnDos, FORM='FORMATTED', STATUS='REPLACE')

          write (6, '(/,a,i3,a,a)', advance='no')                       &
               'Writing spectral function ', J, ' to: ', trim(fnSpc)
          write (6, '(/,a,i3,a,a)') 'Writing DOS ', J, ' to: ',         &
               trim(fnDos)

!         Allocate buffers arrays.
          allocate (buffEn(NTenerg_div))
          allocate (buffSpc(NTenerg_div,nspin))
          allocate (buffDos(NTenerg_div,nspin))

#ifdef MASTER_SLAVE
       else
          allocate (buffSpc(NTenerg_div,nspin))
          allocate (buffDos(NTenerg_div,nspin))
#endif

       endif

!      Write to the output files (energies in eV (from CODATA - 2012)).
#ifdef MASTER_SLAVE
!      Reduce the reults to the IOnode and write
       call MPI_Reduce (spctrl(1,1,J), buffSpc, NTenerg_div*nspin,      &
                        MPI_Double_Precision, MPI_Sum, 0,               &
                        MPI_Comm_MyWorld, MPIerror)
       call MPI_Reduce (dos(1,1,J), buffDos, NTenerg_div*nspin,         &
                        MPI_Double_Precision, MPI_Sum, 0,               &
                        MPI_Comm_MyWorld, MPIerror)

       if (IOnode) then
          spctrl(:,:,J) = buffSpc(:,:)
          dos(:,:,J)    = buffDos(:,:)
#elif defined MPI
       do n = 0,Nodes-1
          if (Node == n .and. IOnode) then
#endif
             do e = 1,NTenerg_div

!               'Ei' in eV (from CODATA - 2012).
                write (iuSpc, '(/,e17.8e3)', advance='no')              &
                     Ei(e) * 13.60569253D0
                write (iuDos, '(/,e17.8e3)', advance='no')              &
                     Ei(e) * 13.60569253D0

                do s = 1,nspin
!                  'spctrl' and 'dos' in 1/eV (from CODATA - 2012).
                   write (iuSpc, '(e17.8e3)', advance='no')             &
                        spctrl(e,s,J) / (13.60569253D0 * pi)
                   write (iuDos, '(e17.8e3)', advance='no')             &
                        dos(e,s,J) / (13.60569253D0 * pi)
                enddo
             enddo
#ifdef MASTER_SLAVE
       end if ! IOnode
#elif defined MPI
          elseif (Node == n) then
             call MPI_Send (Ei, NTenerg_div, MPI_Double_Precision,      &
                            0, 1, MPI_Comm_world, MPIerror)
             call MPI_Send (spctrl(1,1,J), NTenerg_div*nspin,           &
                            MPI_Double_Precision, 0, 2,                 &
                            MPI_Comm_world, MPIerror)
             call MPI_Send (dos(1,1,J), NTenerg_div*nspin,              &
                            MPI_Double_Precision, 0, 3,                 &
                            MPI_Comm_world, MPIerror)
          elseif (Node == 0) then
             call MPI_Recv (buffEn, NTenerg_div, MPI_Double_Precision,  &
                            n, 1, MPI_Comm_world, MPIstatus, MPIerror)
             call MPI_Recv (buffSpc, NTenerg_div*nspin,                 &
                            MPI_Double_Precision, n, 2,                 &
                            MPI_Comm_world, MPIstatus, MPIerror)
             call MPI_Recv (buffDos, NTenerg_div*nspin,                 &
                            MPI_Double_Precision, n, 3,                 &
                            MPI_Comm_world, MPIstatus, MPIerror)
          endif
          if (n /= 0) then
             call MPI_Barrier (MPI_Comm_world, MPIerror)
             if (Node == 0) then
                do e = 1,NTenerg_div

!                  'Ei' in eV (from CODATA - 2012).
                   write (iuSpc, '(/,e17.8e3)', advance='no')           &
                        buffEn(e) * 13.60569253D0
                   write (iuDos, '(/,e17.8e3)', advance='no')           &
                        buffEn(e) * 13.60569253D0

                   do s = 1,nspin

!                     'spctrl' and 'dos' in 1/eV (from CODATA - 2012).
                      write (iuSpc, '(e17.8e3)', advance='no')          &
                           buffSpc(e,s) / (13.60569253D0 * pi)
                      write (iuDos, '(e17.8e3)', advance='no')          &
                           buffDos(e,s) / (13.60569253D0 * pi)
                   enddo
                enddo
             endif
          endif
       enddo ! n = 0,Nodes-1
#endif

!      Close files and free buffers memory.
       if (IONode) then

          call IOclose (iuSpc)
          call IOclose (iuDos)

          deallocate (buffEn)
          deallocate (buffSpc)
          deallocate (buffDos)

#ifdef MASTER_SLAVE
       else
          deallocate (buffSpc)
          deallocate (buffDos)
#endif

       endif

    enddo ! do J = 1,nunitseph+1


  end subroutine writespectral


!  *******************************************************************  !
!                             writecurrent                              !
!  *******************************************************************  !
!  Description: writes calculated currents to output file as follows    !
!                                                                       !
!    - 'slabel_ExVxI.CUR' contains 'Ei Vbias Iel Isymm Iasymm Itot'     !
!                                                                       !
!  where 'Itot = Iel+Isymm+Iasymm'.                                     !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Node                : Actual node (MPI_Comm_rank)            !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  integer nspin               : Number of spin components              !
!  integer label_length        : Length of system label                 !
!  character(label_length) slabel : System Label (for output files)     !
!  character(60) directory     : Working directory                      !
!  integer NIVP                : Number of bias potential points        !
!  real*8 VInitial             : Initial value of the bias potential    !
!  real*8 dV                   : Bias potential step                    !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  real*8 Ei(NTenerg_div)      : Energy grid points                     !
!  TYPE(calcCurr) allcurr(NTenerg_div,nspin,NIVP) : [real] calculated   !
!                                                   currents            !
!  *******************************************************************  !
  subroutine writecurrent

!
!   Modules
!
#ifdef MASTER_SLAVE
    use parallel,        only: IOnode, Node, Nodes, MPI_Comm_MyWorld
#else
    use parallel,        only: IOnode, Node, Nodes
#endif
    use idsrdr_options,  only: nspin, label_length, slabel, directory,  &
                               NIVP, VInitial, dV
    use idsrdr_engrid,   only: NTenerg_div, Ei
#ifdef MASTER_SLAVE
    use idsrdr_current,  only: calcCurr, allcurr, sumCalcCurr
#else
    use idsrdr_current,  only: calcCurr, allcurr
#endif
    use idsrdr_io,       only: IOassign, IOclose
    use idsrdr_string,   only: STRpaste

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: e, s, v, iuExVxI
    real(8) :: Vbias
    real(8), dimension(:), allocatable :: buffEn
    TYPE(calcCurr), allocatable, dimension (:,:,:) :: buffCurr
    character(len=11) :: suffix
    character(len=label_length+70) :: fExVxI
#ifdef MPI
#ifndef MASTER_SLAVE
    integer :: n
    integer, dimension(MPI_Status_Size) :: MPIstatus
#else
    integer :: MPIop
#endif
    integer :: MPIerror, MPIcalcCurr
    integer :: blocklens(1) ! # of elements in each block
    integer :: blockdispl(1) ! byte displacement of each block
    integer :: oldtypes(1) ! type of elements in each block
#endif

    if (IOnode) then

       write (6, '(/,a)') 'Writing calculated currents to *.CUR files'

!      Set file's names.
       suffix = '_ExVxI.CUR'
       call STRpaste (slabel, suffix, fExVxI)
       call STRpaste (directory, fExVxI, fExVxI)

!      Open them.
       call IOassign (iuExVxI)
       open (iuExVxI, FILE=fExVxI, FORM='FORMATTED', STATUS='REPLACE')

!      Allocate buffers arrays.
       allocate (buffEn(NTenerg_div))
       allocate (buffCurr(NTenerg_div,nspin,NIVP))

#ifdef MASTER_SLAVE
    else
       allocate (buffCurr(NTenerg_div,nspin,NIVP))
#endif

    endif

#ifdef MPI
!   Set the description of 'calcCurr' type (3 doubles).
    blocklens(1) = 3
    blockdispl(1) = 0
    oldtypes(1) = MPI_Double_Precision

!   Define structured type and commit it.
    call MPI_Type_create_struct (1, blocklens, blockdispl, oldtypes,    &
                                 MPIcalcCurr, MPIerror)
    call MPI_Type_commit (MPIcalcCurr, MPIerror)
#endif

!   Write to the output files.
#ifdef MASTER_SLAVE
!   Define a "MPI_SUM" operator for type 'calcCurr'.
    call MPI_Op_create (sumCalcCurr, .true., MPIop, MPIerror)

!   Reduce the reults to the IOnode and write
    call MPI_Reduce (allcurr(1,1,1), buffCurr,                          &
                     NTenerg_div*nspin*NIVP, MPIcalcCurr, MPIop,        &
                     0, MPI_Comm_MyWorld, MPIerror)

    if (IOnode) then
       allcurr(:,:,:) = buffCurr(:,:,:)
#elif defined MPI
    do n = 0,Nodes-1
       if (Node == n .and. IOnode) then
#endif
          do e = 1,NTenerg_div
             do s = 1,nspin

                Vbias = VInitial

                do v = 1,NIVP

!                  'Ei' and 'Vbias' in eV (from CODATA - 2012).
                   write (iuExVxI, '(e17.8e3,e17.8e3,e17.8e3,'      //  &
                        'e17.8e3,e17.8e3,e17.8e3)')                     &
                        Ei(e)*13.60569253D0, Vbias*13.60569253D0,       &
                        allcurr(e,s,v)%el, allcurr(e,s,v)%isymm,        &
                        allcurr(e,s,v)%iasymm, allcurr(e,s,v)%el        &
                        + allcurr(e,s,v)%isymm + allcurr(e,s,v)%iasymm

                   Vbias = Vbias + dV

                enddo

                write (iuExVxI, '(/)', advance='no')

             enddo
          enddo
#ifdef MASTER_SLAVE
    end if ! IOnode
#elif defined MPI
       elseif (Node == n) then
          call MPI_Send (Ei, NTenerg_div, MPI_Double_Precision,         &
                         0, 1, MPI_Comm_world, MPIerror)
          call MPI_Send (allcurr(1,1,1), NTenerg_div*nspin*NIVP,        &
                         MPIcalcCurr, 0, 2, MPI_Comm_world, MPIerror)
       elseif (Node == 0) then
          call MPI_Recv (buffEn, NTenerg_div, MPI_Double_Precision,     &
                         n, 1, MPI_Comm_world, MPIstatus, MPIerror)
          call MPI_Recv (buffCurr, NTenerg_div*nspin*NIVP, MPIcalcCurr, &
                         n, 2, MPI_Comm_world, MPIstatus, MPIerror)
       endif
       if (n /= 0) then
          call MPI_Barrier (MPI_Comm_world, MPIerror)
          if (Node == 0) then
             do e = 1,NTenerg_div
                do s = 1,nspin

                   Vbias = VInitial

                   do v = 1,NIVP

!                     'Ei' and 'Vbias' in eV (from CODATA - 2012).
                      write (iuExVxI, '(e17.8e3,e17.8e3,e17.8e3,'   //  &
                           'e17.8e3,e17.8e3,e17.8e3)')                  &
                           buffEn(e)*13.60569253D0, Vbias*13.60569253D0,&
                           buffCurr(e,s,v)%el, buffCurr(e,s,v)%isymm,   &
                           buffCurr(e,s,v)%iasymm, buffCurr(e,s,v)%el   &
                           + buffCurr(e,s,v)%isymm                      &
                           + buffCurr(e,s,v)%iasymm

                      Vbias = Vbias + dV

                   enddo

                   write (iuExVxI, '(/)', advance='no')

                enddo
             enddo
          endif
       endif
    enddo ! n = 0,Nodes-1
#endif

!   Close files and free buffers memory.
    if (IOnode) then

       call IOclose (iuExVxI)

       deallocate (buffEn)
       deallocate (buffCurr)

#ifdef MASTER_SLAVE
    else
       deallocate (buffCurr)
       call MPI_Op_free (MPIop, MPIerror)
#endif

    endif

!   Free MPI data type.
#ifdef MPI
    call MPI_Type_free (MPIcalcCurr, MPIerror)
#endif


  end subroutine writecurrent


!  *******************************************************************  !
!                              writepower                               !
!  *******************************************************************  !
!  Description: writes calculated dissipated powers to output file as   !
!                                                                       !
!    - 'slabel_ExVxP.PWR' contains 'Ei Vbias Power_mode Power_tot'      !
!                                                                       !
!  where 'Power_tot = Sum_mode Power_mode'.                             !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Node                : Actual node (MPI_Comm_rank)            !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  integer nspin               : Number of spin components              !
!  integer label_length        : Length of system label                 !
!  character(label_length) slabel : System Label (for output files)     !
!  character(60) directory     : Working directory                      !
!  integer NIVP                : Number of bias potential points        !
!  real*8 VInitial             : Initial value of the bias potential    !
!  real*8 dV                   : Bias potential step                    !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  real*8 Ei(NTenerg_div)      : Energy grid points                     !
!  integer nModes(neph)        : Number of vibrational modes            !
!  integer ephIdx(ntypeunits+2): Unit index (those with e-ph)           !
!  integer nunitseph           : Number of units with eph               !
!  integer eph_type(nunitseph) : Units types with eph                   !
!  TYPE(phonon) phPwr(NTenerg_div,nspin,NIVP,nunitseph) : [real]        !
!                                                    calculated powers  !
!  *******************************************************************  !
  subroutine writepower

!
!   Modules
!
#ifdef MASTER_SLAVE
    use parallel,        only: IOnode, Node, Nodes, MPI_Comm_MyWorld
#else
    use parallel,        only: IOnode, Node, Nodes
#endif
    use idsrdr_options,  only: nspin, label_length, slabel, directory,  &
                               NIVP, VInitial, dV
    use idsrdr_engrid,   only: NTenerg_div, Ei
    use idsrdr_ephcoupl, only: nModes, ephIdx
    use idsrdr_units,    only: nunitseph, eph_type
    use idsrdr_power,    only: phPwr
    use idsrdr_io,       only: IOassign, IOclose
    use idsrdr_string,   only: STRconcat, STRpaste

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: e, s, v, u, w, idx, iuExVxP
    real(8) :: Vbias
    real(8), dimension(:), allocatable :: buffEn
    real(8), dimension (:,:,:,:), allocatable :: buffPwr
    character(len=15) :: suffix
    character(len=label_length+70) :: fExVxP

#ifdef MPI
#ifndef MASTER_SLAVE
    integer :: n
    integer, dimension(MPI_Status_Size) :: MPIstatus
#endif
    integer :: MPIerror
#endif

    if (IOnode .and. nunitseph > 0) write (6, '(/,a)')                  &
         'Writing calculated dissipated powers to *.PWR files'

    do u = 1,nunitseph ! over units with e-ph intereaction

       idx = ephIdx(eph_type(u))

       if (IOnode) then

!         Set file's names.
          write (suffix,'(i3)') u
          call STRconcat (suffix, '.PWR', suffix)
          call STRpaste ('_ExVxP_', suffix, suffix)
          call STRpaste (slabel, suffix, fExVxP)
          call STRpaste (directory, fExVxP, fExVxP)

!         Open them.
          call IOassign (iuExVxP)
          open (iuExVxP, FILE=fExVxP, FORM='FORMATTED', STATUS='REPLACE')

!         Allocate buffers arrays.
          allocate (buffEn(NTenerg_div))
          allocate (buffPwr(NTenerg_div,nspin,NIVP,nModes(idx)))

#ifdef MASTER_SLAVE
       else
          allocate (buffPwr(NTenerg_div,nspin,NIVP,nModes(idx)))
#endif

       endif

!      Write to the output files.
#ifdef MASTER_SLAVE
!      Reduce the reults to the IOnode and write
       call MPI_Reduce (phPwr(u)%P(1,1,1,1), buffPwr,                   &
                        NTenerg_div*nspin*NIVP*nModes(idx),             &
                        MPI_Double_Precision, MPI_Sum,                  &
                        0, MPI_Comm_MyWorld, MPIerror)

       if (IOnode) then
          phPwr(u)%P(:,:,:,:) = buffPwr(:,:,:,:)
#elif defined MPI
       do n = 0,Nodes-1
          if (Node == n .and. IONode) then
#endif
             do e = 1,NTenerg_div
                do s = 1,nspin

                   Vbias = VInitial

                   do v = 1,NIVP

!                     'phPwr' in eV/s and, 'Ei' and 'Vbias'
!                     in eV (from CODATA - 2012).
                      write (iuExVxP, '(/,e17.8e3,e17.8e3)',            &
                           advance='no') Ei(e)*13.60569253D0,           &
                           Vbias*13.60569253D0
                      do w = 1,nModes(idx)
                         write (iuExVxP, '(e17.8e3)', advance='no')     &
                              phPwr(u)%P(e,s,v,w)*13.60569253D0
                      enddo
                      write (iuExVxP, '(e17.8e3)')                      &
                           SUM(phPwr(u)%P(e,s,v,1:nModes(idx)))         &
                           *13.60569253D0

                      Vbias = Vbias + dV

                   enddo

                   write (iuExVxP, '(/)', advance='no')

                enddo
             enddo
#ifdef MASTER_SLAVE
       end if ! IOnode
#elif defined MPI
          elseif (Node == n) then
             call MPI_Send (Ei, NTenerg_div, MPI_Double_Precision,      &
                            0, 1, MPI_Comm_world, MPIerror)
             call MPI_Send (phPwr(u)%P(1,1,1,1),                        &
                            NTenerg_div*nspin*NIVP*nModes(idx),         &
                            MPI_Double_Precision, 0, 2,                 &
                            MPI_Comm_world, MPIerror)
          elseif (Node == 0) then
             call MPI_Recv (buffEn, NTenerg_div, MPI_Double_Precision,  &
                            n, 1, MPI_Comm_world, MPIstatus, MPIerror)
             call MPI_Recv (buffPwr,                                    &
                            NTenerg_div*nspin*NIVP*nModes(idx),         &
                            MPI_Double_Precision, n, 2,                 &
                            MPI_Comm_world, MPIstatus, MPIerror)
          endif
          if (n /= 0) then
             call MPI_Barrier (MPI_Comm_world, MPIerror)
             if (Node == 0) then
                do e = 1,NTenerg_div
                   do s = 1,nspin

                      Vbias = VInitial

                      do v = 1,NIVP

!                        'phPwr' in eV/s and, 'Ei' and 'Vbias'
!                        in eV (from CODATA - 2012).
                         write (iuExVxP, '(/,e17.8e3,e17.8e3)',         &
                              advance='no') Vbias*13.60569253D0,        &
                              buffEn(e)*13.60569253D0
                         do w = 1,nModes(idx)
                            write (iuExVxP, '(e17.8e3)', advance='no')  &
                                 buffPwr(e,s,v,w)*13.60569253D0
                         enddo
                         write (iuExVxP, '(e17.8e3)')                   &
                              SUM(buffPwr(e,s,v,1:nModes(idx)))         &
                              *13.60569253D0

                         Vbias = Vbias + dV

                      enddo

                      write (iuExVxP, '(/)', advance='no')

                   enddo
                enddo
             endif
          endif
       enddo ! n = 0,Nodes-1
#endif

!      Close files and free buffers memory.
       if (IOnode) then

          call IOclose (iuExVxP)

          deallocate (buffEn)
          deallocate (buffPwr)

#ifdef MASTER_SLAVE
       else
          deallocate (buffPwr)
#endif

       endif

    enddo ! do u = 1,nunitseph


  end subroutine writepower


!  *******************************************************************  !
!                               writedIdV                               !
!  *******************************************************************  !
!  Description: writes calculated differential conductances ('dI/dV')   !
!  to output file as follows                                            !
!                                                                       !
!    - 'slabel_ExVxdI.dIdV' contains 'Ei Vbias dIdVel dIdVsymm          !
!    dIdVasymm dIdVtot'                                                 !
!                                                                       !
!  where 'dIdVtot = dIdVel+dIdVsymm+dIdVasymm'.                         !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Node                : Actual node (MPI_Comm_rank)            !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  integer nspin               : Number of spin components              !
!  integer label_length        : Length of system label                 !
!  character(label_length) slabel : System Label (for output files)     !
!  character(60) directory     : Working directory                      !
!  integer NIVP                : Number of bias potential points        !
!  real*8 VInitial             : Initial value of the bias potential    !
!  real*8 dV                   : Bias potential step                    !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  real*8 Ei(NTenerg_div)      : Energy grid points                     !
!  TYPE(alldIdV) dIdV(NTenerg_div,nspin,NIVP) : [real] calculated       !
!                                               dI/dV                   !
!  *******************************************************************  !
  subroutine writedIdV

!
!   Modules
!
#ifdef MASTER_SLAVE
    use parallel,        only: IOnode, Node, Nodes, MPI_Comm_MyWorld
#else
    use parallel,        only: IOnode, Node, Nodes
#endif
    use idsrdr_options,  only: nspin, label_length, slabel, directory,  &
                               NIVP, VInitial, dV
    use idsrdr_engrid,   only: NTenerg_div, Ei
#ifdef MASTER_SLAVE
    use idsrdr_conduct,  only: alldIdV, dIdV, sumAlldIdV
#else
    use idsrdr_conduct,  only: alldIdV, dIdV
#endif
    use idsrdr_io,       only: IOassign, IOclose
    use idsrdr_string,   only: STRpaste

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: e, s, v, iuExVxdI
    real(8) :: Vbias
    real(8), dimension(:), allocatable :: buffEn
    TYPE(alldIdV), allocatable, dimension (:,:,:) :: buffdIdV
    character(len=13) :: suffix
    character(len=label_length+70) :: fExVxdI
#ifdef MPI
#ifndef MASTER_SLAVE
    integer :: n
    integer, dimension(MPI_Status_Size) :: MPIstatus
#else
    integer :: MPIop
#endif
    integer :: MPIerror, MPIalldIdV
    integer :: blocklens(1) ! # of elements in each block
    integer :: blockdispl(1) ! byte displacement of each block
    integer :: oldtypes(1) ! type of elements in each block
#endif

    if (IOnode) then

       write (6, '(/,a)') 'Writing calculated differential '       //   &
            'conductance to *.dIdV files'

!      Set file's names.
       suffix = '_ExVxdI.dIdV'
       call STRpaste (slabel, suffix, fExVxdI)
       call STRpaste (directory, fExVxdI, fExVxdI)

!      Open them.
       call IOassign (iuExVxdI)
       open (iuExVxdI, FILE=fExVxdI, FORM='FORMATTED', STATUS='REPLACE')

!      Allocate buffers arrays.
       allocate (buffEn(NTenerg_div))
       allocate (buffdIdV(NTenerg_div,nspin,NIVP))

#ifdef MASTER_SLAVE
    else
       allocate (buffdIdV(NTenerg_div,nspin,NIVP))
#endif

    endif

#ifdef MPI
!   Set the description of 'alldIdV' type (3 doubles).
    blocklens(1) = 3
    blockdispl(1) = 0
    oldtypes(1) = MPI_Double_Precision

!   Define structured type and commit it.
    call MPI_Type_create_struct (1, blocklens, blockdispl, oldtypes,    &
                                 MPIalldIdV, MPIerror)
    call MPI_Type_commit (MPIalldIdV, MPIerror)
#endif

!   Write to the output files.
#ifdef MASTER_SLAVE
!   Define a "MPI_SUM" operator for type 'calcCurr'.
    call MPI_Op_create (sumAlldIdV, .true., MPIop, MPIerror)

!   Reduce the reults to the IOnode and write
    call MPI_Reduce (dIdV(1,1,1), buffdIdV,                             &
                     NTenerg_div*nspin*NIVP, MPIalldIdV, MPIop,         &
                     0, MPI_Comm_MyWorld, MPIerror)

    if (IOnode) then
       dIdV(:,:,:) = buffdIdV(:,:,:)
#elif defined MPI
    do n = 0,Nodes-1
       if (Node == n .and. IOnode) then
#endif
          do e = 1,NTenerg_div
             do s = 1,nspin

                Vbias = VInitial

                do v = 1,NIVP

!                  'Ei' and 'Vbias' in eV, and 'dIdV'
!                  in A/eV (from CODATA - 2012).
                   write (iuExVxdI, '(e17.8e3,e17.8e3,e17.8e3,'     //  &
                        'e17.8e3,e17.8e3,e17.8e3)')                     &
                        Ei(e)*13.60569253D0, Vbias*13.60569253D0,       &
                        dIdV(e,s,v)%el/13.60569253D0,                   &
                        dIdV(e,s,v)%symm/13.60569253D0,                 &
                        dIdV(e,s,v)%asymm/13.60569253D0,                &
                        (dIdV(e,s,v)%el + dIdV(e,s,v)%symm              &
                        + dIdV(e,s,v)%asymm)/13.60569253D0

                   Vbias = Vbias + dV

                enddo

                write (iuExVxdI, '(/)', advance='no')

             enddo
          enddo
#ifdef MASTER_SLAVE
    end if ! IOnode
#elif defined MPI
       elseif (Node == n) then
          call MPI_Send (Ei, NTenerg_div, MPI_Double_Precision,         &
                         0, 1, MPI_Comm_world, MPIerror)
          call MPI_Send (dIdV(1,1,1), NTenerg_div*nspin*NIVP,           &
                         MPIalldIdV, 0, 2, MPI_Comm_world, MPIerror)
       elseif (Node == 0) then
          call MPI_Recv (buffEn, NTenerg_div, MPI_Double_Precision,     &
                         n, 1, MPI_Comm_world, MPIstatus, MPIerror)
          call MPI_Recv (buffdIdV, NTenerg_div*nspin*NIVP, MPIalldIdV,  &
                         n, 2, MPI_Comm_world, MPIstatus, MPIerror)
       endif
       if (n /= 0) then
          call MPI_Barrier (MPI_Comm_world, MPIerror)
          if (Node == 0) then
             do e = 1,NTenerg_div
                do s = 1,nspin

                   Vbias = VInitial

                   do v = 1,NIVP

!                     'Ei' and 'Vbias' in eV, and 'dIdV'
!                     in A/eV (from CODATA - 2012).
                      write (iuExVxdI, '(e17.8e3,e17.8e3,e17.8e3,' //   &
                           'e17.8e3,e17.8e3,e17.8e3)')                  &
                           buffEn(e)*13.60569253D0, Vbias*13.60569253D0,&
                           buffdIdV(e,s,v)%el/13.60569253D0,            &
                           buffdIdV(e,s,v)%symm/13.60569253D0,          &
                           buffdIdV(e,s,v)%asymm/13.60569253D0,         &
                           (buffdIdV(e,s,v)%el + buffdIdV(e,s,v)%symm   &
                           + buffdIdV(e,s,v)%asymm)/13.60569253D0

                      Vbias = Vbias + dV

                   enddo

                   write (iuExVxdI, '(/)', advance='no')

                enddo
             enddo
          endif
       endif
    enddo ! n = 0,Nodes-1
#endif

!   Close files and free buffers memory.
    if (IOnode) then

       call IOclose (iuExVxdI)

       deallocate (buffEn)
       deallocate (buffdIdV)

#ifdef MASTER_SLAVE
    else
       deallocate (buffdIdV)
       call MPI_Op_free (MPIop, MPIerror)
#endif

    endif

!   Free MPI data type.
#ifdef MPI
    call MPI_Type_free (MPIalldIdV, MPIerror)
#endif


  end subroutine writedIdV


!  *******************************************************************  !
!                              writed2IdV2                              !
!  *******************************************************************  !
!  Description: writes calculated derivative of the differential        !
!  conductances ('d2I/dV2') to output file as follows                   !
!                                                                       !
!    - 'slabel_ExVxdI.d2IdV2' contains 'Ei Vbias d2IdV2el d2IdV2symm    !
!    d2IdV2asymm d2IdV2tot'                                             !
!                                                                       !
!  where 'd2IdV2tot = d2IdV2el+d2IdV2symm+d2IdV2asymm'.                 !
!                                                                       !
!  Written by Pedro Brandimarte, Jan 2014.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    January 2014                                    !
!  *********************** INPUT FROM MODULES ************************  !
!  logical IOnode              : True if it is the I/O node             !
!  integer Node                : Actual node (MPI_Comm_rank)            !
!  integer Nodes               : Total number of nodes (MPI_Comm_size)  !
!  integer nspin               : Number of spin components              !
!  integer label_length        : Length of system label                 !
!  character(label_length) slabel : System Label (for output files)     !
!  character(60) directory     : Working directory                      !
!  integer NIVP                : Number of bias potential points        !
!  real*8 VInitial             : Initial value of the bias potential    !
!  real*8 dV                   : Bias potential step                    !
!  integer NTenerg_div         : Number of energy grid points per node  !
!  real*8 Ei(NTenerg_div)      : Energy grid points                     !
!  TYPE(alldIdV) d2IdV2(NTenerg_div,nspin,NIVP) : [real] calculated     !
!                                                 d2I/dV2               !
!  *******************************************************************  !
  subroutine writed2IdV2

!
!   Modules
!
#ifdef MASTER_SLAVE
    use parallel,        only: IOnode, Node, Nodes, MPI_Comm_MyWorld
#else
    use parallel,        only: IOnode, Node, Nodes
#endif
    use idsrdr_options,  only: nspin, label_length, slabel, directory,  &
                               NIVP, VInitial, dV
    use idsrdr_engrid,   only: NTenerg_div, Ei
#ifdef MASTER_SLAVE
    use idsrdr_conduct,  only: alldIdV, d2IdV2, sumAlldIdV
#else
    use idsrdr_conduct,  only: alldIdV, d2IdV2
#endif
    use idsrdr_io,       only: IOassign, IOclose
    use idsrdr_string,   only: STRpaste

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: e, s, v, iuExVxd2I
    real(8) :: Vbias
    real(8), dimension(:), allocatable :: buffEn
    real(8), parameter :: Ry2 = 13.60569253D0*13.60569253D0
    TYPE(alldIdV), allocatable, dimension (:,:,:) :: buffd2IdV2
    character(len=16) :: suffix
    character(len=label_length+70) :: fExVxd2I

#ifdef MPI
#ifndef MASTER_SLAVE
    integer :: n
    integer, dimension(MPI_Status_Size) :: MPIstatus
#else
    integer :: MPIop
#endif
    integer :: MPIerror, MPIalld2IdV2
    integer :: blocklens(1) ! # of elements in each block
    integer :: blockdispl(1) ! byte displacement of each block
    integer :: oldtypes(1) ! type of elements in each block
#endif

    if (IOnode) then

       write (6, '(/,a)') 'Writing calculated differential '       //   &
            'conductance to *.d2IdV2 files'

!      Set file's names.
       suffix = '_ExVxd2I.d2IdV2'
       call STRpaste (slabel, suffix, fExVxd2I)
       call STRpaste (directory, fExVxd2I, fExVxd2I)

!      Open them.
       call IOassign (iuExVxd2I)
       open (iuExVxd2I, FILE=fExVxd2I, FORM='FORMATTED',                &
            STATUS='REPLACE')

!      Allocate buffers arrays.
       allocate (buffEn(NTenerg_div))
       allocate (buffd2IdV2(NTenerg_div,nspin,NIVP))

#ifdef MASTER_SLAVE
    else
       allocate (buffd2IdV2(NTenerg_div,nspin,NIVP))
#endif

    endif

#ifdef MPI
!   Set the description of 'alldIdV' type (3 doubles).
    blocklens(1) = 3
    blockdispl(1) = 0
    oldtypes(1) = MPI_Double_Precision

!   Define structured type and commit it.
    call MPI_Type_create_struct (1, blocklens, blockdispl, oldtypes,    &
                                 MPIalld2IdV2, MPIerror)
    call MPI_Type_commit (MPIalld2IdV2, MPIerror)
#endif

!   Write to the output files.
#ifdef MASTER_SLAVE
!   Define a "MPI_SUM" operator for type 'calcCurr'.
    call MPI_Op_create (sumAlldIdV, .true., MPIop, MPIerror)

!   Reduce the reults to the IOnode and write
    call MPI_Reduce (d2IdV2(1,1,1), buffd2IdV2,                         &
                     NTenerg_div*nspin*NIVP, MPIalld2IdV2, MPIop,       &
                     0, MPI_Comm_MyWorld, MPIerror)

    if (IOnode) then
       d2IdV2(:,:,:) = buffd2IdV2(:,:,:)
#elif defined MPI
    do n = 0,Nodes-1
       if (Node == n .and. IOnode) then
#endif
          do e = 1,NTenerg_div
             do s = 1,nspin

                Vbias = VInitial

                do v = 1,NIVP

!                  'Ei' and 'Vbias' in eV, and 'd2I/dV2'
!                  in A/V^2 (from CODATA - 2012).
                   write (iuExVxd2I, '(e17.8e3,e17.8e3,e17.8e3,'    //  &
                        'e17.8e3,e17.8e3,e17.8e3)') Ei(e)*13.60569253D0,&
                        Vbias*13.60569253D0, d2IdV2(e,s,v)%el/Ry2,      &
                        d2IdV2(e,s,v)%symm/Ry2, d2IdV2(e,s,v)%asymm/Ry2,&
                        (d2IdV2(e,s,v)%el + d2IdV2(e,s,v)%symm          &
                        + d2IdV2(e,s,v)%asymm)/Ry2

                   Vbias = Vbias + dV

                enddo

                write (iuExVxd2I, '(/)', advance='no')

             enddo
          enddo
#ifdef MASTER_SLAVE
    end if ! IOnode
#elif defined MPI
       elseif (Node == n) then
          call MPI_Send (Ei, NTenerg_div, MPI_Double_Precision,         &
                         0, 1, MPI_Comm_world, MPIerror)
          call MPI_Send (d2IdV2(1,1,1), NTenerg_div*nspin*NIVP,         &
                         MPIalld2IdV2, 0, 2, MPI_Comm_world, MPIerror)
       elseif (Node == 0) then
          call MPI_Recv (buffEn, NTenerg_div, MPI_Double_Precision,     &
                         n, 1, MPI_Comm_world, MPIstatus, MPIerror)
          call MPI_Recv (buffd2IdV2, NTenerg_div*nspin*NIVP,            &
                         MPIalld2IdV2, n, 2, MPI_Comm_world,            &
                         MPIstatus, MPIerror)
       endif
       if (n /= 0) then
          call MPI_Barrier (MPI_Comm_world, MPIerror)
          if (Node == 0) then
             do e = 1,NTenerg_div
                do s = 1,nspin

                   Vbias = VInitial

                   do v = 1,NIVP

!                     'Ei' and 'Vbias' in eV, and 'd2I/dV2'
!                     in A/V^2 (from CODATA - 2012).
                      write (iuExVxd2I,                                 &
                           '(e17.8e3,e17.8e3,e17.8e3,e17.8e3,'     //   &
                           'e17.8e3,e17.8e3)') buffEn(e)*13.60569253D0, &
                           Vbias*13.60569253D0,                         &
                           buffd2IdV2(e,s,v)%el/Ry2,                    &
                           buffd2IdV2(e,s,v)%symm/Ry2,                  &
                           buffd2IdV2(e,s,v)%asymm/Ry2,                 &
                           (buffd2IdV2(e,s,v)%el                        &
                           + buffd2IdV2(e,s,v)%symm                     &
                           + buffd2IdV2(e,s,v)%asymm)/Ry2

                      Vbias = Vbias + dV

                   enddo

                   write (iuExVxd2I, '(/)', advance='no')

                enddo
             enddo
          endif
       endif
    enddo ! n = 0,Nodes-1
#endif

!   Close files and free buffers memory.
    if (IOnode) then

       call IOclose (iuExVxd2I)

       deallocate (buffEn)
       deallocate (buffd2IdV2)

#ifdef MASTER_SLAVE
    else
       deallocate (buffd2IdV2)
       call MPI_Op_free (MPIop, MPIerror)
#endif

    endif

!   Free MPI data type.
#ifdef MPI
    call MPI_Type_free (MPIalld2IdV2, MPIerror)
#endif


  end subroutine writed2IdV2


!  *******************************************************************  !


END MODULE idsrdr_out

