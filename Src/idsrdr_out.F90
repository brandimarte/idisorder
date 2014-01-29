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
!  *******************************************************************  !
  subroutine output

!
!   Modules
!
    use parallel,        only: IOnode

    if (IOnode) write (6,'(/,28("*"),a,29("*"))')                       &
            ' Writing output files '

!   Write spectral function and DOS to files.
    call writespectral

!   Write calculated currents to files.
    call writecurrent

!   Write calculated dissipated powers to files.
    call writepower

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
    use parallel,        only: IOnode, Node, Nodes
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
    integer :: J, n, e, s, iuSpc, iuDos
!!$    real(8), parameter :: pi = 3.14159265358979323846264338327950241D0
    real(8), dimension(:), allocatable :: buffEn
    real(8), dimension(:,:), allocatable :: buffSpc, buffDos
    character(len=10) :: suffix
    character(len=label_length+70) :: fnSpc, fnDos
#ifdef MPI
    integer :: MPIerror
    integer, dimension(MPI_Status_Size) :: MPIstatus
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

!               OBS.: Multiply 'Ei' by '13.60569253D0'
!               for writing in eV (from CODATA - 2012).
                write (iuSpc, '(/,e17.8e3)', advance='no') Ei(e) !Ry
                write (iuDos, '(/,e17.8e3)', advance='no') Ei(e) !Ry

                do s = 1,nspin
!                  OBS.: Divide 'spctrl' and 'dos' by
!                  '(13.60569253D0 * pi' for writing in eV
!                  (from CODATA - 2012).
                   write (iuSpc, '(e17.8e3)', advance='no')             &
                        spctrl(e,s,J) !Ry
                   write (iuDos, '(e17.8e3)', advance='no')             &
                        dos(e,s,J) !Ry
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

!                  OBS.: Multiply 'Ei' by '13.60569253D0'
!                  for writing in eV (from CODATA - 2012).
                   write (iuSpc, '(/,e17.8e3)', advance='no')           &
                        buffEn(e) !Ry
                   write (iuDos, '(/,e17.8e3)', advance='no')           &
                        buffEn(e) !Ry

                   do s = 1,nspin

!                     OBS.: Divide 'spctrl' and 'dos' by
!                     '(13.60569253D0 * pi' for writing in eV
!                     (from CODATA - 2012).
                      write (iuSpc, '(e17.8e3)', advance='no')          &
                           buffSpc(e,s) !Ry
                      write (iuDos, '(e17.8e3)', advance='no')          &
                           buffDos(e,s) !Ry
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
!  Description: writes calculated currents to output files as follows   !
!                                                                       !
!    - 'slabel_VxI.CUR' contains 'Vbias Iel Isymm Iasymm Itotal'        !
!    - 'slabel_ExVxI.CUR' contains 'Ei Vbias Iel Isymm Iasymm Itotal'   !
!    - 'slabel_ExVxIel.CUR' contains 'Ei Vbias Iel'                     !
!    - 'slabel_ExVxIsymm.CUR' contains 'Ei Vbias Isymm'                 !
!    - 'slabel_ExVxIasymm.CUR' contains 'Ei Vbias Iasymm'               !
!    - 'slabel_ExVxItot.CUR' contains 'Ei Vbias Itotal'                 !
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
    use parallel,        only: IOnode, Node, Nodes
    use idsrdr_options,  only: nspin, label_length, slabel, directory,  &
                               NIVP, VInitial, dV
    use idsrdr_engrid,   only: NTenerg_div, Ei
    use idsrdr_current,  only: calcCurr, allcurr
    use idsrdr_io,       only: IOassign, IOclose
    use idsrdr_string,   only: STRpaste

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: n, e, s, v, iuVxI, iuExVxI, iuExVxIel, iuExVxIsy,        &
               iuExVxIasy, iuExVxItot
    real(8) :: Vbias
    real(8), dimension(:), allocatable :: buffEn
    TYPE(calcCurr), allocatable, dimension (:,:,:) :: buffCurr
    character(len=16) :: suffix
    character(len=label_length+70) :: fVxI, fExVxI, fExVxIel,           &
                                      fExVxIsy, fExVxIasy, fExVxItot
#ifdef MPI
    integer :: MPIerror, MPIcalcCurr
    integer, dimension(MPI_Status_Size) :: MPIstatus
    integer :: blocklens(1) ! # of elements in each block
    integer :: blockdispl(1) ! byte displacement of each block
    integer :: oldtypes(1) ! type of elements in each block
#endif

    if (IOnode) then

       write (6, '(/,a)') 'Writing calculated currents to *.CUR files'

!      Set file's names.
       suffix = '_VxI.CUR'
       call STRpaste (slabel, suffix, fVxI)
       call STRpaste (directory, fVxI, fVxI)
       suffix = '_ExVxI.CUR'
       call STRpaste (slabel, suffix, fExVxI)
       call STRpaste (directory, fExVxI, fExVxI)
       suffix = '_ExVxIel.CUR'
       call STRpaste (slabel, suffix, fExVxIel)
       call STRpaste (directory, fExVxIel, fExVxIel)
       suffix = '_ExVxIsy.CUR'
       call STRpaste (slabel, suffix, fExVxIsy)
       call STRpaste (directory, fExVxIsy, fExVxIsy)
       suffix = '_ExVxIasy.CUR'
       call STRpaste (slabel, suffix, fExVxIasy)
       call STRpaste (directory, fExVxIasy, fExVxIasy)
       suffix = '_ExVxItot.CUR'
       call STRpaste (slabel, suffix, fExVxItot)
       call STRpaste (directory, fExVxItot, fExVxItot)

!      Open them.
       call IOassign (iuVxI)
       open (iuVxI, FILE=fVxI, FORM='FORMATTED', STATUS='REPLACE')
       call IOassign (iuExVxI)
       open (iuExVxI, FILE=fExVxI, FORM='FORMATTED', STATUS='REPLACE')
       call IOassign (iuExVxIel)
       open (iuExVxIel, FILE=fExVxIel, FORM='FORMATTED',                &
            STATUS='REPLACE')
       call IOassign (iuExVxIsy)
       open (iuExVxIsy, FILE=fExVxIsy, FORM='FORMATTED',            &
            STATUS='REPLACE')
       call IOassign (iuExVxIasy)
       open (iuExVxIasy, FILE=fExVxIasy, FORM='FORMATTED',          &
            STATUS='REPLACE')
       call IOassign (iuExVxItot)
       open (iuExVxItot, FILE=fExVxItot, FORM='FORMATTED',              &
            STATUS='REPLACE')

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
!   Reduce the reults to the IOnode and write
    call MPI_Reduce (allcurr(1,1,1), buffCurr,                          &
                     NTenerg_div*nspin*NIVP, MPIcalcCurr, MPI_Sum,      &
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

!                  OBS.: Multiply 'Ei' and 'Vbias' by '13.60569253D0'
!                  for writing in eV (from CODATA - 2012).
                   write (iuVxI, '(e17.8e3,e17.8e3,e17.8e3,'        //  &
                        'e17.8e3,e17.8e3)') Vbias, allcurr(e,s,v)%el,   &
                        allcurr(e,s,v)%isymm, allcurr(e,s,v)%iasymm,    &
                        allcurr(e,s,v)%el + allcurr(e,s,v)%isymm +      &
                        allcurr(e,s,v)%iasymm
                   write (iuExVxI,                                      &
                        '(e17.8e3,e17.8e3,e17.8e3,'                 //  &
                        'e17.8e3,e17.8e3,e17.8e3)') Ei(e), Vbias,       &
                        allcurr(e,s,v)%el, allcurr(e,s,v)%isymm,        &
                        allcurr(e,s,v)%iasymm, allcurr(e,s,v)%el        &
                        + allcurr(e,s,v)%isymm + allcurr(e,s,v)%iasymm
                   write (iuExVxIel, '(e17.8e3,e17.8e3,e17.8e3)')       &
                        Ei(e), Vbias, allcurr(e,s,v)%el
                   write (iuExVxIsy, '(e17.8e3,e17.8e3,e17.8e3)')       &
                        Ei(e), Vbias, allcurr(e,s,v)%isymm
                   write (iuExVxIasy, '(e17.8e3,e17.8e3,e17.8e3)')      &
                        Ei(e), Vbias, allcurr(e,s,v)%iasymm
                   write (iuExVxItot, '(e17.8e3,e17.8e3,e17.8e3)')      &
                        Ei(e), Vbias, allcurr(e,s,v)%el                 &
                        + allcurr(e,s,v)%isymm + allcurr(e,s,v)%iasymm

                   Vbias = Vbias + dV

                enddo
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

!                     OBS.: Multiply 'Ei' and 'Vbias' by '13.60569253D0'
!                     for writing in eV (from CODATA - 2012).
                      write (iuVxI, '(e17.8e3,e17.8e3,e17.8e3,'     //  &
                           'e17.8e3,e17.8e3)') Vbias,                   &
                           buffCurr(e,s,v)%el, buffCurr(e,s,v)%isymm,   &
                           buffCurr(e,s,v)%iasymm, buffCurr(e,s,v)%el   &
                           + buffCurr(e,s,v)%isymm                      &
                           + buffCurr(e,s,v)%iasymm
                      write (iuExVxI,                                   &
                           '(e17.8e3,e17.8e3,e17.8e3,e17.8e3,'      //  &
                           'e17.8e3,e17.8e3)') buffEn(e), Vbias,        &
                           buffCurr(e,s,v)%el, buffCurr(e,s,v)%isymm,   &
                           buffCurr(e,s,v)%iasymm, buffCurr(e,s,v)%el   &
                           + buffCurr(e,s,v)%isymm                      &
                           + buffCurr(e,s,v)%iasymm
                      write (iuExVxIel, '(e17.8e3,e17.8e3,e17.8e3)')    &
                           buffEn(e), Vbias, buffCurr(e,s,v)%el
                      write (iuExVxIsy, '(e17.8e3,e17.8e3,e17.8e3)')    &
                           buffEn(e), Vbias, buffCurr(e,s,v)%isymm
                      write (iuExVxIasy, '(e17.8e3,e17.8e3,e17.8e3)')   &
                           buffEn(e), Vbias, buffCurr(e,s,v)%iasymm
                      write (iuExVxItot, '(e17.8e3,e17.8e3,e17.8e3)')   &
                           buffEn(e), Vbias, buffCurr(e,s,v)%el         &
                           + buffCurr(e,s,v)%isymm                      &
                           + buffCurr(e,s,v)%iasymm

                      Vbias = Vbias + dV

                   enddo
                enddo
             enddo
          endif
       endif
    enddo ! n = 0,Nodes-1
#endif

!   Close files and free buffers memory.
    if (IOnode) then

       call IOclose (iuVxI)
       call IOclose (iuExVxI)
       call IOclose (iuExVxIel)
       call IOclose (iuExVxIsy)
       call IOclose (iuExVxIasy)
       call IOclose (iuExVxItot)

       deallocate (buffEn)
       deallocate (buffCurr)

#ifdef MASTER_SLAVE
    else
       deallocate (buffCurr)
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
!  Description: writes calculated dissipated powers to output files     !
!  as follows                                                           !
!                                                                       !
!    - 'slabel_VxPtot.PWR' contains 'Vbias Power_tot'                   !
!    - 'slabel_ExPtot.PWR' contains 'Ei Power_tot'                      !
!    - 'slabel_ExVxPtot.PWR' contains 'Ei Vbias Power_tot'              !
!    - 'slabel_VxP.PWR' contains 'Vbias Power_mode'                     !
!    - 'slabel_EXP.PWR' contains 'Ei Power_mode'                        !
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
    use parallel,        only: IOnode, Node, Nodes
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
    integer :: n, e, s, v, u, w, idx,                                   &
               iuVxPtot, iuExPtot, iuExVxPtot, iuVxP, iuExP
    real(8) :: Vbias
    real(8), dimension(:), allocatable :: buffEn
    real(8), dimension (:,:,:,:), allocatable :: buffPwr
    character(len=8) :: suffix
    character(len=18) :: prefix
    character(len=label_length+70) :: fVxPtot, fExPtot, fExVxPtot,      &
                                      fVxP, fExP

#ifdef MPI
    integer :: MPIerror
    integer, dimension(MPI_Status_Size) :: MPIstatus
#endif

    do u = 1,nunitseph ! over units with e-ph intereaction

       idx = ephIdx(eph_type(u))

       if (IOnode) then

          write (6, '(/,a)') 'Writing calculated dissipated powers ' // &
               'to *.PWR files'

!         Set file's names.
          write (suffix,'(i3)') u
          call STRconcat (suffix, '.PWR', suffix)
          call STRpaste ('_VxPtot_', suffix, prefix)
          call STRpaste (slabel, prefix, fVxPtot)
          call STRpaste (directory, fVxPtot, fVxPtot)
          call STRpaste ('_ExPtot_', suffix, prefix)
          call STRpaste (slabel, prefix, fExPtot)
          call STRpaste (directory, fExPtot, fExPtot)
          call STRpaste ('_ExVxPtot_', suffix, prefix)
          call STRpaste (slabel, prefix, fExVxPtot)
          call STRpaste (directory, fExVxPtot, fExVxPtot)
          call STRpaste ('_VxP_', suffix, prefix)
          call STRpaste (slabel, prefix, fVxP)
          call STRpaste (directory, fVxP, fVxP)
          call STRpaste ('_ExP_', suffix, prefix)
          call STRpaste (slabel, prefix, fExP)
          call STRpaste (directory, fExP, fExP)

!         Open them.
          call IOassign (iuVxPtot)
          open (iuVxPtot, FILE=fVxPtot, FORM='FORMATTED',               &
               STATUS='REPLACE')
          call IOassign (iuExPtot)
          open (iuExPtot, FILE=fExPtot, FORM='FORMATTED',               &
               STATUS='REPLACE')
          call IOassign (iuExVxPtot)
          open (iuExVxPtot, FILE=fExVxPtot, FORM='FORMATTED',           &
               STATUS='REPLACE')
          call IOassign (iuVxP)
          open (iuVxP, FILE=fVxP, FORM='FORMATTED', STATUS='REPLACE')
          call IOassign (iuExP)
          open (iuExP, FILE=fExP, FORM='FORMATTED', STATUS='REPLACE')

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

!                     OBS.: Multiply 'Ei' and 'Vbias' by '13.60569253D0'
!                     for writing in eV (from CODATA - 2012).
                      write (iuVxPtot, '(e17.8e3,e17.8e3)') Vbias,      &
                           SUM(phPwr(u)%P(e,s,v,1:nModes(idx)))
                      write (iuExPtot, '(e17.8e3,e17.8e3)') Ei(e),      &
                           SUM(phPwr(u)%P(e,s,v,1:nModes(idx)))
                      write (iuExVxPtot, '(e17.8e3,e17.8e3,e17.8e3)')   &
                           Ei(e), Vbias,                                &
                           SUM(phPwr(u)%P(e,s,v,1:nModes(idx)))

                      write (iuVxP, '(/,e17.8e3)', advance='no') Vbias
                      write (iuExP, '(/,e17.8e3)', advance='no') Ei(e)
                      do w = 1,nModes(idx)
                         write (iuVxP, '(e17.8e3)', advance='no')       &
                              phPwr(u)%P(e,s,v,w)
                         write (iuExP, '(e17.8e3)', advance='no')       &
                              phPwr(u)%P(e,s,v,w)
                      enddo

                      Vbias = Vbias + dV

                   enddo
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

!                     OBS.: Multiply 'Ei' and 'Vbias' by '13.60569253D0'
!                     for writing in eV (from CODATA - 2012).
                         write (iuVxPtot, '(e17.8e3,e17.8e3)') Vbias,   &
                              SUM(buffPwr(e,s,v,1:nModes(idx)))
                         write (iuExPtot, '(e17.8e3,e17.8e3)')          &
                              buffEn(e),                                &
                              SUM(buffPwr(e,s,v,1:nModes(idx)))
                         write (iuExVxPtot, '(e17.8e3,e17.8e3,e17.8e3)')&
                              buffEn(e), Vbias,                         &
                              SUM(buffPwr(e,s,v,1:nModes(idx)))

                         write (iuVxP, '(/,e17.8e3)', advance='no') Vbias
                         write (iuExP, '(/,e17.8e3)', advance='no')     &
                              buffEn(e)
                         do w = 1,nModes(idx)
                            write (iuVxP, '(e17.8e3)', advance='no')    &
                                 buffPwr(e,s,v,w)
                            write (iuExP, '(e17.8e3)', advance='no')    &
                                 buffPwr(e,s,v,w)
                         enddo

                         Vbias = Vbias + dV

                      enddo
                   enddo
                enddo
             endif
          endif
       enddo ! n = 0,Nodes-1
#endif

!      Close files and free buffers memory.
       if (IOnode) then

          call IOclose (iuVxPtot)
          call IOclose (iuExPtot)
          call IOclose (iuExVxPtot)
          call IOclose (iuVxP)
          call IOclose (iuExP)

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
!  to output files as follows                                           !
!                                                                       !
!    - 'slabel_VxdI.dIdV' contains 'Vbias dIdVel dIdVsymm dIdVasymm     !
!    dIdVtotal'                                                         !
!    - 'slabel_ExVxdI.dIdV' contains 'Ei Vbias dIdVel dIdVsymm          !
!    dIdVasymm dIdVtotal'                                               !
!    - 'slabel_ExVxdIel.dIdV' contains 'Ei Vbias dIdVel'                !
!    - 'slabel_ExVxdIsy.dIdV' contains 'Ei Vbias dIdVsymm'              !
!    - 'slabel_ExVxdIasy.dIdV' contains 'Ei Vbias dIdVasymm'            !
!    - 'slabel_ExVxdItot.dIdV' contains 'Ei Vbias dIdVtotal'            !
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
    use parallel,        only: IOnode, Node, Nodes
    use idsrdr_options,  only: nspin, label_length, slabel, directory,  &
                               NIVP, VInitial, dV
    use idsrdr_engrid,   only: NTenerg_div, Ei
    use idsrdr_conduct,  only: alldIdV, dIdV
    use idsrdr_io,       only: IOassign, IOclose
    use idsrdr_string,   only: STRpaste

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: n, e, s, v, iuVxdI, iuExVxdI, iuExVxdIel, iuExVxdIsy,    &
               iuExVxdIasy, iuExVxdItot
    real(8) :: Vbias
    real(8), dimension(:), allocatable :: buffEn
    TYPE(alldIdV), allocatable, dimension (:,:,:) :: buffdIdV
    character(len=16) :: suffix
    character(len=label_length+70) :: fVxdI, fExVxdI, fExVxdIel,        &
                                      fExVxdIsy, fExVxdIasy, fExVxdItot
#ifdef MPI
    integer :: MPIerror, MPIalldIdV
    integer, dimension(MPI_Status_Size) :: MPIstatus
    integer :: blocklens(1) ! # of elements in each block
    integer :: blockdispl(1) ! byte displacement of each block
    integer :: oldtypes(1) ! type of elements in each block
#endif

    if (IOnode) then

       write (6, '(/,a)') 'Writing calculated differential '       //   &
            'conductance to *.dIdV files'

!      Set file's names.
       suffix = '_VxdI.dIdV'
       call STRpaste (slabel, suffix, fVxdI)
       call STRpaste (directory, fVxdI, fVxdI)
       suffix = '_ExVxdI.dIdV'
       call STRpaste (slabel, suffix, fExVxdI)
       call STRpaste (directory, fExVxdI, fExVxdI)
       suffix = '_ExVxdIel.dIdV'
       call STRpaste (slabel, suffix, fExVxdIel)
       call STRpaste (directory, fExVxdIel, fExVxdIel)
       suffix = '_ExVxdIsy.dIdV'
       call STRpaste (slabel, suffix, fExVxdIsy)
       call STRpaste (directory, fExVxdIsy, fExVxdIsy)
       suffix = '_ExVxdIasy.dIdV'
       call STRpaste (slabel, suffix, fExVxdIasy)
       call STRpaste (directory, fExVxdIasy, fExVxdIasy)
       suffix = '_ExVxdItot.dIdV'
       call STRpaste (slabel, suffix, fExVxdItot)
       call STRpaste (directory, fExVxdItot, fExVxdItot)

!      Open them.
       call IOassign (iuVxdI)
       open (iuVxdI, FILE=fVxdI, FORM='FORMATTED', STATUS='REPLACE')
       call IOassign (iuExVxdI)
       open (iuExVxdI, FILE=fExVxdI, FORM='FORMATTED', STATUS='REPLACE')
       call IOassign (iuExVxdIel)
       open (iuExVxdIel, FILE=fExVxdIel, FORM='FORMATTED',              &
            STATUS='REPLACE')
       call IOassign (iuExVxdIsy)
       open (iuExVxdIsy, FILE=fExVxdIsy, FORM='FORMATTED',              &
            STATUS='REPLACE')
       call IOassign (iuExVxdIasy)
       open (iuExVxdIasy, FILE=fExVxdIasy, FORM='FORMATTED',            &
            STATUS='REPLACE')
       call IOassign (iuExVxdItot)
       open (iuExVxdItot, FILE=fExVxdItot, FORM='FORMATTED',            &
            STATUS='REPLACE')

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
!   Reduce the reults to the IOnode and write
    call MPI_Reduce (dIdV(1,1,1), buffdIdV,                             &
                     NTenerg_div*nspin*NIVP, MPIalldIdV, MPI_Sum,       &
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

!                  OBS.: Multiply 'Ei' and 'Vbias' by '13.60569253D0'
!                  for writing in eV (from CODATA - 2012).
                   write (iuVxdI, '(e17.8e3,e17.8e3,e17.8e3,'       //  &
                        'e17.8e3,e17.8e3)') Vbias, dIdV(e,s,v)%el,      &
                        dIdV(e,s,v)%symm, dIdV(e,s,v)%asymm,            &
                        dIdV(e,s,v)%el + dIdV(e,s,v)%symm               &
                        + dIdV(e,s,v)%asymm
                   write (iuExVxdI,                                     &
                        '(e17.8e3,e17.8e3,e17.8e3,'                 //  &
                        'e17.8e3,e17.8e3,e17.8e3)') Ei(e), Vbias,       &
                        dIdV(e,s,v)%el, dIdV(e,s,v)%symm,               &
                        dIdV(e,s,v)%asymm, dIdV(e,s,v)%el               &
                        + dIdV(e,s,v)%symm + dIdV(e,s,v)%asymm
                   write (iuExVxdIel, '(e17.8e3,e17.8e3,e17.8e3)')      &
                        Ei(e), Vbias, dIdV(e,s,v)%el
                   write (iuExVxdIsy, '(e17.8e3,e17.8e3,e17.8e3)')      &
                        Ei(e), Vbias, dIdV(e,s,v)%symm
                   write (iuExVxdIasy, '(e17.8e3,e17.8e3,e17.8e3)')     &
                        Ei(e), Vbias, dIdV(e,s,v)%asymm
                   write (iuExVxdItot, '(e17.8e3,e17.8e3,e17.8e3)')     &
                        Ei(e), Vbias, dIdV(e,s,v)%el                    &
                        + dIdV(e,s,v)%symm + dIdV(e,s,v)%asymm

                   Vbias = Vbias + dV

                enddo
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

!                     OBS.: Multiply 'Ei' and 'Vbias' by '13.60569253D0'
!                     for writing in eV (from CODATA - 2012).
                      write (iuVxdI, '(e17.8e3,e17.8e3,e17.8e3,'    //  &
                           'e17.8e3,e17.8e3)') Vbias,                   &
                           buffdIdV(e,s,v)%el, buffdIdV(e,s,v)%symm,    &
                           buffdIdV(e,s,v)%asymm, buffdIdV(e,s,v)%el    &
                           + buffdIdV(e,s,v)%symm + buffdIdV(e,s,v)%asymm
                      write (iuExVxdI,                                  &
                           '(e17.8e3,e17.8e3,e17.8e3,e17.8e3,'     //   &
                           'e17.8e3,e17.8e3)') buffEn(e), Vbias,        &
                           buffdIdV(e,s,v)%el, buffdIdV(e,s,v)%symm,    &
                           buffdIdV(e,s,v)%asymm, buffdIdV(e,s,v)%el    &
                           + buffdIdV(e,s,v)%symm + buffdIdV(e,s,v)%asymm
                      write (iuExVxdIel, '(e17.8e3,e17.8e3,e17.8e3)')   &
                           buffEn(e), Vbias, buffdIdV(e,s,v)%el
                      write (iuExVxdIsy, '(e17.8e3,e17.8e3,e17.8e3)')   &
                           buffEn(e), Vbias, buffdIdV(e,s,v)%symm
                      write (iuExVxdIasy, '(e17.8e3,e17.8e3,e17.8e3)')  &
                           buffEn(e), Vbias, buffdIdV(e,s,v)%asymm
                      write (iuExVxdItot, '(e17.8e3,e17.8e3,e17.8e3)')  &
                           buffEn(e), Vbias, buffdIdV(e,s,v)%el         &
                           + buffdIdV(e,s,v)%symm + buffdIdV(e,s,v)%asymm

                      Vbias = Vbias + dV

                   enddo
                enddo
             enddo
          endif
       endif
    enddo ! n = 0,Nodes-1
#endif

!   Close files and free buffers memory.
    if (IOnode) then

       call IOclose (iuVxdI)
       call IOclose (iuExVxdI)
       call IOclose (iuExVxdIel)
       call IOclose (iuExVxdIsy)
       call IOclose (iuExVxdIasy)
       call IOclose (iuExVxdItot)

       deallocate (buffEn)
       deallocate (buffdIdV)

#ifdef MASTER_SLAVE
    else
       deallocate (buffdIdV)
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
!  conductances ('d2I/dV2') to output files as follows                  !
!                                                                       !
!    - 'slabel_VxdI.d2IdV2' contains 'Vbias d2IdV2el d2IdV2symm         !
!    d2IdV2asymm d2IdV2total'                                           !
!    - 'slabel_ExVxdI.d2IdV2' contains 'Ei Vbias d2IdV2el d2IdV2symm    !
!    d2IdV2asymm d2IdV2total'                                           !
!    - 'slabel_ExVxdIel.d2IdV2' contains 'Ei Vbias d2IdV2el'            !
!    - 'slabel_ExVxdIsy.d2IdV2' contains 'Ei Vbias d2IdV2symm'          !
!    - 'slabel_ExVxdIasy.d2IdV2' contains 'Ei Vbias d2IdV2asymm'        !
!    - 'slabel_ExVxdItot.d2IdV2' contains 'Ei Vbias d2IdV2total'        !
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
    use parallel,        only: IOnode, Node, Nodes
    use idsrdr_options,  only: nspin, label_length, slabel, directory,  &
                               NIVP, VInitial, dV
    use idsrdr_engrid,   only: NTenerg_div, Ei
    use idsrdr_conduct,  only: alldIdV, d2IdV2
    use idsrdr_io,       only: IOassign, IOclose
    use idsrdr_string,   only: STRpaste

#ifdef MPI
    include "mpif.h"
#endif

!   Local variables.
    integer :: n, e, s, v, iuVxd2I, iuExVxd2I, iuExVxd2Iel,             &
               iuExVxd2Isy, iuExVxd2Iasy, iuExVxd2Itot
    real(8) :: Vbias
    real(8), dimension(:), allocatable :: buffEn
    TYPE(alldIdV), allocatable, dimension (:,:,:) :: buffd2IdV2
    character(len=19) :: suffix
    character(len=label_length+70) :: fVxd2I, fExVxd2I, fExVxd2Iel,     &
                                      fExVxd2Isy, fExVxd2Iasy,          &
                                      fExVxd2Itot
#ifdef MPI
    integer :: MPIerror, MPIalld2IdV2
    integer, dimension(MPI_Status_Size) :: MPIstatus
    integer :: blocklens(1) ! # of elements in each block
    integer :: blockdispl(1) ! byte displacement of each block
    integer :: oldtypes(1) ! type of elements in each block
#endif

    if (IOnode) then

       write (6, '(/,a)') 'Writing calculated differential '       //   &
            'conductance to *.d2IdV2 files'

!      Set file's names.
       suffix = '_Vxd2I.d2IdV2'
       call STRpaste (slabel, suffix, fVxd2I)
       call STRpaste (directory, fVxd2I, fVxd2I)
       suffix = '_ExVxd2I.d2IdV2'
       call STRpaste (slabel, suffix, fExVxd2I)
       call STRpaste (directory, fExVxd2I, fExVxd2I)
       suffix = '_ExVxd2Iel.d2IdV2'
       call STRpaste (slabel, suffix, fExVxd2Iel)
       call STRpaste (directory, fExVxd2Iel, fExVxd2Iel)
       suffix = '_ExVxd2Isy.d2IdV2'
       call STRpaste (slabel, suffix, fExVxd2Isy)
       call STRpaste (directory, fExVxd2Isy, fExVxd2Isy)
       suffix = '_ExVxd2Iasy.d2IdV2'
       call STRpaste (slabel, suffix, fExVxd2Iasy)
       call STRpaste (directory, fExVxd2Iasy, fExVxd2Iasy)
       suffix = '_ExVxd2Itot.d2IdV2'
       call STRpaste (slabel, suffix, fExVxd2Itot)
       call STRpaste (directory, fExVxd2Itot, fExVxd2Itot)

!      Open them.
       call IOassign (iuVxd2I)
       open (iuVxd2I, FILE=fVxd2I, FORM='FORMATTED', STATUS='REPLACE')
       call IOassign (iuExVxd2I)
       open (iuExVxd2I, FILE=fExVxd2I, FORM='FORMATTED',                &
            STATUS='REPLACE')
       call IOassign (iuExVxd2Iel)
       open (iuExVxd2Iel, FILE=fExVxd2Iel, FORM='FORMATTED',            &
            STATUS='REPLACE')
       call IOassign (iuExVxd2Isy)
       open (iuExVxd2Isy, FILE=fExVxd2Isy, FORM='FORMATTED',            &
            STATUS='REPLACE')
       call IOassign (iuExVxd2Iasy)
       open (iuExVxd2Iasy, FILE=fExVxd2Iasy, FORM='FORMATTED',          &
            STATUS='REPLACE')
       call IOassign (iuExVxd2Itot)
       open (iuExVxd2Itot, FILE=fExVxd2Itot, FORM='FORMATTED',          &
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
!   Reduce the reults to the IOnode and write
    call MPI_Reduce (d2IdV2(1,1,1), buffd2IdV2,                         &
                     NTenerg_div*nspin*NIVP, MPIalld2IdV2, MPI_Sum,     &
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

!                  OBS.: Multiply 'Ei' and 'Vbias' by '13.60569253D0'
!                  for writing in eV (from CODATA - 2012).
                   write (iuVxd2I, '(e17.8e3,e17.8e3,e17.8e3,'      //  &
                        'e17.8e3,e17.8e3)') Vbias, d2IdV2(e,s,v)%el,    &
                        d2IdV2(e,s,v)%symm, d2IdV2(e,s,v)%asymm,        &
                        d2IdV2(e,s,v)%el + d2IdV2(e,s,v)%symm           &
                        + d2IdV2(e,s,v)%asymm
                   write (iuExVxd2I,                                    &
                        '(e17.8e3,e17.8e3,e17.8e3,'                 //  &
                        'e17.8e3,e17.8e3,e17.8e3)') Ei(e), Vbias,       &
                        d2IdV2(e,s,v)%el, d2IdV2(e,s,v)%symm,           &
                        d2IdV2(e,s,v)%asymm, d2IdV2(e,s,v)%el           &
                        + d2IdV2(e,s,v)%symm + d2IdV2(e,s,v)%asymm
                   write (iuExVxd2Iel, '(e17.8e3,e17.8e3,e17.8e3)')     &
                        Ei(e), Vbias, d2IdV2(e,s,v)%el
                   write (iuExVxd2Isy, '(e17.8e3,e17.8e3,e17.8e3)')     &
                        Ei(e), Vbias, d2IdV2(e,s,v)%symm
                   write (iuExVxd2Iasy, '(e17.8e3,e17.8e3,e17.8e3)')    &
                        Ei(e), Vbias, d2IdV2(e,s,v)%asymm
                   write (iuExVxd2Itot, '(e17.8e3,e17.8e3,e17.8e3)')    &
                        Ei(e), Vbias, d2IdV2(e,s,v)%el                  &
                        + d2IdV2(e,s,v)%symm + d2IdV2(e,s,v)%asymm

                   Vbias = Vbias + dV

                enddo
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

!                     OBS.: Multiply 'Ei' and 'Vbias' by '13.60569253D0'
!                     for writing in eV (from CODATA - 2012).
                      write (iuVxd2I, '(e17.8e3,e17.8e3,e17.8e3,'   //  &
                           'e17.8e3,e17.8e3)') Vbias,                   &
                           buffd2IdV2(e,s,v)%el,                        &
                           buffd2IdV2(e,s,v)%symm,                      &
                           buffd2IdV2(e,s,v)%asymm,                     &
                           buffd2IdV2(e,s,v)%el                         &
                           + buffd2IdV2(e,s,v)%symm                     &
                           + buffd2IdV2(e,s,v)%asymm
                      write (iuExVxd2I,                                 &
                           '(e17.8e3,e17.8e3,e17.8e3,e17.8e3,'     //   &
                           'e17.8e3,e17.8e3)') buffEn(e), Vbias,        &
                           buffd2IdV2(e,s,v)%el,                        &
                           buffd2IdV2(e,s,v)%symm,                      &
                           buffd2IdV2(e,s,v)%asymm,                     &
                           buffd2IdV2(e,s,v)%el                         &
                           + buffd2IdV2(e,s,v)%symm                     &
                           + buffd2IdV2(e,s,v)%asymm
                      write (iuExVxd2Iel, '(e17.8e3,e17.8e3,e17.8e3)')  &
                           buffEn(e), Vbias, buffd2IdV2(e,s,v)%el
                      write (iuExVxd2Isy, '(e17.8e3,e17.8e3,e17.8e3)')  &
                           buffEn(e), Vbias, buffd2IdV2(e,s,v)%symm
                      write (iuExVxd2Iasy, '(e17.8e3,e17.8e3,e17.8e3)') &
                           buffEn(e), Vbias, buffd2IdV2(e,s,v)%asymm
                      write (iuExVxd2Itot, '(e17.8e3,e17.8e3,e17.8e3)') &
                           buffEn(e), Vbias, buffd2IdV2(e,s,v)%el       &
                           + buffd2IdV2(e,s,v)%symm                     &
                           + buffd2IdV2(e,s,v)%asymm

                      Vbias = Vbias + dV

                   enddo
                enddo
             enddo
          endif
       endif
    enddo ! n = 0,Nodes-1
#endif

!   Close files and free buffers memory.
    if (IOnode) then

       call IOclose (iuVxd2I)
       call IOclose (iuExVxd2I)
       call IOclose (iuExVxd2Iel)
       call IOclose (iuExVxd2Isy)
       call IOclose (iuExVxd2Iasy)
       call IOclose (iuExVxd2Itot)

       deallocate (buffEn)
       deallocate (buffd2IdV2)

#ifdef MASTER_SLAVE
    else
       deallocate (buffd2IdV2)
#endif

    endif

!   Free MPI data type.
#ifdef MPI
    call MPI_Type_free (MPIalld2IdV2, MPIerror)
#endif


  end subroutine writed2IdV2


!  *******************************************************************  !


END MODULE idsrdr_out

