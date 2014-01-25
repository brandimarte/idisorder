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

  implicit none

  PUBLIC  :: output
  PRIVATE :: writespectral, writecurrent


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
    real(8), parameter :: pi = 3.1415926535897932384626433832795028841D0
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

       endif

!      Write to the output files (energies in eV (from CODATA - 2012)).
#ifdef MPI
       do n = 0,Nodes-1
          if (Node == n .and. Node == 0) then
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
#ifdef MPI
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
       enddo
#endif

!      Close files and free buffers memory.
       if (Node == 0) then

          call IOclose (iuSpc)
          call IOclose (iuDos)

          deallocate (buffEn)
          deallocate (buffSpc)
          deallocate (buffDos)

       endif

    enddo ! do J = 1,nunitseph


  end subroutine writespectral


!  *******************************************************************  !
!                             writecurrent                              !
!  *******************************************************************  !
!  Description: writes calculated currents to output files as follows   !
!                                                                       !
!    - 'slabel_VxI.CUR' contains 'Vbias Iel Isymm Iasymm Itotal'        !
!    - 'slabel_ExI.CUR' contains 'Ei Iel Isymm Iasymm Itotal'           !
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
    integer :: i, n, e, s, iuExVxI, iuExVxIel, iuExVxIsymm,             &
               iuExVxIasymm, iuExVxItot
    real(8) :: Vbias
    real(8), parameter :: pi = 3.1415926535897932384626433832795028841D0
    real(8), dimension(:), allocatable :: buffEn
    TYPE(calcCurr), allocatable, dimension (:,:,:) :: buffCurr
    character(len=16) :: suffix
    character(len=label_length+70) :: fExVxI, fExVxIel, fExVxIsymm,     &
                                      fExVxIasymm, fExVxItot
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
       suffix = '_ExVxI.CUR'
       call STRpaste (slabel, suffix, fExVxI)
       call STRpaste (directory, fExVxI, fExVxI)
       suffix = '_ExVxIel.CUR'
       call STRpaste (slabel, suffix, fExVxIel)
       call STRpaste (directory, fExVxIel, fExVxIel)
       suffix = '_ExVxIsymm.CUR'
       call STRpaste (slabel, suffix, fExVxIsymm)
       call STRpaste (directory, fExVxIsymm, fExVxIsymm)
       suffix = '_ExVxIasymm.CUR'
       call STRpaste (slabel, suffix, fExVxIasymm)
       call STRpaste (directory, fExVxIasymm, fExVxIasymm)
       suffix = '_ExVxItot.CUR'
       call STRpaste (slabel, suffix, fExVxItot)
       call STRpaste (directory, fExVxItot, fExVxItot)

!      Open them.
       call IOassign (iuExVxI)
       open (iuExVxI, FILE=fExVxI, FORM='FORMATTED', STATUS='REPLACE')
       call IOassign (iuExVxIel)
       open (iuExVxIel, FILE=fExVxIel, FORM='FORMATTED',                &
            STATUS='REPLACE')
       call IOassign (iuExVxIsymm)
       open (iuExVxIsymm, FILE=fExVxIsymm, FORM='FORMATTED',            &
            STATUS='REPLACE')
       call IOassign (iuExVxIasymm)
       open (iuExVxIasymm, FILE=fExVxIasymm, FORM='FORMATTED',          &
            STATUS='REPLACE')
       call IOassign (iuExVxItot)
       open (iuExVxItot, FILE=fExVxItot, FORM='FORMATTED',              &
            STATUS='REPLACE')

!      Allocate buffers arrays.
       allocate (buffEn(NTenerg_div))
       allocate (buffCurr(NTenerg_div,nspin,NIVP))

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
#ifdef MPI
    do n = 0,Nodes-1
       if (Node == n .and. Node == 0) then
#endif
          do e = 1,NTenerg_div
             do s = 1,nspin

                Vbias = VInitial

                do i = 1,NIVP

!                  OBS.: Multiply 'Ei' and 'Vbias' by '13.60569253D0'
!                  for writing in eV (from CODATA - 2012).
                   write (iuExVxI,                                      &
                        '(e17.8e3,e17.8e3,e17.8e3,'                 //  &
                        'e17.8e3,e17.8e3,e17.8e3)') Ei(e), Vbias,       &
                        allcurr(e,s,i)%el, allcurr(e,s,i)%isymm,        &
                        allcurr(e,s,i)%iasymm, allcurr(e,s,i)%el        &
                        + allcurr(e,s,i)%isymm + allcurr(e,s,i)%iasymm
                   write (iuExVxIel,                                    &
                        '(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')            &
                        Ei(e), Vbias, allcurr(e,s,i)%el
                   write (iuExVxIsymm,                                  &
                        '(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')            &
                        Ei(e), Vbias, allcurr(e,s,i)%isymm
                   write (iuExVxIasymm,                                 &
                        '(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')            &
                        Ei(e), Vbias, allcurr(e,s,i)%iasymm
                   write (iuExVxItot,                                   &
                        '(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')            &
                        Ei(e), Vbias, allcurr(e,s,i)%el                 &
                        + allcurr(e,s,i)%isymm + allcurr(e,s,i)%iasymm

                   Vbias = Vbias + dV

                enddo
             enddo
          enddo
#ifdef MPI
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

                   do i = 1,NIVP

!                     OBS.: Multiply 'Ei' and 'Vbias' by '13.60569253D0'
!                     for writing in eV (from CODATA - 2012).
                      write (iuExVxI,                                   &
                           '(e17.8e3,e17.8e3,e17.8e3,e17.8e3,'    //    &
                           'e17.8e3,e17.8e3)') buffEn(e), Vbias,        &
                           buffCurr(e,s,i)%el, buffCurr(e,s,i)%isymm,   &
                           buffCurr(e,s,i)%iasymm, buffCurr(e,s,i)%el   &
                           + buffCurr(e,s,i)%isymm                      &
                           + buffCurr(e,s,i)%iasymm
                      write (iuExVxIel,                                 &
                           '(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')         &
                           buffEn(e), Vbias, buffCurr(e,s,i)%el
                      write (iuExVxIsymm,                               &
                           '(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')         &
                           buffEn(e), Vbias, buffCurr(e,s,i)%isymm
                      write (iuExVxIasymm,                              &
                           '(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')         &
                           buffEn(e), Vbias, buffCurr(e,s,i)%iasymm
                      write (iuExVxItot,                                &
                           '(e17.8e3,e17.8e3,e17.8e3,e17.8e3)')         &
                           buffEn(e), Vbias, buffCurr(e,s,i)%el         &
                           + buffCurr(e,s,i)%isymm                      &
                           + buffCurr(e,s,i)%iasymm

                      Vbias = Vbias + dV

                   enddo
                enddo
             enddo
          endif
       endif
    enddo
#endif

!   Close files and free buffers memory.
    if (Node == 0) then

       call IOclose (iuExVxI)
       call IOclose (iuExVxIel)
       call IOclose (iuExVxIsymm)
       call IOclose (iuExVxIasymm)
       call IOclose (iuExVxItot)

       deallocate (buffEn)
       deallocate (buffCurr)

    endif

!   Free MPI data type.
#ifdef MPI
    call MPI_Type_free (MPIcalcCurr, MPIerror)
#endif


  end subroutine writecurrent


!  *******************************************************************  !


END MODULE idsrdr_out

