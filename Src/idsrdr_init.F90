!  *******************************************************************  !
!  I-Disorder Fortran Code                                              !
!                                                                       !
!  Written by Alexandre Reily Rocha and Pedro Brandimarte, 2007-2013    !
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
!                          MODULE idsrdr_init                           !
!  *******************************************************************  !
!  Description: electron transport in disordered systems with           !
!  electron-phonon interaction.                                         !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !

MODULE idsrdr_init

!
! Modules
!
  use parallel,        only:
  use fdf

  implicit none

  PUBLIC ! default is public

! System Label (to name output files).
  integer, parameter :: label_length = 60
  character(len=label_length), save :: slabel


CONTAINS


!  *******************************************************************  !
!                                 init                                  !
!  *******************************************************************  !
!  Description: returns the address at 'B_Fdat' of a given 'key' or     !
!  returns 'KEY_MISSING' (through the global variable 'Koderr').        !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !
  subroutine init

!
! Modules
!
    use parallel,        only: Node, Nodes, IOnode

    include "mpif.h"

!   Local variables.
    external :: timer
#ifdef MPI
    integer :: MPIerror ! Return error code in MPI routines
#endif


!   Initialise MPI and set processor number.
#ifdef MPI
    call MPI_Init (MPIerror)
    call MPI_Comm_Rank (MPI_Comm_World, Node, MPIerror)
    call MPI_Comm_Size (MPI_Comm_World, Nodes, MPIerror)
#endif

    IOnode = Node .eq. 0

!   Print version information.
    if (IOnode) then
#ifdef MPI
       if (Nodes.gt.1) then
          write(6,'(/,a,i4,a)')                                         &
               '* Running on ', Nodes, ' nodes in parallel'
       else
          write(6,'(/,a,i4,a)') '* Running in serial mode with MPI'
       endif
#else
       write(6,'(/,a,i4,a)') '* Running in serial mode'
#endif
    endif

!   Start time counter.
    call timer ('i-disorder', 0)

!   Initialize some variables.
!!$    kpoint = 0.d0
!!$    nsc = 1
    slabel = ""

!   Initialise read.
    call initread

!   Read simulation data.
!!$    call read_options


  end subroutine init


!  *******************************************************************  !
!                               initread                                !
!  *******************************************************************  !
!  Description: initialise the reading of the data for I-Disorder.      !
!                                                                       !
!  Uses FDF (Flexible Data Fromat) package of J.M.Soler and A.Garcia.   !
!                                                                       !
!  Taken from 'reinit' subroutine from SIESTA package.                  !
!                                                                       !
!  Written by Pedro Brandimarte, Sep 2013.                              !
!  Instituto de Fisica                                                  !
!  Universidade de Sao Paulo                                            !
!  e-mail: brandimarte@gmail.com                                        !
!  ***************************** HISTORY *****************************  !
!  Original version:    September 2013                                  !
!  *******************************************************************  !
  subroutine initread

!
! Modules
!
    use parallel,        only: IOnode
    use fdf

    include "mpif.h"

!   Local variables.
    character string*20
    character filein*20, fileout*20, line*150 
    character slabel_default*59
    integer :: count, length, lun, lun_tmp
    logical :: debug_input, file_exists
    external :: io_assign, io_close
#ifdef MPI
    integer :: MPIerror
#endif

!   Print welcome and presentation.
    if (IOnode) then
       write (6,'(/a)')                                                 &
            '                           ***************************     '
       write (6,'(a)')                                                  &
            '                           *  WELCOME TO I-DISORDER  *     '
       write (6,'(a)')                                                  &
            '                           ***************************     '

!      Dump data file to output file and generate scratch file
!      for FDF to read from (except if INPUT_DEBUG exists).
       inquire (file='INPUT_DEBUG', exist=debug_input)
       if (debug_input) then
          write (6,'(a)') 'WARNING: ' //                                &
               'I-Disorder is reading its input from file INPUT_DEBUG'
           
          call io_assign(lun)
          filein = 'INPUT_DEBUG'
          open (lun, file='INPUT_DEBUG', form='formatted', status='old')
          rewind (lun)
       else
          write (6,'(/a)') 'initread: Reading from standard input'
          lun = 5
          call io_assign (lun_tmp)
          do
             call system_clock (count)
             write (string,*) count
             filein = 'INPUT_TMP.' // adjustl(string)
             inquire (file=filein, exist=file_exists)
             if (.not. file_exists) exit
          end do
          open (lun_tmp, file=filein, form='formatted', status='replace')
          rewind (lun_tmp)
       endif

       write (6,'(a,23(1h*),a,28(1h*))') '***',                         &
            ' Dump of input data file '

10     continue
       read (lun, err=20, end=20, fmt='(a)') line
       call chrlen (line, 0, length)
       if (length .ne. 0) then
          write(6,'(a)') line(1:length)
          if (.not. debug_input) write (lun_tmp,'(a)') line(1:length)
       endif
       goto 10
20     continue

       write (6,'(a,23(1h*),a,29(1h*))') '***',                         &
            ' End of input data file '

!      Choose proper file for fdf processing.
       if (debug_input) then
          call io_close (lun)
       else
          call io_close (lun_tmp)
       endif

!      Set up fdf.
       fileout = 'fdf.log'
       call fdf_init (filein, fileout)

!      Defile System Label (short name to label files).
       slabel_default = 'i-disorder'
       slabel = fdf_string ('SystemLabel', slabel_default)
!!$       label_length = size (slabel)
       write (6,'(a,a)') 'initread: System Label: ', slabel
       write (6,'(a,71(1h-))') 'initread: '

    endif

!   Global broadcast of system label.
#ifdef MPI
    call MPI_Bcast (slabel, label_length, MPI_Character, 0,             &
         MPI_Comm_World, MPIerror)
#endif


  end subroutine initread


!  *******************************************************************  !


END MODULE idsrdr_init
